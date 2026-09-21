#!/usr/bin/env python3
"""Prepare member-0 S/Pi scan cards in an available export; never integrate.

The command-line entry point holds the shared MG5 lock. A future serialized
launcher may call configure while holding that lock itself. Tests use private
card/topology fixtures and do not acquire or modify a physical export.
"""
import argparse
import copy
import contextlib
import fcntl
import io
import json
import math
from pathlib import Path
import shutil

from campaign import ROOT, STUDY, GF, MW, MZ, digest, now, save, setup_module
from parameter_variation_inputs import validate


def write_parameters(process, scenario, benchmark, width):
    from models.check_param_card import ParamCard
    from models import import_ufo, model_reader
    path = process/'Cards/param_card.dat'
    card = ParamCard(str(path))
    if 'decaymass' in card:
        raise ValueError('This bounded parameter scan requires massless decays')
    alpha = next(row['alpha_s'] for row in scenario['coupling_checkpoints'] if row['Q_GeV'] == MZ)
    mt = scenario['mt_GeV']
    parameters = {('sminputs', 1): benchmark['inputs']['alpha_Gmu_inverse'],
        ('sminputs', 2): GF, ('sminputs', 3): alpha,
        ('mass', 6): mt, ('yukawa', 6): mt, ('mass', 23): MZ,
        ('mass', 5): 0., ('yukawa', 5): 0., ('mass', 25): 125., ('mass', 15): 1.777,
        ('decay', 6): width['gamma_nlo_pdf'], ('decay', 24): scenario['fixed_W_width_GeV'],
        ('decay', 23): 2.4952, ('decay', 25): .00407, ('decay', 15): 0.}
    for (block, code), value in parameters.items():
        try:
            card[block].get((code,)).value = value
        except KeyError:
            if value != 0.:
                raise
    model = model_reader.ModelReader(import_ufo.import_model('loop_sm-no_b_mass'))
    card.update_dependent(model, None, 20)
    model.set_parameters_and_couplings(card)
    if (abs(float(card['mass'].get((24,)).value)-MW) > 1.e-10
            or abs(model.get_mass(24)-MW) > 1.e-10
            or model.get_mass(6) != mt
            or float(card['yukawa'].get((6,)).value) != mt):
        raise ValueError('Dependent model masses or top Yukawa disagree with the scan')
    card.write(str(path), precision=16)
    return dict(mt_GeV=mt, top_yukawa_GeV=mt, model_internal_MW_GeV=model.get_mass(24),
                alpha_s_at_MZ=alpha, fixed_W_width_GeV=scenario['fixed_W_width_GeV'],
                top_width_lo_GeV=width['gamma_lo'], top_width_nlo_GeV=width['gamma_nlo_pdf'])


def configure(process, inputs_path, scenario_name, variant, seed, mode='onshell', accuracy=.015):
    from madgraph.various.banner import RunCardNLO
    process, inputs_path = Path(process).resolve(), Path(inputs_path).resolve()
    record = json.loads(inputs_path.read_text())
    if record.get('status') != 'parameter inputs verified; no variation integrations':
        raise ValueError('Require verified parameter-variation inputs')
    benchmark, _ = validate(record)
    matches = [row for row in record['scenarios'] if row['name'] == scenario_name]
    if len(matches) != 1 or variant not in ('S', 'Pi') or mode not in ('onshell', 'top-bw', 'all-bw'):
        raise ValueError('Unknown parameter scenario or S/Pi configuration')
    if not isinstance(seed, int) or not 0 < seed < 900000000 or not math.isfinite(accuracy) or not 0. < accuracy < 1.:
        raise ValueError('Invalid scan seed or initial accuracy')
    scenario = copy.deepcopy(matches[0])
    setup = setup_module()
    if setup.export_w_treatment(process)[0] != mode:
        raise ValueError('Export topology differs from requested W treatment')
    source = process/'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90'
    if digest(source) != digest(ROOT/'Template/fNLO/FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90'):
        raise ValueError('Parameter comparisons require the unchanged benchmark analysis')
    width = next(row for row in scenario['top_widths'] if row['scale_factor'] == 1.
                 and row['w_treatment'] == ('onshell' if mode == 'onshell' else 'bw'))
    sampling = 'w-current' if mode == 'all-bw' else 'flat'
    base_name = '%s_%s_core-w-ht-half_separate_%d' % (variant, mode, seed)
    if sampling == 'w-current':
        base_name += '_wcurrent'
    name = base_name+'__'+scenario_name+'_member0'
    archive, base = (process/'study_cards'/value for value in (name, base_name))
    if any(path.exists() for path in (archive, base, process/'Events'/name, process/'Events'/base_name)):
        raise ValueError('Existing scan cards or run; never overwrite or restart implicitly')
    archive.mkdir(parents=True)
    shutil.copy2(process/'Cards/param_card.dat', archive/'param_card.before.dat')
    physical = write_parameters(process, scenario, benchmark, width)
    setup_stdout = io.StringIO()
    with contextlib.redirect_stdout(setup_stdout):
        setup.configure(argparse.Namespace(process_dir=str(process), variant=variant,
            top_width_lo=width['gamma_lo'], top_width_nlo=width['gamma_nlo_pdf'],
            w_treatment=mode, top_width_w_treatment=width['w_treatment'],
            decay_bottom_mass=0., top_width_bottom_mass=0.,
            width_source=str(inputs_path)+'; scenario='+scenario_name+'; sha256='+digest(inputs_path),
            pdf_id=scenario['PDF']['id'], production_scale='core-w-ht-half', production_sampling=sampling,
            decay_scales='separate', ecm=scenario['ecm_GeV'], seed=seed,
            points=1000, grid_points=1000, iterations=3, accuracy=accuracy, job_seconds=1.))
    (archive/'base_setup_stdout.txt').write_text(setup_stdout.getvalue())
    # Family spreads are central-member cross-checks, with no Hessian/replica
    # error inferred from these runs. Keep every signed physical scale weight.
    run = RunCardNLO(str(process/'Cards/run_card.dat'))
    run.set('reweight_pdf', False, user=True, raiseerror=True)
    run.set('qes_over_ref', 1., user=True, raiseerror=True)
    run.write(str(process/'Cards/run_card.dat'), template=str(process/'Cards/run_card_default.dat'))
    manifest = json.loads((base/'manifest.json').read_text())
    manifest['settings'].update(reweight_pdf=False, qes_over_ref=1.)
    for filename in ('param_card.dat', 'run_card.dat', 'decay_card.dat', 'FO_analyse_card.dat'):
        shutil.copy2(process/'Cards'/filename, archive/filename)
        manifest['hashes'][filename] = digest(archive/filename)
    manifest.update(run_name=name, parameter_scenario=scenario,
        parameter_variation_inputs=str(inputs_path), parameter_variation_inputs_sha256=digest(inputs_path),
        physical_parameters=physical, pdf_error_weights=False,
        retained_weights='Nominal member 0 and all 81 signed physical scale points; no PDF uncertainty ensemble.',
        base_configuration=str(base/'manifest.json'), base_configuration_sha256=digest(base/'manifest.json'),
        prior_param_card_sha256=digest(archive/'param_card.before.dat'),
        script_sha256=digest(Path(__file__)), configured_utc=now(),
        status='parameter-scan cards configured; no integration performed')
    save(archive/'manifest.json', manifest)
    return name, archive


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('process', type=Path)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--scenario', required=True)
    parser.add_argument('--variant', choices=('S', 'Pi'), required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--mode', choices=('onshell', 'top-bw', 'all-bw'), default='onshell')
    parser.add_argument('--accuracy', type=float, default=.015)
    args = parser.parse_args()
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        name, archive = configure(args.process, args.inputs, args.scenario, args.variant,
                                  args.seed, args.mode, args.accuracy)
    print('Final parameter-scan cards:', archive)
    print('No integration launched; main S/Pi convergence must determine the scan allocation.')
