#!/usr/bin/env python3
"""Prepare ttW studies with on-shell tops and on-shell/finite-width Ws; never run."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

VARIANTS = {
    'LO': ('LO', 'LO', 'ADDITIVE'),
    'P': ('NLO', 'LO', 'ADDITIVE'),
    'D': ('LO', 'NLO', 'ADDITIVE'),
    'S': ('NLO', 'NLO', 'ADDITIVE'),
    'PiD': ('LO', 'NLO', 'MULTIPLICATIVE'),
    'Pi': ('NLO', 'NLO', 'MULTIPLICATIVE'),
}
W_TREATMENTS = ('onshell', 'top-bw', 'all-bw')
PRODUCTION_SCALES = ('core-w-ht-half', 'core-ht-half', 'fixed')


def process_commands(charge, flavours, output, corrected='both', real_only=False,
                     partonic=False, w_treatment='onshell', decay_bottom_mass=0.):
    """One export per charge and ordered (Wt,Wtbar,Wassoc) flavour assignment."""
    output = str(Path(output).expanduser().absolute())
    if not re.fullmatch(r'[\w./+\-]+', output):
        raise ValueError('MadGraph output path must contain no spaces or shell metacharacters')
    if w_treatment not in W_TREATMENTS:
        raise ValueError('Unknown W treatment: ' + w_treatment)
    if not math.isfinite(decay_bottom_mass) or decay_bottom_mass < 0.:
        raise ValueError('Decay bottom mass must be finite and nonnegative')
    lp = {'e': 'e+ ve', 'mu': 'mu+ vm'}
    lm = {'e': 'e- ve~', 'mu': 'mu- vm~'}
    order = '[real=QCD]' if real_only else '[QCD]'
    torder = ' ' + order if corrected in ('both', 't') else ''
    aorder = ' ' + order if corrected in ('both', 'tbar') else ''
    w = 'w+' if charge == 'plus' else 'w-'
    initial = ('u d~' if charge == 'plus' else 'd u~') if partonic else 'p p'
    assoc = (lp if charge == 'plus' else lm)[flavours[2]]
    if w_treatment == 'onshell':
        decays = ('(t > w+ b QED=1%s, w+ > %s), '
                  '(t~ > w- b~ QED=1%s, w- > %s)' %
                  (torder, lp[flavours[0]], aorder, lm[flavours[1]]))
    else:
        # The W is internal to each complete three-body decay block. Do not
        # introduce a separate W decay node or an additional branching factor.
        decays = ('(t > b %s QED=2%s), (t~ > b~ %s QED=2%s)' %
                  (lp[flavours[0]], torder, lm[flavours[1]], aorder))
    core = '%s QCD=2 QED=1' % w
    if w_treatment == 'all-bw':
        # No W resonance selector: integrate the lepton-pair virtuality over
        # the physical phase space, not an on-shell/resonance-window decay.
        core = '%s QCD=2 QED=2' % assoc
    else:
        decays += ', %s > %s' % (w, assoc)
    return '\n'.join([
        '# ttW W treatment: ' + w_treatment + '; tops remain on shell',
        'import model loop_sm-no_b_mass',
        'set decay_bottom_mass %.16g' % decay_bottom_mass,
        'define p = g u c d s b u~ c~ d~ s~ b~',
        'generate %s > t t~ %s %s, %s' % (initial, core, order, decays),
        'output fNLO ' + output,
        '',
    ])


def export_decay_bottom_mass(process, param):
    """Read the decay-only input, checking the generated model separation."""
    path = process / 'Cards/decay_mass_scheme.json'
    if not path.is_file():
        if 'decaymass' in param:
            raise ValueError('DECAYMASS without a generated mass scheme; re-export')
        return 0., None
    scheme = json.loads(path.read_text())
    if (scheme.get('format') != 1 or scheme.get('production_model') != 'loop_sm-no_b_mass'
            or scheme.get('decay_model') != 'loop_sm'
            or scheme.get('parameter') != ['decaymass', 5]
            or scheme.get('alpha_s_flavours') != 5):
        raise ValueError('Unsupported exported decay mass scheme')
    mass = float(param['decaymass'].get((5,)).value)
    if not math.isfinite(mass) or mass <= 0.:
        raise ValueError('A massive-decay export requires positive DECAYMASS(5); re-export for zero')
    source = process / 'Source/MODEL/get_mass_width_fcts.f'
    if not source.is_file() or 'GET_DECAY_MASS_FROM_ID=DC_MDL_MB' not in source.read_text().upper():
        raise ValueError('Missing generated decay-only mass lookup; re-export')
    return mass, path


def export_w_treatment(process):
    """Infer the actual treatment from exported trees, not directory names.

    Inspect every production subprocess, including the Born core: two direct
    three-body tops alone do not distinguish ttW from a dileptonic ttbar run.
    """
    paths = sorted(process.glob('SubProcesses/P*/decay_chain_info.dat'))
    if not paths:
        raise ValueError('Missing exported decay topology; re-export this ttW process')
    treatments = set()
    for path in paths:
        records = [line.split() for line in path.read_text().splitlines() if line.strip()]
        nodes = {int(r[1]): (int(r[2]), int(r[3]), r[7:])
                 for r in records if r[0] == 'NODE'}
        leaves = {int(r[1]): int(r[3]) for r in records if r[0] == 'DECAY_LEAF'}
        born = {r[1] for r in records if r[0] == 'CONTEXT' and r[2] == 'BORN'}
        if len(born) != 1:
            raise ValueError('Expected one Born core in ' + str(path))
        final = sorted(int(r[3]) for r in records
                       if r[0] == 'CORE_LEG' and r[1] in born and r[4] == 'F')
        tops = {pdg: node for node, (parent, pdg, _) in nodes.items()
                if abs(pdg) == 6 and parent == 0}
        if set(tops) != {6, -6}:
            raise ValueError('Expected two on-shell top nodes in ' + str(path))
        ws = [(parent, pdg) for parent, pdg, _ in nodes.values() if abs(pdg) == 24]
        associated = [pdg for parent, pdg in ws if parent == 0]
        nested = [(parent, pdg) for parent, pdg in ws if parent != 0]
        if len(associated) == 1 and sorted(final) == sorted([6, -6] + associated):
            if set(nested) == {(tops[6], 24), (tops[-6], -24)} and len(nodes) == 5:
                treatment = 'onshell'
            elif not nested and len(nodes) == 3:
                treatment = 'top-bw'
            else:
                raise ValueError('Unsupported mixed W decay tree in ' + str(path))
        elif not ws and len(nodes) == 2 and len(final) == 4:
            leptons = sorted(pdg for pdg in final if abs(pdg) != 6)
            if leptons not in ([-11, 12], [-13, 14], [-12, 11], [-14, 13]):
                raise ValueError('Expected an associated e/mu charged current in ' + str(path))
            treatment = 'all-bw'
        else:
            raise ValueError('Not a supported trilepton ttW topology: ' + str(path))
        for parent, pdg, children in nodes.values():
            if len(children) % 2:
                raise ValueError('Malformed decay children in ' + str(path))
            daughters = sorted(nodes[int(ref)][1] if kind == 'NODE' else leaves[int(ref)]
                               for kind, ref in zip(children[::2], children[1::2]))
            sign = 1 if pdg > 0 else -1
            leptons = [[-11*sign, 12*sign], [-13*sign, 14*sign]]
            if abs(pdg) == 24:
                expected = [sorted(pair) for pair in leptons]
            elif abs(pdg) == 6 and treatment == 'onshell':
                expected = [sorted([5*sign, 24*sign])]
            elif abs(pdg) == 6:
                expected = [sorted([5*sign] + pair) for pair in leptons]
            else:
                expected = []
            if daughters not in expected:
                raise ValueError('Expected leptonic W/top decays in ' + str(path))
        treatments.add(treatment)
    if len(treatments) != 1:
        raise ValueError('W treatments differ between exported subprocesses')
    return treatments.pop(), paths


def mg5_root():
    for parent in Path(__file__).resolve().parents:
        if (parent / 'madgraph/fks/fks_decay.py').is_file():
            return parent
    raise ValueError('Run this script from the source checkout, not its exported copy')


def configure(args):
    root = mg5_root()
    sys.path.insert(0, str(root))
    from madgraph.fks import fks_decay
    from madgraph.various.banner import RunCardNLO
    from models.check_param_card import ParamCard

    process = Path(args.process_dir).resolve()
    cards = process / 'Cards'
    for name in ('param_card.dat', 'run_card.dat', 'run_card_default.dat', 'decay_card.dat'):
        if not (cards / name).is_file():
            raise ValueError('Missing fNLO bundle card: ' + str(cards / name))
    if not (process / 'SubProcesses/decay_chain_parameters.f90').is_file():
        raise ValueError('This is not a current fNLO export')
    if 'DECAY_SCALE_GROUPING' not in (process / 'SubProcesses/decay_chain_parameters.f90').read_text():
        raise ValueError('Re-export from the updated checkout to include signed-PDG scale axes')
    grouped_w = args.production_scale == 'core-w-ht-half'
    scale_sources = {name: process / 'SubProcesses' / name for name in
                     ('setscales.f90', 'decay_chain_scales.f90', 'decay_chain_parameters.f90')}
    if grouped_w:
        for name, marker in [('setscales.f90', 'w_system_core_ht_half'),
                             ('decay_chain_scales.f90', 'function w_system_core_ht_half'),
                             ('decay_chain_parameters.f90', 'PRODUCTION_SCALE_GROUPING')]:
            if not scale_sources[name].is_file() or marker not in scale_sources[name].read_text():
                raise ValueError('Re-export from the updated checkout to include W-system CORE HT/2')
    source = process / 'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90'
    if not source.is_file():
        raise ValueError('Re-export from the updated checkout to include the product analysis')
    treatment, topology_paths = export_w_treatment(process)
    if treatment != args.w_treatment:
        raise ValueError('Export has W treatment %s, not %s; generate a separate export' %
                         (treatment, args.w_treatment))
    production_sampling = getattr(args, 'production_sampling', 'flat')
    if production_sampling not in ('flat', 'w-current'):
        raise ValueError('Production sampling must be flat or w-current')
    sampling_sources = {name: process / 'SubProcesses' / name for name in (
        'phase_space_kinematics.f90', 'factorized_block_kinematics.f90',
        'decay_chain_parameters.f90', 'decay_chain_kinematics.f90',
        'nlo_decay_kinematics.f90')}
    if production_sampling == 'w-current':
        if treatment != 'all-bw':
            raise ValueError('W-current sampling requires an all-bw production core')
        for name, marker in (
                ('factorized_block_kinematics.f90', 'generate_factorized_current_nbody'),
                ('decay_chain_parameters.f90', 'PRODUCTION_PHASE_SPACE_SAMPLING'),
                ('decay_chain_kinematics.f90', 'call generate_current_nbody'),
                ('nlo_decay_kinematics.f90', 'call generate_current_nbody')):
            path = sampling_sources[name]
            if not path.is_file() or marker not in path.read_text():
                raise ValueError('Re-export to include W-current sampling in both core paths')
    internal_width_paths = [path.with_name('decay_internal_widths.json')
                            for path in topology_paths]
    if treatment != 'onshell':
        for path in internal_width_paths:
            if not path.is_file():
                raise ValueError('Re-export to preserve internal decay W widths at initialization')
            data = json.loads(path.read_text())
            if data.get('format') != 1 or 24 not in data.get('pdgs', []):
                raise ValueError('Export does not certify finite internal decay W widths')
    expected_width_treatment = 'onshell' if treatment == 'onshell' else 'bw'
    if args.top_width_w_treatment != expected_width_treatment:
        raise ValueError('%s requires top total widths calculated with %s Ws' %
                         (treatment, expected_width_treatment))
    if not args.width_source.strip():
        raise ValueError('Record the matched top-width calculation in --width-source')
    param = ParamCard(str(cards / 'param_card.dat'))
    mt = float(param['mass'].get((6,)).value)
    mw = float(param['mass'].get((24,)).value)
    ww = float(param['decay'].get((24,)).value)
    mb = float(param['mass'].get((5,)).value)
    if mb != 0.:
        raise ValueError('This setup uses loop_sm-no_b_mass; production bottom mass must be zero')
    decay_mb, mass_scheme_path = export_decay_bottom_mass(process, param)
    requested_mb = getattr(args, 'decay_bottom_mass', 0.)
    width_mb = getattr(args, 'top_width_bottom_mass', None)
    if not math.isfinite(requested_mb) or requested_mb < 0. or not math.isclose(
            decay_mb, requested_mb, rel_tol=1.e-12, abs_tol=1.e-12):
        raise ValueError('Requested decay bottom mass does not match the export/card')
    if width_mb is None and decay_mb == 0.:
        width_mb = 0.  # Backwards-compatible massless baseline.
    if width_mb is None or not math.isfinite(width_mb) or not math.isclose(
            decay_mb, width_mb, rel_tol=1.e-12, abs_tol=1.e-12):
        raise ValueError('Supply --top-width-bottom-mass matching the decay matrix elements')
    if decay_mb >= mt - (mw if treatment == 'onshell' else 0.):
        raise ValueError('Decay bottom mass closes the top-decay phase space')
    for value in (mt, mw, ww, args.top_width_lo, args.top_width_nlo, args.ecm):
        if not math.isfinite(value) or value <= 0.:
            raise ValueError('Masses, widths and collider energy must be positive and finite')
    if min(args.seed, args.points, args.grid_points, args.iterations) < 1:
        raise ValueError('Seed, integration counts and iterations must be positive')
    accuracy = getattr(args, 'accuracy', -1.)
    if not math.isfinite(accuracy) or not (accuracy == -1. or 0. < accuracy < 1.):
        raise ValueError('Integration accuracy must be -1 (fixed counts) or between zero and one')
    job_time = getattr(args, 'job_seconds', 0.)
    if not math.isfinite(job_time) or job_time < 0.:
        raise ValueError('Job target time must be finite and nonnegative')
    run_name = '%s_%s_%s_%s_%d' % (args.variant, treatment, args.production_scale,
                                  args.decay_scales, args.seed)
    if production_sampling == 'w-current':
        run_name += '_wcurrent'
    if decay_mb:
        run_name += '_mb%s' % ('%.12g' % decay_mb).replace('.', 'p')
    archive = process / 'study_cards' / run_name
    if archive.exists() or (process / 'Events' / run_name).exists():
        raise ValueError('Run or card archive already exists: ' + run_name)
    production, decay, combination = VARIANTS[args.variant]
    fixed = args.production_scale == 'fixed'
    settings = dict(
        req_acc_fo=accuracy, npoints_fo_grid=args.grid_points, niters_fo_grid=1,
        npoints_fo=args.points, niters_fo=args.iterations, iseed=args.seed,
        fo_job_target_time=job_time,
        lpp1=1, lpp2=1, ebeam1=args.ecm / 2., ebeam2=args.ecm / 2.,
        fixed_ren_scale=fixed, fixed_fac_scale=fixed, fixed_qes_scale=fixed,
        mur_ref_fixed=mt+mw/2., muf_ref_fixed=mt+mw/2., qes_ref_fixed=mt+mw/2.,
        dynamical_scale_choice=3, mur_over_ref=1., muf_over_ref=1.,
        reweight_scale=True, rw_rscale=[1., .5, 2.], rw_fscale=[1., .5, 2.],
        reweight_pdf=args.pdf_id is not None,
        pdlabel='lhapdf' if args.pdf_id is not None else 'nn23nlo',
        maxjetflavor=5,
        cut_decays=False, ptj=0., etaj=-1., ptl=0., etal=-1.,
        drll=0., drll_sf=0., mll=0., mll_sf=0., ptgmin=0.,
        pt_min_pdg={}, pt_max_pdg={}, mxx_min_pdg={}, custom_fcts=[],
    )
    if args.pdf_id is not None:
        settings['lhaid'] = args.pdf_id
    run = RunCardNLO(str(cards / 'run_card.dat'))
    for key, value in settings.items():
        # Persist explicit settings even when absent from the default card
        # (e.g. the hidden maxjetflavor input in an fNLO export).
        run.set(key, value, user=True, raiseerror=True)
    widths, references = {6: args.top_width_lo}, {6: mt}
    if treatment != 'all-bw':
        widths[24], references[24] = ww, mw
    decay_text = fks_decay.decay_card_text(
        widths, references,
        nlo_width_pdgs={6}, nlo_widths={6: args.top_width_nlo},
        nlo_decay_combination=combination, production_order=production, decay_order=decay,
        decay_scale_variation_mode='INDEPENDENT', decay_scale_factors=(1., .5, 2.),
        decay_scale_grouping='SIGNED_PDG' if args.decay_scales == 'separate' else 'SPECIES',
        decay_width_scale_modes={6: 'AUTO'}, production_scale_momenta='CORE',
        production_scale_grouping='W_SYSTEM' if grouped_w else 'NONE',
        production_phase_space_sampling='W_CURRENT' if production_sampling == 'w-current' else 'FLAT',
        production_sampling_mass=mw if production_sampling == 'w-current' else None,
        production_sampling_width=ww if production_sampling == 'w-current' else None)
    # Preserve the input cards before replacing any of them.
    archive.mkdir(parents=True)
    (archive / 'before').mkdir()
    for name in ('run_card.dat', 'decay_card.dat', 'FO_analyse_card.dat'):
        if (cards / name).is_file():
            shutil.copy2(cards / name, archive / 'before' / name)
    run.write(str(cards / 'run_card.dat'), template=str(cards / 'run_card_default.dat'))
    (cards / 'decay_card.dat').write_text(decay_text)
    (cards / 'FO_analyse_card.dat').write_text(
        'FO_ANALYSIS_FORMAT=HwU\nFO_ANALYSE=analysis_HwU_pp_ttxw_product.o\n'
        'FO_EXTRALIBS=\nFO_EXTRAPATHS=\nFO_INCLUDEPATHS=\n')
    for name in ('run_card.dat', 'decay_card.dat', 'FO_analyse_card.dat', 'param_card.dat'):
        shutil.copy2(cards / name, archive / name)
    if mass_scheme_path is not None:
        shutil.copy2(mass_scheme_path, archive / mass_scheme_path.name)
    hashes = {name: hashlib.sha256((archive / name).read_bytes()).hexdigest()
              for name in ('run_card.dat', 'decay_card.dat', 'FO_analyse_card.dat', 'param_card.dat')}
    hashes['analysis_source'] = hashlib.sha256(source.read_bytes()).hexdigest()
    topology_hashes = {str(path.relative_to(process)): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in topology_paths}
    manifest = dict(variant=args.variant, run_name=run_name, settings=settings,
                    w_treatment=treatment, tops_on_shell=True,
                    production_bottom_mass=mb, decay_bottom_mass=decay_mb,
                    bottom_mass_scheme='on-shell' if decay_mb else 'massless',
                    alpha_s_flavours=5, top_width_bottom_mass=width_mb,
                    decay_mass_scheme_hash=(hashlib.sha256(mass_scheme_path.read_bytes()).hexdigest()
                                            if mass_scheme_path else None),
                    top_width_lo=args.top_width_lo, top_width_nlo=args.top_width_nlo,
                    top_width_w_treatment=args.top_width_w_treatment,
                    width_source=args.width_source, top_width_reference_scale=mt,
                    w_width=ww, production_scale=args.production_scale,
                    production_sampling=production_sampling,
                    production_sampling_parameters=(dict(mass_GeV=mw, width_GeV=ww)
                                                     if production_sampling == 'w-current' else None),
                    sampling_runtime_hashes={name: hashlib.sha256(path.read_bytes()).hexdigest()
                                             for name, path in sampling_sources.items() if path.is_file()},
                    production_scale_grouping='W_SYSTEM' if grouped_w else 'NONE',
                    production_scale_objects=('t, tbar, associated W system, radiation' if grouped_w else
                                              'native production-core particles' if not fixed else
                                              'fixed mt + MW/2'),
                    production_w_virtuality='q = p_lepton + p_neutrino; no pole-mass projection'
                                            if grouped_w else None,
                    scale_runtime_hashes={name: hashlib.sha256(path.read_bytes()).hexdigest()
                                          for name, path in scale_sources.items() if path.is_file()},
                    production_core_objects=('t, tbar, associated lepton, associated neutrino, radiation'
                                             if treatment == 'all-bw' else
                                             't, tbar, associated on-shell W, radiation'),
                    w_width_prescription=('NWA normalization only' if treatment == 'onshell' else
                                          'fixed-width internal propagators with real masses'),
                    topology_hashes=topology_hashes,
                    internal_width_hashes={str(path.relative_to(process)):
                                           hashlib.sha256(path.read_bytes()).hexdigest()
                                           for path in internal_width_paths if path.is_file()},
                    decay_scale_grouping=args.decay_scales, hashes=hashes,
                    git_head=subprocess.check_output(
                        ['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip(),
                    git_status=subprocess.check_output(
                        ['git', '-C', str(root), 'status', '--short'], text=True),
                    status='configured; no integration performed')
    (archive / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    if treatment == 'all-bw' and args.production_scale == 'core-ht-half':
        print('Note: native dynamic CORE HT/2 resolves the associated lepton and neutrino. '
              'Comparison with an on-shell associated W also changes the scale definition; '
              'see the production-scale study in ttw_product_study.md.', file=sys.stderr)
    print('Cards archived in ' + str(archive))
    print('Next, from %s: bin/calculate_xsect NLO -f -n %s' % (process, run_name))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    cmd = sub.add_parser('commands', help='print commands; does not generate or launch')
    cmd.add_argument('--charge', choices=['plus', 'minus'], required=True)
    cmd.add_argument('--flavours', choices=['e', 'mu'], nargs=3, default=['e', 'e', 'mu'])
    cmd.add_argument('--output', required=True)
    cmd.add_argument('--corrected', choices=['both', 't', 'tbar', 'neither'], default='both')
    cmd.add_argument('--w-treatment', choices=W_TREATMENTS, default='onshell',
                     help='on-shell Ws, BW Ws in top decays only, or BW Ws everywhere')
    cmd.add_argument('--decay-bottom-mass', type=float, default=0.,
                     help='generation-time on-shell bottom mass in decays only, in GeV')
    cmd.add_argument('--real-only', action='store_true', help='generation smoke test only, not physical NLO')
    cmd.add_argument('--partonic', action='store_true', help='one initial channel for a generation smoke test')
    config = sub.add_parser('configure', help='update exported cards, with a snapshot; does not run')
    config.add_argument('--process-dir', required=True)
    config.add_argument('--variant', choices=list(VARIANTS), required=True)
    config.add_argument('--top-width-lo', type=float, required=True)
    config.add_argument('--top-width-nlo', type=float, required=True)
    config.add_argument('--w-treatment', choices=W_TREATMENTS, default='onshell',
                        help='must match the generated topology; cannot change it through cards')
    config.add_argument('--top-width-w-treatment', choices=['onshell', 'bw'], required=True,
                        help='explicit W convention of BOTH supplied physical total top widths')
    config.add_argument('--decay-bottom-mass', type=float, default=0.,
                        help='must match the generated DECAYMASS(5), or zero for massless decays')
    config.add_argument('--top-width-bottom-mass', type=float,
                        help='bottom mass used in BOTH total widths; required for massive decays')
    config.add_argument('--width-source', required=True, help='matched width calculation and input record')
    config.add_argument('--pdf-id', type=int, help='LHAPDF ID; absent selects bundled nn23nlo for pilots')
    config.add_argument('--production-scale', choices=PRODUCTION_SCALES, default='core-w-ht-half',
                        help='default groups the associated W current; core-ht-half is the native diagnostic')
    config.add_argument('--production-sampling', choices=['flat', 'w-current'], default='flat',
                        help='opt-in all-bw current proposal; does not change propagators or physical widths')
    config.add_argument('--decay-scales', choices=['separate', 'shared'], default='separate',
                        help='separate t/tbar axes (81 points), or shared (27); with NLO decays')
    config.add_argument('--ecm', type=float, default=13000.)
    config.add_argument('--seed', type=int, default=31701)
    config.add_argument('--points', type=int, default=1000)
    config.add_argument('--grid-points', type=int, default=200)
    config.add_argument('--iterations', type=int, default=3)
    config.add_argument('--accuracy', type=float, default=-1.,
                        help='positive total-rate target enables adaptive integration; -1 uses fixed counts')
    config.add_argument('--job-seconds', type=float, default=0.,
                        help='optional target CPU seconds per refinement split; 0 keeps default scheduling')
    args = parser.parse_args()
    try:
        if args.action == 'commands':
            print(process_commands(args.charge, args.flavours, args.output,
                                   args.corrected, args.real_only, args.partonic,
                                   args.w_treatment, args.decay_bottom_mass), end='')
        else:
            configure(args)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(2, str(error) + '\n')


if __name__ == '__main__':
    main()
