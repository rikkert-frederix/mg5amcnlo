#!/usr/bin/env python3
"""Generate ttW commands or configure a fresh fNLO export; never launch a run."""
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


def process_commands(charge, flavours, output, corrected='both', real_only=False,
                     partonic=False):
    """One export per charge and ordered (Wt,Wtbar,Wassoc) flavour assignment."""
    output = str(Path(output).expanduser().absolute())
    if not re.fullmatch(r'[\w./+\-]+', output):
        raise ValueError('MadGraph output path must contain no spaces or shell metacharacters')
    lp = {'e': 'e+ ve', 'mu': 'mu+ vm'}
    lm = {'e': 'e- ve~', 'mu': 'mu- vm~'}
    order = '[real=QCD]' if real_only else '[QCD]'
    torder = ' ' + order if corrected in ('both', 't') else ''
    aorder = ' ' + order if corrected in ('both', 'tbar') else ''
    w = 'w+' if charge == 'plus' else 'w-'
    initial = ('u d~' if charge == 'plus' else 'd u~') if partonic else 'p p'
    assoc = (lp if charge == 'plus' else lm)[flavours[2]]
    return '\n'.join([
        'import model loop_sm-no_b_mass',
        'define p = g u c d s b u~ c~ d~ s~ b~',
        'generate %s > t t~ %s QCD=2 QED=1 %s, '
        '(t > w+ b QED=1%s, w+ > %s), '
        '(t~ > w- b~ QED=1%s, w- > %s), %s > %s' %
        (initial, w, order, torder, lp[flavours[0]], aorder, lm[flavours[1]], w, assoc),
        'output fNLO ' + output,
        '',
    ])


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
    source = process / 'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90'
    if not source.is_file():
        raise ValueError('Re-export from the updated checkout to include the product analysis')
    param = ParamCard(str(cards / 'param_card.dat'))
    mt = float(param['mass'].get((6,)).value)
    mw = float(param['mass'].get((24,)).value)
    ww = float(param['decay'].get((24,)).value)
    mb = float(param['mass'].get((5,)).value)
    if mb != 0.:
        raise ValueError('This setup uses loop_sm-no_b_mass; bottom mass must be zero')
    for value in (mt, mw, ww, args.top_width_lo, args.top_width_nlo, args.ecm):
        if not math.isfinite(value) or value <= 0.:
            raise ValueError('Masses, widths and collider energy must be positive and finite')
    if min(args.seed, args.points, args.grid_points, args.iterations) < 1:
        raise ValueError('Seed, integration counts and iterations must be positive')
    run_name = '%s_%s_%s_%d' % (args.variant, args.production_scale, args.decay_scales, args.seed)
    archive = process / 'study_cards' / run_name
    if archive.exists() or (process / 'Events' / run_name).exists():
        raise ValueError('Run or card archive already exists: ' + run_name)
    production, decay, combination = VARIANTS[args.variant]
    fixed = args.production_scale == 'fixed'
    settings = dict(
        req_acc_fo=-1., npoints_fo_grid=args.grid_points, niters_fo_grid=1,
        npoints_fo=args.points, niters_fo=args.iterations, iseed=args.seed,
        lpp1=1, lpp2=1, ebeam1=args.ecm / 2., ebeam2=args.ecm / 2.,
        fixed_ren_scale=fixed, fixed_fac_scale=fixed, fixed_qes_scale=fixed,
        mur_ref_fixed=mt+mw/2., muf_ref_fixed=mt+mw/2., qes_ref_fixed=mt+mw/2.,
        dynamical_scale_choice=3, mur_over_ref=1., muf_over_ref=1.,
        reweight_scale=True, rw_rscale=[1., .5, 2.], rw_fscale=[1., .5, 2.],
        reweight_pdf=args.pdf_id is not None,
        pdlabel='lhapdf' if args.pdf_id is not None else 'nn23nlo',
        cut_decays=False, ptj=0., etaj=-1., ptl=0., etal=-1.,
        drll=0., drll_sf=0., mll=0., mll_sf=0., ptgmin=0.,
        pt_min_pdg={}, pt_max_pdg={}, mxx_min_pdg={}, custom_fcts=[],
    )
    if args.pdf_id is not None:
        settings['lhaid'] = args.pdf_id
    run = RunCardNLO(str(cards / 'run_card.dat'))
    for key, value in settings.items():
        run[key] = value
    decay_text = fks_decay.decay_card_text(
        {6: args.top_width_lo, 24: ww}, {6: mt, 24: mw},
        nlo_width_pdgs={6}, nlo_widths={6: args.top_width_nlo},
        nlo_decay_combination=combination, production_order=production, decay_order=decay,
        decay_scale_variation_mode='INDEPENDENT', decay_scale_factors=(1., .5, 2.),
        decay_scale_grouping='SIGNED_PDG' if args.decay_scales == 'separate' else 'SPECIES',
        decay_width_scale_modes={6: 'AUTO'}, production_scale_momenta='CORE')
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
    hashes = {name: hashlib.sha256((archive / name).read_bytes()).hexdigest()
              for name in ('run_card.dat', 'decay_card.dat', 'FO_analyse_card.dat', 'param_card.dat')}
    hashes['analysis_source'] = hashlib.sha256(source.read_bytes()).hexdigest()
    manifest = dict(variant=args.variant, run_name=run_name, settings=settings,
                    top_width_lo=args.top_width_lo, top_width_nlo=args.top_width_nlo,
                    width_source=args.width_source, top_width_reference_scale=mt,
                    w_width=ww, production_scale=args.production_scale,
                    decay_scale_grouping=args.decay_scales, hashes=hashes,
                    git_head=subprocess.check_output(
                        ['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip(),
                    git_status=subprocess.check_output(
                        ['git', '-C', str(root), 'status', '--short'], text=True),
                    status='configured; no integration performed')
    (archive / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
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
    cmd.add_argument('--real-only', action='store_true', help='generation smoke test only, not physical NLO')
    cmd.add_argument('--partonic', action='store_true', help='one initial channel for a generation smoke test')
    config = sub.add_parser('configure', help='update exported cards, with a snapshot; does not run')
    config.add_argument('--process-dir', required=True)
    config.add_argument('--variant', choices=list(VARIANTS), required=True)
    config.add_argument('--top-width-lo', type=float, required=True)
    config.add_argument('--top-width-nlo', type=float, required=True)
    config.add_argument('--width-source', required=True, help='matched width calculation and input record')
    config.add_argument('--pdf-id', type=int, help='LHAPDF ID; absent selects bundled nn23nlo for pilots')
    config.add_argument('--production-scale', choices=['core-ht-half', 'fixed'], default='core-ht-half')
    config.add_argument('--decay-scales', choices=['separate', 'shared'], default='separate',
                        help='separate t/tbar axes (81 points), or shared (27); with NLO decays')
    config.add_argument('--ecm', type=float, default=13000.)
    config.add_argument('--seed', type=int, default=31701)
    config.add_argument('--points', type=int, default=1000)
    config.add_argument('--grid-points', type=int, default=200)
    config.add_argument('--iterations', type=int, default=3)
    args = parser.parse_args()
    try:
        if args.action == 'commands':
            print(process_commands(args.charge, args.flavours, args.output,
                                   args.corrected, args.real_only, args.partonic), end='')
        else:
            configure(args)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(2, str(error) + '\n')


if __name__ == '__main__':
    main()
