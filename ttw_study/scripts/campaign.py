#!/usr/bin/env python3
"""Reproducible preflight and serialized generation for the full ttW study."""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / 'ttw_study'
WIDTH = ROOT.parent / 'top-decay-virtual-cdr/width'
PDF_NAME = 'NNPDF40_nlo_as_01180'
SPEC = ROOT / 'Template/fNLO/FixedOrderAnalysis/ttw_product_study.md'
GF, MT, MW, MZ = 1.1663787e-5, 172.5, 80.385, 91.1876
CORES = 64
RUNTIME_SOURCES = (
    'madgraph/fks/fks_decay.py', 'models/check_param_card.py',
    'madgraph/iolibs/export_fks.py', 'madgraph/iolibs/export_spin_density.py',
    'madgraph/iolibs/template_files/loop_optimized/CT_interface.inc',
    'madgraph/interface/amcatnlo_run_interface.py', 'madgraph/various/banner.py',
    'madgraph/madevent/sum_html.py',
    'Template/fNLO/Source/ranmar.f90',
    'Template/fNLO/SubProcesses/phase_space_kinematics.f90',
    'Template/fNLO/SubProcesses/factorized_block_kinematics.f90',
    'Template/fNLO/SubProcesses/decay_chain_kinematics.f90',
    'Template/fNLO/SubProcesses/nlo_decay_kinematics.f90',
    'Template/fNLO/SubProcesses/multiplicative_nlo_decay.f90',
    'Template/fNLO/SubProcesses/driver_mintFO.f90',
    'Template/fNLO/SubProcesses/makefile_fks_dir',
    'Template/fNLO/SubProcesses/decay_chain_parameters.f90',
    'Template/fNLO/SubProcesses/decay_chain_scales.f90',
    'Template/fNLO/SubProcesses/setscales.f90',
    'Template/fNLO/SubProcesses/fks_weights.f90',
    'Template/fNLO/SubProcesses/spin_density_weight_lines.f90',
    'Template/fNLO/SubProcesses/top_decay_virtual_dispatch.f90',
    'Template/fNLO/SubProcesses/ajob_template',
    'Template/fNLO/FixedOrderAnalysis/HwU.f90',
    'Template/fNLO/FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90',
    'Template/fNLO/SubProcesses/mint_module.f90',
    'Template/fNLO/FixedOrderAnalysis/ttw_product_setup.py',
    'ttw_study/scripts/campaign.py')


def now():
    return datetime.now(timezone.utc).isoformat()


def output(command, cwd=ROOT):
    return subprocess.check_output(list(map(str, command)), cwd=cwd, text=True).strip()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def lhapdf_environment(base=None):
    """Put the frozen local PDF data before any inherited LHAPDF path."""
    data = ROOT/'LHAPDF/share/LHAPDF'
    if not (data/'lhapdf.conf').is_file() or not (data/'pdfsets.index').is_file():
        raise ValueError('Local LHAPDF data directory is incomplete: '+str(data))
    environment = dict(os.environ if base is None else base)
    inherited = [path for path in environment.get('LHAPDF_DATA_PATH','').split(':') if path]
    environment['LHAPDF_DATA_PATH'] = ':'.join(dict.fromkeys([str(data),*inherited]))
    return environment


def calculator(command):
    raw = output(command, WIDTH)
    values = {}
    for line in raw.splitlines():
        if '=' not in line:
            continue
        key, value = line.split('=', 1)
        try:
            values[key.strip()] = float(value.strip())
        except ValueError:
            values[key.strip()] = value.strip()
    return {'command': list(map(str, command)), 'stdout': raw, 'values': values}


def preflight(args):
    import lhapdf
    lhapdf.setVerbosity(0)
    pdf = lhapdf.mkPDF(PDF_NAME, 0)
    pdf_dir = ROOT / 'LHAPDF/share/LHAPDF' / PDF_NAME
    metadata = pdf_dir / (PDF_NAME + '.info')
    expected_hash = '9362949bda8c0ae1d6ba37542ee31ddbae010369d126ac45ab68ce8b306212a1'
    if digest(metadata) != expected_hash:
        raise ValueError('Official PDF metadata differs from the frozen benchmark')
    checkpoints = [(MW, .1202908028622), (MT/2, .1190016017860),
                   (MZ, .1180021695134), (MT, .1076703477820),
                   (MT*2, .0983392832531)]
    coupling = []
    for scale, expected in checkpoints:
        actual = pdf.alphasQ(scale)
        if not math.isclose(actual, expected, rel_tol=0., abs_tol=5e-13):
            raise ValueError('PDF alpha-s checkpoint failed at %s' % scale)
        coupling.append({'Q_GeV': scale, 'alpha_s': actual, 'benchmark': expected})
    members = []
    for member in range(101):
        path = pdf_dir / ('%s_%04d.dat' % (PDF_NAME, member))
        replica = lhapdf.mkPDF(PDF_NAME, member)
        sample = replica.xfxQ(21, .1, MT)
        if not math.isfinite(sample):
            raise ValueError('Nonfinite PDF member %s' % member)
        members.append({'member': member, 'sha256': digest(path),
                        'xg_at_x0p1_Qmt': sample})
    w = calculator([WIDTH / 'w_decay_width', MW, MW, 5, .118])
    v = w['values']
    w0 = v['Gamma_W_LO [GeV]']
    wc = v['Delta_Gamma_QCD_NLO [GeV]'] / v['alpha_s(muR)']
    ww = w0 + wc * pdf.alphasQ(MW)
    leptonic = GF * MW**3/(6*math.pi*math.sqrt(2))
    assert abs(ww-leptonic*(3+6*(1+pdf.alphasQ(MW)/math.pi))) < 1e-13
    w.update(gamma_lo=w0, qcd_coefficient=wc, gamma_nlo_pdf=ww,
             gamma_lepton=leptonic, branching_e=leptonic/ww)
    reference = {(0., 'onshell'): (1.4806285092, 1.3535485212),
                 (0., 'bw'): (1.4576010900, 1.3324848297),
                 (4.8, 'onshell'): (1.4765338515, 1.3513387258),
                 (4.8, 'bw'): (1.4535389808, 1.3302875304)}
    tops = []
    for mb in (0., 4.8):
        for mode in ('onshell', 'bw'):
            for factor in (.5, 1., 2.):
                q = MT*factor
                row = calculator([WIDTH/'top_decay_width', MT, MW, mb, q,
                                  0. if mode == 'onshell' else ww, .118, 1., 5])
                v = row['values']
                g0 = v['Gamma_LO [GeV]']
                coefficient = v['Delta_Gamma_QCD^NLO [GeV]']/v['alpha_s(muR)']
                gnlo = g0+coefficient*pdf.alphasQ(q)
                if factor == 1.:
                    for actual, expected in zip((g0, gnlo), reference[mb, mode]):
                        assert abs(actual-expected) < 5.1e-11, (actual, expected)
                row.update(mb_GeV=mb, w_treatment=mode, muR_GeV=q,
                           scale_factor=factor, gamma_lo=g0, qcd_coefficient=coefficient,
                           alpha_s_pdf=pdf.alphasQ(q), gamma_nlo_pdf=gnlo)
                tops.append(row)
    alpha = math.sqrt(2)*GF*MW**2*(1-MW**2/MZ**2)/math.pi
    data = {'created_utc': now(), 'specification_sha256': digest(SPEC),
            'git_head': output(['git', 'rev-parse', 'HEAD']),
            'width_git_head': output(['git', 'rev-parse', 'HEAD'], WIDTH),
            'width_sources': {p.name: digest(p) for p in WIDTH.glob('*.f90')},
            'compiler': output(['gfortran', '--version']).splitlines()[0],
            'cxx_compiler': output(['g++', '--version']).splitlines()[0],
            'python': sys.version, 'lhapdf': lhapdf.version(),
            'fastjet': output(['fastjet-config', '--version']),
            'PDF': {'name': PDF_NAME, 'id': 331700, 'DataVersion': 1,
                    'metadata_sha256': digest(metadata), 'members': members},
            'inputs': {'mt': MT, 'mw': MW, 'mz': MZ, 'GF': GF,
                       'alpha_Gmu_inverse': 1/alpha, 'sin2_theta_W': 1-MW**2/MZ**2,
                       'Gamma_W_fixed': ww, 'production_mb': 0., 'nf': 5},
            'coupling_checkpoints': coupling, 'W_width': w, 'top_widths': tops,
            'max_MG5_cores': CORES, 'status': 'preflight passed; no integration'}
    save(STUDY/'inputs/benchmark.json', data)
    (STUDY/'inputs/preflight_source_changes.patch').write_text(output(['git', 'diff']))
    print('Preflight passed: 101 PDFs, five couplings, fixed W and 12 top-width scale rows.')
    print('Gamma_W = %.13f GeV; 1/alpha_Gmu = %.13f' % (ww, 1/alpha))
    for row in tops:
        if row['scale_factor'] == 1.:
            print('mb=%s %s: Gamma_LO=%.13f, Gamma_NLO_PDF=%.13f GeV' % (
                row['mb_GeV'], row['w_treatment'], row['gamma_lo'], row['gamma_nlo_pdf']))


def setup_module():
    spec = importlib.util.spec_from_file_location(
        'ttw_setup', ROOT/'Template/fNLO/FixedOrderAnalysis/ttw_product_setup.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def generate(args):
    if not (STUDY/'inputs/benchmark.json').exists():
        raise ValueError('Run the physical-input preflight first')
    name = 'TTW%s_%s_%s_mb%s_%s' % (
        args.charge, args.w_treatment, ''.join(args.flavours),
        str(args.decay_bottom_mass).replace('.', 'p'), args.corrected)
    if args.tag:
        if not args.tag.replace('_', '').isalnum():
            raise ValueError('Export tags may contain letters, digits and underscores')
        name += '_'+args.tag
    process = STUDY/'processes'/name
    if process.exists():
        raise ValueError('Output already exists; inspect its state before further work: '+str(process))
    card = STUDY/'inputs'/(name+'.mg5')
    preamble = '\n'.join([
        'set auto_update 0', 'set automatic_html_opening False',
        'set notification_center False', 'set nb_core %d' % CORES,
        'set run_mode 2', 'set low_mem_multicore_nlo_generation False',
        'set lhapdf '+str(ROOT/'LHAPDF/bin/lhapdf-config'),
        'set fastjet '+str(STUDY/'local/fastjet/bin/fastjet-config'), ''])
    card.write_text(preamble+setup_module().process_commands(
        args.charge, args.flavours, str(process), args.corrected,
        w_treatment=args.w_treatment, decay_bottom_mass=args.decay_bottom_mass))
    log = STUDY/'logs'/(name+'_generation.log')
    metadata = {'name': name, 'started_utc': now(), 'parent_pid': os.getpid(),
                'command_card': str(card), 'log': str(log), 'max_cores': CORES,
                'status': 'running', 'git_head': output(['git', 'rev-parse', 'HEAD']),
                'args': {k:v for k,v in vars(args).items() if k != 'func'}}
    source_patch = STUDY/'inputs'/(name+'_source.patch')
    source_patch.write_text(output(['git', 'diff']))
    metadata['source_patch_sha256'] = digest(source_patch)
    record = STUDY/'inputs'/(name+'_generation.json')
    # A single campaign lock serializes all MG5 launches. Actual process/exit
    # status, not this record, is authoritative when resuming an interrupted turn.
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        start = time.monotonic()
        environment = lhapdf_environment()
        metadata['environment'] = {'LHAPDF_DATA_PATH':environment['LHAPDF_DATA_PATH']}
        with log.open('w') as stream:
            command = ['taskset', '-c', '0-63', sys.executable, str(ROOT/'bin/mg5_aMC'), str(card)]
            child = subprocess.Popen(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                                     env=environment)
            metadata['child_pid'] = child.pid
            save(record, metadata)
            result = child.wait()
        metadata.update(returncode=result, elapsed_s=time.monotonic()-start,
                        finished_utc=now(), status='finished' if result == 0 else 'failed')
        if not (process/'Cards/decay_card.dat').exists():
            metadata['status'] = 'failed: missing completed export'
        save(record, metadata)
    print(json.dumps(metadata, indent=2))
    if metadata['status'] != 'finished':
        raise RuntimeError('Generation failed; see '+str(log))


def benchmark_param_card(process, benchmark, w_treatment):
    from models.check_param_card import ParamCard
    from models import import_ufo, model_reader
    path = process/'Cards/param_card.dat'
    original = process/'study_cards/original_param_card.dat'
    original.parent.mkdir(exist_ok=True)
    if not original.exists():
        shutil.copy2(path, original)
    card = ParamCard(str(path))
    mb = float(card['decaymass'].get((5,)).value) if 'decaymass' in card else 0.
    width_mode = 'onshell' if w_treatment == 'onshell' else 'bw'
    width = next(row for row in benchmark['top_widths']
                 if row['mb_GeV'] == mb and row['w_treatment'] == width_mode
                 and row['scale_factor'] == 1.)
    params = {('sminputs', 1): benchmark['inputs']['alpha_Gmu_inverse'],
              ('sminputs', 2): GF,
              ('sminputs', 3): next(p['alpha_s'] for p in benchmark['coupling_checkpoints']
                                    if p['Q_GeV'] == MZ),
              ('mass', 6): MT, ('yukawa', 6): MT, ('mass', 23): MZ,
              ('mass', 5): 0., ('yukawa', 5): 0.,
              ('mass', 25): 125., ('mass', 15): 1.777,
              ('decay', 6): width['gamma_nlo_pdf'],
              ('decay', 24): benchmark['W_width']['gamma_nlo_pdf'],
              ('decay', 23): 2.4952, ('decay', 25): .00407, ('decay', 15): 0.}
    for (block, code), value in params.items():
        try:
            card[block].get((code,)).value = value
        except KeyError:
            # Restricted zero Yukawas are absent from fresh cards; retain the
            # restriction instead of inventing an external model parameter.
            if value != 0.:
                raise
    # Recompute internal MW from the actual restricted model's weak inputs.
    # The decay-only mass is independent of these production parameters.
    model = model_reader.ModelReader(import_ufo.import_model('loop_sm-no_b_mass'))
    card.update_dependent(model, None, 20)
    model.set_parameters_and_couplings(card)
    assert abs(float(card['mass'].get((24,)).value)-MW) < 1e-10
    assert abs(model.get_mass(24)-MW) < 1e-10
    card.write(str(path), precision=16)
    save(process/'study_cards/benchmark_parameters.json', {
        'created_utc': now(), 'param_card_sha256': digest(path),
        'model_internal_MW': model.get_mass(24),
        'benchmark_sha256': digest(STUDY/'inputs/benchmark.json'),
        'decay_mb': mb, 'top_width': width})
    return mb, width


def run(args):
    process = Path(args.process_dir).resolve()
    if process.parent != STUDY/'processes':
        raise ValueError('Run only a process generated within this campaign')
    benchmark = json.loads((STUDY/'inputs/benchmark.json').read_text())
    setup = setup_module()
    treatment, _ = setup.export_w_treatment(process)
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        mb, width = benchmark_param_card(process, benchmark, treatment)
        settings = argparse.Namespace(
            process_dir=str(process), variant=args.variant,
            top_width_lo=width['gamma_lo'], top_width_nlo=width['gamma_nlo_pdf'],
            w_treatment=treatment, top_width_w_treatment=width['w_treatment'],
            decay_bottom_mass=mb, top_width_bottom_mass=mb,
            width_source=str(STUDY/'inputs/benchmark.json')+'; sha256='+
                         digest(STUDY/'inputs/benchmark.json'),
            pdf_id=331700, production_scale=args.production_scale,
            production_sampling=getattr(args, 'production_sampling', 'flat'),
            decay_scales='separate', ecm=args.ecm, seed=args.seed,
            points=args.points, grid_points=args.grid_points, iterations=args.iterations,
            accuracy=args.accuracy, job_seconds=getattr(args, 'job_seconds', 60.))
        setup.configure(settings)
        run_name = '%s_%s_%s_separate_%d' % (
            args.variant, treatment, args.production_scale, args.seed)
        if settings.production_sampling == 'w-current':
            run_name += '_wcurrent'
        if mb:
            run_name += '_mb%s' % ('%.12g' % mb).replace('.', 'p')
        archive = process/'study_cards'/run_name
        grid_reference = getattr(args,'grid_reference',None)
        grid_provenance = None
        if grid_reference:
            grid_provenance = import_grid(process,archive,Path(grid_reference))
        log = STUDY/'logs'/(process.name+'_'+run_name+'.log')
        command = ['taskset', '-c', '0-63', sys.executable,
                   str(process/'bin/calculate_xsect'), 'NLO', '-f', '-n', run_name,
                   '--multicore', '--nb_core=64']
        if grid_reference:
            command.append('--only_generation')
        record = {'started_utc': now(), 'parent_pid': os.getpid(), 'command': command,
                  'log': str(log), 'status': 'running', 'max_cores': CORES,
                  'source_hashes': {name: digest(ROOT/name) for name in RUNTIME_SOURCES},
                  'export_generation': json.loads((STUDY/'inputs'/
                                                  (process.name+'_generation.json')).read_text()),
                  'pilot': not args.main_run,
                  'grid_reference': grid_provenance,
                  'requested_total_rate_accuracy': args.accuracy}
        start = time.monotonic()
        environment = lhapdf_environment()
        record['environment'] = {'LHAPDF_DATA_PATH':environment['LHAPDF_DATA_PATH']}
        with log.open('w') as stream:
            child = subprocess.Popen(command, cwd=process, stdout=stream, stderr=subprocess.STDOUT,
                                     env=environment)
            record['child_pid'] = child.pid
            save(archive/'execution.json', record)
            code = child.wait()
        hwu = sorted((process/'Events'/run_name).glob('*.HwU'))
        record.update(finished_utc=now(), elapsed_s=time.monotonic()-start,
                      returncode=code, outputs={str(p): digest(p) for p in hwu},
                      status='finished' if code == 0 and hwu else 'failed')
        save(archive/'execution.json', record)
    print(json.dumps(record, indent=2))
    if record['status'] != 'finished':
        raise RuntimeError('Integration failed or lacks HwU output; see '+str(log))


def import_grid(process, archive, reference_manifest):
    """Transfer trained grids after a scheduling-only fix, never final weights.

    The native --only_generation path rewrites saved job directory roots and
    performs an independent refinement. Require identical physical cards and
    topology; the source failed run and its trained grids remain unchanged.
    """
    reference_manifest = reference_manifest.resolve()
    old_archive = reference_manifest.parent
    old_process = old_archive.parents[1]
    if old_process.parent != STUDY/'processes' or old_process == process:
        raise ValueError('Grid recovery requires another study export')
    previous = json.loads(reference_manifest.read_text())
    current = json.loads((archive/'manifest.json').read_text())
    execution = json.loads((old_archive/'execution.json').read_text())
    for key in ('variant','w_treatment','production_scale','decay_bottom_mass',
                'decay_scale_grouping','topology_hashes','internal_width_hashes'):
        if current[key] != previous[key]:
            raise ValueError('Grid physical/topology mismatch: '+key)
    if current.get('production_sampling', 'flat') != previous.get('production_sampling', 'flat'):
        raise ValueError('Grid sampling-map mismatch')
    if current.get('production_sampling_parameters') != previous.get('production_sampling_parameters'):
        raise ValueError('Grid sampling-parameter mismatch')
    for filename in ('param_card.dat','decay_card.dat','FO_analyse_card.dat'):
        expected = previous['hashes'][filename]
        if digest(old_archive/filename) != expected or digest(archive/filename) != expected:
            raise ValueError('Grid input-card mismatch: '+filename)
    excluded = {'iseed','req_acc_fo','npoints_fo','niters_fo','npoints_fo_grid',
                'niters_fo_grid','fo_job_target_time'}
    physical = lambda manifest: {key:value for key,value in manifest['settings'].items() if key not in excluded}
    if physical(previous) != physical(current):
        raise ValueError('Grid production input mismatch')
    scheduling = {'madgraph/interface/amcatnlo_run_interface.py','ttw_study/scripts/campaign.py',
                  'Template/fNLO/SubProcesses/ajob_template'}
    for source,checksum in execution['source_hashes'].items():
        if source not in scheduling and digest(ROOT/source) != checksum:
            raise ValueError('Grid physics source changed: '+source)
    for source in ('FixedOrderAnalysis/HwU.f90',
                   'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90','SubProcesses/mint_module.f90',
                   'SubProcesses/phase_space_kinematics.f90',
                   'SubProcesses/factorized_block_kinematics.f90',
                   'SubProcesses/decay_chain_parameters.f90',
                   'SubProcesses/decay_chain_kinematics.f90',
                   'SubProcesses/nlo_decay_kinematics.f90'):
        if digest(old_process/source) != digest(process/source):
            raise ValueError('Exported grid runtime changed: '+source)
    paths = [Path('SubProcesses/job_status.pkl')]
    for grid in sorted(old_process.glob('SubProcesses/P*/all_G*/res_0.dat')):
        directory = grid.parent
        if directory.name != 'all_G1':
            raise ValueError('Grid recovery currently supports one unsplit channel group per subprocess')
        required = {'mint_grids','grid.MC_integer','res_0.dat','res.dat','results.dat',
                    'contribution_results_0.dat','contribution_results.dat','log.txt'}
        regular = {path.name:path for path in directory.iterdir()
                   if path.is_file() and not path.is_symlink()}
        if not required <= regular.keys() or any(directory.glob('res_[1-9]*.dat')):
            raise ValueError('Recovery requires a complete training-only directory')
        # The native resume path also reads HTML bookkeeping (results.dat).
        # Copy every regular training artifact except plots, not symlinks to
        # old inputs. Worker startup regenerates links in the fresh export.
        paths.extend(path.relative_to(old_process) for name,path in sorted(regular.items())
                     if not name.endswith(('.HwU','.top','.lhe','.lhe.gz')))
    if len(paths) < 6:
        raise ValueError('No completed adaptive grids to recover')
    for relative in paths:
        if not (old_process/relative).is_file() or (process/relative).exists():
            raise ValueError('Missing source or existing destination grid: '+str(relative))
    checksums = {}
    for relative in paths:
        destination = process/relative
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(old_process/relative,destination)
        checksums[str(relative)] = digest(destination)
        assert checksums[str(relative)] == digest(old_process/relative)
    return dict(manifest=str(reference_manifest),manifest_sha256=digest(reference_manifest),
                files=checksums,mode='native saved-grid resume with independent refinement; no final weights reused')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('preflight')
    p.set_defaults(func=preflight)
    p = sub.add_parser('generate')
    p.add_argument('--charge', choices=('plus', 'minus'), required=True)
    p.add_argument('--w-treatment', choices=('onshell', 'top-bw', 'all-bw'), default='onshell')
    p.add_argument('--flavours', nargs=3, choices=('e', 'mu'), default=['e','e','mu'])
    p.add_argument('--corrected', choices=('both', 't', 'tbar', 'neither'), default='both')
    p.add_argument('--decay-bottom-mass', type=float, default=0.)
    p.add_argument('--tag', default='')
    p.set_defaults(func=generate)
    p = sub.add_parser('run')
    p.add_argument('--process-dir', required=True)
    p.add_argument('--variant', choices=('LO','P','D','S','PiD','Pi'), required=True)
    p.add_argument('--seed', type=int, required=True)
    p.add_argument('--grid-reference',type=Path,
                   help='recover a compatible failed run grid into a fresh export after a scheduler-only fix')
    p.add_argument('--points', type=int, default=1000)
    p.add_argument('--grid-points', type=int, default=200)
    p.add_argument('--iterations', type=int, default=3)
    p.add_argument('--accuracy', type=float, default=-1.)
    p.add_argument('--job-seconds', type=float, default=60.)
    p.add_argument('--main-run', action='store_true',
                   help='label as a main-campaign candidate; validation is still required')
    p.add_argument('--production-scale', choices=('core-w-ht-half','core-ht-half','fixed'),
                   default='core-w-ht-half')
    p.add_argument('--production-sampling', choices=('flat','w-current'), default='flat',
                   help='all-bw integration proposal only; fresh grids are required')
    p.add_argument('--ecm', type=float, default=13000.)
    p.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
