#!/usr/bin/env python3
"""Bounded, central-weight diagnostic integrations in an isolated ttW export.

This is a tail/stability pilot for one beam ordering, not a cross-section
campaign. Workers start from a saved proposal and adapt during fixed-count
iterations; a fresh training arm can start without grids. Limits apply to
each execute invocation. Original exports and campaign records are never modified.
"""
import argparse
import asyncio
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time

from campaign import ROOT, STUDY, lhapdf_environment, now

DEFAULT_TAG = 'tailpilot_20261014'
ORIGINAL = STUDY/'processes/TTWplus_onshell_muemu_mb0p0_both_precision_onshell_v3'
SUBPROCESS = 'P0_dxu_ttxwp_t_wpb_wp_mupvm_tx_wmbx_wm_emvex_wp_mupvm'
RUNTIME_FILES = ('FKS_params.dat', 'configs_and_props_info.dat',
    'decay_chain_info.dat', 'leshouche_info.dat', 'nlo_contribution_info.dat',
    'nlo_decay_info_2.dat', 'nlo_decay_info_3.dat', 'param_card.dat', 'symfact.dat')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def record_path(tag):
    if not re.fullmatch(r'[A-Za-z0-9_]+', tag):
        raise ValueError('Tag must contain only letters, numbers and underscores')
    return STUDY/'inputs'/('tail_pilot_'+tag+'.json')


def ignore_copy(directory, names):
    excluded = {'Events', 'HTML', 'study_cards', 'ME5_debug', 'ME5_debug.old'}
    return [name for name in names if name in excluded or name.startswith('all_G')]


def central_run_cards(process):
    """Keep readable cards consistent with compiled central-only controls."""
    for name in ('run_card.dat', 'run_card_default.dat'):
        path = process/'Cards'/name
        if not path.is_file():
            continue
        text = path.read_text()
        for field in ('reweight_scale', 'reweight_pdf'):
            text, count = re.subn(r'^([ \t]*)(?:True|False|\.true\.|\.false\.)'
                r'([ \t]*=[ \t]*'+field+r'\b)', r'\1False\2', text,
                flags=re.MULTILINE|re.IGNORECASE)
            if count == 0:
                text += '\nFalse = '+field+' ! Central-only diagnostic pilot.\n'
        path.write_text(text)


def prepare(args):
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        return prepare_isolated(args)


def prepare_isolated(args):
    path = record_path(args.tag)
    process = STUDY/'processes'/('TTWplus_onshell_muemu_mb0p0_both_'+args.tag)
    if path.exists() or process.exists():
        raise ValueError('Pilot output already exists: '+str(process))
    original_subprocess = ORIGINAL/'SubProcesses'/SUBPROCESS
    source_grid = original_subprocess/'all_G1'
    shutil.copytree(ORIGINAL, process, symlinks=True, ignore=ignore_copy)
    proposal = process/'pilot/proposal'
    proposal.mkdir(parents=True)
    provenance = {}
    for name in ('mint_grids', 'grid.MC_integer'):
        source = source_grid/name
        shutil.copy2(source, proposal/name)
        provenance[name] = dict(source=str(source), sha256=digest(proposal/name),
            original_mtime_ns=source.stat().st_mtime_ns, bytes=source.stat().st_size)
    # These inputs are compiled into setrun. The clone keeps the same physical
    # scales, PDF member and matrix elements; only variation work is disabled.
    include = process/'Source/run_card.inc'
    text = include.read_text()
    for field in ('LSCALEVAR', 'LPDFVAR'):
        text, count = re.subn(r'('+field+r'\(1\)\s*=\s*)\.TRUE\.',
                             r'\1.FALSE.', text, flags=re.IGNORECASE)
        if count != 1:
            raise ValueError('Expected one active '+field+' control in clone')
    include.write_text(text)
    central_run_cards(process)
    card = process/'Cards/decay_card.dat'
    text = re.sub(r'^\s*INDEPENDENT\s*=\s*decay_scale_variation_mode',
                  'NONE = decay_scale_variation_mode', card.read_text(), flags=re.MULTILINE)
    card.write_text(text)
    from tail_pilot_instrument import instrument
    instrumentation = instrument(process)
    data = dict(schema='bounded_ttw_tail_pilot_v1', tag=args.tag, created_utc=now(),
        status='building', process=str(process), subprocess=SUBPROCESS,
        original_process=str(ORIGINAL), proposal=str(proposal), grid_provenance=provenance,
        proposal_description='Snapshot of the original export all_G1 grids at preparation; '
            'grid_provenance identifies the copied files. Default workers share this snapshot.',
        physics='W+ mu,e,mu; m_b=0; all Ws and tops on shell; production and top decays NLO; '
            'central scales and central NNPDF40_nlo_as_01180 member only; dxu beam ordering only.',
        instrumentation=instrumentation, executions=[], build_log=str(process/'pilot/build.log'))
    save(path, data)
    start = time.monotonic()
    with (process/'pilot/build.log').open('w') as stream:
        # setrun_bridge reads run_card.inc into libgeneric.a; rebuilding only
        # the subprocess would silently retain the old reweight controls.
        result = subprocess.run(['make', '-j8'], cwd=process/'Source',
            stdout=stream, stderr=subprocess.STDOUT, env=lhapdf_environment(), timeout=1200)
        if result.returncode == 0:
            result = subprocess.run(['make', '-j8', 'madevent_mintFO'],
                cwd=process/'SubProcesses'/SUBPROCESS, stdout=stream, stderr=subprocess.STDOUT,
                env=lhapdf_environment(), timeout=1200)
    data.update(status='prepared' if result.returncode == 0 else 'build failed',
        build_returncode=result.returncode, build_elapsed_s=time.monotonic()-start)
    save(path, data)
    print(path, flush=True)
    if result.returncode:
        raise RuntimeError('Build failed: '+data['build_log'])


def worker_input(points, iterations=1, restart=-1):
    return '\n'.join((f'NPOINTS = {points}', f'NITERATIONS = {iterations}', 'ACCURACY = -1',
        'ADAPT_GRID = 2', 'MULTICHANNEL = 1', 'SUM_HELICITY = 1', 'NCHANS = 4',
        'CHANNEL = 1.1 1.2 2.1 2.2', 'SPLIT = 1', 'WGT_MULT= 1.0',
        'RUN_MODE = all', f'RESTART = {restart}', ''))


def set_virtual_fraction(worker, fraction):
    """Change a worker-local sampling proposal/floor, never shared inputs."""
    grid = worker/'mint_grids'
    lines = grid.read_text().splitlines()
    previous = []
    for index, line in enumerate(lines):
        words = line.split()
        if len(words) == 3 and words[0] == 'AVE':
            previous.append(float(words[1].replace('D', 'e').replace('d', 'e')))
            lines[index] = 'AVE %.17g %s' % (fraction, words[2])
    if len(previous) != 4:
        raise ValueError('Expected four terminal virtual-fraction grid entries')
    grid.write_text('\n'.join(lines)+'\n')
    path = worker/'FKS_params.dat'
    text = path.read_text()
    value = ('%.17e' % fraction).replace('e', 'd')
    for field in ('VirtualFraction', 'MinVirtualFraction'):
        text, count = re.subn(r'(^#'+field+r'[ \t]*\n)[^\n]+',
            lambda match: match.group(1)+value, text, flags=re.MULTILINE)
        if count != 1:
            raise ValueError('Expected one active '+field+' parameter')
    if path.is_symlink():
        path.unlink()
    path.write_text(text)
    return previous


def make_worker(process, proposal, label, variant, replica, seed, points,
                iterations=1, restart=-1, virtual_fraction=None):
    parent = process/'SubProcesses'/SUBPROCESS
    worker = parent/f'tail_{label}_{variant}_{replica:02d}_{seed}'
    worker.mkdir()
    for name in RUNTIME_FILES:
        (worker/name).symlink_to('../'+name)
    for name in ('mint_grids', 'grid.MC_integer'):
        shutil.copy2(proposal/name, worker/name)
    text = (process/'Cards/decay_card.dat').read_text()
    text, count = re.subn(r'^\s*(?:ADDITIVE|MULTIPLICATIVE)\s*=\s*nlo_decay_combination',
        ('ADDITIVE' if variant == 'S' else 'MULTIPLICATIVE')+' = nlo_decay_combination',
        text, flags=re.MULTILINE)
    if count != 1:
        raise ValueError('Missing decay combination control')
    (worker/'decay_card.dat').write_text(text)
    (worker/'randinit').write_text(f'r={seed}\n')
    (worker/'input_app.txt').write_text(worker_input(points, iterations, restart))
    previous = set_virtual_fraction(worker, virtual_fraction) if virtual_fraction is not None else None
    return dict(variant=variant, replica=replica, seed=seed, points=points,
        iterations=iterations, restart=restart, proposal=str(proposal),
        virtual_fraction=virtual_fraction, previous_grid_virtual_fractions=previous,
        refinement_stage=1, split=1, weight_multiplier=1., directory=str(worker),
        status='pending')


def runtime_provenance(process):
    """Small per-execution identity; no previous-campaign/archive scans."""
    binary = process/'SubProcesses'/SUBPROCESS/'madevent_mintFO'
    helper = Path(__file__).with_name('tail_pilot_instrument.py')
    sources = ('SubProcesses/driver_mintFO.f90', 'SubProcesses/driver_mintFO_bridge.f',
               'SubProcesses/fks_weights.f90', 'Source/run_card.inc')
    helper_hash = digest(helper) if helper.is_file() else None
    return dict(binary=str(binary), binary_sha256=digest(binary),
        instrumentation_version='tail_pilot_instrument.py:'+str(helper_hash),
        instrumentation_helper_sha256=helper_hash,
        source_sha256={name:digest(process/name) for name in sources if (process/name).is_file()})


def write_worker_manifest(job, runtime, grid_provenance):
    worker = Path(job['directory'])
    data = dict(schema='ttw_tail_pilot_worker_v1', created_utc=now(),
        inputs={key:value for key,value in job.items() if key != 'status'},
        runtime=runtime, grid_provenance=grid_provenance,
        input_sha256={name:digest(worker/name) for name in ('input_app.txt', 'randinit',
            'decay_card.dat', 'mint_grids', 'grid.MC_integer', 'FKS_params.dat')})
    path = worker/'pilot_manifest.json'
    save(path, data)
    job['manifest'] = str(path)


async def run_workers(path, data, execution, parallel, timeout):
    semaphore = asyncio.Semaphore(parallel)
    active = {}
    environment = lhapdf_environment()
    environment.update(MG5AMC_REFINEMENT_STAGE='1', OMP_NUM_THREADS='1',
        OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    process = Path(data['process'])
    binary = process/'SubProcesses'/SUBPROCESS/'madevent_mintFO'

    async def launch(job):
        async with semaphore:
            start = time.monotonic()
            worker = Path(job['directory'])
            with (worker/'input_app.txt').open('r') as source, (worker/'log.txt').open('w') as log:
                child = await asyncio.create_subprocess_exec(str(binary), cwd=worker,
                    stdin=source, stdout=log, stderr=asyncio.subprocess.STDOUT,
                    env=environment, start_new_session=True)
                active[child.pid] = child
                job.update(status='running', child_pid=child.pid, started_utc=now())
                save(path, data)
                print(f"Started {job['variant']} seed={job['seed']} pid={child.pid}", flush=True)
                try:
                    code = await asyncio.wait_for(child.wait(), timeout)
                    job['status'] = 'finished' if code == 0 and (worker/'res.dat').is_file() else 'failed'
                except asyncio.TimeoutError:
                    os.killpg(child.pid, signal.SIGTERM)
                    try:
                        await asyncio.wait_for(child.wait(), 5)
                    except asyncio.TimeoutError:
                        os.killpg(child.pid, signal.SIGKILL)
                        await child.wait()
                    code = child.returncode
                    job['status'] = 'timed out'
                finally:
                    active.pop(child.pid, None)
                job.update(returncode=code, elapsed_s=time.monotonic()-start, finished_utc=now())
                save(path, data)
                print(f"{job['status']}: {job['variant']} seed={job['seed']} "
                      f"elapsed={job['elapsed_s']:.1f}s", flush=True)

    try:
        await asyncio.gather(*(launch(job) for job in execution['jobs']))
    finally:
        for pid, child in active.items():
            try:
                os.killpg(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            await child.wait()
        for job in execution['jobs']:
            if job['status'] == 'running':
                job.update(status='interrupted', finished_utc=now())
        save(path, data)


def execute(args):
    path = record_path(args.tag)
    variants = ('S', 'Pi') if args.variant == 'both' else (args.variant,)
    total = len(variants)*args.replicas
    if not re.fullmatch(r'[A-Za-z0-9_]+', args.label):
        raise ValueError('Invalid execution label')
    if (not 1 <= total <= 16 or not 1 <= args.points <= 1000000
            or not 1 <= args.iterations <= 4 or total*args.points*args.iterations > 4000000):
        raise ValueError('Pilot cap: <=16 workers, <=1M points/iteration, <=4 iterations, <=4M requested calls')
    if args.virtual_fraction is not None and (not math.isfinite(args.virtual_fraction)
            or not 0 < args.virtual_fraction <= 1):
        raise ValueError('Virtual sampling fraction must be in (0,1]')
    if not 1 <= args.workers <= 16 or not 1 <= args.timeout <= 1200:
        raise ValueError('Pilot cap: <=16 concurrent workers and <=1200 seconds/worker')
    if args.seed_start < 1 or args.seed_start+total >= 900000000:
        raise ValueError('Invalid seed range')
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        data = json.loads(path.read_text())
        if data['status'] not in ('prepared', 'checkpoint'):
            raise ValueError('Pilot must be prepared or at a completed checkpoint')
        if any(run['label'] == args.label for run in data['executions']):
            raise ValueError('Execution label already exists; choose a fresh label')
        used_seeds = {job['seed'] for run in data['executions'] for job in run['jobs']}
        requested_seeds = set(range(args.seed_start, args.seed_start+total))
        if used_seeds & requested_seeds and not args.replay:
            raise ValueError('Seeds already used by this pilot; choose fresh seeds or explicit --replay')
        execution = dict(label=args.label, started_utc=now(), status='running',
            parent_pid=os.getpid(), concurrency=args.workers, timeout_s=args.timeout,
            requested_calls=total*args.points*args.iterations,
            same_stream_replay=bool(args.replay), jobs=[])
        process = Path(data['process'])
        proposal = (args.proposal or Path(data['proposal'])).resolve()
        grid_provenance = {name:dict(source=str(proposal/name), sha256=digest(proposal/name))
            for name in ('mint_grids', 'grid.MC_integer')}
        execution.update(proposal=str(proposal), grid_provenance=grid_provenance,
            virtual_fraction=args.virtual_fraction, iterations=args.iterations, restart=args.restart)
        execution['runtime'] = runtime_provenance(process)
        for variant in variants:
            for replica in range(args.replicas):
                seed = args.seed_start+len(execution['jobs'])
                job = make_worker(process, proposal, args.label, variant,
                    replica, seed, args.points, args.iterations, args.restart, args.virtual_fraction)
                write_worker_manifest(job, execution['runtime'], grid_provenance)
                execution['jobs'].append(job)
        data['executions'].append(execution)
        data['status'] = 'running'
        save(path, data)
        start = time.monotonic()
        try:
            asyncio.run(run_workers(path, data, execution, args.workers, args.timeout))
            execution['status'] = ('finished' if all(j['status'] == 'finished'
                for j in execution['jobs']) else 'failed')
        except BaseException as error:
            execution.update(status='interrupted', error=repr(error))
            raise
        finally:
            execution.update(elapsed_s=time.monotonic()-start, finished_utc=now())
            data['status'] = 'checkpoint'
            save(path, data)
        if execution['status'] != 'finished':
            raise RuntimeError('Diagnostic worker failed; inspect its local log')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag', default=DEFAULT_TAG)
    sub = parser.add_subparsers(dest='action', required=True)
    sub.add_parser('prepare').set_defaults(func=prepare)
    command = sub.add_parser('execute')
    command.add_argument('--label', default='initial')
    command.add_argument('--points', type=int, default=100000)
    command.add_argument('--replicas', type=int, default=8)
    command.add_argument('--variant', choices=('S', 'Pi', 'both'), default='both')
    command.add_argument('--seed-start', type=int, default=910001)
    command.add_argument('--workers', type=int, default=8)
    command.add_argument('--timeout', type=int, default=1200)
    command.add_argument('--proposal', type=Path,
        help='experimental initial folder containing mint_grids and grid.MC_integer')
    command.add_argument('--virtual-fraction', type=float,
        help='experimental initial virtual fraction and minimum floor, in worker-local files')
    command.add_argument('--iterations', type=int, default=1)
    command.add_argument('--restart', type=int, choices=(0, -1), default=-1,
        help='-1 reads the supplied proposal; 0 starts a fresh adaptive training grid')
    command.add_argument('--replay', action='store_true',
        help='label intentional same-seed/grid replays; these are not independent samples')
    command.set_defaults(func=execute)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
