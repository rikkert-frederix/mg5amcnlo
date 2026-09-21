#!/usr/bin/env python3
"""Compile changed kernels separately; never mutate an archived export."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import shutil
from campaign import ROOT, STUDY, digest, save, now


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('process', type=Path)
    args = parser.parse_args()
    process = args.process.resolve()
    if process.parent != STUDY/'processes':
        raise ValueError('Expected a local study export')
    build = Path(tempfile.mkdtemp(prefix='born_alignment_', dir=STUDY/'local'))
    template = ROOT/'Template/fNLO/SubProcesses'
    records = []
    for subdir in sorted(process.glob('SubProcesses/P*')):
        if not (subdir/'nlo_decay_info_2.dat').is_file():
            continue
        local = build/subdir.name
        local.mkdir()
        flags = ['gfortran', '-O0', '-g', '-fcheck=all', '-ffree-line-length-none',
                 '-ffunction-sections', '-I'+str(local), '-I'+str(subdir),
                 '-I'+str(process/'lib')]
        sources = [template/(n+'.f90') for n in
                   ('nlo_decay_kinematics', 'spin_density_matrix_results',
                    'multiplicative_nlo_decay', 'driver_mintFO')]
        source = STUDY/'scripts/check_born_alignment.f90'
        for path in sources+[source,Path(__file__).resolve()]:
            shutil.copy2(path,local/path.name)
        result = subprocess.run(flags+['-c', *[str(local/p.name) for p in sources]], cwd=local,
                                text=True, capture_output=True)
        (local/'compile.log').write_text(result.stdout+result.stderr)
        result.check_returncode()
        names = ('fnlo_process_common', 'process_dimensions', 'process_dimensions_bridge',
                 'fks_metadata', 'fks_metadata_bridge', 'nlo_contribution_bundle',
                 'nlo_decay_metadata', 'decay_chain_metadata', 'decay_chain_parameters',
                 'factorized_phase_space', 'factorized_block_kinematics',
                 'phase_space_kinematics', 'boostwdir2', 'decay_chain_kinematics')
        objects = [subdir/(n+'.o') for n in names]
        libraries = [process/'lib'/('lib'+n+'.a') for n in ('model','generic','pdf','dhelas')]
        result = subprocess.run(flags+[str(local/source.name),str(local/'nlo_decay_kinematics.o'),
            str(local/'multiplicative_nlo_decay.o'), str(local/'spin_density_matrix_results.o'),
            *map(str,objects), '-Wl,--gc-sections', '-Wl,--start-group',
            *map(str,libraries), '-Wl,--end-group', '-L'+str(ROOT/'LHAPDF/lib'),
            '-lLHAPDF','-lstdc++','-o',str(local/'check')], cwd=local,
            text=True,capture_output=True)
        (local/'link.log').write_text(result.stdout+result.stderr)
        if result.returncode:
            print(result.stderr)
        result.check_returncode()
        result = subprocess.run(['taskset','-c','0',str(local/'check')],cwd=subdir,
                                text=True,capture_output=True)
        (local/'run.log').write_text(result.stdout+result.stderr)
        print(result.stdout+result.stderr)
        result.check_returncode()
        if 'BORN_ALIGNMENT_PASS' not in result.stdout:
            raise ValueError('Missing positive runtime sentinel')
        records.append(dict(subprocess=subdir.name,output=result.stdout,
            source_snapshots={str(local/p.name):digest(local/p.name) for p in sources+[source,Path(__file__).resolve()]},
            linked_objects={str(p):digest(p) for p in objects},
            metadata={str(p):digest(p) for p in subdir.glob('*decay*info*.dat')}))
        for mode, expected in [('bad_point','not aligned with its snapshot'),
                               ('bad_order','Born alignment permutation is invalid')]:
            rejected = subprocess.run(['taskset','-c','0',str(local/'check'),mode],
                cwd=subdir,text=True,capture_output=True)
            (local/(mode+'.log')).write_text(rejected.stdout+rejected.stderr)
            if rejected.returncode == 0 or expected not in rejected.stdout:
                raise ValueError('Alignment guard did not reject '+mode)
    if not records:
        raise ValueError('No matching subprocesses')
    save(build/'result.json',dict(created_utc=now(),process=str(process),records=records,
        sources={str(p):digest(p) for p in sources+[source]},
        limitation='Deterministic Born-map checks, not a full MC product audit.'))
    print(build)


if __name__ == '__main__':
    main()
