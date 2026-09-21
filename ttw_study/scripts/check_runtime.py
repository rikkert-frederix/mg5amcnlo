#!/usr/bin/env python3
"""Validate couplings through the actual generated process's compiled libraries."""
import argparse
import json
import math
from pathlib import Path
import subprocess

from campaign import ROOT, STUDY, digest, now, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('process_dir', type=Path)
    args = parser.parse_args()
    process = args.process_dir.resolve()
    if process.parent != STUDY/'processes':
        raise ValueError('Not a study process')
    archive = STUDY/'inputs'/(process.name+'_runtime_coupling.json')
    if archive.exists():
        raise ValueError('Runtime check already archived')
    build = STUDY/'local/runtime_checks'/process.name
    build.mkdir(parents=True, exist_ok=False)
    source = STUDY/'scripts/check_runtime_coupling.f90'
    libraries = [process/'lib/libgeneric.a', process/'lib/libpdf.a']
    command = ['gfortran', '-O0', '-g', '-fcheck=all',
               '-ffpe-trap=invalid,zero,overflow', '-I'+str(process/'lib'),
               '-I'+str(process/'Source'), str(source), '-Wl,--start-group',
               *map(str, libraries), '-Wl,--end-group',
               '-L'+str(ROOT/'LHAPDF/lib'), '-lLHAPDF', '-lstdc++',
               '-o', str(build/'check_coupling')]
    compile_result = subprocess.run(command, cwd=build, text=True, capture_output=True)
    (build/'build.log').write_text(compile_result.stdout+compile_result.stderr)
    compile_result.check_returncode()
    result = subprocess.run(['taskset', '-c', '0-63', str(build/'check_coupling')],
                            cwd=process, text=True, capture_output=True)
    (build/'run.log').write_text(result.stdout+result.stderr)
    result.check_returncode()
    measured = [list(map(float, line.split()[1:])) for line in result.stdout.splitlines()
                if line.startswith('COUPLING')]
    expected = json.loads((STUDY/'inputs/benchmark.json').read_text())['coupling_checkpoints']
    assert len(measured) == len(expected) == 5
    rows = []
    for (q, value), reference in zip(measured, expected):
        assert q == reference['Q_GeV']
        assert math.isclose(value, reference['alpha_s'], rel_tol=0., abs_tol=5.e-14)
        rows.append(dict(Q_GeV=q, runtime_alpha_s=value, pdf_alpha_s=reference['alpha_s']))
    executables = sorted(process.glob('SubProcesses/P*/madevent_mintFO'))
    assert executables
    record = dict(created_utc=now(), status='passed', process=str(process),
                  compile_command=command, stdout=result.stdout, coupling_checkpoints=rows,
                  source_sha256=digest(source),
                  compiled_library_sha256={str(p): digest(p) for p in libraries},
                  integration_executable_sha256={str(p): digest(p) for p in executables},
                  note='Same compiled PDF initialization and alphas objects as integration; no generated-source edits.')
    save(archive, record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
