#!/usr/bin/env python3
"""Compile and probe unchanged factorized phase-space kernels, without MC.

This establishes local map support/Jacobians, not a generated narrow-W
cross-section limit or convergence of the production-current integral.
"""
import argparse
from pathlib import Path
import subprocess
import tempfile

from campaign import ROOT, digest, now, save


def run(output, process):
    if output.exists():
        raise ValueError('Refuse to overwrite a phase-space support record')
    fixture=ROOT/'tests/input_files/fks_decay'
    source=ROOT/'Template/fNLO'
    paths=[fixture/'phase_space_test_dimensions.f90',source/'Source/kin_functions.f90',
           source/'SubProcesses/phase_space_kinematics.f90',
           source/'SubProcesses/factorized_block_kinematics.f90',fixture/'bw_support_checks.f90']
    generated={}
    for relative in (() if process is None else ('SubProcesses/phase_space_kinematics.f90',
                     'SubProcesses/factorized_block_kinematics.f90',
                     'SubProcesses/decay_chain_kinematics.f90',
                     'SubProcesses/nlo_decay_kinematics.f90')):
        if digest(process/relative)!=digest(source/relative):
            raise ValueError('Generated map differs from audited source: '+relative)
        generated[relative]=digest(process/relative)
    with tempfile.TemporaryDirectory(prefix='ttw_bw_support_') as directory:
        executable=Path(directory)/'check'
        command=['gfortran','-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
                 '-o',str(executable),*map(str,paths)]
        compilation=subprocess.run(command,cwd=directory,capture_output=True,text=True)
        if compilation.returncode:
            raise ValueError(compilation.stdout+compilation.stderr)
        result=subprocess.run([str(executable)],cwd=directory,capture_output=True,text=True)
        if result.returncode or 'PASS: full-support' not in result.stdout:
            raise ValueError(result.stdout+result.stderr)
    values={}
    for line in result.stdout.splitlines():
        key,*rest=line.split()
        if len(rest)==1:
            values[key]=float(rest[0])
    record=dict(created_utc=now(),status='local support and measure probes passed; generated limits still required',
                process=str(process.resolve()) if process else None,generated_source_sha256=generated,
                source_only=process is None,
                source_sha256={str(path.relative_to(ROOT)):digest(path) for path in paths},
                script_sha256=digest(Path(__file__)),compiler_command=command,
                stdout=result.stdout,values=values,
                support='Associated current at partonic energies 500/1000/2000 GeV; '
                        'top BW current at mb=0/4.8 GeV and GammaW factors 1/0.05; '
                        'interior points near both physical endpoints and away from/on the W pole.',
                limitations='Finite point checks plus source inspection, not exhaustive end-to-end '
                            'event generation. The zero-measure endpoints themselves are not sampled. '
                            'No matrix elements, PDFs, cuts, resonance-selector test, '
                            'cross-section narrow-W limit or convergence certification.')
    save(output,record)
    print(result.stdout,end='')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    target=parser.add_mutually_exclusive_group(required=True)
    target.add_argument('--process',type=Path)
    target.add_argument('--source-only',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.output,args.process)
