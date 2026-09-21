#!/usr/bin/env python3
"""Preserve completed, still-present worker histograms before reconfiguration.

This archive enables a separate audit of split-level covariance. Merely
retaining splits does not certify them as independent integration replicas.
Never harvest a process with a live or subsequently reconfigured run.
"""
import argparse
import json
from pathlib import Path
import pickle
import re
import tarfile

from campaign import STUDY, digest, now, save
from rng_history import history_paths


def harvest(audit_path,output):
    report = json.loads(audit_path.read_text())
    hwu = Path(report['path'])
    process = hwu.parents[2]
    if process.parent != STUDY/'processes' or report['execution']['status'] != 'finished':
        raise ValueError('Require a completed study execution')
    if digest(hwu) != report['output_sha256']:
        raise ValueError('Final histogram checksum changed')
    if output.exists() or output.with_suffix('.json').exists():
        raise ValueError('Worker archive already exists')
    # This is the trusted, locally generated MG5 job-status pickle from the
    # verified study process, never a downloaded or user-supplied pickle.
    status = process/'SubProcesses/job_status.pkl'
    jobs = pickle.loads(status.read_bytes())
    seed = report['manifest']['settings']['iseed']
    copied = []
    files = []
    for job in jobs:
        directory = Path(job['dirname']).resolve()
        if not directory.is_relative_to(process/'SubProcesses') or job['split'] <= 0:
            raise ValueError('Expected completed local split-job directories')
        log = (directory/'log.txt').read_text()
        seeds = re.findall(r'with seed\s+(\d+)',log)
        if seeds != [str(seed)]:
            raise ValueError('Worker no longer belongs to this run: '+str(directory))
        if 'Time spent in Total :' not in log or not (directory/'MADatNLO.HwU').is_file():
            raise ValueError('Incomplete worker output')
        relative = directory.relative_to(process)
        selected = [directory/name for name in ('MADatNLO.HwU','input_app.txt','log.txt','res.dat')]
        checksums = {path.name:digest(path) for path in selected}
        files.extend(selected)
        copied.append(dict(directory=str(relative),job=job,files=checksums))
    history=history_paths(process,hwu.parent.name)
    history_files={str(p.relative_to(process)):digest(p) for p in history}
    files.extend(history)
    exclusions=sorted(p for p in (hwu.parent/'split_outliers').rglob('*') if p.is_file())
    outlier_files={str(p.relative_to(process)):digest(p) for p in exclusions}
    files.extend(exclusions)
    output.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(output,'w:gz',compresslevel=1) as archive:
        for i,path in enumerate(files):
            archive.add(path,arcname=str(path.relative_to(process)),recursive=False)
            if i % 64 == 63:
                print('Archived %d/%d worker files' % (i+1,len(files)),flush=True)
    record = dict(created_utc=now(),audit=str(audit_path.resolve()),audit_sha256=digest(audit_path),
                  run=hwu.parent.name,seed=seed,worker_count=len(copied),workers=copied,
                  history_files=history_files,split_outlier_files=outlier_files,
                  job_status_sha256=digest(status),archive_sha256=digest(output),
                  status='raw split outputs preserved; covariance interpretation not yet certified')
    save(output.with_suffix('.json'),record)
    print(output,flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('audit',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    harvest(args.audit,args.output)
