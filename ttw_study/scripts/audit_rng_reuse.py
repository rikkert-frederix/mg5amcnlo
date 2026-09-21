#!/usr/bin/env python3
"""Freeze completed-run stage logs and diagnose exact RANMAR seed reuse.

This accepts legacy generators for diagnosis only; it does not certify a
statistical estimator or change any existing run/audit record.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
import tarfile

from campaign import STUDY, digest, now, save
from rng_history import STAGE_FILE, STAGE_ROWS, history_paths


def audit(audit_path,output):
    report=json.loads(audit_path.read_text())
    hwu=Path(report['path'])
    process=hwu.parents[2]
    if (process.parent!=STUDY/'processes' or report['execution']['status']!='finished'
            or digest(hwu)!=report['output_sha256']):
        raise ValueError('Require an unchanged completed study histogram/audit')
    if output.exists() or output.with_suffix('.json').exists():
        raise ValueError('Refuse to overwrite an RNG diagnostic')
    paths=history_paths(process,hwu.parent.name)
    streams=defaultdict(list)
    counts={}
    seed=report['manifest']['settings']['iseed']
    for path in paths:
        match=STAGE_FILE.fullmatch(str(path.relative_to(process)))
        if not match:
            continue
        stage=int(match[1])
        rows=STAGE_ROWS.findall(path.read_text())
        if not rows:
            raise ValueError('Empty native stage log')
        counts[str(stage)]=len(rows)
        for directory,log in rows:
            pairs=re.findall(r'Ranmar initialization seeds\s+(\d+)\s+(\d+)',log)
            if (len(pairs)!=1 or re.findall(r'with seed\s+(\d+)',log)!=[str(seed)]
                    or 'Time spent in Total :' not in log):
                raise ValueError('Incomplete or mismatched worker seed evidence')
            streams[tuple(map(int,pairs[0]))].append(dict(stage=stage,directory=directory))
    collisions=[dict(pair=pair,uses=uses) for pair,uses in streams.items() if len(uses)>1]
    output.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(output,'w:gz',compresslevel=1) as archive:
        for path in paths:
            archive.add(path,arcname=str(path.relative_to(process)),recursive=False)
    record=dict(created_utc=now(),audit=str(audit_path.resolve()),audit_sha256=digest(audit_path),
                archive=str(output.resolve()),archive_sha256=digest(output),seed=seed,
                files={str(p.relative_to(process)):digest(p) for p in paths},
                stage_worker_counts=counts,distinct_initialization_pairs=len(streams),
                repeated_pairs=len(collisions),collisions=collisions,
                status='completed-run RNG reuse diagnostic; no statistical coverage certified')
    save(output.with_suffix('.json'),record)
    print('%d repeated initialization pairs among %d stage-worker executions' %
          (len(collisions),sum(counts.values())),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('audit',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    audit(args.audit,args.output)
