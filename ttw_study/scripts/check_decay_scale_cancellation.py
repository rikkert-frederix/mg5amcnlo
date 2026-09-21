#!/usr/bin/env python3
"""Archive the paired inclusive decay-scale test without a stable-W rescaling.

Applicable to the full associated current as well as an explicit on-shell W.
This uses the existing matched BW-width proof and joint-batch implementation.
"""
import argparse
from pathlib import Path

from campaign import digest, now, save
from inclusive_report import check_joint_decay_scales


def run(batches,output):
    if output.exists():
        raise ValueError('Refuse to overwrite a paired decay-scale diagnostic')
    report=check_joint_decay_scales(batches)
    pulls=[abs(p) for p in report['pulls'] if p is not None]
    report.update(created_utc=now(),script_sha256=digest(Path(__file__)),
                  tested_nontrivial_points=len(pulls),maximum_absolute_conditional_pull=max(pulls),
                  status='paired inclusive decay-scale pilot; convergence and NLO normalization remain separate')
    save(output,report)
    print(len(pulls),'nontrivial points; maximum absolute conditional pull',max(pulls),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('batches',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    run(args.batches,args.output)
