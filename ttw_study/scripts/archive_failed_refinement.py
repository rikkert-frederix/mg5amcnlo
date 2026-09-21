#!/usr/bin/env python3
"""Preserve a frozen unsuccessful adaptive run, never a physics estimate.

Call only after stopping its exact calculator process tree. Native stage
rollover may already have replaced some worker histograms: report that loss
explicitly and retain every available per-stage summary and log.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import pickle
import re
import tarfile

from campaign import ROOT, STUDY, digest, now, save


def archive(process, name, output):
    process=process.resolve()
    if process.parent!=STUDY/'processes':
        raise ValueError('Only local study processes may be archived')
    if output.exists() or output.with_suffix('.json').exists():
        raise ValueError('Do not overwrite failed-run diagnostics')
    card_archive=process/'study_cards'/name
    execution=json.loads((card_archive/'execution.json').read_text())
    pid=execution['child_pid']
    proc=Path('/proc')/str(pid)
    if str(process/'bin/calculate_xsect') not in (proc/'cmdline').read_bytes().decode().replace('\0',' '):
        raise ValueError('Calculator PID does not match the exact study export')
    state=(proc/'status').read_text()
    if not re.search(r'^State:\s+T\s',state,re.M):
        raise ValueError('Freeze the calculator before copying its state')
    status=process/'SubProcesses/job_status.pkl'
    # This pickle is the locally generated state of the verified process.
    jobs=pickle.loads(status.read_bytes())
    if not jobs or any(j.get('niters_done')!=1 or j.get('split',0)<=0 for j in jobs):
        raise ValueError('Expected the completed first-refinement split ledger')
    paths={status}
    for directory in (process/'Cards',card_archive):
        paths.update(p for p in directory.rglob('*') if p.is_file() and not p.is_symlink())
    paths.update(p for p in (process/'SubProcesses').glob('*.f90') if not p.is_symlink())
    for directory in process.glob('SubProcesses/P*/all_G1*'):
        if directory.is_dir() and not directory.is_symlink():
            paths.update(p for p in directory.iterdir() if p.is_file() and not p.is_symlink())
    summaries=[]
    for job in jobs:
        directory=Path(job['dirname']).resolve()
        if not directory.is_relative_to(process/'SubProcesses'):
            raise ValueError('Job directory outside the verified export')
        requested=(directory/'input_app.txt').read_text()
        npoints=int(re.search(r'^NPOINTS\s*=\s*(\d+)',requested,re.M)[1])
        summaries.append(dict(directory=str(directory.relative_to(process)),
                              completed_npoints=job['npoints_done'],next_requested_npoints=npoints,
                              result=job['result'],error=job['error'],
                              absolute_result=job['resultABS'],absolute_error=job['errorABS'],
                              contribution_results=job.get('contribution_results'),
                              surviving_histogram=(directory/'MADatNLO.HwU').is_file(),
                              completed_stage_log=(directory/'log_MINT1.txt').is_file(),
                              completed_stage_summary=(directory/'res_1.dat').is_file()))
    paths=sorted(paths)
    checksums={str(p.relative_to(process)):digest(p) for p in paths}
    output.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(output,'w:gz',compresslevel=1) as stream:
        for p in paths:
            stream.add(p,arcname=str(p.relative_to(process)),recursive=False)
    if any(digest(process/p)!=sha for p,sha in checksums.items()):
        raise ValueError('State changed during the frozen snapshot')
    result=sum(j['result'] for j in jobs)
    error=math.sqrt(sum(j['error']**2 for j in jobs))
    record=dict(created_utc=now(),process=str(process),run=name,calculator_pid=pid,
                status='frozen nonconverged adaptive diagnostic; excluded from physics',
                reason='First refinement has order-one relative error; controller requested 100x points',
                source_execution=execution,worker_count=len(jobs),workers=summaries,
                completed_points=sum(j['npoints_done'] for j in jobs),
                next_requested_points=sum(j['next_requested_npoints'] for j in summaries),
                completed_stage_result_pb=result,completed_stage_error_pb=error,
                relative_error=abs(error/result),
                largest_variance_workers=sorted(summaries,key=lambda j:j['error']**2,reverse=True)[:10],
                surviving_histograms=sum(j['surviving_histogram'] for j in summaries),
                completed_stage_logs=sum(j['completed_stage_log'] for j in summaries),
                completed_stage_summaries=sum(j['completed_stage_summary'] for j in summaries),
                archive=str(output.resolve()),archive_sha256=digest(output),files=checksums,
                script_sha256=digest(Path(__file__)),
                limitations='Not a completed MG5 run or a final HwU result. Native second-stage startup '
                            'has replaced some first-stage histograms; no complete covariance or '
                            'fiducial prediction may be reconstructed from this partial snapshot.')
    save(output.with_suffix('.json'),record)
    print(json.dumps({k:record[k] for k in ('status','completed_points','next_requested_points',
          'completed_stage_result_pb','completed_stage_error_pb','relative_error','surviving_histograms',
          'completed_stage_logs','completed_stage_summaries','archive_sha256')},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--process',type=Path,required=True)
    parser.add_argument('--run',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    archive(args.process,args.run,args.output)
