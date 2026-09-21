"""Audit actual worker seeds, including the samples used to learn proposals.

Distinct initialization pairs are a necessary stream-separation check, not a
proof of PRNG independence, rare-tail convergence or adaptive-stop coverage.
"""
import hashlib
from pathlib import Path
import re

from campaign import ROOT, digest

RNG_SOURCE='Source/ranmar.f90'
MINT_SOURCE='SubProcesses/mint_module.f90'
STAGE_FILE=re.compile(r'Events/[^/]+/alllogs_(\d+)\.html')
STAGE_ROWS=re.compile(r'<a name=(/[^<>\s]+)></a>.*?<PRE>\n(.*?)\n</PRE>',re.S)


def history_paths(process,run):
    paths=[process/RNG_SOURCE,process/MINT_SOURCE]
    logs=sorted((process/'Events'/run).glob('alllogs_*.html'))
    if not logs or any(not STAGE_FILE.fullmatch(str(p.relative_to(process))) for p in logs):
        raise ValueError('Missing or malformed native per-stage log collection')
    return paths+logs


def refinement_steps(log):
    steps=[int(s) for s in re.findall(r'Refining results, step\s+(-?\d+)',log)]
    if not steps:
        raise ValueError('Missing native refinement history for the batch-independence audit')
    return steps


def validate_refinement_history(log):
    """Fail closed on legacy archives without verified per-stage evidence."""
    steps=refinement_steps(log)
    if max(steps)>1 or len(set(steps))>1:
        raise ValueError('Repeated refinement rounds reuse native random streams; '
                         'independent conditional batch covariance is not certified')


def validate_history(log,seed,raw_files,final_logs=None):
    steps=refinement_steps(log)
    if 'Setting up grids' not in log or min(steps)<1:
        raise ValueError('Stage audit currently requires fresh in-run training, not imported grids')
    for source in (RNG_SOURCE,MINT_SOURCE):
        if (source not in raw_files or hashlib.sha256(raw_files[source]).hexdigest()
                !=digest(ROOT/'Template/fNLO'/source)):
            raise ValueError('Per-stage source evidence differs from the audited RNG/estimator')
    expected={0,*steps}
    stage_files={int(match[1]):raw for name,raw in raw_files.items()
                 if (match:=STAGE_FILE.fullmatch(name))}
    if set(stage_files)!=expected or len(stage_files)!=len(raw_files)-2:
        raise ValueError('Incomplete or duplicate per-stage log evidence')
    seen={}
    last={}
    counts={}
    for stage,raw in sorted(stage_files.items()):
        text=raw.decode()
        rows=STAGE_ROWS.findall(text)
        if not rows or len(rows)!=text.count('<PRE>') or len(rows)!=text.count('<a name='):
            raise ValueError('Malformed native stage log collection')
        directories=set()
        for directory,worker_log in rows:
            directory='SubProcesses'+directory
            if directory in directories:
                raise ValueError('Duplicate worker in one native stage')
            directories.add(directory)
            base=re.findall(r'with seed\s+(\d+)',worker_log)
            pairs=re.findall(r'Ranmar initialization seeds\s+(\d+)\s+(\d+)',worker_log)
            namespace=re.findall(r'Ranmar stream namespace refinement-v1 stage (\d+) split (\d+)',worker_log)
            if (base!=[str(seed)] or len(pairs)!=1 or len(namespace)!=1
                    or int(namespace[0][0])!=stage or 'Time spent in Total :' not in worker_log):
                raise ValueError('Missing, incomplete or inconsistent stage-specific worker stream')
            suffix=re.fullmatch(r'(?:all|born)_G\d+(?:_(\d+))?',Path(directory).name)
            if not suffix or int(suffix[1] or 0)!=int(namespace[0][1]):
                raise ValueError('Stage log split namespace disagrees with worker directory')
            pair=tuple(map(int,pairs[0]))
            if not (0<=pair[0]<31329 and 0<=pair[1]<30082):
                raise ValueError('Invalid RANMAR initialization pair')
            if pair in seen:
                raise ValueError('Random-stream collision across workers/stages: %s and %s' %
                                 (seen[pair],(stage,directory)))
            seen[pair]=(stage,directory)
            last[directory]=dict(stage=stage,pair=pair,
                                 log_sha256=hashlib.sha256(worker_log.encode()).hexdigest())
        counts[str(stage)]=len(rows)
    for directory,checksum in (final_logs or {}).items():
        if directory not in last or last[directory]['stage']==0 or last[directory]['log_sha256']!=checksum:
            raise ValueError('Final worker is not its last complete archived refinement')
    return dict(scheme='refinement-v1',stage_worker_counts=counts,
                distinct_initialization_pairs=len(seen),final_workers_checked=len(final_logs or {}),
                limitations='Actual initialization pairs checked across fresh training and all refinements. '
                            'Not proof of PRNG independence, retraining, rare-tail convergence, '
                            'or unbiased adaptive stopping.')
