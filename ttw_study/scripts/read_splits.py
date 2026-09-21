#!/usr/bin/env python3
"""Audit archived equal-count split batches and preserve joint weight vectors.

The batches are independent conditional on the trained importance proposal,
with separate subprocess/channel strata. They are not independent full runs
with retraining. Their use still requires rare-tail/convergence diagnostics.
"""
import argparse
from collections import defaultdict
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tarfile

import numpy as np

from campaign import ROOT, digest, now, save
from load_results import canonical_columns, load
from madgraph.various.histograms import HwUList
from rng_history import validate_history, validate_refinement_history


def validate_jobs(workers):
    groups = defaultdict(list)
    for index, worker in enumerate(workers):
        job = worker['job']
        if job['mint_mode'] != -1 or job['niters'] != 1 or job['niters_done'] != 1:
            raise ValueError('Require exactly one complete iteration on a trained grid')
        if job['npoints_done'] != job['npoints'] or job['npoints'] <= 0:
            raise ValueError('Require the full fixed batch count without extra sampling')
        groups[(job['p_dir'],job['channel'])].append(index)
    for indices in groups.values():
        if len(indices) < 5:
            raise ValueError('Require at least five independent batches per stratum')
        jobs = [workers[i]['job'] for i in indices]
        splits = {job['split'] for job in jobs}
        metadata=jobs[0].get('split_outlier')
        if any(job.get('split_outlier')!=metadata for job in jobs):
            raise ValueError('Inconsistent split exclusion metadata')
        original_count=metadata['original_count'] if metadata else len(jobs)
        excluded=set(metadata['excluded_splits']) if metadata else set()
        if metadata and (original_count!=len(jobs)+1 or len(excluded)!=1 or not metadata.get('audit')):
            raise ValueError('Invalid single-split exclusion metadata')
        if len(splits)!=len(jobs) or splits & excluded or splits | excluded != set(range(1,original_count+1)):
            raise ValueError('Missing or duplicate split index')
        for key in ('npoints','configs','nchans','accuracy','run_mode','wgt_frac'):
            if any(job[key] != jobs[0][key] for job in jobs):
                raise ValueError('Unequal batch settings: '+key)
        if any(job['wgt_mult'] != 1./original_count or not math.isclose(
                job.get('split_result_scale',1.),original_count/len(jobs),rel_tol=1.e-12)
                for job in jobs):
            raise ValueError('Unexpected split weight normalization')
    if not groups:
        raise ValueError('Empty batch archive')
    return groups


def validate_exclusion_evidence(workers,evidence,process):
    """A missing split must be justified by the archived, checksum-verified audit."""
    for worker in workers:
        job=worker['job']
        metadata=job.get('split_outlier')
        if not metadata:
            continue
        audit=str(Path(metadata['audit']).relative_to(process))
        if audit not in evidence:
            raise ValueError('Missing archived split exclusion audit')
        matches=[g for g in evidence[audit]['groups'] if
                 (g['p_dir'],g['channel'],g['run_mode'])==(job['p_dir'],job['channel'],job['run_mode'])]
        if len(matches)!=1:
            raise ValueError('Ambiguous split exclusion audit')
        group=matches[0]
        raw=[r for r in group['raw_jobs'] if r['split']==job['split']]
        if (group['original_count']!=metadata['original_count'] or
                [group['excluded_split']]!=metadata['excluded_splits'] or len(raw)!=1 or
                not math.isclose(group['retained_scale'],job['split_result_scale'],rel_tol=1.e-12)):
            raise ValueError('Split exclusion audit differs from retained job')
        for key in ('result','error','resultABS','errorABS'):
            if not math.isclose(raw[0][key]*group['retained_scale'],job[key],rel_tol=1.e-12,abs_tol=1.e-30):
                raise ValueError('Retained result differs from split exclusion audit')


def split_ensembles(contributions, strata):
    """Convert already 1/N-weighted split contributions to stratum replicas.

    Sum the stratum MEANS for a cross section. Delete a batch from only one
    stratum at a time for nonlinear observables (independent_jackknife_many).
    """
    contributions = np.asarray(contributions,dtype=float)
    strata = np.asarray(strata)
    if len(contributions) != len(strata) or not np.isfinite(contributions).all():
        raise ValueError('Inconsistent or nonfinite batch vectors')
    result = []
    for label in sorted(set(strata)):
        selected = contributions[strata == label]
        if len(selected) < 5:
            raise ValueError('Too few batches in a stratum')
        result.append(len(selected)*selected)
    return result


def read(record_path, output):
    record_path, output = Path(record_path).resolve(), Path(output).resolve()
    if output.exists() or output.with_suffix('.json').exists():
        raise ValueError('Refuse to overwrite a batch audit')
    record = json.loads(record_path.read_text())
    archive_path = record_path.with_suffix('.gz')
    if digest(archive_path) != record['archive_sha256']:
        raise ValueError('Archive checksum mismatch')
    if digest(record['audit']) != record['audit_sha256']:
        raise ValueError('Combined-output audit changed')
    reference = load(record['audit'])
    run_log=Path(reference['report']['execution']['log']).read_text()
    if not record.get('history_files'):
        validate_refinement_history(run_log)
    workers = record['workers']
    groups = validate_jobs(workers)
    expected = {}
    directories = {worker['directory']:i for i,worker in enumerate(workers)}
    if len(directories) != len(workers):
        raise ValueError('Repeated worker directory')
    for worker in workers:
        for filename,checksum in worker['files'].items():
            expected[worker['directory']+'/'+filename] = checksum
    expected.update(record.get('history_files',{}))
    expected.update(record.get('split_outlier_files',{}))
    history={}
    exclusions={}
    vectors, errors, seeds, initial_channels = {}, {}, {}, {}
    # Streaming avoids random seeks/repeated decompression of this large file.
    with tarfile.open(archive_path,'r|gz') as archive:
        for member in archive:
            if not member.isfile() or member.name not in expected:
                raise ValueError('Unexpected archive entry: '+member.name)
            raw = archive.extractfile(member).read()
            if hashlib.sha256(raw).hexdigest() != expected.pop(member.name):
                raise ValueError('Worker checksum mismatch: '+member.name)
            if member.name in record.get('history_files',{}):
                history[member.name]=raw
                continue
            if member.name in record.get('split_outlier_files',{}):
                if member.name.endswith('/audit.json'):
                    exclusions[member.name]=json.loads(raw)
                continue
            directory,filename = member.name.rsplit('/',1)
            index = directories[directory]
            if filename == 'log.txt':
                log = raw.decode()
                found = re.findall(r'Ranmar initialization seeds\s+(\d+)\s+(\d+)',log)
                base = re.findall(r'with seed\s+(\d+)',log)
                if len(found) != 1 or base != [str(record['seed'])]:
                    raise ValueError('Unverified worker random stream')
                seeds[index] = tuple(map(int,found[0]))
                initial_channels[index] = re.findall(r'^channel.*$',log.split('------- iteration')[0],re.M)
                if not initial_channels[index]:
                    raise ValueError('Missing initial channel/proposal diagnostics')
            elif filename == 'MADatNLO.HwU':
                histograms = HwUList(io.StringIO(raw.decode()),raw_labels=True)
                if [hist.title for hist in histograms] != reference['titles']:
                    raise ValueError('Worker histogram layout changed')
                columns = list(histograms[0].bins[0].wgts)
                indices = canonical_columns(columns,reference['report']['variant'])
                rows = [row for hist in histograms for row in hist.bins]
                if not np.array_equal(np.asarray([row.boundaries for row in rows]),reference['edges']):
                    raise ValueError('Worker bin edges changed')
                data = np.asarray([[row.wgts[label] for label in columns] for row in rows])
                if not np.isfinite(data).all():
                    raise ValueError('Nonfinite worker weight')
                factor=workers[index]['job'].get('split_result_scale',1.)
                vectors[index] = factor*data[:,indices]
                errors[index] = factor*data[:,columns.index('dy')]
                if len(vectors) % 8 == 0:
                    print('Parsed %d/%d batch vectors' % (len(vectors),len(workers)),flush=True)
    if expected or len(vectors) != len(workers) or len(seeds) != len(workers):
        raise ValueError('Incomplete worker archive')
    if len(set(seeds.values())) != len(workers):
        raise ValueError('Random-stream collision between workers')
    process = Path(reference['report']['path']).parents[2]
    validate_exclusion_evidence(workers,exclusions,process)
    stage_audit=validate_history(run_log,record['seed'],history,
        {w['directory']:w['files']['log.txt'] for w in workers}) if history else None
    strata = np.empty(len(workers),dtype=int)
    for label,indices in enumerate(groups.values()):
        if any(initial_channels[i] != initial_channels[indices[0]] for i in indices):
            raise ValueError('Initial proposal/channel diagnostics differ within a stratum')
        strata[indices] = label
    values = np.asarray([vectors[i] for i in range(len(workers))])
    error = np.asarray([errors[i] for i in range(len(workers))])
    total = values.sum(axis=0)
    difference = total-reference['values']
    tolerance = 3.e-6*(np.abs(values).sum(axis=0)+np.abs(reference['values']))+1.e-20
    if not (np.abs(difference) <= tolerance).all():
        raise ValueError('Worker sum does not reproduce every final weight/bin')
    error_sum = np.sqrt(np.sum(error**2,axis=0))
    if not np.allclose(error_sum,reference['errors'],rtol=3.e-6,atol=1.e-20):
        raise ValueError('Independent worker error sum disagrees with final HwU')
    process = Path(reference['report']['path']).parents[2]
    sources = {}
    for relative in ('SubProcesses/mint_module.f90','Source/ranmar.f90'):
        if digest(process/relative) != digest(ROOT/'Template/fNLO'/relative):
            raise ValueError('Current source audit does not match generated estimator')
        sources[relative] = digest(process/relative)
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output,contributions=values,individual_errors=error,
                        strata=strata,random_seeds=np.asarray([seeds[i] for i in range(len(workers))]),
                        edges=reference['edges'],offsets=reference['offsets'],
                        titles=reference['titles'],weights=reference['weights'])
    save(output.with_suffix('.json'),dict(
        created_utc=now(),input=str(record_path),input_sha256=digest(record_path),
        arrays_sha256=digest(output),estimator_sources=sources,
        refinement_rng_audit=stage_audit,
        split_outlier_audits=exclusions,
        final_split_filtering_applied=any('split_outlier' in w['job'] for w in workers),
        worker_count=len(workers),strata={str(key):len(indices) for key,indices in groups.items()},
        maximum_weight_sum_difference_pb=float(np.abs(difference).max()),
        maximum_error_sum_difference_pb=float(np.abs(error_sum-reference['errors']).max()),
        status='equal-count conditional split-batch audit passed; full-run convergence not certified',
        convention='Per-stratum equal-weight batches; independent stratum variances add. '
                   'Distinct native RANMAR seeds and common trained proposals checked. '
                   'Not a claim of independent retraining or asymptotic coverage in sparse tails. '
                   'If split filtering was applied, errors condition on retention and omit selection bias.'))
    print(output,flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    read(args.record,args.output)
