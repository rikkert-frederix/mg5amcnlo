#!/usr/bin/env python3
"""Native-current batches with explicitly audited deterministic MINT rounding.

Unlike stopping on a measured error, initialize_even_random_numbers chooses
the effective count before sampling. Preserve both requested and actual counts.
The general study reader's exact-requested-count guard is not changed.
"""
import argparse
import copy
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

import numpy as np

from audit_native_current import native_vectors
from campaign import digest, now, save
from load_results import load
from madgraph.various.histograms import HwUList
from read_splits import validate_jobs
from rng_history import validate_history


def effective_count(requested,dimensions):
    if requested<2 or dimensions<1:
        raise ValueError('Unsupported native count/dimension')
    ng=int((requested/2.)**(1./dimensions))
    if ng<1:
        raise ValueError('Invalid native sampling-cell count')
    cells=ng**dimensions
    return max(requested//cells,2)*cells


def count_evidence(log,job):
    configuration=re.findall(r'about to integrate\s+(\d+)\s+(\d+)\s+(\d+)',log)
    if len(configuration)!=1:
        raise ValueError('Missing native integration dimensions/count')
    dimensions,requested,iterations=map(int,configuration[0])
    actual=effective_count(requested,dimensions)
    updates=[tuple(map(int,p)) for p in re.findall(r'Update # PS points:\s+(\d+)\s+-->\s+(\d+)',log)]
    expected_updates=([(requested,actual)],) if requested!=actual else \
                     ([],[(requested,requested)])
    if (requested!=job['npoints'] or actual!=job['npoints_done'] or iterations!=1
            or len(re.findall(r'------- iteration\s+\d+',log))!=1
            or updates not in expected_updates):
        raise ValueError('Native draw count does not match its deterministic pre-sampling rule')
    if updates and log.index('Update # PS points:')>log.index('Ranmar initialization seeds'):
        raise ValueError('Native point-count change was not recorded before sampling')
    return dict(requested=requested,effective=actual,dimensions=dimensions,
                logged_noop_update=updates==[(requested,requested)],
                policy='MINT initialize_even_random_numbers; selected before the first draw')


def read(record_path,output,*,variant='LO',kind='current'):
    if output.exists() or output.with_suffix('.json').exists():
        raise ValueError('Refuse to overwrite native current batches')
    record=json.loads(record_path.read_text())
    if (digest(record_path.with_suffix('.gz'))!=record['archive_sha256']
            or digest(record['audit'])!=record['audit_sha256']):
        raise ValueError('Native current raw archive or audit changed')
    final=load(record['audit'])
    if variant not in ('LO','P') or kind not in ('current','stable'):
        raise ValueError('Unknown native inclusive order or process')
    if variant=='LO' and kind=='current':
        vectors=native_vectors
        scope='Undecayed stable-top plus associated current only; no fiducial decay observables'
    else:
        from native_inclusive import native_vectors as inclusive_vectors, native_scope
        vectors=lambda histograms: inclusive_vectors(histograms,variant)
        scope=native_scope(kind)
    if final['report'].get('scope')!=scope or final['report']['variant']!=variant:
        raise ValueError('Native count policy does not match the audited process/order')
    workers=record['workers']
    directories={w['directory']:i for i,w in enumerate(workers)}
    if len(directories)!=len(workers):
        raise ValueError('Duplicate native worker directory')
    expected={w['directory']+'/'+name:sha for w in workers for name,sha in w['files'].items()}
    expected.update(record['history_files'])
    logs,history,values,errors,seeds,initial={},{},{},{},{},{}
    with tarfile.open(record_path.with_suffix('.gz'),'r|gz') as archive:
        for member in archive:
            if not member.isfile():
                raise ValueError('Unexpected native worker archive entry')
            raw=archive.extractfile(member).read()
            if hashlib.sha256(raw).hexdigest()!=expected.pop(member.name,None):
                raise ValueError('Changed or unexpected native worker file')
            if member.name in record['history_files']:
                history[member.name]=raw
                continue
            directory,name=member.name.rsplit('/',1)
            index=directories[directory]
            if name=='log.txt':
                log=raw.decode()
                logs[index]=log
                pair=re.findall(r'Ranmar initialization seeds\s+(\d+)\s+(\d+)',log)
                if len(pair)!=1 or re.findall(r'with seed\s+(\d+)',log)!=[str(record['seed'])]:
                    raise ValueError('Unverified native worker stream')
                seeds[index]=tuple(map(int,pair[0]))
                initial[index]=re.findall(r'^channel.*$',log.split('------- iteration')[0],re.M)
                if not initial[index]:
                    raise ValueError('Missing native trained-channel diagnostics')
            elif name=='MADatNLO.HwU':
                row=vectors(HwUList(io.StringIO(raw.decode()),raw_labels=True))
                if row['titles']!=final['titles'] or not np.array_equal(row['edges'],final['edges']):
                    raise ValueError('Native worker histogram layout changed')
                values[index]=row['values'][:,row['indices']]
                errors[index]=row['values'][:,row['columns'].index('dy')]
    if expected or set(logs)!=set(range(len(workers))) or set(values)!=set(logs):
        raise ValueError('Incomplete native worker archive')
    counts={i:count_evidence(logs[i],w['job']) for i,w in enumerate(workers)}
    # Reuse structural/equal-weight checks on an in-memory validation view.
    # Raw job records and their requested/effective counts are never rewritten.
    effective_workers=copy.deepcopy(workers)
    for i,w in enumerate(effective_workers):
        w['job']['npoints']=counts[i]['effective']
    groups=validate_jobs(effective_workers)
    strata=np.empty(len(workers),dtype=int)
    for label,indices in enumerate(groups.values()):
        if any(counts[i]!=counts[indices[0]] or initial[i]!=initial[indices[0]] for i in indices):
            raise ValueError('Unequal native count/proposal policy within a stratum')
        strata[indices]=label
    proof=validate_history(Path(final['report']['execution']['log']).read_text(),record['seed'],history,
                           {w['directory']:w['files']['log.txt'] for w in workers})
    array=np.asarray([values[i] for i in range(len(workers))])
    error=np.asarray([errors[i] for i in range(len(workers))])
    tolerance=3.e-6*(np.abs(array).sum(axis=0)+np.abs(final['values']))+1.e-20
    if (not np.all(np.abs(array.sum(axis=0)-final['values'])<=tolerance)
            or not np.allclose(np.sqrt((error**2).sum(axis=0)),final['errors'],rtol=3.e-6,atol=1.e-20)):
        raise ValueError('Native batch sums do not reproduce every final weight/error')
    np.savez_compressed(output,contributions=array,individual_errors=error,strata=strata,
                        random_seeds=np.asarray([seeds[i] for i in range(len(workers))]),
                        titles=final['titles'],offsets=final['offsets'],edges=final['edges'],weights=final['weights'])
    save(output.with_suffix('.json'),dict(created_utc=now(),input=str(record_path.resolve()),
        input_sha256=digest(record_path),arrays_sha256=digest(output),worker_count=len(workers),
        refinement_rng_audit=proof,native_variant=variant,native_kind=kind,
        count_policy={workers[i]['directory']:counts[i] for i in counts},
        script_sha256=digest(Path(__file__)),
        status='equal-count conditional native-current batches; deterministic effective-count audit passed',
        limitations='Count chosen before sampling, with original requested/effective counts retained. '
                    'No independent retraining or rare-tail/adaptive-stop coverage certified.'))
    print('Audited',len(workers),'native batches and all 183 canonical weights',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--variant',choices=('LO','P'),default='LO')
    parser.add_argument('--kind',choices=('stable','current'),default='current')
    args=parser.parse_args()
    read(args.record,args.output,variant=args.variant,kind=args.kind)
