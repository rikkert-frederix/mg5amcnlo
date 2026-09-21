#!/usr/bin/env python3
"""Collect complete, independently retrained main vectors without flavour sums.

The resulting groups feed full_flavour_statistics.py. Within-run conditional
variances are diagnostics for retraining spread; they must not be added again
to a variance estimated from independently retrained complete run vectors.
"""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from load_results import assert_matched_physics, assert_same_layout, load
from run_main_campaign import cases
from run_width_mass_pilots import verified_pairs

DONE='main integrations finished; flavour sums and convergence analysis required'


def collapse_batches(contributions,strata,final_values):
    values=np.asarray(contributions,dtype=float)
    labels=np.asarray(strata)
    final=np.asarray(final_values,dtype=float)
    if (values.ndim!=3 or values.shape[1:]!=final.shape or labels.shape!=(len(values),) or
            not np.isfinite(values).all() or not np.isfinite(final).all()):
        raise ValueError('Invalid full-run batch vectors or final histogram')
    total=values.sum(axis=0)
    tolerance=3.e-6*(np.abs(values).sum(axis=0)+np.abs(final))+1.e-20
    if not (np.abs(total-final)<=tolerance).all():
        raise ValueError('Full-run vector differs from its audited final histogram')
    variance=np.zeros_like(total)
    for label in sorted(set(labels)):
        rows=values[labels==label]
        if len(rows)<5:
            raise ValueError('Too few conditional batches in a main stratum')
        # Archived contributions already contain their 1/N stratum weight.
        variance+=len(rows)*np.var(rows,axis=0,ddof=1)
    return total,variance


def inventory(queue):
    if queue['status']!=DONE or not queue['cases']:
        raise ValueError('Require the complete main integration campaign')
    replicas=max(case['replica'] for case in queue['cases'])+1
    expected=cases(replicas,queue['cases'][0]['seed'])
    if queue['cases']!=expected or [j['case'] for j in queue['jobs']]!=expected:
        raise ValueError('Missing or reordered full-flavour/prescription/retraining inventory')
    return replicas


def collect(queue_path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite main replica vectors')
    queue=json.loads(queue_path.read_text())
    count=inventory(queue)
    if queue['benchmark_sha256']!=digest(STUDY/'inputs/benchmark.json'):
        raise ValueError('Main benchmark changed')
    groups,rows,streams={},[],set()
    reference,layout=None,None
    for job in queue['jobs']:
        case=job['case']
        result=load(job['audit'])
        report=result['report']
        manifest,execution=report['manifest'],report['execution']
        generation=execution['export_generation']['args']
        expected=dict(charge=case['charge'],flavours=case['flavours'],w_treatment='onshell',
                      decay_bottom_mass=0.,corrected='both',tag=queue['tag'])
        if (job['variant']!=case['variant'] or report['variant']!=case['variant'] or
                manifest['settings']['iseed']!=case['seed'] or execution['pilot'] or
                execution.get('grid_reference') is not None or
                manifest['production_scale']!='core-w-ht-half' or
                manifest['decay_scale_grouping']!='separate' or
                any(generation[k]!=v for k,v in expected.items())):
            raise ValueError('Main physical assignment or independent-retraining convention changed')
        if any(queue['source_hashes'].get(p)!=sha for p,sha in execution['source_hashes'].items()):
            raise ValueError('Main run did not use the frozen runtime sources')
        projected=dict(result,report=copy.deepcopy(report))
        projected['report']['execution']['export_generation']['args'].update(
            charge='plus',flavours=['e','e','e'])
        # Only the two explicitly summed/disjoint labels are projected out.
        # PDFs, all masses/widths/scales/cuts and the active decay axes match.
        if reference is not None:
            assert_matched_physics([reference,projected])
            assert_same_layout([layout,result])
        else:
            reference=projected
            layout={key:result[key] for key in ('titles','edges','offsets','weights')}
        pairs=verified_pairs(job)
        if pairs & streams:
            raise ValueError('Main retrainings/flavours share actual stage initializations')
        streams.update(pairs)
        for h,title in enumerate(result['titles']):
            if ('W-' if case['charge']=='plus' else 'W+') in title:
                a,b=map(int,result['offsets'][h:h+2])
                if np.any(result['values'][a:b]):
                    raise ValueError('Single-charge main run populated the other charge')
        with np.load(job['batches']) as batches:
            if (list(batches['titles'])!=result['titles'] or
                    not np.array_equal(batches['edges'],result['edges']) or
                    not np.array_equal(batches['offsets'],result['offsets'])):
                raise ValueError('Main batch layout differs from its final audit')
            mean,conditional_variance=collapse_batches(batches['contributions'],batches['strata'],result['values'])
        key=case['charge']+'__'+case['variant']+'__'+''.join(case['flavours'])
        group=groups.setdefault(key,dict(charge=case['charge'],variant=case['variant'],flavours=case['flavours'],
            replicas=[],vectors=[],conditional_variances=[]))
        group['replicas'].append(case['replica'])
        group['vectors'].append(mean)
        group['conditional_variances'].append(conditional_variance)
        rows.append(dict(case=case,audit=job['audit'],audit_sha256=digest(job['audit']),
            batches=job['batches'],batches_sha256=digest(job['batches']),stage_initializations=len(pairs)))
        print('Collected complete run',len(rows),'/',len(queue['jobs']),job['run'],flush=True)
    arrays=dict(layout)
    definitions={}
    for key,group in groups.items():
        if group['replicas']!=list(range(count)):
            raise ValueError('Missing or repeated independently retrained estimate in a flavour group')
        arrays[key]=np.asarray(group.pop('vectors'))
        arrays[key+'__conditional_variances']=np.asarray(group.pop('conditional_variances'))
        definitions[key]=group
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'),**arrays)
    save(output,dict(created_utc=now(),status='complete independent main vectors collected; flavour reduction and convergence required',
        input=str(queue_path.resolve()),input_sha256=digest(queue_path),
        benchmark_sha256=queue['benchmark_sha256'],arrays_sha256=digest(output.with_suffix('.npz')),
        script_sha256=digest(Path(__file__)),groups=definitions,samples=rows,
        distinct_stage_initializations=len(streams),cross_run_rng_overlap=0,
        convention='Each array row is one complete independently retrained cross-section vector. '
            'Average within a flavour group and sum all eight disjoint flavour means at identical weights. '
            'No arbitrary pairing across independent flavours, prescriptions or charges. '
            'Within-run conditional variances are diagnostic only; never add them again to retraining variance.',
        limitations='This archive supplies data for convergence/precision tests; it does not certify them. '
            'Fiducial precision, sparse-bin coverage, common final bins, scale/PDF reductions and '
            'BW/massive/parameter robustness remain required.'))
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('queue',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    collect(args.queue,args.output)
