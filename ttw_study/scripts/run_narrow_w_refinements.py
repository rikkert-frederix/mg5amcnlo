#!/usr/bin/env python3
"""Five independent strict-NLO on-shell retrainings for the narrow-W reference."""
import argparse
import json
from pathlib import Path
import time

from audit_virtuals import audit as audit_virtuals
from campaign import ROOT, STUDY, RUNTIME_SOURCES, digest, generate, now, save
from narrow_w_inputs import COMMON_SCALE
from read_splits import read
from run_narrow_w_pilots import run_case
from run_normalization_checks import DONE as PREVIOUS_DONE
from run_width_mass_pilots import verified_pairs

DONE='narrow-W strict reference retrainings finished; combine and inspect convergence'


def cases(seed_start=87001,replicas=5):
    if not isinstance(replicas,int) or replicas<5:
        raise ValueError('At least five independently retrained reference estimates are required')
    if not isinstance(seed_start,int) or not 0<seed_start<=900000000-replicas:
        raise ValueError('Invalid narrow-reference seed range')
    return [dict(variant='S',w_treatment='onshell',width_factor=1.,decay_bottom_mass=0.,
        production_sampling='flat',seed=seed_start+i,replica=i) for i in range(replicas)]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',required=True,type=Path)
    parser.add_argument('--width-inputs',required=True,type=Path)
    parser.add_argument('--accuracy',type=float,default=.01)
    parser.add_argument('--seed-start',type=int,default=87001)
    parser.add_argument('--replicas',type=int,default=5)
    args=parser.parse_args()
    if not args.tag.replace('_','').isalnum() or not 0<args.accuracy<1:
        raise ValueError('Invalid narrow-reference settings')
    args.width_inputs=args.width_inputs.resolve()
    inputs=json.loads(args.width_inputs.read_text())
    if (inputs['benchmark_sha256']!=digest(STUDY/'inputs/benchmark.json')
            or inputs['script_sha256']!=digest(STUDY/'scripts/narrow_w_inputs.py')
            or inputs['common_scale_GeV']!=COMMON_SCALE):
        raise ValueError('Common-scale narrow-W input convention changed')
    selected=cases(args.seed_start,args.replicas)
    destination=STUDY/'inputs'/('narrow_w_refinement_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Existing narrow-reference queue; do not restart implicitly')
    files=[Path(__file__).resolve(),args.width_inputs,STUDY/'inputs/benchmark.json',
        *[STUDY/'scripts'/name for name in ('run_narrow_w_pilots.py','narrow_w_inputs.py','run_normalization_checks.py',
            'run_inclusive.py','audit_results.py','audit_virtuals.py','harvest_splits.py','read_splits.py',
            'rng_history.py','run_width_mass_pilots.py','load_results.py')]]
    sources={p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    sources.update({str(p.relative_to(ROOT)):digest(p) for p in files})
    record=dict(created_utc=now(),status='waiting for absolute normalization integrations',
        source_hashes=sources,after_queue=str(args.after_queue.resolve()),width_inputs=str(args.width_inputs),
        max_cores=64,cases=selected,jobs=[],accuracy=args.accuracy,grid_points=1000,
        scope='Five independent S on-shell reference retrainings at the narrow-W common fixed scales. '
              'No grid imports. Repeated estimates will be averaged, not summed as disjoint flavours. '
              'Retain the original lower-statistics reference for an independent convergence comparison.',
        precision='Initial 1% total-rate target per replica; actual between-retraining fiducial/shape '
                  'uncertainty, not this request, determines limiting-test sensitivity.')
    save(destination,record)
    def unchanged():
        if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
            raise ValueError('Frozen narrow-reference input/source changed')
    try:
        while True:
            previous=json.loads(args.after_queue.read_text())
            if previous['status']=='stopped':
                raise ValueError('Absolute normalization predecessor stopped')
            if previous['status']==PREVIOUS_DONE:
                break
            time.sleep(15)
        unchanged()
        process=STUDY/'processes'/('TTWplus_onshell_eemu_mb0p0_both_'+args.tag)
        record.update(status='running strict narrow-reference retrainings',
                      predecessor_sha256=digest(args.after_queue),process=str(process))
        save(destination,record)
        generate(argparse.Namespace(charge='plus',w_treatment='onshell',flavours=['e','e','mu'],
                                    corrected='both',decay_bottom_mass=0.,tag=args.tag))
        streams=set()
        for case in selected:
            unchanged()
            record['current_case']=case
            save(destination,record)
            job=run_case(process,case,inputs,args.width_inputs,args.accuracy,sources,grid_points=1000)
            batches=STUDY/'results'/(process.name+'_'+job['run']+'_batches.npz')
            read(Path(job['workers']).with_suffix('.json'),batches)
            job['batches']=str(batches)
            virtuals=STUDY/'results'/(process.name+'_'+job['run']+'_virtuals.json')
            audit_virtuals(Path(job['workers']).with_suffix('.json'),virtuals)
            job['analytic_virtual_checks']=str(virtuals)
            pairs=verified_pairs(job)
            if streams & pairs:
                raise ValueError('Narrow-reference retrainings reuse actual stage streams')
            streams.update(pairs)
            record['jobs'].append(job)
            record['distinct_stage_pairs']=len(streams)
            save(destination,record)
        record.update(status=DONE,current_case=None,finished_utc=now())
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__=='__main__':
    main()
