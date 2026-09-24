#!/usr/bin/env python3
"""Fresh independent Pi retrainings at QES/2 and 2*QES after pilot queues."""
import argparse
import json
from pathlib import Path

from audit_virtuals import audit as audit_virtuals
from campaign import ROOT, STUDY, RUNTIME_SOURCES, digest, generate, now, save
from production_invariance import compare_subprocess, subprocess_inventory
from read_splits import read as read_batches
from run_scale_validation import cases as technical_cases, compare, reference_controls, run_case
from run_width_mass_pilots import verified_pairs

TAG='qes_pi_retrain_v1'
SEEDS=(93001,93002)
DONE='independent Pi QES retraining pair finished; inspect cancellation and convergence'
DEPENDENCIES={
    'scale_validation_queue_tech_v4.json':'scale validation pilots finished; inspect identities and convergence',
    'normalization_queue_absolute_v5.json':'absolute normalization integrations finished; inspect joint comparisons and convergence',
    'small_mass_queue_smallmass_v5.json':'generated small-mass pilots finished; inspect continuity and statistical sensitivity',
    'rng_pilot_queue_alignment_controls_v1.json':'stage-audited S/Pi pilots finished; inspect physics and convergence',
}


def allocation():
    selected=[case for case in technical_cases() if case['name'] in ('qes_half','qes_two')]
    if ([case['name'] for case in selected]!=['qes_half','qes_two'] or
            [case['qes'] for case in selected]!=[.5,2.] or len(set(SEEDS))!=2):
        raise ValueError('Unexpected QES retraining allocation')
    return [dict(case=case,variant='Pi',seed=seed) for case,seed in zip(selected,SEEDS)]


def source_identity(reference, candidate, destination):
    if destination.exists():
        raise ValueError('Existing QES source-transition audit')
    old,new=subprocess_inventory(reference),subprocess_inventory(candidate)
    if old.keys()!=new.keys() or len(old)!=2:
        raise ValueError('QES production subprocesses changed')
    proof={str(key):dict(reference_directory=str(old[key]),candidate_directory=str(new[key]),
        **compare_subprocess(old[key],new[key],allow_kinematic_update=True)) for key in sorted(old)}
    for row in proof.values():
        if not row.get('kinematic_matrix_transition'):
            raise ValueError('Expected an explicit generic loop-helper source transition')
    save(destination,dict(created_utc=now(),status='production amplitudes and FKS ownership preserved with explicit generic loop-helper transition',
        reference=str(reference),candidate=str(candidate),subprocesses=proof,
        limitation='The fresh analytic decay-reference retry changes validation code but leaves the analytic physics provider unchanged. This source audit does not replace numerical convergence checks.'))


def run():
    destination=STUDY/'inputs'/('qes_pi_retraining_queue_'+TAG+'.json')
    process=STUDY/'processes'/('TTWplus_onshell_eemu_mb0p0_both_'+TAG)
    if destination.exists() or process.exists():
        raise ValueError('QES retraining queue or export already exists; inspect before further action')
    dependencies={}
    for name,status in DEPENDENCIES.items():
        path=STUDY/'inputs'/name
        record=json.loads(path.read_text())
        if record['status']!=status:
            raise ValueError('Incomplete QES predecessor: '+name)
        dependencies[str(path)]=digest(path)
    benchmark=STUDY/'inputs/benchmark.json'
    earlier=STUDY/'inputs/qes_pi_retraining_plan_v1.json'
    if json.loads(earlier.read_text())['new_integrations_started']!=0:
        raise ValueError('Earlier QES plan already launched integrations')
    dependencies[str(earlier)]=digest(earlier)
    dependencies[str(benchmark)]=digest(benchmark)
    source_names=list(RUNTIME_SOURCES)+['ttw_study/scripts/'+name for name in (
        'run_qes_pi_retraining.py','run_scale_validation.py','run_inclusive.py',
        'read_splits.py','audit_results.py','harvest_splits.py','audit_virtuals.py',
        'run_width_mass_pilots.py','production_invariance.py','joint_comparison.py')]
    sources={name:digest(ROOT/name) for name in source_names}
    selected=allocation()
    record=dict(created_utc=now(),status='preparing independent Pi QES retrainings',
        max_cores=64,tag=TAG,source_hashes=sources,dependencies=dependencies,
        cases=selected,jobs=[],current_case=None,process=str(process),accuracy=.01,
        predecessor_plan=str(earlier),
        seed_note='Seeds 71401/71402 from the older unlaunched plan were subsequently used by central-scale runs; 93001/93002 are fresh.',
        scope='Two fresh on-shell W+ e,e,mu Pi grids at QES/2 and 2*QES, with all 81 scale and 101 PDF weights. '
              'Compare both with the same archived central Pi control and inspect their direct difference.',
        limitations='Conditional learned-grid pilot; no QES cancellation or full-flavour convergence conclusion is automatic.')
    save(destination,record)
    def unchanged():
        if (any(digest(ROOT/name)!=sha for name,sha in sources.items()) or
                any(digest(path)!=sha for path,sha in dependencies.items())):
            raise ValueError('Frozen QES source or predecessor changed')
    try:
        unchanged()
        controls_path=STUDY/'inputs/rng_pilot_queue_alignment_controls_v1.json'
        references,streams=reference_controls(controls_path)
        old_technical=json.loads((STUDY/'inputs/scale_validation_queue_tech_v4.json').read_text())
        old_seeds={int(json.loads(Path(job['audit']).read_text())['manifest']['settings']['iseed'])
                   for job in old_technical['jobs']}
        old_seeds.update(int(json.loads(Path(job['audit']).read_text())['manifest']['settings']['iseed'])
                         for job in json.loads(controls_path.read_text())['jobs'])
        for job in old_technical['jobs']:
            pairs=verified_pairs(job)
            if streams & pairs:
                raise ValueError('Original QES and reference streams overlap')
            streams.update(pairs)
        if any(case['seed'] in old_seeds for case in selected):
            raise ValueError('New QES seed already used by a technical pilot')
        record['earlier_distinct_stage_pairs']=len(streams)
        record['reference_batches']=str(references['onshell','Pi'])
        unchanged()
        record.update(status='generating fresh Pi QES export',current_case=selected[0])
        save(destination,record)
        generate(argparse.Namespace(charge='plus',flavours=['e','e','mu'],
            w_treatment='onshell',decay_bottom_mass=0.,corrected='both',tag=TAG))
        reference_process=STUDY/'processes/TTWplus_onshell_eemu_mb0p0_both_alignment_controls_v1'
        identity=STUDY/'inputs'/('qes_pi_retraining_source_identity_'+TAG+'.json')
        source_identity(reference_process,process,identity)
        record['source_identity']=str(identity)
        record['source_identity_sha256']=digest(identity)
        for row in selected:
            unchanged()
            case,seed=row['case'],row['seed']
            record.update(status='running independent Pi QES retraining',current_case=row)
            save(destination,record)
            job=run_case(process,case,'Pi',seed,.01,sources,None)
            job.update(case=case,variant='Pi',seed=seed)
            batches=STUDY/'results'/('%s_%s_Pi_batches.npz'%(TAG,case['name']))
            read_batches(Path(job['workers']).with_suffix('.json'),batches)
            job['batches']=str(batches)
            pairs=verified_pairs(job)
            if streams & pairs:
                raise ValueError('New QES stage streams overlap earlier study controls')
            streams.update(pairs)
            virtual=STUDY/'results'/('%s_%s_Pi_virtuals.json'%(TAG,case['name']))
            audit_virtuals(Path(job['workers']).with_suffix('.json'),virtual)
            job['analytic_virtual_checks']=str(virtual)
            comparison=STUDY/'results'/('%s_%s_Pi_vs_central.json'%(TAG,case['name']))
            compare(batches,references['onshell','Pi'],case,comparison)
            job['comparison']=str(comparison)
            job['finished_utc']=now()
            record['jobs'].append(job)
            record['distinct_stage_pairs_including_predecessors']=len(streams)
            save(destination,record)
        record.update(status=DONE,current_case=None,finished_utc=now())
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__=='__main__':
    run()
