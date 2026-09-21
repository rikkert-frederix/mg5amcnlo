#!/usr/bin/env python3
"""One fresh massive Pi pilot validating the Born-order fix end to end."""
import argparse
import json
from pathlib import Path
from campaign import ROOT, STUDY, RUNTIME_SOURCES, benchmark_param_card, digest, now, run, save
from production_invariance import run as production_identity
from check_phase_space_support import run as support
from run_sampler_pilots import run_settings
from audit_results import audit
from harvest_splits import harvest
from audit_virtuals import audit as virtuals
from read_splits import read
from check_decay_scale_cancellation import run as scale_check


def main():
    process = STUDY/'processes/TTWplus_onshell_eemu_mb4p8_both_alignment_v1'
    reference = STUDY/'processes/TTWplus_onshell_eemu_mb0p0_both_rng_v1'
    destination = STUDY/'inputs/alignment_pilot_v1.json'
    if destination.exists():
        raise ValueError('Never overwrite or silently restart this pilot')
    case = dict(w_treatment='onshell',variant='Pi',production_sampling='flat',seed=61001)
    name = 'Pi_onshell_core-w-ht-half_separate_61001_mb4p8'
    record = dict(created_utc=now(),status='preparing',process=str(process),run=name,case=case,
        source_hashes={p:digest(ROOT/p) for p in RUNTIME_SOURCES},max_cores=64,
        limitation='One massive W+ e,e,mu technical Pi pilot; full study remains incomplete.')
    save(destination,record)
    try:
        benchmark_param_card(process,json.loads((STUDY/'inputs/benchmark.json').read_text()),'onshell')
        identity = STUDY/'inputs/alignment_v1_production_identity.json'
        production_identity(reference,process,identity)
        support_path = STUDY/'inputs/alignment_v1_phase_space_support.json'
        support(support_path,process)
        record.update(status='integrating',production_identity=str(identity),support=str(support_path))
        save(destination,record)
        run(run_settings(process,case,.03))
        audit(process/'Events'/name/'MADatNLO.HwU')
        audit_path = STUDY/'results/audits'/(process.name+'_'+name+'.json')
        workers = STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
        harvest(audit_path,workers)
        vpath = STUDY/'results/alignment_v1_virtual_checks.json'
        virtuals(workers.with_suffix('.json'),vpath)
        batches = STUDY/'results/alignment_v1_batches.npz'
        read(workers.with_suffix('.json'),batches)
        scales = STUDY/'results/alignment_v1_inclusive_decay_scales.json'
        scale_check(batches,scales)
        record.update(status='completed technical pilot; inspect precision and remaining study',
            audit=str(audit_path),workers=str(workers),batches=str(batches),
            virtual_checks=str(vpath),inclusive_decay_scales=str(scales),finished_utc=now())
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__ == '__main__':
    main()
