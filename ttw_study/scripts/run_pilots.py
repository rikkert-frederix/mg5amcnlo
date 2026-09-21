#!/usr/bin/env python3
"""Serialized six-prescription, two-charge pilot queue with immutable records.

This is a pilot campaign, not the main statistical sample. It waits for an
existing study launch, checks its output, and stops at the first failed job.
"""
import argparse
import fcntl
import json
from pathlib import Path
import sys

from campaign import ROOT, STUDY, RUNTIME_SOURCES, digest, generate, now, run, save
from audit_results import audit
from load_results import load


def compatible_reference(path,charge,variant,process):
    """Reuse only audited physics after an explicitly scheduling-only change."""
    result = load(path)
    report = result['report']
    metadata = report['manifest']
    generation = report['execution']['export_generation']['args']
    required = dict(charge=charge,flavours=['e','e','mu'],w_treatment='onshell',
                    decay_bottom_mass=0.,corrected='both')
    if any(generation[key] != value for key,value in required.items()):
        raise ValueError('Reference charge/flavour/process mismatch')
    if report['variant'] != variant or metadata['settings']['lhaid'] != 331700:
        raise ValueError('Reference prescription or PDF mismatch')
    if metadata['production_scale'] != 'core-w-ht-half' or metadata['decay_scale_grouping'] != 'separate':
        raise ValueError('Reference scale definition mismatch')
    if not metadata['width_source'].endswith(digest(STUDY/'inputs/benchmark.json')):
        raise ValueError('Reference benchmark/width source changed')
    allowed = {'madgraph/interface/amcatnlo_run_interface.py','ttw_study/scripts/campaign.py',
               'Template/fNLO/SubProcesses/ajob_template'}
    for source,checksum in report['execution']['source_hashes'].items():
        if source not in allowed and digest(ROOT/source) != checksum:
            raise ValueError('Reference physics source changed: '+source)
    hwu = Path(report['path'])
    old_process = hwu.parents[2]
    for source in ('FixedOrderAnalysis/HwU.f90',
                   'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90','SubProcesses/mint_module.f90'):
        if digest(old_process/source) != digest(process/source):
            raise ValueError('Reference generated physics differs: '+source)
    return old_process,hwu.parent.name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag', default='r3')
    parser.add_argument('--accuracy', type=float, default=.05)
    parser.add_argument('--job-seconds', type=float, default=1.,
                        help='grid-based split target; refinement PDF work is absent from grid timing')
    parser.add_argument('--charges', nargs='+', choices=('plus', 'minus'), default=['plus','minus'])
    parser.add_argument('--previous-plus-strict', type=Path,
                        help='reuse an already completed compatible W+ S HwU after a scheduling-only update')
    parser.add_argument('--plus-grid-reference',type=Path,
                        help='recover compatible W+ S training grids into this fresh export')
    parser.add_argument('--reference-queue',type=Path,
                        help='reuse completed compatible audits from an earlier stopped queue')
    args = parser.parse_args()
    record_path = STUDY/'inputs'/('onshell_pilot_queue_'+args.tag+'.json')
    if record_path.exists():
        raise ValueError('Queue already recorded; inspect its progress before resuming')
    frozen = {name: digest(ROOT/name) for name in RUNTIME_SOURCES}
    references = {}
    if args.reference_queue:
        for job in json.loads(args.reference_queue.read_text())['jobs']:
            report = json.loads(Path(job['audit']).read_text())
            key = (report['execution']['export_generation']['args']['charge'],report['variant'])
            if key in references:
                raise ValueError('Duplicate reference prescription')
            references[key] = Path(job['audit'])
    record = dict(started_utc=now(), status='waiting for existing study launch',
                  source_hashes=frozen, max_cores=64, jobs=[], current_job=None,
                  accuracy=args.accuracy,job_seconds=args.job_seconds,
                  queue_script_sha256=digest(Path(__file__)))
    if args.reference_queue:
        record['reference_queue'] = dict(path=str(args.reference_queue.resolve()),sha256=digest(args.reference_queue))
    save(record_path, record)
    try:
        with (STUDY/'mg5.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
        for charge in args.charges:
            process = STUDY/'processes'/('TTW%s_onshell_eemu_mb0p0_both_%s' % (charge,args.tag))
            if not process.exists():
                generate(argparse.Namespace(charge=charge, w_treatment='onshell',
                                            flavours=['e','e','mu'], corrected='both',
                                            decay_bottom_mass=0., tag=args.tag))
            base = 32000 if charge == 'plus' else 33000
            variants = [(variant,base+100*i+1) for i,variant in
                        enumerate(('S','Pi','LO','P','D','PiD'))]
            for variant, seed in variants:
                if frozen != {name:digest(ROOT/name) for name in RUNTIME_SOURCES}:
                    raise RuntimeError('Source changed during queue; export a new version before continuing')
                name = '%s_onshell_core-w-ht-half_separate_%d' % (variant,seed)
                job_process = process
                if (charge,variant) in references:
                    job_process,name = compatible_reference(references[charge,variant],charge,variant,process)
                if charge == 'plus' and variant == 'S' and args.previous_plus_strict:
                    previous = args.previous_plus_strict.resolve()
                    job_process, name = previous.parents[2], previous.parent.name
                    metadata = json.loads((job_process/'study_cards'/name/'manifest.json').read_text())
                    generation = json.loads((STUDY/'inputs'/(job_process.name+'_generation.json')).read_text())
                    assert generation['args']['charge'] == 'plus'
                    assert metadata['variant'] == 'S' and metadata['w_treatment'] == 'onshell'
                    assert metadata['decay_bottom_mass'] == 0.
                    assert metadata['settings']['lhaid'] == 331700
                    for source in ('FixedOrderAnalysis/HwU.f90',
                                   'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90',
                                   'SubProcesses/mint_module.f90'):
                        assert digest(job_process/source) == digest(process/source)
                record.update(status='running', current_job=dict(process=str(job_process), run=name))
                save(record_path, record)
                execution = job_process/'study_cards'/name/'execution.json'
                if execution.exists():
                    saved = json.loads(execution.read_text())
                    if saved['status'] != 'finished':
                        raise RuntimeError('Existing job did not finish: '+str(execution))
                else:
                    run(argparse.Namespace(process_dir=str(job_process), variant=variant, seed=seed,
                                           accuracy=args.accuracy, points=1000, grid_points=200,
                                           job_seconds=args.job_seconds,
                                           grid_reference=(args.plus_grid_reference if charge == 'plus' and variant == 'S' else None),
                                           iterations=3, production_scale='core-w-ht-half',
                                           ecm=13000., main_run=False))
                hwu = job_process/'Events'/name/'MADatNLO.HwU'
                audit_path = STUDY/'results/audits'/(job_process.name+'_'+name+'.json')
                if not audit_path.exists():
                    audit(hwu)
                record['jobs'].append(dict(process=str(job_process), run=name, audit=str(audit_path),
                                           output_sha256=digest(hwu), finished_utc=now()))
                save(record_path, record)
        record.update(status='pilot queue finished; main precision not certified',
                      current_job=None, finished_utc=now())
    except BaseException as error:
        record.update(status='stopped', error=repr(error), stopped_utc=now())
        save(record_path, record)
        raise
    save(record_path, record)


if __name__ == '__main__':
    main()
