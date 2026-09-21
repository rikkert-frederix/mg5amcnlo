#!/usr/bin/env python3
"""Archive inspection of completed ttW/ttbar pilots, without releasing main MC."""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from compare_current_batches import all_stage_pairs, batch_record
from pilot_report import clean


def run(output):
    if output.exists():
        raise ValueError('Refuse to overwrite a reference inspection')
    inputs,streams,counts={},set(),{}
    def record(path):
        path=Path(path).resolve()
        inputs[str(path)]=digest(path)
        return json.loads(path.read_text())
    def verify_batch(path):
        path=Path(path).resolve()
        if str(path) in counts:
            return
        metadata,worker= batch_record(path)
        audit_path=Path(worker['audit'])
        audit=record(audit_path)
        if (digest(path)!=metadata['arrays_sha256'] or
                digest(metadata['input'])!=metadata['input_sha256'] or
                digest(audit_path)!=worker['audit_sha256'] or
                digest(audit['path'])!=audit['output_sha256']):
            raise ValueError('Reference batch/archive/final-output evidence changed')
        process=Path(audit['path']).parents[2]
        archive=process/'study_cards'/Path(audit['path']).parent.name
        hashes=audit['manifest'].get('files',audit['manifest'].get('hashes',{}))
        for filename,sha in hashes.items():
            if filename.endswith('.dat') and digest(archive/filename)!=sha:
                raise ValueError('Reference input card archive changed')
        pairs=all_stage_pairs(process,worker)
        if len(pairs)!=metadata['refinement_rng_audit']['distinct_initialization_pairs'] or streams & pairs:
            raise ValueError('Missing or overlapping stage streams in reference inspection')
        streams.update(pairs)
        counts[str(path)]=len(pairs)
        inputs[str(path)]=digest(path)
        record(path.with_suffix('.json'))
        record(metadata['input'])
    comparisons={}
    for charge in ('plus','minus'):
        path=STUDY/'results'/('reference_2005_'+charge+'_ref2005_v5_joint.json')
        joint=record(path)
        original=record(joint['input'])
        if digest(joint['input'])!=joint['input_sha256'] or digest(path.with_suffix('.npz'))!=joint['arrays_sha256']:
            raise ValueError('Joint ttW reference comparison changed')
        inputs[str(path.with_suffix('.npz'))]=joint['arrays_sha256']
        for row in joint['batch_inputs']:
            if digest(row['path'])!=row['sha256']:
                raise ValueError('Joint ttW input changed')
            verify_batch(row['path'])
        for artifact in original['artifacts'].values():
            if 'analytic_virtual_checks' in artifact:
                virtual=record(artifact['analytic_virtual_checks'])
                if (virtual['status']!='all archived top/antitop analytic virtual checks passed' or
                        virtual['maximum_relative_difference']>virtual['tolerance']):
                    raise ValueError('ttW analytic validation is incomplete')
        comparisons[charge]=joint['comparisons']
    repeat=record(STUDY/'results/reference_2005_plus_P_v4_v5_retraining.json')
    for path,sha in repeat['inputs'].items():
        if digest(path)!=sha:
            raise ValueError('Independent P-repeat evidence changed')
        verify_batch(path)
    ttbar_path=STUDY/'results/ttbar_reference_ref1901_v6.json'
    ttbar=record(ttbar_path)
    ttbar_queue=record(STUDY/'inputs/ttbar_reference_queue_ref1901_v6.json')
    preflight=record(STUDY/'inputs/ttbar_reference_preflight.json')
    if ttbar_queue['status']!='ttbar reference pilot finished; inspect comparison and convergence':
        raise ValueError('ttbar reference has not completed')
    if ttbar['reference_preflight_sha256']!=digest(STUDY/'inputs/ttbar_reference_preflight.json'):
        raise ValueError('ttbar reference conventions changed')
    for variant,path in ttbar['batches'].items():
        verify_batch(path)
    for job in ttbar_queue['jobs']:
        audit=record(job['audit'])
        if audit['manifest']['physics']!=preflight['physics']:
            raise ValueError('ttbar physical inputs do not match the supplied reference')
        if job['variant']!='LO':
            virtual=record(job['analytic_virtual_checks'])
            if virtual['status']!='all archived top/antitop analytic virtual checks passed':
                raise ValueError('Missing ttbar analytic validation')
    curves={}
    for label,result in ttbar['comparisons'].items():
        error=np.asarray(result['mc_error'])
        value=np.asarray(result['value'])
        reference=np.asarray(result['reference'])
        curves[label]=dict(maximum_central_residual_over_combined_mc_error=result['max_abs_central_bin_pull'],
            nominal_relative_mc_error_range=[np.min(error[:,0]/value[:,0]),np.max(error[:,0]/value[:,0])],
            maximum_seven_scale_residual_over_our_mc_error=np.max(np.abs((value-reference)/error)))
    save(output,clean(dict(created_utc=now(),
        status='reference pilots inspected; direct-scale and absolute-normalization evidence still required',
        inputs=inputs,script_sha256=digest(Path(__file__)),
        benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),ttw=comparisons,ttbar=curves,
        plus_P_independent_repeat=repeat['nominal_fiducial'],stage_stream_counts=counts,
        total_distinct_initialization_pairs=len(streams),cross_reference_rng_overlap=0,
        assessment='Both ttW NWA central values and production-only estimates are compatible with the '
            'rounded literature values at the measured conditional precision. Both ttbar normalization '
            'conventions have central-bin residuals below 1.7 combined MC errors. The independent W+ P '
            'estimates differ by about two conditional errors; retain both and continue retraining checks.',
        limitations='ttW published MC errors are unavailable. ttbar conditional shape errors range up to '
            'about 9%; this is pilot calibration, not percent-level shape agreement. No global goodness-of-fit '
            'is inferred from correlated bins or scales. The main full-flavour campaign remains unreleased; '
            'direct scale/QES checks, absolute normalization, and the original convergence/scope ledger remain required.')))
    print(output,';',len(streams),'distinct reference stage initializations')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    run(parser.parse_args().output)
