#!/usr/bin/env python3
"""Inspect completed direct-scale identities with immutable provenance.

Partial inspections take an immutable snapshot of the live queue. Central
production-scale changes are physical diagnostics, not equality tests; this
report audits their completed samples but does not certify their agreement.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from campaign import ROOT, STUDY, digest, now, save
from compare_current_batches import all_stage_pairs, batch_record
from pilot_report import RATES, clean
from replica_statistics import ratio
from run_scale_validation import DONE, cases, overlapping_points


def residual_summary(values, errors):
    pulls=ratio(np.asarray(values),np.asarray(errors))
    finite=np.isfinite(pulls)
    return dict(finite_comparisons=int(finite.sum()),
                maximum_absolute_residual=(float(np.max(np.abs(pulls[finite]))) if finite.any() else None),
                comparisons_above_three=int(np.sum(np.abs(pulls[finite])>3.)))


def summarize_identity(comparison, arrays):
    points=comparison['reference_points']
    expected=[row[2] for row in overlapping_points(comparison['case']['factors'])]
    if points!=expected:
        raise ValueError('Direct comparison has incorrect physical scale coordinates')
    values,errors=arrays['rate_values'],arrays['rate_mc_errors']
    if values.shape!=(4,9,len(points)) or errors.shape!=values.shape:
        raise ValueError('Incomplete rate or matching-coordinate arrays')
    if (not np.isfinite(values[:3]).all() or not np.isfinite(errors[:3]).all()
            or np.any(errors[np.isfinite(errors)]<0.)):
        raise ValueError('Invalid absolute rate comparison or MC uncertainty')
    rates={}
    for i,label in enumerate(RATES):
        supplied=comparison['rates'][label]
        for name,expected_row in (('direct_pb',values[0,i]),('reference_pb',values[1,i]),
                ('difference_pb',values[2,i]),('difference_mc_error_pb',errors[2,i])):
            if not np.array_equal(np.asarray(supplied[name]),expected_row):
                raise ValueError('Rate summary differs from its weight archive')
        rates[label]=dict(direct_pb=values[0,i,0],reweighted_pb=values[1,i,0],
            difference_pb=values[2,i,0],difference_mc_error_pb=errors[2,i,0],
            ratio=values[3,i,0],ratio_mc_error=errors[3,i,0],
            residual_over_mc_error=ratio(values[2,i,0],errors[2,i,0]))
    spectra={}
    for local,title in zip((3,4,6,11,12,20),comparison['spectra']):
        value,error=arrays['normalized_%d'%local],arrays['normalized_mc_errors_%d'%local]
        if value.shape!=error.shape or value.shape[0]!=4 or value.shape[-1]!=len(points):
            raise ValueError('Incomplete normalized-spectrum comparison')
        nominal=residual_summary(value[2,:,0],error[2,:,0])
        supplied=comparison['spectra'][title]
        if (nominal['finite_comparisons']!=supplied['finite_nominal_pulls'] or
                nominal['maximum_absolute_residual']!=supplied['maximum_absolute_nominal_pull'] or
                nominal['comparisons_above_three']!=supplied['nominal_bins_above_three']):
            raise ValueError('Spectrum summary differs from its weight archive')
        spectra[title]=dict(nominal=nominal,all_matching_weights=residual_summary(value[2],error[2]),
                            normalized_to=supplied['normalized_to'])
    return clean(dict(case=comparison['case'],weight_comparison_entries=len(points),
        distinct_physical_scale_coordinates=len(set(map(tuple,points))),
        nominal_rates=rates,all_rate_weights=residual_summary(values[2],errors[2]),
        normalized_spectra=spectra))


def frozen_source_drift(source_hashes, allow_source_drift=False):
    """Return current hashes that differ from a completed queue's freeze.

    An inspection normally proves that its current runtime still matches the
    one used to generate the queue.  Historical inspection after a later,
    separately documented source repair is useful too, but it must be an
    explicit choice and the divergence must become part of the result.
    """
    drift = {path: dict(frozen_sha256=sha, current_sha256=digest(ROOT/path))
             for path, sha in source_hashes.items()
             if digest(ROOT/path) != sha}
    if drift and not allow_source_drift:
        raise ValueError('Frozen technical source or benchmark changed: '+next(iter(drift)))
    return drift


def run(queue_path, output, allow_partial=False, allow_source_drift=False):
    snapshot_path=output.with_suffix('.queue.json')
    if output.exists() or snapshot_path.exists():
        raise ValueError('Refuse to overwrite technical inspection evidence')
    queue_bytes=queue_path.read_bytes()
    queue=json.loads(queue_bytes)
    if queue['cases']!=cases():
        raise ValueError('Unexpected technical case inventory')
    complete=queue['status']==DONE
    if not complete and not allow_partial:
        raise ValueError('Require a complete technical queue, or an explicit partial inspection')
    source_drift = frozen_source_drift(queue['source_hashes'], allow_source_drift)
    inputs,streams,counts={},{},{}

    def record(path):
        path=Path(path).resolve()
        inputs[str(path)]=digest(path)
        return json.loads(path.read_text())

    def verify_batch(path):
        path=Path(path).resolve()
        if str(path) in counts:
            return streams[str(path)]
        metadata,worker=batch_record(path)
        audit=record(worker['audit'])
        if (digest(path)!=metadata['arrays_sha256'] or
                digest(metadata['input'])!=metadata['input_sha256'] or
                digest(worker['audit'])!=worker['audit_sha256'] or
                digest(audit['path'])!=audit['output_sha256'] or
                digest(Path(worker['audit']).with_suffix('.npz'))!=audit['arrays_sha256']):
            raise ValueError('Technical batch/archive/final-output evidence changed')
        process=Path(audit['path']).parents[2]
        archive=process/'study_cards'/Path(audit['path']).parent.name
        for filename,sha in audit['manifest']['hashes'].items():
            if filename.endswith('.dat') and digest(archive/filename)!=sha:
                raise ValueError('Technical input card archive changed')
        hook=audit['manifest'].get('custom_user_hook')
        if hook and digest(archive/Path(hook['source']).name)!=hook['sha256']:
            raise ValueError('Archived direct-scale hook changed')
        pairs=all_stage_pairs(process,worker)
        if len(pairs)!=metadata['refinement_rng_audit']['distinct_initialization_pairs']:
            raise ValueError('Missing technical stage streams')
        if any(pairs & previous for previous in streams.values()):
            raise ValueError('Different technical samples share actual random streams')
        streams[str(path)]=pairs
        counts[str(path)]=len(pairs)
        inputs[str(path)]=digest(path)
        record(path.with_suffix('.json'))
        record(metadata['input'])
        return pairs

    controls_path=Path(queue['validation_queue'])
    controls=record(controls_path)
    for path,sha in queue['reference_queues_sha256'].items():
        if digest(path)!=sha:
            raise ValueError('Technical predecessor/reference queue changed')
        record(path)
    references={}
    for job in controls['jobs']:
        audit=record(job['audit'])
        key=audit['manifest']['w_treatment'],job['variant']
        if key in references:
            raise ValueError('Duplicated reference control')
        references[key]=Path(job['batches']).resolve()
        verify_batch(job['batches'])
    if set(references)!={(mode,variant) for mode in ('onshell','all-bw') for variant in ('S','Pi')}:
        raise ValueError('Incomplete reference controls')
    identities,diagnostics,pending={}, {}, []
    seen=set()
    for job in queue['jobs']:
        audit=record(job['audit'])
        case=audit['manifest']['technical_test']
        key=case['name'],job['variant']
        if key in seen or case not in cases() or job['variant'] not in ('S','Pi'):
            raise ValueError('Unexpected or duplicated technical job')
        seen.add(key)
        identity=not case['name'].startswith('central_')
        required={'batches','analytic_virtual_checks'} | ({'comparison'} if identity else set())
        if not required.issubset(job):
            pending.append(list(key))
            continue
        direct_path=Path(job['batches']).resolve()
        pairs=verify_batch(direct_path)
        virtual=record(job['analytic_virtual_checks'])
        metadata=json.loads(direct_path.with_suffix('.json').read_text())
        if (virtual['status']!='all archived top/antitop analytic virtual checks passed' or
                virtual['checked_points']<=0 or virtual['maximum_relative_difference']>1.e-8 or
                virtual['tolerance']!=1.e-8 or
                Path(virtual['input']).resolve()!=Path(metadata['input']).resolve() or
                virtual['input_sha256']!=metadata['input_sha256']):
            raise ValueError('Technical analytic virtual validation is incomplete')
        validation={field:virtual[field] for field in ('checked_points','maximum_relative_difference','tolerance','status')}
        if identity:
            path=Path(job['comparison'])
            comparison=record(path)
            reference_path=references[case['mode'],job['variant']]
            if (comparison['case']!=case or set(comparison['inputs'])!={str(direct_path),str(reference_path)}
                    or any(digest(p)!=sha for p,sha in comparison['inputs'].items())
                    or digest(path.with_suffix('.npz'))!=comparison['arrays_sha256']):
                raise ValueError('Technical comparison refers to changed or mismatched evidence')
            inputs[str(path.with_suffix('.npz').resolve())]=comparison['arrays_sha256']
            with np.load(path.with_suffix('.npz')) as arrays:
                summary=summarize_identity(comparison,arrays)
            summary.update(variant=job['variant'],analytic_virtual_checks=validation,
                stage_stream_counts=[len(pairs),len(streams[str(reference_path)])],
                cross_sample_stage_overlap=0)
            identities['__'.join(key)]=summary
        else:
            diagnostics['__'.join(key)]=dict(case=case,variant=job['variant'],batches=str(direct_path),
                stage_initializations=len(pairs),analytic_virtual_checks=validation,
                interpretation='A changed central production scale is a physical comparison, not an equality test.')
    expected={(case['name'],variant) for case in cases() for variant in ('S','Pi')}
    pending.extend([list(key) for key in sorted(expected-seen)])
    if complete and (pending or len(identities)!=8 or len(diagnostics)!=4):
        raise ValueError('Complete queue has incomplete technical evidence')
    if not identities:
        raise ValueError('No complete direct-scale identity is available')
    output.parent.mkdir(parents=True,exist_ok=True)
    snapshot_path.write_bytes(queue_bytes)
    inputs[str(snapshot_path.resolve())]=digest(snapshot_path)
    result=dict(created_utc=now(),
        status=('complete technical queue inspected; assess residuals and central-scale diagnostics' if complete
                else 'partial technical inspection; remaining cases required'),
        original_queue=str(queue_path.resolve()),queue_status_at_snapshot=queue['status'],
        queue_snapshot_sha256=hashlib.sha256(queue_bytes).hexdigest(),
        inputs=inputs,source_hashes={str(Path(__file__).resolve()):digest(Path(__file__)),
            **{str(STUDY/'scripts'/name):digest(STUDY/'scripts'/name) for name in
               ('compare_current_batches.py','pilot_report.py','rng_history.py')},
            **{str((ROOT/path).resolve()):sha for path,sha in queue['source_hashes'].items()}},
        benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),
        current_source_drift=source_drift,
        identities=identities,central_scale_samples=diagnostics,pending_cases=pending,
        stage_stream_counts=counts,total_distinct_initialization_pairs=sum(counts.values()),
        cross_sample_stage_overlap=0,
        convention='Matching direct/reweighted coordinates are compared before any scale envelope. '
            'The nominal entry duplicates one physical grid coordinate. Errors are conditional on '
            'trained proposals, with independent run/beam-stratum deletion and within-run correlations.',
        limitations='Correlated bins and scale points are not independent goodness-of-fit trials. '
            'Structural zero sectors have undefined pulls. These single-assignment pilots do not '
            'establish independent-retraining coverage or final fiducial/tail precision. '
            'This inspection alone does not release the main campaign.')
    save(output,clean(result))
    print(output,';',len(identities),'identities inspected;',len(pending),'cases pending',flush=True)
    return clean(result)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queue',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--allow-partial',action='store_true')
    parser.add_argument('--allow-source-drift',action='store_true',
                        help='inspect immutable completed evidence after a recorded later source change')
    args=parser.parse_args()
    run(args.queue,args.output,args.allow_partial,args.allow_source_drift)
