#!/usr/bin/env python3
"""Summarize a completed small-mass pilot without treating correlated pulls as trials."""
import argparse
import json
import math
from pathlib import Path

from campaign import digest, now, save


def inspect(queue_path, report_path, output):
    if output.exists():
        raise ValueError('Refuse to overwrite a small-mass inspection')
    queue = json.loads(queue_path.read_text())
    report = json.loads(report_path.read_text())
    if (len(queue['jobs']) != 12 or not queue['status'].startswith('generated small-mass pilots finished')
            or report['inputs'][str(queue_path.resolve())] != digest(queue_path)
            or report['arrays_sha256'] != digest(report_path.with_suffix('.npz'))
            or report['cross_sample_rng_overlap'] != 0
            or report['distinct_initialization_pairs'] != queue['distinct_stage_pairs']):
        raise ValueError('Incomplete or changed small-mass result evidence')
    modes = {}
    for mode, row in report['modes'].items():
        rates, shapes = [], []
        for title, labels in row['rates'].items():
            for rate, entries in labels.items():
                for contrast in ('S_mass_over_zero', 'Pi_mass_over_zero'):
                    for index, mass in enumerate((1., .1)):
                        value = entries[contrast]['value'][index]
                        error = entries[contrast]['mc_error'][index]
                        if value is None or error is None or error <= 0:
                            continue
                        if not math.isfinite(value) or not math.isfinite(error):
                            raise ValueError('Nonfinite small-mass rate comparison')
                        rates.append(dict(abs_pull=abs((value-1)/error), title=title,
                            rate=rate, contrast=contrast, mass_GeV=mass,
                            ratio=value, conditional_mc_error=error))
        for local, spectrum in row['spectra'].items():
            for index, mass in enumerate((1., .1)):
                for variant_index, variant in enumerate(('S', 'Pi')):
                    pulls = spectrum['nominal_mass_ratios_residuals_over_error'][index][variant_index]
                    for bin_index, pull in enumerate(pulls):
                        if pull is None:
                            continue
                        if not math.isfinite(pull):
                            raise ValueError('Nonfinite small-mass shape comparison')
                        shapes.append(dict(abs_pull=abs(pull), signed_pull=pull,
                            local=local, title=spectrum['title'], mass_GeV=mass,
                            variant=variant, bin_index_zero_based=bin_index))
        if not rates or not shapes:
            raise ValueError('Missing small-mass contrasts')
        nominal = row['rates']['R04_b25 W+ rates: input/leptons/masses/1b/2b/2b_0j']
        fiducial = {
            name: {contrast:dict(ratios=nominal[name][contrast]['value'],
                conditional_mc_errors=nominal[name][contrast]['mc_error'])
                for contrast in ('S_mass_over_zero', 'Pi_mass_over_zero',
                                 'mass_change_in_Pi_over_S')}
            for name in ('fiducial_1b', 'fiducial_2b')}
        modes[mode] = dict(stage_stream_counts=row['stage_stream_counts'],
            rate_entries=len(rates), maximum_nominal_rate_pull=max(rates,key=lambda x:x['abs_pull']),
            rate_entries_above_three=sum(x['abs_pull']>3 for x in rates),
            normalized_shape_entries=len(shapes),
            maximum_nominal_normalized_shape_pull=max(shapes,key=lambda x:x['abs_pull']),
            normalized_shape_entries_above_three=sum(x['abs_pull']>3 for x in shapes),
            nominal_fiducial_ratios=fiducial,
            production_identity_inputs=row['production_identity_inputs'])
    exclusions = []
    for job in queue['jobs']:
        audit = json.loads(Path(job['audit']).read_text())
        for path, sha in audit['split_outlier_audits'].items():
            if digest(path) != sha:
                raise ValueError('Split-exclusion evidence changed')
            detail = json.loads(Path(path).read_text())
            exclusions.append(dict(case=job['case'], audit=path, sha256=sha,
                groups=[{key:group[key] for key in ('channel','excluded_split',
                    'original_count','variance_fraction','deviation','retained_scale')}
                    for group in detail['groups']]))
    save(output,dict(created_utc=now(),
        status='pilot small-mass continuity compatible at conditional precision; further precision untested',
        inputs={str(p.resolve()):digest(p) for p in
                (queue_path,report_path,report_path.with_suffix('.npz'))},
        cases=12,distinct_initialization_pairs=report['distinct_initialization_pairs'],
        cross_sample_rng_overlap=0,modes=modes,split_exclusions=exclusions,
        interpretation='The finite-mass points bound current sensitivity but do not prove the mathematical limit. '
            'Rates and bins share events; individual residuals are not independent global-significance trials.',
        limitations='One ordered W+ flavour assignment and trained-grid conditional errors. '
            'The opt-in data-dependent split exclusion omits selection bias. '
            'These pilots do not certify percent-level shapes or full-flavour robustness.',
        script_sha256=digest(Path(__file__))))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('queue',type=Path)
    parser.add_argument('report',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    inspect(args.queue,args.report,args.output)
