#!/usr/bin/env python3
"""Matched ttW reference rates with independently checked common-scale batches."""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import digest, now, save
from compare_current_batches import all_stage_pairs, batch_record
from pilot_report import clean
from read_splits import split_ensembles
from reference_batches import POINTS
from replica_statistics import independent_jackknife_many
from run_literature_reference import unexpanded


def run(path, output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite matched-reference statistics')
    original = json.loads(path.read_text())
    streams, counts, ensembles, records = set(), [], [], []
    for variant in ('LO', 'P', 'S'):
        artifact = original['artifacts'][variant]
        batch = Path(artifact['batches'])
        meta, record = batch_record(batch)
        if (digest(batch) != meta['arrays_sha256'] or digest(meta['input']) != meta['input_sha256']
                or digest(record['audit']) != record['audit_sha256']):
            raise ValueError('Changed reference batch evidence')
        audit = json.loads(Path(record['audit']).read_text())
        if (audit['variant'] != variant or audit['manifest']['physics'] != original['physics']
                or digest(audit['path']) != audit['output_sha256']):
            raise ValueError('Reference physics or final output changed')
        # A verified completed LO/P variant may be retained while S gets a
        # fresh export. Its stage evidence belongs to its original process.
        pairs = all_stage_pairs(Path(audit['path']).parents[2], record)
        if pairs & streams:
            raise ValueError('Actual reference stage streams overlap between variants')
        streams.update(pairs)
        counts.append(len(pairs))
        with np.load(batch) as data:
            groups = split_ensembles(data['contributions'], data['strata'])
        ensembles.append(groups)
        records.append(dict(path=str(batch), sha256=meta['arrays_sha256']))
    boundaries = np.cumsum([0]+[len(groups) for groups in ensembles])
    def transform(*means):
        born, production, strict = [sum(means[a:b]) for a,b in zip(boundaries[:-1],boundaries[1:])]
        return np.asarray([production, unexpanded(strict, born)])
    result = independent_jackknife_many([row for groups in ensembles for row in groups], transform)
    charge = 'plus' if original['physics'] and 'TTWplus_' in original['process'] else 'minus'
    slot = 1 if charge == 'plus' else 3
    seven = [i+1 for i,(r,f) in enumerate(POINTS) if r/f in (.5,1.,2.)]
    comparisons = {}
    for i, label in enumerate(('P', 'unexpanded_NWA')):
        values = result['value'][i,slot]*1.e6
        errors = result['mc_error'][i,slot]*1.e6
        target = original['comparisons'][label]['published_rounded_ab']
        comparisons[label] = dict(value_ab=values[0], joint_mc_error_ab=errors[0],
            published_rounded_ab=target, residual_over_our_joint_mc_error=(values[0]-target)/errors[0],
            seven_point_envelope_ab=[min(values[seven]),max(values[seven])],
            common_scale_values_ab=values, common_scale_joint_mc_errors_ab=errors)
    np.savez_compressed(output.with_suffix('.npz'), **result)
    save(output, clean(dict(created_utc=now(), input=str(path), input_sha256=digest(path),
        batch_inputs=records, stage_stream_counts=counts, cross_variant_rng_overlap=0,
        comparisons=comparisons, scale_points=['nominal']+POINTS,
        arrays_sha256=digest(output.with_suffix('.npz')), script_sha256=digest(Path(__file__)),
        convention='Independent LO/P/S and beam strata; top widths fixed under common scale variation. '
                   'All rates and weight covariance are in the arrays. Seven-point envelope excludes '
                   'the two anti-correlated production corners; MC errors are separate.',
        limitations='Conditional trained-grid MC, not retraining coverage. Published MC errors '
                    'are unavailable; residual/our-error is not a reference-inclusive statistical pull.',
        status='matched ttW reference joint statistics evaluated; inspect residuals and convergence')))
    print(json.dumps(clean(comparisons), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.report, args.output)
