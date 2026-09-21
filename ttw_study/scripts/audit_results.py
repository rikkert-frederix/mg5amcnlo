#!/usr/bin/env python3
"""Archive weight-level pilot checks; never promote a pilot to a physics result."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

import numpy as np

from campaign import ROOT, STUDY, digest, now, save

sys.path.insert(0, str(ROOT))
from madgraph.various.histograms import HwUList

RATES = ('input', 'leptons', 'mass_requirements', 'fiducial_1b', 'fiducial_2b',
         '2b_0extra', '2b_1extra', '2b_2extra', '2b_3plus_extra')


def audit(path, diagnostic=False):
    path = path.resolve()
    process = path.parents[2]
    run_name = path.parent.name
    archive = process/'study_cards'/run_name
    manifest = json.loads((archive/'manifest.json').read_text())
    execution = json.loads((archive/'execution.json').read_text())
    if execution['status'] != 'finished' or execution['outputs'][str(path)] != digest(path):
        raise ValueError('Output differs from completed execution record')
    corrected_accumulator = digest(process/'FixedOrderAnalysis/HwU.f90') == digest(
        ROOT/'Template/fNLO/FixedOrderAnalysis/HwU.f90')
    if not corrected_accumulator and not diagnostic:
        raise ValueError('Older histogram source; explicitly request diagnostic-only audit')
    histograms = HwUList(str(path), raw_labels=True)
    if len(histograms) != 210:
        raise ValueError('Expected all 210 histograms')
    columns = list(histograms[0].bins[0].wgts)
    scale_labels = [label for label in columns if 'muR=' in label]
    pdf_labels = [label for label in columns if label.startswith('PDF=')]
    variant = manifest['variant']
    expected_scale_count = 9 if variant in ('LO', 'P') else 81
    assert len(scale_labels) == expected_scale_count, (variant, len(scale_labels))
    assert sorted(int(re.search(r'PDF=(\d+)', label)[1]) for label in pdf_labels) == list(range(331700,331801))
    scale_central = [label for label in scale_labels if all(
        float(value) == 1. for value in re.findall(r'(?:muR|muF|d-?6)=\s*([\d.]+)', label))]
    assert len(scale_central) == 1
    central_i, error_i = columns.index('central value'), columns.index('dy')
    physical_columns = [i for i, label in enumerate(columns) if label == 'central value'
                        or label in scale_labels or label in pdf_labels]
    rates, arrays, titles, edges, offsets = {}, [], [], [], [0]
    nominal_mismatch = partition_mismatch = 0.
    strict_extra_max = 0.
    for hist in histograms:
        assert all(list(row.wgts) == columns for row in hist.bins)
        data = np.asarray([[row.wgts[label] for label in columns] for row in hist.bins])
        assert np.isfinite(data).all()
        arrays.append(data)
        titles.append(hist.title)
        edges.extend([row.boundaries for row in hist.bins])
        offsets.append(offsets[-1]+len(hist.bins))
        scale_nominal = data[:, columns.index(scale_central[0])]
        nominal = data[:, central_i]
        delta = np.abs(scale_nominal-nominal)
        tolerance = 2.e-6*np.maximum(np.abs(nominal), np.abs(scale_nominal))+1.e-20
        assert (delta <= tolerance).all(), hist.title
        nominal_mismatch = max(nominal_mismatch, float(delta.max()))
        if ' rates:' not in hist.title:
            continue
        assert len(data) == 9
        total = data[4, physical_columns]
        partition = data[5:9, physical_columns].sum(axis=0)
        delta = np.abs(total-partition)
        tolerance = 3.e-6*(np.abs(total)+np.abs(data[5:9, physical_columns]).sum(axis=0))+1.e-20
        if corrected_accumulator:
            assert (delta <= tolerance).all(), ('Exclusive jet partition failed', hist.title)
        partition_mismatch = max(partition_mismatch, float(delta.max()))
        if variant in ('LO','P','D','S'):
            strict_extra_max = max(strict_extra_max, float(np.abs(data[7:9,physical_columns]).max()))
        rates[hist.title] = {name: dict(value_pb=float(row[central_i]),
                                        mc_error_pb=float(row[error_i]),
                                        relative_mc_error=(float(abs(row[error_i]/row[central_i]))
                                                           if row[central_i] else None),
                                        weights={columns[i]: float(row[i]) for i in physical_columns})
                            for name, row in zip(RATES, data)}
    if variant in ('LO','P','D','S'):
        assert strict_extra_max == 0.
    target = STUDY/'results/audits'/(process.name+'_'+run_name)
    if target.with_suffix('.json').exists() or target.with_suffix('.npz').exists():
        raise ValueError('Audit already archived: '+str(target))
    target.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(str(target)+'.npz', values=np.concatenate(arrays),
                        edges=np.asarray(edges), offsets=np.asarray(offsets),
                        titles=np.asarray(titles), columns=np.asarray(columns))
    record = dict(created_utc=now(), path=str(path), output_sha256=digest(path),
                  status='diagnostic only' if diagnostic else 'technical checks passed; precision not certified',
                  variant=variant, manifest=manifest, execution=execution,
                  histogram_count=len(histograms), scale_count=len(scale_labels), pdf_count=len(pdf_labels),
                  corrected_accumulator=corrected_accumulator,
                  maximum_nominal_weight_mismatch_pb=nominal_mismatch,
                  maximum_jet_partition_mismatch_pb=partition_mismatch,
                  strict_2plus_extra_maximum_pb=strict_extra_max,
                  rates=rates, arrays_sha256=digest(str(target)+'.npz'))
    exclusions=sorted((path.parent/'split_outliers').glob('step_*/audit.json'))
    record['split_outlier_audits']={str(p):digest(p) for p in exclusions}
    if exclusions:
        record['split_outlier_limitation']='Data-dependent split exclusion occurred; retained-sample errors omit selection bias.'
    save(target.with_suffix('.json'), record)
    print(str(target.with_suffix('.json')))
    for title, bins in rates.items():
        if title.startswith('R04_b25'):
            print(title)
            for name in RATES:
                row = bins[name]
                print('  %-20s %.8e +/- %.3e pb' % (name, row['value_pb'], row['mc_error_pb']))
    print('Jet-partition maximum absolute mismatch: %.3e pb' % partition_mismatch)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', type=Path, nargs='+')
    parser.add_argument('--diagnostic', action='store_true')
    args = parser.parse_args()
    for path in args.paths:
        audit(path, args.diagnostic)


if __name__ == '__main__':
    main()
