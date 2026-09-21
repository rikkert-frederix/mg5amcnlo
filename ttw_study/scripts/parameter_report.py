#!/usr/bin/env python3
"""Rates and selected spectra from independently retrained parameter scans.

All scale comparisons are point matched. Shared baselines and cut-bank
correlations are preserved in joint nominal covariance factors. This report
does not certify precision or turn a PDF-family spread into a PDF error.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from campaign import ROOT, STUDY, digest, now, save
from covariance_storage import nominal_copy
from joint_report import derived_rates, normalization_bins
from load_results import POINTS, WEIGHTS
from main_binning import apply_plan, layout_signature
import parameter_charge_statistics as charge_statistics
from parameter_replicas import read as read_vectors
from parameter_statistics import LABELS, comparison_keys, estimate, select
from parameter_variation_results import SCALE_WEIGHTS
from pilot_report import RATES, clean, scales
from replica_statistics import ratio

CONFIGS = ('R04_b25', 'R04_b30', 'R04_b40', 'R03_b25', 'R05_b25')
DERIVED = ('acceptance_1b', 'acceptance_2b', 'two_b_given_one_b',
           'veto_0extra', 'fraction_1extra', 'fraction_2extra', 'fraction_3plus_extra')
DEFAULT_SPECTRA = (3, 4, 6, 11, 12, 20)


def apply_common_binning(layout, ensembles, conditional, plan):
    if layout['weights'] != SCALE_WEIGHTS:
        raise ValueError('Require the complete canonical parameter scale grid')
    compatible = dict(layout, weights=WEIGHTS)
    if plan.get('layout_sha256') == layout_signature(compatible):
        adjusted = dict(plan, layout_sha256=layout_signature(layout))
        origin = 'main 183-weight layout; identical bins/titles and common 82 coordinates verified'
    elif plan.get('layout_sha256') == layout_signature(layout):
        adjusted, origin = plan, 'exact 82-weight parameter layout'
    else:
        raise ValueError('Common main/scan bin plan differs in geometry or weight coordinates')
    updated, values, diagnostics, metadata = apply_plan(layout, ensembles, conditional, adjusted)
    metadata['plan_layout_compatibility'] = origin
    return updated, values, diagnostics, metadata


def rate_banks(layout):
    banks = {}
    for h, title in enumerate(layout['titles']):
        if ' rates:' not in title:
            continue
        config = title.split()[0]
        charge = 'plus' if 'W+' in title else 'minus' if 'W-' in title else None
        a, b = map(int, layout['offsets'][h:h+2])
        if (config, charge) in banks or b-a != 9:
            raise ValueError('Duplicated or incomplete parameter rate bank')
        banks[config, charge] = h
    if set(banks) != {(config, charge) for config in CONFIGS for charge in ('plus', 'minus')}:
        raise ValueError('Require all five cut configurations for both charges')
    return banks


def features(ensembles, layout, banks, configs, local=None, nominal=False):
    result = {}
    for key, values in ensembles.items():
        rows = []
        for config in configs:
            h = banks[config, key[3]]
            rate_start = int(layout['offsets'][h])
            if local is None:
                row = values[:, rate_start:rate_start+9]
            else:
                a, b = map(int, layout['offsets'][h+local-1:h+local+1])
                parents = [rate_start+i for i in normalization_bins(local)]
                row = np.concatenate([values[:, a:b], values[:, parents].sum(axis=1)[:, None, :]], axis=1)
            rows.append(row[:, :, :1] if nominal else row)
        result[key] = np.concatenate(rows, axis=1)
    return result


def one_comparison(result):
    return {key: value[:, 0] if key == 'covariance_factor' else value[0]
            for key, value in result.items() if key != 'ensemble_order'}


def describe(value, error, bias):
    if value.shape != (82,):
        raise ValueError('Require nominal and all 81 scale values')
    points = {'central': clean(value[0])}
    points.update({point: clean(value[i+1]) for i, point in enumerate(POINTS)})
    return dict(scale=scales.describe(points), nominal_retraining_mc_error=clean(error[0]),
                nominal_jackknife_bias_estimate=clean(bias[0]))


def summaries(result, bin_index):
    rows = {}
    for i, label in enumerate(LABELS):
        value, error, bias = (result[key][i, bin_index] for key in
                              ('value', 'mc_error', 'nonlinear_bias_estimate'))
        row = describe(value, error, bias)
        if label in LABELS[:4]:
            row['nominal_relative_mc_error'] = clean(ratio(error[0], value[0]))
        if label in ('product_shift_change', 'product_ratio_change'):
            row['nominal_change_over_mc_error'] = clean(ratio(value[0], error[0]))
            row['change_exceeds_three_mc_errors'] = bool(np.isfinite(error[0]) and error[0] > 0.
                                                         and abs(value[0]) > 3.*error[0])
        rows[label] = row
    return rows


def normalize(row):
    return ratio(row[:-1], row[-1])


def charge_summaries(result, bin_index, quantity_labels):
    rows = {}
    for q, quantity in enumerate(quantity_labels):
        rows[quantity] = {}
        for i, label in enumerate(charge_statistics.LABELS):
            value, error, bias = (result[key][q, i, bin_index] for key in
                                  ('value', 'mc_error', 'nonlinear_bias_estimate'))
            summary = describe(value, error, bias)
            if label in charge_statistics.LABELS[4:]:
                summary.update(nominal_change_over_mc_error=clean(ratio(value[0], error[0])),
                    change_exceeds_three_mc_errors=bool(np.isfinite(error[0]) and error[0] > 0.
                                                        and abs(value[0]) > 3.*error[0]))
            rows[quantity][label] = summary
    return rows


def precision(result):
    report = {}
    for i, label in enumerate(LABELS[:4]):
        value, error = result['value'][i, :, 0], result['mc_error'][i, :, 0]
        selected = np.isfinite(value)&np.isfinite(error)&(value > 0.)&(error > 0.)
        relative = error[selected]/value[selected]
        report[label] = dict(positive_nonzero_error_bins=int(selected.sum()),
            bins_below_two_percent=int(np.sum(relative < .02)),
            relative_error_quantiles=clean(np.quantile(relative, [.1, .5, .9])) if relative.size else None,
            bins_with_undefined_error=int(np.sum(~np.isfinite(error))),
            bins_with_nonpositive_central=int(np.sum(value <= 0.)))
    return report


def run(path, output, binning_path=None, locals=DEFAULT_SPECTRA, full_flavour=True):
    path, output = Path(path).resolve(), Path(output).resolve()
    directory = output.parent/(output.stem+'_spectra')
    if output.exists() or output.with_suffix('.npz').exists() or directory.exists():
        raise ValueError('Refuse to overwrite parameter reductions or partial spectrum archives')
    locals = tuple(locals)
    if not locals or len(set(locals)) != len(locals) or any(i not in range(2, 22) for i in locals):
        raise ValueError('Select distinct spectrum IDs 2--21')
    source, layout, ensembles, conditional = read_vectors(path)
    comparisons = source['comparisons']
    for comparison in comparisons:
        comparison_keys(comparison, full_flavour)
    charge_pairs, unpaired_charges = charge_statistics.pair_comparisons(comparisons, full_flavour)
    # Reject stale/unreferenced groups before any output is created.
    selected = select(ensembles, comparisons, full_flavour)
    if set(selected) != set(ensembles):
        raise ValueError('Unreferenced parameter retraining group')
    binning = None
    if binning_path is not None:
        content = Path(binning_path).read_bytes()
        plan = json.loads(content)
        layout, ensembles, conditional, binning = apply_common_binning(layout, ensembles, conditional, plan)
        binning.update(path=str(Path(binning_path).resolve()), sha256=hashlib.sha256(content).hexdigest(), plan=plan)
    banks = rate_banks(layout)
    if charge_pairs:
        for config in CONFIGS:
            for local in locals:
                charge_edges = []
                for charge in ('plus', 'minus'):
                    h = banks[config, charge]+local-1
                    a, b = map(int, layout['offsets'][h:h+2])
                    charge_edges.append(layout['edges'][a:b])
                if not np.array_equal(*charge_edges):
                    raise ValueError('Paired charge spectra have different bin boundaries')
    directory.mkdir(parents=True)
    arrays = dict(layout, contrast_labels=LABELS, charge_contrast_labels=charge_statistics.LABELS,
        charge_quantity_labels=charge_statistics.CHARGE_LABELS,
        charge_shape_quantity_labels=charge_statistics.CHARGE_SHAPE_LABELS,
        charge_derived_quantity_labels=charge_statistics.CHARGE_DERIVED_LABELS)
    report = dict(created_utc=now(), status='parameter observables reduced; precision and sensitivity inspection required',
        input=str(path), input_sha256=digest(path), input_arrays_sha256=source['arrays_sha256'],
        benchmark_sha256=source['benchmark_sha256'], comparisons=comparisons,
        flavour_scope='all eight ordered direct e/mu assignments' if full_flavour else 'explicit selected flavour sums',
        full_flavour=full_flavour, local_histogram_ids=locals, binning=binning,
        rates={}, derived_rates={}, spectra={}, shared_nominal_covariance={}, retraining_diagnostics={},
        charge_observables=dict(pairs=charge_pairs, unpaired=unpaired_charges, rates={}, derived_rates={}, spectra={},
            convention='Sum raw flavour and charge cross sections before combined normalization. '
                       'Charge ratios/asymmetries use matching scale coordinates. Store absolute '
                       'S/Pi and parameter changes, including at zero or negative asymmetry. '
                       'The asymmetry of an acceptance or normalized shape is labelled separately '
                       'from the rate charge asymmetry A_W.'),
        convention='Average independent retrainings within each flavour, then sum actual flavours. '
                    'Normalize and compare at matching scale coordinates before envelopes. Shared '
                    'baseline groups enter joint nominal covariance once. Conditional noise is not added again.',
        limitations='Precision flags do not prove coverage or adaptive-stop stability. Retain original '
                    'continuous ranges without folded overflow. PDF member 0 only: PDF uncertainty '
                    'belongs to the main member-wise analysis; PDF-family spreads and +/-1 GeV top-mass '
                    'scans are not additional Gaussian uncertainties. Physical parameter and full-study '
                    'conclusions require the selected integrations and their convergence inspection.')
    def keep(prefix, result, joint=False):
        arrays[prefix] = result['value']
        arrays[prefix+'_mc_errors'] = result['mc_error']
        arrays[prefix+'_covariance_factor_nominal'] = nominal_copy(result['covariance_factor'])
        arrays[prefix+'_nonlinear_bias_nominal'] = nominal_copy(result['nonlinear_bias_estimate'])
        if joint:
            report['shared_nominal_covariance'][prefix] = dict(
                comparisons=comparisons, contrast_labels=LABELS,
                ensemble_order=clean(result['ensemble_order']),
                factor_array=prefix+'_covariance_factor_nominal', shape=arrays[prefix+'_covariance_factor_nominal'].shape)
    def joint(prefix, values, transform=None):
        result = estimate(values, comparisons, transform, full_flavour)
        keep(prefix, result, joint=True)
    nominal_rates = features(ensembles, layout, banks, CONFIGS, nominal=True)
    joint('joint_all_rates_nominal', nominal_rates)
    joint('joint_all_derived_nominal', nominal_rates,
          lambda row: np.asarray([derived_rates(block) for block in row.reshape(len(CONFIGS), 9, 1)]))
    report['shared_nominal_covariance']['joint_all_rates_nominal']['feature_labels'] = [
        config+'/'+rate for config in CONFIGS for rate in RATES]
    report['shared_nominal_covariance']['joint_all_derived_nominal'].update(configs=CONFIGS, derived_labels=DERIVED)
    del nominal_rates
    for c, comparison in enumerate(comparisons):
        subset = select(ensembles, [comparison], full_flavour)
        identifier = '%s__%s__%s' % (comparison['scenario'], comparison['w_treatment'], comparison['charge'])
        report['rates'][identifier], report['derived_rates'][identifier], report['spectra'][identifier] = {}, {}, {}
        for config in CONFIGS:
            data = features(subset, layout, banks, [config])
            result = one_comparison(estimate(data, [comparison], full_flavour=full_flavour))
            prefix = 'c%d_%s_rates' % (c, config)
            keep(prefix, result)
            report['rates'][identifier][config] = {name: summaries(result, i) for i, name in enumerate(RATES)}
            result = one_comparison(estimate(data, [comparison], derived_rates, full_flavour))
            keep(prefix+'_derived', result)
            report['derived_rates'][identifier][config] = {name: summaries(result, i) for i, name in enumerate(DERIVED)}
            del result, data
            for local in locals:
                data = features(subset, layout, banks, [config], local)
                absolute = one_comparison(estimate({key: value[:, :-1] for key, value in data.items()},
                                                    [comparison], full_flavour=full_flavour))
                normalized = one_comparison(estimate(data, [comparison], normalize, full_flavour))
                coverage = one_comparison(estimate(data, [comparison], lambda row: normalize(row).sum(axis=0)[None], full_flavour))
                prefix = 'c%d_%s_h%d' % (c, config, local)
                keep(prefix+'_absolute', absolute)
                keep(prefix+'_normalized', normalized)
                h = banks[config, comparison['charge']]+local-1
                a, b = map(int, layout['offsets'][h:h+2])
                detail = dict(title=layout['titles'][h], comparison=comparison, edges=layout['edges'][a:b],
                    normalization_rate_bins=[RATES[i] for i in normalization_bins(local)],
                    units='Absolute cross sections per bin in pb; normalized entries are bin fractions.',
                    absolute=[summaries(absolute, i) for i in range(b-a)],
                    normalized=[summaries(normalized, i) for i in range(b-a)],
                    measured_range_coverage=summaries(coverage, 0), precision=precision(normalized))
                target = directory/(prefix+'.json.gz')
                target.write_bytes(gzip.compress(json.dumps(clean(detail), sort_keys=True, allow_nan=False).encode(), mtime=0))
                report['spectra'][identifier][config+'_h'+str(local)] = dict(path=str(target), sha256=digest(target),
                    precision=detail['precision'])
                del absolute, normalized, coverage, detail, data
        print('Reduced parameter rates and spectra:', identifier, flush=True)
    # Cross-scenario/charge covariance uses only the nominal coordinate here;
    # complete retraining vectors retain all cross-weight information.
    for config in CONFIGS:
        for local in locals:
            data = features(ensembles, layout, banks, [config], local, nominal=True)
            joint('joint_%s_h%d_absolute_nominal' % (config, local),
                  {key: value[:, :-1] for key, value in data.items()})
            joint('joint_%s_h%d_normalized_nominal' % (config, local), data, normalize)
    def keep_charge(prefix, result, pairs=None, quantity_labels=None):
        keep(prefix, result)
        if pairs is not None:
            report['shared_nominal_covariance'][prefix] = dict(
                pairs=pairs, quantity_labels=quantity_labels, contrast_labels=charge_statistics.LABELS,
                ensemble_order=clean(result['ensemble_order']),
                factor_array=prefix+'_covariance_factor_nominal', shape=arrays[prefix+'_covariance_factor_nominal'].shape)
    if charge_pairs:
        paired = [comparison for pair in charge_pairs for comparison in charge_statistics.paired_comparisons(pair)]
        charge_values = select(ensembles, paired, full_flavour)
        nominal_rates = features(charge_values, layout, banks, CONFIGS, nominal=True)
        result = charge_statistics.estimate(nominal_rates, charge_pairs, full_flavour=full_flavour)
        keep_charge('joint_charge_all_rates_nominal', result, charge_pairs, charge_statistics.CHARGE_LABELS)
        report['shared_nominal_covariance']['joint_charge_all_rates_nominal']['feature_labels'] = [
            config+'/'+rate for config in CONFIGS for rate in RATES]
        result = charge_statistics.estimate(nominal_rates, charge_pairs,
            lambda row: np.asarray([derived_rates(block) for block in row.reshape(len(CONFIGS), 9, 1)]), full_flavour)
        keep_charge('joint_charge_all_derived_nominal', result, charge_pairs, charge_statistics.CHARGE_DERIVED_LABELS)
        report['shared_nominal_covariance']['joint_charge_all_derived_nominal'].update(configs=CONFIGS, derived_labels=DERIVED)
        del result, nominal_rates
        for c, pair in enumerate(charge_pairs):
            subset = select(charge_values, charge_statistics.paired_comparisons(pair), full_flavour)
            identifier = '%s__%s' % (pair['scenario'], pair['w_treatment'])
            section = report['charge_observables']
            section['rates'][identifier], section['derived_rates'][identifier], section['spectra'][identifier] = {}, {}, {}
            for config in CONFIGS:
                data = features(subset, layout, banks, [config])
                result = one_comparison(charge_statistics.estimate(data, [pair], full_flavour=full_flavour))
                prefix = 'charge_c%d_%s_rates' % (c, config)
                keep_charge(prefix, result)
                section['rates'][identifier][config] = {name: charge_summaries(result, i, charge_statistics.CHARGE_LABELS)
                                                         for i, name in enumerate(RATES)}
                result = one_comparison(charge_statistics.estimate(data, [pair], derived_rates, full_flavour))
                keep_charge(prefix+'_derived', result)
                section['derived_rates'][identifier][config] = {name: charge_summaries(result, i, charge_statistics.CHARGE_DERIVED_LABELS)
                                                                 for i, name in enumerate(DERIVED)}
                del result, data
                for local in locals:
                    data = features(subset, layout, banks, [config], local)
                    absolute = one_comparison(charge_statistics.estimate({key: value[:, :-1] for key, value in data.items()},
                                                                        [pair], full_flavour=full_flavour))
                    normalized = one_comparison(charge_statistics.estimate(data, [pair], normalize, full_flavour))
                    coverage = one_comparison(charge_statistics.estimate(data, [pair],
                        lambda row: normalize(row).sum(axis=0)[None], full_flavour))
                    prefix = 'charge_c%d_%s_h%d' % (c, config, local)
                    keep_charge(prefix+'_absolute', absolute)
                    keep_charge(prefix+'_normalized', normalized)
                    hp, hm = (banks[config, charge]+local-1 for charge in ('plus', 'minus'))
                    a, b = map(int, layout['offsets'][hp:hp+2])
                    detail = dict(pair=pair, titles=[layout['titles'][h] for h in (hp, hm)], edges=layout['edges'][a:b],
                        normalization_rate_bins=[RATES[i] for i in normalization_bins(local)],
                        units='Absolute charge sums and their changes are pb per bin; charge ratios, asymmetries, '
                              'normalized quantities and their changes are dimensionless.',
                        absolute=[charge_summaries(absolute, i, charge_statistics.CHARGE_LABELS) for i in range(b-a)],
                        normalized=[charge_summaries(normalized, i, charge_statistics.CHARGE_SHAPE_LABELS) for i in range(b-a)],
                        measured_range_coverage=charge_summaries(coverage, 0, charge_statistics.CHARGE_DERIVED_LABELS))
                    target = directory/(prefix+'.json.gz')
                    target.write_bytes(gzip.compress(json.dumps(clean(detail), sort_keys=True, allow_nan=False).encode(), mtime=0))
                    section['spectra'][identifier][config+'_h'+str(local)] = dict(path=str(target), sha256=digest(target))
                    del data, absolute, normalized, coverage, detail
            print('Reduced parameter charge observables:', identifier, flush=True)
        for config in CONFIGS:
            for local in locals:
                data = features(charge_values, layout, banks, [config], local, nominal=True)
                result = charge_statistics.estimate({key: value[:, :-1] for key, value in data.items()},
                                                    charge_pairs, full_flavour=full_flavour)
                keep_charge('joint_charge_%s_h%d_absolute_nominal' % (config, local), result,
                            charge_pairs, charge_statistics.CHARGE_LABELS)
                result = charge_statistics.estimate(data, charge_pairs, normalize, full_flavour)
                keep_charge('joint_charge_%s_h%d_normalized_nominal' % (config, local), result,
                            charge_pairs, charge_statistics.CHARGE_SHAPE_LABELS)
                del result, data
    for key, values in ensembles.items():
        h = banks['R04_b25', key[3]]
        start = int(layout['offsets'][h])
        indices = [start+i for i in (0, 3, 4)]
        count = len(values)
        between = np.var(values[:, indices, 0], axis=0, ddof=1)/count
        within = conditional[key][:, indices].sum(axis=0)/count**2
        name = '__'.join(key[:-1]+(''.join(key[-1]),))
        report['retraining_diagnostics'][name] = dict(retrainings=count,
            rate_labels=['input', 'fiducial_1b', 'fiducial_2b'], between_retraining_variance_of_mean=between,
            conditional_variance_of_combination=within, variance_ratio=ratio(between, within),
            interpretation='Diagnostic comparison only; no variance addition or automatic error rescaling.')
    np.savez_compressed(output.with_suffix('.npz'), **arrays)
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),
        source_hashes={name:digest(STUDY/'scripts'/name) for name in ('parameter_report.py', 'parameter_statistics.py', 'parameter_charge_statistics.py',
            'parameter_replicas.py', 'parameter_variation_results.py', 'covariance_storage.py', 'main_binning.py',
            'rebinning.py', 'joint_report.py', 'pilot_report.py', 'replica_statistics.py')})
    report['source_hashes']['Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py'] = digest(
        ROOT/'Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py')
    save(output, clean(report))
    print(output, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--binning', type=Path)
    parser.add_argument('--histograms', type=int, nargs='+', default=DEFAULT_SPECTRA)
    parser.add_argument('--selected-flavours', action='store_true', help='Explicitly label selected flavour subsets instead of requiring all eight')
    args = parser.parse_args()
    run(args.input, args.output, args.binning, args.histograms, not args.selected_flavours)
