#!/usr/bin/env python3
"""Full-flavour W/mass observables from complete independent retrainings.

Physical states shared by comparisons occur once in the joint covariance.
Scale/PDF reductions follow all flavour sums and nonlinear observables.
Spectrum archives are written separately to bound retained output memory.
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
from load_results import WEIGHTS
from main_binning import apply_plan
from pilot_report import RATES, clean, describe as describe_weights
from replica_statistics import ratio
from robustness_replicas import read as read_vectors
import robustness_statistics as statistics

CONFIGS = ('R04_b25', 'R04_b30', 'R04_b40', 'R03_b25', 'R05_b25')
DERIVED = ('acceptance_1b', 'acceptance_2b', 'two_b_given_one_b',
           'veto_0extra', 'fraction_1extra', 'fraction_2extra', 'fraction_3plus_extra')
DEFAULT_SPECTRA = (3, 4, 6, 10, 11, 12, 20)
MIGRATIONS = tuple(config+'/'+rate+'_over_'+CONFIGS[0]
                   for config in CONFIGS[1:] for rate in RATES[3:5])


def rate_banks(layout):
    offsets = np.asarray(layout['offsets'])
    if (layout['weights'] != WEIGHTS or len(layout['titles']) != 210
            or offsets.shape != (211,) or offsets[0] != 0
            or offsets[-1] != len(layout['edges']) or np.any(np.diff(offsets) <= 0)):
        raise ValueError('Require the complete ten-bank, 183-weight robustness layout')
    banks = {}
    for h, title in enumerate(layout['titles']):
        if ' rates:' not in title:
            continue
        config = title.split()[0]
        charge = 'plus' if 'W+' in title else 'minus' if 'W-' in title else None
        if h % 21 or offsets[h+1]-offsets[h] != 9 or (config, charge) in banks:
            raise ValueError('Duplicated, misplaced or incomplete robustness rate bank')
        banks[config, charge] = h
    if set(banks) != {(config, charge) for config in CONFIGS for charge in ('plus', 'minus')}:
        raise ValueError('Require all five cut configurations for both charges')
    return banks


def features(ensembles, layout, banks, configs, local=None, nominal=False):
    selected = {}
    for key, values in ensembles.items():
        mode, mass, charge, variant, flavours = key
        source = values[..., :1] if nominal else values
        rows = []
        for config in configs:
            h = banks[config, charge]
            start = int(layout['offsets'][h])
            if local is None:
                row = source[:, start:start+9]
            else:
                a, b = map(int, layout['offsets'][h+local-1:h+local+1])
                parents = [start+i for i in normalization_bins(local)]
                row = np.concatenate([source[:, a:b], source[:, parents].sum(axis=1)[:, None, :]], axis=1)
            rows.append(row)
        selected[key] = np.concatenate(rows, axis=1)
    return selected


def normalize(row):
    return ratio(row[:-1], row[-1])


def all_derived(row):
    return np.concatenate([derived_rates(block) for block in row.reshape(len(CONFIGS), 9, row.shape[-1])])


def migrations(row):
    rates = row.reshape(len(CONFIGS), 9, row.shape[-1])[:, 3:5]
    return ratio(rates[1:], rates[0]).reshape(len(MIGRATIONS), row.shape[-1])


def one_comparison(result):
    return {key: value[:, 0] if key == 'covariance_factor' else value[0]
            for key, value in result.items() if key != 'ensemble_order'}


def describe(value, error, bias, relative=False, change=False):
    row = describe_weights(value)
    row.update(nominal_retraining_mc_error=error[0], nominal_jackknife_bias_estimate=bias[0])
    if relative:
        row['nominal_relative_mc_error'] = ratio(error[0], abs(value[0]))
    if change:
        row.update(nominal_change_over_mc_error=ratio(value[0], error[0]),
            change_exceeds_three_mc_errors=bool(np.isfinite(error[0]) and error[0] > 0.
                                                and abs(value[0]) > 3.*error[0]))
    return row


def summaries(result, index, quantity_labels=None):
    if quantity_labels is not None:
        labels = statistics.CHARGE_CONTRAST_LABELS
        return {quantity: {label: describe(*(result[key][q, i, index] for key in
                    ('value', 'mc_error', 'nonlinear_bias_estimate')), change=i >= 4)
                for i, label in enumerate(labels)} for q, quantity in enumerate(quantity_labels)}
    return {label: describe(*(result[key][i, index] for key in
                ('value', 'mc_error', 'nonlinear_bias_estimate')), relative=i < 4,
                change=label in ('product_shift_change', 'product_ratio_change'))
            for i, label in enumerate(statistics.LABELS)}


def precision(result):
    report = {}
    for i, label in enumerate(statistics.LABELS[:4]):
        value, error = result['value'][i, :, 0], result['mc_error'][i, :, 0]
        selected = np.isfinite(value)&np.isfinite(error)&(value > 0.)&(error > 0.)
        relative = error[selected]/value[selected]
        report[label] = dict(positive_nonzero_error_bins=int(selected.sum()),
            bins_below_two_percent=int(np.sum(relative < .02)),
            relative_error_quantiles=np.quantile(relative, [.1, .5, .9]) if relative.size else None,
            bins_with_undefined_error=int(np.sum(~np.isfinite(error))),
            bins_with_nonpositive_central=int(np.sum(value <= 0.)))
    return report


def keep(arrays, prefix, result):
    arrays[prefix] = result['value']
    arrays[prefix+'_mc_errors'] = result['mc_error']
    arrays[prefix+'_covariance_factor_nominal'] = nominal_copy(result['covariance_factor'])
    arrays[prefix+'_nonlinear_bias_nominal'] = nominal_copy(result['nonlinear_bias_estimate'])


def archive_arrays(path, arrays):
    np.savez_compressed(path, **arrays)
    return dict(path=str(path), sha256=digest(path))


def archive_spectrum(directory, prefix, detail, arrays):
    array_record = archive_arrays(directory/(prefix+'.npz'), arrays)
    detail['arrays'] = array_record
    path = directory/(prefix+'.json.gz')
    path.write_bytes(gzip.compress(json.dumps(clean(detail), sort_keys=True, allow_nan=False).encode(), mtime=0))
    record = dict(path=str(path), sha256=digest(path), arrays=array_record)
    if 'precision' in detail:
        record['precision'] = detail['precision']
    return record


def run(path, output, binning_path=None, locals=DEFAULT_SPECTRA, full_flavour=True):
    path, output = Path(path).resolve(), Path(output).resolve()
    directory = output.parent/(output.stem+'_spectra')
    if output.exists() or output.with_suffix('.npz').exists() or directory.exists():
        raise ValueError('Refuse to overwrite robustness reductions or partial spectrum archives')
    locals = tuple(locals)
    if not locals or len(set(locals)) != len(locals) or any(i not in range(2, 22) for i in locals):
        raise ValueError('Select distinct spectrum IDs 2--21')
    source, layout, ensembles, conditional = read_vectors(path)
    comparisons = source['comparisons']
    for comparison in comparisons:
        statistics.comparison_keys(comparison, full_flavour)
    charge_pairs, unpaired = statistics.pair_comparisons(comparisons, full_flavour)
    if set(statistics.select(ensembles, comparisons, full_flavour)) != set(ensembles):
        raise ValueError('Unreferenced robustness retraining group')
    banks = rate_banks(layout)
    binning = None
    if binning_path is not None:
        content = Path(binning_path).read_bytes()
        plan = json.loads(content)
        layout, ensembles, conditional, binning = apply_plan(layout, ensembles, conditional, plan)
        binning.update(path=str(Path(binning_path).resolve()), sha256=hashlib.sha256(content).hexdigest(), plan=plan)
        banks = rate_banks(layout)
    if charge_pairs:
        for config in CONFIGS:
            for local in locals:
                edges = []
                for charge in ('plus', 'minus'):
                    h = banks[config, charge]+local-1
                    a, b = map(int, layout['offsets'][h:h+2])
                    edges.append(layout['edges'][a:b])
                if not np.array_equal(*edges):
                    raise ValueError('Paired charge spectra have different bin boundaries')
    directory.mkdir(parents=True)
    arrays = dict(layout, contrast_labels=statistics.LABELS,
        charge_contrast_labels=statistics.CHARGE_CONTRAST_LABELS,
        charge_quantity_labels=statistics.CHARGE_LABELS,
        charge_shape_quantity_labels=statistics.CHARGE_SHAPE_LABELS,
        charge_derived_quantity_labels=statistics.CHARGE_DERIVED_LABELS)
    def sections():
        return dict(rates={}, derived_rates={}, migrations={}, spectra={})
    report = dict(created_utc=now(),
        status='robustness observables reduced; precision and sensitivity inspection required',
        input=str(path), input_sha256=digest(path), input_arrays_sha256=source['arrays_sha256'],
        benchmark_sha256=source['benchmark_sha256'], comparisons=comparisons,
        flavour_scope='all eight ordered direct e/mu assignments' if full_flavour else 'explicit selected flavour sums',
        full_flavour=full_flavour, local_histogram_ids=locals, binning=binning,
        shared_nominal_covariance={}, retraining_diagnostics={}, **sections(),
        charge_observables=dict(pairs=charge_pairs, unpaired=unpaired, **sections(),
            convention='Sum raw flavour and charge cross sections before combined normalization. '
                       'Charge ratios/asymmetries use matching scale/PDF coordinates. Absolute '
                       'S/Pi and W/mass changes remain defined at zero or negative asymmetry. '
                       'Acceptance, migration and shape asymmetries are distinct from rate A_W.'),
        convention='Average complete independent retrainings within each flavour, then sum actual flavours. '
                   'Normalize and compare at matching scale/PDF coordinates before uncertainty reductions. '
                   'Shared physical states enter joint covariance once, including opposite signs for an '
                   'intermediate target/reference. Conditional noise is not added to retraining variance.',
        limitations='Precision flags do not prove coverage or adaptive-stop stability. Common bins retain '
                    'the original measured ranges and exclude folded overflow. PDF uncertainty, scale '
                    'envelopes and MC errors remain separate. W/mass conclusions require the real '
                    'campaign and convergence/sensitivity inspection. Full cross-weight correlations '
                    'remain available in the independent input vectors; stored covariance factors '
                    'use the nominal coordinate. Spectrum archives are separate from rate arrays.')
    for charged, items in ((False, comparisons), (True, charge_pairs)):
        if not items:
            continue
        selected_comparisons = ([row for pair in items for row in statistics.paired_comparisons(pair)]
                                if charged else items)
        selected = statistics.select(ensembles, selected_comparisons, full_flavour)
        estimator = statistics.estimate_charges if charged else statistics.estimate
        labels = statistics.CHARGE_CONTRAST_LABELS if charged else statistics.LABELS
        quantities = statistics.CHARGE_LABELS if charged else None
        shape_quantities = statistics.CHARGE_SHAPE_LABELS if charged else None
        derived_quantities = statistics.CHARGE_DERIVED_LABELS if charged else None
        section = report['charge_observables'] if charged else report
        stem = 'charge_' if charged else ''

        def joint(prefix, values, transform=None, quantity_labels=None, destination=None, feature_labels=None):
            result = estimator(values, items, transform, full_flavour)
            keep(arrays if destination is None else destination, prefix, result)
            metadata = dict(contrast_labels=labels, ensemble_order=clean(result['ensemble_order']),
                factor_array=prefix+'_covariance_factor_nominal',
                shape=result['covariance_factor'].shape[:-1],
                feature_labels=feature_labels, quantity_labels=quantity_labels)
            metadata['pairs' if charged else 'comparisons'] = items
            report['shared_nominal_covariance'][prefix] = metadata

        nominal = features(selected, layout, banks, CONFIGS, nominal=True)
        for name, transform, names, qlabels in (
                ('rates', None, [config+'/'+rate for config in CONFIGS for rate in RATES], quantities),
                ('derived', all_derived, [config+'/'+rate for config in CONFIGS for rate in DERIVED], derived_quantities),
                ('migrations', migrations, MIGRATIONS, derived_quantities)):
            joint('joint_'+stem+'all_'+name+'_nominal', nominal, transform, qlabels, feature_labels=names)
        del nominal
        for c, comparison in enumerate(items):
            members = statistics.paired_comparisons(comparison) if charged else [comparison]
            subset = statistics.select(selected, members, full_flavour)
            identifier = comparison['name'] if charged else comparison['name']+'__'+comparison['charge']
            for name in ('rates', 'derived_rates', 'spectra'):
                section[name][identifier] = {}
            all_rates = features(subset, layout, banks, CONFIGS)
            result = one_comparison(estimator(all_rates, [comparison], migrations, full_flavour))
            keep(arrays, '%sc%d_migrations' % (stem, c), result)
            section['migrations'][identifier] = {name: summaries(result, i, derived_quantities)
                                                 for i, name in enumerate(MIGRATIONS)}
            del all_rates, result
            for config in CONFIGS:
                data = features(subset, layout, banks, [config])
                prefix = '%sc%d_%s_rates' % (stem, c, config)
                result = one_comparison(estimator(data, [comparison], full_flavour=full_flavour))
                keep(arrays, prefix, result)
                section['rates'][identifier][config] = {name: summaries(result, i, quantities)
                                                        for i, name in enumerate(RATES)}
                result = one_comparison(estimator(data, [comparison], derived_rates, full_flavour))
                keep(arrays, prefix+'_derived', result)
                section['derived_rates'][identifier][config] = {name: summaries(result, i, derived_quantities)
                                                                for i, name in enumerate(DERIVED)}
                del data, result
                for local in locals:
                    data = features(subset, layout, banks, [config], local)
                    absolute = one_comparison(estimator({key: value[:, :-1] for key, value in data.items()},
                                                        [comparison], full_flavour=full_flavour))
                    normalized = one_comparison(estimator(data, [comparison], normalize, full_flavour))
                    coverage = one_comparison(estimator(data, [comparison],
                        lambda row: normalize(row).sum(axis=0)[None], full_flavour))
                    prefix = '%sc%d_%s_h%d' % (stem, c, config, local)
                    charges = ('plus', 'minus') if charged else (comparison['charge'],)
                    histograms = [banks[config, charge]+local-1 for charge in charges]
                    a, b = map(int, layout['offsets'][histograms[0]:histograms[0]+2])
                    spectrum_arrays = dict(edges=layout['edges'][a:b], weights=WEIGHTS, contrast_labels=labels)
                    if charged:
                        spectrum_arrays.update(absolute_quantity_labels=quantities,
                            normalized_quantity_labels=shape_quantities, coverage_quantity_labels=derived_quantities)
                    for name, result in (('absolute', absolute), ('normalized', normalized),
                                         ('measured_range_coverage', coverage)):
                        keep(spectrum_arrays, name, result)
                    detail = dict(edges=layout['edges'][a:b], titles=[layout['titles'][h] for h in histograms],
                        normalization_rate_bins=[RATES[i] for i in normalization_bins(local)],
                        units=('Absolute charge sums and their changes are pb per bin; charge ratios/asymmetries, '
                               'normalized quantities and their changes are dimensionless.' if charged else
                               'Absolute cross sections are pb per bin; normalized quantities are parent-rate bin fractions.'),
                        absolute=[summaries(absolute, i, quantities) for i in range(b-a)],
                        normalized=[summaries(normalized, i, shape_quantities) for i in range(b-a)],
                        measured_range_coverage=summaries(coverage, 0, derived_quantities))
                    detail['pair' if charged else 'comparison'] = comparison
                    detail['precision'] = (dict(combined=precision({key: normalized[key][0]
                        for key in ('value', 'mc_error')})) if charged else precision(normalized))
                    section['spectra'][identifier][config+'_h'+str(local)] = archive_spectrum(
                        directory, prefix, detail, spectrum_arrays)
                    del absolute, normalized, coverage, result, spectrum_arrays, detail, data
            print('Reduced robustness charge observables:' if charged else 'Reduced robustness observables:',
                  identifier, flush=True)
        # Keep each observable's joint nominal factors in its own archive.
        # Complete input vectors retain arbitrary cross-observable/weight covariance.
        for config in CONFIGS:
            for local in locals:
                data = features(selected, layout, banks, [config], local, nominal=True)
                prefix = 'joint_'+stem+config+'_h'+str(local)
                joint_arrays = dict(contrast_labels=labels)
                for charge in ('plus', 'minus'):
                    h = banks[config, charge]+local-1
                    a, b = map(int, layout['offsets'][h:h+2])
                    joint_arrays['edges_'+charge] = layout['edges'][a:b]
                for kind, values, transform, qlabels in (
                        ('absolute', {key: value[:, :-1] for key, value in data.items()}, None, quantities),
                        ('normalized', data, normalize, shape_quantities)):
                    joint(prefix+'_'+kind+'_nominal', values, transform, qlabels, destination=joint_arrays)
                location = archive_arrays(directory/(prefix+'.npz'), joint_arrays)
                for kind in ('absolute', 'normalized'):
                    report['shared_nominal_covariance'][prefix+'_'+kind+'_nominal']['arrays'] = location
                del data, values, joint_arrays
    for key, values in ensembles.items():
        mode, mass, charge, variant, flavours = key
        start = int(layout['offsets'][banks['R04_b25', charge]])
        indices = [start+i for i in (0, 3, 4)]
        count = len(values)
        between = np.var(values[:, indices, 0], axis=0, ddof=1)/count
        within = conditional[key][:, indices].sum(axis=0)/count**2
        name = '__'.join((mode, str(mass), charge, variant, ''.join(flavours)))
        report['retraining_diagnostics'][name] = dict(retrainings=count,
            rate_labels=['input', 'fiducial_1b', 'fiducial_2b'], between_retraining_variance_of_mean=between,
            conditional_variance_of_combination=within, variance_ratio=ratio(between, within),
            interpretation='Diagnostic comparison only; no variance addition or automatic error rescaling.')
    location = archive_arrays(output.with_suffix('.npz'), arrays)
    for metadata in report['shared_nominal_covariance'].values():
        metadata.setdefault('arrays', location)
    sources = ('robustness_report.py', 'robustness_statistics.py', 'robustness_replicas.py',
        'run_robustness_campaign.py', 'parameter_statistics.py', 'parameter_charge_statistics.py',
        'full_flavour_statistics.py', 'load_results.py', 'covariance_storage.py', 'main_binning.py',
        'rebinning.py', 'joint_report.py', 'pilot_report.py', 'pdf_statistics.py', 'replica_statistics.py')
    report.update(arrays_sha256=location['sha256'],
        source_hashes={name: digest(STUDY/'scripts'/name) for name in sources})
    scale_source = 'Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py'
    report['source_hashes'][scale_source] = digest(ROOT/scale_source)
    save(output, clean(report))
    print(output, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--binning', type=Path)
    parser.add_argument('--histograms', type=int, nargs='+', default=DEFAULT_SPECTRA)
    args = parser.parse_args()
    run(args.input, args.output, args.binning, args.histograms)
