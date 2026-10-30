#!/usr/bin/env python3
"""Compare complete fixed-count production hard-split and baseline ensembles."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {key: clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def validate_inventory(record):
    jobs = record['jobs']
    expected = {(variant, scale, replica, group)
                for variant in ('S', 'Pi') for scale in (0.1, 0.3)
                for replica in range(4) for group in (0, 1)}
    actual = [(job['variant'], job['scale'], job['replica'], job['group'])
              for job in jobs]
    if len(jobs) != 32 or set(actual) != expected or len(set(actual)) != 32:
        raise ValueError('Expected 32 distinct S/Pi × h=0.1/0.3 × replica0..3 × group0/1 jobs')
    if any(job.get('status') != 'finished' or job.get('returncode') != 0
           for job in jobs):
        raise ValueError('All 32 fixed-count jobs must finish before analysis')
    if len({job['seed'] for job in jobs}) != 32:
        raise ValueError('Independent jobs require distinct seeds')
    if any(job['iterations'] != 1 for job in jobs):
        raise ValueError('This reduction requires one fixed-count iteration per job')
    return jobs


def native_result(directory, metadata=None):
    """Read total first rows; subsequent rows describe individual channels."""
    directory = Path(directory)
    res = directory/'res.dat'
    results = directory/'results.dat'
    if res.is_file():
        row = list(map(float, res.read_text().splitlines()[0].split()))
        value, error, cpu = row[2], row[3], row[6]
    else:
        row = list(map(float, results.read_text().splitlines()[0].split()))
        value, error = row[9], row[1]
        match = re.search(r'Time spent in Total\s*:\s*(\S+)',
                          (directory/'log.txt').read_text())
        if not match:
            raise ValueError('Missing native CPU total: ' + str(directory))
        cpu = float(match[1])
    if metadata:
        for name, actual in [('native_value_pb', value),
                             ('native_error_pb', error), ('native_cpu_s', cpu)]:
            if name in metadata and not math.isclose(metadata[name], actual,
                                                     rel_tol=1e-10, abs_tol=1e-16):
                raise ValueError('Cached native total differs: ' + str(directory) + ' ' + name)
    if not all(math.isfinite(item) for item in (value, error, cpu)) or error < 0 or cpu < 0:
        raise ValueError('Invalid native total: ' + str(directory))
    return dict(value_pb=value, error_pb=error, cpu_s=cpu, directory=str(directory))


def estimate(values, variances):
    """Arithmetic independent-replica mean, native SEM and scatter SEM."""
    values, variances = np.asarray(values), np.asarray(variances)
    if values.shape != variances.shape or len(values) < 2:
        raise ValueError('Matching values/variances and at least two replicas required')
    n = len(values)
    return dict(mean=np.mean(values, axis=0),
                native_variance=np.sum(variances, axis=0)/n**2,
                scatter_variance=np.var(values, axis=0, ddof=1)/n)


def rate_summary(rows, training_cpu=0):
    values = np.array([row['value_pb'] for row in rows])
    variances = np.array([row['error_pb']**2 for row in rows])
    result = estimate(values, variances)
    variance = float(result['native_variance'])
    scatter_variance = float(result['scatter_variance'])
    cpu = sum(row['cpu_s'] for row in rows)
    return dict(replicas=len(rows), mean_pb=float(result['mean']),
        native_sem_pb=math.sqrt(variance), scatter_sem_pb=math.sqrt(scatter_variance),
        native_variance_times_cpu_pb2_s=variance*cpu,
        scatter_variance_times_cpu_pb2_s=scatter_variance*cpu,
        total_native_cpu_s=cpu, training_cpu_s=training_cpu,
        native_variance_times_cpu_with_training_pb2_s=variance*(cpu+training_cpu),
        largest_replica_native_variance_fraction=float(max(variances)/sum(variances))
            if sum(variances) else None,
        values_pb=values.tolist(), errors_pb=np.sqrt(variances).tolist(), samples=rows)


def read_histograms(directory):
    from madgraph.various.histograms import HwUList
    path = Path(directory)/'MADatNLO.HwU'
    histograms = HwUList(str(path), raw_labels=True)
    titles, edges, values, variances, offsets = [], [], [], [], [0]
    for histogram in histograms:
        titles.append(histogram.title.strip())
        for row in histogram.bins:
            labels = {key.strip().lower(): key for key in row.wgts}
            central = next((labels[key] for key in ('central value', 'central')
                            if key in labels), None)
            error = next((labels[key] for key in ('dy', 'stat_error')
                          if key in labels), None)
            if central is None or error is None:
                raise ValueError('Missing central value/dy in ' + str(path))
            edges.append(row.boundaries)
            values.append(row.wgts[central])
            variances.append(row.wgts[error]**2)
        offsets.append(len(values))
    return dict(titles=titles, edges=np.array(edges), offsets=np.array(offsets)), \
        np.array(values), np.array(variances), str(path)


def same_layout(reference, candidate):
    return (reference['titles'] == candidate['titles'] and
            np.array_equal(reference['edges'], candidate['edges']) and
            np.array_equal(reference['offsets'], candidate['offsets']))


def contrast(first, second):
    difference = first['mean_pb'] - second['mean_pb']
    error = math.hypot(first['native_sem_pb'], second['native_sem_pb'])
    return dict(difference_pb=difference, native_error_pb=error,
                descriptive_pull=difference/error if error else None,
                first_over_second=first['mean_pb']/second['mean_pb']
                    if second['mean_pb'] else None)


def differential_closure(layout, first, second):
    reports = []
    for index, title in enumerate(layout['titles']):
        if 'rates:' in title:
            continue
        start, end = layout['offsets'][index:index+2]
        a, b = first['mean'][start:end], second['mean'][start:end]
        # A population selection is descriptive, not a statistical acceptance cut.
        mask = ((a > 0) & (b > 0) &
                (a >= 0.005*np.sum(np.maximum(a, 0))) &
                (b >= 0.005*np.sum(np.maximum(b, 0))))
        errors = np.sqrt(first['native_variance'][start:end] +
                         second['native_variance'][start:end])
        mask &= errors > 0
        if not np.any(mask):
            continue
        pulls = (a[mask] - b[mask])/errors[mask]
        positions = np.flatnonzero(mask)
        maximum = int(np.argmax(np.abs(pulls)))
        reports.append(dict(title=title, populated_bins=len(pulls),
            median_absolute_descriptive_pull=float(np.median(np.abs(pulls))),
            p90_absolute_descriptive_pull=float(np.quantile(np.abs(pulls), 0.9)),
            maximum_absolute_descriptive_pull=float(np.max(np.abs(pulls))),
            maximum_bin_edges=layout['edges'][start+positions[maximum]].tolist(),
            median_relative_difference=float(np.median((a[mask]-b[mask])/b[mask]))))
    return reports


def make_plots(directory, rates, spectra, layout, legacy_pi_baseline=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    figure, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    labels = ['baseline', 'h=0.1', 'h=0.3']
    for variant, offset, color in [('S', -0.05, 'tab:blue'), ('Pi', 0.05, 'tab:orange')]:
        positions = [1, 2] if variant == 'Pi' and legacy_pi_baseline else [0, 1, 2]
        keys = [variant+'_'+['baseline', 'h10', 'h30'][i] for i in positions]
        rows = [rates[key] for key in keys]
        axes[0, 0].errorbar(np.array(positions)+offset,
            [row['mean_pb']*1000 for row in rows],
            yerr=[row['native_sem_pb']*1000 for row in rows], fmt='o',
            color=color, capsize=3, label=variant+' native SEM')
        axes[0, 1].plot(np.array(positions)+offset,
            [row['native_variance_times_cpu_pb2_s'] for row in rows],
            'o-', color=color, label=variant)
        if any(row['training_cpu_s'] for row in rows):
            axes[0, 1].plot(np.array(positions)+offset,
                [row['native_variance_times_cpu_with_training_pb2_s'] for row in rows],
                ':', color=color, label=variant+' including new training')
    axes[0, 0].set_ylabel('One-beam inclusive cross section [fb]')
    axes[0, 0].set_title('Fixed-count arithmetic means')
    axes[0, 1].set_ylabel(r'Native variance × CPU [pb$^2$ s]')
    axes[0, 1].set_yscale('log')
    axes[0, 1].set_title('Smaller values indicate lower measured cost')
    for ax in axes[0]:
        ax.set_xticks(range(3), labels)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.2)
    examples = [('S', 'R04_b25 W+ 1b HTjets', r'$H_T^{jets}$ [GeV]'),
                ('Pi', 'R04_b25 W+ 2b subleading bjet pT', r'$p_T(b_2)$ [GeV]')]
    missing = []
    for ax, (variant, title, xlabel) in zip(axes[1], examples):
        if title not in layout['titles']:
            ax.set_visible(False)
            missing.append(title)
            continue
        index = layout['titles'].index(title)
        start, end = layout['offsets'][index:index+2]
        boundaries = layout['edges'][start:end]
        centers, widths = np.mean(boundaries, axis=1), np.diff(boundaries, axis=1)[:, 0]
        positive_values = []
        for suffix, label, color in [('baseline', 'baseline', '0.4'),
                                     ('h10', 'h=0.1', 'tab:blue'),
                                     ('h30', 'h=0.3', 'tab:orange')]:
            if variant == 'Pi' and suffix == 'baseline' and legacy_pi_baseline:
                continue
            row = spectra[variant+'_'+suffix]
            values = row['mean'][start:end]*1000/widths
            errors = np.sqrt(row['native_variance'][start:end])*1000/widths
            mask = values > 0
            positive_values.append(values[mask])
            ax.errorbar(centers[mask], values[mask], yerr=errors[mask],
                        fmt='.', ms=3, color=color, alpha=0.85, label=label)
            ax.fill_between(centers, np.where(mask, np.maximum(values-errors, 1e-15), np.nan),
                            np.where(mask, values+errors, np.nan),
                            color=color, alpha=0.12)
        ax.set_yscale('log')
        positive = np.concatenate(positive_values)
        if positive.size:
            maximum = float(np.max(positive))
            ax.set_ylim(max(float(np.min(positive))*0.4, maximum*1e-5),
                        maximum*3)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$d\sigma/dx$ [fb/GeV]')
        ax.set_title(variant+': '+title)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.2)
    subtitle = ('Legacy Pi reference excluded: production-real spin frames corrected'
                if legacy_pi_baseline else
                'Native MC errors; four split replicas versus eight baseline replicas')
    figure.suptitle('Production hard-split diagnostic: central weights, one beam ordering\n'+subtitle)
    figure.tight_layout()
    paths = []
    for extension in ('pdf', 'png'):
        path = directory/('production_hard_split.'+extension)
        figure.savefig(path, dpi=160)
        paths.append(str(path.resolve()))
    plt.close(figure)
    return dict(files=paths, missing_examples=missing,
                positive_bins_only_for_log_plots=True,
                omitted_legacy_pi_baseline=legacy_pi_baseline)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path,
        default=Path('ttw_study/checks/hard_split_fixed_count_20261014.json'))
    parser.add_argument('--baseline', type=Path,
        default=Path('ttw_study/results/tail_pilot_final_20261014.json'))
    parser.add_argument('--training-record', type=Path)
    parser.add_argument('--output', type=Path,
        default=Path('ttw_study/results/production_hard_split_20261014.json'))
    parser.add_argument('--plot-dir', type=Path)
    parser.add_argument('--legacy-pi-baseline', action='store_true',
        help='retain old Pi data as historical only; omit its baseline pulls, efficiency ratios and plot points')
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.npz').exists():
        raise ValueError('Choose fresh output filenames')
    record, baseline = [json.loads(path.read_text()) for path in (args.record, args.baseline)]
    legacy_pi_baseline = bool(args.legacy_pi_baseline or
                              record.get('legacy_pi_baseline_incomparable', False))
    jobs = validate_inventory(record)
    training_cost = {}
    if args.training_record:
        for job in json.loads(args.training_record.read_text())['jobs']:
            if job.get('status') != 'finished' or job.get('returncode') != 0:
                raise ValueError('Training cost record contains unfinished jobs')
            key = (job['variant'], job['scale'])
            training_cost[key] = training_cost.get(key, 0) + native_result(job['directory'])['cpu_s']
    layout = None
    paths, grouped = [], {}

    def load(directory):
        nonlocal layout
        candidate, values, variances, path = read_histograms(directory)
        if layout is None:
            layout = candidate
        elif not same_layout(layout, candidate):
            raise ValueError('Histogram layout differs: ' + path)
        paths.append(path)
        return values, variances

    for job in jobs:
        key = (job['variant'], job['scale'], job['replica'])
        values, variances = load(job['directory'])
        grouped.setdefault(key, []).append((native_result(job['directory'], job), values, variances))
    rates, spectra, vectors = {}, {}, {}
    for variant in ('S', 'Pi'):
        for scale in (0.1, 0.3):
            rows, all_values, all_variances = [], [], []
            for replica in range(4):
                parts = grouped[(variant, scale, replica)]
                rows.append(dict(value_pb=sum(part[0]['value_pb'] for part in parts),
                    error_pb=math.sqrt(sum(part[0]['error_pb']**2 for part in parts)),
                    cpu_s=sum(part[0]['cpu_s'] for part in parts),
                    replica=replica, groups=[part[0] for part in parts]))
                all_values.append(sum(part[1] for part in parts))
                all_variances.append(sum(part[2] for part in parts))
            label = variant+'_h%d' % round(scale*100)
            rates[label] = rate_summary(rows, training_cost.get((variant, scale), 0))
            spectra[label] = estimate(all_values, all_variances)
            vectors[label+'_replica_values_pb'] = np.array(all_values)
            vectors[label+'_replica_native_variances_pb2'] = np.array(all_variances)
        selected = [worker for worker in baseline['workers']
            if worker['metadata'].get('campaign_job', {}).get('execution_label') == 'initial'
            and worker['metadata']['campaign_job']['variant'] == variant
            and not worker['metadata']['campaign_job'].get('same_stream_replay', False)]
        if len(selected) != 8:
            raise ValueError('Expected eight initial baseline workers per prescription')
        rows, all_values, all_variances = [], [], []
        for worker in sorted(selected, key=lambda row: row['metadata']['campaign_job']['replica']):
            directory = worker['metadata']['campaign_job']['directory']
            rows.append(native_result(directory))
            values, variances = load(directory)
            all_values.append(values)
            all_variances.append(variances)
        label = variant+'_baseline'
        rates[label] = rate_summary(rows)
        rates[label]['comparable_reference'] = not (variant == 'Pi' and legacy_pi_baseline)
        if not rates[label]['comparable_reference']:
            rates[label]['reference_limitation'] = (
                'Legacy Pi used production-real and underlying-Born decay spin densities '
                'without consistent frame alignment. It is not a normalization or '
                'pure sampling-efficiency reference for the corrected integrand.')
        spectra[label] = estimate(all_values, all_variances)
        vectors[label+'_replica_values_pb'] = np.array(all_values)
        vectors[label+'_replica_native_variances_pb2'] = np.array(all_variances)
    closure = {}
    for variant in ('S', 'Pi'):
        closure[variant] = dict(h10_minus_h30=contrast(rates[variant+'_h10'], rates[variant+'_h30']),
            differential_h10_minus_h30=differential_closure(
                layout, spectra[variant+'_h10'], spectra[variant+'_h30']))
        comparable = rates[variant+'_baseline']['comparable_reference']
        closure[variant]['baseline_comparable'] = comparable
        for suffix in ('h10', 'h30'):
            closure[variant][suffix+'_minus_baseline'] = (
                contrast(rates[variant+'_'+suffix], rates[variant+'_baseline'])
                if comparable else None)
        if not comparable:
            closure[variant]['baseline_comparison_unavailable_reason'] = \
                rates[variant+'_baseline']['reference_limitation']
        baseline_cost = rates[variant+'_baseline']['native_variance_times_cpu_pb2_s']
        for label in ('h10', 'h30'):
            row = rates[variant+'_'+label]
            row['native_variance_cpu_ratio_to_baseline'] = (
                row['native_variance_times_cpu_pb2_s']/baseline_cost if comparable else None)
            row['native_variance_cpu_with_training_ratio_to_baseline'] = \
                (row['native_variance_times_cpu_with_training_pb2_s']/baseline_cost
                 if comparable else None)
    for label, row in spectra.items():
        vectors.update({label+'_'+key: value for key, value in row.items()})
    vectors.update(edges=layout['edges'], offsets=layout['offsets'], titles=np.array(layout['titles']))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    npz_path = args.output.with_suffix('.npz')
    np.savez_compressed(npz_path, **vectors)
    figures = make_plots(args.plot_dir, rates, spectra, layout, legacy_pi_baseline) if args.plot_dir else {}
    manifest = dict(schema='production_hard_split_diagnostic_v1',
        created_utc=datetime.now(timezone.utc).isoformat(),
        record=str(args.record.resolve()), record_sha256=digest(args.record),
        baseline=str(args.baseline.resolve()), baseline_sha256=digest(args.baseline),
        training_record=str(args.training_record.resolve()) if args.training_record else None,
        legacy_pi_baseline_incomparable=legacy_pi_baseline,
        script_sha256=digest(__file__), scope=record.get('scope'),
        fixed_count_jobs=32, rates=rates, comparisons=closure,
        bin_population_definition='Both positive bin yields >=0.5% of each histogram positive-yield sum; rates histograms excluded.',
        histogram_paths=paths, complete_bin_vectors=str(npz_path.resolve()),
        complete_bin_vectors_sha256=digest(npz_path), figures=figures,
        limitations=[
            'Central weights only, one W+ ordered lepton channel and one production beam ordering; no PDF/scale bands or full-flavour prediction.',
            'Four split replicas versus eight baseline replicas provide noisy heavy-tail variance estimates; descriptive pulls are not Gaussian pass/fail tests.',
            'Bin errors are native MC errors conditional on trained proposals. Bins are correlated; no independent-bin goodness-of-fit interpretation is applied.',
            'Training cost, when requested, is added only for the new split proposals; baseline historical training cost is not available.',
            'All signed bins are retained in NPZ; nonpositive bins are omitted only from logarithmic example plots.'])
    args.output.write_text(json.dumps(clean(manifest), indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: {name: row[name] for name in
        ('mean_pb', 'native_sem_pb', 'scatter_sem_pb', 'total_native_cpu_s',
         'native_variance_times_cpu_pb2_s')} for key, row in rates.items()}, indent=2))


if __name__ == '__main__':
    main()
