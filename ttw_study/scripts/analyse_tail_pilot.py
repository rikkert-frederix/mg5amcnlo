#!/usr/bin/env python3
"""Summarise bounded ttW weight traces without altering or discarding samples.

Example: python3 analyse_tail_pilot.py PILOT --output report.json --plot-dir figures
Explicit replay pairs can be specified as --replay-pair WORKER REPLAY_WORKER.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

TYPE_NAMES = {1: 'real', 2: 'Born', 3: 'integrated', 4: 'soft',
              5: 'collinear', 6: 'soft_collinear', 14: 'virtual', 15: 'average_Born'}
NUMBER = r'[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[EeDd][-+]?\d+)?'


def number(value):
    return float(value.replace('D', 'E').replace('d', 'e'))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    return value


def accurate_sum(values):
    try:
        return math.fsum(values)
    except (ValueError, OverflowError):
        return float('nan')


def conditioning(values):
    net = accurate_sum(values)
    absolute = accurate_sum(abs(v) for v in values)
    return dict(sum=net, sum_absolute=absolute,
                cancellation_ratio=absolute / abs(net) if net else None,
                exact_zero_with_nonzero_terms=(net == 0 and absolute > 0))


def point_diagnostics(point, production_component=2):
    lines = point['lines']
    leaves = point.get('leaves', [])
    values = [r['central'] for r in lines] if lines else [r['central'] for r in leaves]
    result = conditioning(values)
    result['decomposition'] = 'additive_weight_lines' if lines else 'multiplicative_leaves' if leaves else 'absent'
    result['native_sum'] = point['net']
    result['native_vs_fsum_absolute'] = abs(point['net'] - result['sum'])
    result['native_vs_fsum_relative'] = (result['native_vs_fsum_absolute'] / abs(result['sum'])
                                       if result['sum'] else None)
    result['native_vs_fsum_scaled_by_sum_absolute'] = (
        result['native_vs_fsum_absolute'] / result['sum_absolute'] if result['sum_absolute'] else None)
    result['reported_sum_absolute'] = point['sumabs']
    result['component_type'] = []
    groups = {}
    for row in lines:
        groups.setdefault((row['component'], row['type']), []).append(row['central'])
    for (component, kind), members in sorted(groups.items()):
        result['component_type'].append(dict(component=component, type=kind,
            type_name=TYPE_NAMES.get(kind, str(kind)), **conditioning(members)))
    components = {}
    for row in lines:
        components.setdefault(row['component'], []).append(row['central'])
    result['components'] = [dict(component=c, **conditioning(v)) for c, v in sorted(components.items())]
    if result['components']:
        result['dominant_component'] = max(result['components'], key=lambda r: abs(r['sum']))['component']
    production = [r for r in lines if r['component'] == production_component]
    result['production'] = conditioning([r['central'] for r in production])
    result['production_local_subtraction'] = conditioning(
        [r['central'] for r in production if r['type'] in (1, 4, 5, 6)])
    result['largest_absolute_line'] = max(lines, key=lambda r: abs(r['central'])) if lines else None
    return result


def parse_trace(path, production_component=2):
    """Parse complete POINT/END records; preserve optional summary/moment rows."""
    points, summaries, moments, nonfinite, extras, current = [], [], [], [], [], None
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        fields = raw.split()
        if not fields or fields[0].startswith('#'):
            continue
        tag, values = fields[0], fields[1:]
        try:
            if tag == 'POINT':
                if len(values) < 6:
                    raise ValueError('POINT requires ordinal,channel,config,ndim,net,sumabs')
                current = dict(ordinal=int(values[0]), channel=int(values[1]), config=int(values[2]),
                               ndim=int(values[3]), net=number(values[4]), sumabs=number(values[5]), lines=[], leaves=[])
            elif tag == 'VOLUME' and current is not None:
                current['volume'] = number(values[0])
            elif tag == 'X' and current is not None:
                current['x'] = [number(v) for v in values]
                if len(current['x']) != current['ndim']:
                    raise ValueError('X length does not match ndim')
            elif tag == 'LINE' and current is not None:
                if len(values) < 12:
                    raise ValueError('LINE requires four integers and eight floats')
                current['lines'].append(dict(index=int(values[0]), type=int(values[1]),
                    component=int(values[2]), fks=int(values[3]), central=number(values[4]),
                    raw=[number(v) for v in values[5:8]], bjorken_x=[number(v) for v in values[8:10]],
                    scale_r=number(values[10]), scale_f=number(values[11]),
                    extra=[number(v) for v in values[12:]]))
            elif tag == 'LEAF' and current is not None:
                current['leaves'].append(dict(mask=int(values[0]), central=number(values[1])))
            elif tag == 'PRODUCTION' and current is not None:
                if len(values) < 7:
                    raise ValueError('PRODUCTION requires ifks,invprob,xi,1-y,xiHat,xiMax,shat')
                current['production_coordinates'] = dict(fks=int(values[0]), inverse_probability=number(values[1]),
                    xi=number(values[2]), one_minus_y=number(values[3]), xi_hat=number(values[4]),
                    xi_max=number(values[5]), shat=number(values[6]))
            elif tag == 'END' and current is not None:
                current['diagnostics'] = point_diagnostics(current, production_component)
                points.append(current)
                current = None
            elif tag == 'SUMMARY':
                if len(values) >= 7:
                    row = dict(points=int(values[0]), sum=number(values[1]), sum_squares=number(values[2]),
                               sum_absolute=number(values[3]), positive=int(values[4]), negative=int(values[5]),
                               nonfinite=int(values[6]))
                    if len(values) > 7:
                        row['extra'] = values[7:]
                    summaries.append(row)
                else:
                    extras.append(dict(tag=tag, values=values))
            elif tag in ('MOMENT', 'TYPEMOMENT'):
                offset = 1 if tag == 'TYPEMOMENT' else 0
                if len(values) >= 5 + offset:
                    row = dict(component=int(values[0]), points=int(values[1 + offset]),
                               sum=number(values[2 + offset]), sum_squares=number(values[3 + offset]),
                               sum_absolute=number(values[4 + offset]))
                    if offset:
                        row['type'] = int(values[1])
                    moments.append(row)
                else:
                    extras.append(dict(tag=tag, values=values))
            elif tag == 'NONFINITE':
                nonfinite.append(values)
            else:
                extras.append(dict(tag=tag, values=values))
        except (ValueError, OverflowError) as error:
            raise ValueError(f'{path}:{lineno}: {error}') from error
    return dict(points=points, summaries=summaries, component_moments=moments,
                nonfinite_records=nonfinite, extra_records=extras,
                incomplete_trailing_point=current is not None)


def native_log(directory):
    result = {}
    candidates = [directory / name for name in ('log.txt', 'logfile', 'stdout.log', 'worker.log', 'output.log')]
    for path in candidates:
        if not path.is_file():
            continue
        content = path.read_text(errors='replace')
        final = re.findall(r'Final result:\s*(' + NUMBER + r')\s*\+/-\s*(' + NUMBER + r')', content)
        iterations = re.findall(r'Maximum number of iterations is:\s*(\d+)', content)
        virtual = re.findall(r'^accumulated results V\s+1\s*=\s*(' + NUMBER + r')\s*\+/-\s*(' + NUMBER + r')',
                             content, re.MULTILINE)
        timing = re.findall(r'Time spent in ([^:\n]+)\s*:\s*(' + NUMBER + r')', content)
        if final:
            result.update(value_pb=number(final[-1][0]), error_pb=number(final[-1][1]))
        if iterations:
            result['maximum_iterations'] = int(iterations[-1])
        if virtual:
            result['production_virtual'] = dict(value_pb=number(virtual[-1][0]), error_pb=number(virtual[-1][1]),
                limitation='Production virtual residual estimate V 1; values/errors are rounded in the native text log.')
        if timing:
            result['cpu_seconds'] = {name.strip(): number(value) for name, value in timing}
        if final or iterations:
            result.update(path=str(path), sha256=sha256(path))
            break
    return result


def worker_report(path, root, production_component=2):
    parsed = parse_trace(path, production_component)
    points = parsed['points']
    metadata = {}
    for name in ('manifest.json', 'worker.json', 'metadata.json'):
        file = path.parent / name
        if file.is_file():
            metadata[name] = json.loads(file.read_text())
    row = dict(worker=str(path.parent.relative_to(root)), trace=str(path), sha256=sha256(path),
               metadata=metadata, native=native_log(path.parent), **parsed)
    nonfinite_path = path.parent / 'nonfinite.dat'
    if nonfinite_path.is_file():
        row['raw_nonfinite'] = dict(path=str(nonfinite_path), sha256=sha256(nonfinite_path),
            records=[line.split()[1:] for line in nonfinite_path.read_text().splitlines()
                     if line.split() and line.split()[0] == 'NONFINITE'],
            limitation='Raw nonfinite log is capped; retained entry count is a lower bound if the cap is reached.')
    finite = [p for p in points if math.isfinite(p['net'])]
    row['saved_record_count'] = len(points)
    row['largest_saved_absolute_weight'] = max((abs(p['net']) for p in finite), default=None)
    row['maximum_saved_sum_absolute'] = max((p['diagnostics']['sum_absolute'] for p in points), default=None)
    row['maximum_saved_native_vs_fsum_absolute'] = max(
        (p['diagnostics']['native_vs_fsum_absolute'] for p in points), default=None)
    row['maximum_saved_native_vs_fsum_relative'] = max(
        (p['diagnostics']['native_vs_fsum_relative'] for p in points
         if p['diagnostics']['native_vs_fsum_relative'] is not None), default=None)
    if parsed['summaries']:
        summary = parsed['summaries'][-1]
        n = summary['points']
        if n > 0:
            mean = summary['sum'] / n
            central_moment = max(0., summary['sum_squares'] - summary['sum'] * mean)
            row['all_point_moments'] = dict(mean=mean, centered_sum_squares=central_moment,
                independent_point_standard_error=math.sqrt(central_moment / (n * (n - 1))) if n > 1 else None,
                summary=summary)
            if finite:
                largest = max(finite, key=lambda p: abs(p['net']))
                row['largest_saved_raw_second_moment_fraction'] = (
                    largest['net'] ** 2 / summary['sum_squares'] if summary['sum_squares'] else None)
                row['largest_saved_centered_second_moment_fraction'] = (
                    (largest['net'] - mean) ** 2 / central_moment if central_moment else None)
                unique = {p['ordinal']: p for p in finite}
                row['all_saved_raw_second_moment_fraction_lower_bound'] = (
                    accurate_sum(p['net'] ** 2 for p in unique.values()) / summary['sum_squares']
                    if summary['sum_squares'] else None)
    return row


def replay_comparison(first, second):
    """Compare retained ordinals only; full stream equality needs all-point hashes."""
    left = {p['ordinal']: p for p in first['points']}
    right = {p['ordinal']: p for p in second['points']}
    common = sorted(left.keys() & right.keys())
    mismatches, differences = [], []
    for ordinal in common:
        a, b = left[ordinal], right[ordinal]
        # Omit derived diagnostics; compare the saved raw data.
        comparable = ('channel', 'config', 'ndim', 'net', 'sumabs', 'volume', 'x', 'lines', 'leaves', 'production_coordinates')
        changed = [key for key in comparable if a.get(key) != b.get(key)]
        if changed:
            mismatches.append(dict(ordinal=ordinal, changed=changed))
        differences.append(abs(a['net'] - b['net']))
    return dict(first=first['worker'], second=second['worker'], common_saved_ordinals=len(common),
        first_only_saved_ordinals=sorted(left.keys() - right.keys()),
        second_only_saved_ordinals=sorted(right.keys() - left.keys()),
        mismatches=mismatches, saved_records_equal=bool(common) and not mismatches
        and left.keys() == right.keys(), maximum_absolute_net_difference=max(differences, default=None),
        summaries_equal=bool(first['summaries']) and first['summaries'] == second['summaries'],
        limitation='Only retained records and available summary moments are compared; this is not a higher precision test.')


def independent_ensembles(rows):
    """Compare native worker errors with independent-worker scatter by allocation."""
    groups = {}
    for row in rows:
        job = row['metadata'].get('campaign_job', {})
        label = job.get('execution_label', '')
        role = ('smoke' if 'smoke' in label.lower() else 'replay' if job.get('same_stream_replay')
                else 'training' if any(token in label.lower() for token in ('train', 'warmup'))
                else 'independent' if job else 'unknown')
        row['role'] = role
        native = row['native']
        if role != 'independent' or 'value_pb' not in native or 'error_pb' not in native:
            continue
        key = (label, job.get('variant'), job.get('points'), job.get('iterations', 1))
        groups.setdefault(key, []).append(row)
    result = []
    for (label, variant, points, iterations), members in sorted(groups.items(), key=lambda item: str(item[0])):
        n = len(members)
        values = [r['native']['value_pb'] for r in members]
        mean = accurate_sum(values) / n
        errors = [r['native']['error_pb'] for r in members]
        native_error = math.sqrt(accurate_sum(e * e for e in errors)) / n
        scatter_error = math.sqrt(accurate_sum((v - mean) ** 2 for v in values) / (n * (n - 1))) if n > 1 else None
        variance = accurate_sum(e * e for e in errors)
        ensemble = dict(execution_label=label, variant=variant, requested_points=points,
            requested_iterations=iterations,
            workers=n, seeds=[r['metadata']['campaign_job'].get('seed') for r in members],
            native_mean_pb=mean, mean_native_error_pb=native_error,
            independent_scatter_error_pb=scatter_error,
            scatter_to_native_error_ratio=scatter_error / native_error if scatter_error is not None and native_error else None,
            largest_worker_native_variance_fraction=max((e * e for e in errors), default=0.) / variance if variance else None,
            worker_values_pb=values, worker_errors_pb=errors,
            total_elapsed_s=accurate_sum(r['metadata']['campaign_job'].get('elapsed_s', 0.) for r in members),
            limitations='The scatter and variance estimates have limited precision for a small number of independent workers and heavy tails; replay, smoke and training workers are excluded.')
        cpu = [r['native'].get('cpu_seconds', {}).get('Total') for r in members]
        if all(value is not None for value in cpu):
            ensemble['total_native_cpu_s'] = accurate_sum(cpu)
            ensemble['native_variance_times_cpu_pb2_s'] = native_error ** 2 * ensemble['total_native_cpu_s']
        virtual = [r['native'].get('production_virtual') for r in members]
        if all(value is not None for value in virtual):
            virtual_error = math.sqrt(accurate_sum(v['error_pb'] ** 2 for v in virtual)) / n
            ensemble['production_virtual_mean_pb'] = accurate_sum(v['value_pb'] for v in virtual) / n
            ensemble['production_virtual_error_pb'] = virtual_error
            ensemble['production_virtual_worker_errors_pb'] = [v['error_pb'] for v in virtual]
            if 'total_native_cpu_s' in ensemble:
                ensemble['production_virtual_variance_times_cpu_pb2_s'] = virtual_error ** 2 * ensemble['total_native_cpu_s']
        result.append(ensemble)
    return result


def make_report(root, production_component=2, replay_pairs=(), campaign_record=None):
    root = Path(root).resolve()
    rows = [worker_report(path, root, production_component) for path in sorted(root.rglob('record_breakers.dat'))]
    index = {r['worker']: r for r in rows}
    campaign = json.loads(Path(campaign_record).read_text()) if campaign_record else None
    if campaign:
        by_path = {str(Path(r['trace']).parent.resolve()): r for r in rows}
        earlier = {}
        for execution in campaign.get('executions', []):
            for job in execution.get('jobs', []):
                row = by_path.get(str(Path(job['directory']).resolve()))
                if row is None:
                    continue
                row['metadata']['campaign_job'] = dict(job, execution_label=execution.get('label'),
                    same_stream_replay=execution.get('same_stream_replay', False))
                key = (job.get('variant'), job.get('seed'))
                if execution.get('same_stream_replay') and key in earlier:
                    row['metadata']['campaign_job']['replay_of'] = earlier[key]
                else:
                    earlier.setdefault(key, row['worker'])
    comparisons = []
    for first, second in replay_pairs:
        if first not in index or second not in index:
            raise ValueError(f'Replay worker absent: {first}, {second}')
        comparisons.append(replay_comparison(index[first], index[second]))
    for row in rows:
        for metadata in row['metadata'].values():
            target = metadata.get('replay_of') if isinstance(metadata, dict) else None
            if target in index and not any(c['first'] == target and c['second'] == row['worker'] for c in comparisons):
                comparisons.append(replay_comparison(index[target], row))
    manifests = {}
    for name in ('manifest.json', 'pilot_manifest.json', 'launch.json'):
        path = root / name
        if path.is_file():
            manifests[name] = json.loads(path.read_text())
    ensembles = independent_ensembles(rows)
    return clean(dict(format=1, pilot_root=str(root), production_component=production_component,
        workers=rows, worker_count=len(rows), pilot_metadata=manifests,
        independent_ensembles=ensembles,
        campaign_record=str(campaign_record) if campaign_record else None,
        campaign_record_sha256=sha256(Path(campaign_record)) if campaign_record else None,
        replays=comparisons, script_sha256=sha256(Path(__file__)),
        limitations=[
            'Traces retain record breakers, possibly capped; they are not an unbiased sample of the tail distribution.',
            'The largest saved weight is a lower bound on the global maximum unless an all-point maximum is separately supplied.',
            'Saved-weight second-moment shares are lower bounds on contributions from all tail points; all-point summaries are needed.',
            'All-point moments can span adaptive iterations; their simple standard error need not equal the final native MINT estimate.',
            'Cancellation ratios quantify conditioning, not the accuracy of the individual matrix-element/subtraction terms.',
            'Repeatability and accurate final summation do not establish accuracy of double-precision amplitudes.',
            'Nonfinite JSON quantities are represented as null; zero-denominator cancellation is separately flagged.'
        ]))


def plot_report(report, output, execution_labels=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if output.exists():
        raise ValueError(f'Refuse to overwrite plot directory: {output}')
    output.mkdir(parents=True)
    rows = [row for row in report['workers'] if row.get('role') == 'independent'
            and row.get('native', {}).get('value_pb')]
    available_labels = {row['metadata'].get('campaign_job', {}).get('execution_label') for row in rows}
    if execution_labels is None and len(available_labels) > 1 and 'initial' in available_labels:
        execution_labels = ['initial']
    if execution_labels:
        rows = [row for row in rows if row['metadata'].get('campaign_job', {}).get('execution_label') in execution_labels]
    if not rows:
        raise ValueError('No completed independent workers available for diagnostic figure')
    with plt.rc_context({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
                         'legend.fontsize': 9, 'pdf.fonttype': 42}):
        fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.7))
        fig.subplots_adjust(left=.10, right=.98, bottom=.12, top=.91, hspace=.38, wspace=.30)
        max_replica = 0
        for variant, color, marker, shift in (('S', '#0072B2', 'o', -.18), ('Pi', '#D55E00', '^', .18)):
            peers = [row for row in rows if row['metadata'].get('campaign_job', {}).get('variant') == variant]
            peers.sort(key=lambda row: row['metadata']['campaign_job'].get('seed', 0))
            replicas = [int(row['metadata']['campaign_job'].get('replica', i)) + 1 for i, row in enumerate(peers)]
            max_replica = max([max_replica] + replicas)
            relative_errors = [100 * row['native']['error_pb'] / abs(row['native']['value_pb']) for row in peers]
            axes[0, 0].plot(replicas, relative_errors, marker=marker, color=color, lw=1.0,
                            ms=5, label=r'$\Pi$' if variant == 'Pi' else 'S')
            available = [(i, row['largest_saved_centered_second_moment_fraction']) for i, row in zip(replicas, peers)
                         if row.get('largest_saved_centered_second_moment_fraction') is not None]
            axes[0, 1].bar([i + shift for i, _ in available], [100 * v for _, v in available],
                           width=.34, color=color, alpha=.85)
            for row in peers:
                saved = row['points']
                mean = abs(row['native']['value_pb'])
                axes[1, 0].plot([p['ordinal'] for p in saved], [abs(p['net']) / mean for p in saved],
                                marker=marker, color=color, alpha=.55, ms=2.5, lw=.7)
                selected = [p for p in saved if p.get('production_coordinates', {}).get('xi', 0) > 0
                            and p['diagnostics']['cancellation_ratio']]
                axes[1, 1].scatter([p['production_coordinates']['xi'] for p in selected],
                                   [p['diagnostics']['cancellation_ratio'] for p in selected],
                                   marker=marker, color=color, alpha=.6, s=14, linewidths=0)
        axes[0, 0].set(xlabel='Independent worker index', ylabel='Native relative integration error [%]',
                       title='A  Quoted error at fixed sample count')
        axes[0, 0].legend(frameon=False)
        axes[0, 1].set(xlabel='Independent worker index', ylabel='Largest saved point / total variance [%]',
                       ylim=(0, 100), title='B  Concentration of estimated variance')
        for axis in axes[0]:
            axis.set_xticks(range(1, max_replica + 1))
        axes[1, 0].set(xlabel='Evaluation ordinal', ylabel=r'Saved $|w|$ / worker mean',
                       xscale='log', yscale='log', title='C  Growth of record net weights')
        axes[1, 1].set(xlabel=r'Production radiation $\xi$',
                       ylabel=r'$\sum_i |w_i|/|\sum_i w_i|$',
                       xscale='log', yscale='log', title='D  Cancellation in saved points')
        for axis in axes.flat:
            axis.grid(alpha=.15, which='major')
            axis.set_axisbelow(True)
        fig.suptitle(r'Bounded $t\bar t W^+$ diagnostic: central weights, $m_b=0$, on-shell Ws', fontsize=12)
        allocation_text = ', '.join(sorted({row['metadata']['campaign_job']['execution_label'] for row in rows}))
        fig.text(.10, .035, 'Allocation: '+allocation_text+'. Smoke, training and replay jobs excluded; saved points are record breakers, not a tail histogram.\n'
                 'Panel D sums additive weight lines for S and product leaves for Pi. Cancellation is conditioning, not an accuracy measurement.',
                 fontsize=8, va='bottom')
        fig.savefig(output / 'tail_diagnostics.pdf')
        fig.savefig(output / 'tail_diagnostics.png', dpi=180)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pilot_root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot-dir', type=Path)
    parser.add_argument('--plot-execution-label', action='append',
                        help='Allocation to include in plots (repeatable); default is initial if several are available')
    parser.add_argument('--production-component', type=int, default=2)
    parser.add_argument('--campaign-record', type=Path,
                        help='Optional runner JSON used to attach worker seeds, variants and explicit replay labels')
    parser.add_argument('--replay-pair', nargs=2, action='append', default=[], metavar=('WORKER', 'REPLAY'))
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError(f'Refuse to overwrite report: {args.output}')
    report = make_report(args.pilot_root, args.production_component, args.replay_pair, args.campaign_record)
    if not report['worker_count']:
        raise ValueError('No record_breakers.dat files found')
    if args.plot_dir:
        plot_report(report, args.plot_dir, args.plot_execution_label)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(f'Saved {report["worker_count"]} worker reports to {args.output}')


if __name__ == '__main__':
    main()
