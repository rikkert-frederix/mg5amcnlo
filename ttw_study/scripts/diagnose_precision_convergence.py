#!/usr/bin/env python3
"""Read native logs only: precision timing, refinement work and weight tails.

This does not run MG5, inspect old worker archives, reject samples or alter
campaign state. Log moments locate suspect workers, not individual events.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np

NUMBER = r'[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[EeDd][-+]?\d+)?'
ROWS = re.compile(r'<a name=(/[^<>\s]+)></a>.*?<PRE>\n(.*?)\n</PRE>', re.S)


def number(value):
    return float(value.replace('D', 'E').replace('d', 'e'))


def worker(directory, log):
    final = re.findall(r'Final result:\s*('+NUMBER+r')\s*\+/-\s*('+NUMBER+r')', log)
    if not final:
        return None
    value, error = map(number, final[-1])
    iterations = re.findall(r'Maximum number of iterations is:\s*(\d+)', log)
    points = re.findall(r'Number of phase-space points per iteration:\s*(\d+)', log)
    channels = [dict(channel=int(c), config=int(config), samples=int(n), accepted=int(a),
                     absolute=number(absolute), signed=number(signed))
                for c, config, n, a, absolute, signed in re.findall(
                    r'^channel\s+(\d+)\s*:\s*(\d+)\s+[TF]\s+(\d+)\s+(\d+)\s+('
                    +NUMBER+r')\s+('+NUMBER+r')', log, re.M)]
    virtual = re.findall(r'^accumulated results Virtual\s*=\s*('+NUMBER+r')\s*\+/-\s*('+NUMBER+r')', log, re.M)
    components = [dict(component=int(i), value=number(v), error=number(e))
                  for i, v, e in re.findall(r'^accumulated results C\s+(\d+)\s*=\s*('
                      +NUMBER+r')\s*\+/-\s*('+NUMBER+r')', log, re.M)]
    stability = {}
    for label in ('Stable PS point', 'Unstable PS point (and rescued)',
                  'Exceptional PS point (unstable and not rescued)', 'Quadruple precision used'):
        counts = re.findall(re.escape(label)+r':\s*(\d+)', log)
        stability[label] = int(counts[-1]) if counts else None
    cpu = re.findall(r'Time spent in ([^:\n]+):\s*('+NUMBER+r')', log)
    return dict(directory=directory, value_pb=value, error_pb=error,
                requested_points_per_iteration=int(points[-1]) if points else None,
                maximum_iterations=int(iterations[-1]) if iterations else None,
                sampled_channel_calls=sum(c['samples'] for c in channels),
                channels=channels, components=components, stability=stability,
                virtual_error_pb=number(virtual[-1][1]) if virtual else None,
                cpu_seconds={name.strip():number(t) for name, t in cpu})


def summarize(rows):
    variance = sum(r['error_pb']**2 for r in rows)
    ordered = sorted(rows, key=lambda r:r['error_pb'], reverse=True)
    groups = {}
    for r in rows:
        group = r['directory'].rsplit('_', 1)[0]
        groups.setdefault(group, []).append(r)
    result = dict(workers=len(rows), sampled_channel_calls=sum(r['sampled_channel_calls'] for r in rows),
                  sum_pb=sum(r['value_pb'] for r in rows), error_pb=variance**.5,
                  top_one_variance_fraction=ordered[0]['error_pb']**2/variance if variance else None,
                  top_four_variance_fraction=sum(r['error_pb']**2 for r in ordered[:4])/variance if variance else None,
                  worst_workers=[])
    result['strata'] = {name:dict(workers=len(peers),
        sampled_channel_calls=sum(p['sampled_channel_calls'] for p in peers),
        sum_pb=sum(p['value_pb'] for p in peers),
        error_pb=sum(p['error_pb']**2 for p in peers)**.5)
        for name, peers in groups.items()}
    for r in ordered[:4]:
        peers = groups[r['directory'].rsplit('_', 1)[0]]
        median = float(np.median([p['value_pb'] for p in peers]))
        result['worst_workers'].append(dict(r, variance_fraction=r['error_pb']**2/variance,
            peer_median_pb=median, displacement_in_own_errors=(r['value_pb']-median)/r['error_pb']))
    return result


def run_record(process, run, case, execution):
    event = process/'Events'/run
    history = []; last = {}; completed_calls = 0
    for path in sorted(event.glob('alllogs_*.html'), key=lambda p:int(p.stem.split('_')[-1])):
        stage = int(path.stem.split('_')[-1])
        rows = [w for directory, log in ROWS.findall(path.read_text()) if (w:=worker(directory, log))]
        if not rows:
            continue
        completed_calls += sum(r['sampled_channel_calls'] for r in rows)
        if stage:
            last.update({r['directory']:dict(r, stage=stage) for r in rows})
        history.append(dict(stage=stage, **summarize(rows),
            retained_results=summarize(list(last.values())) if last else None))
    e = json.loads(execution.read_text()) if execution.exists() else {}
    output = dict(case=case, process=str(process), run=run, execution_status=e.get('status'),
        elapsed_hours=e['elapsed_s']/3600 if 'elapsed_s' in e else None, history=history,
        completed_logged_channel_calls=completed_calls,
        limitations='Only complete logged stages are included; interrupted uncollected workers are omitted. '
        'Different refinements replace earlier samples, and may include different strata. '
        'Worker moments do not identify an individual point or distinguish physical tails from roundoff.')
    if last:
        output['last_logged_worker_results'] = summarize(list(last.values()))
    components = sorted(event.glob('contribution_results_*.txt'), key=lambda p:int(p.stem.split('_')[-1]))
    if components:
        output['final_component_summary'] = [dict(component=int(i), name=name, value_pb=number(v), error_pb=number(err))
            for i, name, v, err in re.findall(r'^COMPONENT (\d+) (\S+) ('+NUMBER+r') ('+NUMBER+r')', components[-1].read_text(), re.M)]
    return output


def plot(results, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if output.exists():
        raise ValueError('Refuse to overwrite plot directory')
    output.mkdir(parents=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for variant, marker in (('S', 'o'), ('Pi', '^')):
        rows = [r for r in results if r['case']['variant']==variant and r['execution_status']=='finished']
        x = [r['completed_logged_channel_calls']/1e6 for r in rows]
        y = [100*r['last_logged_worker_results']['error_pb']/abs(r['last_logged_worker_results']['sum_pb']) for r in rows]
        axes[0, 0].scatter(x, y, marker=marker, label=variant)
        for r, a, b in zip(rows, x, y):
            axes[0, 0].annotate(str(r['case']['seed']), (a, b), xytext=(3, 4), textcoords='offset points', fontsize=7)
    axes[0, 0].axhline(.5, color='grey', ls='--', label='Requested 0.5%')
    axes[0, 0].set(xscale='log', xlabel='Logged channel calls, all stages [millions]', ylabel='Final quoted relative error [%]')
    axes[0, 0].margins(x=.2, y=.15)
    axes[0, 0].legend(fontsize=8)
    latest = max((r for r in results if r['case']['variant']=='S' and r['execution_status']=='finished'),
                 key=lambda r:r['elapsed_hours'])
    label = 'S%d'%latest['case']['seed']
    rows = [s for s in latest['history'] if s['stage']]
    stages = [s['stage'] for s in rows]
    axes[0, 1].plot(stages, [1e3*s['retained_results']['error_pb'] for s in rows], 'o-', label='Retained workers')
    axes[0, 1].set(xticks=stages, xlabel=label+' refinement stage', ylabel='Quoted cross-section error [fb]')
    axes[0, 1].margins(x=.15, y=.2)
    for s in rows:
        axes[0, 1].annotate('%.1f M new calls'%(s['sampled_channel_calls']/1e6),
            (s['stage'], 1e3*s['retained_results']['error_pb']), xytext=(0, 8), textcoords='offset points', ha='center', fontsize=8)
    variance = latest['last_logged_worker_results']
    top = variance['worst_workers']
    fractions = [100*r['variance_fraction'] for r in top]
    labels = [('dxu' if '_dxu_' in r['directory'] else 'udx')+' / '+r['directory'].split('/')[-1] for r in top]
    axes[1, 0].bar(range(5), fractions+[100-sum(fractions)])
    axes[1, 0].set(xticks=range(5), xticklabels=labels+['Other workers'], ylabel='Share of final variance [%]', title=label+': final worker influence')
    axes[1, 0].tick_params(axis='x', rotation=20, labelsize=8)
    components = latest['final_component_summary']
    labels = ['Born', 'Production NLO', 'Top decay NLO', 'Antitop decay NLO', 'Width term']
    axes[1, 1].bar(range(5), [1e3*c['error_pb'] for c in components])
    axes[1, 1].set(xticks=range(5), xticklabels=labels, yscale='log', ylabel='Component quoted error [fb]', title=label+': component errors')
    axes[1, 1].tick_params(axis='x', rotation=20, labelsize=8)
    fig.suptitle('ttW on-shell precision campaign: sampling tails (no samples removed)', fontsize=12, y=1.02)
    fig.savefig(output/'convergence.pdf', bbox_inches='tight')
    fig.savefig(output/'convergence.png', dpi=160, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot-dir', type=Path)
    args = parser.parse_args()
    record = json.loads(args.campaign.read_text())
    jobs = list(record['jobs'])
    case = record.get('current_case')
    if case and case not in [j['case'] for j in jobs]:
        study = args.campaign.resolve().parents[1]
        process = study/'processes'/('TTW%s_onshell_%s_mb0p0_both_%s'%(
            case['charge'], ''.join(case['flavours']), record['tag']))
        run = '%s_onshell_core-w-ht-half_separate_%d'%(case['variant'], case['seed'])
        jobs.append(dict(case=case, process=str(process), run=run, execution=str(process/'study_cards'/run/'execution.json')))
    results = []
    for j in jobs:
        result = run_record(Path(j['process']), j['run'], j['case'], Path(j['execution']))
        results.append(result)
        stages = result['history']
        print(j['case']['seed'], j['case']['variant'], 'hours='+str(result['elapsed_hours']),
              'stages='+','.join(str(s['stage']) for s in stages), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise ValueError('Refuse to overwrite diagnostic')
    args.output.write_text(json.dumps(dict(campaign=str(args.campaign), status=record['status'], runs=results), indent=2)+'\n')
    if args.plot_dir:
        plot(results, args.plot_dir)


if __name__ == '__main__':
    main()
