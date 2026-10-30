#!/usr/bin/env python3
"""Figures and LaTeX tables for the current on-shell ttW scale comparison."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np

from collect_decay_scale_report import DEFAULT_OUTPUT, digest, save_json, paired_ratio_error

COLORS = {'S': '#0072B2', 'Pi': '#D55E00'}
LABELS = {'S': r'$S$', 'Pi': r'$\Pi$'}
WIDTHS = {'auto': 'Scale-dependent widths', 'fixed': 'Fixed widths'}
RATE_NAMES = ['inclusive rate', 'fiducial baseline', 'fiducial 0b rate',
              'fiducial 1b rate', 'fiducial 2b rate']
RATE_LABELS = ['Inclusive', 'Fiducial baseline', 'Exactly 0 b jets',
               'Exactly 1 b jet', 'Exactly 2 b jets']
SPECTRA = [
    ('lepton_pt', [
        (2, 'e+ pT', 12., 240., r'$p_T(e^+)$ [GeV]'),
        (2, 'e- pT', 12., 240., r'$p_T(e^-)$ [GeV]')]),
    ('muon_eta', [
        (2, 'mu pT', 12., 240., r'$p_T(\mu^+)$ [GeV]'),
        (2, 'e+ eta', -2.6, 2.6, r'$\eta(e^+)$')]),
    ('bottom_pt', [
        (2, 'b jet pT', 20., 320., r'$p_T(j_b)$ [GeV]'),
        (2, 'bbar jet pT', 20., 320., r'$p_T(j_{\bar b})$ [GeV]')]),
    ('lepton_bottom', [
        (2, 'e+ b mass', 0., 300., r'$m(e^+,j_b)$ [GeV]'),
        (2, 'mu b delta R', .4, 4.4, r'$\Delta R(\mu^+,j_b)$')]),
    ('trilepton_angles', [
        (2, 'trilepton mass', 0., 600., r'$m(3\ell)$ [GeV]'),
        (2, 'same-sign delta phi', 0., np.pi, r'$|\Delta\phi_{\mathrm{SS}}|$')]),
    ('low_tag_categories', [
        (0, 'mu pT', 12., 180., r'$p_T(\mu^+)$ [GeV]'),
        (1, 'mu pT', 12., 180., r'$p_T(\mu^+)$ [GeV]')]),
]


class Data:
    def __init__(self, directory):
        self.directory = directory
        self.report = json.loads((directory/'results.json').read_text())
        if digest(directory/'collected.npz') != self.report['collected_sha256']:
            raise ValueError('Collected array archive changed')
        self.arrays = np.load(directory/'collected.npz')
        self.titles = list(self.arrays['titles'])
        self.offsets = self.arrays['offsets']
        self.edges = self.arrays['edges']
        self.points = [tuple(p) for p in self.arrays['points']]
        self.groups = {key: [self.points.index(tuple(p)) for p in values]
                       for key, values in self.report['scale_groups'].items()}
        self.figures = []
        self.selection = []
        self.replicated = 'replica_blocks' in self.report
        self.rate_names = RATE_NAMES[1:] if self.replicated else RATE_NAMES
        self.rate_labels = RATE_LABELS[1:] if self.replicated else RATE_LABELS

    def spectrum(self, key, title, low, high):
        h = self.titles.index(title)
        s = slice(self.offsets[h], self.offsets[h+1])
        edges = self.edges[s]
        # Include complete native bins only. No statistical covariance is lost
        # through unrecorded rebinning; the displayed ranges are predeclared.
        use = (edges[:, 0] >= low-1e-6) & (edges[:, 1] <= high+1e-6)
        out = {field: self.arrays[key+'__'+field][s][use]
               for field in ('central', 'error', 'scales')}
        out['edges'] = edges[use]
        if self.replicated:
            out['replica_central'] = self.arrays[key+'__replica_central'][:,s][:,use]
            out['empirical_error'] = self.arrays[key+'__empirical_error'][s][use]
        return out

    def envelope(self, values, group):
        return (np.min(values[..., self.groups[group]], axis=-1),
                np.max(values[..., self.groups[group]], axis=-1))

    def save(self, fig, name):
        fig_dir = self.directory/'figures'
        fig_dir.mkdir(exist_ok=True)
        for ext in ('pdf', 'png'):
            fig.savefig(fig_dir/(name+'.'+ext), dpi=175,
                        bbox_inches='tight', facecolor='white')
        plt.close(fig)
        self.figures.append(dict(name=name, pdf_sha256=digest(fig_dir/(name+'.pdf')),
                                 png_sha256=digest(fig_dir/(name+'.png'))))


def axes_style(ax):
    ax.grid(alpha=.18)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8)


def rates_figure(data):
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 4.8), sharey=True)
    for ax, width in zip(axes, WIDTHS):
        for variant, delta in [('S', -.13), ('Pi', .13)]:
            for i, name in enumerate(data.rate_names):
                row = data.report['rates'][variant+'_'+width][name]
                reference = data.report['rates']['S_'+width][name]['central_pb']
                y = i+delta
                broad, narrow = row['combined63'], row['correlated7']
                ax.plot(np.array([broad['low_pb'], broad['high_pb']])/reference,
                        [y, y], color=COLORS[variant], alpha=.25, lw=10, solid_capstyle='butt')
                ax.plot(np.array([narrow['low_pb'], narrow['high_pb']])/reference,
                        [y, y], color=COLORS[variant], lw=2.0)
                ax.errorbar(row['central_pb']/reference, y,
                    xerr=row['mc_error_pb']/reference, fmt='o', ms=3.5,
                    color=COLORS[variant], capsize=2, zorder=5)
        ax.axvline(1, color='0.3', lw=.8, ls=':')
        ax.set_title(WIDTHS[width], fontsize=10)
        ax.set_xlabel(r'Cross section / $S$ central', fontsize=9)
        ax.set_xlim(.77, 1.43)
        ax.set_yticks(range(len(data.rate_names)), data.rate_labels)
        axes_style(ax)
    axes[0].invert_yaxis()
    handles = [Line2D([], [], color=COLORS[v], marker='o', label=LABELS[v]) for v in COLORS]
    handles += [Line2D([], [], color='0.35', alpha=.25, lw=9, label='Independent: 63 points'),
                Line2D([], [], color='0.35', lw=2, label='Correlated: 7 points')]
    fig.legend(handles=handles, loc='lower center', ncol=2, fontsize=9,
               bbox_to_anchor=(.5, -.01), frameon=False)
    fig.suptitle(r'$pp\to t\bar t W^+\to b\bar b\,e^+e^-\mu^++3\nu$, 13 TeV', fontsize=11)
    fig.tight_layout(rect=(0, .14, 1, .94))
    data.save(fig, '01_rates')


def envelope_figure(data):
    groups = ['production7', 'decay3', 'decay9', 'common3',
              'correlated7', 'combined21', 'combined63']
    names = ['Prod.\n7', 'Decay\n3', 'Decay\n9', 'Common\n3',
             'Corr.\n7', 'Indep.\n21', 'Indep.\n63']
    fig, axes = plt.subplots(2, 2, figsize=(7.3, 6), sharex=True, sharey=True)
    for col, width in enumerate(WIDTHS):
        for ax, rate in zip(axes[:, col], ('inclusive rate', 'fiducial 2b rate')):
            for variant, delta in [('S', -.1), ('Pi', .1)]:
                row = data.report['rates'][variant+'_'+width][rate]
                low = [-row[g]['down_percent'] for g in groups]
                high = [row[g]['up_percent'] for g in groups]
                xx = np.arange(len(groups))+delta
                ax.vlines(xx, low, high, color=COLORS[variant], lw=4, alpha=.8)
                ax.scatter(xx, np.zeros(len(xx)), s=14, color=COLORS[variant], zorder=3)
            ax.axhline(0, color='0.3', lw=.7)
            ax.set_title(('Varying widths' if width=='auto' else 'Fixed widths')+
                         '; '+(('baseline' if data.replicated else 'inclusive')
                                 if rate=='inclusive rate' else '2b'), fontsize=10)
            ax.set_xticks(range(len(groups)), names, fontsize=8)
            ax.set_ylim(-16, 19)
            axes_style(ax)
    for ax in axes[:, 0]:
        ax.set_ylabel('Scale-envelope deviation [%]')
    fig.legend(handles=[Line2D([], [], color=COLORS[v], lw=4, label=LABELS[v]) for v in COLORS],
               loc='upper center', ncol=2, frameon=False, bbox_to_anchor=(.5, 1.01))
    fig.tight_layout(rect=(0, 0, 1, .96))
    data.save(fig, '02_scale_prescriptions')


def response_figure(data):
    fig, axes = plt.subplots(2, 2, figsize=(7.3, 6.8), sharex=True, sharey=True)
    factors = (.5, 1., 2.)
    matrices = {}
    for row, variant in enumerate(COLORS):
        for col, width in enumerate(WIDTHS):
            record = data.report['rates'][variant+'_'+width]['inclusive rate']
            points = record['scale_points']
            matrix = np.array([[points[','.join(map(str, (p, p, d, d)))]
                                for p in factors] for d in factors])
            matrix = 100*(matrix/record['central_pb']-1)
            matrices[variant+'_'+width] = matrix
            ax = axes[row, col]
            im = ax.imshow(matrix, origin='lower', cmap='RdBu_r', vmin=-16, vmax=16)
            for i, j in np.ndindex(matrix.shape):
                ax.text(j, i, f'{matrix[i,j]:+.2f}%', ha='center', va='center',
                        fontsize=11, color='white' if abs(matrix[i,j])>10 else 'black')
            ax.set_xticks(range(3), ['1/2', '1', '2'])
            ax.set_yticks(range(3), ['1/2', '1', '2'])
            ax.set_title(LABELS[variant]+': '+WIDTHS[width], fontsize=11)
            if row == 1:
                ax.set_xlabel(r'Production factor $\xi_R=\xi_F$')
            if col == 0:
                ax.set_ylabel(r'Decay factor $\xi_t=\xi_{\bar t}$')
    fig.subplots_adjust(left=.09, right=.87, bottom=.08, top=.91, hspace=.28, wspace=.24)
    cax = fig.add_axes([.9, .18, .02, .58])
    label = 'Baseline-selected' if data.replicated else 'Inclusive'
    fig.colorbar(im, cax=cax, label=label+' rate deviation [%]')
    fig.suptitle(label+' production / decay scale response', fontsize=13)
    data.save(fig, '03_scale_response')
    data.report['inclusive_response_matrices_percent'] = matrices


def step_band(ax, edges, low, high, color, alpha):
    x = np.r_[edges[:, 0], edges[-1, 1]]
    ax.fill_between(x, np.r_[low, low[-1]], np.r_[high, high[-1]],
                    step='post', color=color, alpha=alpha, lw=0)


def point_matched_ratio(numerator, denominator):
    # The same scale vector is used in numerator and denominator. Unknown
    # MC covariance does not enter deterministic scale envelopes.
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan),
                     where=denominator>0)


def spectrum_figure(data, filename, specs):
    fig = plt.figure(figsize=(7.3, 7.8))
    grid = fig.add_gridspec(2, 2, hspace=.5, wspace=.3)
    for row, (category, suffix, low, high, xlabel) in enumerate(specs):
        title = f'fiducial {category}b {suffix}'
        ratio_min, ratio_max = 1., 1.
        for width in WIDTHS:
            p = data.spectrum('Pi_'+width, title, low, high)
            s = data.spectrum('S_'+width, title, low, high)
            a, b = data.envelope(point_matched_ratio(p['scales'], s['scales']), 'combined63')
            ratio_min = min(ratio_min, float(np.nanmin(a)))
            ratio_max = max(ratio_max, float(np.nanmax(b)))
            if data.replicated:
                mc = paired_ratio_error(p['replica_central'], s['replica_central'])
                c = point_matched_ratio(p['central'], s['central'])
                ratio_min = min(ratio_min, float(np.nanmin(c-mc)))
                ratio_max = max(ratio_max, float(np.nanmax(c+mc)))
        margin = max(.04, .08*(ratio_max-ratio_min))
        for col, width in enumerate(WIDTHS):
            cell = grid[row, col].subgridspec(2, 1, height_ratios=[2.4, 1], hspace=.06)
            upper = fig.add_subplot(cell[0])
            lower = fig.add_subplot(cell[1], sharex=upper)
            snapshots = {}
            for variant in COLORS:
                key = variant+'_'+width
                d = data.spectrum(key, title, low, high)
                snapshots[variant] = d
                edges = d['edges']
                center = edges.mean(axis=1)
                density = 1000/(edges[:, 1]-edges[:, 0])
                for group, alpha in [('combined63', .16), ('correlated7', .27)]:
                    a, b = data.envelope(d['scales'], group)
                    step_band(upper, edges, a*density, b*density, COLORS[variant], alpha)
                upper.stairs(d['central']*density, np.r_[edges[:, 0], edges[-1, 1]],
                             color=COLORS[variant], lw=1.25)
                upper.errorbar(center, d['central']*density, yerr=d['error']*density,
                               fmt='.', ms=2, color=COLORS[variant], lw=.65, capsize=0)
                data.selection.append(dict(run=key, title=title, low=low, high=high,
                    native_bins=len(edges), nonpositive_bins=int((d['central']<=0).sum()),
                    median_relative_mc_error=float(np.median(
                        d['error'][d['central']>0]/d['central'][d['central']>0]))))
            numerator, denominator = snapshots['Pi'], snapshots['S']
            ratio = point_matched_ratio(numerator['scales'], denominator['scales'])
            for group, alpha in [('combined63', .2), ('correlated7', .4)]:
                a, b = data.envelope(ratio, group)
                step_band(lower, edges, a, b, COLORS['Pi'], alpha)
            central_ratio = point_matched_ratio(numerator['central'], denominator['central'])
            lower.plot(center, central_ratio, '.-', ms=2.5, lw=.8, color=COLORS['Pi'])
            if data.replicated:
                ratio_mc = paired_ratio_error(numerator['replica_central'],
                                              denominator['replica_central'])
                lower.errorbar(center, central_ratio, yerr=ratio_mc, fmt='none',
                               color=COLORS['Pi'], lw=.65, capsize=1.2)
            lower.axhline(1, color='0.3', lw=.7, ls=':')
            lower.set_ylim(min(.9, ratio_min-margin), max(1.1, ratio_max+margin))
            lower.yaxis.set_major_locator(MaxNLocator(4))
            lower.set_ylabel(r'$\Pi/S$', fontsize=9)
            lower.set_xlabel(xlabel, fontsize=10)
            upper.set_ylabel('fb / GeV' if '[GeV]' in xlabel else 'fb / unit', fontsize=9)
            upper.set_title(f'{category}b; '+('Varying widths' if width=='auto' else 'Fixed widths'), fontsize=10)
            upper.tick_params(labelbottom=False)
            upper.set_xlim(edges[0, 0], edges[-1, 1])
            # Linear axes show signed noisy bins instead of silently dropping them.
            upper.axhline(0, color='0.3', lw=.5)
            axes_style(upper)
            axes_style(lower)
    fig.subplots_adjust(left=.1, right=.985, top=.92, bottom=.065)
    fig.legend(handles=[Line2D([], [], color=COLORS[v], lw=1.5, label=LABELS[v]) for v in COLORS]+
        [Line2D([], [], color='0.4', lw=8, alpha=.2, label='I63 (light)'),
         Line2D([], [], color='0.4', lw=8, alpha=.45, label='C7 (dark)')],
        loc='upper center', ncol=4, fontsize=9, frameon=False, bbox_to_anchor=(.5, .985))
    data.save(fig, filename)


def replica_figure(data):
    fig, axes = plt.subplots(2, 2, figsize=(7.3, 5.8), sharex=True)
    count = data.report['replica_blocks']['count']
    xx = np.arange(1, count+1)
    for row, variant in enumerate(COLORS):
        for col, width in enumerate(WIDTHS):
            record = data.report['rates'][variant+'_'+width]['fiducial baseline']
            reference = record['central_pb']
            values = 100*(np.asarray(record['replica_central_pb'])/reference-1)
            errors = 100*np.asarray(record['replica_native_error_pb'])/reference
            native = 100*record['mc_error_pb']/reference
            empirical = 100*record['empirical_mc_error_pb']/reference
            ax = axes[row,col]
            ax.axhspan(-native, native, color='0.65', alpha=.3)
            ax.axhline(empirical, color='0.3', ls='--', lw=.8)
            ax.axhline(-empirical, color='0.3', ls='--', lw=.8)
            ax.axhline(0, color='0.4', lw=.5)
            ax.errorbar(xx, values, yerr=errors, fmt='o', color=COLORS[variant],
                        ms=3, capsize=2, lw=.8)
            ax.set_title(LABELS[variant]+': '+WIDTHS[width], fontsize=10)
            ax.set_xticks(xx)
            axes_style(ax)
            if row==1: ax.set_xlabel('Independent replica block')
            if col==0: ax.set_ylabel('Baseline deviation from mean [%]')
    fig.legend(handles=[Line2D([],[],color='0.65',alpha=.4,lw=8,label='Mean: propagated native MC error'),
                        Line2D([],[],color='0.3',ls='--',label='Mean: empirical replica SE')],
               loc='upper center', ncol=2, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0,0,1,.94))
    data.save(fig,'10_replica_precision')


def scale_tex(row, group):
    entry = row[group]
    return r'$^{+%.2f}_{-%.2f}\%%$' % (entry['up_percent'], entry['down_percent'])


def write_tables(data):
    out = data.directory
    columns = 'Baseline & $0b$ & $1b$ & $2b$' if data.replicated else 'Inclusive & Baseline & $0b$ & $1b$ & $2b$'
    text = [r'\begin{tabular}{ll'+len(data.rate_names)*'r'+'}', r'\toprule',
            'Width & Approx. & '+columns+r' \\', r'\midrule']
    for width in WIDTHS:
        for variant in COLORS:
            rows = data.report['rates'][variant+'_'+width]
            values = [r'$%.6f\pm%.6f$' % (1000*rows[n]['central_pb'], 1000*rows[n]['mc_error_pb'])
                      for n in data.rate_names]
            text.append(('Varying' if width=='auto' else 'Fixed')+' & '+LABELS[variant]+' & '+
                        ' & '.join(values)+r' \\')
    text += [r'\bottomrule', r'\end{tabular}']
    (out/'rates_table.tex').write_text('\n'.join(text)+'\n')
    text = [r'\begin{tabular}{lllrrrr}', r'\toprule',
            r'Width & Approx. & Selection & C7 & I21 & I63 & Full81 \\', r'\midrule']
    for width in WIDTHS:
        for variant in COLORS:
            for name, label in [('inclusive rate', 'Inclusive'), ('fiducial baseline', 'Baseline'),
                                ('fiducial 0b rate', '$0b$'), ('fiducial 1b rate', '$1b$'),
                                ('fiducial 2b rate', '$2b$')]:
                row = data.report['rates'][variant+'_'+width][name]
                if name not in data.rate_names:
                    continue
                text.append(('Varying' if width=='auto' else 'Fixed')+' & '+LABELS[variant]+' & '+label+' & '+
                    ' & '.join(scale_tex(row,g) for g in ('correlated7','combined21','combined63','all81_diagnostic'))+r' \\')
    text += [r'\bottomrule', r'\end{tabular}']
    (out/'scales_table.tex').write_text('\n'.join(text)+'\n')
    text = [r'\begin{tabular}{ll'+len(data.rate_names)*'r'+'}', r'\toprule',
            'Width & Quantity & '+columns+r' \\', r'\midrule']
    for width in WIDTHS:
        rows = data.report['rate_comparisons'][width]
        names = ['inclusive rate','fiducial cutflow','fiducial 0b rate','fiducial 1b rate','fiducial 2b rate']
        if data.replicated:
            names = names[1:]
        text.append(('Varying' if width=='auto' else 'Fixed')+r' & $100(\Pi/S-1)$ & '+
            ' & '.join('$%+.2f$'% (100*(rows[n]['central_ratio']-1)) for n in names)+r' \\')
        if data.replicated:
            text.append(r' & Paired MC error on $\Pi/S$ & '+
                ' & '.join('$%.5f$'% rows[n]['paired_mc_error'] for n in names)+r' \\')
        for group,label in [('correlated7','C7'),('combined63','I63')]:
            text.append(' & '+label+' ratio range & '+' & '.join(
                r'$[%.4f,%.4f]$' % tuple(rows[n]['scale_envelopes'][group]) for n in names)+r' \\')
    text += [r'\bottomrule', r'\end{tabular}']
    (out/'ratios_table.tex').write_text('\n'.join(text)+'\n')
    text = [r'\begin{tabular}{rlrrrr'+('r' if data.replicated else '')+'}', r'\toprule',
            r'$m(e^+,j_b)$ [GeV] & Width & $\Pi/S$ & I63 ratio range & MC $S$ [\%] & MC $\Pi$ [\%]'+
            (r' & MC $\Pi/S$' if data.replicated else '')+r' \\',
            r'\midrule']
    for low, high in [(90.,120.), (120.,150.), (150.,180.)]:
        for width in WIDTHS:
            s = data.spectrum('S_'+width, 'fiducial 2b e+ b mass', low, high)
            p = data.spectrum('Pi_'+width, 'fiducial 2b e+ b mass', low, high)
            if len(s['central']) != 1:
                raise ValueError('Mass table requires the stated native bin')
            a,b = data.envelope(point_matched_ratio(p['scales'],s['scales']), 'combined63')
            line = r'$%g$--$%g$ & %s & $%.4f$ & $[%.4f,%.4f]$ & $%.2f$ & $%.2f$' % (
                low, high, 'Varying' if width=='auto' else 'Fixed',
                p['central'][0]/s['central'][0], a[0], b[0],
                100*s['error'][0]/s['central'][0], 100*p['error'][0]/p['central'][0])
            if data.replicated:
                line += ' & $%.5f$'% paired_ratio_error(p['replica_central'],s['replica_central'])[0]
            text.append(line+r' \\')
    text += [r'\bottomrule',r'\end{tabular}']
    (out/'mass_ratio_table.tex').write_text('\n'.join(text)+'\n')
    if data.replicated:
        text = [r'\begin{tabular}{lllrrr}',r'\toprule',
                r'Width & Approx. & Selection & Native MC [\%] & Replica SE [\%] & Replica/native \\',r'\midrule']
        for width in WIDTHS:
            for variant in COLORS:
                for name,label in zip(data.rate_names,['Baseline','$0b$','$1b$','$2b$']):
                    row = data.report['rates'][variant+'_'+width][name]
                    text.append(('%s & %s & %s & $%.3f$ & $%.3f$ & $%.2f$' % (
                        'Varying' if width=='auto' else 'Fixed',LABELS[variant],label,
                        100*row['relative_mc_error'],
                        100*row['empirical_mc_error_pb']/abs(row['central_pb']),row['empirical_over_native']))+r' \\')
        text += [r'\bottomrule',r'\end{tabular}']
        (out/'precision_table.tex').write_text('\n'.join(text)+'\n')
    with (out/'histogram_precision.csv').open('w') as stream:
        columns = list(data.report['precision_inventory'][0])
        writer = csv.DictWriter(stream, columns)
        writer.writeheader()
        writer.writerows(data.report['precision_inventory'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
        'axes.labelsize':9,'savefig.dpi':175,'pdf.fonttype':42})
    data = Data(args.output.resolve())
    rates_figure(data)
    envelope_figure(data)
    response_figure(data)
    for number, (name, specs) in enumerate(SPECTRA, 4):
        spectrum_figure(data, '%02d_%s'%(number,name), specs)
    if data.replicated:
        replica_figure(data)
    write_tables(data)
    save_json(data.directory/'figure_manifest.json', dict(figures=data.figures,
        selections=data.selection, renderer_sha256=digest(__file__),
        convention='All bins are native; absolute density bands and point-matched ratio bands. '
            'Light I63 and dark C7. '+
            ('Absolute bars propagate native MC errors across replicas; ratio bars are paired-replica jackknife SEs.'
             if data.replicated else 'Absolute bars are per-run MC errors; ratio panels carry scale envelopes only because cross-run covariance is unavailable.')))
    save_json(data.directory/'results.json', data.report)
    print('Wrote',len(data.figures),'figures and LaTeX tables to',data.directory)


if __name__=='__main__':
    main()
