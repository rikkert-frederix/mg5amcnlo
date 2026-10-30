#!/usr/bin/env python3
"""Plot completed bounded pilot arms from their compact numerical report."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    if args.output_dir.exists():raise ValueError('Refuse to overwrite figures')
    report=json.loads(args.report.read_text())
    names=('initial','virtual_floor_005','S_proposal')
    labels=('Baseline proposal','5% virtual sampling','S-trained proposal')
    rows=[next(e for e in report['independent_ensembles'] if e['execution_label']==name and e['variant']=='S') for name in names]
    x=np.arange(3);colors=('#0072b2','#e69f00','#009e73')
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    axes[0,0].bar(x-.18,[e['mean_native_error_pb']*1e6 for e in rows],.36,label='Native errors combined')
    axes[0,0].bar(x+.18,[e['independent_scatter_error_pb']*1e6 for e in rows],.36,label='Between-worker spread')
    axes[0,0].set(ylabel=r'Error of ensemble mean [$10^{-6}$ pb]',title='Eight independent S workers per arm')
    axes[0,0].legend(fontsize=8)
    for ax,key,title in ((axes[0,1],'native_variance_times_cpu_pb2_s','Full-integral variance × CPU'),
                         (axes[1,0],'production_virtual_variance_times_cpu_pb2_s','Virtual-residual variance × CPU')):
        values=[e[key]/rows[0][key] for e in rows]
        ax.bar(x,values,color=colors)
        ax.axhline(1,color='grey',ls='--',lw=1)
        ax.set(ylabel='Relative to baseline; lower is better',title=title)
        for i,v in enumerate(values):ax.text(i,v,'%.3g'%v,ha='center',va='bottom')
        ax.margins(y=.18)
    for i,name in enumerate(names):
        workers=[w for w in report['workers'] if w.get('role')=='independent'
                 and w['metadata']['campaign_job']['execution_label']==name
                 and w['metadata']['campaign_job']['variant']=='S']
        # Role is attached by the ensemble collector; tolerate report variants.
        if not workers:
            workers=[w for w in report['workers'] if w['metadata'].get('campaign_job',{}).get('execution_label')==name
                     and w['metadata']['campaign_job']['variant']=='S']
        values=[100*w['largest_saved_centered_second_moment_fraction'] for w in workers]
        axes[1,1].scatter(np.full(len(values),i)+np.linspace(-.12,.12,len(values)),values,color=colors[i])
    axes[1,1].set(ylabel='Largest point / worker variance [%]',ylim=(0,105),title='Concentration persists in every arm')
    for ax in axes.flat:
        ax.set_xticks(x,labels,rotation=12,ha='right',fontsize=9)
    fig.suptitle('Central-only ttW diagnostic: sampler controls',fontsize=14)
    fig.supxlabel('Each arm uses 8 × 100,000 points. CPU excludes proposal training.\n'
                  'Variance estimates fluctuate strongly; these are pilot comparisons, not performance guarantees.',fontsize=9)
    args.output_dir.mkdir(parents=True)
    fig.savefig(args.output_dir/'sampler_controls.pdf',bbox_inches='tight')
    fig.savefig(args.output_dir/'sampler_controls.png',dpi=160,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':main()
