#!/usr/bin/env python3
"""Nominal W-treatment diagnostic; conditional single-flavour pilots only."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from campaign import STUDY,digest,now,save


def run(paths,prefix):
    prefix=Path(prefix).resolve()
    outputs=[prefix.with_suffix(suffix) for suffix in ('.pdf','.png','.json')]
    if any(path.exists() for path in outputs):
        raise ValueError('Refuse to overwrite W-treatment pilot figure artifacts')
    reports={}
    for path in map(Path,paths):
        report=json.loads(path.read_text())
        if (report['status']!='conditional three-W-treatment pilot comparison; full-flavour retraining required'
                or report['flavours']!=['e','e','mu'] or report['charge'] in reports
                or digest(path.with_suffix('.npz'))!=report['arrays_sha256']
                or digest(STUDY/'inputs/benchmark.json')!=report['benchmark_sha256']):
            raise ValueError('Require both unchanged matching single-assignment W-treatment pilots')
        reports[report['charge']]=report
    if set(reports)!={'plus','minus'}:
        raise ValueError('Require both charges')
    comparisons=('top_W','associated_W','all_W')
    ticklabels=('top BW /\non-shell','all BW /\ntop BW','all BW /\non-shell')
    points={}
    fig,axes=plt.subplots(1,2,figsize=(8.2,3.8),sharey=True)
    for ax,rate,title in zip(axes,('fiducial_1b','fiducial_2b'),('At least one b jet','At least two b jets')):
        points[rate]={}
        for charge,label,color,marker,offset in (('plus',r'$W^+$','#2865ad','o',-.07),
                                                 ('minus',r'$W^-$','#c94b32','s',.07)):
            rows=[reports[charge]['rates']['R04_b25'][rate][name]['product_ratio_double_ratio'] for name in comparisons]
            value=np.array([row['scale']['central'] for row in rows])
            error=np.array([row['nominal_conditional_mc_error'] for row in rows])
            if not np.isfinite(value).all() or not np.isfinite(error).all() or (error<=0.).any():
                raise ValueError('Invalid nominal double-ratio values or conditional errors')
            points[rate][charge]=dict(value=value.tolist(),conditional_mc_error=error.tolist())
            ax.errorbar(np.arange(3)+offset,value,yerr=error,fmt=marker,color=color,
                        label=label,capsize=3,markersize=4.5,linewidth=1.1)
        ax.axhline(1.,color='0.5',linestyle='--',linewidth=.8,zorder=0)
        ax.set_xticks(range(3),ticklabels)
        ax.set_title(title,fontsize=11)
        ax.grid(axis='y',alpha=.2)
        ax.set_xlim(-.35,2.35)
    axes[0].set_ylabel(r'$(\Pi/S)_{\rm target}/(\Pi/S)_{\rm reference}$')
    fig.suptitle(r'W-treatment sensitivity; ordered $(e,e,\mu)$; common W-system $H_T/2$',fontsize=11,y=.98)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(.5,.91))
    fig.text(.5,.02,'Pilot diagnostics: conditional MC errors. Treatment ratios share samples.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.065,1,.90))
    prefix.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(outputs[0]);fig.savefig(outputs[1],dpi=180);plt.close(fig)
    save(outputs[2],dict(created_utc=now(),status='single-assignment conditional pilot figure; not a full-flavour prediction',
        inputs={str(Path(path).resolve()):digest(path) for path in paths},comparison_order=comparisons,
        points=points,script_sha256=digest(Path(__file__)),outputs={str(path):digest(path) for path in outputs[:2]}))
    print(outputs[0])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reports',nargs=2,type=Path,required=True)
    parser.add_argument('--output-prefix',type=Path,required=True)
    args=parser.parse_args()
    run(args.reports,args.output_prefix)
