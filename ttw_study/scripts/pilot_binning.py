#!/usr/bin/env python3
"""Assess candidate bins for four priority shapes from audited joint batches.

Candidates are NOT frozen publication bins. Their estimates come from only
one ordered flavour/prescription. The required other flavours and strict
prediction must be checked before any binning/convergence certification.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from campaign import digest, now, save
from joint_report import estimate, normalization_bins
from pilot_report import clean
from rebinning import rebin
from replica_statistics import ratio

CANDIDATES = {
    4: [0,60,100,140,180,240,320,440,600],
    11:[0,20,30,40,50,60,80,100,140,200,400],
    12:[0,30,50,70,90,110,130,150,160,180,220,280,400],
    20:[0,20,40,60,80,100,150,200,300,400],
}
LABELS = {4:r'$H_T^{\ell}$ [GeV]',11:r'$p_T^{b_2}$ [GeV]',
          12:r'$m_{b\ell}^{\rm minimax}$ [GeV]',20:r'$p_T^{\ell_1}$ [GeV]'}


def run(path,output):
    if output.exists() or output.with_suffix('.pdf').exists():
        raise ValueError('Candidate study already exists')
    metadata=json.loads(path.with_suffix('.json').read_text())
    if digest(path) != metadata['arrays_sha256']:
        raise ValueError('Joint batch archive checksum mismatch')
    with np.load(path) as archive:
        values,strata,offsets,edges = [archive[key] for key in ('contributions','strata','offsets','edges')]
        titles=list(archive['titles'])
    start_h=next(h for h,title in enumerate(titles) if title.startswith('R04_b25')
                 and ' rates:' in title and np.any(values[:,offsets[h],0]))
    rate_start=int(offsets[start_h])
    report=dict(created_utc=now(),input=str(path.resolve()),input_sha256=digest(path),
                status='candidate pilot bins, not frozen; incomplete flavours and prescriptions',
                notes='All bins are retained; raw batch weights are rebinned before jackknife. '
                      'MC cost multipliers assume inverse-square-root scaling and must be remeasured. '
                      'No total-rate error substitutes for a shape-bin error. Overflows remain outside '
                      'the original finite histogram range and the rates supply normalization.',
                spectra={})
    fig,axes=plt.subplots(2,2,figsize=(9.,6.6),constrained_layout=True)
    for (local,requested),ax in zip(CANDIDATES.items(),axes.flat):
        h=start_h+local-1
        first,last=map(int,offsets[h:h+2])
        bins=rebin(values[:,first:last,:],edges[first:last],requested)
        denominator=values[:,[rate_start+i for i in normalization_bins(local)],:].sum(axis=1)
        joint=np.concatenate([bins,denominator[:,None,:]],axis=1)
        result=estimate(joint,strata,lambda total: ratio(total[:-1],total[-1]))
        coverage=estimate(joint,strata,lambda total: ratio(total[:-1].sum(axis=0),total[-1]))
        central,error=result['value'][:,0],result['mc_error'][:,0]
        relative=ratio(error,central)
        cost=(relative/.02)**2
        report['spectra'][titles[h]]=dict(edges=requested,normalized_bin_integrals=central,
            joint_mc_errors=error,relative_errors=relative,cost_multiplier_to_2percent=cost,
            central_coverage=coverage['value'][0],central_coverage_error=coverage['mc_error'][0],
            all_weights_normalized=result['value'],all_weights_joint_errors=result['mc_error'])
        widths=np.diff(requested)
        centers=np.asarray(requested[:-1])+widths/2
        ax.errorbar(centers,central/widths,yerr=error/widths,xerr=widths/2,
                    fmt='o',markersize=3.5,linewidth=1.,color='#286b8e',capsize=2)
        ax.set_xlabel(LABELS[local])
        ax.set_ylabel(r'$(1/\sigma_{\rm region})\,d\sigma/dx$ [GeV$^{-1}$]')
        ax.set_xlim(requested[0],requested[-1])
        ax.set_ylim(bottom=min(0.,1.1*float(np.nanmin((central-error)/widths))))
        ax.grid(axis='y',alpha=.2)
        ax.set_title(r'In-range fraction %.3f $\pm$ %.3f' % (coverage['value'][0],coverage['mc_error'][0]),fontsize=9)
    fig.suptitle('Candidate binning: single-flavour pilot; conditional MC errors only',fontsize=12)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output.with_suffix('.pdf'))
    fig.savefig(output.with_suffix('.png'),dpi=150)
    plt.close(fig)
    save(output,clean(report))
    print(output,flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('batches',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.batches,args.output)
