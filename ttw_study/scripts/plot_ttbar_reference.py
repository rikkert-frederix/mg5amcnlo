#!/usr/bin/env python3
"""Standalone matched ttbar angular calibration with explicit MC residuals."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from campaign import STUDY, digest, now, save


def run(report_path,output):
    if any(output.with_suffix(ext).exists() for ext in ('.pdf','.png','.json')):
        raise ValueError('Refuse to overwrite a ttbar calibration figure')
    report=json.loads(report_path.read_text())
    preflight_path=STUDY/'inputs/ttbar_reference_preflight.json'
    if digest(preflight_path)!=report['reference_preflight_sha256']:
        raise ValueError('ttbar reference inputs changed')
    preflight=json.loads(preflight_path.read_text())
    inputs={str(report_path.resolve()):digest(report_path),str(preflight_path):digest(preflight_path)}
    for path in map(Path,report['batches'].values()):
        metadata=json.loads(path.with_suffix('.json').read_text())
        if digest(path)!=metadata['arrays_sha256']:
            raise ValueError('ttbar reference batch evidence changed')
        inputs[str(path)]=digest(path)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'ps.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(9.2,6.),sharex='col',
        gridspec_kw={'height_ratios':[3,1.35]},layout='constrained')
    for col,(label,expanded) in enumerate((('NLO_normalized',False),('NLO_expanded_normalized',True))):
        name='norm-incl-'+('exp-' if expanded else '')+'dPhill-NLO-NNPDF-mt.dat'
        path=STUDY/'references/data-expanded-ratios'/name
        if digest(path)!=preflight['reference']['curves'][name]['sha256']:
            raise ValueError('Published curve changed')
        inputs[str(path)]=digest(path)
        curve=np.loadtxt(path)
        result=report['comparisons'][label]
        central=np.asarray(result['value'])[:,0]
        error=np.asarray(result['mc_error'])[:,0]
        target=np.asarray(result['reference'])[:,0]
        target_error=np.asarray(result['reference_central_mc_error'])
        if not np.array_equal(curve[:,1],target) or not np.array_equal(curve[:,8],target_error):
            raise ValueError('Comparison does not match the supplied curve')
        x=curve[:,0]
        axes[0,col].errorbar(x,target,yerr=target_error,color='0.2',marker='s',markersize=3,
            linewidth=1.,capsize=2,label='Supplied NLO reference')
        axes[0,col].errorbar(x,central,yerr=error,color='#1565a3',marker='o',markersize=4,
            linestyle='none',capsize=3,label='Generated strict NLO (S)')
        axes[0,col].set_title('Expanded normalized distribution' if expanded else 'NLO differential / NLO total')
        axes[0,col].set_ylim(.68,1.42)
        combined=np.hypot(error,target_error)
        pull=(central-target)/combined
        axes[1,col].axhspan(-1,1,color='0.94')
        axes[1,col].axhline(0,color='0.45',linewidth=.8)
        axes[1,col].plot(x,pull,'o',color='#1565a3',markersize=4)
        axes[1,col].set_ylim(-3,3)
        axes[1,col].set_yticks([-2,0,2])
        axes[1,col].set_xlabel(r'$\Delta\phi(\ell^+,\ell^-)/\pi$')
        axes[1,col].set_xlim(0,1)
        for row in (0,1):axes[row,col].grid(axis='y',alpha=.2)
    axes[0,0].set_ylabel(r'Normalized density in $\Delta\phi/\pi$')
    axes[1,0].set_ylabel('Residual /\ncombined MC error')
    axes[0,0].legend(frameon=False,loc='upper left',fontsize=9)
    fig.suptitle(r'Inclusive dilepton $t\bar t$ calibration; $\mu_R=\mu_F=m_t$, NNPDF3.1 NLO',fontsize=12)
    fig.supxlabel('Pilot MC errors conditional on trained grids; correlated bin residuals are not a global fit test.',fontsize=9)
    output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ('.pdf','.png'):fig.savefig(output.with_suffix(ext),dpi=180)
    plt.close(fig)
    save(output.with_suffix('.json'),dict(created_utc=now(),inputs=inputs,script_sha256=digest(Path(__file__)),
        outputs={ext:digest(output.with_suffix(ext)) for ext in ('.pdf','.png')},
        scope='Matched central-scale NLO pilot calibration only. Error bars show MC uncertainty, not scale bands. '
            'Each residual uses the generated conditional error and supplied reference MC error in quadrature; '
            'the reference inter-bin covariance is unavailable. The expanded label refers to the normalized '
            'distribution, not to an unexpanded total top-width denominator.'))
    print(output.with_suffix('.pdf'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.report,args.output)
