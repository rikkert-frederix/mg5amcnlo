#!/usr/bin/env python3
"""Standalone generated narrow-W diagnostic figure and reproducibility record."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from campaign import STUDY, digest, now, save


def main():
    paths = [STUDY/'results'/('narrow_v2_'+v+'_generated_limit.json') for v in ('LO','S','Pi')]
    reports = [json.loads(p.read_text()) for p in paths]
    for path, report in zip(paths, reports):
        if digest(path.with_suffix('.npz')) != report['arrays_sha256']:
            raise ValueError('Changed narrow-W arrays')
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.6), sharex=True)
    for col, report in enumerate(reports):
        rates = next(v for k,v in report['rates'].items() if k.startswith('R04_b25'))
        for row, key in enumerate(('fiducial_1b','fiducial_2b')):
            ax = axes[row,col]
            for mode, color, marker, offset in [('top-bw','#2266aa','o',.96),('all-bw','#bc5035','s',1.04)]:
                selected = [i for i,s in enumerate(report['samples']) if s['case']['w_treatment']==mode]
                x = np.asarray([report['samples'][i]['case']['width_factor'] for i in selected])
                values = np.asarray(rates[key]['ratio'])[selected]
                errors = np.asarray(rates[key]['ratio_mc_error'])[selected]
                ax.errorbar(x*offset, values, yerr=errors, color=color, marker=marker,
                            linestyle='none', capsize=3, markersize=4, label=mode)
            ax.axhline(1., color='0.4', linewidth=.8)
            ax.set_xscale('log')
            ax.set_xlim(.036,1.4)
            ax.set_xticks([.05,.2,1.], ['0.05','0.2','1'])
            ax.grid(axis='y',alpha=.2)
            if row==0:
                ax.set_title({'LO':'LO','S':'Strict NLO (S)','Pi':r'Product ($\Pi$)'}[report['variant']])
            if col==0:
                ax.set_ylabel(('One-b' if row==0 else 'Two-b')+' rate ratio')
            if row==1:
                ax.set_xlabel(r'$\epsilon=\Gamma_W/\Gamma_W^{\rm ref}$')
    axes[0,0].legend(frameon=False,fontsize=9)
    fig.suptitle(r'Generated $\epsilon^3\sigma_{\rm BW}/\sigma_{\rm on\ shell}$; '
                 r'$W^+$, ordered $e,e,\mu$; common fixed scales', fontsize=12)
    fig.text(.5,.015,'Pilot diagnostics; conditional Monte Carlo errors. Width points share the on-shell reference.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.045,1,.94])
    destination=STUDY/'figures/generated_narrow_w'
    destination.parent.mkdir(exist_ok=True)
    for extension in ('.pdf','.png'):
        fig.savefig(destination.with_suffix(extension),dpi=170)
    plt.close(fig)
    save(destination.with_suffix('.json'),dict(created_utc=now(),
        inputs={str(p):digest(p) for p in paths},script_sha256=digest(Path(__file__)),
        outputs={ext:digest(destination.with_suffix(ext)) for ext in ('.pdf','.png')},
        scope='One-assignment technical limit diagnostic, not main phenomenology. '
        'Horizontal offsets distinguish modes; all use exactly epsilon=1,0.2,0.05. '
        'Each panel has its own vertical scale; errors are conditional joint MC, not scale/PDF bands.'))


if __name__=='__main__':
    main()
