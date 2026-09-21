#!/usr/bin/env python3
"""Rate/acceptance and scale-response tables from complete main reductions."""
import argparse
import math
from pathlib import Path

import numpy as np

from campaign import digest, now, save
from main_figures import MainData, numeric, row_values

VARIANTS=('LO','P','D','S','PiD','Pi')
LABELS=dict(LO=r'$\mathrm{LO}$',P='$P$',D='$D$',S='$S$',PiD=r'$\Pi_D$',Pi=r'$\Pi$')


def estimate_tex(row,factor=1.):
    value,error=row_values(row)[:2]*factor
    if not np.isfinite(value) or not np.isfinite(error):return r'$\mathrm{undefined}$'
    if error<0.:raise ValueError('Negative MC error in a main table')
    if error==0.:return '$%s\\pm0$'%format(value,'.6g')
    decimals=1-int(math.floor(math.log10(error)))
    value,error=round(value,decimals),round(error,decimals)
    digits=max(0,decimals)
    return '$%s\\pm%s$'%(format(value,'.%df'%digits),format(error,'.%df'%digits))


def band_tex(row,band,relative=False):
    central,_,low,high=row_values(row,band)
    if not np.isfinite([low,high]).all() or (relative and (not np.isfinite(central) or central<=0.)):
        return r'$\mathrm{undefined}$'
    if relative:
        low,high=100*(np.array([low,high])/central-1.)
        return '$[%+.1f,%+.1f]$'%(low,high)
    return '$[%.3f,%.3f]$'%(low,high)


def scalar_tex(value,percent=False,signed=False):
    value=numeric(value)
    if not np.isfinite(value):return r'$\mathrm{undefined}$'
    return '$%s$'%format(value*(100 if percent else 1),('+.2f' if signed else '.2f') if percent else '.3g')


def table(rows,columns):
    return '\\begin{tabular}{'+columns+'}\n\\toprule\n'+'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}\n'


def line(values):
    return ' & '.join(values)+r' \\'


def captions(data,text):
    warning=(r'\textbf{SYNTHETIC LAYOUT TEST; no physics predictions.} '
             if data.scope=='synthetic_plot_layout_only' else '')
    return '\\caption{'+warning+text+'}\n'


def texts(data):
    rates=[line(['Prescription',r'$\sigma^+_{\geq1b}$',r'$\sigma^-_{\geq1b}$',
                 r'$\sigma^+_{\geq2b}$',r'$\sigma^-_{\geq2b}$']),r'\midrule']
    efficiencies=[line(['Prescription',r'$A^+_{1b}$',r'$A^-_{1b}$',r'$A^+_{2b}$',r'$A^-_{2b}$',
                        r'$R_{\rm charge}^{2b}$',r'$A_{\rm charge}^{2b}$']),r'\midrule']
    charge_title=data.rate_title().replace('W+','W+/W-')
    charge_rows=data.report['charge_rates'][charge_title]['fiducial_2b']
    for variant in VARIANTS:
        rates.append(line([LABELS[variant]]+[estimate_tex(data.rates(selection,charge)[variant],1000.)
            for selection in ('fiducial_1b','fiducial_2b') for charge in ('plus','minus')]))
        efficiencies.append(line([LABELS[variant]]+[estimate_tex(data.rates(selection,charge,derived=True)[variant])
            for selection in ('acceptance_1b','acceptance_2b') for charge in ('plus','minus')]+
            [estimate_tex(charge_rows[variant][quantity]) for quantity in ('plus_over_minus','asymmetry')]))
    prefix='\\begin{table}[htbp]\n\\centering\n\\small\n\\setlength{\\tabcolsep}{3pt}\n'
    rate_text=prefix+table(rates,'lrrrr')+'\\par\\medskip\n'+table(efficiencies,'lrrrrrr')
    rate_text+=captions(data,r'Full direct $e/\mu$ flavour sum at 13 TeV with on-shell Ws and massless decay bottoms. '
        r'Cross sections are in fb; displayed errors are independent-retraining MC errors only. '
        r'$A_{1b,2b}=\sigma_{1b,2b}/\sigma_{\rm input}$, '
        r'$R_{\rm charge}=\sigma^+/\sigma^-$ and '
        r'$A_{\rm charge}=(\sigma^+-\sigma^-)/(\sigma^++\sigma^-)$. '
        r'All prescriptions use the same NLO PDF; P retains its LO top-width convention.')
    rate_text+='\\label{tab:mainrates}\n\\end{table}\n'
    shifts=[line(['Region',r'$\Delta_\Pi$ [fb]',r'$R_\Pi$',r'Prod. 7',r'Decay 9',r'Combined 21',r'Combined 63',
                  r'$\delta_{\rm PDF}R_\Pi$']),r'\midrule']
    for charge,sign in (('plus','+'),('minus','-')):
        for selection,region in (('fiducial_1b','1b'),('fiducial_2b','2b')):
            rows=data.rates(selection,charge)
            ratio=rows['Pi_over_S']
            shifts.append(line([r'$W^{%s},\,\geq%s$'%(sign,region),estimate_tex(rows['Pi_minus_S'],1000.),
                estimate_tex(ratio)]+[band_tex(ratio,band) for band in ('production7','decay9','combined21','combined63')]+
                [scalar_tex(ratio['pdf']['error_symmetric'])]))
    responses=[line(['Charge','Prescription',r'MC [\%]',r'Prod. 7 [\%]',r'Decay 9 [\%]',r'Combined 63 [\%]',
                     r'$\rho_P$ [\%]',r'$\rho_D$ [\%]']),r'\midrule']
    for charge,sign in (('plus','+'),('minus','-')):
        for variant in ('P','S','Pi'):
            row=data.rates('fiducial_2b',charge)[variant]
            value,error=row_values(row)[:2]
            mc=error/value if value>0. else np.nan
            slopes=row['scale']['log_responses']
            responses.append(line(['$W^{%s}$'%sign,LABELS[variant],scalar_tex(mc,percent=True)]+
                [band_tex(row,band,relative=True) for band in ('production7','decay9','combined63')]+
                [scalar_tex(slopes[name]['relative'],percent=True,signed=True) for name in ('production','decay')]))
    response_text=prefix+table(shifts,'lrrrrrrr')+'\\par\\medskip\n'+table(responses,'llrrrrrr')
    response_text+=captions(data,r'Point-matched $\Delta_\Pi=\Pi-S$ and $R_\Pi=\Pi/S$ before envelope/PDF reduction. '
        r'Top: nominal values with MC errors, scale-only intervals for $R_\Pi$, and its separate one-standard-deviation '
        r'NNPDF replica error. Bottom: two-b P/S/$\Pi$ MC precision and scale-envelope deviations from each nominal '
        r'cross section. The logarithmic responses are '
        r'$\rho_X=[\sigma(\xi_X=2)-\sigma(\xi_X=1/2)]/[\sigma_0\ln4]$; '
        r'production varies $\mu_R=\mu_F$ together and decay varies top/antitop together. '
        r'These scale ranges and responses are not Gaussian errors or evidence of NNLO accuracy.')
    response_text+='\\label{tab:mainresponses}\n\\end{table}\n'
    return {'main_rates.tex':rate_text,'main_responses.tex':response_text}


def render(data,output):
    if output.exists():raise ValueError('Refuse to overwrite main table output')
    content=texts(data)
    output.mkdir(parents=True)
    outputs={}
    for name,value in content.items():
        path=output/name
        path.write_text(value)
        outputs[str(path.resolve())]=digest(path)
    save(output/'manifest.json',dict(created_utc=now(),data_scope=data.scope,
        status=('synthetic table-layout test; no physics predictions' if data.scope=='synthetic_plot_layout_only' else
                'main numerical tables rendered; precision and paper interpretation require inspection'),
        inputs=data.inputs,outputs=outputs,source_hashes={str(Path(__file__).resolve()):digest(Path(__file__)),
            str(Path(__file__).with_name('main_figures.py').resolve()):digest(Path(__file__).with_name('main_figures.py'))},
        convention='Rates converted from pb to fb. Ratios/acceptances and all band coordinates use the '
            'complete full-run reduction with retained covariance. MC, PDF and scale quantities remain separate.'))
    print(output,flush=True)


def run(path,output):
    if output.exists():raise ValueError('Refuse to overwrite main table output')
    render(MainData(path),output)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args()
    run(args.report,args.output_dir)
