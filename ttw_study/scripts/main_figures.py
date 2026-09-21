#!/usr/bin/env python3
"""Draft paper figures 2--7 from complete full-flavour main reductions.

Bands are scale envelopes of the plotted observable; bars are independent-
retraining MC errors. PDF uncertainty remains separate in the source report.
Figure creation is not a numerical-convergence certification.
"""
import argparse
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

from campaign import ROOT, STUDY, digest, now, save
from load_results import WEIGHTS
from main_replicas import inventory
from pilot_report import clean, scales

COLORS=dict(S='#0072B2',Pi='#D55E00',P='#009E73',PiD='#CC79A7')
NAMES=dict(S='S',Pi=r'$\Pi$',P='P',PiD=r'$\Pi_D$')
BANDS=('production7','top3','antitop3','decay3','decay9','combined21','combined63','common3','all81_diagnostic')
BAND_LABELS=('Prod.\n7','Top\n3','Antitop\n3','Shared decay\n3','Independent decay\n9',
             'Combined\n21','Combined\n63','Common\n3','Full grid\n81')
XLABELS={3:r'$H_T^j$ [GeV]',4:r'$H_T^\ell$ [GeV]',5:r'$\Delta R(\ell_1,b_1)$',
         6:r'$|\Delta\phi_{\rm SS}|$ [rad]',7:r'$|\Delta\eta_{\rm SS}|$',8:r'$S_T$ [GeV]',
         11:r'$p_{T,b_2}$ [GeV]',12:r'$m_{\ell b}^{\rm minimax}$ [GeV]',
         16:r'$p_{T,j_2}^{\rm extra}$ [GeV]'}


def numeric(value):
    return float(value) if value is not None else np.nan


def row_values(row,band='combined63'):
    scale=row['scale']
    points=scale.get('points',{})
    if set(points)!=set(WEIGHTS[1:82]):
        raise ValueError('A plotted observable lacks the complete signed scale grid')
    selected=[points[','.join('%g'%v for v in point)] for point in sorted(scales.GROUPS[band])]
    finite=[numeric(value) for value in selected if value is not None and np.isfinite(value)]
    expected=[min(finite),max(finite)] if len(finite)==len(selected) else None
    reported=scale[band]
    if (reported['npoints']!=len(selected) or reported['nvalid']!=len(finite) or
            reported['envelope']!=expected):
        raise ValueError('Plotted band differs from its declared matching scale points')
    low,high=expected if expected is not None else (np.nan,np.nan)
    return np.asarray([numeric(scale['central']),numeric(row['nominal_retraining_mc_error']),low,high])


def series(detail,label,kind='normalized',band='combined63',density=True):
    edges=np.asarray(detail['edges'],dtype=float)
    if (edges.ndim!=2 or edges.shape[1]!=2 or len(edges)!=len(detail[kind]) or
            np.any(edges[:,1]<=edges[:,0]) or not np.array_equal(edges[1:,0],edges[:-1,1])):
        raise ValueError('Figure bins must be complete, ordered and contiguous')
    values=np.asarray([row_values(row[label],band) for row in detail[kind]])
    if density:
        values=values/np.diff(edges,axis=1)
    return edges,values


class MainData:
    def __init__(self,path):
        self.path=path.resolve()
        self.report=json.loads(self.path.read_text())
        report=self.report
        if report['status']!='full-flavour main observables reduced; inspect precision and convergence':
            raise ValueError('Require a complete full-flavour main reduction')
        source_path=Path(report['input'])
        source=json.loads(source_path.read_text())
        if source['status']!='complete independent main vectors collected; flavour reduction and convergence required':
            raise ValueError('Main figures require independently retrained complete source vectors')
        if (digest(source_path)!=report['input_sha256'] or
                digest(source_path.with_suffix('.npz'))!=report['input_arrays_sha256'] or
                report['input_arrays_sha256']!=source['arrays_sha256'] or
                digest(source['input'])!=source['input_sha256'] or
                digest(self.path.with_suffix('.npz'))!=report['arrays_sha256'] or
                digest(STUDY/'inputs/benchmark.json')!=report['benchmark_sha256']):
            raise ValueError('Main figure input, replica, campaign or benchmark evidence changed')
        count=inventory(json.loads(Path(source['input']).read_text()))
        if (len(source['groups'])!=96 or len(source['samples'])!=96*count or
                any(group['replicas']!=list(range(count)) for group in source['groups'].values())):
            raise ValueError('Main figures require every ordered flavour and independent retraining')
        self.inputs={str(p.resolve()):digest(p) for p in
            (self.path,self.path.with_suffix('.npz'),source_path,source_path.with_suffix('.npz'),
             Path(source['input']),STUDY/'inputs/benchmark.json')}
        with np.load(self.path.with_suffix('.npz')) as arrays:
            self.titles=list(arrays['titles'])
            self.offsets=arrays['offsets'].copy()
            if list(arrays['weights'])!=WEIGHTS or len(self.titles)!=210:
                raise ValueError('Main figure inventory is not the complete canonical histogram bank')
        self.banner=r'13 TeV; on-shell W; $m_b^{\rm decay}=0$; full direct $e/\mu$ flavour sum'
        self.scope='complete_full_flavour_main_reduction'

    def rate_title(self,charge='plus',config='R04_b25'):
        prefix=config+(' W+' if charge=='plus' else ' W-')+' rates:'
        found=[title for title in self.report['rates'] if title.startswith(prefix)]
        if len(found)!=1:raise ValueError('Missing or ambiguous plotted rate bank')
        return found[0]

    def rates(self,label,charge='plus',config='R04_b25',derived=False):
        return self.report['derived_rates' if derived else 'rates'][self.rate_title(charge,config)][label]

    def histogram(self,local,charge='plus',charge_comparison=False):
        h=self.titles.index(self.rate_title(charge))+local-1
        title=self.titles[h]
        if charge_comparison:title=title.replace('W+','W+/W-')
        entry=self.report['charge_spectra' if charge_comparison else 'spectra'][title]
        path=Path(entry['path'])
        if digest(path)!=entry['sha256']:
            raise ValueError('Plotted spectrum changed after the main reduction')
        self.inputs[str(path.resolve())]=entry['sha256']
        detail=json.loads(gzip.decompress(path.read_bytes()))
        if detail['title']!=title:
            raise ValueError('Plotted spectrum title differs from its inventory')
        return detail


def band_curve(ax,edges,values,color,label,linestyle='-',fill=True):
    x=edges.mean(axis=1)
    boundaries=np.r_[edges[:,0],edges[-1,1]]
    central,error,low,high=values.T
    # Matplotlib 3.6 still inserts zero into autoscaling for stairs even
    # with baseline=None. An explicit post-step line has no invented base.
    ax.step(boundaries,np.r_[central,central[-1]],where='post',color=color,
            linewidth=1.35,linestyle=linestyle,label=label)
    if fill:
        ax.fill_between(boundaries,np.r_[low,low[-1]],np.r_[high,high[-1]],step='post',color=color,alpha=.16)
    valid=np.isfinite(central)&np.isfinite(error)&(error>=0.)
    ax.errorbar(x[valid],central[valid],yerr=error[valid],fmt='.',color=color,markersize=2.8,
                linewidth=.65,capsize=1.2)
    ax.set_xlim(boundaries[0],boundaries[-1])
    ax.grid(axis='y',alpha=.18)


def shape_figure(data,locals_,charge_panel=False):
    fig,axes=plt.subplots(2,len(locals_),figsize=(7.8,5.7),sharex='col',
        gridspec_kw={'height_ratios':[3,1.5]},layout='constrained',squeeze=False)
    for col,local in enumerate(locals_):
        detail=data.histogram(local)
        for variant in ('P','S','Pi'):
            edges,values=series(detail,variant)
            band_curve(axes[0,col],edges,values,COLORS[variant],NAMES[variant],
                       linestyle='--' if variant=='P' else '-',fill=variant!='P')
        parent='/'.join(detail['normalization_rate_bins'])
        axes[0,col].set_title('W+; '+parent.replace('fiducial_','')+'; '+XLABELS[local])
        axes[0,col].set_ylabel(r'$1/\sigma_{\rm fid}\;d\sigma/dX$')
        axes[0,col].ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useMathText=True)
        axes[0,col].yaxis.set_major_locator(MaxNLocator(nbins=4))
        if charge_panel:
            charge=data.histogram(local,charge_comparison=True)
            selected=dict(charge,normalized=[{variant:row[variant]['plus_shape_over_minus_shape']
                for variant in ('S','Pi')} for row in charge['normalized']])
            for variant in ('S','Pi'):
                edges,values=series(selected,variant,density=False)
                band_curve(axes[1,col],edges,values,COLORS[variant],NAMES[variant])
            axes[1,col].set_ylabel(r'$f_+(X)/f_-(X)$')
        else:
            edges,values=series(detail,'Pi_over_S',density=False)
            band_curve(axes[1,col],edges,values,COLORS['Pi'],r'$\Pi/S$: 63 points')
            _,shared=series(detail,'Pi_over_S',band='combined21',density=False)
            x=np.r_[edges[:,0],edges[-1,1]]
            for boundary in (2,3):
                axes[1,col].step(x,np.r_[shared[:,boundary],shared[-1,boundary]],where='post',
                    color=COLORS['Pi'],linestyle=':',linewidth=.8,
                    label='21-point edges' if boundary==2 else None)
            axes[1,col].set_ylabel(r'Normalized $\Pi/S$')
        axes[1,col].axhline(1.,color='0.5',linewidth=.65)
        axes[1,col].set_xlabel(XLABELS[local])
        axes[1,col].xaxis.set_major_locator(MaxNLocator(nbins=4))
    axes[0,0].legend(frameon=False,ncol=3,fontsize=9)
    axes[1,0].legend(frameon=False,fontsize=8)
    fig.suptitle(data.banner,fontsize=11)
    fig.supxlabel('S/Pi: scale-only 63-point bands; P: nominal control. Bars: retraining MC errors.\nMatched parent-rate normalization.',fontsize=8)
    return fig


def scale_figure(data):
    fig,axes=plt.subplots(2,1,figsize=(7.8,6.4),sharex=True,layout='constrained')
    x=np.arange(len(BANDS))
    for ax,charge in zip(axes,('plus','minus')):
        rows=data.rates('fiducial_2b',charge)
        for variant,offset in (('S',-.12),('Pi',.12)):
            values=np.asarray([row_values(rows[variant],band) for band in BANDS])*1000.
            central,error,low,high=values.T
            ax.vlines(x+offset,low,high,color=COLORS[variant],linewidth=5.,alpha=.4)
            ax.errorbar(x+offset,central,yerr=error,fmt='o',color=COLORS[variant],
                        markersize=3,capsize=2,linewidth=.9,label=NAMES[variant])
        ax.set_xticks(x,BAND_LABELS,fontsize=8,rotation=30,ha='right')
        ax.set_title('W'+('+' if charge=='plus' else '-')+'; two-b fiducial rate')
        ax.set_ylabel(r'$\sigma_{\geq2b}$ [fb]')
        ax.grid(axis='y',alpha=.2)
    axes[0].legend(frameon=False,ncol=2)
    fig.suptitle(data.banner,fontsize=11)
    fig.supxlabel('Thick ranges: the stated scale-point envelopes.\nMarkers and thin bars: nominal values and retraining MC errors.',fontsize=8)
    return fig


def extra_jet_figure(data):
    fig,axes=plt.subplots(1,2,figsize=(7.8,4.3),layout='constrained')
    x=np.arange(4)
    lower=0.
    for variant,offset in (('S',-.2),('PiD',0.),('Pi',.2)):
        values=np.asarray([row_values(data.rates(label)[variant]) for label in
            ('2b_0extra','2b_1extra','2b_2extra','2b_3plus_extra')])*1000.
        central,error,low,high=values.T
        lower=min(lower,float(np.nanmin(np.r_[low,central-error])))
        axes[0].vlines(x+offset,low,high,color=COLORS[variant],linewidth=5.,alpha=.4)
        axes[0].errorbar(x+offset,central,yerr=error,fmt='o',markersize=4,color=COLORS[variant],
                        capsize=2,linewidth=.8,label=NAMES[variant])
    axes[0].set_xticks(x,['0','1','2',r'$\geq3$'])
    axes[0].set_xlabel('Extra jets in the two-b selection')
    axes[0].set_ylabel('Cross section [fb]')
    axes[0].set_yscale('symlog',linthresh=.01)
    axes[0].set_ylim(bottom=1.2*lower if lower<0. else -.002)
    axes[0].grid(axis='y',alpha=.2)
    axes[0].legend(frameon=False,ncol=3,fontsize=9)
    detail=data.histogram(16)
    strict=series(detail,'S',kind='absolute')[1]
    if np.any(strict[:,[0,2,3]]!=0.):
        raise ValueError('Strict NLO must vanish in the second-extra-jet spectrum')
    for variant in ('PiD','Pi'):
        edges,values=series(detail,variant,kind='absolute')
        band_curve(axes[1],edges,values*1000.,COLORS[variant],NAMES[variant])
    axes[1].axhline(0.,color=COLORS['S'],linewidth=1.,label='S = 0')
    axes[1].set_xlabel(XLABELS[16])
    axes[1].set_ylabel(r'$d\sigma/dp_T$ [fb/GeV]')
    axes[1].legend(frameon=False,fontsize=9)
    fig.suptitle(data.banner+'; W+',fontsize=11)
    fig.supxlabel('63-point scale bands and retraining MC bars. Product-only sectors contain partial higher-order terms.\nNo ratio to zero S is formed.',fontsize=8)
    return fig


def migration_figure(data):
    fig,axes=plt.subplots(2,2,figsize=(7.8,6.),layout='constrained')
    choices=((np.array([25.,30.,40.]),('R04_b25','R04_b30','R04_b40'),r'$b$-jet threshold [GeV]; $R=0.4$'),
             (np.array([.3,.4,.5]),('R03_b25','R04_b25','R05_b25'),r'Jet radius $R$; $p_{T,b}\geq25$ GeV'))
    for col,(x,configs,xlabel) in enumerate(choices):
        for row,(quantity,ylabel) in enumerate((('two_b_given_one_b',r'$\sigma_{\geq2b}/\sigma_{\geq1b}$'),
                                               ('veto_0extra',r'$\sigma_{2b,0j}/\sigma_{\geq2b}$'))):
            ax=axes[row,col]
            for variant in ('S','Pi'):
                values=np.asarray([row_values(data.rates(quantity,config=config,derived=True)[variant]) for config in configs])
                central,error,low,high=values.T
                ax.fill_between(x,low,high,color=COLORS[variant],alpha=.16)
                ax.errorbar(x,central,yerr=error,fmt='o-',markersize=4,linewidth=1.,capsize=2,
                            color=COLORS[variant],label=NAMES[variant])
            ax.set_xlabel(xlabel)
            ax.set_xticks(x)
            ax.set_ylabel(ylabel)
            ax.grid(axis='y',alpha=.2)
    axes[0,0].legend(frameon=False,ncol=2)
    fig.suptitle(data.banner+'; W+',fontsize=11)
    fig.supxlabel('63-point scale bands; bars: retraining MC errors.\nCut selections and the two panels share events; their points are correlated.',fontsize=8)
    return fig


def render(data,output):
    if output.exists():raise ValueError('Refuse to overwrite a main figure directory')
    output.mkdir(parents=True)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'ps.fonttype':42})
    specifications=(('figure_2_scales',lambda:scale_figure(data)),
        ('figure_3_decay',lambda:shape_figure(data,(11,12))),
        ('figure_4_activity',lambda:shape_figure(data,(3,8,4))),
        ('figure_5_angles_charge',lambda:shape_figure(data,(6,7,5),charge_panel=True)),
        ('figure_6_extra_jets',lambda:extra_jet_figure(data)),
        ('figure_7_migrations',lambda:migration_figure(data)))
    assets={}
    for name,draw in specifications:
        fig=draw()
        files={}
        for suffix in ('.pdf','.png'):
            path=output/(name+suffix)
            fig.savefig(path,dpi=170,bbox_inches='tight',pad_inches=.06)
            files[str(path.resolve())]=digest(path)
        plt.close(fig)
        assets[name]=files
        print('Rendered',name,flush=True)
    save(output/'manifest.json',clean(dict(created_utc=now(),
        status=('synthetic figure-layout test; no physics predictions' if data.scope=='synthetic_plot_layout_only' else
            'main figure drafts rendered; numerical convergence and paper selection require inspection'),
        data_scope=data.scope,banner=data.banner,
        inputs=data.inputs,assets=assets,source_hashes={str(p):digest(p) for p in
            (Path(__file__).resolve(),ROOT/'Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py')},
        scale_points={band:sorted(scales.GROUPS[band]) for band in BANDS},
        conventions='Scale bands are recomputed from each plotted observable at the selected weight coordinates. '
            'MC bars come from the full-run reduction. Density conversion divides bin integrals and their '
            'MC errors by the bin width. Normalized shapes retain their parent-rate denominator. '
            'PDF uncertainty is separate and is not included in these bands.',
        limitations='These figures do not by themselves certify bin precision, retraining coverage or '
            'the interpretation of a resolved shift. The complete publication still requires validation '
            'and W/mass robustness figures, tables, variations, and the specification audit.')))


def run(path,output):
    if output.exists():raise ValueError('Refuse to overwrite a main figure directory')
    render(MainData(path),output)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args()
    run(args.report,args.output_dir)
