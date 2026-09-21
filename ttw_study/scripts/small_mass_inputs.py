#!/usr/bin/env python3
"""Matched total widths and current normalization for generated mb -> 0 tests."""
import argparse
import copy
import json
import math
from pathlib import Path

from bw_branching import partial_width
from campaign import MT, MW, STUDY, WIDTH, calculator, digest, now, save

MASSES=(0.,1.,.1,4.8)
FACTORS=(.5,1.,2.)


def width_row(mass,mode,factor,alpha,benchmark,scans):
    """Reuse unchanged calculator outputs at mt; calculate only missing rows."""
    known=[row for row in benchmark['top_widths'] if row['mb_GeV']==mass
           and row['w_treatment']==mode and row['scale_factor']==factor]
    source='frozen physical benchmark'
    if not known and factor==1.:
        ww=0. if mode=='onshell' else benchmark['W_width']['gamma_nlo_pdf']
        known=[row for row in scans['small_bottom_mass']
               if row['mb_GeV']==mass and row['W_width_GeV']==ww]
        source='archived width-only continuity scan at mt'
    if len(known)>1:
        raise ValueError('Ambiguous archived small-mass width')
    if known:
        row=copy.deepcopy(known[0])
    else:
        ww=0. if mode=='onshell' else benchmark['W_width']['gamma_nlo_pdf']
        row=calculator([WIDTH/'top_decay_width',MT,MW,mass,MT*factor,ww,.118,1.,5])
        values=row['values']
        row.update(gamma_lo=values['Gamma_LO [GeV]'],
                   qcd_coefficient=values['Delta_Gamma_QCD^NLO [GeV]']/values['alpha_s(muR)'])
        source='new width-calculator evaluation at the requested scale'
    row.update(mb_GeV=mass,w_treatment=mode,scale_factor=factor,muR_GeV=MT*factor,
               alpha_s_pdf=alpha,gamma_nlo_pdf=row['gamma_lo']+row['qcd_coefficient']*alpha,
               reuse_source=source)
    if not all(math.isfinite(row[k]) for k in ('gamma_lo','qcd_coefficient','gamma_nlo_pdf')):
        raise ValueError('Nonfinite matched small-mass width')
    if min(row['gamma_lo'],row['gamma_nlo_pdf'])<=0.:
        raise ValueError('Nonpositive matched total width')
    return row


def validate(record,benchmark):
    expected={(m,w,f) for m in MASSES for w in ('onshell','bw') for f in FACTORS}
    rows=record['top_widths']
    if len(rows)!=len(expected) or {(r['mb_GeV'],r['w_treatment'],r['scale_factor']) for r in rows}!=expected:
        raise ValueError('Incomplete or repeated small-mass width coordinates')
    for row in rows:
        if (row['muR_GeV']!=MT*row['scale_factor'] or
                not math.isclose(row['gamma_nlo_pdf'],row['gamma_lo']+
                                 row['qcd_coefficient']*row['alpha_s_pdf'],rel_tol=1.e-13)):
            raise ValueError('Wrong scale or coefficient-level coupling matching')
    for mode in ('onshell','bw'):
        for factor in FACTORS:
            selected={r['mb_GeV']:r for r in rows if r['w_treatment']==mode and r['scale_factor']==factor}
            for key in ('gamma_lo','gamma_nlo_pdf'):
                shifts={m:abs(selected[m][key]/selected[0.][key]-1.) for m in (4.8,1.,.1)}
                if not shifts[.1]<shifts[1.]<shifts[4.8] or shifts[.1]>1.e-5:
                    raise ValueError('Width-only small-mass continuity failed')
    if record['fixed_W_width_GeV']!=benchmark['W_width']['gamma_nlo_pdf']:
        raise ValueError('The physical W width must stay fixed')


def run(output):
    import lhapdf
    if output.exists():
        raise ValueError('Refuse to overwrite matched small-mass inputs')
    base_path=STUDY/'inputs/benchmark.json'
    scan_path=STUDY/'widths/continuity_scans.json'
    benchmark=json.loads(base_path.read_text())
    scans=json.loads(scan_path.read_text())
    if (scans['benchmark_sha256']!=digest(base_path)
            or scans['executable_sha256']!=digest(WIDTH/'top_decay_width')
            or any(digest(WIDTH/p)!=sha for p,sha in benchmark['width_sources'].items())):
        raise ValueError('Archived calculator or physical benchmark changed')
    lhapdf.setVerbosity(0)
    pdf=lhapdf.mkPDF('NNPDF40_nlo_as_01180',0)
    rows=[width_row(m,w,f,pdf.alphasQ(MT*f),benchmark,scans)
          for m in MASSES for w in ('onshell','bw') for f in FACTORS]
    branching=benchmark['W_width']['branching_e']
    ww=benchmark['W_width']['gamma_nlo_pdf']
    current=[]
    for mass in MASSES:
        partial,error=partial_width(MT,MW,mass,ww)
        total=next(r['gamma_lo'] for r in rows if r['mb_GeV']==mass
                   and r['w_treatment']=='bw' and r['scale_factor']==1.)
        if not math.isclose(partial/total,branching,rel_tol=2.e-10):
            raise ValueError('Matched small-mass current normalization failed')
        current.append(dict(mb_GeV=mass,partial_LO_GeV=partial,quadrature_error_GeV=error,
                            total_LO_GeV=total,integrated_branching=partial/total))
    record=dict(created_utc=now(),status='small-mass width/current inputs checked; no generated continuity result',
        benchmark_sha256=digest(base_path),archived_width_scan=str(scan_path),
        archived_width_scan_sha256=digest(scan_path),width_executable_sha256=digest(WIDTH/'top_decay_width'),
        width_sources=benchmark['width_sources'],top_widths=rows,current_normalization=current,
        masses_GeV=MASSES,scale_factors=FACTORS,fixed_W_width_GeV=ww,
        script_sha256=digest(Path(__file__)),
        convention='Five-flavour PDF coupling for every scale; production mb=0, decay on-shell mb varied. '
                   'Fixed weak inputs, fixed physical W width, matched total LO/NLO widths. '
                   'No additional branching factor on generated BW currents.',
        limitations='Deterministic widths and LO current integrals only. Generated virtual/pole/subtraction '
                    'and IR-safe rate/shape continuity require independent MC with these exact inputs.')
    validate(record,benchmark)
    save(output,record)
    print('Prepared',len(rows),'matched width rows; current normalization verified for',len(current),'masses.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
