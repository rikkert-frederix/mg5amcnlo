#!/usr/bin/env python3
"""Matched fixed-common-scale inputs for a generated narrow-W limit test.

Only the external W width is narrowed; weak couplings and top masses stay
fixed. Each selected leptonic W density therefore has a 1/epsilon pole
normalization. Compare epsilon^3 times the three-W result to the physical-
width on-shell reference. This is a technical limit, not a width uncertainty.
"""
import argparse
import json
import math
from pathlib import Path

from bw_branching import partial_width
from campaign import MT, MW, STUDY, WIDTH, calculator, digest, now, save

COMMON_SCALE=MT+MW/2.
FACTORS=(1.,.2,.05)


def rescaling(factor):
    if not math.isfinite(factor) or not 0.<factor<=1.:
        raise ValueError('Narrow-W factor must be in (0,1]')
    return factor**3


def width_row(mass,factor,onshell,alpha):
    benchmark=json.loads((STUDY/'inputs/benchmark.json').read_text())
    physical=benchmark['W_width']['gamma_nlo_pdf']
    ww=physical*factor
    row=calculator([WIDTH/'top_decay_width',MT,MW,mass,COMMON_SCALE,
                    0. if onshell else ww,.118,1.,5])
    values=row['values']
    lo=values['Gamma_LO [GeV]']
    coefficient=values['Delta_Gamma_QCD^NLO [GeV]']/values['alpha_s(muR)']
    row.update(decay_bottom_mass_GeV=mass,width_factor=factor,w_width_GeV=ww,
               width_mode='onshell' if onshell else 'bw',mu_reference_GeV=COMMON_SCALE,
               gamma_lo=lo,qcd_coefficient=coefficient,alpha_s_pdf=alpha,gamma_nlo=lo+coefficient*alpha)
    if not all(math.isfinite(row[k]) and row[k]>0. for k in ('gamma_lo','gamma_nlo','w_width_GeV')):
        raise ValueError('Nonpositive/nonfinite matched width')
    return row


def run(output):
    if output.exists():
        raise ValueError('Refuse to overwrite narrow-W inputs')
    import lhapdf
    lhapdf.setVerbosity(0)
    pdf=lhapdf.mkPDF('NNPDF40_nlo_as_01180',0)
    alpha=pdf.alphasQ(COMMON_SCALE)
    benchmark_path=STUDY/'inputs/benchmark.json'
    benchmark=json.loads(benchmark_path.read_text())
    for name,sha in benchmark['width_sources'].items():
        if digest(WIDTH/name)!=sha:
            raise ValueError('Width-calculator source changed since the benchmark')
    baseline=benchmark['W_width']['branching_e']
    rows=[]
    for mass in (0.,4.8):
        on=width_row(mass,1.,True,alpha)
        rows.append(on)
        for factor in FACTORS:
            row=width_row(mass,factor,False,alpha)
            partial,error=partial_width(MT,MW,mass,row['w_width_GeV'])
            scaled_branching=factor*partial/row['gamma_lo']
            if not math.isclose(scaled_branching,baseline,rel_tol=2.e-10,abs_tol=0.):
                raise ValueError('Fixed-coupling current does not have the assumed pole normalization')
            row.update(leptonic_partial_LO_GeV=partial,quadrature_error_GeV=error,
                       rescaled_top_branching=scaled_branching,
                       relative_LO_width_change=row['gamma_lo']/on['gamma_lo']-1.,
                       relative_NLO_width_change=row['gamma_nlo']/on['gamma_nlo']-1.,
                       cross_section_rescaling=rescaling(factor))
            rows.append(row)
        selected=[row for row in rows if row['decay_bottom_mass_GeV']==mass and row['width_mode']=='bw']
        for key in ('relative_LO_width_change','relative_NLO_width_change'):
            shifts=[abs(row[key]) for row in selected]
            if not all(b<a for a,b in zip(shifts,shifts[1:])) or shifts[-1]>.001:
                raise ValueError('Matched width narrowing does not approach the on-shell reference')
    record=dict(created_utc=now(),status='matched common-scale narrow-W inputs and pole normalization checked; no MC',
                benchmark_sha256=digest(benchmark_path),common_scale_GeV=COMMON_SCALE,
                alpha_s_pdf=alpha,scale_checkpoints={str(f):pdf.alphasQ(f*COMMON_SCALE) for f in (.5,1.,2.)},
                width_factors=FACTORS,rows=rows,script_sha256=digest(Path(__file__)),
                width_executable_sha256=digest(WIDTH/'top_decay_width'),width_sources=benchmark['width_sources'],
                branching_reference=baseline,
                convention='Fixed weak inputs, no top off-shellness. Compare epsilon^3 times all-bw/top-bw '
                           'cross sections to an on-shell result with physical GammaW. This is NOT an '
                           'additional W branching factor and does not modify generated event weights. '
                           'Top total LO/NLO widths are recomputed with the same W prescription and '
                           'five-flavour PDF coupling at common muR=muF=QES=mu_t=mu_antitop=mt+MW/2.',
                source='https://arxiv.org/pdf/1204.1513, Eqs. (2.8)-(2.9); '
                       'epsilon^3 normalization is our fixed-coupling three-leptonic-current derivation',
                limitations='Independent deterministic current/width checks only. Full-support generated '
                            'LO/S/Pi rates and shapes, statistical narrowing convergence, and massive '
                            'small-mass continuity are separate outstanding tests.')
    save(output,record)
    print('Prepared',len(rows),'common-scale width rows; no generated narrow-W limit claimed.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.output)
