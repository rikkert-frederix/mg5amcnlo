#!/usr/bin/env python3
"""Fixed-coupling leptonic-current normalization of the top-W BW convolution.

This is a deterministic width/current check, not an additional W branching
factor to multiply into generated three-body decay events.
"""
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from campaign import GF, MT, MW, STUDY, WIDTH, calculator, digest, now, save


def shape(s, mt, mb):
    x, beta2 = s/mt**2, (mb/mt)**2
    lam = (1-x-beta2)**2-4*x*beta2
    if s <= 0. or lam <= 0.:
        return 0.
    return math.sqrt(lam)*((1-beta2)**2+x*(1+beta2)-2*x*x)


def current_density(s, mt, mw, mb, gw, gf=GF):
    """dGamma(t -> b l nu)/ds from fixed-g two-body factors, at LO."""
    if s <= 0.:
        return 0.
    # g^2=4 sqrt(2) GF MW^2 is fixed; replacing MW by sqrt(s) in the
    # kinematics does NOT redefine the electroweak coupling or GF.
    top_two_body = gf*mt**3/(8*math.pi*math.sqrt(2))*mw**2/s*shape(s,mt,mb)
    w_leptonic = gf*mw**2*math.sqrt(s)/(6*math.pi*math.sqrt(2))
    return top_two_body*math.sqrt(s)*w_leptonic/(math.pi*((s-mw**2)**2+(mw*gw)**2))


def convolution_density(s, mt, mw, mb, gw, gf=GF):
    kernel = gf*mt**3/(8*math.pi*math.sqrt(2))*shape(s,mt,mb)
    return kernel*mw*gw/(math.pi*((s-mw**2)**2+(mw*gw)**2))


def partial_width(mt, mw, mb, gw):
    a = math.atan(-mw/gw)
    b = math.atan(((mt-mb)**2-mw**2)/(mw*gw))
    def integrand(theta):
        s = mw**2+mw*gw*math.tan(theta)
        jacobian = mw*gw/math.cos(theta)**2
        return current_density(s,mt,mw,mb,gw)*jacobian
    return quad(integrand,a,b,epsabs=2.e-13,epsrel=2.e-12,limit=200)


def main():
    import json
    output = STUDY/'inputs/bw_branching_normalization.json'
    if output.exists():
        raise ValueError('Refuse to overwrite the BW current-normalization check')
    benchmark = json.loads((STUDY/'inputs/benchmark.json').read_text())
    reference_gw = benchmark['W_width']['gamma_nlo_pdf']
    gamma_lep = GF*MW**3/(6*math.pi*math.sqrt(2))
    records = []
    for mb in (0.,4.8):
        for width_factor in (1.,.2,.05):
            gw = reference_gw*width_factor
            branching = gamma_lep/gw
            ratios = [current_density(s,MT,MW,mb,gw)/convolution_density(s,MT,MW,mb,gw)
                      for s in (1.e-5,MW**2/4,MW**2,.99*(MT-mb)**2)]
            if not np.allclose(ratios,branching,rtol=2.e-14,atol=0.):
                raise ValueError('Pointwise fixed-current normalization identity failed')
            value,error = partial_width(MT,MW,mb,gw)
            width = calculator([WIDTH/'top_decay_width',MT,MW,mb,MT,gw,.118,1.,5])
            total_lo = width['values']['Gamma_LO [GeV]']
            residual = value-branching*total_lo
            if abs(residual)>2.e-10*abs(value):
                raise ValueError('Integrated leptonic current disagrees with matched BW width')
            records.append(dict(mb_GeV=mb,W_width_factor=width_factor,W_width_GeV=gw,
                                leptonic_partial_LO_GeV=value,quadrature_error_GeV=error,
                                total_BW_LO_GeV=total_lo,integrated_branching=value/total_lo,
                                expected_branching=branching,residual_GeV=residual,
                                calculator=width))
    save(output,dict(created_utc=now(),status='LO fixed-current pointwise and full-support integral checks passed',
                     benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),
                     script_sha256=digest(Path(__file__)),records=records,
                     source='Campbell-Ellis arXiv:1204.1513 Eqs. 2.2 and 2.8; direct fixed-g current factorization',
                     NLO_argument='Top QCD corrections replace the transverse-current kernel while retaining '
                                  'the same fixed weak-coupling factors. Thus Gamma_leptonic_i=B_l*Gamma_BW_i '
                                  'for i=0,1 with the matched convolution and a fixed external GammaW. '
                                  'This is an algebraic argument, not an independent NLO integration.',
                     scope='Stable-ttW branching comparison extends to top-bw with an on-shell associated W. '
                           'It does not extend to all-bw production. No extra branching factor is applied to events. '
                           'The generator narrow-W and virtuality-support tests remain separate.'))
    print(output,flush=True)


if __name__ == '__main__':
    main()
