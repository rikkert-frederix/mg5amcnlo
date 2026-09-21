#!/usr/bin/env python3
"""Independent width-calculator continuity tests (not matrix-element validation)."""
import json
import math

from campaign import WIDTH, STUDY, MT, MW, calculator, digest, now, save


def main():
    benchmark = json.loads((STUDY/'inputs/benchmark.json').read_text())
    ww = benchmark['W_width']['gamma_nlo_pdf']
    alpha = next(row['alpha_s'] for row in benchmark['coupling_checkpoints']
                 if row['Q_GeV'] == MT)
    path = STUDY/'widths/continuity_scans.json'
    if path.exists():
        raise ValueError('Scan archive already exists')

    def evaluate(mb, width):
        row = calculator([WIDTH/'top_decay_width', MT, MW, mb, MT, width, .118, 1., 5])
        values = row['values']
        g0 = values['Gamma_LO [GeV]']
        coefficient = values['Delta_Gamma_QCD^NLO [GeV]']/values['alpha_s(muR)']
        row.update(mb_GeV=mb, W_width_GeV=width, gamma_lo=g0,
                   qcd_coefficient=coefficient, gamma_nlo_pdf=g0+coefficient*alpha)
        assert all(math.isfinite(row[key]) for key in ('gamma_lo','qcd_coefficient','gamma_nlo_pdf'))
        return row

    narrow, small_mass = [], []
    for mb in (0., 4.8):
        reference = evaluate(mb, 0.)
        for factor in (1., .5, .1, .03, .01, .003):
            row = evaluate(mb, ww*factor)
            row.update(width_factor=factor,
                       relative_LO_change=row['gamma_lo']/reference['gamma_lo']-1.,
                       relative_NLO_change=row['gamma_nlo_pdf']/reference['gamma_nlo_pdf']-1.)
            narrow.append(row)
        for key in ('relative_LO_change','relative_NLO_change'):
            shifts = [abs(row[key]) for row in narrow if row['mb_GeV'] == mb]
            assert all(b < a for a,b in zip(shifts,shifts[1:])), shifts
            assert shifts[-1] < 1.e-4, shifts[-1]
    for width in (0., ww):
        reference = evaluate(0., width)
        for mb in (4.8, 2., 1., .1, .01, .001):
            row = evaluate(mb, width)
            row.update(relative_LO_change=row['gamma_lo']/reference['gamma_lo']-1.,
                       relative_NLO_change=row['gamma_nlo_pdf']/reference['gamma_nlo_pdf']-1.)
            small_mass.append(row)
        assert abs(small_mass[-1]['relative_LO_change']) < 1.e-7
        assert abs(small_mass[-1]['relative_NLO_change']) < 1.e-7
    save(path, dict(created_utc=now(), status='width continuity passed',
                    benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),
                    executable_sha256=digest(WIDTH/'top_decay_width'),
                    alpha_s_pdf=alpha, narrow_W=narrow, small_bottom_mass=small_mass,
                    limitation='These test the independent width calculator only; generated decay matrix elements, subtraction and phase space need separate checks.'))
    print('Passed width continuity: 12 narrow-W and 12 small-bottom-mass rows.')


if __name__ == '__main__':
    main()
