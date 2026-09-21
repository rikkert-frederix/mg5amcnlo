#!/usr/bin/env python3
"""Verify bounded parameter-scan inputs without launching integrations.

The frozen benchmark and original variation preflight remain unchanged.
Every scenario contains matched on-shell/BW widths at all three decay scales.
PDF-family and alpha-s comparisons use member 0; PDF errors belong to the
separate, fully correlated main replica calculation.
"""
import argparse
import copy
import json
import math
from pathlib import Path

from campaign import ROOT, STUDY, WIDTH, MT, MW, MZ, digest, now, save

SCENARIOS = (
    ('baseline', 'baseline', 'NNPDF40_nlo_as_01180', 331700, MT, 13000.),
    ('ct18', 'pdf_family', 'CT18NLO', 14400, MT, 13000.),
    ('msht20', 'pdf_family', 'MSHT20nlo_as118', 27100, MT, 13000.),
    ('as117', 'alpha_s', 'NNPDF40_nlo_as_01170', 333900, MT, 13000.),
    ('as119', 'alpha_s', 'NNPDF40_nlo_as_01190', 334100, MT, 13000.),
    ('mt171p5', 'top_mass', 'NNPDF40_nlo_as_01180', 331700, 171.5, 13000.),
    ('mt173p5', 'top_mass', 'NNPDF40_nlo_as_01180', 331700, 173.5, 13000.),
    ('ecm13600', 'beam_energy', 'NNPDF40_nlo_as_01180', 331700, MT, 13600.),
)
FACTORS = (.5, 1., 2.)
WIDTH_MODES = ('onshell', 'bw')


def close(a, b):
    return math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=2.e-12, abs_tol=1.e-14)


def source_rows(definition, benchmark, preflight):
    name, family, pdf_name, pdf_id, mt, energy = definition
    if family in ('pdf_family', 'alpha_s'):
        matches = [row for row in preflight['PDFs'] if row['name'] == pdf_name and row['id'] == pdf_id]
        if len(matches) != 1:
            raise ValueError('Missing or ambiguous installed alternative PDF: '+name)
        rows = matches[0]['top_widths']
        source = 'parameter_variations_preflight.PDFs.'+pdf_name+'.top_widths'
    elif family == 'top_mass':
        rows = [row for row in preflight['top_mass_scan'] if row['mt_GeV'] == mt]
        source = 'parameter_variations_preflight.top_mass_scan.mt='+str(mt)
    else:
        rows = benchmark['top_widths']
        source = 'benchmark.top_widths'
    return [copy.deepcopy(row) for row in rows if row['mb_GeV'] == 0.], source


def validate_scenario(scenario, definition, benchmark, preflight):
    name, family, pdf_name, pdf_id, mt, energy = definition
    expected = dict(name=name, family=family, mt_GeV=mt, top_yukawa_GeV=mt,
                    ecm_GeV=energy, beam_energy_GeV=energy/2., mw_GeV=MW,
                    production_mb_GeV=0., decay_mb_GeV=0., alpha_s_flavours=5,
                    fixed_W_width_GeV=benchmark['W_width']['gamma_nlo_pdf'],
                    top_width_reference_scale_GeV=mt, decay_scale_reference_GeV=mt,
                    production_scale='core-w-ht-half', qes_over_ref=1.)
    if any(scenario.get(key) != value for key, value in expected.items()):
        raise ValueError('Unmatched parameter-scan physical inputs: '+name)
    if (scenario['PDF']['name'], scenario['PDF']['id'], scenario['PDF']['member']) != (pdf_name, pdf_id, 0):
        raise ValueError('Wrong PDF family or member: '+name)
    rows, source = source_rows(definition, benchmark, preflight)
    if scenario['width_source'] != source or scenario['top_widths'] != rows:
        raise ValueError('Width inputs differ from their frozen calculation: '+name)
    coordinates = {(row['w_treatment'], row['scale_factor']) for row in rows}
    if len(rows) != 6 or coordinates != {(w, f) for w in WIDTH_MODES for f in FACTORS}:
        raise ValueError('Incomplete or repeated decay-width coordinates: '+name)
    alphas = {row['Q_GeV']: row['alpha_s'] for row in scenario['coupling_checkpoints']}
    if len(alphas) != 5 or set(alphas) != {MW, mt/2., MZ, mt, 2.*mt}:
        raise ValueError('Incomplete PDF-coupling checkpoints: '+name)
    if not all(math.isfinite(a) and 0. < a < 1. for a in alphas.values()):
        raise ValueError('Invalid PDF coupling: '+name)
    for row in rows:
        factor, mode = row['scale_factor'], row['w_treatment']
        if (row['muR_GeV'] != mt*factor or not close(row['alpha_s_pdf'], alphas[mt*factor])
                or not close(row['gamma_nlo_pdf'], row['gamma_lo']+row['qcd_coefficient']*row['alpha_s_pdf'])
                or min(row['gamma_lo'], row['gamma_nlo_pdf']) <= 0.):
            raise ValueError('Wrong scale or coefficient-level width matching: '+name)
        if family in ('pdf_family', 'alpha_s'):
            base = next(r for r in benchmark['top_widths'] if r['mb_GeV'] == 0.
                        and r['scale_factor'] == factor and r['w_treatment'] == mode)
            if any(row[key] != base[key] for key in ('gamma_lo', 'qcd_coefficient')):
                raise ValueError('A PDF change modified the coupling-independent width coefficient')
        if family == 'top_mass':
            command = row['command']
            ww = 0. if mode == 'onshell' else benchmark['W_width']['gamma_nlo_pdf']
            if [float(v) for v in command[1:]] != [mt, MW, 0., mt*factor, ww, .118, 1., 5.]:
                raise ValueError('Top-mass width calculator inputs differ from the scenario')
            values = row['values']
            if (not close(row['gamma_lo'], values['Gamma_LO [GeV]'])
                    or not close(row['qcd_coefficient'], values['Delta_Gamma_QCD^NLO [GeV]']/values['alpha_s(muR)'])):
                raise ValueError('Top-mass width coefficient differs from calculator output')
    for mode in WIDTH_MODES:
        selected = [row for row in rows if row['w_treatment'] == mode]
        for key in ('gamma_lo', 'qcd_coefficient'):
            if not all(close(row[key], selected[0][key]) for row in selected):
                raise ValueError('NLO width running cannot be reproduced from one coefficient')


def pdf_record(definition, benchmark, preflight):
    _, _, name, pdf_id, _, _ = definition
    if pdf_id == benchmark['PDF']['id']:
        base = benchmark['PDF']
        member_hash = next(row['sha256'] for row in base['members'] if row['member'] == 0)
    else:
        base = next(row for row in preflight['PDFs'] if row['id'] == pdf_id and row['name'] == name)
        member_hash = base['member_sha256'][name+'_0000.dat']
    directory = ROOT/'LHAPDF/share/LHAPDF'/name
    files = {str(directory/(name+'.info')): base['metadata_sha256'],
             str(directory/(name+'_0000.dat')): member_hash}
    if any(digest(path) != sha for path, sha in files.items()):
        raise ValueError('Installed PDF member 0 or metadata changed: '+name)
    return dict(name=name, id=pdf_id, member=0, DataVersion=base['DataVersion'], files=files)


def validate(record):
    benchmark_path = STUDY/'inputs/benchmark.json'
    preflight_path = STUDY/'inputs/parameter_variations_preflight.json'
    if (record['benchmark_sha256'] != digest(benchmark_path)
            or record['preflight_sha256'] != digest(preflight_path)):
        raise ValueError('Frozen physical inputs or variation preflight changed')
    benchmark = json.loads(benchmark_path.read_text())
    preflight = json.loads(preflight_path.read_text())
    if preflight['benchmark_sha256'] != digest(benchmark_path):
        raise ValueError('Variation preflight belongs to a different benchmark')
    if (len(record['scenarios']) != len(SCENARIOS)
            or [row['name'] for row in record['scenarios']] != [row[0] for row in SCENARIOS]):
        raise ValueError('Missing, repeated or reordered parameter scenarios')
    for scenario, definition in zip(record['scenarios'], SCENARIOS):
        validate_scenario(scenario, definition, benchmark, preflight)
        if scenario['PDF'] != pdf_record(definition, benchmark, preflight):
            raise ValueError('PDF provenance differs from installed frozen member 0')
    return benchmark, preflight


def run(output):
    import lhapdf
    if output.exists():
        raise ValueError('Refuse to overwrite parameter-variation inputs')
    benchmark_path = STUDY/'inputs/benchmark.json'
    preflight_path = STUDY/'inputs/parameter_variations_preflight.json'
    benchmark = json.loads(benchmark_path.read_text())
    preflight = json.loads(preflight_path.read_text())
    if preflight['benchmark_sha256'] != digest(benchmark_path):
        raise ValueError('Variation preflight belongs to a different benchmark')
    if any(digest(WIDTH/name) != sha for name, sha in benchmark['width_sources'].items()):
        raise ValueError('Width-calculator source differs from the archived calculations')
    lhapdf.setVerbosity(0)
    scenarios = []
    for definition in SCENARIOS:
        name, family, pdf_name, pdf_id, mt, energy = definition
        pdf_info = pdf_record(definition, benchmark, preflight)
        pdf = lhapdf.mkPDF(pdf_name, 0)
        rows, source = source_rows(definition, benchmark, preflight)
        scenario = dict(name=name, family=family, mt_GeV=mt, top_yukawa_GeV=mt,
            ecm_GeV=energy, beam_energy_GeV=energy/2., mw_GeV=MW,
            production_mb_GeV=0., decay_mb_GeV=0., alpha_s_flavours=5,
            fixed_W_width_GeV=benchmark['W_width']['gamma_nlo_pdf'],
            top_width_reference_scale_GeV=mt, decay_scale_reference_GeV=mt,
            production_scale='core-w-ht-half', qes_over_ref=1., PDF=pdf_info,
            top_widths=rows, width_source=source,
            coupling_checkpoints=[dict(Q_GeV=q, alpha_s=pdf.alphasQ(q))
                                  for q in (MW, mt/2., MZ, mt, 2.*mt)])
        validate_scenario(scenario, definition, benchmark, preflight)
        scenarios.append(scenario)
    record = dict(created_utc=now(), status='parameter inputs verified; no variation integrations',
        benchmark_sha256=digest(benchmark_path), preflight_sha256=digest(preflight_path),
        scenarios=scenarios, script_sha256=digest(Path(__file__)),
        width_sources=benchmark['width_sources'], width_executable_sha256=digest(WIDTH/'top_decay_width'),
        lhapdf_version=lhapdf.version(),
        scope='Baseline plus seven bounded member-0 PDF-family, alpha-s, top-mass and beam-energy scenarios. '
              'Matched massless-decay widths support on-shell, top-BW and all-BW exports.',
        conventions='W width and weak inputs remain fixed. Top mass, Yukawa, phase-space mass, '
                    'decay-scale reference and width coefficient move together. Every PDF or alpha-s '
                    'scenario requires a fresh numerator integration with the same coupling as its widths.',
        remaining='Choose S/Pi allocation and observables after main S/Pi convergence. Sum actual ordered '
                  'flavour cross sections for any inclusive e/mu prediction. These inputs are not MC results.')
    validate(record)
    save(output, record)
    print('Verified eight scenarios and 48 matched width rows; no integrations launched.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output)
