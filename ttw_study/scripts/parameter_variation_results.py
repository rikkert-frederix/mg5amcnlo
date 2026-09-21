#!/usr/bin/env python3
"""Audited member-0 parameter scans with all 81 signed physical scale points.

The main 101-member reader remains unchanged. The optional PDF IDs accepted
by vectors are only an explicit projection interface for existing benchmark
controls; a parameter-scan audit always requires zero PDF-replica columns.
"""
import argparse
import io
import json
import math
from pathlib import Path
import re

import numpy as np

from campaign import ROOT, STUDY, GF, MW, MZ, digest, now, save
from load_results import POINTS, WEIGHTS, parameter_signature
from parameter_variation_inputs import validate
from stage_batch_reader import read as read_stage_batches

SCALE_WEIGHTS = WEIGHTS[:82]
NUMBER = r'[0-9.eE+\-]+'


def columns(labels, expected_pdf_ids=()):
    """Order top before antitop by signed label, independent of file order."""
    labels = list(labels)
    if len(set(labels)) != len(labels) or labels.count('central value') != 1 or labels.count('dy') != 1:
        raise ValueError('Duplicate labels or missing nominal/error column')
    if len(set(expected_pdf_ids)) != len(expected_pdf_ids):
        raise ValueError('Duplicate declared PDF member')
    scales, pdfs = {}, set()
    for i, label in enumerate(labels):
        if label in ('central value', 'dy') or re.fullmatch(r'delta_mu_(?:cen|min|max) -?\d+ @aux', label):
            continue
        if expected_pdf_ids and re.fullmatch(r'delta_pdf_(?:cen|min|max) [A-Za-z0-9_.\-]+ @aux', label):
            continue
        pdf = re.fullmatch(r'PDF=\s*(\d+)(?: [A-Za-z0-9_.\-]+)?', label)
        if pdf:
            member = int(pdf[1])
            if member in pdfs:
                raise ValueError('Duplicate PDF member')
            pdfs.add(member)
            continue
        pattern = r'dyn=3\s+muR=\s*('+NUMBER+r')\s+muF=\s*('+NUMBER+r')((?:\s+d-?6=\s*'+NUMBER+r'){2})'
        match = re.fullmatch(pattern, label)
        if not match:
            raise ValueError('Unknown parameter-scan weight: '+label)
        decay = re.findall(r'd(-?6)=\s*('+NUMBER+r')', match[3])
        if {int(pdg) for pdg, _ in decay} != {6, -6}:
            raise ValueError('Require separate top and antitop scale axes')
        decay = {int(pdg): float(value) for pdg, value in decay}
        point = (float(match[1]), float(match[2]), decay[6], decay[-6])
        if point not in POINTS or point in scales:
            raise ValueError('Unsupported or duplicate physical scale point')
        scales[point] = i
    if pdfs != set(expected_pdf_ids):
        raise ValueError('Unexpected or incomplete declared PDF columns')
    if set(scales) != set(POINTS):
        raise ValueError('Require all 81 independent signed scale points')
    return [labels.index('central value')]+[scales[point] for point in POINTS]


def vectors(raw, expected_pdf_ids=()):
    from madgraph.various.histograms import HwUList
    histograms = HwUList(io.StringIO(raw.decode()), raw_labels=True)
    if len(histograms) != 210 or not histograms[0].bins:
        raise ValueError('Require all 210 benchmark measurement histograms')
    labels = list(histograms[0].bins[0].wgts)
    selected = columns(labels, expected_pdf_ids)
    bins = [row for hist in histograms for row in hist.bins]
    if any(list(row.wgts) != labels for row in bins):
        raise ValueError('Histogram weight labels differ between bins')
    all_values = np.asarray([[row.wgts[label] for label in labels] for row in bins])
    edges = np.asarray([row.boundaries for row in bins])
    if not np.isfinite(all_values).all() or not np.isfinite(edges).all() or np.any(edges[:, 1] <= edges[:, 0]):
        raise ValueError('Nonfinite histogram or invalid bin boundaries')
    values, errors = all_values[:, selected], all_values[:, labels.index('dy')]
    if np.any(errors < 0.):
        raise ValueError('Negative nominal MC error')
    central = 1+POINTS.index((1., 1., 1., 1.))
    if not np.allclose(values[:, 0], values[:, central], rtol=3.e-6, atol=1.e-20):
        raise ValueError('Nominal differs from the central scale-coordinate weight')
    return dict(values=values, errors=errors, titles=[hist.title for hist in histograms],
        offsets=np.r_[0, np.cumsum([len(hist.bins) for hist in histograms])], edges=edges,
        weights=list(SCALE_WEIGHTS))


def decay_settings(path):
    result = {}
    for line in path.read_text().splitlines():
        line = line.split('!', 1)[0].strip()
        if not line or line.startswith('#'):
            continue
        if '=' not in line:
            raise ValueError('Malformed decay-card entry')
        value, key = (part.strip() for part in line.split('=', 1))
        if key in result:
            raise ValueError('Duplicate decay-card entry')
        result[key] = value
    return result


def verify_inputs(archive, manifest):
    from models.check_param_card import ParamCard
    from madgraph.various.banner import RunCardNLO
    path = Path(manifest['parameter_variation_inputs'])
    if digest(path) != manifest['parameter_variation_inputs_sha256']:
        raise ValueError('Parameter-scan input record changed')
    record = json.loads(path.read_text())
    benchmark, _ = validate(record)
    scenario = manifest['parameter_scenario']
    if scenario not in record['scenarios'] or manifest['variant'] not in ('S', 'Pi'):
        raise ValueError('Unverified parameter scenario or perturbative prescription')
    for name, sha in manifest['hashes'].items():
        if name.endswith('.dat') and digest(archive/name) != sha:
            raise ValueError('Archived scan card changed: '+name)
    if manifest['hashes']['analysis_source'] != digest(ROOT/'Template/fNLO/FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90'):
        raise ValueError('Parameter scan uses a different measurement')
    mode = manifest['w_treatment']
    if mode not in ('onshell', 'top-bw', 'all-bw'):
        raise ValueError('Unsupported W treatment')
    width = next(row for row in scenario['top_widths'] if row['scale_factor'] == 1.
                 and row['w_treatment'] == ('onshell' if mode == 'onshell' else 'bw'))
    expected = dict(top_width_lo=width['gamma_lo'], top_width_nlo=width['gamma_nlo_pdf'],
        top_width_reference_scale=scenario['mt_GeV'], top_width_w_treatment=width['w_treatment'],
        top_width_bottom_mass=0., production_bottom_mass=0., decay_bottom_mass=0.,
        w_width=scenario['fixed_W_width_GeV'], production_scale='core-w-ht-half',
        production_scale_grouping='W_SYSTEM', decay_scale_grouping='separate', pdf_error_weights=False)
    if any(manifest.get(key) != value for key, value in expected.items()):
        raise ValueError('Inconsistent parameter-scan widths, scales or mass convention')
    param = ParamCard(str(archive/'param_card.dat'))
    alpha = next(row['alpha_s'] for row in scenario['coupling_checkpoints'] if row['Q_GeV'] == MZ)
    physical = {('mass', 6): scenario['mt_GeV'], ('yukawa', 6): scenario['mt_GeV'],
        ('mass', 5): 0., ('mass', 23): MZ, ('mass', 24): MW,
        ('sminputs', 1): benchmark['inputs']['alpha_Gmu_inverse'], ('sminputs', 2): GF,
        ('sminputs', 3): alpha, ('decay', 6): width['gamma_nlo_pdf'],
        ('decay', 24): scenario['fixed_W_width_GeV']}
    if 'decaymass' in param:
        raise ValueError('Unexpected massive decay model in the parameter scan')
    for (block, code), value in physical.items():
        if not math.isclose(float(param[block].get((code,)).value), value, rel_tol=2.e-13, abs_tol=1.e-14):
            raise ValueError('Actual parameter card differs from the scan: %s(%d)' % (block, code))
    settings = manifest['settings']
    required = dict(lhaid=scenario['PDF']['id'], pdlabel='lhapdf', reweight_pdf=False,
        reweight_scale=True, rw_rscale=[1., .5, 2.], rw_fscale=[1., .5, 2.],
        ebeam1=scenario['beam_energy_GeV'], ebeam2=scenario['beam_energy_GeV'],
        fixed_ren_scale=False, fixed_fac_scale=False, fixed_qes_scale=False,
        mur_over_ref=1., muf_over_ref=1., qes_over_ref=1., dynamical_scale_choice=3, maxjetflavor=5)
    if any(settings.get(key) != value for key, value in required.items()):
        raise ValueError('Wrong beam, coupling or physical scale settings')
    run = RunCardNLO(str(archive/'run_card.dat'))
    for key, expected_value in settings.items():
        actual = run[key]
        if isinstance(actual, list) and len(actual) == 1 and not isinstance(expected_value, list):
            actual = actual[0]
        if actual != expected_value:
            raise ValueError('Actual run card differs from archived settings: '+key)
    decay = decay_settings(archive/'decay_card.dat')
    required_decay = dict(production_order='NLO', decay_order='NLO',
        nlo_decay_combination='ADDITIVE' if manifest['variant'] == 'S' else 'MULTIPLICATIVE',
        decay_scale_variation_mode='INDEPENDENT', decay_scale_grouping='SIGNED_PDG',
        production_ren_scale_momenta='CORE', production_scale_grouping='W_SYSTEM')
    required_decay.update({'decay_width_scale_mode(6)': 'AUTO', 'decay_dynamical_scale_choice(6)': '0'})
    if any(decay.get(key) != value for key, value in required_decay.items()):
        raise ValueError('Unexpected actual decay prescription or scale convention')
    for key, value in {'lo_decay_width(6)': width['gamma_lo'], 'nlo_decay_width(6)': width['gamma_nlo_pdf'],
                       'decay_ren_scale(6)': scenario['mt_GeV']}.items():
        if float(decay[key]) != value:
            raise ValueError('Actual decay width or reference differs from matched inputs')
    if tuple(float(v) for v in decay['decay_scale_factors'].split(',')) != (1., .5, 2.):
        raise ValueError('Unexpected signed decay scale-factor list')
    if ('lo_decay_width(24)' in decay) != (mode != 'all-bw'):
        raise ValueError('Wrong number of W branching denominators')
    if mode != 'all-bw' and float(decay['lo_decay_width(24)']) != scenario['fixed_W_width_GeV']:
        raise ValueError('Associated-W branching denominator changed')
    return scenario


def audit(path, output):
    path, output = Path(path).resolve(), Path(output).resolve()
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a parameter-scan audit')
    process, name = path.parents[2], path.parent.name
    archive = process/'study_cards'/name
    manifest = json.loads((archive/'manifest.json').read_text())
    execution = json.loads((archive/'execution.json').read_text())
    if (manifest['run_name'] != name or execution['status'] != 'finished'
            or execution['outputs'].get(str(path)) != digest(path)):
        raise ValueError('Missing or changed completed parameter integration')
    verify_inputs(archive, manifest)
    if digest(process/'FixedOrderAnalysis/HwU.f90') != digest(ROOT/'Template/fNLO/FixedOrderAnalysis/HwU.f90'):
        raise ValueError('Unverified histogram accumulator')
    result = vectors(path.read_bytes())
    rates, partition_max, strict_max = {}, 0., 0.
    for index, title in enumerate(result['titles']):
        if ' rates:' not in title:
            continue
        lo, hi = result['offsets'][index:index+2]
        row = result['values'][lo:hi]
        if len(row) != 9:
            raise ValueError('Incomplete rate bank')
        delta = np.abs(row[4]-row[5:9].sum(axis=0))
        if not (delta <= 3.e-6*(np.abs(row[4])+np.abs(row[5:9]).sum(axis=0))+1.e-20).all():
            raise ValueError('Extra-jet sectors do not partition the two-b rate')
        partition_max = max(partition_max, float(delta.max()))
        if manifest['variant'] == 'S':
            strict_max = max(strict_max, float(np.abs(row[7:9]).max()))
        rates[title] = dict(nominal_pb=row[:, 0].tolist(), mc_error_pb=result['errors'][lo:hi].tolist())
    if len(rates) != 10 or strict_max != 0.:
        raise ValueError('Missing rate banks or forbidden strict-NLO extra jets')
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'), **result)
    report = dict(created_utc=now(), path=str(path), output_sha256=digest(path),
        manifest=manifest, execution=execution, variant=manifest['variant'],
        status='parameter-scan histogram checks passed; precision not certified',
        corrected_accumulator=True, histogram_count=210, scale_count=81, pdf_count=0,
        maximum_jet_partition_mismatch_pb=partition_max, strict_2plus_extra_maximum_pb=strict_max,
        rates=rates, arrays_sha256=digest(output.with_suffix('.npz')), script_sha256=digest(Path(__file__)))
    save(output, report)
    return report


def load(path):
    path = Path(path).resolve()
    report = json.loads(path.read_text())
    if (report['status'] != 'parameter-scan histogram checks passed; precision not certified'
            or not report['corrected_accumulator'] or report['pdf_count'] != 0 or report['scale_count'] != 81
            or digest(report['path']) != report['output_sha256']
            or digest(path.with_suffix('.npz')) != report['arrays_sha256']):
        raise ValueError('Changed or incomplete parameter-scan histogram audit')
    hwu = Path(report['path'])
    archive = hwu.parents[2]/'study_cards'/hwu.parent.name
    verify_inputs(archive, report['manifest'])
    with np.load(path.with_suffix('.npz')) as arrays:
        result = {key: arrays[key].copy() for key in ('values', 'errors', 'offsets', 'edges')}
        result.update(titles=arrays['titles'].tolist(), weights=arrays['weights'].tolist())
    if result['weights'] != SCALE_WEIGHTS or not np.isfinite(result['values']).all():
        raise ValueError('Invalid canonical scan weights')
    result.update(report=report, parameter_signature=parameter_signature(archive/'param_card.dat'))
    return result


def read_batches(record_path, output):
    record = json.loads(Path(record_path).read_text())
    reference = load(record['audit'])
    read_stage_batches(record_path, output, reference, vectors)


def load_batches(path):
    path = Path(path).resolve()
    metadata = json.loads(path.with_suffix('.json').read_text())
    if (not metadata['status'].startswith('equal-count conditional')
            or digest(path) != metadata['arrays_sha256']
            or digest(metadata['input']) != metadata['input_sha256']):
        raise ValueError('Changed or unsuccessful scan batch audit')
    record = json.loads(Path(metadata['input']).read_text())
    if digest(record['audit']) != record['audit_sha256']:
        raise ValueError('Scan final-output audit changed')
    result = load(record['audit'])
    with np.load(path) as data:
        if data['weights'].tolist() != result['weights']:
            raise ValueError('Batch and combined scale-coordinate ordering differ')
        result.update(contributions=data['contributions'].copy(), strata=data['strata'].copy(),
                      seeds=data['random_seeds'].copy())
    result['batch_audit'] = metadata
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    for action in ('audit', 'batches'):
        command = sub.add_parser(action)
        command.add_argument('input', type=Path)
        command.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    (audit if args.action == 'audit' else read_batches)(args.input, args.output)
