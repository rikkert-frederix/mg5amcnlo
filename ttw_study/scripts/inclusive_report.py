#!/usr/bin/env python3
"""Single-assignment branching normalization with an on-shell associated W.

Stable LO supplies LO/D/PiD; stable NLO supplies P/S/Pi. Comparisons are
point-matched before scale/PDF reduction. Only the nominal HwU MC error is
available from a combined stable run; no varied-weight significance is
invented. An optional audited joint-batch input separately tests cancellation
of the decay-scale dependence, retaining correlations between weights.
Top-bw additionally requires the fixed-current/full-support width check;
the absolute stable-production identity never applies to all-bw.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from joint_report import estimate
from load_results import POINTS, WEIGHTS, canonical_columns, load
from pilot_report import clean, describe
from replica_statistics import ratio

STATISTICAL = {'iseed','req_acc_fo','npoints_fo','niters_fo',
               'npoints_fo_grid','niters_fo_grid','fo_job_target_time',
               'fo_job_min_splits'}
STABLE_VARIANT = dict(LO='LO',D='LO',PiD='LO',P='P',S='P',Pi='P')


def production_parameters(path):
    """Compare numeric model inputs, excluding only documented decay inputs."""
    from models.check_param_card import ParamCard
    card = ParamCard(str(path))
    return {name:dict(scale=block.scale,parameters=sorted(
                (tuple(p.lhacode),float(p.value)) for p in block
                if not (name == 'decay' and tuple(p.lhacode) in ((6,),(24,)))))
            for name,block in card.items() if name != 'decaymass'}


def load_stable(execution_path):
    from madgraph.various.histograms import HwUList
    from models.check_param_card import ParamCard
    execution_path = Path(execution_path).resolve()
    archive = execution_path.parent
    execution = json.loads(execution_path.read_text())
    manifest = json.loads((archive/'manifest.json').read_text())
    if execution['status'] != 'finished' or execution['returncode'] != 0:
        raise ValueError('Stable reference did not finish successfully')
    if execution['variant'] not in ('LO','P') or execution['variant'] != manifest['variant']:
        raise ValueError('Unexpected stable-production prescription')
    for name,checksum in manifest['files'].items():
        if digest(archive/name) != checksum:
            raise ValueError('Stable input-card checksum mismatch: '+name)
    reference = Path(manifest['reference_manifest'])
    if digest(reference) != manifest['reference_manifest_sha256']:
        raise ValueError('Stable comparison reference changed')
    reference_manifest = json.loads(reference.read_text())
    if digest(STUDY/'inputs/benchmark.json') != manifest['benchmark_sha256']:
        raise ValueError('Stable benchmark changed')
    if not manifest['no_decay_branching_factor_applied']:
        raise ValueError('Stable reference already contains a branching factor')
    card = ParamCard(str(archive/'param_card.dat'))
    if any(card['decay'].get((pdg,)).value != 0. for pdg in (6,24)):
        raise ValueError('Stable on-shell connectors are not widthless')
    hwu = archive.parents[1]/'Events'/archive.name/'MADatNLO.HwU'
    if digest(hwu) != execution['output_sha256']:
        raise ValueError('Stable output checksum mismatch')
    totals = [h for h in HwUList(str(hwu),raw_labels=True) if h.title.strip() == 'total rate']
    if len(totals) != 1 or totals[0].bins[0].boundaries != (.5,1.5):
        raise ValueError('Unexpected stable inclusive histogram layout')
    row = totals[0].bins[0]
    columns = list(row.wgts)
    values = np.asarray([row.wgts[columns[i]] for i in canonical_columns(columns,execution['variant'])])
    if not np.isfinite(values).all() or row.wgts['dy'] < 0.:
        raise ValueError('Invalid stable reference weights/errors')
    for index in (1+POINTS.index((1.,1.,1.,1.)),82):
        if not np.isclose(values[0],values[index],rtol=3.e-6,atol=1.e-20):
            raise ValueError('Stable nominal duplicate mismatch')
    return dict(values=values,error=row.wgts['dy'],manifest=manifest,
                execution=execution,archive=archive,path=execution_path,
                reference=reference_manifest,
                parameters=production_parameters(archive/'param_card.dat'))


def input_bin(result):
    """All cut configurations share the preselection input, not a fiducial rate."""
    active = [int(result['offsets'][i]) for i,title in enumerate(result['titles'])
              if ' rates:' in title and result['values'][result['offsets'][i],0] != 0.]
    if len(active) != 5:
        raise ValueError('Require five input-rate copies for one active charge')
    if not all(np.array_equal(result['values'][i],result['values'][active[0]]) for i in active):
        raise ValueError('Cut configurations disagree before selection')
    return active[0]


def validate_bw_normalization(manifest, branching):
    """Require the independently checked current/convolution normalization."""
    path = STUDY/'inputs/bw_branching_normalization.json'
    record = json.loads(path.read_text())
    benchmark_path = STUDY/'inputs/benchmark.json'
    if (record['benchmark_sha256'] != digest(benchmark_path) or
            record['script_sha256'] != digest(STUDY/'scripts/bw_branching.py')):
        raise ValueError('BW current-normalization evidence changed')
    benchmark = json.loads(benchmark_path.read_text())
    mb = manifest['decay_bottom_mass']
    matching = [r for r in record['records'] if r['mb_GeV']==mb and r['W_width_factor']==1.]
    if len(matching)!=1 or not np.isclose(matching[0]['integrated_branching'],branching,rtol=2.e-10,atol=0.):
        raise ValueError('Missing or mismatched integrated BW leptonic density')
    widths = [r for r in benchmark['top_widths']
              if r['mb_GeV']==mb and r['w_treatment']=='bw' and r['scale_factor']==1.]
    if (len(widths)!=1 or manifest['top_width_w_treatment']!='bw' or
            manifest['top_width_lo']!=widths[0]['gamma_lo'] or
            manifest['top_width_nlo']!=widths[0]['gamma_nlo_pdf'] or
            manifest['w_width']!=matching[0]['W_width_GeV']):
        raise ValueError('BW numerator and denominator conventions do not match')
    return dict(path=str(path),sha256=digest(path),
                convention='Integrated fixed-coupling three-body density; no extra branching factor applied to generated events')


def compare(result,stable,branching):
    manifest = result['report']['manifest']
    generation = result['report']['execution']['export_generation']['args']
    variant = result['report']['variant']
    if (manifest['w_treatment'] not in ('onshell','top-bw') or
            manifest['production_scale'] != 'core-w-ht-half' or
            any(flavour not in ('e','mu') for flavour in generation['flavours']) or
            len(generation['flavours']) != 3):
        raise ValueError('This branching check requires an on-shell associated W and three direct e/mu decays')
    bw_evidence = (validate_bw_normalization(manifest,branching)
                   if manifest['w_treatment']=='top-bw' else None)
    if (STABLE_VARIANT[variant] != stable['execution']['variant'] or
            generation['charge'] != stable['execution']['charge']):
        raise ValueError('Wrong stable perturbative order or charge')
    physical = lambda settings: {k:v for k,v in settings.items() if k not in STATISTICAL}
    if physical(manifest['settings']) != physical(stable['manifest']['settings']):
        raise ValueError('Stable and decayed physical run settings differ')
    if manifest['settings']['iseed'] == stable['manifest']['settings']['iseed']:
        raise ValueError('Independent-error comparison requires distinct seeds')
    if (manifest['w_width'] != stable['reference']['w_width'] or
            branching != stable['execution']['branching_e']):
        raise ValueError('W normalization convention differs from the stable benchmark')
    archive = Path(result['report']['path']).parents[2]/'study_cards'/manifest['run_name']
    if production_parameters(archive/'param_card.dat') != stable['parameters']:
        raise ValueError('Stable and decayed production parameters differ')
    if not 0. < branching < 1.:
        raise ValueError('Invalid one-flavour W branching fraction')
    factor = branching**3
    index = input_bin(result)
    actual = result['values'][index]
    expected = factor*stable['values']
    delta = actual-expected
    error = np.hypot(result['errors'][index],factor*stable['error'])
    return clean(dict(variant=variant,charge=generation['charge'],flavours=generation['flavours'],
                      w_treatment=manifest['w_treatment'],bw_current_normalization=bw_evidence,
                      stable_variant=stable['execution']['variant'],branching_factor=factor,
                      decayed=describe(actual),expected=describe(expected),ratio=describe(ratio(actual,expected)),
                      weight_labels=WEIGHTS,weight_residuals_pb=delta,
                      central_residual_pb=delta[0],central_mc_error_pb=error,
                      central_pull=ratio(delta[0],error),
                      note='Nominal independent-run pull only. The 81 scale and 101 PDF residuals '
                           'do not have certified MC errors from these combined files.'))


def decay_scale_residuals(values):
    """At each production point subtract its paired central-decay weight."""
    indices = [1+POINTS.index((r,f,1.,1.)) for r,f,_,_ in POINTS]
    return values[...,1:82]-values[...,indices]


def check_joint_decay_scales(batch_path):
    batch_path = Path(batch_path).resolve()
    metadata = json.loads(batch_path.with_suffix('.json').read_text())
    if (digest(batch_path) != metadata['arrays_sha256'] or
            not metadata['status'].startswith('equal-count conditional')):
        raise ValueError('Require an unchanged successful joint-batch audit')
    record_path = Path(metadata['input'])
    if digest(record_path) != metadata['input_sha256']:
        raise ValueError('Joint-batch source record changed')
    record = json.loads(record_path.read_text())
    if digest(record['audit']) != record['audit_sha256']:
        raise ValueError('Joint-batch combined audit changed')
    result = load(record['audit'])
    manifest = result['report']['manifest']
    mode = manifest['w_treatment']
    if mode not in ('onshell','top-bw','all-bw'):
        raise ValueError('Unknown W treatment')
    if mode != 'onshell':
        benchmark = json.loads((STUDY/'inputs/benchmark.json').read_text())
        validate_bw_normalization(manifest,benchmark['W_width']['branching_e'])
    index = input_bin(result)
    with np.load(batch_path) as data:
        reduced = estimate(data['contributions'][:,index,:],data['strata'],decay_scale_residuals)
    return clean(dict(input=str(batch_path),input_sha256=digest(batch_path),
                      variant=result['report']['variant'],w_treatment=mode,points=POINTS,
                      residuals_pb=reduced['value'],conditional_joint_mc_errors_pb=reduced['mc_error'],
                      pulls=ratio(reduced['value'],reduced['mc_error']),
                      note='Same-event decay-weight differences; independent subprocess-stratum variances add. '
                           'The nine central-decay identities are exactly zero and have undefined pulls. '
                           'Conditional on training; full-run convergence remains required.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('audits',nargs='+',type=Path)
    parser.add_argument('--stable',nargs='+',required=True,type=Path)
    parser.add_argument('--joint-batches',nargs='*',default=[],type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Refuse to overwrite a branching-normalization report')
    benchmark = json.loads((STUDY/'inputs/benchmark.json').read_text())
    stable = {}
    for path in args.stable:
        row = load_stable(path)
        key = (row['execution']['charge'],row['execution']['variant'])
        if key in stable:
            raise ValueError('Multiple stable replicas require an explicit replica combiner')
        stable[key] = row
    rows = []
    for path in args.audits:
        result = load(path)
        charge = result['report']['execution']['export_generation']['args']['charge']
        variant = result['report']['variant']
        rows.append(compare(result,stable[(charge,STABLE_VARIANT[variant])],benchmark['W_width']['branching_e']))
    save(args.output,dict(created_utc=now(),status='inclusive normalization pilot; precision not certified',
                         inputs={str(p.resolve()):digest(p) for p in args.audits+args.stable},
                         comparisons=rows,joint_decay_scale_checks=[check_joint_decay_scales(p) for p in args.joint_batches]))
    print(args.output)
