"""Canonical full weight grids for six-variant algebra and replica statistics."""
import itertools
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from campaign import digest

FACTORS = (.5,1.,2.)
POINTS = list(itertools.product(FACTORS,repeat=4))
WEIGHTS = ['central'] + [','.join('%g' % value for value in point) for point in POINTS]
WEIGHTS += ['PDF=%d' % member for member in range(331700,331801)]


def parameter_signature(path):
    """Exact numeric SLHA contents, ignoring comments and record ordering.

    Raw hashes still certify the archived cards. MG5 can add whitespace to
    QNUMBERS comments when round-tripping the same physical parameters.
    """
    from models.check_param_card import ParamCard

    def block_signature(block):
        parameters = []
        for parameter in block:
            value = float(parameter.value)
            if not math.isfinite(value):
                raise ValueError('Unresolved or nonfinite physical parameter')
            parameters.append((tuple(parameter.lhacode),value))
        return dict(scale=block.scale,parameters=sorted(parameters),
                    decay_table={str(pdg):block_signature(table) for pdg,table in block.decay_table.items()})
    card = ParamCard(str(path))
    numeric = {name:block_signature(block) for name,block in card.items()}
    return hashlib.sha256(json.dumps(numeric,sort_keys=True,allow_nan=False).encode()).hexdigest()


def canonical_columns(columns, variant):
    """Broadcast only verified LO/P production weights on inactive decay axes."""
    production_only = variant in ('LO','P')
    scales, pdfs = {}, {}
    for i,label in enumerate(columns):
        if label.startswith('PDF='):
            match = re.match(r'PDF=\s*(\d+)',label)
            if not match:
                raise ValueError('Unrecognized PDF weight label: '+label)
            member = int(match[1])
            if member in pdfs:
                raise ValueError('Duplicate PDF member')
            pdfs[member] = i
        match = re.search(r'muR=\s*([\d.]+)\s+muF=\s*([\d.]+)',label)
        if not match:
            continue
        point = tuple(map(float,match.groups()))
        decay = [(int(pdg),float(factor)) for pdg,factor in
                 re.findall(r'(?:^|\s)d(-?\d+)=\s*([\d.]+)',label)]
        if production_only:
            if decay:
                raise ValueError('LO/P unexpectedly contains decay-scale axes')
        else:
            if len(decay) != 2 or {pdg for pdg,_ in decay} != {6,-6}:
                raise ValueError('NLO decays require both signed scale axes')
            decay = dict(decay)
            point += (decay[6],decay[-6])
        if point in scales:
            raise ValueError('Duplicate scale point')
        scales[point] = i
    expected = set(itertools.product(FACTORS,repeat=2 if production_only else 4))
    if set(scales) != expected:
        raise ValueError('Incomplete scale grid')
    if set(pdfs) != set(range(331700,331801)):
        raise ValueError('Incomplete benchmark PDF grid')
    indices = [columns.index('central value')]
    indices += [scales[point[:2] if production_only else point] for point in POINTS]
    indices += [pdfs[member] for member in range(331700,331801)]
    return indices


def load(audit_path):
    """Read a frozen successful audit and canonicalize its compressed arrays."""
    audit_path = Path(audit_path).resolve()
    report = json.loads(audit_path.read_text())
    if not report['corrected_accumulator'] or report['status'] == 'diagnostic only':
        raise ValueError('Diagnostic/old-accumulator output cannot enter physics comparisons')
    hwu_path = Path(report['path'])
    if digest(hwu_path) != report['output_sha256']:
        raise ValueError('Original histogram checksum mismatch')
    archive = hwu_path.parents[2]/'study_cards'/hwu_path.parent.name
    for name, expected in report['manifest']['hashes'].items():
        if name.endswith('.dat') and digest(archive/name) != expected:
            raise ValueError('Archived input card checksum mismatch: '+name)
    hook = report['manifest'].get('custom_user_hook')
    if hook and digest(archive/Path(hook['source']).name) != hook['sha256']:
        raise ValueError('Archived custom user-hook checksum mismatch')
    array_path = audit_path.with_suffix('.npz')
    if digest(array_path) != report['arrays_sha256']:
        raise ValueError('Weight archive checksum mismatch')
    with np.load(array_path) as data:
        columns = list(data['columns'])
        indices = canonical_columns(columns,report['variant'])
        values = data['values'][:,indices]
        errors = data['values'][:,columns.index('dy')]
        result = dict(report=report, values=values, errors=errors,
                      parameter_signature=parameter_signature(archive/'param_card.dat'),
                      edges=data['edges'].copy(),offsets=data['offsets'].copy(),
                      titles=list(data['titles']),weights=WEIGHTS,
                      canonicalization=('LO/P broadcast on inactive signed decay axes'
                                        if report['variant'] in ('LO','P') else 'signed-axis ordering'))
    assert np.isfinite(values).all()
    # The PDF nominal is an independent duplicate in the output format.
    for index in (1+POINTS.index((1.,1.,1.,1.)),82):
        assert np.allclose(values[:,0],values[:,index],rtol=3.e-6,atol=1.e-20)
    return result


def assert_same_layout(results):
    reference = results[0]
    for result in results[1:]:
        if result['titles'] != reference['titles'] or not np.array_equal(
                result['edges'],reference['edges']) or not np.array_equal(
                result['offsets'],reference['offsets']):
            raise ValueError('Histogram layouts differ')


def assert_matched_physics(results):
    """No algebraic identity may mix different flavours, widths, PDFs or cuts."""
    statistical = {'iseed','req_acc_fo','npoints_fo','niters_fo',
                   'npoints_fo_grid','niters_fo_grid','fo_job_target_time'}
    signatures = []
    for result in results:
        report = result['report']
        manifest = report['manifest']
        generation = report['execution']['export_generation']['args']
        signatures.append(dict(
            process={key:generation[key] for key in
                     ('charge','flavours','w_treatment','decay_bottom_mass','corrected')},
            settings={key:value for key,value in manifest['settings'].items() if key not in statistical},
            widths={key:manifest[key] for key in
                    ('top_width_lo','top_width_nlo','w_width','top_width_reference_scale',
                     'top_width_w_treatment','top_width_bottom_mass')},
            analysis=manifest['hashes']['analysis_source'],
            parameters=result.get('parameter_signature',manifest['hashes']['param_card.dat']),
            production_scale=manifest['production_scale'],
            physical_scale_coordinates=manifest.get('physical_scale_coordinates'),
            custom_user_hook_sha256=(manifest.get('custom_user_hook') or {}).get('sha256'),
            decay_scale_grouping=manifest['decay_scale_grouping']))
    if any(signature != signatures[0] for signature in signatures[1:]):
        raise ValueError('Physical inputs or selected charge/flavours differ')


def strict_identity(results):
    """S-P-D+LO; central MC errors assume independently seeded integrations."""
    if set(results) != {'LO','P','D','S'}:
        raise ValueError('Four variants required for strict-NLO identity')
    ordered = [results[name] for name in ('LO','P','D','S')]
    if any(results[name]['report']['variant'] != name for name in results):
        raise ValueError('Variant label disagrees with archived output')
    assert_same_layout(ordered)
    assert_matched_physics(ordered)
    seeds = [row['report']['manifest']['settings']['iseed'] for row in ordered]
    if len(set(seeds)) != 4:
        raise ValueError('Independence not established: repeated seeds')
    delta = results['S']['values']-results['P']['values']-results['D']['values']+results['LO']['values']
    error = np.sqrt(sum(row['errors']**2 for row in ordered))
    pulls = np.divide(delta[:,0],error,out=np.full(error.shape,np.nan),where=error>0.)
    return dict(difference=delta,central_mc_error=error,central_pull=pulls,
                note='Per-bin pulls are correlated; this is not a global chi-squared test. No varied-weight MC errors inferred.')
