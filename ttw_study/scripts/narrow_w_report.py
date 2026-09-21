#!/usr/bin/env python3
"""Generated narrow-W rates and shapes with a shared on-shell reference.

All scale/PDF coordinates are retained. Coarse diagnostic bins sum raw
worker vectors before normalization; they are not frozen publication bins.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from inclusive_report import STATISTICAL, production_parameters
from load_results import assert_same_layout, load
from narrow_w_inputs import COMMON_SCALE
from narrow_w_statistics import estimate
from pilot_report import RATES, clean
from rebinning import rebin
from run_narrow_w_pilots import DONE, cases, decay_card, width_for
from run_width_mass_pilots import verified_pairs

SPECTRA = (3, 4, 6, 11, 12, 20)


def validate_job(job, result, inputs, inputs_path):
    """Check physical cards, rather than inferring a limit from run names."""
    from models.check_param_card import ParamCard
    manifest = result['report']['manifest']
    case = job['case']
    width = width_for(case, inputs)
    archive = Path(job['process'])/'study_cards'/job['run']
    generation = result['report']['execution']['export_generation']['args']
    expected = dict(charge='plus', flavours=['e', 'e', 'mu'], corrected='both',
                    decay_bottom_mass=0., w_treatment=case['w_treatment'])
    if any(generation[k] != v for k, v in expected.items()):
        raise ValueError('Narrow-W samples require matched ordered production/decay flavours')
    if (manifest['narrow_W_test'] != case or manifest['variant'] != case['variant']
            or manifest['w_treatment'] != case['w_treatment']
            or manifest['top_width_reference_scale'] != COMMON_SCALE
            or manifest['production_scale'] != 'fixed'
            or manifest['physical_scale_coordinates'] != dict(
                common_fixed_scale_GeV=COMMON_SCALE, weights_relative_to_common_scale=True,
                top_width_mode='AUTO')
            or manifest['narrow_W_inputs_sha256'] != digest(inputs_path)
            or manifest['generated_event_weight_rescaling'] != 1.
            or not np.isclose(manifest['comparison_cross_section_rescaling'],
                              case['width_factor']**3, rtol=1.e-14)):
        raise ValueError('Incorrect narrow-W comparison or scale convention')
    for key, value in [('top_width_lo', width['gamma_lo']),
                       ('top_width_nlo', width['gamma_nlo']),
                       ('w_width', width['w_width_GeV'])]:
        if not np.isclose(manifest[key], value, rtol=1.e-13, atol=0.):
            raise ValueError('Unmatched narrow-W width: '+key)
    if (archive/'decay_card.dat').read_text() != decay_card(case, width):
        raise ValueError('Archived decay card differs from the matched common-scale convention')
    card = ParamCard(str(archive/'param_card.dat'))
    for pdg, value in [(6, width['gamma_nlo']), (24, width['w_width_GeV'])]:
        if not np.isclose(card['decay'].get((pdg,)).value, value, rtol=1.e-13, atol=0.):
            raise ValueError('Archived internal width differs from the matched input')
    if card['mass'].get((5,)).value != 0. or 'decaymass' in card:
        raise ValueError('Unexpected bottom mass in massless narrow-W test')
    settings = manifest['settings']
    for key in ('mur_ref_fixed', 'muf_ref_fixed', 'qes_ref_fixed'):
        if settings[key] != COMMON_SCALE:
            raise ValueError('Production reference scales are not common')
    for key in ('fixed_ren_scale', 'fixed_fac_scale', 'fixed_qes_scale'):
        if settings[key] is not True:
            raise ValueError('Narrow-W test requires fixed scales')
    return dict(parameters=production_parameters(archive/'param_card.dat'),
                settings={k: v for k, v in settings.items() if k not in STATISTICAL},
                analysis=manifest['hashes']['analysis_source'],
                scale_runtime=manifest['scale_runtime_hashes'])


def selected_vectors(job, result):
    """Keep only needed observables in memory after reading one full archive."""
    with np.load(job['batches']) as data:
        all_vectors = data['contributions']
        strata = data['strata'].copy()
    rates, spectra, edges = {}, {}, {}
    for h, title in enumerate(result['titles']):
        if 'W+' in title and ' rates:' in title:
            a = int(result['offsets'][h])
            rates[title] = all_vectors[:, a:a+9].copy()
    for local in SPECTRA:
        h = local-1
        a, b = map(int, result['offsets'][h:h+2])
        original = result['edges'][a:b]
        # Four adjacent original bins; retain the entire original range.
        requested = np.r_[original[::4, 0], original[-1, 1]]
        values = rebin(all_vectors[:, a:b], original, requested, axis=1)
        parent = 4 if local in (11, 12) else 3
        spectra[local] = np.concatenate([values, all_vectors[:, parent:parent+1]], axis=1)
        edges[local] = np.column_stack([requested[:-1], requested[1:]])
    return dict(rates=rates, spectra=spectra, edges=edges, strata=strata)


def run(queue_path, variant, output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a generated narrow-W comparison')
    queue = json.loads(queue_path.read_text())
    if queue['status'] != DONE or [j['case'] for j in queue['jobs']] != cases():
        raise ValueError('Require the complete, original 21-case narrow-W queue')
    jobs = [j for j in queue['jobs'] if j['variant'] == variant]
    if len(jobs) != 7:
        raise ValueError('Require one reference and six BW samples')
    inputs_path = Path(queue['width_inputs'])
    inputs = json.loads(inputs_path.read_text())
    if inputs['benchmark_sha256'] != digest(STUDY/'inputs/benchmark.json'):
        raise ValueError('Width input benchmark changed')
    selected, source_rows, streams, counts = [], [], set(), []
    signature, layout = None, None
    for job in jobs:
        result = load(job['audit'])
        current = validate_job(job, result, inputs, inputs_path)
        if signature is not None and current != signature:
            raise ValueError('Narrow-W samples differ beyond the declared W treatment/width')
        if layout is not None:
            assert_same_layout([layout, result])
        signature = current
        layout = {k: result[k] for k in ('titles', 'offsets', 'edges')}
        pairs = verified_pairs(job)
        if pairs & streams:
            raise ValueError('Shared random streams between narrow-W samples')
        streams.update(pairs)
        counts.append(len(pairs))
        selected.append(selected_vectors(job, result))
        source_rows.append(dict(case=job['case'], audit=job['audit'], audit_sha256=digest(job['audit']),
                                batches=job['batches'], batches_sha256=digest(job['batches'])))
        print('Loaded and matched', variant, job['case']['w_treatment'],
              job['case']['width_factor'], flush=True)
    factors = [j['case']['width_factor'] for j in jobs]
    strata = [row['strata'] for row in selected]
    report = dict(created_utc=now(), variant=variant, queue_sha256=digest(queue_path),
                  width_inputs_sha256=digest(inputs_path), samples=source_rows,
                  stage_stream_counts=counts, distinct_initialization_pairs=len(streams),
                  cross_sample_rng_overlap=0, rates={}, spectra={},
                  status='generated narrow-W comparison evaluated; precision and convergence must be inspected',
                  convention='epsilon cubed applied once to cross sections; shared on-shell reference. '
                  'Independent run/beam-stratum deletion retains bin/parent and scale/PDF correlations. '
                  'Arrays retain nominal, 81 scales and 101 PDFs; covariance factors stored for nominal.',
                  limitations='One W+ e,e,mu assignment, fixed common scales, conditional learned-grid errors. '
                  'Coarse bins are diagnostic; agreement at finite epsilon is not a proof of the exact '
                  'limit, full-flavour precision, independent retraining, or all sparse-tail convergence.')
    arrays = {}
    for index, title in enumerate(selected[0]['rates']):
        result = estimate([r['rates'][title] for r in selected], strata, factors)
        prefix = 'rates_%d' % index
        arrays[prefix] = result['value']
        arrays[prefix+'_mc_errors'] = result['mc_error']
        arrays[prefix+'_covariance_factor_nominal'] = result['covariance_factor'][..., 0]
        report['rates'][title] = {label: dict(
            rescaled_pb=result['value'][0, :, i, 0], rescaled_mc_error_pb=result['mc_error'][0, :, i, 0],
            difference_pb=result['value'][1, :, i, 0], difference_mc_error_pb=result['mc_error'][1, :, i, 0],
            ratio=result['value'][2, :, i, 0], ratio_mc_error=result['mc_error'][2, :, i, 0])
            for i, label in enumerate(RATES)}
        del result
    for local in SPECTRA:
        vectors = [r['spectra'][local] for r in selected]
        for normalized in (False, True):
            result = estimate(vectors if normalized else [v[:, :-1] for v in vectors],
                              strata, factors, normalized=normalized)
            prefix = ('normalized_' if normalized else 'absolute_')+str(local)
            arrays[prefix] = result['value']
            arrays[prefix+'_mc_errors'] = result['mc_error']
            arrays[prefix+'_covariance_factor_nominal'] = result['covariance_factor'][..., 0]
            if normalized:
                pulls = np.divide(result['value'][1, :, :, 0], result['mc_error'][1, :, :, 0],
                                  out=np.full_like(result['value'][1, :, :, 0], np.nan),
                                  where=result['mc_error'][1, :, :, 0] > 0.)
                report['spectra'][str(local)] = dict(title=layout['titles'][local-1],
                    parent_region='2b' if local in (11, 12) else '1b', nominal_normalized_pulls=pulls,
                    maximum_absolute_nominal_pull_by_sample=[
                        float(np.nanmax(np.abs(p))) if np.isfinite(p).any() else None for p in pulls],
                    edges=selected[0]['edges'][local])
            del result
        arrays['edges_'+str(local)] = selected[0]['edges'][local]
    np.savez_compressed(output.with_suffix('.npz'), **arrays)
    report['arrays_sha256'] = digest(output.with_suffix('.npz'))
    report['script_sha256'] = digest(Path(__file__))
    save(output, clean(report))
    print(output, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('queue', type=Path)
    parser.add_argument('--variant', required=True, choices=('LO', 'S', 'Pi'))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    run(args.queue, args.variant, args.output)
