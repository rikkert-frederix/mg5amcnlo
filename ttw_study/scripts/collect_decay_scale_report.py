#!/usr/bin/env python3
"""Collect the four current ttW runs and compare width/scale prescriptions.

This is a read-only reduction of completed integrations. Native histogram
bins are retained: no cross-bin or cross-run Monte Carlo covariance is
invented. Archived RunMaterial cards, rather than current cards, define runs.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import re
import tarfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT/'ttw_study/reports/decay_scales_20261019'
RUNS = {
    'S_auto': ('2', 'additive_v2'),
    'S_fixed': ('2', 'additive_v2_FixedWidth'),
    'Pi_auto': ('3', 'multiplicative_v2'),
    'Pi_fixed': ('3', 'multiplicative_v2_FixedWidth'),
}
CAMPAIGN_DIRECTORIES = {
    'S_auto': 'ttw_additive',
    'S_fixed': 'ttw_additive_fixedwidth',
    'Pi_auto': 'ttw_multiplicative',
    'Pi_fixed': 'ttw_multiplicative_fixedwidth',
}
spec = importlib.util.spec_from_file_location('ttw_scale_points',
    ROOT/'Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py')
scale_tools = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scale_tools)
GROUPS = dict(scale_tools.GROUPS)
GROUPS['correlated7'] = {(r, f, r, r) for r, f in scale_tools.SEVEN}
RATE_TITLES = ['inclusive rate', 'fiducial cutflow',
               'fiducial 0b rate', 'fiducial 1b rate', 'fiducial 2b rate']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_value(value):
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    return value


def save_json(path, value):
    path.write_text(json.dumps(json_value(value), indent=2, allow_nan=False)+'\n')


def read_hwu(path):
    titles, arrays, rows = [], [], []
    with path.open() as stream:
        columns = [x.strip() for x in next(stream).split('&')[1:]]
        for line in stream:
            if line.startswith('<histogram>'):
                title = re.search(r'"(.*?)\s*\|', line)[1].strip()
                rows = []
            elif line.startswith('<\\histogram>'):
                titles.append(title)
                arrays.append(np.asarray(rows))
            elif line.strip() and not line.startswith('#'):
                values = np.fromstring(line, sep=' ')
                if len(values) != len(columns):
                    raise ValueError(f'Column count in {path}: {len(values)}')
                rows.append(values)
    if len(titles) != len(set(titles)):
        raise ValueError('Duplicate histogram titles')
    offsets = np.r_[0, np.cumsum([len(a) for a in arrays])]
    values = np.concatenate(arrays)
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite histogram entry')
    return dict(titles=titles, columns=columns, offsets=offsets,
                edges=values[:, :2], central=values[:, 2],
                error=values[:, 3], scales=values[:, 4:])


def card_values(text):
    out = {}
    for line in text.splitlines():
        line = line.split('!')[0].split('#')[0].strip()
        if '=' in line:
            value, key = line.split('=', 1)
            out[key.strip().lower()] = value.strip()
    return out


def read_run(key, output, directory=None):
    replica = directory is not None
    if directory is None:
        proc, run = RUNS[key]
        directory = ROOT/f'PROC_fNLO_loop_sm-no_b_mass_{proc}/Events/{run}'
    directory = Path(directory)
    run = directory.name
    path = directory/'MADatNLO.HwU'
    data = read_hwu(path)
    if len(data['titles']) != 581:
        raise ValueError('Require the complete new fiducial analysis')
    points = [scale_tools.scale_point(label) for label in data['columns'][4:]]
    if len(points) != 81 or set(points) != scale_tools.POINTS:
        raise ValueError('Missing, duplicate or mislabelled scale coordinates')
    index = points.index((1., 1., 1., 1.))
    mismatch = np.max(np.abs(data['central']-data['scales'][:, index]))
    if not np.allclose(data['central'], data['scales'][:, index], rtol=3e-6, atol=1e-18):
        raise ValueError('Nominal and central scale weights disagree')
    input_dir = output/'inputs'/key
    if replica:
        input_dir = input_dir/run
    input_dir.mkdir(parents=True, exist_ok=True)
    archived = {}
    with tarfile.open(directory/'RunMaterial.tar.gz') as archive:
        for name in ('run_card.dat', 'param_card.dat', 'decay_card.dat', 'FO_analyse_card.dat'):
            raw = archive.extractfile('RunMaterial/'+name).read()
            (input_dir/name).write_bytes(raw)
            archived[name] = raw.decode()
    rc, dc = card_values(archived['run_card.dat']), card_values(archived['decay_card.dat'])
    expected_mode = 'AUTO' if key.endswith('auto') else 'EXPLICIT'
    if dc['decay_width_scale_mode(6)'] != expected_mode:
        raise ValueError('Run label disagrees with archived width mode')
    if dc['nlo_decay_combination'] != ('ADDITIVE' if key.startswith('S') else 'MULTIPLICATIVE'):
        raise ValueError('Run label disagrees with archived NLO approximation')
    if key.endswith('fixed'):
        for pdg, orders in ((6, ('lo', 'nlo')), (24, ('lo',))):
            for order, factor in itertools.product(orders, ('.5', '2')):
                candidates = [v for k, v in dc.items() if re.fullmatch(
                    order+r'_decay_width_variation\('+str(pdg)+r',\s*[^)]+\)', k)
                    and float(k.split(',')[1][:-1]) == float(factor)]
                if len(candidates) != 1 or float(candidates[0]) != float(dc[f'{order}_decay_width({pdg})']):
                    raise ValueError('Explicit widths are not fixed at the central values')
    rng = []
    for log in sorted(directory.glob('alllogs_*.html')):
        text = log.read_text()
        rng.extend(tuple(map(int, pair)) for pair in re.findall(
            r'Ranmar initialization seeds\s+(\d+)\s+(\d+)', text))
    metadata = dict(key=key, run=run, path=str(path.relative_to(ROOT)),
        hwu_sha256=digest(path), archive_sha256=digest(directory/'RunMaterial.tar.gz'),
        banner_sha256=digest(next(directory.glob('*banner.txt'))),
        archived_cards={name: digest(input_dir/name) for name in archived},
        run_settings={k: rc.get(k) for k in ('ebeam1', 'ebeam2', 'lhaid',
            'dynamical_scale_choice', 'req_acc_fo', 'iseed', 'reweight_scale', 'reweight_pdf')},
        decay_settings=dc, histogram_count=len(data['titles']),
        scale_points=points, nominal_scale_max_difference_pb=float(mismatch),
        rng_pairs=sorted(set(rng)))
    data['metadata'] = metadata
    return data


def average_replicas(central, error, scales):
    """Equal-weight mean; propagate native errors without variance weighting."""
    count = len(central)
    return dict(central=np.mean(central, axis=0),
                error=np.sqrt(np.sum(np.square(error), axis=0))/count,
                empirical_error=np.std(central, axis=0, ddof=1)/np.sqrt(count),
                scales=np.mean(scales, axis=0))


def paired_ratio_error(numerator, denominator):
    """Delete-one paired-replica jackknife SE of the ratio of ensemble means."""
    numerator, denominator = np.asarray(numerator), np.asarray(denominator)
    count = len(numerator)
    if count < 2 or numerator.shape != denominator.shape:
        raise ValueError('Need matching independent replica blocks')
    den = denominator.sum(axis=0)-denominator
    leave_one_out = np.divide(numerator.sum(axis=0)-numerator, den,
        out=np.full_like(numerator, np.nan, dtype=float), where=den > 0)
    return np.sqrt((count-1)/count*np.sum(
        (leave_one_out-leave_one_out.mean(axis=0))**2, axis=0))


def read_campaign(key, output, campaign_root):
    process = campaign_root/CAMPAIGN_DIRECTORIES[key]
    directories = sorted(p for p in (process/'Events').glob('run_*')
                         if (p/'MADatNLO.HwU').is_file())
    if len(directories) < 2:
        raise ValueError('Need at least two complete campaign replicas: '+str(process))
    replicas = []
    for directory in directories:
        print('Collecting', key, directory.name, flush=True)
        replicas.append(read_run(key, output, directory))
    ref = replicas[0]
    for data in replicas[1:]:
        for field in ('titles', 'columns'):
            if data[field] != ref[field]:
                raise ValueError('Replica histogram layouts differ')
        if not np.array_equal(data['edges'], ref['edges']):
            raise ValueError('Replica bin boundaries differ')
        for field in ('archived_cards', 'run_settings', 'decay_settings'):
            if data['metadata'][field] != ref['metadata'][field]:
                raise ValueError('Replica input settings differ: '+field)
    central = np.stack([d['central'] for d in replicas])
    errors = np.stack([d['error'] for d in replicas])
    averaged = average_replicas(central, errors, np.stack([d['scales'] for d in replicas]))
    data = {k: ref[k] for k in ('titles', 'columns', 'offsets', 'edges')}
    data.update(averaged, replica_central=central, replica_error=errors)
    metadata = dict(ref['metadata'])
    metadata.update(run='equal-weight mean of '+str(len(replicas))+' replicas',
        path=str(process.relative_to(ROOT)), replica_count=len(replicas),
        replica_ids=[p.name for p in directories],
        replicas=[d['metadata'] for d in replicas],
        rng_pairs=sorted(set.union(*(set(map(tuple,d['metadata']['rng_pairs'])) for d in replicas))),
        aggregation='equal-weight mean of completed integrations; native MC errors propagated in quadrature',
        source_snapshots={})
    del metadata['hwu_sha256'], metadata['archive_sha256'], metadata['banner_sha256']
    for relative in ('FixedOrderAnalysis/analysis_HwU_pp_ttxw_trilepton.f90',
                     'SubProcesses/cuts.f90'):
        source = process/relative
        target = output/'inputs'/key/'source'/source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source.read_bytes())
        metadata['source_snapshots'][relative] = digest(source)
    data['metadata'] = metadata
    return data


def hist_slice(data, title):
    i = data['titles'].index(title)
    return slice(int(data['offsets'][i]), int(data['offsets'][i+1]))


def band(values, points, group):
    selected = [points.index(point) for point in sorted(GROUPS[group])]
    return np.min(values[..., selected], axis=-1), np.max(values[..., selected], axis=-1)


def describe_vector(central, error, values, points):
    result = dict(central_pb=central, mc_error_pb=error,
                  relative_mc_error=error/abs(central) if central else None,
                  scale_points={','.join(map(str, point)): float(values[i]) for i, point in enumerate(points)})
    for name in GROUPS:
        low, high = band(values, points, name)
        result[name] = dict(low_pb=low, high_pb=high, npoints=len(GROUPS[name]),
            down_percent=100*(1-low/central) if central else None,
            up_percent=100*(high/central-1) if central else None)
    return result


def collect(output, campaign_root=None):
    output.mkdir(parents=True, exist_ok=True)
    runs = {key: (read_campaign(key, output, campaign_root) if campaign_root else
                  read_run(key, output)) for key in RUNS}
    reference = runs['S_auto']
    points = reference['metadata']['scale_points']
    for key, data in runs.items():
        if (data['titles'] != reference['titles'] or
            data['columns'] != reference['columns'] or
            not np.array_equal(data['edges'], reference['edges'])):
            raise ValueError('Runs have different histogram or scale layouts')
        if data['metadata']['archived_cards']['param_card.dat'] != reference['metadata']['archived_cards']['param_card.dat']:
            raise ValueError('Archived masses, couplings or reference widths differ')
        for field in ('ebeam1', 'ebeam2', 'lhaid', 'dynamical_scale_choice'):
            if data['metadata']['run_settings'][field] != reference['metadata']['run_settings'][field]:
                raise ValueError('Unmatched physics setting: '+field)
        for field, value in reference['metadata']['decay_settings'].items():
            if field == 'nlo_decay_combination' or field.startswith('decay_width_scale_mode'):
                continue
            if data['metadata']['decay_settings'].get(field) != value:
                raise ValueError('Unmatched decay or production setting: '+field)
    report = dict(schema='ttw_current_width_scale_report_v1', scope=
        '13 TeV W+ e+nu_e e-nubar_e mu+nu_mu; one direct-lepton assignment, not a full-flavour or charge sum',
        source_runs={key: d['metadata'] for key, d in runs.items()},
        scale_groups={key: sorted(value) for key, value in GROUPS.items()},
        rates={}, rate_comparisons={}, rng_overlap={}, precision_inventory=[])
    if campaign_root:
        report['schema'] = 'ttw_replicated_width_scale_report_v2'
        ids = reference['metadata']['replica_ids']
        if any(d['metadata']['replica_ids'] != ids for d in runs.values()):
            raise ValueError('Campaign replica blocks do not match')
        blocks = [set.union(*(set(map(tuple, d['metadata']['replicas'][i]['rng_pairs']))
                    for d in runs.values())) for i in range(len(ids))]
        overlaps = [dict(a=ids[i], b=ids[j], shared_pairs=len(a & b))
                    for (i,a),(j,b) in itertools.combinations(enumerate(blocks),2) if a & b]
        if overlaps:
            raise ValueError('Random streams overlap between replica blocks: '+str(overlaps))
        report['replica_blocks'] = dict(ids=ids, count=len(ids),
            cross_block_rng_overlaps=overlaps,
            convention='Matching run indices form paired blocks across S/Pi and width modes; blocks use distinct logged RNG seed pairs.')
    for key, data in runs.items():
        rates = {}
        for title in RATE_TITLES:
            s = hist_slice(data, title)
            i = s.stop-1 if title == 'fiducial cutflow' else s.start
            name = 'fiducial baseline' if title == 'fiducial cutflow' else title
            rates[name] = describe_vector(data['central'][i], data['error'][i], data['scales'][i], points)
            if campaign_root:
                rates[name].update(empirical_mc_error_pb=float(data['empirical_error'][i]),
                    empirical_over_native=float(data['empirical_error'][i]/data['error'][i]),
                    replica_central_pb=data['replica_central'][:,i],
                    replica_native_error_pb=data['replica_error'][:,i])
        report['rates'][key] = rates
        tags = sum(data['scales'][hist_slice(data, f'fiducial {n}b rate').start] for n in range(3))
        base = data['scales'][hist_slice(data, 'fiducial cutflow').stop-1]
        residual = np.max(np.abs(tags-base))
        if not np.allclose(tags, base, rtol=3e-6, atol=1e-17):
            raise ValueError('Exclusive category rates do not close on the common baseline')
        data['metadata']['category_closure_max_pb'] = float(residual)
        for title in data['titles']:
            s = hist_slice(data, title)
            c, e = data['central'][s], data['error'][s]
            threshold = .01*np.max(c)
            use = (c > max(threshold, 0.)) & (e > 0)
            report['precision_inventory'].append(dict(run=key, title=title,
                populated_bins=int(use.sum()), negative_bins=int((c < 0).sum()),
                median_relative_mc_error=float(np.median(e[use]/c[use])) if use.any() else None,
                p90_relative_mc_error=float(np.quantile(e[use]/c[use], .9)) if use.any() else None))
            if campaign_root:
                report['precision_inventory'][-1].update(
                    median_empirical_relative_mc_error=float(np.median(
                        data['empirical_error'][s][use]/c[use])) if use.any() else None)
    for width in ('auto', 'fixed'):
        numerator, denominator = runs['Pi_'+width], runs['S_'+width]
        rows = {}
        for title in RATE_TITLES:
            s = hist_slice(reference, title)
            i = s.stop-1 if title == 'fiducial cutflow' else s.start
            ratios = numerator['scales'][i]/denominator['scales'][i]
            rows[title] = dict(central_ratio=float(
                    numerator['central'][i]/denominator['central'][i]),
                scale_envelopes={name: band(ratios, points, name) for name in GROUPS})
            if campaign_root:
                rows[title]['paired_mc_error'] = float(paired_ratio_error(
                    numerator['replica_central'][:, i], denominator['replica_central'][:, i]))
        report['rate_comparisons'][width] = rows
    for a, b in itertools.combinations(RUNS, 2):
        aa, bb = (set(map(tuple, runs[x]['metadata']['rng_pairs'])) for x in (a, b))
        report['rng_overlap'][a+' / '+b] = len(aa & bb)
    arrays = {key+'__'+field: data[field] for key, data in runs.items()
              for field in ('central', 'error', 'scales')}
    if campaign_root:
        arrays.update({key+'__'+field: data[field] for key, data in runs.items()
            for field in ('empirical_error','replica_central','replica_error')})
        report['paired_width_comparisons'] = {}
        for variant in ('S','Pi'):
            a, b = runs[variant+'_fixed'], runs[variant+'_auto']
            report['paired_width_comparisons'][variant] = {}
            for title in RATE_TITLES:
                s = hist_slice(reference, title)
                i = s.stop-1 if title == 'fiducial cutflow' else s.start
                differences = a['replica_central'][:,i]-b['replica_central'][:,i]
                report['paired_width_comparisons'][variant][title] = dict(
                    difference_pb=np.mean(differences),
                    paired_mc_error_pb=np.std(differences, ddof=1)/np.sqrt(len(differences)),
                    max_replica_abs_difference_pb=np.max(np.abs(differences)))
    arrays.update(edges=reference['edges'], offsets=reference['offsets'],
                  titles=np.asarray(reference['titles']), points=np.asarray(points))
    np.savez_compressed(output/'collected.npz', **arrays)
    report['collected_sha256'] = digest(output/'collected.npz')
    report['collection_script_sha256'] = digest(__file__)
    save_json(output/'results.json', report)
    return runs, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--campaign-root', type=Path,
        help='Collect every Events/run_* integration from the four ttw_* directories under this root')
    args = parser.parse_args()
    _, report = collect(args.output.resolve(),
                        args.campaign_root.resolve() if args.campaign_root else None)
    for key, rows in report['rates'].items():
        print(key)
        for name, row in rows.items():
            print(' ', name, '%.8f +/- %.8f fb' % (1000*row['central_pb'], 1000*row['mc_error_pb']),
                  'C7 / I63:', [(name, round(row[name]['down_percent'], 3),
                    round(row[name]['up_percent'], 3)) for name in ('correlated7', 'combined63')])
    print('RNG overlaps:', report['rng_overlap'])


if __name__ == '__main__':
    main()
