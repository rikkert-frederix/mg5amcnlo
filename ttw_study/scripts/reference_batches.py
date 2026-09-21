#!/usr/bin/env python3
"""Small joint-batch reader for common-scale, no-PDF reference calculations."""
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

import numpy as np

from campaign import digest, now, save
from madgraph.various.histograms import HwUList
from read_splits import validate_jobs
from rng_history import validate_history, validate_refinement_history

FACTORS = (1., .5, 2.)
POINTS = [(r, f) for r in FACTORS for f in FACTORS]


def columns(labels):
    mapping = {}
    for i, label in enumerate(labels):
        if label in ('central value', 'dy') or re.fullmatch(r'delta_mu_(?:cen|min|max) -?\d+ @aux', label):
            continue
        match = re.search(r'muR=\s*([0-9.eE+-]+)\s+muF=\s*([0-9.eE+-]+)', label)
        if not match:
            raise ValueError('Unexpected reference weight: '+label)
        point = tuple(map(float, match.groups()))
        if point in mapping or point not in POINTS:
            raise ValueError('Duplicate or unsupported reference scale point')
        for pdg, factor in re.findall(r'\bd(-?\d+)=\s*([0-9.eE+-]+)', label):
            if abs(int(pdg)) != 6 or float(factor) != point[0]:
                raise ValueError('Reference decay numerator does not follow production muR')
        mapping[point] = i
    if set(mapping) != set(POINTS) or labels.count('central value') != 1 or labels.count('dy') != 1:
        raise ValueError('Reference requires nominal, error and nine common-scale weights')
    return [labels.index('central value')]+[mapping[p] for p in POINTS]


def vectors(raw):
    histograms = HwUList(io.StringIO(raw.decode()), raw_labels=True)
    labels = list(histograms[0].bins[0].wgts)
    selected = columns(labels)
    rows = [b for h in histograms for b in h.bins]
    data = np.asarray([[b.wgts[label] for label in labels] for b in rows])
    if not np.isfinite(data).all():
        raise ValueError('Nonfinite reference histogram')
    return dict(values=data[:,selected], errors=data[:,labels.index('dy')],
                titles=[h.title for h in histograms],
                offsets=np.r_[0,np.cumsum([len(h.bins) for h in histograms])],
                edges=np.asarray([b.boundaries for b in rows]))


def read(record_path, output):
    if output.exists() or output.with_suffix('.json').exists():
        raise ValueError('Refuse to overwrite reference batches')
    record = json.loads(record_path.read_text())
    archive_path = record_path.with_suffix('.gz')
    if digest(archive_path) != record['archive_sha256']:
        raise ValueError('Reference worker archive changed')
    if digest(record['audit']) != record['audit_sha256']:
        raise ValueError('Reference final-output audit changed')
    audit = json.loads(Path(record['audit']).read_text())
    run_log=Path(audit['execution']['log']).read_text()
    if not record.get('history_files'):
        validate_refinement_history(run_log)
    if digest(audit['path']) != audit['output_sha256']:
        raise ValueError('Final reference histogram changed')
    final = vectors(Path(audit['path']).read_bytes())
    workers = record['workers']
    groups = validate_jobs(workers)
    index = {w['directory']: i for i, w in enumerate(workers)}
    expected = {w['directory']+'/'+name: checksum for w in workers
                for name, checksum in w['files'].items()}
    expected.update(record.get('history_files',{}))
    history={}
    arrays, errors, streams, initial = {}, {}, {}, {}
    with tarfile.open(archive_path, 'r|gz') as archive:
        for member in archive:
            if not member.isfile():
                raise ValueError('Unexpected worker archive entry')
            raw = archive.extractfile(member).read()
            if hashlib.sha256(raw).hexdigest() != expected.pop(member.name, None):
                raise ValueError('Changed or unexpected reference worker file')
            if member.name in record.get('history_files',{}):
                history[member.name]=raw
                continue
            directory, name = member.name.rsplit('/',1)
            i = index[directory]
            if name == 'log.txt':
                log = raw.decode()
                seeds = re.findall(r'Ranmar initialization seeds\s+(\d+)\s+(\d+)', log)
                if len(seeds) != 1 or re.findall(r'with seed\s+(\d+)', log) != [str(record['seed'])]:
                    raise ValueError('Unexpected reference worker random stream')
                streams[i] = tuple(map(int, seeds[0]))
                initial[i] = re.findall(r'^channel.*$', log.split('------- iteration')[0], re.M)
                if not initial[i]:
                    raise ValueError('Missing initial grid/channel diagnostics')
            elif name == 'MADatNLO.HwU':
                row = vectors(raw)
                if row['titles'] != final['titles'] or not np.array_equal(row['edges'], final['edges']):
                    raise ValueError('Worker reference layout changed')
                arrays[i], errors[i] = row['values'], row['errors']
    if expected or len(arrays) != len(workers) or len(streams) != len(workers):
        raise ValueError('Incomplete reference worker archive')
    if len(set(streams.values())) != len(workers):
        raise ValueError('Reference random-stream collision')
    stage_audit=validate_history(run_log,record['seed'],history,
        {w['directory']:w['files']['log.txt'] for w in workers}) if history else None
    strata = np.empty(len(workers), dtype=int)
    for label, indices in enumerate(groups.values()):
        if any(initial[i] != initial[indices[0]] for i in indices):
            raise ValueError('Reference proposal differs within a stratum')
        strata[indices] = label
    values = np.asarray([arrays[i] for i in range(len(workers))])
    errors = np.asarray([errors[i] for i in range(len(workers))])
    tolerance = 3.e-6*(np.abs(values).sum(axis=0)+np.abs(final['values']))+1.e-20
    if not (np.abs(values.sum(axis=0)-final['values']) <= tolerance).all():
        raise ValueError('Reference worker sum does not reproduce final weights')
    if not np.allclose(np.sqrt((errors**2).sum(axis=0)), final['errors'],rtol=3.e-6,atol=1.e-20):
        raise ValueError('Reference worker error sum does not reproduce final errors')
    np.savez_compressed(output, contributions=values, individual_errors=errors, strata=strata,
                        random_seeds=np.asarray([streams[i] for i in range(len(workers))]),
                        titles=final['titles'], offsets=final['offsets'], edges=final['edges'])
    save(output.with_suffix('.json'),dict(created_utc=now(), input=str(record_path.resolve()),
                                         input_sha256=digest(record_path), arrays_sha256=digest(output),
                                         worker_count=len(workers), scale_points=POINTS,
                                         refinement_rng_audit=stage_audit,
                                         status='equal-count conditional reference batches audited; full-run convergence not certified'))
    return values, strata
