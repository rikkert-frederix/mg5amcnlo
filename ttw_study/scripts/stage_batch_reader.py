"""Joint equal-count worker vectors with an explicit histogram reader.

This adapter lets new weight layouts use the established estimator and
complete training/refinement stream audit without modifying frozen queues.
It never turns conditional split batches into independent retrainings.
"""
import hashlib
import json
from pathlib import Path
import re
import tarfile

import numpy as np

from campaign import ROOT, digest, now, save
from read_splits import validate_jobs
from rng_history import validate_history


def read(record_path, output, reference, vector_reader):
    record_path, output = Path(record_path).resolve(), Path(output).resolve()
    if output.exists() or output.with_suffix('.json').exists():
        raise ValueError('Refuse to overwrite a joint batch archive')
    record = json.loads(record_path.read_text())
    archive_path = record_path.with_suffix('.gz')
    if (digest(archive_path) != record['archive_sha256']
            or digest(record['audit']) != record['audit_sha256']):
        raise ValueError('Worker archive or final-output audit changed')
    audit = reference['report']
    if audit != json.loads(Path(record['audit']).read_text()):
        raise ValueError('Reference vectors belong to a different output audit')
    if digest(audit['path']) != audit['output_sha256']:
        raise ValueError('Original combined histogram changed')
    if not record.get('history_files'):
        raise ValueError('Require complete fresh-training and refinement stream evidence')
    workers = record['workers']
    if record['worker_count'] != len(workers):
        raise ValueError('Worker count differs from the archive inventory')
    groups = validate_jobs(workers)
    directories = {worker['directory']: i for i, worker in enumerate(workers)}
    if len(directories) != len(workers):
        raise ValueError('Repeated worker directory')
    expected = {worker['directory']+'/'+name: sha for worker in workers
                for name, sha in worker['files'].items()}
    if set(expected) & set(record['history_files']):
        raise ValueError('History entries overlap final-worker files')
    expected.update(record['history_files'])
    history, vectors, errors, seeds, initial = {}, {}, {}, {}, {}
    with tarfile.open(archive_path, 'r|gz') as archive:
        for member in archive:
            if not member.isfile() or member.name not in expected:
                raise ValueError('Unexpected, duplicate or non-file archive entry: '+member.name)
            raw = archive.extractfile(member).read()
            if hashlib.sha256(raw).hexdigest() != expected.pop(member.name):
                raise ValueError('Worker file checksum mismatch: '+member.name)
            if member.name in record['history_files']:
                history[member.name] = raw
                continue
            directory, filename = member.name.rsplit('/', 1)
            index = directories[directory]
            if filename == 'log.txt':
                log = raw.decode()
                found = re.findall(r'Ranmar initialization seeds\s+(\d+)\s+(\d+)', log)
                if len(found) != 1 or re.findall(r'with seed\s+(\d+)', log) != [str(record['seed'])]:
                    raise ValueError('Unverified final-worker random stream')
                seeds[index] = tuple(map(int, found[0]))
                initial[index] = re.findall(r'^channel.*$', log.split('------- iteration')[0], re.M)
                if not initial[index]:
                    raise ValueError('Missing initial proposal/channel diagnostics')
            elif filename == 'MADatNLO.HwU':
                row = vector_reader(raw)
                if (row['titles'] != reference['titles'] or row['weights'] != reference['weights']
                        or not np.array_equal(row['edges'], reference['edges'])
                        or not np.array_equal(row['offsets'], reference['offsets'])):
                    raise ValueError('Worker histogram or weight layout changed')
                vectors[index], errors[index] = row['values'], row['errors']
                if len(vectors) % 16 == 0:
                    print('Parsed %d/%d joint batch vectors' % (len(vectors), len(workers)), flush=True)
    if (expected or len(vectors) != len(workers) or len(errors) != len(workers)
            or len(seeds) != len(workers) or len(initial) != len(workers)):
        raise ValueError('Incomplete worker archive')
    if len(set(seeds.values())) != len(workers):
        raise ValueError('Repeated final-worker random stream')
    run_log = Path(audit['execution']['log']).read_text()
    stage_audit = validate_history(run_log, record['seed'], history,
        {worker['directory']: worker['files']['log.txt'] for worker in workers})
    strata = np.empty(len(workers), dtype=int)
    for label, indices in enumerate(groups.values()):
        if any(initial[i] != initial[indices[0]] for i in indices):
            raise ValueError('Initial proposals differ within one stratum')
        strata[indices] = label
    values = np.asarray([vectors[i] for i in range(len(workers))])
    errors = np.asarray([errors[i] for i in range(len(workers))])
    difference = values.sum(axis=0)-reference['values']
    tolerance = 3.e-6*(np.abs(values).sum(axis=0)+np.abs(reference['values']))+1.e-20
    if not (np.abs(difference) <= tolerance).all():
        raise ValueError('Worker sum does not reproduce every combined weight/bin')
    error_sum = np.sqrt(np.sum(errors**2, axis=0))
    if not np.allclose(error_sum, reference['errors'], rtol=3.e-6, atol=1.e-20):
        raise ValueError('Worker errors do not reproduce combined nominal errors')
    process = Path(audit['path']).parents[2]
    sources = {}
    for relative in ('SubProcesses/mint_module.f90', 'Source/ranmar.f90'):
        sources[relative] = digest(process/relative)
        if sources[relative] != digest(ROOT/'Template/fNLO'/relative):
            raise ValueError('Generated estimator/RNG source differs from its audit')
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, contributions=values, individual_errors=errors,
        strata=strata, random_seeds=np.asarray([seeds[i] for i in range(len(workers))]),
        edges=reference['edges'], offsets=reference['offsets'],
        titles=reference['titles'], weights=reference['weights'])
    save(output.with_suffix('.json'), dict(created_utc=now(), input=str(record_path),
        input_sha256=digest(record_path), arrays_sha256=digest(output), estimator_sources=sources,
        refinement_rng_audit=stage_audit, worker_count=len(workers),
        strata={str(key): len(indices) for key, indices in groups.items()},
        maximum_weight_sum_difference_pb=float(np.abs(difference).max()),
        maximum_error_sum_difference_pb=float(np.abs(error_sum-reference['errors']).max()),
        reader_sha256=digest(Path(__file__)),
        status='equal-count conditional joint batch audit passed; full-run convergence not certified',
        convention='Keep all selected weight/bin correlations and independent strata. The final '
                   'batches share trained proposals; no independent-retraining or rare-tail coverage claim.'))
    print(output, flush=True)
