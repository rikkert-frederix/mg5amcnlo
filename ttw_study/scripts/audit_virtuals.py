#!/usr/bin/env python3
"""Verify analytic top/antitop virtual checks from immutable worker archives."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import tarfile

from campaign import digest, now, save


CHECK = re.compile(r'analytic top-decay virtual validation (\d+)/3 for contribution (\d+); '
                   r'relative difference\s+(\S+)')
DONE = re.compile(r'analytic top-decay virtual provider for contribution (\d+) validated;')


def check_log(log, contributions=(2, 3), tolerance=1.e-8):
    points = {}
    for point, contribution, raw in CHECK.findall(log):
        key = (int(contribution), int(point))
        value = float(raw.replace('D', 'E'))
        if key in points or not math.isfinite(value) or not 0. <= value <= tolerance:
            raise ValueError('Duplicate or failed analytic virtual validation')
        points[key] = value
    expected = {(c, i) for c in contributions for i in (1, 2, 3)}
    if set(points) != expected or sorted(map(int, DONE.findall(log))) != sorted(contributions):
        raise ValueError('Incomplete top/antitop analytic virtual validation')
    if 'Time spent in Total :' not in log:
        raise ValueError('Worker did not complete')
    return {str(c): max(points[c, i] for i in (1, 2, 3)) for c in contributions}


def audit(record_path, output):
    if output.exists():
        raise ValueError('Refuse to overwrite a virtual validation audit')
    record = json.loads(record_path.read_text())
    archive_path = record_path.with_suffix('.gz')
    if digest(archive_path) != record['archive_sha256']:
        raise ValueError('Worker archive checksum changed')
    if digest(record['audit']) != record['audit_sha256']:
        raise ValueError('Combined-output audit changed')
    expected = {w['directory']+'/log.txt': w['files']['log.txt'] for w in record['workers']}
    if len(expected) != record['worker_count']:
        raise ValueError('Duplicate worker directories')
    rows = {}
    with tarfile.open(archive_path, 'r|gz') as archive:
        for member in archive:
            if not member.name.endswith('/log.txt'):
                continue
            raw = archive.extractfile(member).read()
            if hashlib.sha256(raw).hexdigest() != expected.pop(member.name, None):
                raise ValueError('Unexpected, repeated or changed worker log')
            rows[member.name] = check_log(raw.decode())
    if expected:
        raise ValueError('Missing worker logs')
    result = dict(created_utc=now(), input=str(record_path.resolve()),
                  input_sha256=digest(record_path), worker_count=len(rows),
                  checked_points=6*len(rows), tolerance=1.e-8,
                  maximum_relative_difference=max(v for r in rows.values() for v in r.values()),
                  workers=rows, status='all archived top/antitop analytic virtual checks passed',
                  limitations='Three native validation points per corrected branch and worker; '
                              'not a narrow-width, full-virtuality-support or convergence test.')
    save(output, result)
    print('%d checks; maximum relative difference %.6g' %
          (result['checked_points'], result['maximum_relative_difference']), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    audit(args.record, args.output)
