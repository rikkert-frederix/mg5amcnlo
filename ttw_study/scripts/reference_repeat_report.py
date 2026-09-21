#!/usr/bin/env python3
"""Compare independent matched-reference retrainings without pooling shared data."""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import digest, now, save
from compare_current_batches import all_stage_pairs, batch_record
from inclusive_report import STATISTICAL
from joint_comparison import paired_variants
from pilot_report import clean
from reference_batches import POINTS, vectors
from replica_statistics import ratio


def read(path):
    metadata,record=batch_record(path)
    audit_path=Path(record['audit'])
    audit=json.loads(audit_path.read_text())
    if (digest(path)!=metadata['arrays_sha256'] or digest(metadata['input'])!=metadata['input_sha256']
            or digest(audit_path)!=record['audit_sha256']
            or digest(audit['path'])!=audit['output_sha256']):
        raise ValueError('Reference repeat input changed')
    hwu=Path(audit['path'])
    archive=hwu.parents[2]/'study_cards'/hwu.parent.name
    if any(digest(archive/name)!=sha for name,sha in audit['manifest']['files'].items()):
        raise ValueError('Reference repeat card archive changed')
    final=vectors(hwu.read_bytes())
    streams=all_stage_pairs(hwu.parents[2],record)
    if len(streams)!=metadata['refinement_rng_audit']['distinct_initialization_pairs']:
        raise ValueError('Reference repeat stage inventory differs from its audit')
    with np.load(path) as data:
        result=dict(contributions=data['contributions'].copy(),strata=data['strata'].copy())
        if list(data['titles'])!=final['titles'] or not np.array_equal(data['edges'],final['edges']):
            raise ValueError('Reference repeat layout changed')
    result.update(audit=audit,streams=streams,**{k:final[k] for k in ('titles','edges','offsets')})
    return result


def compare(first_path,second_path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a reference convergence comparison')
    a,b=read(first_path),read(second_path)
    manifests=[row['audit']['manifest'] for row in (a,b)]
    signature=lambda m:{**{k:m[k] for k in ('variant','charge','flavours','physics','production_scale_GeV',
        'decay_numerator_scale_GeV','scale_variation')},
        'settings':{k:v for k,v in m['settings'].items() if k not in STATISTICAL}}
    if signature(manifests[0])!=signature(manifests[1]):
        raise ValueError('Reference retrainings do not have matched physical inputs')
    if a['titles']!=b['titles'] or not np.array_equal(a['edges'],b['edges']):
        raise ValueError('Reference retrainings have different observables')
    if a['streams']&b['streams']:
        raise ValueError('Reference estimates share actual stage streams; not independent retrainings')
    transform=lambda x,y:np.asarray([x,y,y-x,ratio(y,x)])
    result=paired_variants(a['contributions'],b['contributions'],a['strata'],b['strata'],transform)
    pulls=ratio(result['value'][2],result['mc_error'][2])
    slot=1 if manifests[0]['charge']=='plus' else 3
    nominal=result['value'][:,slot,0]
    errors=result['mc_error'][:,slot,0]
    np.savez_compressed(output.with_suffix('.npz'),**result,conditional_pulls=pulls,
                        titles=a['titles'],edges=a['edges'],offsets=a['offsets'])
    report=dict(created_utc=now(),status='independent matched reference repeats compared; retraining coverage not certified',
        inputs={str(p.resolve()):digest(p) for p in (first_path,second_path)},
        arrays_sha256=digest(output.with_suffix('.npz')),script_sha256=digest(Path(__file__)),
        variant=manifests[0]['variant'],charge=manifests[0]['charge'],flavours=manifests[0]['flavours'],
        stage_stream_counts=[len(a['streams']),len(b['streams'])],cross_run_stream_overlap=0,
        nominal_fiducial=dict(first_ab=nominal[0]*1.e6,first_mc_error_ab=errors[0]*1.e6,
            second_ab=nominal[1]*1.e6,second_mc_error_ab=errors[1]*1.e6,
            difference_ab=nominal[2]*1.e6,difference_mc_error_ab=errors[2]*1.e6,
            ratio=nominal[3],ratio_mc_error=errors[3],conditional_pull=pulls[slot,0]),
        common_scale_points=['nominal']+POINTS,fiducial_conditional_pulls=pulls[slot],
        maximum_absolute_fiducial_scale_pull=np.max(np.abs(pulls[slot,1:])),
        convention='Second minus first and second/first from independent run/beam-stratum ensembles. '
                   'All common-scale coordinates retained. The samples are neither summed as flavours '
                   'nor treated as arbitrarily paired events. Both original estimates remain available.',
        limitations='Two independent retrainings do not establish stable between-retraining covariance. '
                    'The diagnostic denominator uses each run\'s conditional batch uncertainty; '
                    'correlated scale pulls are not independent trials or a global goodness-of-fit test.')
    save(output,clean(report))
    print(json.dumps(clean(report['nominal_fiducial']),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('first',type=Path)
    parser.add_argument('second',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    compare(args.first,args.second,args.output)
