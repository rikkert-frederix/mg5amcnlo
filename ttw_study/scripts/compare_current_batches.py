#!/usr/bin/env python3
"""Joint-weight LO current normalization, after the independent physics audit."""
import argparse
import json
from pathlib import Path
import re

import numpy as np

from campaign import ROOT, digest, now, save
from inclusive_report import input_bin
from joint_comparison import read, paired_variants
from pilot_report import clean
from replica_statistics import ratio
from rng_history import STAGE_FILE


def batch_record(path):
    metadata=json.loads(path.with_suffix('.json').read_text())
    record=json.loads(Path(metadata['input']).read_text())
    if not metadata.get('refinement_rng_audit'):
        raise ValueError('Both current comparisons require full staged RNG evidence')
    return metadata,record


def all_stage_pairs(process,record):
    pairs=[]
    for name,sha in record['history_files'].items():
        path=process/name
        if digest(path)!=sha:
            raise ValueError('Current-comparison stage evidence changed')
        if STAGE_FILE.fullmatch(name):
            pairs.extend(tuple(map(int,p)) for p in re.findall(
                r'Ranmar initialization seeds\s+(\d+)\s+(\d+)',path.read_text()))
    if not pairs or len(pairs)!=len(set(pairs)):
        raise ValueError('Missing or repeated current-comparison stage streams')
    return set(pairs)


def contrast(native,decayed,branching_factor):
    expected=native*branching_factor
    return np.asarray([expected,decayed,decayed-expected,ratio(decayed,expected)])


def absolute_inputs(inputs):
    result={}
    for name,sha in inputs.items():
        path=Path(name)
        path=path.resolve() if path.is_absolute() else (ROOT/path).resolve()
        if str(path) in result:
            raise ValueError('Duplicate physical input path')
        result[str(path)]=sha
    return result


def compare(native_path,decayed_path,physical_report,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite the joint current comparison')
    native,decayed=read(native_path),read(decayed_path)
    nm,nr=batch_record(native_path)
    dm,dr=batch_record(decayed_path)
    proof=json.loads(physical_report.read_text())
    native_inputs=absolute_inputs(native['report']['original_inputs'])
    physical_inputs=absolute_inputs(proof['inputs'])
    native_execution=[p for p in native_inputs if Path(p).name=='execution.json']
    if (len(native_execution)!=1 or set(physical_inputs)!={native_execution[0],str(Path(dr['audit']).resolve())}
            or any(digest(p)!=sha for p,sha in physical_inputs.items())
            or any(digest(p)!=sha for p,sha in native_inputs.items())
            or proof['associated_branching_factor_applied']):
        raise ValueError('Batches do not match the independently audited physical normalization')
    factor=proof['branching_factor']
    target=input_bin(decayed)
    if (native['report']['scope']!='Undecayed stable-top plus associated current only; no fiducial decay observables'
            or not np.array_equal(decayed['values'][target],proof['decayed_pb'])
            or not np.array_equal(native['values'][0]*factor,proof['expected_pb'])):
        raise ValueError('Native/decayed canonical weights differ from the physical comparison')
    processes=[Path(r['report']['path']).parents[2] for r in (native,decayed)]
    streams=[all_stage_pairs(p,r) for p,r in zip(processes,(nr,dr))]
    if streams[0]&streams[1]:
        raise ValueError('Random-stream collision between independent current calculations')
    result=paired_variants(native['contributions'][:,0,:],decayed['contributions'][:,target,:],
                           native['strata'],decayed['strata'],lambda a,b:contrast(a,b,factor))
    pull=ratio(result['value'][2],result['mc_error'][2])
    np.savez_compressed(output.with_suffix('.npz'),values=result['value'],mc_errors=result['mc_error'],
                        conditional_pulls=pull,weights=native['weights'])
    record=dict(created_utc=now(),status='independent LO current joint-weight pilot; full-run convergence not certified',
                inputs={str(p.resolve()):digest(p) for p in (native_path,decayed_path,physical_report)},
                script_sha256=digest(Path(__file__)),arrays_sha256=digest(output.with_suffix('.npz')),
                branching_factor=factor,associated_branching_factor_applied=False,
                stage_stream_counts=[len(s) for s in streams],cross_run_stream_collisions=0,
                nominal=dict(expected_pb=result['value'][0,0],decayed_pb=result['value'][1,0],
                             difference_pb=result['value'][2,0],difference_mc_error_pb=result['mc_error'][2,0],
                             ratio=result['value'][3,0],ratio_mc_error=result['mc_error'][3,0],
                             conditional_pull=pull[0]),
                scale_maximum_absolute_conditional_pull=np.max(np.abs(pull[1:82])),
                pdf_maximum_absolute_conditional_pull=np.max(np.abs(pull[82:])),
                conditional_pulls=pull,weight_labels=native['weights'],
                limitations='One ordered all-bw LO assignment at matched fixed production scales. '
                            'Independent stratum batches and paired weights; conditional on learned proposals. '
                            'Correlated pulls are not a global goodness-of-fit test or sparse-tail coverage proof. '
                            'Not a narrow-W, NLO, dynamic-scale or full-flavour validation. '
                            'Native requested/effective cell-rounded counts remain separately archived.')
    save(output,clean(record))
    print('Joint current nominal pull:',pull[0],'; maximum scale/PDF pulls:',
          record['scale_maximum_absolute_conditional_pull'],record['pdf_maximum_absolute_conditional_pull'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native',required=True,type=Path)
    parser.add_argument('--decayed',required=True,type=Path)
    parser.add_argument('--physical-report',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    compare(args.native,args.decayed,args.physical_report,args.output)
