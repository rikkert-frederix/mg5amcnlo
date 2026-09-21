#!/usr/bin/env python3
"""Archive the native template's current rates for joint-weight batch analysis.

This is a native inclusive-current audit, not a 210-histogram decay analysis.
The original native manifest/execution are never changed.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import ROOT, STUDY, digest, now, save
from harvest_splits import harvest
from load_results import canonical_columns
from madgraph.various.histograms import HwUList
from rng_history import history_paths, validate_history


def native_vectors(histograms):
    if ([h.title.strip() for h in histograms]!=['total rate','total rate Born']
            or any(len(h.bins)!=5 for h in histograms)):
        raise ValueError('Expected the two five-bin native template histograms')
    rows=[row for h in histograms for row in h.bins]
    edges=np.asarray([row.boundaries for row in rows])
    if not np.array_equal(edges,[(i+.5,i+1.5) for _ in range(2) for i in range(5)]):
        raise ValueError('Native template bin boundaries changed')
    columns=list(rows[0].wgts)
    indices=canonical_columns(columns,'LO')
    values=np.asarray([[row.wgts[c] for c in columns] for row in rows])
    if (not np.isfinite(values).all() or np.any(values[:,columns.index('dy')]<0.)
            or np.any(values[[1,2,3,4,6,7,8,9]]) or not np.array_equal(values[0],values[5])):
        raise ValueError('Expected identical LO total/Born rates and empty unused template bins')
    return dict(values=values,columns=columns,indices=indices,edges=edges,
                titles=[h.title for h in histograms],offsets=[0,5,10])


def audit(execution_path,output,workers):
    execution_path=execution_path.resolve()
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a native-current audit')
    execution=json.loads(execution_path.read_text())
    manifest_path=execution_path.parent/'manifest.json'
    manifest=json.loads(manifest_path.read_text())
    hwu=Path(execution['output'])
    process=hwu.parents[2]
    if (process.parent!=STUDY/'processes' or execution['status']!='finished'
            or manifest['variant']!='LO' or digest(hwu)!=execution['output_sha256']
            or not manifest['no_branching_factor_applied']
            or (process/'Cards/decay_card.dat').exists()):
        raise ValueError('Require an unchanged completed undecayed LO current reference')
    for name,sha in manifest['files'].items():
        if digest(execution_path.parent/name)!=sha:
            raise ValueError('Native current card changed')
    if digest(process/'FixedOrderAnalysis/HwU.f90')!=digest(ROOT/'Template/fNLO/FixedOrderAnalysis/HwU.f90'):
        raise ValueError('Native current accumulator differs from the corrected runtime')
    histograms=HwUList(str(hwu),raw_labels=True)
    vectors=native_vectors(histograms)
    if (not np.array_equal(vectors['values'][0,vectors['indices']],execution['canonical_weights_pb'])
            or vectors['values'][0,vectors['columns'].index('dy')]!=execution['nominal_mc_error_pb']):
        raise ValueError('Native current weights/error differ from the frozen execution')
    history={str(p.relative_to(process)):p.read_bytes() for p in history_paths(process,hwu.parent.name)}
    proof=validate_history(Path(execution['log']).read_text(),manifest['settings']['iseed'],history)
    if (proof!=execution['rng_history_audit'] or
            any(digest(process/name)!=sha for name,sha in execution['rng_history_hashes'].items())):
        raise ValueError('Native current RNG history changed')
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'),**{k:v for k,v in vectors.items() if k!='indices'})
    record=dict(created_utc=now(),path=str(hwu),output_sha256=digest(hwu),variant='LO',
                status='native inclusive-current technical audit passed; convergence not certified',
                corrected_accumulator=True,histogram_count=2,scale_count=9,pdf_count=101,
                arrays_sha256=digest(output.with_suffix('.npz')),
                manifest=dict(manifest,hashes=manifest['files']),execution=execution,
                original_inputs={str(p):digest(p) for p in (execution_path,manifest_path)},
                scope='Undecayed stable-top plus associated current only; no fiducial decay observables')
    save(output,record)
    harvest(output,workers)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('execution',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--workers',required=True,type=Path)
    args=parser.parse_args()
    audit(args.execution,args.output,args.workers)
