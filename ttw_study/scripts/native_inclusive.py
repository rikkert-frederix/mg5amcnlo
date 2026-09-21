#!/usr/bin/env python3
"""Audited LO/NLO native stable-ttW and associated-current total rates.

The NLO total includes real, virtual and subtraction contributions. The
separately booked Born histogram is retained, never substituted for it.
"""
import json
from pathlib import Path

import numpy as np

from campaign import ROOT, STUDY, digest, now, save
from harvest_splits import harvest
from load_results import canonical_columns
from rng_history import history_paths, validate_history


def native_scope(kind):
    if kind=='current':
        return 'Undecayed stable-top plus associated current only; no fiducial decay observables'
    if kind=='stable':
        return 'Undecayed stable t,tbar,W only; no branching factors or fiducial decay observables'
    raise ValueError('Unknown native inclusive process')


def native_vectors(histograms,variant):
    if variant not in ('LO','P'):
        raise ValueError('Native total rates support LO or production NLO only')
    if ([h.title.strip() for h in histograms]!=['total rate','total rate Born']
            or any(len(h.bins)!=5 for h in histograms)):
        raise ValueError('Expected the two five-bin native template histograms')
    rows=[row for h in histograms for row in h.bins]
    edges=np.asarray([row.boundaries for row in rows])
    if not np.array_equal(edges,[(i+.5,i+1.5) for _ in range(2) for i in range(5)]):
        raise ValueError('Native template bin boundaries changed')
    columns=list(rows[0].wgts)
    if any(list(row.wgts)!=columns for row in rows):
        raise ValueError('Native weight labels differ between bins')
    indices=canonical_columns(columns,variant)
    values=np.asarray([[row.wgts[c] for c in columns] for row in rows])
    if (not np.isfinite(values).all() or np.any(values[:,columns.index('dy')]<0.)
            or np.any(values[[1,2,3,4,6,7,8,9]])):
        raise ValueError('Nonfinite weights, negative errors or nonempty unused native bins')
    if variant=='LO' and not np.array_equal(values[0],values[5]):
        raise ValueError('LO total and Born histograms differ')
    return dict(values=values,columns=columns,indices=indices,edges=edges,
                titles=[h.title for h in histograms],offsets=[0,5,10])


def audit(execution_path,output,workers):
    from madgraph.various.histograms import HwUList
    from models.check_param_card import ParamCard
    execution_path=execution_path.resolve()
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a native inclusive audit')
    execution=json.loads(execution_path.read_text())
    manifest_path=execution_path.parent/'manifest.json'
    manifest=json.loads(manifest_path.read_text())
    variant,kind=manifest['variant'],manifest['native_kind']
    scope=native_scope(kind)
    hwu=Path(execution['output'])
    process=hwu.parents[2]
    if (process.parent!=STUDY/'processes' or execution['status']!='finished'
            or execution['returncode']!=0 or execution['variant']!=variant
            or digest(hwu)!=execution['output_sha256']
            or not manifest['no_branching_factor_applied']
            or (process/'Cards/decay_card.dat').exists()):
        raise ValueError('Require a completed undecayed native inclusive reference')
    for name,sha in manifest['files'].items():
        if digest(execution_path.parent/name)!=sha:
            raise ValueError('Native inclusive input card changed')
    if (digest(STUDY/'inputs/benchmark.json')!=manifest['benchmark_sha256']
            or digest(manifest['reference_manifest'])!=manifest['reference_manifest_sha256']):
        raise ValueError('Native benchmark/reference changed')
    for relative in ('FixedOrderAnalysis/HwU.f90','FixedOrderAnalysis/analysis_HwU_template.f90'):
        if digest(process/relative)!=manifest['source_hashes']['Template/fNLO/'+relative]:
            raise ValueError('Native histogram code differs from frozen source')
    card=ParamCard(str(execution_path.parent/'param_card.dat'))
    benchmark=json.loads((STUDY/'inputs/benchmark.json').read_text())
    expected_w=benchmark['W_width']['gamma_nlo_pdf'] if kind=='current' else 0.
    if card['decay'].get((6,)).value!=0. or card['decay'].get((24,)).value!=expected_w:
        raise ValueError('Native external/internal width convention differs')
    vectors=native_vectors(HwUList(str(hwu),raw_labels=True),variant)
    if (not np.array_equal(vectors['values'][0,vectors['indices']],execution['canonical_weights_pb'])
            or vectors['values'][0,vectors['columns'].index('dy')]!=execution['nominal_mc_error_pb']):
        raise ValueError('Native total weights/error differ from frozen execution')
    history={str(p.relative_to(process)):p.read_bytes() for p in history_paths(process,hwu.parent.name)}
    proof=validate_history(Path(execution['log']).read_text(),manifest['settings']['iseed'],history)
    if proof!=execution['rng_history_audit']:
        raise ValueError('Native stage history differs from frozen execution')
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'),**{k:v for k,v in vectors.items() if k!='indices'})
    save(output,dict(created_utc=now(),path=str(hwu),output_sha256=digest(hwu),variant=variant,
        status='native inclusive technical audit passed; convergence not certified',
        corrected_accumulator=True,histogram_count=2,scale_count=9,pdf_count=101,
        arrays_sha256=digest(output.with_suffix('.npz')),
        manifest=dict(manifest,hashes=manifest['files']),execution=execution,
        original_inputs={str(p):digest(p) for p in (execution_path,manifest_path)},scope=scope,
        script_sha256=digest(Path(__file__))))
    harvest(output,workers)
