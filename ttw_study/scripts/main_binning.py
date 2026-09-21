"""Apply a common bin plan to complete independently retrained main vectors.

Conditional per-bin variances do not contain off-diagonal covariance, so
they cannot supply merged-bin variances. Keep rate diagnostics unchanged
and mark merged conditional diagonals unavailable. Physics errors are
recomputed from the rebinned full-run vectors by the main statistics code.
"""
import hashlib
import json

import numpy as np

from rebinning import spans

CONTINUOUS_IDS=tuple(local for local in range(2,22) if local not in (2,14))


def layout_signature(layout):
    record={key:np.asarray(layout[key]).tolist() for key in ('edges','offsets','titles','weights')}
    return hashlib.sha256(json.dumps(record,sort_keys=True,allow_nan=False,separators=(',',':')).encode()).hexdigest()


def apply_plan(layout,vectors,conditional,plan):
    if plan.get('schema')!='main_common_binning_v1' or plan.get('layout_sha256')!=layout_signature(layout):
        raise ValueError('Bin plan does not match the original main histogram layout')
    definitions=plan.get('histograms')
    if (not isinstance(definitions,dict) or not definitions or
            any(str(int(key))!=key or int(key) not in CONTINUOUS_IDS for key in definitions)):
        raise ValueError('Specify continuous local IDs 3--21 excluding 14; rates and categorical jet bins cannot be merged')
    titles,offsets,edges=layout['titles'],np.asarray(layout['offsets']),np.asarray(layout['edges'])
    if (len(titles)!=210 or offsets.shape!=(211,) or offsets[0]!=0 or
            offsets[-1]!=len(edges) or np.any(np.diff(offsets)<=0)):
        raise ValueError('Require the full ten-bank main histogram layout')
    if set(vectors)!=set(conditional):
        raise ValueError('Main vector and conditional-diagnostic groups differ')
    requested={int(key):np.asarray(value,dtype=float) for key,value in definitions.items()}
    intervals,parts,new_offsets,original_by_local,changes=[],[],[0],{},[]
    for h,title in enumerate(titles):
        local=h%21+1
        a,b=map(int,offsets[h:h+2])
        original=edges[a:b]
        if local==1 and (' rates:' not in title or b-a!=9):
            raise ValueError('A main bank is missing its nine unmodified parent rates')
        if local in requested:
            if local in original_by_local and not np.array_equal(original,original_by_local[local]):
                raise ValueError('A common observable has different original bins across cuts/charges')
            original_by_local[local]=original
            selected=spans(original,requested[local])
            parts.append(np.column_stack([requested[local][:-1],requested[local][1:]]))
            intervals.extend((a+first,a+last) for first,last in selected)
            changes.append(dict(title=title,local_id=local,original_bins=b-a,retained_bins=len(selected)))
        else:
            parts.append(original)
            intervals.extend((i,i+1) for i in range(a,b))
        new_offsets.append(len(intervals))
    if not any(last-first>1 for first,last in intervals):
        raise ValueError('Bin plan does not merge any original bin')
    rebinned,diagnostics={},{}
    for key,value in vectors.items():
        value=np.asarray(value)
        variance=np.asarray(conditional[key])
        if value.ndim!=3 or value.shape[1]!=len(edges) or variance.shape!=value.shape[:2]:
            raise ValueError('Full-run vectors or conditional diagonals have an inconsistent bin layout')
        target=np.empty((len(value),len(intervals),value.shape[2]))
        diagonal=np.full((len(value),len(intervals)),np.nan)
        for j,(first,last) in enumerate(intervals):
            if last-first==1:
                target[:,j]=value[:,first]
                diagonal[:,j]=variance[:,first]
            else:
                target[:,j]=value[:,first:last].sum(axis=1)
        rebinned[key]=target
        diagnostics[key]=diagonal
    updated=dict(layout,edges=np.concatenate(parts),offsets=np.asarray(new_offsets,dtype=int))
    return updated,rebinned,diagnostics,dict(original_bins=len(edges),retained_bins=len(intervals),
        common_across_all_cuts_charges_flavours_and_prescriptions=True,changed_histograms=changes,
        original_range_retained=True,continuous_overflows_folded=False,
        conditional_diagonal_policy='Parent rates and unmerged bins are unchanged. Merged-bin conditional '
            'variances are unavailable without within-run off-diagonal covariance and are marked NaN. '
            'All reported physics errors use the rebinned independently retrained vectors.',
        interpretation='This applies a declared common plan; it does not freeze bins or certify precision.')
