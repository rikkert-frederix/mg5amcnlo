#!/usr/bin/env python3
"""Absolute native/decayed normalization with matched scale/PDF covariance."""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from compare_current_batches import all_stage_pairs, batch_record, contrast
from inclusive_report import STATISTICAL, STABLE_VARIANT, input_bin, production_parameters, validate_bw_normalization
from joint_comparison import paired_variants, read
from native_inclusive import native_scope
from pilot_report import clean
from replica_statistics import ratio


def physical_settings(settings):
    result={k:v for k,v in settings.items() if k not in STATISTICAL}
    # Older standard manifests omit this unit-valued default; direct QES
    # tests state it explicitly. No nonunit QES change is ignored.
    result.setdefault('qes_over_ref',1.)
    return result


def match(native,decayed):
    nr,dr=native['report'],decayed['report']
    nm,dm=nr['manifest'],dr['manifest']
    generation=dr['execution']['export_generation']['args']
    kind=nm['native_kind']
    if (nr['scope']!=native_scope(kind) or nr['variant']!=STABLE_VARIANT[dr['variant']]
            or nm['charge']!=generation['charge'] or not nm['no_branching_factor_applied']
            or generation['corrected']!='both' or len(generation['flavours'])!=3
            or any(f not in ('e','mu') for f in generation['flavours'])):
        raise ValueError('Native and decayed process/order/normalization mismatch')
    expected_scale='fixed' if kind=='current' else 'core-w-ht-half'
    expected_modes=('all-bw',) if kind=='current' else ('onshell','top-bw')
    if (dm['w_treatment'] not in expected_modes or dm['production_scale']!=expected_scale
            or nm['production_scale']!=expected_scale
            or (kind=='current' and generation['flavours'][2]!=nm['associated_flavour'])
            or dm.get('custom_user_hook')):
        raise ValueError('Unmatched W treatment, associated flavour, scale or user hook')
    if physical_settings(nm['settings'])!=physical_settings(dm['settings']):
        raise ValueError('Native and decayed physical run settings differ')
    coordinates=dm.get('physical_scale_coordinates')
    if coordinates and (coordinates['numerator_factors']!=[1.]*4 or coordinates['qes_over_ref']!=1.):
        raise ValueError('Use a central technical control, not a direct varied-scale run')
    archives=[Path(r['path']).parents[2]/'study_cards'/r['manifest']['run_name'] for r in (nr,dr)]
    if production_parameters(archives[0]/'param_card.dat')!=production_parameters(archives[1]/'param_card.dat'):
        raise ValueError('Native and decayed production parameters differ')
    benchmark_path=STUDY/'inputs/benchmark.json'
    benchmark=json.loads(benchmark_path.read_text())
    if digest(benchmark_path)!=nm['benchmark_sha256']:
        raise ValueError('Native normalization benchmark changed')
    if dm['w_width']!=benchmark['W_width']['gamma_nlo_pdf']:
        raise ValueError('Unmatched physical W width')
    branching=benchmark['W_width']['branching_e']
    if dm['w_treatment']=='onshell':
        widths=[r for r in benchmark['top_widths'] if r['mb_GeV']==dm['decay_bottom_mass']
                and r['w_treatment']=='onshell' and r['scale_factor']==1.]
        if (len(widths)!=1 or dm['top_width_w_treatment']!='onshell'
                or dm['top_width_lo']!=widths[0]['gamma_lo']
                or dm['top_width_nlo']!=widths[0]['gamma_nlo_pdf']):
            raise ValueError('Unmatched on-shell top width normalization')
        width_proof=dict(path=str(benchmark_path),sha256=digest(benchmark_path))
    else:
        width_proof=validate_bw_normalization(dm,branching)
    if kind=='current' and nm['internal_W_width_GeV']!=dm['w_width']:
        raise ValueError('Native associated current has a different physical W width')
    return dict(native_kind=kind,charge=nm['charge'],decayed_variant=dr['variant'],
        native_variant=nr['variant'],w_treatment=dm['w_treatment'],flavours=generation['flavours'],
        branching_factor=branching**(2 if kind=='current' else 3),
        associated_branching_factor_applied=kind=='stable',width_normalization=width_proof,
        production_scale=expected_scale)


def compare(native_path,decayed_path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite absolute normalization evidence')
    native,decayed=read(native_path),read(decayed_path)
    physics=match(native,decayed)
    streams=[]
    for path,result in ((native_path,native),(decayed_path,decayed)):
        metadata,record=batch_record(path)
        pairs=all_stage_pairs(Path(result['report']['path']).parents[2],record)
        if len(pairs)!=metadata['refinement_rng_audit']['distinct_initialization_pairs']:
            raise ValueError('Native/decayed stage inventory differs from the batch audit')
        streams.append(pairs)
    if streams[0]&streams[1]:
        raise ValueError('Native and decayed calculations share an actual random stream')
    target=input_bin(decayed)
    result=paired_variants(native['contributions'][:,0,:],decayed['contributions'][:,target,:],
        native['strata'],decayed['strata'],lambda n,d: contrast(n,d,physics['branching_factor']))
    pulls=ratio(result['value'][2],result['mc_error'][2])
    np.savez_compressed(output.with_suffix('.npz'),**result,conditional_pulls=pulls,weights=native['weights'])
    record=dict(created_utc=now(),status='absolute normalization evaluated; inspect convergence and residuals',
        inputs={str(p.resolve()):digest(p) for p in (native_path,decayed_path)},physics=physics,
        stage_stream_counts=[len(s) for s in streams],cross_run_stream_collisions=0,
        arrays_sha256=digest(output.with_suffix('.npz')),script_sha256=digest(Path(__file__)),
        nominal=dict(expected_pb=result['value'][0,0],decayed_pb=result['value'][1,0],
            difference_pb=result['value'][2,0],difference_mc_error_pb=result['mc_error'][2,0],
            ratio=result['value'][3,0],ratio_mc_error=result['mc_error'][3,0],conditional_pull=pulls[0]),
        scale_maximum_absolute_conditional_pull=np.nanmax(np.abs(pulls[1:82])),
        pdf_maximum_absolute_conditional_pull=np.nanmax(np.abs(pulls[82:])),
        conditional_pulls=pulls,weight_labels=native['weights'],
        convention='Native total (including real/virtual/subtraction at NLO) times the matched branching '
                    'density versus the preselection decayed rate. All signed decay-scale coordinates '
                    'and 101 PDF members retained. Independent native/decayed stratum jackknife; '
                    'no extra associated-W branching factor for full currents.',
        limitations='Conditional on the trained integration proposals. Correlated scale/PDF pulls '
                    'are not independent trials or a global goodness-of-fit test. Full-flavour sums, '
                    'independent retraining and fiducial/tail precision remain separate checks.')
    save(output,clean(record))
    print(physics['charge'],physics['w_treatment'],physics['decayed_variant'],
          'ratio',record['nominal']['ratio'],'+/-',record['nominal']['ratio_mc_error'],flush=True)
    return clean(record)


def queue_report(queue_path,controls_path,companions_path,output):
    from run_normalization_checks import DONE
    from run_rng_pilots import DONE as CONTROLS_DONE
    from run_width_mass_pilots import DONE as COMPANIONS_DONE
    if output.exists():
        raise ValueError('Existing normalization summary')
    queue,controls,companions=[json.loads(p.read_text()) for p in (queue_path,controls_path,companions_path)]
    if (queue['status']!=DONE or len(queue['jobs'])!=14 or controls['status']!=CONTROLS_DONE
            or companions['status']!=COMPANIONS_DONE):
        raise ValueError('Normalization and companion queues must be complete')
    technical_path=Path(queue['after_queue'])
    if digest(technical_path)!=queue['predecessor_sha256']:
        raise ValueError('Normalization predecessor changed')
    technical=json.loads(technical_path.read_text())
    natives={}
    candidates=[]
    for job in queue['jobs']:
        case=job['case']
        if case['kind']=='decayed':
            candidates.append(job)
        else:
            key=case['kind'],case['charge'],case['variant']
            if key in natives:
                raise ValueError('Duplicate native reference; replicas require explicit combination')
            natives[key]=job
    candidates += [j for j in technical['jobs'] if
                   json.loads(Path(j['audit']).read_text())['manifest']['technical_test']['name']=='central_fixed']
    candidates += [j for source in (controls,companions) for j in source['jobs']
                   if j['variant'] in ('S','Pi')]
    reports=[]
    seen=set()
    for job in candidates:
        audit=json.loads(Path(job['audit']).read_text())
        manifest=audit['manifest']
        if (manifest['decay_bottom_mass']!=0. or
                (manifest['w_treatment']=='all-bw' and manifest['production_scale']!='fixed')):
            continue
        charge=audit['execution']['export_generation']['args']['charge']
        key=charge,manifest['w_treatment'],job['variant']
        if key in seen:
            raise ValueError('Duplicate decayed normalization reference')
        seen.add(key)
        kind='current' if manifest['w_treatment']=='all-bw' else 'stable'
        native=natives[(kind,charge,STABLE_VARIANT[job['variant']])]
        target=output.with_name(output.stem+'_%s_%s_%s.json'%key)
        report=compare(Path(native['batches']),Path(job['batches']),target)
        reports.append(dict(path=str(target),sha256=digest(target),physics=report['physics'],nominal=report['nominal']))
    expected={(charge,mode,variant) for charge in ('plus','minus')
              for mode in ('onshell','top-bw','all-bw')
              for variant in (('LO','P','S','Pi') if mode=='all-bw' else ('S','Pi'))}
    if seen!=expected:
        raise ValueError('Incomplete absolute normalization coverage')
    save(output,dict(created_utc=now(),status='sixteen absolute normalization comparisons evaluated; inspect agreement',
        inputs={str(p.resolve()):digest(p) for p in (queue_path,controls_path,companions_path,technical_path)},
        comparisons=reports,script_sha256=digest(Path(__file__)),
        note='Each pair retains covariance across all weights. Comparisons sharing a native reference '
             'are correlated; no combined chi-squared or ensemble-wide significance is inferred.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='action',required=True)
    p=sub.add_parser('pair')
    p.add_argument('--native',required=True,type=Path)
    p.add_argument('--decayed',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p=sub.add_parser('queue')
    p.add_argument('queue',type=Path)
    p.add_argument('--controls',required=True,type=Path)
    p.add_argument('--companions',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    if args.action=='pair':
        compare(args.native,args.decayed,args.output)
    else:
        queue_report(args.queue,args.controls,args.companions,args.output)
