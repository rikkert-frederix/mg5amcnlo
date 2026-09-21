#!/usr/bin/env python3
"""Matched fixed-mt, inclusive dilepton ttbar spin-correlation calibration.

Separate NNPDF3.1/reference inputs, not the main ttW benchmark. The reference
files labelled 'exp' expand the normalized distribution; unlabelled files
divide strict NLO differential by strict NLO total. Neither label refers to
an unexpanded top-width denominator. Keep these two expansions distinct.
"""
import argparse
import fcntl
import json
import math
from pathlib import Path
import shutil
import sys
import time

import numpy as np

from campaign import ROOT, STUDY, RUNTIME_SOURCES, benchmark_param_card, digest, now, save
from audit_virtuals import audit as audit_virtuals
from compare_current_batches import all_stage_pairs, batch_record
from harvest_splits import harvest
from pilot_report import clean
from read_splits import split_ensembles
from reference_batches import POINTS, read, vectors
from replica_statistics import independent_jackknife_many, ratio
from run_inclusive import launch

ANALYSIS = 'analysis_HwU_pp_ttx_leptons_reference.f90'


def run_settings(baseline, physics, seed, accuracy, job_seconds):
    if (not isinstance(seed,int) or not 0<seed<=900000000 or
            not math.isfinite(job_seconds) or job_seconds<=0.):
        raise ValueError('Invalid ttbar sampling settings')
    settings = dict(baseline['settings'])
    settings.update(lhaid=303400,iseed=seed,req_acc_fo=accuracy,fo_job_target_time=job_seconds,
                    fixed_ren_scale=True,fixed_fac_scale=True,fixed_qes_scale=True,
                    mur_ref_fixed=physics['mt'],muf_ref_fixed=physics['mt'],qes_ref_fixed=physics['mt'],
                    reweight_scale=True,reweight_pdf=False)
    return settings


def decay_card(variant, physics):
    from madgraph.fks.fks_decay import decay_card_text
    if variant not in ('LO','S','Pi'):
        raise ValueError('ttbar reference uses LO auxiliary, S and Pi')
    g0, gn = physics['top_width_LO'], physics['top_width_NLO']
    return decay_card_text(
        {6:g0,24:physics['W_width']},{6:physics['mt'],24:physics['MW']},
        nlo_width_pdgs={6},nlo_widths={6:gn},
        production_order='LO' if variant == 'LO' else 'NLO',
        decay_order='LO' if variant == 'LO' else 'NLO',
        nlo_decay_combination='MULTIPLICATIVE' if variant == 'Pi' else 'ADDITIVE',
        decay_scale_variation_mode='CORRELATED',decay_scale_factors=(1.,.5,2.),
        decay_width_scale_modes={6:'EXPLICIT'},
        lo_width_variations={(6,f):g0 for f in (.5,2.)},
        nlo_width_variations={(6,f):gn for f in (.5,2.)},production_scale_momenta='CORE')


def normalized(total):
    """Bin densities in DeltaPhi/pi, normalized to the independent rates bin."""
    return ratio(total[2:12],total[0])/.1


def expanded_normalized(born, strict):
    return (ratio(strict[2:12],born[0])-
            ratio(born[2:12],born[0])*ratio(strict[0]-born[0],born[0]))/.1


def check_histogram(row):
    if row['offsets'].tolist() != [0,2,12,22]:
        raise ValueError('Unexpected ttbar reference histogram layout')
    value = row['values']
    for first,last,rate in ((2,12,0),(12,22,1)):
        tolerance = 3.e-6*(np.abs(value[first:last]).sum(axis=0)+np.abs(value[rate]))+1.e-20
        if not (np.abs(value[first:last].sum(axis=0)-value[rate]) <= tolerance).all():
            raise ValueError('Angular integral does not equal its normalization rate')
    if not np.allclose(value[:,0],value[:,1],rtol=2.e-7,atol=1.e-20):
        raise ValueError('Nominal/common central scale identity failed')


def parameters(process, preflight):
    import lhapdf
    from models import import_ufo, model_reader
    from models.check_param_card import ParamCard
    physics = preflight['physics']
    name = preflight['PDF']['name']
    directory = ROOT/'LHAPDF/share/LHAPDF'/name
    if digest(directory/(name+'.info')) != preflight['PDF']['metadata_sha256']:
        raise ValueError('ttbar reference PDF metadata changed')
    for filename,checksum in preflight['PDF']['member_sha256'].items():
        if digest(directory/filename) != checksum:
            raise ValueError('ttbar reference PDF changed')
    benchmark_param_card(process,json.loads((STUDY/'inputs/benchmark.json').read_text()),'onshell')
    card = ParamCard(str(process/'Cards/param_card.dat'))
    gf,mw,mz = (physics[k] for k in ('GF','MW','MZ'))
    alpha_inverse = 1./(math.sqrt(2)*gf*mw**2/math.pi*(1-mw**2/mz**2))
    pdf = lhapdf.mkPDF(name,0)
    for block,pdg,value in (('sminputs',1,alpha_inverse),('sminputs',2,gf),
                            ('sminputs',3,pdf.alphasQ(mz)),
                            ('decay',6,physics['top_width_NLO']),('decay',24,physics['W_width'])):
        card[block].get((pdg,)).value = value
    model = model_reader.ModelReader(import_ufo.import_model('loop_sm-no_b_mass'))
    card.update_dependent(model,None,20)
    model.set_parameters_and_couplings(card)
    if abs(model.get_mass(24)-mw)>1.e-10 or card['mass'].get((5,)).value != 0.:
        raise ValueError('ttbar reference model masses do not match')
    card.write(str(process/'Cards/param_card.dat'),precision=16)


def compare(batches, preflight, output):
    ensembles, streams = {}, set()
    stage_counts = {}
    for variant,path in batches.items():
        metadata, record = batch_record(path)
        audit = json.loads(Path(record['audit']).read_text())
        if (digest(path) != metadata['arrays_sha256'] or
                digest(metadata['input']) != metadata['input_sha256'] or
                digest(record['audit']) != record['audit_sha256'] or
                digest(audit['path']) != audit['output_sha256']):
            raise ValueError('ttbar reference batches changed')
        ensemble_streams = all_stage_pairs(Path(audit['path']).parents[2],record)
        if (len(ensemble_streams) != metadata['refinement_rng_audit']['distinct_initialization_pairs'] or
                streams.intersection(ensemble_streams)):
            raise ValueError('Missing or overlapping ttbar training/refinement streams')
        streams.update(ensemble_streams)
        stage_counts[variant] = len(ensemble_streams)
        with np.load(path) as data:
            ensembles[variant] = split_ensembles(data['contributions'],data['strata'])
    reports = {}
    order = [(1,1),(2,2),(.5,.5),(1,.5),(.5,1),(2,1),(1,2)]
    indices = [1+POINTS.index(point) for point in order]
    nborn = len(ensembles['LO'])
    for label, expanded in (('NLO_normalized',False),('NLO_expanded_normalized',True)):
        if expanded:
            estimate = independent_jackknife_many(
                ensembles['LO']+ensembles['S'],
                lambda *means: expanded_normalized(sum(means[:nborn]),sum(means[nborn:])))
        else:
            estimate = independent_jackknife_many(ensembles['S'],lambda *means: normalized(sum(means)))
        filename = 'norm-incl-'+('exp-' if expanded else '')+'dPhill-NLO-NNPDF-mt.dat'
        path = STUDY/'references/data-expanded-ratios'/filename
        if digest(path) != preflight['reference']['curves'][filename]['sha256']:
            raise ValueError('Published reference curve changed')
        target = np.loadtxt(path)
        value,error = estimate['value'][:,indices], estimate['mc_error'][:,indices]
        residual = value-target[:,1:8]
        pulls = residual[:,0]/np.hypot(error[:,0],target[:,8])
        reports[label] = dict(value=value,mc_error=error,reference=target[:,1:8],
                              reference_central_mc_error=target[:,8],residual=residual,
                              central_bin_pulls=pulls,max_abs_central_bin_pull=np.max(np.abs(pulls)),
                              interpretation='per-bin MC compatibility, not a global chi-squared test; reference covariance unavailable')
    nstrict = len(ensembles['S'])
    contrast = independent_jackknife_many(
        ensembles['S']+ensembles['Pi'],lambda *means:
        normalized(sum(means[nstrict:]))-normalized(sum(means[:nstrict])))
    result = dict(created_utc=now(),status='matched ttbar pilot compared; inspect precision and residuals',
                  reference_preflight_sha256=digest(STUDY/'inputs/ttbar_reference_preflight.json'),
                  batches={v:str(p) for v,p in batches.items()},scale_points=order,
                  stage_stream_counts=stage_counts,cross_variant_stage_overlap=0,
                  comparisons=reports,product_minus_strict_normalized=dict(
                      value=contrast['value'][:,indices],mc_error=contrast['mc_error'][:,indices]),
                  limitations='One massless e+mu- assignment, flavour-independent inclusive observable. '
                              'Joint errors are conditional on trained grids. No NNLO coefficient inferred '
                              'from the NNLO-PDF reference curves or from this NLO-PDF Pi-S difference.')
    save(output,clean(result))


def run(args, sources, record, destination):
    from madgraph.various.banner import RunCardNLO
    process = STUDY/'processes'/('TTbar_onshell_epmum_both_'+args.tag)
    if process.exists():
        raise ValueError('Existing ttbar export; inspect before resuming')
    preflight = json.loads((STUDY/'inputs/ttbar_reference_preflight.json').read_text())
    physics = preflight['physics']
    baseline = json.loads(args.reference_manifest.read_text())
    batches = {}
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        command_card = STUDY/'inputs'/(process.name+'.mg5')
        command_card.write_text('\n'.join([
            'set auto_update 0','set automatic_html_opening False','set notification_center False',
            'set nb_core 64','set run_mode 2','set low_mem_multicore_nlo_generation False',
            'set lhapdf '+str(ROOT/'LHAPDF/bin/lhapdf-config'),
            'set fastjet '+str(STUDY/'local/fastjet/bin/fastjet-config'),
            'import model loop_sm-no_b_mass','set decay_bottom_mass 0',
            'define p = g u c d s b u~ c~ d~ s~ b~',
            'generate p p > t t~ QCD=2 QED=0 [QCD], '
            '(t > w+ b QED=1 [QCD], w+ > e+ ve), '
            '(t~ > w- b~ QED=1 [QCD], w- > mu- vm~)',
            'output fNLO '+str(process),'']))
        launch(['taskset','-c','0-63',sys.executable,str(ROOT/'bin/mg5_aMC'),str(command_card)],
               ROOT,STUDY/'logs'/(process.name+'_generation.log'),
               STUDY/'inputs'/(process.name+'_generation.json'),dict(source_hashes=sources,max_cores=64))
        if not (process/'Cards/decay_card.dat').is_file():
            raise ValueError('ttbar export did not complete')
        if digest(process/'FixedOrderAnalysis'/ANALYSIS) != sources['Template/fNLO/FixedOrderAnalysis/'+ANALYSIS]:
            raise ValueError('ttbar export analysis does not match')
        parameters(process,preflight)
        for i,variant in enumerate(('LO','S','Pi')):
            if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
                raise ValueError('Frozen ttbar reference source changed')
            seed = args.seed_start+i
            name = '%s_ref1901_mt_%d' % (variant,seed)
            record.update(status='running',current_run=name)
            save(destination,record)
            archive = process/'study_cards'/name
            archive.mkdir()
            settings = run_settings(baseline,physics,seed,args.accuracy,args.job_seconds)
            card = RunCardNLO(str(process/'Cards/run_card.dat'))
            for key,value in settings.items():
                card.set(key,value,user=True,raiseerror=True)
            card.write(str(process/'Cards/run_card.dat'),template=str(process/'Cards/run_card_default.dat'))
            (process/'Cards/decay_card.dat').write_text(decay_card(variant,physics))
            (process/'Cards/FO_analyse_card.dat').write_text(
                'FO_ANALYSIS_FORMAT=HwU\nFO_ANALYSE='+ANALYSIS.replace('.f90','.o')+'\n'
                'FO_EXTRALIBS=\nFO_EXTRAPATHS=\nFO_INCLUDEPATHS=\n')
            for filename in ('run_card.dat','param_card.dat','decay_card.dat','FO_analyse_card.dat'):
                shutil.copy2(process/'Cards'/filename,archive/filename)
            manifest = dict(variant=variant,settings=settings,physics=physics,flavours=['e','mu'],
                            source_hashes=sources,max_cores=64,PDF=preflight['PDF'],
                            hashes={p.name:digest(p) for p in archive.iterdir()},
                            base_manifest_sha256=digest(args.reference_manifest),
                            reference_preflight_sha256=digest(STUDY/'inputs/ttbar_reference_preflight.json'))
            save(archive/'manifest.json',manifest)
            execution = dict(variant=variant,pilot=True,max_cores=64,source_hashes=sources)
            launch(['taskset','-c','0-63',sys.executable,str(process/'bin/calculate_xsect'),
                    'NLO','-f','-n',name,'--multicore','--nb_core=64'],process,
                   STUDY/'logs'/(process.name+'_'+name+'.log'),archive/'execution.json',execution)
            hwu = process/'Events'/name/'MADatNLO.HwU'
            if not hwu.is_file():
                execution.update(status='failed: no final histogram',missing_output=str(hwu))
                save(archive/'execution.json',execution)
                raise ValueError('ttbar calculator returned without its final histogram')
            row = vectors(hwu.read_bytes())
            check_histogram(row)
            audit_path = STUDY/'results/audits'/(process.name+'_'+name+'.json')
            save(audit_path,dict(created_utc=now(),path=str(hwu),output_sha256=digest(hwu),
                                 variant=variant,manifest=manifest,execution=execution,
                                 status='reference layout/normalization identities passed; no literature agreement certified',
                                 inclusive_pb=float(row['values'][0,0]),inclusive_mc_error_pb=float(row['errors'][0])))
            worker_archive = STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
            harvest(audit_path,worker_archive)
            batch_path = STUDY/'results'/(process.name+'_'+name+'_batches.npz')
            read(worker_archive.with_suffix('.json'),batch_path)
            batches[variant] = batch_path
            job = dict(variant=variant,audit=str(audit_path),batches=str(batch_path),
                       process=str(process),run=name,workers=str(worker_archive))
            if variant != 'LO':
                virtuals = STUDY/'results'/(process.name+'_'+name+'_virtual_checks.json')
                audit_virtuals(worker_archive.with_suffix('.json'),virtuals)
                job['analytic_virtual_checks'] = str(virtuals)
            record['jobs'].append(job)
            save(destination,record)
    compare(batches,preflight,STUDY/'results'/('ttbar_reference_'+args.tag+'.json'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',required=True,type=Path)
    parser.add_argument('--reference-manifest',required=True,type=Path)
    parser.add_argument('--accuracy',type=float,default=.003)
    parser.add_argument('--seed-start',type=int,default=49501)
    parser.add_argument('--job-seconds',type=float,default=.05,
                        help='Short batches for fast ttbar LO channels; the reader still requires >=5 per stratum')
    args = parser.parse_args()
    if (not args.tag.replace('_','').isalnum() or not 0.<args.accuracy<1. or
            not 0<args.seed_start<=900000000-3 or
            not math.isfinite(args.job_seconds) or args.job_seconds<=0.):
        raise ValueError('Invalid tag or accuracy')
    destination = STUDY/'inputs'/('ttbar_reference_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('ttbar queue already exists')
    sources = {p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    for p in (Path(__file__).resolve(),STUDY/'scripts/reference_batches.py',
              STUDY/'scripts/run_inclusive.py',STUDY/'scripts/harvest_splits.py',
              STUDY/'scripts/audit_virtuals.py',STUDY/'scripts/compare_current_batches.py',
              STUDY/'scripts/read_splits.py',
              STUDY/'scripts/replica_statistics.py',STUDY/'scripts/rng_history.py',
              STUDY/'inputs/ttbar_reference_preflight.json',
              ROOT/'Template/fNLO/FixedOrderAnalysis'/ANALYSIS,
              ROOT/'Template/fNLO/Source/ranmar.f90'):
        sources[str(p.relative_to(ROOT))] = digest(p)
    record = dict(created_utc=now(),status='waiting for ttW literature queue',source_hashes=sources,
                  after_queue=str(args.after_queue.resolve()),max_cores=64,jobs=[],
                  seed_start=args.seed_start,job_seconds=args.job_seconds,accuracy=args.accuracy)
    save(destination,record)
    try:
        while True:
            previous = json.loads(args.after_queue.read_text())
            if previous['status'] == 'stopped':
                raise ValueError('Preceding literature queue stopped')
            if previous['status'] == 'reference pilots finished; agreement not automatically certified':
                break
            time.sleep(15)
        if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
            raise ValueError('Frozen ttbar source changed before launch')
        run(args,sources,record,destination)
        record.update(status='ttbar reference pilot finished; inspect comparison and convergence',finished_utc=now())
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)
