#!/usr/bin/env python3
"""Separate matched NWA validation configuration for arXiv:2005.09427.

This does not modify the main benchmark. LO is an auxiliary Born sample
using the NLO PDF, not the paper's LO-PDF prediction. P can be compared
directly to NWA_LOdecay. S and the auxiliary LO reconstruct the paper's
additive NLO numerator with its unexpanded, fixed NLO-width denominator.
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

from campaign import ROOT, STUDY, RUNTIME_SOURCES, benchmark_param_card, digest, generate, now, save
from audit_virtuals import audit as audit_virtuals
from harvest_splits import harvest
from reference_batches import read as read_batches
from run_inclusive import launch

MT,MW,MZ,GF = 172.5,80.385,91.1876,1.166378e-5
G0,GN,GW = 1.48063,1.35355,2.09767
MU = MT+MW/2.
ANALYSIS = 'analysis_HwU_pp_ttxw_reference.f90'
TARGETS_AB = {'plus':{'P':127.,'unexpanded_NWA':123.},
              'minus':{'P':69.8,'unexpanded_NWA':68.}}


def decay_card(variant):
    from madgraph.fks.fks_decay import decay_card_text
    if variant not in ('LO','P','S'):
        raise ValueError('Reference uses only LO auxiliary, P and S')
    return decay_card_text(
        {6:G0,24:GW},{6:MU,24:MW},nlo_width_pdgs={6},nlo_widths={6:GN},
        production_order='LO' if variant == 'LO' else 'NLO',
        decay_order='NLO' if variant == 'S' else 'LO',nlo_decay_combination='ADDITIVE',
        decay_scale_variation_mode='CORRELATED',decay_scale_factors=(1.,.5,2.),
        decay_width_scale_modes={6:'EXPLICIT'},
        lo_width_variations={(6,f):G0 for f in (.5,2.)},
        nlo_width_variations={(6,f):GN for f in (.5,2.)},production_scale_momenta='CORE')


def unexpanded(strict,born):
    """Undo only the additive width counterterm, then use the NLO widths."""
    g = (GN-G0)/G0
    return (np.asarray(strict)+2*g*np.asarray(born))/(1+g)**2


def central_reference_comparison(results,charge):
    # The paper rounds Table 4/5 rates to 0.1 ab and does not list MC errors.
    # Residual/our-error is not a reference-inclusive statistical pull.
    g = (GN-G0)/G0
    strict,born,production = (results[v] for v in ('S','LO','P'))
    value = float(unexpanded(strict['fiducial_pb']['central value'],born['fiducial_pb']['central value']))
    error = math.hypot(strict['fiducial_pb']['dy'],2*g*born['fiducial_pb']['dy'])/(1+g)**2
    rows = {}
    for label,actual,uncertainty in (
            ('P',production['fiducial_pb']['central value'],production['fiducial_pb']['dy']),
            ('unexpanded_NWA',value,error)):
        target = TARGETS_AB[charge][label]
        rows[label] = dict(value_ab=actual*1e6,mc_error_ab=uncertainty*1e6,
                           published_rounded_ab=target,relative_residual=actual*1e6/target-1.,
                           residual_over_our_mc_error=((actual*1e6-target)/(uncertainty*1e6)
                                                       if uncertainty > 0 else None),
                           reference_mc_error='not supplied; not included in the diagnostic denominator')
    return rows


def parameters(process):
    import lhapdf
    from models import import_ufo,model_reader
    from models.check_param_card import ParamCard
    pdf_record = next(row for row in json.loads(
        (STUDY/'inputs/parameter_variations_preflight.json').read_text())['PDFs']
        if row['name'] == 'NNPDF30_nlo_as_0118')
    pdf_directory = ROOT/'LHAPDF/share/LHAPDF/NNPDF30_nlo_as_0118'
    if digest(pdf_directory/'NNPDF30_nlo_as_0118.info') != pdf_record['metadata_sha256']:
        raise ValueError('Reference PDF metadata changed since installation audit')
    for name,checksum in pdf_record['member_sha256'].items():
        if digest(pdf_directory/name) != checksum:
            raise ValueError('Reference PDF member changed: '+name)
    benchmark_param_card(process,json.loads((STUDY/'inputs/benchmark.json').read_text()),'onshell')
    card = ParamCard(str(process/'Cards/param_card.dat'))
    alpha_inverse = 1./(math.sqrt(2)*GF*MW**2/math.pi*(1-MW**2/MZ**2))
    pdf = lhapdf.mkPDF('NNPDF30_nlo_as_0118',0)
    for block,pdg,value in (('sminputs',1,alpha_inverse),('sminputs',2,GF),
                            ('sminputs',3,pdf.alphasQ(MZ)),('decay',6,GN),
                            ('decay',24,GW),('decay',23,2.50775)):
        card[block].get((pdg,)).value = value
    model = model_reader.ModelReader(import_ufo.import_model('loop_sm-no_b_mass'))
    card.update_dependent(model,None,20)
    model.set_parameters_and_couplings(card)
    if abs(model.get_mass(24)-MW)>1.e-10:
        raise ValueError('Reference weak inputs did not reproduce MW')
    card.write(str(process/'Cards/param_card.dat'),precision=16)
    return dict(GF=GF,mt=MT,MW=MW,MZ=MZ,alpha_Gmu_inverse=alpha_inverse,CKM='diagonal',
                top_width_lo=G0,top_width_nlo=GN,W_width=GW,Z_width=2.50775,
                alpha_s_at_mt=pdf.alphasQ(MT),alpha_s_at_mu0=pdf.alphasQ(MU),
                PDF='NNPDF30_nlo_as_0118',PDF_id=260000,PDF_member=0,
                PDF_DataVersion=pdf_record['DataVersion'],
                PDF_metadata_sha256=pdf_record['metadata_sha256'],
                PDF_member0_sha256=pdf_record['member_sha256']['NNPDF30_nlo_as_0118_0000.dat'],
                reference='arXiv:2005.09427, Sec. 3, Tables 4/5 and footnote 2')


def reusable_reference(variant,charge,physics,settings,artifact):
    """Use a completed, unchanged variant without relabelling it as new MC."""
    from inclusive_report import STATISTICAL
    audit_path=Path(artifact['audit'])
    audit=json.loads(audit_path.read_text())
    manifest,execution=audit['manifest'],audit['execution']
    hwu=Path(audit['path'])
    archive=hwu.parents[2]/'study_cards'/hwu.parent.name
    batch=Path(artifact['batches'])
    metadata=json.loads(batch.with_suffix('.json').read_text())
    worker_record=Path(metadata['input'])
    workers=json.loads(worker_record.read_text())
    physical=lambda s:{k:v for k,v in s.items() if k not in STATISTICAL}
    if (audit['variant']!=variant or manifest['charge']!=charge or manifest['physics']!=physics
            or physical(manifest['settings'])!=physical(settings)
            or execution['status']!='finished' or execution['returncode']!=0
            or digest(hwu)!=audit['output_sha256'] or digest(hwu)!=execution['output_sha256']
            or digest(batch)!=metadata['arrays_sha256']
            or digest(worker_record)!=metadata['input_sha256']
            or digest(audit_path)!=workers['audit_sha256']
            or digest(artifact['workers'])!=workers['archive_sha256']
            or not metadata.get('refinement_rng_audit')
            or any(digest(archive/name)!=sha for name,sha in manifest['files'].items())):
        raise ValueError('Reused reference is incomplete, changed or physically mismatched')
    if variant=='S' and not Path(artifact['analytic_virtual_checks']).is_file():
        raise ValueError('Missing analytic virtual audit for reused strict reference')
    return execution,str(archive/'execution.json')


def run_charge(charge,tag,reference_manifest,accuracy,seed,sources,reused=None):
    from madgraph.various.banner import RunCardNLO
    from madgraph.various.histograms import HwUList
    flavours = ['e','mu','e'] if charge == 'plus' else ['mu','e','e']
    generate(argparse.Namespace(charge=charge,w_treatment='onshell',flavours=flavours,
                                corrected='both',decay_bottom_mass=0.,tag=tag))
    process = STUDY/'processes'/('TTW%s_onshell_%s_mb0p0_both_%s' % (charge,''.join(flavours),tag))
    analysis = process/'FixedOrderAnalysis'/ANALYSIS
    if digest(analysis) != sources['Template/fNLO/FixedOrderAnalysis/'+ANALYSIS]:
        raise ValueError('Fresh export lacks the unchanged reference analysis')
    reference = json.loads(reference_manifest.read_text())
    if reference['w_treatment'] != 'onshell':
        raise ValueError('Reference requires an on-shell starting configuration')
    results, artifacts, runs = {}, {}, {}
    reused=reused or {}
    if set(reused)-{'LO','P','S'}:
        raise ValueError('Unknown reused reference prescription')
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        physics = parameters(process)
        save(process/'study_cards/reference_2005_parameters.json',physics)
        for index,variant in enumerate(('LO','P','S')):
            if any(digest(ROOT/path)!=checksum for path,checksum in sources.items()):
                raise ValueError('Reference source changed; prepare a fresh queue')
            settings = dict(reference['settings'])
            settings.update(lhaid=260000,iseed=seed+index,req_acc_fo=accuracy,fo_job_target_time=1.,
                            fixed_ren_scale=True,fixed_fac_scale=True,fixed_qes_scale=True,
                            mur_ref_fixed=MU,muf_ref_fixed=MU,qes_ref_fixed=MU,
                            reweight_pdf=False,reweight_scale=True)
            if variant in reused:
                results[variant],runs[variant]=reusable_reference(variant,charge,physics,settings,reused[variant])
                artifacts[variant]=dict(reused[variant],reused_completed_variant=True)
                continue
            name = '%s_ref2005_fixed_%d' % (variant,seed+index)
            archive = process/'study_cards'/name
            archive.mkdir()
            card = RunCardNLO(str(process/'Cards/run_card.dat'))
            for key,value in settings.items():
                card.set(key,value,user=True,raiseerror=True)
            card.write(str(process/'Cards/run_card.dat'),template=str(process/'Cards/run_card_default.dat'))
            (process/'Cards/decay_card.dat').write_text(decay_card(variant))
            (process/'Cards/FO_analyse_card.dat').write_text(
                'FO_ANALYSIS_FORMAT=HwU\nFO_ANALYSE='+ANALYSIS.replace('.f90','.o')+'\n'
                'FO_EXTRALIBS=\nFO_EXTRAPATHS=\nFO_INCLUDEPATHS=\n')
            for filename in ('run_card.dat','param_card.dat','decay_card.dat','FO_analyse_card.dat'):
                shutil.copy2(process/'Cards'/filename,archive/filename)
            manifest = dict(variant=variant,settings=settings,physics=physics,flavours=flavours,charge=charge,
                            source_hashes=sources,files={p.name:digest(p) for p in archive.iterdir()},
                            base_manifest=str(reference_manifest),base_manifest_sha256=digest(reference_manifest),
                            prescription='strict additive; constant explicit LO/NLO widths under all scales',
                            production_scale_GeV=MU,decay_numerator_scale_GeV=MU,
                            scale_variation='production muR, muF; both decay numerators follow production muR',
                            auxiliary_LO_uses_NLO_PDF=True,max_cores=64)
            save(archive/'manifest.json',manifest)
            execution = dict(variant=variant,charge=charge,pilot=True,max_cores=64,source_hashes=sources)
            launch(['taskset','-c','0-63',sys.executable,str(process/'bin/calculate_xsect'),
                    'NLO','-f','-n',name,'--multicore','--nb_core=64'],process,
                   STUDY/'logs'/(process.name+'_'+name+'.log'),archive/'execution.json',execution)
            hwu = process/'Events'/name/'MADatNLO.HwU'
            if not hwu.is_file():
                execution.update(status='failed: no final histogram',missing_output=str(hwu))
                save(archive/'execution.json',execution)
                raise ValueError('Reference calculator returned without its final histogram; inspect worker log')
            histograms = HwUList(str(hwu),raw_labels=True)
            if len(histograms) != 2 or any(len(h.bins)!=2 for h in histograms):
                raise ValueError('Unexpected reference histogram layout')
            slot = 0 if charge == 'plus' else 1
            if any(value != 0. for row in histograms[1-slot].bins for value in row.wgts.values()):
                raise ValueError('Reference histogram has the wrong charge')
            if any(not math.isfinite(value) for h in histograms for row in h.bins for value in row.wgts.values()):
                raise ValueError('Nonfinite reference weight')
            execution.update(output_sha256=digest(hwu),input_pb=dict(histograms[slot].bins[0].wgts),
                             fiducial_pb=dict(histograms[slot].bins[1].wgts))
            save(archive/'execution.json',execution)
            results[variant] = execution
            runs[variant] = str(archive/'execution.json')
            audit_path = STUDY/'results/audits'/(process.name+'_'+name+'.json')
            save(audit_path,dict(created_utc=now(),path=str(hwu),output_sha256=digest(hwu),
                variant=variant,manifest=manifest,execution=execution,
                status='reference layout and finite weights passed; literature agreement not certified'))
            workers = STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
            harvest(audit_path,workers)
            batches = STUDY/'results'/(process.name+'_'+name+'_batches.npz')
            read_batches(workers.with_suffix('.json'),batches)
            artifacts[variant] = dict(audit=str(audit_path),workers=str(workers),batches=str(batches))
            if variant == 'S':
                virtuals = STUDY/'results'/(process.name+'_'+name+'_virtual_checks.json')
                audit_virtuals(workers.with_suffix('.json'),virtuals)
                artifacts[variant]['analytic_virtual_checks'] = str(virtuals)
    report = dict(created_utc=now(),status='matched reference pilot; inspect deviations before validation',
                  process=str(process),charge=charge,physics=physics,comparisons=central_reference_comparison(results,charge),
                  artifacts=artifacts,
                  runs=runs,
                  unexpanded_conversion='[S + 2*(Gamma_NLO/Gamma_LO-1)*LO]/(Gamma_NLO/Gamma_LO)^2',
                  note='Auxiliary LO uses NLO PDF. No comparison to LO-PDF table entries. '
                       'Published MC errors are unavailable, so residual/our-error is not a combined pull.')
    save(STUDY/'results'/('reference_2005_'+charge+'_'+tag+'.json'),report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',type=Path)
    parser.add_argument('--reference-manifest',required=True,type=Path)
    parser.add_argument('--accuracy',type=float,default=.01)
    parser.add_argument('--seed-start',type=int,default=47001)
    parser.add_argument('--reuse-artifacts',type=Path,
                        help='Explicit audited per-charge/variant artifacts to retain without rerunning')
    args = parser.parse_args()
    if (not args.tag.replace('_','').isalnum() or not 0.<args.accuracy<1.
            or not 0 < args.seed_start < 900000000-102):
        raise ValueError('Invalid tag or accuracy')
    destination = STUDY/'inputs'/('literature_reference_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Existing reference queue; inspect before resuming')
    sources = {path:digest(ROOT/path) for path in RUNTIME_SOURCES}
    for path in (Path(__file__).resolve(),STUDY/'scripts/run_inclusive.py',
                 STUDY/'scripts/inclusive_report.py',
                 STUDY/'scripts/harvest_splits.py',STUDY/'scripts/reference_batches.py',
                 STUDY/'scripts/rng_history.py',STUDY/'scripts/audit_virtuals.py',
                 ROOT/'Template/fNLO/FixedOrderAnalysis'/ANALYSIS):
        sources[str(path.relative_to(ROOT))] = digest(path)
    reuse={}
    if args.reuse_artifacts:
        reuse=json.loads(args.reuse_artifacts.read_text())
        if set(reuse)-{'plus','minus'}:
            raise ValueError('Unknown charge in reference reuse record')
        sources[str(args.reuse_artifacts.resolve().relative_to(ROOT))]=digest(args.reuse_artifacts)
    record = dict(created_utc=now(),status='prepared',source_hashes=sources,max_cores=64,
                  seed_start=args.seed_start,jobs=[])
    save(destination,record)
    try:
        if args.after_queue:
            record.update(status='waiting for validation pilots',after_queue=str(args.after_queue.resolve()))
            save(destination,record)
            while True:
                previous = json.loads(args.after_queue.read_text())
                if previous['status'] == 'stopped':
                    raise ValueError('Preceding validation queue stopped')
                if previous['status'] == 'validation pilots finished; scientific validation still required':
                    break
                time.sleep(15)
        for i,charge in enumerate(('plus','minus')):
            if any(digest(ROOT/path)!=checksum for path,checksum in sources.items()):
                raise ValueError('Frozen reference source changed before launch')
            record.update(status='running',current_charge=charge)
            save(destination,record)
            record['jobs'].append(run_charge(charge,args.tag,args.reference_manifest,args.accuracy,
                                               args.seed_start+100*i,sources,reuse.get(charge)))
            save(destination,record)
        record.update(status='reference pilots finished; agreement not automatically certified',finished_utc=now())
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)
