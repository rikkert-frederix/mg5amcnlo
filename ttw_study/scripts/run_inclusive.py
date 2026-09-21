#!/usr/bin/env python3
"""Stable-ttW LO/NLO normalization reference with the frozen physical inputs.

Run only when the main queue is idle: the common lock prevents overlapping
MG5 launches. No decayed-event acceptance or branching rescaling is applied.
The native HT/2 scale equals W-system CORE HT/2 for an explicit stable W.
"""
import argparse
import fcntl
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

from campaign import (ROOT, STUDY, RUNTIME_SOURCES, benchmark_param_card,
                      digest, lhapdf_environment, now, output, save)


def launch(command, cwd, log, record_path, record):
    start = time.monotonic()
    environment=lhapdf_environment()
    record.update(started_utc=now(),command=command,log=str(log),status='running')
    record['environment']={'LHAPDF_DATA_PATH':environment['LHAPDF_DATA_PATH']}
    with log.open('w') as stream:
        child = subprocess.Popen(command,cwd=cwd,stdout=stream,stderr=subprocess.STDOUT,
                                 env=environment)
        record['pid'] = child.pid
        save(record_path,record)
        code = child.wait()
    record.update(status='finished' if code == 0 else 'failed',returncode=code,
                  finished_utc=now(),elapsed_s=time.monotonic()-start)
    save(record_path,record)
    if code:
        raise RuntimeError('Failed command; see '+str(log))


def run(args):
    from madgraph.various.banner import RunCardNLO
    from models.check_param_card import ParamCard
    from madgraph.various.histograms import HwUList

    if not args.tag.replace('_','').isalnum():
        raise ValueError('Use an alphanumeric/underscore export tag')
    if not 0. < args.accuracy < 1.:
        raise ValueError('Accuracy must be positive and below one')
    process = STUDY/'processes'/('TTW%s_stable_%s' % (args.charge,args.tag))
    if process.exists():
        raise ValueError('Existing reference export: inspect before resuming')
    benchmark = json.loads((STUDY/'inputs/benchmark.json').read_text())
    reference = json.loads(Path(args.reference_manifest).read_text())
    if reference['w_treatment'] != 'onshell' or reference['production_scale'] != 'core-w-ht-half':
        raise ValueError('Stable-W reference currently matches only the on-shell W-system scale')
    sources = {name:digest(ROOT/name) for name in RUNTIME_SOURCES}
    sources['Template/fNLO/FixedOrderAnalysis/analysis_HwU_template.f90'] = digest(
        ROOT/'Template/fNLO/FixedOrderAnalysis/analysis_HwU_template.f90')
    sources['ttw_study/scripts/run_inclusive.py'] = digest(Path(__file__))
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        command_path = STUDY/'inputs'/(process.name+'.mg5')
        command_path.write_text('\n'.join([
            'set auto_update 0','set automatic_html_opening False','set notification_center False',
            'set nb_core 64','set run_mode 2','set low_mem_multicore_nlo_generation False',
            'set lhapdf '+str(ROOT/'LHAPDF/bin/lhapdf-config'),
            'set fastjet '+str(STUDY/'local/fastjet/bin/fastjet-config'),
            'import model loop_sm-no_b_mass','define p = g u c d s b u~ c~ d~ s~ b~',
            'generate p p > t t~ w%s QCD=2 QED=1 [QCD]' % ('+' if args.charge == 'plus' else '-'),
            'output fNLO '+str(process),'']))
        record = dict(charge=args.charge,source_hashes=sources,max_cores=64,
                      git_head=output(['git','rev-parse','HEAD']))
        launch(['taskset','-c','0-63',sys.executable,str(ROOT/'bin/mg5_aMC'),str(command_path)],
               ROOT,STUDY/'logs'/(process.name+'_generation.log'),
               STUDY/'inputs'/(process.name+'_generation.json'),record)
        if not (process/'Cards/run_card.dat').exists() or (process/'Cards/decay_card.dat').exists():
            raise RuntimeError('Missing stable-production export, or unexpected decay card')
        benchmark_param_card(process,benchmark,'onshell')
        card = ParamCard(str(process/'Cards/param_card.dat'))
        # External on-shell t, tbar and W use zero amplitude widths, exactly
        # as the forced connectors in the normalized decay-chain sample.
        for pdg in (6,24):
            card['decay'].get((pdg,)).value = 0.
        card.write(str(process/'Cards/param_card.dat'),precision=16)
        (process/'Cards/FO_analyse_card.dat').write_text(
            'FO_ANALYSIS_FORMAT=HwU\nFO_ANALYSE=analysis_HwU_template.o\n'
            'FO_EXTRALIBS=\nFO_EXTRAPATHS=\nFO_INCLUDEPATHS=\n')
        for index,variant in enumerate(('LO','P')):
            seed = args.seed+index
            name = '%s_stable_core_ht_half_%d' % (variant,seed)
            archive = process/'study_cards'/name
            archive.mkdir(parents=True)
            settings = dict(reference['settings'])
            settings.update(iseed=seed,req_acc_fo=args.accuracy,fo_job_target_time=1.)
            run_card = RunCardNLO(str(process/'Cards/run_card.dat'))
            for key,value in settings.items():
                run_card.set(key,value,user=True,raiseerror=True)
            run_card.write(str(process/'Cards/run_card.dat'),
                           template=str(process/'Cards/run_card_default.dat'))
            for filename in ('run_card.dat','param_card.dat','FO_analyse_card.dat'):
                shutil.copy2(process/'Cards'/filename,archive/filename)
            manifest = dict(variant=variant,settings=settings,source_hashes=sources,
                            benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),
                            reference_manifest=str(Path(args.reference_manifest).resolve()),
                            reference_manifest_sha256=digest(args.reference_manifest),
                            files={path.name:digest(path) for path in archive.iterdir()},
                            scale_definition='stable t,tbar,W and production radiation mT sum / 2',
                            no_decay_branching_factor_applied=True,
                            external_amplitude_widths_zero=[6,24],max_cores=64)
            save(archive/'manifest.json',manifest)
            execution = dict(variant=variant,charge=args.charge,max_cores=64,
                             source_hashes=sources,pilot=True)
            launch(['taskset','-c','0-63',sys.executable,str(process/'bin/calculate_xsect'),
                    'LO' if variant == 'LO' else 'NLO','-f','-n',name,'--multicore','--nb_core=64'],
                   process,STUDY/'logs'/(process.name+'_'+name+'.log'),archive/'execution.json',execution)
            hwu_path = process/'Events'/name/'MADatNLO.HwU'
            if not hwu_path.exists():
                execution['status'] = 'failed: no histogram output'
                save(archive/'execution.json',execution)
                raise RuntimeError('Stable reference produced no final HwU')
            hwu = HwUList(str(hwu_path),raw_labels=True)
            total = next(hist for hist in hwu if hist.title.strip() == 'total rate')
            execution.update(output_sha256=digest(hwu_path),
                             total_pb=dict(total.bins[0].wgts),
                             branching_e=benchmark['W_width']['branching_e'])
            save(archive/'execution.json',execution)
            print(variant,json.dumps(execution['total_pb']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--charge',choices=('plus','minus'),required=True)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--reference-manifest',required=True)
    parser.add_argument('--accuracy',type=float,default=.02)
    parser.add_argument('--seed',type=int,default=41001)
    run(parser.parse_args())
