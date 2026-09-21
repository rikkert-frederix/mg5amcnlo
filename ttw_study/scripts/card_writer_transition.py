#!/usr/bin/env python3
"""Prove the narrow DECAYMASS writer transition preserves massless cards.

The historical source is reconstructed by reversing exactly one known line,
then checked against the predecessor's frozen hash before it is evaluated.
No old export, card, source record or run is modified.
"""
import argparse
import hashlib
import json
from pathlib import Path

from campaign import ROOT, digest, now, save

SOURCE='models/check_param_card.py'
BEFORE="        elif self.name.startswith('decay'):\n"
AFTER="        elif self.name.startswith('decay_table'):\n"
STATUS='DECAYMASS writer transition verified; massless card serialization unchanged'


def old_source(current,expected_hash):
    if current.count(AFTER)!=1 or BEFORE in current:
        raise ValueError('The current writer is not the specific DECAYMASS header fix')
    previous=current.replace(AFTER,BEFORE,1)
    if hashlib.sha256(previous.encode()).hexdigest()!=expected_hash:
        raise ValueError('Writer change is not exactly the audited one-line correction')
    return previous


def validate(previous,record):
    differences={name:(sha,digest(ROOT/name)) for name,sha in previous['source_hashes'].items()
                 if digest(ROOT/name)!=sha}
    if record['status']!=STATUS or set(differences)!={SOURCE}:
        raise ValueError('Transition is only valid for the single DECAYMASS writer correction')
    if list(differences[SOURCE])!=record['source_transition']:
        raise ValueError('Writer transition source fingerprints differ')
    old_source((ROOT/SOURCE).read_text(),differences[SOURCE][0])
    if any(digest(name)!=sha for name,sha in record['evidence_sha256'].items()):
        raise ValueError('Writer transition evidence changed')
    if len(record['massless_cards'])!=len(previous['jobs']) or not record['massive_roundtrip_passed']:
        raise ValueError('Incomplete writer transition card checks')
    return record


def run(predecessor_path,failed_queue_path,regressions_path,output):
    if output.exists():
        raise ValueError('Refuse to overwrite a writer transition record')
    previous=json.loads(predecessor_path.read_text())
    failed=json.loads(failed_queue_path.read_text())
    if (not previous['status'].startswith('stage-audited S/Pi pilots finished')
            or len(previous['jobs'])!=4 or failed['status']!='stopped' or failed['jobs']):
        raise ValueError('Require completed controls and a stopped pre-integration massive queue')
    process=Path(failed['current_job']['process'])
    if list((process/'study_cards').glob('*/execution.json')):
        raise ValueError('Failed massive export already has an integration record')
    source=(ROOT/SOURCE).read_text()
    before=old_source(source,previous['source_hashes'][SOURCE])
    namespace={'__name__':'ttw_frozen_card_writer','__file__':str(ROOT/SOURCE)}
    exec(compile(before,'<hash-verified frozen card writer>','exec'),namespace)
    from models.check_param_card import ParamCard
    old_card=namespace['ParamCard']
    records=[]
    evidence={str(p.resolve()):digest(p) for p in (predecessor_path,failed_queue_path,regressions_path)}
    if '\nOK\n' not in regressions_path.read_text():
        raise ValueError('Writer/model regressions did not pass')
    for job in previous['jobs']:
        card_path=Path(job['process'])/'study_cards'/job['run']/'param_card.dat'
        old,new=old_card(str(card_path)),ParamCard(str(card_path))
        if 'decaymass' in old or 'decaymass' in new:
            raise ValueError('The preserved control is not massless')
        text=new.write(precision=16)
        if old.write(precision=16)!=text:
            raise ValueError('Massless card serialization changed')
        evidence[str(card_path)]=digest(card_path)
        records.append(dict(path=str(card_path),serialized_sha256=hashlib.sha256(text.encode()).hexdigest()))
    original=process/'study_cards/original_param_card.dat'
    broken=process/'Cards/param_card.dat'
    fixed=ParamCard(ParamCard(str(original)).write(precision=16))
    legacy=old_card(old_card(str(original)).write(precision=16))
    if ('decaymass' in legacy or 'decaymass' in ParamCard(str(broken))
            or fixed['decaymass'].get((5,)).value!=4.8 or fixed['mass'].get((5,)).value!=0.):
        raise ValueError('The archived failure does not reproduce the specific writer bug')
    evidence.update({str(p):digest(p) for p in (original,broken,Path(__file__).resolve())})
    record=dict(created_utc=now(),status=STATUS,source_transition=[previous['source_hashes'][SOURCE],digest(ROOT/SOURCE)],
                evidence_sha256=evidence,massless_cards=records,massive_roundtrip_passed=True,
                limitations='Exact one-line Python writer transition only. No Fortran, physical input or RNG change. '
                            'The failed export remains untouched and has no MC output; fresh massive exports are required.')
    validate(previous,record)
    save(output,record)
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predecessor',required=True,type=Path)
    parser.add_argument('--failed-queue',required=True,type=Path)
    parser.add_argument('--regressions',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    run(args.predecessor,args.failed_queue,args.regressions,args.output)
