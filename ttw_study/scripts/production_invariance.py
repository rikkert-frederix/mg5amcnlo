#!/usr/bin/env python3
"""Read-only structural comparison of massless and massive-decay exports.

Compare the actual production density/virtual source and production-owned
FKS records. Decay code and decay-owned singular regions may change. This
is a source-identity check, not a substitute for pole/soft tests, numerical
amplitude tests, or a small-mass continuity study.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

from campaign import STUDY, digest, now, save
from inclusive_report import production_parameters


def fixed_statements(source):
    """Join generated fixed-form continuations, without altering statements."""
    rows=[]
    for line in source.splitlines():
        if not line.strip() or line[0] in 'cC*!':
            continue
        if len(line)<6:
            raise ValueError('Unexpected short fixed-form statement')
        if line[5] not in ' 0':
            if not rows:
                raise ValueError('Orphan fixed-form continuation')
            rows[-1]+=line[6:]
        else:
            rows.append(line[6:])
    return rows


def production_owner(source):
    rows=[line.split() for line in source.splitlines()]
    if ['FORMAT','3'] not in rows:
        raise ValueError('Require version-3 contribution ownership')
    owners=[row for row in rows if row[:1]==['CONTRIBUTION'] and row[2]=='PRODUCTION']
    if len(owners)!=1 or len(owners[0])!=10:
        raise ValueError('Require exactly one complete production contribution')
    row=owners[0]
    identifier,first,last,representative,virtual,parent,occurrence,node=map(int,[row[1],*row[3:]])
    if (identifier!=1 or first!=1 or last<first or not first<=representative<=last
            or virtual!=1 or (parent,occurrence,node)!=(0,0,0)):
        raise ValueError('Unexpected production FKS owner')
    grids=[row for row in rows if row[:2]==['VIRTUAL_GRID',str(identifier)]]
    if not grids:
        raise ValueError('Production virtual grids missing')
    return dict(row=owners[0],virtual_grids=grids,first=first,last=last)


def production_fks(source,first,last):
    """Select all and only production DATA records, retaining their indices."""
    result={}
    globals_={'FKS_I_D','FKS_J_D','NEED_COLOR_LINKS_D'}
    indexed={'FKS_J_FROM_I_D','PARTICLE_TYPE_D','PDG_TYPE_D'}
    seen={region:set() for region in range(first,last+1)}
    for statement in fixed_statements(source):
        compact=re.sub(r'\s+','',statement).upper()
        if not compact.startswith('DATA'):
            continue
        match=re.fullmatch(r'DATA([^/]+)/([^/]+)/',compact)
        if not match:
            raise ValueError('Unrecognized FKS DATA statement: '+statement)
        key,value=match.groups()
        if key in globals_:
            values=value.split(',')
            if len(values)<last or key in result:
                raise ValueError('Incomplete or duplicated global FKS data')
            result[key]=values[first-1:last]
            continue
        match=re.match(r'\(?([A-Z_]+)\((\d+),',key)
        if not match or match[1] not in indexed:
            raise ValueError('Unknown indexed FKS data: '+key)
        region=int(match[2])
        if first<=region<=last:
            if key in result:
                raise ValueError('Duplicate production FKS data: '+key)
            result[key]=value
            seen[region].add(match[1])
    if not globals_.issubset(result) or any(names!=indexed for names in seen.values()):
        raise ValueError('Production FKS data incomplete')
    return result


def production_files(directory):
    files=list(directory.glob('spin_density_production_*.f'))
    required={'spin_density_production_born.f','spin_density_production_real_1.f',
              'spin_density_production_born_contribution_1_link_1.f',
              'spin_density_production_virtual_contribution_1.f'}
    if not required.issubset({p.name for p in files}):
        raise ValueError('Incomplete production density providers')
    virtual=directory/'VContribution1'
    local=[p for p in virtual.iterdir()
           if p.is_file() and not p.is_symlink() and p.suffix in ('.f','.inc')]
    if not {'born_matrix.f','loop_matrix.f','loop_num.f','helas_calls_uvct_1.f'}.issubset(
            {p.name for p in local}):
        raise ValueError('Incomplete generated production virtual provider')
    flavours=directory/'born_leshouche.inc'
    if not flavours.is_file():
        raise ValueError('Production Born flavour inventory missing')
    # Symlinked coupl.inc also declares private decay parameters, so it is
    # deliberately not a byte-identity target. Coupling/model invariance has
    # separate model-reader regression tests; retain this limitation below.
    return {str(p.relative_to(directory)):p for p in files+local+[flavours]}


def core_key(source):
    """Match ordered physical production legs, not decay-model aliases."""
    rows=[line.split() for line in source.splitlines()]
    if ['FORMAT','4'] not in rows:
        raise ValueError('Require version-4 decay/core metadata')
    born=[row for row in rows if row[:1]==['CONTEXT'] and row[2]=='BORN']
    if len(born)!=1 or len(born[0])!=6:
        raise ValueError('Require one complete Born core context')
    identifier=born[0][1]
    count=int(born[0][4])
    legs=[row for row in rows if row[:2]==['CORE_LEG',identifier]]
    if (len(legs)!=count or any(len(row)!=5 for row in legs)
            or [int(row[2]) for row in legs]!=list(range(1,count+1))
            or [row[4] for row in legs]!=['I','I']+['F']*(count-2)):
        raise ValueError('Malformed or incomplete ordered production core')
    return tuple((int(row[3]),row[4]) for row in legs)


def subprocess_inventory(process):
    result={}
    for path in sorted((process/'SubProcesses').glob('P*')):
        if not path.is_dir():
            continue
        key=core_key((path/'decay_chain_info.dat').read_text())
        if key in result:
            raise ValueError('Ambiguous repeated physical production core')
        result[key]=path
    if not result:
        raise ValueError('Missing production subprocesses')
    return result


def kinematic_helpers(source):
    """Separate only the two generic reduction helpers from generated amplitudes."""
    pattern=r'(?ims)^      SUBROUTINE (\w*BUILD_KINEMATIC_MATRIX)\(.*?^      END[ \t]*$'
    matches=list(re.finditer(pattern,source))
    if (len(matches)!=2 or matches[1].group(1).upper()!=
            matches[0].group(1).upper().replace('BUILD_KINEMATIC_MATRIX','MP_BUILD_KINEMATIC_MATRIX')):
        raise ValueError('Missing unique DP/QP kinematic-matrix helper pair')
    bodies='\n'.join(match.group(0) for match in matches)
    return re.sub(pattern,'',source),bodies


def compare_subprocess(reference,candidate,*,allow_kinematic_update=False):
    left,right=production_files(reference),production_files(candidate)
    if set(left)!=set(right):
        raise ValueError('Production source inventory differs')
    hashes={}
    transitions={}
    for name in sorted(left):
        a,b=digest(left[name]),digest(right[name])
        if a!=b:
            if not allow_kinematic_update or name!='VContribution1/CT_interface.f':
                raise ValueError('Production source changed: '+str(right[name]))
            old_code,old_helpers=kinematic_helpers(left[name].read_text())
            new_code,new_helpers=kinematic_helpers(right[name].read_text())
            if old_code!=new_code:
                raise ValueError('Production source changed outside kinematic helpers: '+str(right[name]))
            transitions[name]=dict(reference_sha256=a,candidate_sha256=b,
                unchanged_source_sha256=hashlib.sha256(old_code.encode()).hexdigest(),
                reference_helpers_sha256=hashlib.sha256(old_helpers.encode()).hexdigest(),
                candidate_helpers_sha256=hashlib.sha256(new_helpers.encode()).hexdigest())
            continue
        hashes[name]=a
    owners=[]
    regions=[]
    for directory in (reference,candidate):
        owner=production_owner((directory/'nlo_contribution_info.dat').read_text())
        owners.append(owner)
        regions.append(production_fks((directory/'fks_info.inc').read_text(),owner['first'],owner['last']))
    if owners[0]!=owners[1] or regions[0]!=regions[1]:
        raise ValueError('Production FKS ownership, grids or region data changed')
    result=dict(production_sources_sha256=hashes,production_ownership=owners[0],
                production_fks_data=regions[0],
                full_fks_source_hashes={str(p/'fks_info.inc'):digest(p/'fks_info.inc')
                                        for p in (reference,candidate)})
    if transitions:
        result['kinematic_matrix_transition']=transitions
    return result


def model_functions(path):
    source=path.read_text()
    marker='      DOUBLE PRECISION FUNCTION GET_DECAY_MASS_FROM_ID(ID)'
    if source.count(marker)!=1:
        raise ValueError('Missing unique decay-only mass getter')
    production,decay=source.split(marker)
    if 'DC_' in production.upper():
        raise ValueError('Production mass/width getter refers to private decay model')
    return production,decay


def run(reference,candidate,output,*,allow_kinematic_update=False):
    reference,candidate=reference.resolve(),candidate.resolve()
    if output.exists() or reference==candidate:
        raise ValueError('Require distinct exports and a new audit destination')
    generations=[]
    cards=[]
    for process in (reference,candidate):
        generation=STUDY/'inputs'/(process.name+'_generation.json')
        metadata=json.loads(generation.read_text())
        if metadata['status']!='finished':
            raise ValueError('Export generation did not finish')
        generations.append(metadata['args'])
        cards.append(process/'Cards/param_card.dat')
    if (generations[0]['decay_bottom_mass']!=0. or generations[1]['decay_bottom_mass']<=0.
            or any(generations[0][key]!=generations[1][key]
                   for key in ('charge','w_treatment','flavours','corrected'))):
        raise ValueError('Require matched zero/positive-mass production channels')
    scheme_path=candidate/'Cards/decay_mass_scheme.json'
    scheme=json.loads(scheme_path.read_text())
    if (scheme.get('production_model')!='loop_sm-no_b_mass'
            or scheme.get('decay_model')!='loop_sm' or scheme.get('alpha_s_flavours')!=5
            or scheme.get('parameter')!=['decaymass',5]
            or scheme.get('bottom_mass_scheme')!='on-shell'):
        raise ValueError('Not the declared 5FS-production/massive-decay scheme')
    from models.check_param_card import ParamCard
    parsed=[ParamCard(str(p)) for p in cards]
    if (any(float(card['mass'].get((5,)).value)!=0. for card in parsed)
            or 'decaymass' in parsed[0]
            or 'decaymass' not in parsed[1]
            or float(parsed[1]['decaymass'].get((5,)).value)!=generations[1]['decay_bottom_mass']):
        raise ValueError('Actual production/decay bottom masses do not match generation')
    if (production_parameters(cards[0])!=production_parameters(cards[1])
            or parsed[0]['decay'].get((24,)).value!=parsed[1]['decay'].get((24,)).value):
        raise ValueError('Production model inputs, including the internal W width, differ')
    getters=[model_functions(p/'Source/MODEL/get_mass_width_fcts.f') for p in (reference,candidate)]
    if getters[0][0]!=getters[1][0] or 'DC_MDL_MB' not in getters[1][1]:
        raise ValueError('Production mass getter changed or massive decay getter absent')
    directories=[subprocess_inventory(process) for process in (reference,candidate)]
    if not directories[0] or directories[0].keys()!=directories[1].keys():
        raise ValueError('Production subprocess inventory differs')
    records={str(key):dict(reference_directory=str(directories[0][key]),
                          candidate_directory=str(directories[1][key]),
                          **compare_subprocess(directories[0][key],directories[1][key],
                                               allow_kinematic_update=allow_kinematic_update))
             for key in sorted(directories[0])}
    report=dict(created_utc=now(),status='production source and FKS-region identity passed',
                reference=str(reference),candidate=str(candidate),scheme=scheme,
                scheme_sha256=digest(scheme_path),param_card_sha256={str(p):digest(p) for p in cards},
                production_getter_sha256=hashlib.sha256(getters[0][0].encode()).hexdigest(),
                subprocesses=records,script_sha256=digest(Path(__file__)),
                limitations='Byte-identical production providers and matching numeric inputs/production FKS data. '
                            'Common coupling declarations may include additional private decay parameters. '
                            'This does not replace model-reader, numerical-amplitude, massive pole/soft, '
                            'or small-mass continuity checks. Top total widths are intentionally allowed to differ.')
    if any('kinematic_matrix_transition' in row for row in records.values()):
        report['status']='production amplitude and FKS-region identity passed with kinematic-helper transition'
        report['limitations']=('The two generic DP/QP kinematic-matrix helpers changed; both versions and '
            'the unchanged remainder are hashed explicitly. Every production amplitude, loop routing, '
            'coupling input and FKS-region check remains required. This is not byte identity of the '
            'entire numerical backend or a numerical convergence certificate. Top widths may differ.')
    save(output,report)
    print('Production source/FKS identity passed for',len(records),'subprocesses:',output,flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',required=True,type=Path)
    parser.add_argument('--candidate',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--allow-kinematic-update',action='store_true',
        help='Audit a change confined to the two generic reduction-matrix helpers; all amplitudes must match')
    args=parser.parse_args()
    run(args.reference,args.candidate,args.output,allow_kinematic_update=args.allow_kinematic_update)
