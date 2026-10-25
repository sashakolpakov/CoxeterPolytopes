#!/usr/bin/env python3
"""SOURCE ONLY: remaining exact CURRENT universal-receiver joins.

Named batches reuse the frozen literal engine. Accepted complete native-profile,
Schur, whole-factor and physical-star propositions remain explicit premises.
No batch proves those premises, S8, protected74 or a global finite57 theorem.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import types

ROOT=Path(__file__).resolve().parents[1]
LOADER='scripts/check_finite57_current_adaptive9_path_binding_20261025.py'
LOADER_SHA='76022d95fcc5ec6e61bf1ccaea65635adda1ab8593d0bf921d75754f14bd139c'
LEDGER='results/finite57_current_universal_receiver_retention_source_ledger_20261025.json'
LEDGER_SHA='68359555db44cbe8d20a8b2d0110debd2ec9ae93d94fb3e2f9af58826422b797'
D11='fresh_affine_d11_leaf_pair_transfer_20261017'
REPAIR='fresh_d12_universal_ideal_output_repairs_20261018'
U='fresh_affine_d12_leaf_pair_transfer_20261017'
LARGE='fresh_large_finite_d_leaf_half_transfer_20261017'
BC12='fresh_order12_bc_near_terminal_half_transfer_20261018'
BC13='fresh_order13_bc_inner_terminal_transfer_20261018'
D13='fresh_order13_d_center_leaf_half_transfer_20261018'
D14='fresh_order14_d_near_fork_half_transfer_20261018'
D15='fresh_order15_bd_basic_fractional_transfer_20261018'
BRIDGE='results/finite57_current_universal_receiver_retention_bridge_20261025.md'
WORKER='scripts/check_finite57_current_universal_receivers_20261025.py'
G4='scripts/check_finite57_g4_lower_cycle_terminal_bindings_20261025.py'
FINITE='scripts/check_finite57_terminal_collectors_syntax_safe_20261025.py'
EXTRA_PINS={
 f'scripts/{D11}.py':'b1297fe9e36bf3795e32cd5b294e0043c05cba97ff2dccbfc868a194180a90f8',
 f'results/{D11}.md':'c91b3c73c66507508dae9ff4fb3bbebdf3fea1dd22eaa72c9992cd1fc5a7d351',
 f'results/{D11}.json':'6ce315f8e0827db3415d9c3ce21cba25b2e13c7ea42dd6f3b51b3b32ff8c61da',
 f'scripts/{REPAIR}.py':'f4591d4897d321b8eba46093ab128c79e89942739af3c17af5eef22023571b83',
 f'results/{REPAIR}.md':'231480a5fd21be324041d7e08cfc6a7577a85513e438a12d36656b15a93c2c12',
 f'results/{REPAIR}.json':'82850cfdfe1ff1a9d56f4f76f5466d02d7539d74670e4172aa3f7c7643a0916d',
 BRIDGE:'309bbd85a766466bcb3bb4437d7a5d50ffe8c614c4892f40f4ccc06d842893eb',
 'results/finite57_current_universal_receiver_bridge_root_acceptance_20261025.json':'6bc292fb08a50c15cfa1cff2ef9e43def44d4febd942904f0cdd1033eac4838e',
 'results/review_finite57_affine_order12_d_current_receivers_source_20261025.md':'603790d2a06379b86396d499c3166a04a90824437fcec89fbca247851beaf1f5',
 'results/finite57_affine_order12_d_current_receiver_root_acceptance_20261025.json':'57e5d59ee226c494310a36efb5bcdc17a8e18b203ebf3892e787861b452d4a7d',
 G4:'f17014594587a68ef766ae245db372550fe89adc02bedf9e44c477a8d79eb242',
 FINITE:'930ffff5601e6096359d3b2f8846fa6d7d050f95fb2aabe449447e0d3d0bf9ca',
}
MAX_OUTPUT=500*1024
BATCHES=('order12-d','five-local','final-quiet')+tuple(f'large-d-{i}' for i in range(4))
ORDER12_MARKS=tuple(product((8,9),(10,11)))
ORDER13_MARKS=tuple(product((9,10),(11,12)))
THEOREM_BATCHES={
 'finite57.current.affine_order12_D_opposite_leaves.whole_face.v1':[(WORKER,'order12-d')],
 'finite57.current.BC12_near_terminal.whole_face.v1':[(WORKER,'five-local'),(G4,'order12-near')],
 'finite57.current.C12_near_companion_zero.whole_face.v1':[(WORKER,'five-local')],
 'finite57.current.BC13_inner_terminal.whole_face.v1':[(WORKER,'five-local'),(FINITE,'universal-bc13')],
 'finite57.current.C13_inner_first_zero.whole_face.v1':[(WORKER,'five-local')],
 'finite57.current.D13_center_leaf.whole_face.v1':[(WORKER,'five-local'),(FINITE,'universal-d13')],
 'finite57.current.D14_near_fork.whole_face.v1':[(WORKER,'five-local'),(FINITE,'universal-d14-near')],
 'finite57.current.D15_basic_fractional.whole_face.v1':[(WORKER,'five-local'),(FINITE,'universal-d15-basic')],
 'finite57.current.largeD_leaf9_distinct_half.whole_face.v1':[(WORKER,f'large-d-{i}') for i in range(4)],
 'finite57.current.final_guard_quiet_receiver.whole_face.v1':[(WORKER,'final-quiet'),
     ('scripts/check_finite57_quiet_b13_current_priority_bridge_20261025.py','quiet-b13')],
}
PRIORITIES={
 ('C',10,0,9):'finite57.current.C10_endpoint_zero.whole_face.v1',
 **{('B',11,8,j):'finite57.current.BC11_zero.whole_face.v1' for j in (9,10)},
 **{('C',11,*p):'finite57.current.BC11_zero.whole_face.v1' for p in ((0,9),(0,10),(1,10))},
 ('B',12,0,9):'finite57.current.B12_center_terminal.whole_face.v1',
 **{('B',12,8,j):'finite57.current.BC12_near_terminal.whole_face.v1' for j in (10,11)},
 ('C',12,1,10):'finite57.current.C12_near_companion_zero.whole_face.v1',
 **{('D',12,i,j):'finite57.current.affine_order12_D_opposite_leaves.whole_face.v1' for i,j in ORDER12_MARKS},
 **{('B',13,8,j):'finite57.current.BC13_inner_terminal.whole_face.v1' for j in (11,12)},
 **{('D',13,*p):'finite57.current.D13_center_leaf.whole_face.v1' for p in ((0,11),(0,12),(8,9),(8,10))},
 **{('D',14,*p):'finite57.current.D14_near_fork.whole_face.v1' for p in ((1,12),(1,13),(8,10),(8,11))},
 **{('D',15,*p):'finite57.current.D15_basic_fractional.whole_face.v1' for p in ((2,13),(2,14),(8,11),(8,12))},
}


def need(ok,message):
    if not ok:raise ValueError(message)


def raw(name):
    p=(ROOT/name).resolve()
    need(p.is_relative_to(ROOT) and p.is_file(),'bounded repository source')
    with p.open('rb') as stream:data=stream.read(4*1024**2+1)
    need(0<len(data)<=4*1024**2,'bounded source read')
    return data


def prepare():
    encoded=raw(LOADER);need(sha256(encoded).hexdigest()==LOADER_SHA,'exact shared current loader')
    a=types.ModuleType('finite57_universal_current_loader');a.__file__=str(ROOT/LOADER)
    exec(compile(encoded,a.__file__,'exec'),a.__dict__)
    r,c,extra=a.prepare()  # Enforces the already-held shared proof lock.
    ledger=raw(LEDGER);need(sha256(ledger).hexdigest()==LEDGER_SHA,'accepted universal source ledger')
    entries=json.loads(ledger,object_pairs_hook=c.unique)['sources']
    pins={p['path']:p['sha256'] for p in entries}
    need(len(pins)==len(entries),'unique universal ledger paths')
    pins.update(EXTRA_PINS)
    for name,digest in pins.items():
        if name in c.PINS:need(c.PINS[name]==digest,'accepted/base pin collision')
        c.PINS[name]=digest;need(c.digest(name)==digest,'universal proof/source/data identity')
        c.READ.add(name)
    return a,r,c,extra,ledger


def domain(c,rows,key,expected,message):
    need(type(rows) is list and len(rows)==len(expected),message)
    observed=[key(p) for p in rows]
    need(len(observed)==len(set(observed)) and set(observed)==set(expected),message)
    return {key(p):p for p in rows}


def context(c,g,white):
    return sha256(c.canonical([g,white]).encode()).hexdigest()


def role(c,g,white,marks):
    physical=c.E.matrix(g);ids=sorted(set(white)|set(marks));out=[]
    for nodes in c.E.components(physical,ids,c.FUEL):
        local=[m for m in marks if m in nodes]
        if not local:continue
        typ,images,ideal=c.E.classify(physical,nodes,local,c.W,c.FUEL)
        out.append(dict(nodes=nodes,marks=local,type=list(typ),affine=ideal,
                        all_native_marked_images=[list(v) for v in images]))
    return out


def current(c,g,white,marks,expected,extra=(),complete=False):
    result=c.E.corner(g,white,marks,c.W,c.FUEL,extra=extra,complete=complete)
    need(F(result['weight'])==F(expected),'exact CURRENT receiver value')
    return dict(actual_corner=marks,current=result,native_roles=role(c,g,white,marks))


def receipt(c,key,g,white,corners,values):
    need(len(corners)==len(values) and len({tuple(sorted(m)) for m in corners})==len(corners),
         'distinct actual receipt corners')
    out=[current(c,g,white,m,v) for m,v in zip(corners,values)]
    total=sum((F(p['current']['weight']) for p in out),F())
    need(total>=1,'CURRENT distinct-corner sum at least one')
    return dict(application_key=key,physical_context_sha256=context(c,g,white),
                current_receivers=out,current_total=str(total))


def assumed(kind,**fields):
    return dict(kind=kind,runtime_verified=False,**fields)


def application(key,theorem,selector,assumptions,**fields):
    return dict(application_key=key,theorem=theorem,exact_source_selector=selector,
        status='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION',
        conclusion='current whole face payment>=1',assumptions=assumptions,**fields)


def order12(c):
    repaired=c.field(REPAIR,('affine_D11','rows'))
    wanted=list(product(ORDER12_MARKS,('left','right'),('isolated4','mate3_1','D5','A1_A3')))
    by=domain(c,repaired,lambda p:(tuple(p['old_marks']),p['spectator_free_extreme_side'],p['branch']),
              wanted,'exact CURRENT32 order12D repair key domain')
    out=[]
    for (old,side,branch),p in by.items():
        a,b=old;c0,d=12,13;quiet,other=(c0,d) if side=='left' else (d,c0)
        near,far=(a,b) if side=='left' else (b,a)
        mate=lambda v:17-v if v<10 else 21-v
        extra={'isolated4':1,'mate3_1':0,'D5':5,'A1_A3':4}[branch]
        edges=[(quiet,mate(near))]
        if branch=='isolated4':edges += [(other,far,'sqrt(2)'),(other,14)]
        elif branch=='mate3_1':edges += [(other,far),(other,mate(far))]
        elif branch=='D5':edges += [(other,far),(14,15),(15,16),(16,17),(16,18),(other,17)]
        else:edges += [(other,far),(other,14),(other,15),(15,16),(16,17)]
        g=c.extend(c.affine('D',12),2+extra,edges)
        white=sorted(set(range(len(g)))-{a,b,c0,d})
        corners=[[far,quiet], [quiet,other] if branch in ('isolated4','mate3_1') else [near,other]]
        need(c.E.matrix(g)==c.E.matrix(p['gram']) and p['white']==white and
             [q['marks'] for q in p['corners']]==corners,'same complete32 physical Grams/whites/corners')
        values=['1/2','1/2' if branch=='isolated4' else '1']
        row=receipt(c,[list(old),side,branch],g,white,corners,values)
        row['exact_source_selector']=[REPAIR,'affine_D11','rows',repaired.index(p)]
        out.append(row)
    both=c.field(D11,('two_A1_A3_spectators','rows'))
    wanted=list(product(ORDER12_MARKS,(8,9),(10,11),(12,13)))
    by=domain(c,both,lambda p:(tuple(p['marks']),*p['leaf_attachments'],p['A1_contacting_black']),
              wanted,'separate32 BOTH-nonempty order12D key domain')
    sealed=[]
    for (old,left,right,small),p in by.items():
        a,b=old;other=25-small
        g=c.extend(c.affine('D',12),6,[(12,left),(13,right),(small,14),(other,15),(15,16),(16,17)])
        white=sorted(set(range(12))-{a,b})+[14,15,16,17]
        need(c.E.matrix(g)==c.E.matrix(p['gram']) and p['white']==white and p['signature']==[16,1,1],
             'same complete separate32 physical Gram/white/signature record')
        need([q['marks'] for q in p['corners']]==[[b,12],[12,13],[13,a]],'original separate32 corner incidence')
        row=receipt(c,[list(old),left,right,small],g,white,[[b,12],[13,a]],['1','1'])
        row['exact_source_selector']=[D11,'two_A1_A3_spectators','rows',both.index(p)]
        sealed.append(row)
    assumptions=[assumed('accepted_complete_all_native_profile_Schur_and_whole_factor_propositions',source=BRIDGE),
                 assumed('independent_current_S5',theorem='finite57.current.S5.whole_face.v1')]
    apps=[application(['D',12,*p],'finite57.current.affine_order12_D_opposite_leaves.whole_face.v1',
            [D11,'removed_pairs'],assumptions,kind='D',order=12,marked_image=list(p),
            finite_receiver_domains=['current32','separate32'],finite_receiver_S5_dependency=False) for p in ORDER12_MARKS]
    need(all(c.W.affine_weight('D',12,*p)==0 for p in ORDER12_MARKS),'order12D is a zero priority, never intrinsic unit')
    return dict(current32=out,separate32=sealed,theorem_applications=apps)


def zero_kernel(c,g,vector):
    physical=c.E.matrix(g);need(len(vector)==len(g),'exact immediate kernel order')
    values=[c.E.number(v) for v in vector]
    need(all((-v).negative() for v in values),'strictly positive original immediate kernel')
    for row in physical:
        total=c.E.Radical.scalar(0)
        for x,y in zip(row,values):total=total+x*y
        need(not total.terms,'exact zero original immediate kernel product')


def local_immediates(c):
    """Only missing immediate selectors; existing finite checks stay separate."""
    configs=(
      (BC12,'immediate_affine_roles',12,'B',(8,10),('end4','next3'),[(8,7,'sqrt(2)',10),(8,6,'1',10)]),
      (BC13,'immediate_affine_roles',13,'B',(8,11),('end4','next3'),[(8,7,'sqrt(2)',11),(8,6,'1',11)]),
      (BC13,'immediate_affine_roles',13,'C',(2,10),('left_next3','right_next3'),[(2,4,'1',10),(10,8,'1',2)]),
      (D13,'immediate_affine_interfaces',13,'D',(0,11),('inner4','next3'),[(0,1,'sqrt(2)',11),(0,2,'1',11)]),
      (D14,'immediate_corners',14,'D',(1,12),('inner4','next3'),[(1,2,'sqrt(2)',12),(1,3,'1',12)]),
      (D15,'D_immediate_corners',15,'D',(2,13),('inner4','next3'),[(2,3,'sqrt(2)',13),(2,4,'1',13)]),
    )
    out=[]
    for stem,selector,n,kind,selected,names,recipes in configs:
        rows=c.field(stem,(selector,))
        chosen=[p for p in rows if p.get('kind',kind)==kind]
        by=domain(c,chosen,lambda p:p['profile'],names,'exact immediate native selector domain')
        for name,(omit,attach,magnitude,far) in zip(names,recipes):
            p=by[name];g=c.extend(c.affine(kind,n),1,[(n,attach,magnitude)])
            white=sorted(set(range(n+1))-{omit,far,n});marks=[far,n]
            actual=c.E.corner(g,white,marks,c.W,c.FUEL)
            physical=c.E.matrix(g);blocks=c.E.components(physical,sorted(set(white)|set(marks)),c.FUEL)
            nodes=next(b for b in blocks if set(marks)<=set(b))
            stored_nodes=p.get('component',p.get('nodes'))
            need(set(stored_nodes)==set(nodes),'immediate full component equals source physical node set')
            sub=[[g[i][j] for j in stored_nodes] for i in stored_nodes]
            if len(p['gram'])==len(g):need(c.E.matrix(p['gram'])==physical,'full original immediate Gram')
            else:need(c.E.matrix(p['gram'])==c.E.matrix(sub),'original labelled immediate sub-Gram')
            kernel=p.get('kernel',p.get('positive_kernel'))
            zero_kernel(c,sub,kernel)
            if 'marked_pairs' in p:
                need(actual['marked_images']==p['marked_pairs'],'ALL source native marked images match immediate selector')
            else:
                need(p['canonical_pair'] in actual['marked_images'],'source labelled canonical immediate image belongs to complete native orbit')
                typ,images,ideal=c.E.classify(physical,nodes,marks,c.W,c.FUEL)
                need(ideal and actual['marked_images']==[list(v) for v in images],
                     'ALL immediate native images, including source canonical orientation')
            need(actual['kind']=='ideal','immediate receiver is actually ideal')
            expected=p.get('current_weight',p.get('retained_literal_weight',p.get('weight')))
            need(F(actual['weight'])==F(expected),'CURRENT immediate literal equals saved exact selected role')
            own=c.W.affine_weight(kind,n,*selected)
            priority=[]
            if not F(expected):
                priority=[assumed('independent_small_ideal_priority',theorem='finite57.current.BC11_zero.whole_face.v1')]
            else:need(own+F(expected)>=1,'selected plus DISTINCT immediate current output pays')
            out.append(dict(application_key=[stem,kind,name],exact_source_selector=[stem,selector,rows.index(p)],
                physical_context_sha256=context(c,g,white),selected_corner=list(selected),
                distinct_actual_corner=marks,selected_intrinsic_current=str(own),current_receiver=actual,
                native_roles=role(c,g,white,marks),assumptions=priority,
                conclusion='independent whole-face priority first' if priority else 'distinct current payments sum>=1'))
    finite_links=[
      dict(source=BC12,selector='whole_finite_corners',worker='scripts/check_finite57_g4_lower_cycle_terminal_bindings_20261025.py',batch='order12-near'),
      *[dict(source=s,selector=sel,worker='scripts/check_finite57_terminal_collectors_syntax_safe_20261025.py',batch=b)
        for s,sel,b in ((BC13,'whole_finite_corners','universal-bc13'),(D13,'finite_corner_representatives','universal-d13'),
                        (D14,'finite_whole_component_corners','universal-d14-near'),(D15,'D_whole_finite_corners','universal-d15-basic'))],
    ]
    domains=[('finite57.current.BC12_near_terminal.whole_face.v1','B',12,[(8,10),(8,11)]),
             ('finite57.current.BC13_inner_terminal.whole_face.v1','B',13,[(8,11),(8,12)]),
             ('finite57.current.D13_center_leaf.whole_face.v1','D',13,[(0,11),(0,12),(8,9),(8,10)]),
             ('finite57.current.D14_near_fork.whole_face.v1','D',14,[(1,12),(1,13),(8,10),(8,11)]),
             ('finite57.current.D15_basic_fractional.whole_face.v1','D',15,[(2,13),(2,14),(8,11),(8,12)])]
    apps=[]
    for theorem,kind,n,pairs in domains:
        for pair in pairs:
            apps.append(application([kind,n,*pair],theorem,[BRIDGE,'Sections2–5'],
                [assumed('accepted_complete_all_native_profile_Schur_and_whole_factor_propositions',source=BRIDGE),
                 assumed('independent_small_ideal_priority',theorem='finite57.current.BC11_zero.whole_face.v1'),
                 assumed('exact_existing_finite_selector_receipts',requirements=finite_links)],
                kind=kind,order=n,marked_image=list(pair)))
    cf=c.field(BC13,('formulas','C'))
    need(cf['resistance']==8 and F(cf['consequences']['no_next3_strict_gap'])==F(1,4) and
         cf['consequences']['next3_payment']=='whole_affine_B11_retained_distance7_unit',
         'independent C13FIRST exact R8 and strictgap1/4 profile propositions')
    c_immediate=[p for p in out if p['application_key'][:2]==[BC13,'C']]
    need(len(c_immediate)==2 and all(p['current_receiver']['forest']==[['B',11]] and
         p['current_receiver']['marked_images']==[[6,9],[6,10]] and p['current_receiver']['weight']=='1' for p in c_immediate),
         'C13FIRST two distinct native next3 choices are wholeB11 distance7 CURRENT units')
    c_first=application(['C',13,2,10],'finite57.current.C13_inner_first_zero.whole_face.v1',
        [BC13,'formulas','C'],[assumed('accepted_complete_C13_native_profile_and_whole_factor_propositions',source=BRIDGE)],
        kind='C',order=13,marked_image=[2,10],resistance=8,no_next3_strict_gap='1/4',
        next3_receiver_keys=[p['application_key'] for p in c_immediate],independent_S8_dependency=False)
    apps.append(c_first)
    c12=c.field(BC12,('formulas','C'))
    need(c12['resistance']==9 and c12['marks']==[1,10] and
         c.W.affine_weight('C',12,1,10)==0 and F(c12['consequences']['uniform_strict_gap'])==F(1,9),
         'C12 companion independent exact R9 and strictgap1/9 profile propositions')
    apps.append(application(['C',12,1,10],'finite57.current.C12_near_companion_zero.whole_face.v1',
        [BC12,'formulas','C'],[assumed('accepted_complete_C12_native_profile_and_whole_factor_propositions',source=BRIDGE)],
        kind='C',order=12,marked_image=[1,10],resistance=9,uniform_strict_gap='1/9',
        proof='Complete C12 costs<=4+4 at R9 leave strictgap1/9; otherwise an actual adjacent short1/2 path pays independently.',
        independent_S8_dependency=False))
    return dict(immediate_receiver_joins=out,required_existing_finite_receipts=finite_links,theorem_applications=apps,
        independent_C_companions=[dict(kind='C',order=12,marked_image=[1,10],
            proof='R9 costs<=4+4 gives gap1/9; an actual short1/2 path pays independently'),
            dict(kind='C',order=13,marked_image=[2,10],
            proof='Next3 gives the bound wholeB11 distance7 current unit; otherwise R8 costs<=3+3 gives gap1/4')],
        finite_selectors_recomputed=False)


def protected_option(a,r,c,p,core,option):
    """OLD protected option only; no generic currentprotected74 premise."""
    marks=option['preferred_corner'];nodes=option['marked_component'];g=core['gram']
    physical=c.E.matrix(g);white=sorted(set(range(len(g)))-set(marks)-set(option['omitted_black']))
    original=next(q for q in p['all_corners'] if q['preferred_corner']==marks)
    need(all(option[k]==v for k,v in original.items()),'OLD option is the SAME original physical corner record')
    need(sorted(tuple(b) for b in c.E.components(physical,sorted(set(white)|set(marks)),c.FUEL))==
         sorted(tuple(b['nodes']) for b in original['blocks']),'SAME complete physical OLD corner component partition')
    adj=a.adjacency(c,g,nodes);path=a.path(adj,marks)
    forbidden={v for v in nodes if v not in option['possible_contacts'] or F(core['gauge'][v])*F(core['inverse'][v][v])>2}
    leaves=[v for v in nodes if len(adj[v])==1 and v in forbidden]
    simple=[v for v in leaves if next(iter(adj[v])) in forbidden and
            len(adj[next(iter(adj[v]))])==2 and adj[v][next(iter(adj[v]))]==3]
    marked=[m for m in marks if m in leaves]
    features=[v for v in nodes if len(adj[v])>=3 or any(x!=3 for x in adj[v].values())]
    near=min((a.distances(adj,[m])[f] for m in marks for f in features),default=999)
    need(len(path)==10 and option['distance']==9 and marked==option['protected_marked_leaves'] and
         option['base_order']==len(nodes) and near==option['feature_distance'] and
         simple==option['protected_simple_leaves'] and option['forced_finite']==bool(simple and len(nodes)>=10) and
         option['ideal_basic_pair']==(near<=5),
         'OLD exact leaf9, permitted-contact seal, native feature and physical path association')
    labels=c.E.labels(physical,nodes);small=[]
    for ideal,kinds in ((False,('A','B','D')),(True,('A','B','C','D'))):
        if ideal and option['forced_finite']:continue
        for order in range(len(nodes),18):
            for kind in kinds:
                target=c.W.classical_graph(kind,order,affine=ideal)
                for mp in r.embeddings(labels,target,c.FUEL,c.E):
                    if any(len(target[mp[nodes.index(v)]])!=len(adj[v]) for v in forbidden):continue
                    pair=tuple(sorted(mp[nodes.index(v)] for v in marks))
                    need(c.W._distances(target,pair[0])[pair[1]]==9 and
                         any(len(target[mp[nodes.index(v)]])==1 for v in marked),'persistent OLD sealed marked leaf9')
                    key=dict(high_index=p['index'],core_key=p['core_key'],base_nodes=nodes,selected_marks=marks,
                             target_kind=kind,target_order=order,affine=ideal,induced_map=mp,marked_image=list(pair))
                    value=c.W.affine_weight(kind,order,*pair) if ideal else c.W.finite_same_weight((kind,order),*pair)
                    priority=None
                    if value<F(1,2):
                        need(ideal and (kind,order,*pair) in PRIORITIES,'OLD protected option has unsupported CURRENT image')
                        priority=PRIORITIES[kind,order,*pair]
                    small.append(dict(key=key,intrinsic_current=str(value),
                        conditional_whole_face_priority=priority,half_available_without_priority=value>=F(1,2)))
    need(len({c.canonical(q['key']) for q in small})==len(small),'unique OLD option native map keys')
    return dict(actual_corner=marks,white=white,physical_context_sha256=context(c,g,white),base_nodes=nodes,
        retained_path=path,protected_marked_leaves=marked,sealed_core_nodes=sorted(forbidden),
        complete_small_native_maps=small,
        large_native_domain=dict(finite_ranks=[18,58],marked_distance=9,marked_leaf=True,current_lower='1/2',
            ideal_orders=[18,58],ideal_feature_distance_at_most=5,ideal_current='1',
            classification_premise='Whole finite/affine classification and OLD contact seals preserve this marked leaf and feature; affineA has no leaf.'),
        generic_current_protected74_used=False)


def large_roles(c):
    """Complete typed large domain, independent of genericprotected74."""
    digest=sha256();counts={'finite':0,'ideal':0}
    for ideal,kinds in ((False,('A','B','D')),(True,('B','C','D'))):
        for order in range(18,59):
            for kind in kinds:
                target=c.W.classical_graph(kind,order,affine=ideal)
                features=c.W._features(target,finite=not ideal)
                ds=[c.W._distances(target,i) for i in range(order)]
                for x,y in combinations(range(order),2):
                    c.E.charge(c.FUEL)
                    if ds[x][y]!=9 or not any(len(target[m])==1 for m in (x,y)):continue
                    if ideal and min((ds[m][f] for m in (x,y) for f in features),default=999)>5:continue
                    value=c.W.affine_weight(kind,order,x,y) if ideal else c.W.finite_same_weight((kind,order),x,y)
                    need(value==1 if ideal else value>=F(1,2),'independent OLD large protected typed lower')
                    counts['ideal' if ideal else 'finite']+=1
                    digest.update(c.canonical([ideal,kind,order,x,y,str(value)]).encode()+b'\n')
    return dict(counts=counts,complete_ordered_native_role_sha256=digest.hexdigest(),ambient=58,
                generic_current_protected74_used=False)


def coupon(r,c,p,core,q):
    need(q['core_key']==p['core_key'] and c.E.matrix(q['core_gram'])==c.E.matrix(core['gram']) and
         q['core_inverse']==core['inverse'] and q['selected_marks']==p['preferred_corner'] and
         q['selected_base_nodes']==p['marked_component'],'exact OLDcoupon high-key/physical/inverse association')
    nodes=q['selected_base_nodes'];marks=q['selected_marks'];physical=c.E.matrix(q['core_gram'])
    need(13<=len(nodes)<=17 and q['small_component_order_bound'] in (1,5),'OLD whole baseD13..17 and proved small bound')
    typ,images,ideal=c.E.classify(physical,nodes,marks,c.W,c.FUEL)
    need(typ==('D',len(nodes)) and not ideal,'original OLDcoupon base is native finiteD')
    maps=q['exposed_D10'];need(len(maps)==2,'both distinct OLD exposedD20 maps')
    need(len({tuple(z['map']) for z in maps})==2,'distinct labelled originalD20 maps')
    paying_pairs={tuple(sorted((z['small_mark'],z['D10_mark']))) for z in maps}
    need(len(paying_pairs)==1,'BOTH exposedD20 maps give the SAME actual paying pair')
    paying_pair=next(iter(paying_pairs));old_leaf=q['old_opposite_leaf'];D10_mark=maps[0]['D10_mark']
    need(old_leaf in p['omitted_black'] and D10_mark in marks and old_leaf!=D10_mark,
         'old opposite leaf is ORIGINAL omitted black; D10 mark is ORIGINAL selected black')
    nonleaf=next(v for v in marks if v!=D10_mark)
    need([i for i in range(len(physical)) if i!=old_leaf and physical[old_leaf][i].terms]==[nonleaf] and
         c.E.labels(physical,[old_leaf,nonleaf])==[[0,3],[3,0]],
         'actual full-core old leaf contacts ONLY the original nonleaf selected mark, simply')
    alternate=[z for z in p['all_corners'] if tuple(sorted(z['preferred_corner']))==paying_pair]
    need(len(alternate)==1 and alternate[0]['marked_component'] is None,
         'coupon pair is exactly ONE ORIGINAL adjacent corner, in separate marked components')
    complete_maps=list(r.embeddings(c.E.labels(physical,nodes),c.W.classical_graph('D',20),c.FUEL,c.E))
    need(len(complete_maps)==2 and {tuple(mp) for mp in complete_maps}=={tuple(z['map']) for z in maps},
         'COMPLETE original inducedD20 map domain, every original node retained')
    roles=[]
    for z in maps:
        mp=z['map'];need(len(mp)==len(nodes) and len(set(mp))==len(mp) and set(mp)<=set(range(20)),
                        'complete inducedD20 map schema')
        target=c.W.classical_graph('D',20);lab=c.E.labels(physical,nodes)
        need(all(lab[i][j]==target[mp[i]].get(mp[j],0) for i in range(len(nodes)) for j in range(len(nodes))),
             'every sourceD20 image preserves ALL native edges and nonedges')
        lookup=dict(zip(mp,nodes));d10=sorted(lookup[i] for i in range(10,20))
        need(d10==z['whole_D10_nodes'] and z['D10_mark'] in marks and
             mp[nodes.index(z['D10_mark'])] in (18,19) and
             z['small_mark']==q['old_opposite_leaf'] and z['small_component_order_bound']==q['small_component_order_bound'],
             'same exposed wholeD10 physicalnodes/shortleaf/smallmark association')
        t,im,aff=c.E.classify(physical,d10,[z['D10_mark']],c.W,c.FUEL)
        need(t==('D',10) and not aff and set(im)=={(8,),(9,)},'ALL exposed nativeD10 shortleaf images')
        need(all(not physical[z['small_mark']][v].terms for v in d10),'exposed physicalD10 orthogonal to old opposite mark')
        roles.append(dict(D20_induced_map=mp,whole_D10_nodes=d10,shortleaf_mark=z['D10_mark'],
                          all_native_marked_images=[list(v) for v in im],small_mark=z['small_mark']))
    # Every possible finite small type and every marked node is covered by
    # the broad CURRENT small/D10 clause; I2(m),m>=5 is one literal class.
    small=[('A',n) for n in range(1,6)]+[('B',n) for n in range(2,6)]+[('D',4),('D',5),('F',4),('H',3),('H',4),('I',2)]
    for t in small:
        for mark in range(t[1]):
            for leaf in (8,9):need(c.W.finite_cross_weight(t,mark,('D',10),leaf)>=F(1,2),'complete whole-small CURRENT half')
    return dict(source_index=q['index'],original_selected_native_images=[list(v) for v in images],exposed_D20_maps=roles,
        actual_paying_corner=alternate[0]['preferred_corner'],original_actual_corner_selector=['all_corners',p['all_corners'].index(alternate[0])],
        small_component_order_bound=q['small_component_order_bound'],current_distinct_halves=['1/2','1/2'],
        whole_small_seal_premise='OLD complete inverse/contact and residual16 proof, including every indirect tail; displayedD20 alone is not an all-rank proof.')


def large(a,r,c,shard):
    cases,cores,_=r.domain(c)
    old=c.field(LARGE,('protected_cases',));coupons=c.field(LARGE,('coupon_cases',));routing=c.field(LARGE,('routing',))
    expected=[p for p in cases if p['method']=='finite_leaf9_or_ideal_feature5']
    by=domain(c,old,lambda p:p['index'],[p['index'] for p in expected],'exact OLDcomplete74 physical high indices')
    need(routing==dict(whole_component_at_most17=2531,adaptive_simple_leaf_neighbor_degree2=9,
        protected_simple_leaf_neighbor_degree2=18,two_protected_leaf9_corners=40,new_small_component_D10_half=16),
        'exact OLD acyclic2614 routing proposition')
    coupon_by=domain(c,coupons,lambda p:p['index'],
        [p['index'] for p in old if p['method']=='new_small_component_D10_half'],'exact OLDcoupon16 index set')
    two=[p for p in old if p['method']=='two_protected_leaf9_corners']
    need(len(two)==40 and len(coupon_by)==16,'complete OLD40+16 split, not genericprotected74')
    keys=[];out=[]
    for p in expected:
        source=by[p['index']]
        if source['method'] not in ('two_protected_leaf9_corners','new_small_component_D10_half'):continue
        key=dict(high_index=p['index'],core_key=p['core_key'],base_nodes=p['marked_component'],selected_marks=p['preferred_corner'])
        keys.append(key)
        if p['index']%4!=shard:continue
        core=cores[p['core_key']]
        options=[protected_option(a,r,c,p,core,q) for q in source['protected_options']]
        need(len({tuple(sorted(q['actual_corner'])) for q in options})==len(options),'distinct ORIGINAL protected option corners')
        coupon_row=None
        if source['method']=='two_protected_leaf9_corners':need(len(options)>=2,'OLD40 has two distinct protected halves')
        else:
            need(len(options)==1,'OLD16 has precisely one protected option')
            coupon_row=coupon(r,c,p,core,coupon_by[p['index']])
            need(set(options[0]['actual_corner'])!=set([coupon_row['exposed_D20_maps'][0]['small_mark'],
                coupon_row['exposed_D20_maps'][0]['shortleaf_mark']]),'coupon actual corner distinct from original protected half')
        assumptions=[assumed('independent_current_S8',theorem='finite57.current.S8.whole_face.v1'),
            assumed('OLD_complete74_split_and_residual16_whole_factor_proof',source=LARGE),
            assumed('independent_small_ideal_priorities_and_exact_selected_zero_exclusions',
                    source='accepted complete original protected-map/native-profile propositions; each small native priority request is keyed below')]
        out.append(application(key,'finite57.current.largeD_leaf9_distinct_half.whole_face.v1',
            [LARGE,'protected_cases',old.index(source)],assumptions,protected_options=options,coupon=coupon_row,
            original_selected_current_half=dict(type='D',rank_range=[20,58],marked_images='(r-11,r-2) or (r-11,r-1)',value='1/2'),
            generic_current_protected74_used=False))
    return dict(shard=shard,shards=4,complete_old40plus16_keys=keys,theorem_applications=out,
                independent_large_native_roles=large_roles(c),old_split_is_not_generic_current_protected74=True)


def quiet_physical(c,p):
    a,b=p['marks'];side=p['nonextreme_side'];g=c.extend(c.affine('D',13),2,[])
    for root,mark,name in ((13,a,p['nonextreme_profile'] if side=='left' else p['extreme_profile']),
                           (14,b,p['extreme_profile'] if side=='left' else p['nonextreme_profile'])):
        node=mark if name in ('isolated','isolated4') else (19-mark if mark<11 else 23-mark) if name in ('mate3','mate4') else (0 if mark<11 else 8)
        doubled=name in ('isolated4','mate4')
        if doubled:g[root][root]='4'
        g[root][node]=g[node][root]='-2' if doubled else '-1'
    return g


def final_quiet(c):
    original=c.field(U,('nonextreme_zero_cases',));repairs=c.field(c.FINAL,('quiet_repairs',))
    need(len(original)==48,'complete original order13 nonextreme48 domain')
    wanted=[3+6*i+j for i in range(8) for j in range(3)]
    by=domain(c,repairs,lambda p:p['index'],wanted,'exact finalguard24 original quiet indices')
    out=[]
    for index,p in by.items():
        q=original[index];old=ORDER13_MARKS[index//12];side=('left','right')[(index%12)//6]
        profile=('isolated4','mate4','center3')[index%3]
        need(tuple(q['marks'])==old and q['nonextreme_side']==side and q['extreme_profile']=='mate3' and
             q['nonextreme_profile']==p['profile']==profile,'exact original ordered quiet-profile source association')
        g=quiet_physical(c,q);white=sorted(set(range(13))-set(old))
        need(c.E.matrix(g)==c.E.matrix(q['gram']) and q['white']==white,'complete actual quiet Gram and original white identity')
        quiet=[14,old[0]] if side=='left' else [old[1],13]
        paying=([old[1],13] if side=='left' else [14,old[0]]) if profile!='center3' else [13,14]
        need(p['quiet_corner']==quiet and p['paying_corner']==paying and set(quiet)!=set(paying),
             'two DISTINCT actual quiet/pay corners from SAME physical star')
        for kinds in product(('A','D'),repeat=3):
            extra=tuple((kind,15) for kind in kinds)
            first=current(c,g,white,quiet,'1/2',extra,True)
            need(first['current']['complete'] and not first['current']['ignored_unsealed_components'] and
                 sum(t[1] for t in first['current']['forest'])==58,'complete guardTRUE finite rank58 quiet corner')
            premises=[assumed('accepted_full_star_zero_energy_seal',source=BRIDGE,
                assertion='R12 nonextreme8+4 gives E0; fullZ inverse SPD and disjoint columns force BOTH entire columns zero. Same three labelled whole A15/D15 factors are orthogonal to the entire star and remain whole at BOTH receivers.'),
                assumed('complete_original_native_profiles_and_actual_factor_fiber_geometry',source=U)]
            if profile=='mate4':
                second=current(c,g,white,paying,'0')
                need(second['current']['kind']=='ideal' and second['current']['forest']==[['B',13]] and
                     second['current']['marked_images']==p['images']==[[10,11],[10,12]],'exact WHOLE current B13zero priority image')
                premises.append(assumed('independent_quiet_B13_current_priority',
                    theorem='finite57.current.quiet_B13_long_orbit.whole_face.v1',
                    worker='scripts/check_finite57_quiet_b13_current_priority_bridge_20261025.py',batch='quiet-b13',
                    consumer_key=[index,*kinds],marked_images=p['images']))
                total=None
            else:
                second=current(c,g,white,paying,'1/2' if profile=='center3' else '1',extra,True)
                need(second['current']['complete'] and not second['current']['ignored_unsealed_components'] and
                     sum(t[1] for t in second['current']['forest'])==58,'same complete actual rank58 paying forest')
                total=str(F(first['current']['weight'])+F(second['current']['weight']))
            out.append(application([index,*kinds],'finite57.current.final_guard_quiet_receiver.whole_face.v1',
                [c.FINAL,'quiet_repairs',repairs.index(p)],premises,
                original_source_selector=[U,'nonextreme_zero_cases',index],
                original_source_context_sha256=context(c,q['gram'],white),
                physical_context_sha256=context(c,g,white),ordered_whole_witnesses=[dict(copy=i,type=[k,15]) for i,k in enumerate(kinds)],
                same_three_witness_copies_at_both_corners=True,quiet_receiver=first,paying_receiver=second,
                distinct_current_total=total,intrinsic_zero_credited_as_unit=False))
    need(len(out)==192 and len({tuple(p['application_key']) for p in out})==192,'exact24 quiet keys times8 ordered whole contexts')
    return dict(theorem_applications=out,quiet24_complete_at_declared_premises=True,
                independent_quiet_B13_runtime_verified=False,finalguard_one_sided_eight_not_recomputed=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=BATCHES,required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args=parser.parse_args();own=str(Path(__file__).resolve().relative_to(ROOT));before=raw(own)
    a,r,c,extra,ledger=prepare()
    if args.batch=='order12-d':result=order12(c)
    elif args.batch=='five-local':result=local_immediates(c)
    elif args.batch=='final-quiet':result=final_quiet(c)
    else:result=large(a,r,c,int(args.batch[-1]))
    for name in c.READ:need(c.digest(name)==c.PINS[name],'consumed universal source/input changed')
    need(raw(own)==before and sha256(raw(LOADER)).hexdigest()==LOADER_SHA and raw(LEDGER)==ledger,
         'universal controller/ledger changed during batch')
    need(raw(a.EXTRA)==extra and sha256(raw(a.ROUTE)).hexdigest()==a.ROUTE_SHA and
         sha256(raw(r.COLLECTOR)).hexdigest()==r.COLLECTOR_SHA and
         sha256(c.raw(c.MANIFEST)).hexdigest()==c.MANIFEST_SHA,'shared frozen loader chain changed')
    payload=json.dumps(dict(status='PASS_NAMED_CURRENT_UNIVERSAL_RECEIVER_JOINS_AT_DECLARED_PREMISES',
        batch=args.batch,source_sha256=sha256(before).hexdigest(),result=result,
        theorem_required_named_receipts=THEOREM_BATCHES,
        consumed_inputs={p:c.PINS[p] for p in sorted(c.READ)},fuel_used=c.FUEL[0],
        scope='Only exact declared finite physical/native-image receiver joins and conditional applications. Complete native profiles/classification, saved inertia/inverse and whole-factor geometry remain accepted proof inputs; no transitive regeneration.',
        independent_S8_verified=False,current_protected74_verified=False,
        current_terminal_bridge_complete=False,global_finite57_verified=False),sort_keys=True)
    need(len(payload.encode())+1<=MAX_OUTPUT,'complete receiver output cap; refusal without truncating keys')
    print(payload)


if __name__=='__main__':main()
