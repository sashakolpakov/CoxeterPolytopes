#!/usr/bin/env python3
"""Whole finite B13 endpoint zero with two sealed existing-coupon repairs.

Complete selected maps and current chosen-output interfaces; no forest DP.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import sympy as sp
import fresh_d17_leaf_a16_context_coupled_transfer_20261018 as S
import fresh_order13_affine_b_singleton_transfer_20261017 as O
import fresh_order13_bc_inner_terminal_transfer_20261018 as BC
import fresh_order10_b_terminal_half_transfer_20261018 as B10
import fresh_small_d15_fork_coupon_provenance_20261018 as V

ROOT=S.ROOT;T=S.T;A=S.A;Q=S.Q;D=S.D;B=S.B;W=S.W
TAG='fresh_b13_endpoint_retention_20261018';PAIR=(0,12)
PROMOTE='fresh_ideal5_weighted_dimension71_20261017'
COUPON='fresh_finite8_d12_b10_coupon_budget_20261017'
FORK='fresh_order13_affine_b_fork_projection_transfer_20261017'
EXTRA=(O.TAG,BC.TAG,B10.TAG,PROMOTE,COUPON,FORK)
STEMS=tuple(dict.fromkeys((S.TAG,*S.STEMS,*EXTRA,*[s for _,ss,_ in V.ROUTES for s in ss])))
enc=S.enc

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return S.actual_base(t)-(64 if t==('B',13)else 0)
actual_cross=S.actual_cross
context=S.context
def value(fs):return S.value(fs)-64*fs.count(('B',13))

def native_images(g,nodes,marks):
    """Full induced B13 isomorphisms, invariant under positive gauges."""
    if len(nodes)!=13:return []
    labels=[[0]*13 for _ in range(13)]
    for i,j in combinations(range(13),2):
        x=sp.sympify(g[nodes[i]][nodes[j]])
        if x==0:continue
        assert x<0
        r=sp.simplify(4*x*x/(sp.sympify(g[nodes[i]][nodes[i]])*sp.sympify(g[nodes[j]][nodes[j]])))
        if r not in(1,2):return []
        labels[i][j]=labels[j][i]=3 if r==1 else 4
    target=T.target(13,'B');maps=T.P.embeddings(labels,target)
    return sorted({tuple(sorted(m[nodes.index(k)]for k in marks))for m in maps})

def corner(g,white,marks):
    ids=sorted(set(white)|set(marks));parts=A.components(g,ids)
    shared=next((p for p in parts if set(marks)<=set(p)),None)
    images=[] if shared is None else native_images(g,shared,marks)
    return dict(marks=marks,images=images,debit=F(-1)if PAIR in images else F())

def support(over):
    target=T.target(13,'B');rows=[];total=0;un=set()
    for c in over['isolation']['cases']:
        if(c.get('component_order_bound')or 0)<13 or len(c['marked_component'])>13:continue
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        maps=T.T.mapped(g,c['marked_component'],c['preferred_corner'],target)
        pairs=sorted({tuple(m['marked_images'])for m in maps});assert PAIR not in pairs
        rows.append(dict(index=c['index'],bound=c['component_order_bound'],nodes=c['marked_component'],marks=c['preferred_corner'],maps=maps,images=pairs))
        total+=len(maps);un.update(pairs)
    assert len(rows)==801 and total==174
    assert Counter(r['bound']for r in rows)=={13:556,14:108,15:87,16:50}
    assert T.P.embeddings(target,target)==[list(range(13))]
    return dict(patterns=801,maps=total,rows=rows,marked_union=sorted(un),changed_selected_maps=0)

def repairs(docs,over):
    old={p['index']:p for p in docs[O.TAG]['conditional_long_orbit']['cases']}
    out=[]
    for index in(166,1744):
        c=over['isolation']['cases'][index]
        p=O.inspect_case(c,over['isolation']['cores'][c['core_key']]);assert enc(p)=={k:old[index][k]for k in p}
        g=A.matrix(p['gram']);ww=p['white'];rows=[];pay=[];hits=[]
        for r in p['corners']:
            audit=corner(g,ww,r['marks']);payment=F();roles=[]
            if audit['debit']:hits.append(r['marks'])
            if r['kind']=='finite':
                bs=r['components'];aa,bb=[next(b for b in bs if m in b['nodes'])for m in r['marks']]
                for sm,la,lm in((aa,bb,r['marks'][1]),(bb,aa,r['marks'][0])):
                    if(sm['type'],la['type'])==('A1','D12'):
                        images=T.mark_pairs(g,la['nodes'],[lm],12,'D')
                        if images==[(10,),(11,)]:
                            payment=F(1);roles.append(dict(small=sm,large=la,large_mark=lm,images=images,weight='1'))
            if payment:pay.append(r['marks'])
            rows.append(dict(**audit,components=r.get('components',[]),roles=roles,payment=payment))
        assert len(hits)==len(pay)==1 and hits[0]!=pay[0]
        assert sum(r['payment']for r in rows)==1
        out.append(dict(index=index,full_old_seal=p,corners=rows,old_endpoint_corner=hits[0],distinct_existing_unit_corner=pay[0],total='1'))
    # The coupon source defines a global intrinsic rule, not a case-specific credit.
    rules=docs[COUPON]['additional_rules']
    assert rules is not None and any('A1 and D12'in r and 'short fork leaf'in r for r in rules)
    assert F(S.actual_cross(('A',1),('D',12)),64)>=2
    return dict(seals=out,existing_coupon_source=COUPON,global_rule_text=rules,no_coupon_increment=True)

def interfaces(docs):
    literals=[];all_b=[]
    for old,(st,case,g,w,cs)in zip(docs[S.TAG]['current_interfaces']['literals'],T.literal_inputs(docs),strict=True):
        assert(old['source'],old['case'])==(st,enc(case));rr=[]
        for prev,(marks,*_)in zip(old['corners'],cs,strict=True):
            assert prev['marks']==marks;r=corner(g,w,marks);assert not r['debit']
            rr.append(dict(**r,old_lower=prev['lower'],new_lower=prev['lower']))
        literals.append(dict(source=st,case=case,corners=rr,total=old['total']))
    assert len(literals)==216
    later=[];count=0
    for st,keys in B.FINITE.items():
        for key in keys:
            for path,p,g,w,bs in B.corners(docs[st][key],key):
                r=corner(A.matrix(g),w,p['marks']);assert not r['debit']
                later.append(dict(source=st,path=path,**r));count+=1
    assert count==564
    newer=[]
    collections=[(S.R.TAG,docs[S.R.TAG]['local_proof']['whole_proofs']),
        (S.H.TAG,docs[S.H.TAG]['local_proof']['whole_proofs']),
        (S.S.TAG,docs[S.S.TAG]['local_proof']['whole_proofs']),
        (S.N.TAG,docs[S.N.TAG]['local_proof']['proofs']),
        (S.R.R.TAG,docs[S.R.R.TAG]['local_seals']),
        (S.H.Z.U.TAG,[docs[S.H.Z.U.TAG]['local_proof']]),
        (S.TAG,docs[S.TAG]['local_proof']['old_ideal17_seals']),
        (S.TAG,docs[S.TAG]['local_proof']['A16_seals'])]
    for st,ps in collections:
        for p in ps:
            ph=p.get('old_complete_seal',p)
            for c in p['corners']:
                r=corner(A.matrix(ph['gram']),ph['white'],c['marks']);assert not r['debit']
                newer.append(dict(source=st,index=p.get('index',p.get('case',p.get('name','A14'))),**r))
    variants=[]
    for p in docs[S.L.C.R.S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            rr=[corner(A.matrix(v['gram']),v['white'],c['marks'])for c in v['corners']]
            assert not any(r['debit']for r in rr)
            variants.append(dict(index=p['index'],extra_contact=v['extra_contact'],corners=rr))
    assert len(variants)==24
    d12=[]
    for st,n in((B.PART,280),(B.REST,52)):
        ps=[(path,p)for path,p in W.walk(docs[st])if F(p['payment'])>0];assert len(ps)==n
        for path,p in ps:
            if any(b['type']=='B13'for b in p['roles']):
                assert list(PAIR)not in p['images'];all_b.append(dict(source=st,path=path,marks=p['marks'],images=p['images'],payment=p['payment']))
        d12.append(dict(source=st,positive_chosen_outputs=n))
    return dict(literals=literals,literal_changes=0,later_corners=later,later_corner_count=count,newer_corners=newer,
        D14_variants=variants,D12_positive_domains=d12,D12_B13_outputs=all_b)

def universal_roles(docs):
    fork=[]
    for p in docs[FORK]['cases']:
        if 'geometry'not in p.get('proof',{}):continue
        old=p['proof']['geometry'];g=A.matrix(old['gram'])
        for c in old['corners']:
            r=corner(g,old['white'],c['marks']);assert not r['debit']
            if r['images']:fork.append(dict(index=p['index'],**r))
    assert len(fork)==2 and all(r['images']==[(0,11)]for r in fork)
    bc=BC.finite_corners();assert BC.enc(bc)==docs[BC.TAG]['whole_finite_corners']
    bcr=[]
    for p in bc:
        for comp in p['components']:
            if set(p['paying_corner'])<=set(comp):
                images=native_images(p['gram'],comp,p['paying_corner']);assert images==[(0,10)]
                bcr.append(dict(profile=p['profile'],images=images,distance=10))
    assert len(bcr)==1
    b10=B10.finite_corners();assert enc(b10)==docs[B10.TAG]['finite_corners']
    joined=[]
    for p in b10:
        if p['role'].startswith('whole_B'):
            n=len(p['whole_components'][0]);assert n in range(10,14)
            joined.append(dict(rank=n,role=p['role'],marks=p['marked_nodes'],distance=9,endpoint12_not_used=True))
    assert len(joined)==4
    return dict(affine_B_fork=fork,BC13_inner_terminal=bcr,B10_joined=joined,
        universal_output_inventory=[dict(route=k,sources=ss,prior_scope=desc)for k,ss,desc in V.ROUTES],
        retained_near_endpoint_pair=[1,12],changed_pair_distance=12,
        current_A1_D12_shortleaf_unit_unchanged=True)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    over=docs[A.OVER]
    sel=support(over);rep=repairs(docs,over);inter=interfaces(docs);univ=universal_roles(docs)
    assert F(S.actual_base(('B',13)),64)==78 and F(actual_base(('B',13)),64)==77
    weights=[dict(pair=p,weight='0'if p==PAIR else'1')for p in combinations(range(13),2)]
    assert sum(F(r['weight'])for r in weights)==77
    paths={f'scripts/{TAG}.py',f'scripts/{V.TAG}.py',f'results/{V.TAG}.md'}|{f'{f}/{s}.{e}'for s in STEMS for f,e in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_B13_ENDPOINT_EXISTING_COUPON_RETENTION',sources={p:digest(p)for p in sorted(paths)},baseline=S.TAG,
        old_B13_cost='78',new_B13_cost='77',changed_pair=PAIR,weighted_pairs=weights,complete_selected_support=sel,
        whole_component_repairs=rep,current_interfaces=inter,universal_retained_roles=univ,
        ideal_weights_changed=False,added_coupons=[],triangle_minimum='2',quadrilateral_minimum='1',
        no_new_forest_count=True,new_dimension_bound=None,
        scope='Whole B13 endpoint unit removed using two full old affine/radical seals and an existing distinct A1/D12 shortleaf unit; complete current chosen-output retention. No forest DP or new bound.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert a.verify or not p.exists()
    out=run()
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
