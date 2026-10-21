#!/usr/bin/env python3
"""Five whole-core exchanges: D17 leaf quarters and guarded A16 five-eighths.

Exact finite accounting is in the separate ambient58 frontier. This source
reconstructs the full physical seals and current chosen-output interfaces.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import fresh_a16_a1_seven_eighth_second_coupon_20261018 as S
import fresh_a13_endpoint_a1_context_retention_20261018 as I
import fresh_ideal7_finite_cross_payments_20261017 as O

ROOT=S.ROOT;H=S.H;R=S.R;N=S.N;C=S.C;L=S.L;B=S.B;W=S.W;T=S.T;A=S.A;Q=S.Q;D=S.D
TAG='fresh_d17_leaf_a16_context_coupled_transfer_20261018'
STEMS=tuple(dict.fromkeys((S.TAG,I.TAG,*S.STEMS,*I.STEMS,O.TAG)))
DOWN={(0,15),(0,16)};CASES=[662,675,828]
enc=S.enc
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return S.actual_base(t)+{('A',17):22,('D',17):-32,('A',13):-64}.get(t,0)
def actual_cross(a,b):
    return S.actual_cross(a,b)+(32 if frozenset((a,b))==frozenset((('A',1),('A',16)))else 0)
def context(fs):
    co=Counter(fs)
    return S.context(fs)+((-32*co['A',16]+64*co['A',13]+(64*co['A',2]*co['A',4]if co['A',10]else 0))if co['A',1]else 0)
def value(fs):return sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)

def updated_lower(previous,ch):
    # The older literal ledger never booked the later D10 endpoint credit.
    # Its actual current5/32 is already included in the base-cost table.
    credit=max((r['old']-F(previous)for r in ch['roles']if r['role']=='A17_endpoints'),default=F())
    credit=max(F(),credit)
    return F(previous)+ch['delta']+credit,credit

def change(g,white,marks):
    ids=sorted(set(white)|set(marks))
    if Q.inertia(A.sub(g,ids))!=[len(ids),0,0]:return dict(delta=F(),roles=[],components=[])
    blocks=D.finite_blocks(g,ids)
    witnesses={k:[b for b in blocks if b['type']==k and not set(marks)&set(b['nodes'])]for k in('A1','A10')}
    aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
    roles=[]
    def add(name,old,new,**kw):roles.append(dict(role=name,old=old,new=new,delta=new-old,**kw))
    if aa is bb:
        kind=aa['type']
        if kind=='D17':
            images=T.mark_pairs(g,aa['nodes'],marks,17,'D')
            if set(images)&DOWN:
                assert set(images)<=DOWN
                add('D17_long_endpoint_shortleaf',F(1,2),F(1,4),images=images)
        elif kind=='A17':
            images=T.mark_pairs(g,aa['nodes'],marks,17,'A')
            if (0,16)in images:
                assert images==[(0,16)]
                add('A17_endpoints',F(5,32),F(1,2),images=images)
        elif kind=='A16':
            images=T.mark_pairs(g,aa['nodes'],marks,16,'A')
            if set(images)&L.ACHANGE:
                assert set(images)<=L.ACHANGE
                add('A16_guarded_five_eighths',F(7,8)if witnesses['A1']else F(),
                    F(5,8)if witnesses['A1']else F(),images=images,whole_A1=witnesses['A1'])
        elif I.endpoint(g,blocks,marks):
            add('inherited_A13_endpoint_A1_guard',F(1),F(1)if witnesses['A1']else F(),whole_A1=witnesses['A1'])
    else:
        for sm,la,m,lm in((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
            if(sm['type'],la['type'])==('A1','A16'):
                images=T.mark_pairs(g,la['nodes'],[lm],16,'A')
                if images==[(0,),(15,)]:add('A1_A16_endpoint_quarter',F(),F(1,4),images=images)
            if(sm['type'],la['type'])==('A2','A4'):
                images=T.mark_pairs(g,la['nodes'],[lm],4,'A')
                if images==[(0,),(3,)]:
                    add('A2_A4_endpoint_A1_A10_guarded_quarter',F(),
                        F(1,4)if witnesses['A1']and witnesses['A10']else F(),images=images,
                        whole_A1=witnesses['A1'],whole_A10=witnesses['A10'])
    assert len(roles)<=1
    return dict(delta=sum((r['delta']for r in roles),F()),roles=roles,components=blocks)

def local(docs):
    over=docs[A.OVER];support=T.selected(over,17)
    assert support['complete_maps']==168 and support['complete_eligible_patterns']==140
    assert not set(map(tuple,support['marked_pair_union']))&DOWN
    frozen=docs[H.E.SUPPORT]['complete_selected_support']['rows']
    assert len(frozen)==140
    for a,b in zip(support['rows'],frozen,strict=True):
        assert a['index']==b['index']and enc(a['maps'])==b['maps']
    inputs={(stem,case):(g,w,cs)for stem,case,g,w,cs in T.literal_inputs(docs)if isinstance(case,int)}
    current={(r['source'],r['case']):r for r in docs[S.TAG]['current_interfaces']['literals']if isinstance(r['case'],int)}
    proofs=[]
    for case in CASES:
        g,white,cs=inputs[O.TAG,case];n=len(g)
        old=next(p for p in docs[O.TAG]['order17']['replacement_cases']if p['case']==case)
        black=old['black'];assert n==19 and all(g[i][i]==2 for i in range(n))
        assert Q.inertia(g)==[17,1,1]and Q.inertia(A.sub(g,white))==[15,0,0]
        assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[16,0,0]for b in black)
        affine=sorted(set(range(n))-{3,13});af=A.sub(g,affine)
        assert Q.inertia(af)==[16,0,1]
        kernels=O.R.nullspace(af);assert len(kernels)==1
        kernel=kernels[0]
        if all(x<0 for x in kernel):kernel=[-x for x in kernel]
        assert all(x>0 for x in kernel)and all(sum(af[i][j]*kernel[j]for j in range(17))==0 for i in range(17))
        radical=list(map(F,old['radical']))
        assert len(radical)==19 and radical[3]*radical[13]<0
        assert all(sum(g[i][j]*radical[j]for j in range(n))==0 for i in range(n))
        assert old['possible_contacts']==[3,13]
        rows=[]
        for prev,(marks,*_)in zip(current[O.TAG,case]['corners'],cs,strict=True):
            assert prev['marks']==marks;ch=change(g,white,marks)
            payment,credit=updated_lower(prev['lower'],ch);assert payment>=0
            rows.append(dict(marks=marks,previous=prev['lower'],**ch,previous_unbooked_credit=credit,payment=payment))
        assert sorted(r['payment']for r in rows)==[0,F(1,4),F(1,4),F(1,2)]
        proofs.append(dict(case=case,gram=g,white=white,black=black,affine_nodes=affine,
            positive_affine_kernel=kernel,full_radical=radical,possible_contacts=[3,13],
            short_path_if_nonorthogonal=[3,'new_white',13],corners=rows,total='1'))
    a16=[]
    for old in docs[S.TAG]['local_proof']['whole_proofs']:
        c=over['isolation']['cases'][old['index']]
        p=L.a_seal(c,over['isolation']['cores'][c['core_key']]);assert enc(p)==old['full_seal']
        rows=[]
        for prev in old['corners']:
            ch=change(A.matrix(p['gram']),p['white'],prev['marks'])
            rows.append(dict(marks=prev['marks'],previous=prev['payment'],**ch,payment=F(prev['payment'])+ch['delta']))
        assert sorted(r['payment']for r in rows)==[0,F(1,8),F(1,4),F(5,8)]
        added=next(r for r in rows if any(q['role'].startswith('A2_A4')for q in r['roles']))
        rr=added['roles'][0];assert len(rr['whole_A1'])==len(rr['whole_A10'])==1
        assert len({tuple(r['marks'])for r in rows if r['payment']})==3
        a16.append(dict(index=old['index'],full_seal=p,gram=p['gram'],white=p['white'],corners=rows,total='1'))
    d10=docs['fresh_order10_d_opposite_fork_fractional_transfer_20261018']
    # Endpoint provenance: the new half retains the genuine D10 five-thirty-seconds output.
    def find_table(x):
        if isinstance(x,dict):
            if 'A17_weighted_pairs'in x:return x['A17_weighted_pairs']
            for v in x.values():
                r=find_table(v)
                if r is not None:return r
        return None
    table=find_table(d10);assert table is not None
    weights={tuple(r['pair']):F(r['weight'])for r in table}
    assert weights[0,16]==F(5,32);weights[0,16]=F(1,2)
    assert sum(weights.values())==F(actual_base(('A',17)),64)==F(453,4)
    assert F(actual_base(('D',17)),64)==F(227,2)
    return dict(complete_D17_patterns=140,complete_D17_maps=168,complete_selected_support=support,
        changed_D17_selected_maps=0,D17_pairs=sorted(DOWN),old_ideal17_seals=proofs,A16_seals=a16,
        A17_weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(weights.items())],
        ordinary_D17_cost='227/2',ordinary_A17_cost='453/4',ordinary_A16_cost='433/4',
        A16_cost_with_A1='219/2',A1_A16_component_pair_increment='1/2',
        new_A2_A4_component_pair_cost='1',new_guard=['whole unmarked A1','whole unmarked A10'],
        witness_multiplicity_not_counted=True)

def interfaces(docs):
    literals=[];negative=[];positive=[]
    for old,(stem,case,g,white,cs)in zip(docs[S.TAG]['current_interfaces']['literals'],T.literal_inputs(docs),strict=True):
        assert(old['source'],old['case'])==(stem,enc(case));rows=[]
        for prev,(marks,*_)in zip(old['corners'],cs,strict=True):
            assert prev['marks']==marks;ch=change(g,white,marks);payment,credit=updated_lower(prev['lower'],ch)
            assert payment>=0
            if ch['delta']<0:
                assert stem==O.TAG and case in CASES
                negative.append(dict(source=stem,case=case,marks=marks,delta=ch['delta']))
            if ch['delta']>0:positive.append(dict(source=stem,case=case,marks=marks,delta=ch['delta']))
            rows.append(dict(marks=marks,previous=prev['lower'],**ch,previous_unbooked_credit=credit,lower=payment))
        assert sum(r['lower']for r in rows)>=1
        literals.append(dict(source=stem,case=case,corners=rows,total=sum(r['lower']for r in rows)))
    assert len(literals)==216 and len(negative)==4
    later=[];count=0
    for stem,keys in B.FINITE.items():
        rows=[]
        for key in keys:
            for path,p,g,white,blocks in B.corners(docs[stem][key],key):
                ch=change(A.matrix(g),white,p['marks']);assert ch['delta']>=0,(stem,path,ch)
                rows.append(dict(path=path,marks=p['marks'],delta=ch['delta'],roles=ch['roles']));count+=1
        later.append(dict(source=stem,corners=rows))
    assert count==564
    newer=[]
    collections=[(R.TAG,docs[R.TAG]['local_proof']['whole_proofs']),
        (H.TAG,docs[H.TAG]['local_proof']['whole_proofs']),
        (S.TAG,docs[S.TAG]['local_proof']['whole_proofs']),
        (N.TAG,docs[N.TAG]['local_proof']['proofs']),
        (R.R.TAG,docs[R.R.TAG]['local_seals']),
        (H.Z.U.TAG,[docs[H.Z.U.TAG]['local_proof']])]
    for stem,proofs in collections:
        for p in proofs:
            physical=p.get('old_complete_seal',p);rows=[]
            for q in p['corners']:
                ch=change(A.matrix(physical['gram']),physical['white'],q['marks'])
                if ch['delta']<0:assert stem in(R.TAG,S.TAG)and p['index']in L.ACASES
                rows.append(dict(marks=q['marks'],delta=ch['delta'],roles=ch['roles']))
            newer.append(dict(source=stem,index=p.get('index',p.get('name','A14')),corners=rows))
    variants=[]
    for p in docs[L.C.R.S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            rows=[]
            for q in v['corners']:
                ch=change(A.matrix(v['gram']),v['white'],q['marks']);assert ch['delta']>=0
                rows.append(dict(marks=q['marks'],delta=ch['delta']))
            variants.append(dict(index=p['index'],extra_contact=v['extra_contact'],corners=rows))
    assert len(variants)==24
    d12=[]
    for stem,expected in((B.PART,280),(B.REST,52)):
        ps=[(path,p)for path,p in W.walk(docs[stem])if F(p['payment'])>0];assert len(ps)==expected
        for _,p in ps:assert not any(b['type']in('A16','D17')for b in p['roles'])
        d12.append(dict(source=stem,positive_outputs=expected,changed_types_absent=True))
    a13=docs[I.TAG]['local_proof'];a13rows=[]
    for p in a13['radical_proofs']:
        ch=change(A.matrix(p['gram']),p['white'],p['selected_marks']);assert ch['delta']==0
        a13rows.append(dict(index=p['index'],delta=ch['delta']))
    p=a13['saturated_two_witness_proof']
    for v in p['variants']:
        ch=change(A.matrix(v['gram']),v['white'],p['selected_marks']);assert ch['delta']==0
        a13rows.append(dict(index=p['index'],extra_contact=v['extra_contact'],delta=ch['delta']))
    return dict(literals=literals,literal_debits=negative,literal_credits=positive,later_packets=later,later_corners=count,
        newer_seals=newer,retained_D14_variants=variants,chosen_D12_outputs=d12,retained_A13_selected_units=a13rows)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    proof=local(docs);inter=interfaces(docs)
    costs=[]
    for t in C.TYPES:
        costs.append(dict(type=t,ordinary=F(actual_base(t),64)))
        for u in C.TYPES:assert actual_cross(t,u)<=C.upper_cross(C.representative(t),C.representative(u))
    assert F(actual_cross(('A',1),('A',16)),64)==F(3,2)
    assert F(actual_cross(('A',4),('A',12)),64)==8
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_D17_LEAF_A16_CONTEXT_COUPLED_TRANSFER',sources={p:digest(p)for p in sorted(paths)},
        baseline=S.TAG,composed_A13=I.TAG,scale=64,local_proof=proof,current_interfaces=inter,actual_component_costs=costs,
        actual_type_pair_domination=len(C.TYPES)**2,ideal_weights_changed=False,triangle_minimum='2',quadrilateral_minimum='1',
        new_finite_count=False,completed_dimension_bound=58,new_dimension_bound=None,
        scope='Five full physical seals and complete current chosen-output retention, composed with the frozen A13 guard. No forest maximum is asserted by this local packet.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
