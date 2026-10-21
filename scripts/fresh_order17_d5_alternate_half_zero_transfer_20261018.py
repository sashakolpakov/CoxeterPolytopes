#!/usr/bin/env python3
"""Four order17 eighths replaced by two distinct whole finite halves.

Complete inherited support and physical seals, current literal interfaces,
and a complete increasing-family bound from the frozen ambient57 frontier.
No new forest optimization or native census is performed.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import fresh_ambient57_payment_frontier_20261018 as C
import fresh_d5_d12_alternate_coupon_rebalance_20261018 as O
import fresh_order17_center_contextual_zero_transfer_20261018 as I

ROOT=C.ROOT;P=C.P;T=O.T;A=O.A;Q=O.Q;D=O.D;S=O.S;M=S.M
TAG='fresh_order17_d5_alternate_half_zero_transfer_20261018'
CHANGE=M.CHANGE;EXPECTED=M.EXPECTED;enc=C.enc
LITERAL=(T.A.OVER,T.A.OLD16,T.A.OLD17,T.A.D13,T.A.D14,T.T.PARTIAL,T.D.TAG,T.T.NEW_D)
STEMS=tuple(dict.fromkeys((C.TAG,P.TAG,O.TAG,S.TAG,M.O.TAG,M.O.OVER,S.B.TAG,
    M.TAG,S.B.H.TAG,I.TAG,I.Z.TAG,I.Z.O.SMALL,C.Z.TAG,C.END.TAG,*LITERAL)))
FINAL_COUNT_MD='6b6a9582260ac280f50e299250d226eec1d2dbca0e444b114d21da5a1eb1d474'

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return C.actual_base(t)
def actual_cross(a,b):
    return C.actual_cross(a,b)+(16 if frozenset((a,b))==frozenset((('D',5),('A',12)))else 0)
def context(fs):
    co=Counter(fs)
    return C.context(fs)+(8*co['A',1]*co['D',5]if co['A',11]and not(co['A',9]and co['A',2])else 0)
def value(fs):
    return sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)

def change(g,white,marks,target):
    ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids));roles=[];images=[];delta=F()
    if sig==[len(ids),0,0]:
        blocks=D.finite_blocks(g,ids);co=Counter(b['type']for b in blocks)
        aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
        if aa is not bb:
            for sm,la,m,n in((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                if sm['type']!='D5' or T.mark_pairs(g,sm['nodes'],[m],5,'D')!=[(0,)]:continue
                if la['type']=='A12':
                    li=T.mark_pairs(g,la['nodes'],[n],12,'A')
                    if li==[(0,),(11,)]:
                        delta=F(1,8);roles.append(dict(role='D5_long_A12_endpoint',small=sm,large=la,images=li,increment=delta))
                elif la['type']=='A1' and co['A11']and not(co['A9']and co['A2']):
                    delta=F(1,8);roles.append(dict(role='A1_D5_long_with_A11',small=sm,large=la,
                        whole_A11=[b for b in blocks if b['type']=='A11'],increment=delta))
    else:
        blocks=[]
        if sig[1]==0:
            shared=next((b for b in A.components(g,ids)if set(marks)<=set(b)),None)
            if shared and len(shared)==17:
                graph=T.T.labels(g,shared)
                if graph is not None:
                    images=sorted({tuple(sorted(mp[shared.index(m)]for m in marks))for mp in M.P.embeddings(graph,target)})
                    if set(images)&CHANGE:
                        assert set(images)<=CHANGE;delta=-F(1,8)
                        roles.append(dict(role='affine_D16_four_eighths',images=images,increment=delta))
    assert len(roles)<=1
    return dict(delta=delta,roles=roles,components=blocks,ideal_images=images)

def local(docs):
    old=O.local_seals(docs);assert enc(old)==docs[O.TAG]['local_proofs']
    prior=docs[I.TAG]['local_proof'];target=prior['target']
    weights={tuple(r['pair']):F(r['weight'])for r in prior['weighted_pairs']}
    assert len(weights)==136 and sum(weights.values())==F(227,2)
    for p in CHANGE:assert weights[p]==F(1,8);weights[p]=F()
    assert sum(weights.values())==113
    _,dist,_=M.P.distance_data(target)
    assert all(dist[a][b]==11 for a,b in CHANGE)
    for auto in M.P.embeddings(target,target):assert {tuple(sorted(auto[i]for i in p))for p in CHANGE}==CHANGE
    for field in('retained_basic_pairs','retained_short8','retained_triangle_pairs'):
        assert all(weights[tuple(p)]==1 for p in prior[field])
    support=docs[M.TAG]['reconstructed_support'];assert support==docs[I.Z.TAG]['reconstructed_support']
    small=docs[I.Z.O.SMALL];patterns=[p for p in small['patterns']if p['order_bound']>=17 and len(p['nodes'])<=17]
    record=next(r for r in small['records']if r['order']==17 and r['kind']=='D')
    assert len(patterns)==140 and len(support)==len(record['case_embeddings'])==71
    assert sum(len(r['embeddings'])for r in support)==488
    hits=[]
    for row,saved in zip(support,record['case_embeddings'],strict=True):
        assert row['index']==saved['case']and row['embeddings']==saved['embeddings']and row['marked_images']==saved['marked_images']
        p=next(p for p in patterns if p['index']==row['index'])
        maps=[mp for mp in row['embeddings']if tuple(sorted(mp[i]for i in p['marks']))in CHANGE]
        if maps:
            assert len(maps)==8 and len(p['nodes'])==17;hits.append(dict(index=row['index'],maps=maps))
    assert [r['index']for r in hits]==EXPECTED
    proofs=[]
    for p in old['proofs']:
        core=p['old_seal'];g=A.matrix(core['gram']);white=core['white'];corners=[]
        for q in p['corners']:
            ch=change(g,white,q['marks'],target);pay=F(q['payment'])+ch['delta'];assert pay>=0
            corners.append(dict(marks=q['marks'],previous=q['payment'],**ch,payment=pay))
        assert sorted(c['delta']for c in corners)==[-F(1,8),F(),F(),F(1,8)]
        assert sorted(c['payment']for c in corners)==[F(),F(),F(1,2),F(1,2)]
        paying=[c for c in corners if c['payment']];assert paying[0]['marks']!=paying[1]['marks']
        if p['index']in(2371,2539):
            inc=next(c for c in corners if c['delta']>0)
            assert Counter(b['type']for b in inc['components'])=={'A1':1,'A11':1,'D5':1}
        proofs.append(dict(index=p['index'],physical_seal=core,corners=corners,total='1'))
    return dict(target=target,changed_pairs=sorted(CHANGE),old_weight='1/8',new_weight='0',old_D_cost='227/2',new_D_cost='113',
        family_totals=dict(A='85',B='111',C='108',D='113'),new_envelope='113/17',weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(weights.items())],
        complete_patterns=140,complete_successful_contexts=71,complete_maps=488,affected_contexts=hits,affected_maps=32,
        physical_proofs=proofs,retained_basic_pairs=prior['retained_basic_pairs'],retained_short8=prior['retained_short8'],
        retained_triangle_pairs=prior['retained_triangle_pairs'],native_search_repeated=False)

def literal(docs,target):
    rows=[];hits=Counter()
    old=docs[P.TAG]['current_interfaces']['literals']
    for prev,(stem,case,g,white,cs)in zip(old,T.literal_inputs(docs),strict=True):
        assert(prev['source'],prev['case'])==(stem,enc(case));corners=[]
        for q,(marks,*_)in zip(prev['corners'],cs,strict=True):
            assert q['marks']==marks;ch=change(g,white,marks,target);pay=F(q['lower'])+ch['delta'];assert pay>=0
            hits.update(r['role']for r in ch['roles']);corners.append(dict(marks=marks,previous=q['lower'],**ch,lower=pay))
        total=sum(c['lower']for c in corners);assert total>=1
        rows.append(dict(source=stem,case=case,corners=corners,total=total))
    assert len(rows)==216 and hits['affine_D16_four_eighths']==4
    assert hits['D5_long_A12_endpoint']==2 and hits['A1_D5_long_with_A11']==2
    return dict(rows=rows,configurations=216,changed_roles=dict(hits),minimum_total=min(r['total']for r in rows))

def finite(previous):
    assert previous['actual_finite']['exact_actual_maximum']=='372'
    actual={(tuple(r['counts']),r['k5'],r['a12']):r for r in previous['actual_finite']['classes']}
    rows=[];refined=set()
    for family,source in(('ordinary','rows'),('quarter','quarter_rows')):
        for r in previous['finite_upper'][source]:
            cs=r['counts'];k=r['k5'];a=r['a12'];rem=r['residual_rank'];key=(tuple(cs),k,a)
            # A12 and A11 are relaxed to their complete representative fibers;
            # every A1 consumes one of the available residual ranks.
            if not k or not(a or cs[2]and rem):continue
            delta=F(k*a,4)+(F(k*rem,8)if cs[2]else 0)
            base=F(r['upper']);method='complete_coarse_upper'
            if base+delta>372:
                base=F(actual[key]['maximum']);method='complete_actual_class_upper';refined.add(key)
            upper=base+delta;assert upper<372
            rows.append(dict(family=family,counts=cs,k5=k,a12=a,residual_rank=rem,old_upper=r['upper'],
                delta_upper=delta,retained_base=base,new_upper=upper,method=method))
    assert len(refined)==3 and max(r['new_upper']for r in rows)==F(739,2)
    witness=[('D',17),('D',17),('D',23)];assert sum(t[1]for t in witness)==57
    assert value(witness)==C.value(witness)==64*372
    return dict(complete_increasing_family_rows=rows,three_refined_keys=sorted(refined),increasing_family_upper='739/2',
        inherited_complete_classes=2992,new_finite_maximum='372',unchanged_actual_attainer=witness,
        new_global_DP=False,above_curvature_classes_recounted=False,old_38_count_not_reasserted=True)

def ideal(previous):
    rows=[]
    assert previous['complete_profiles']==2270 and len(previous['negative_profiles'])==41
    assert F(previous['current_envelopes']['16'])==F(227,34)
    for r in previous['rows']:
        saving=F(r['ranks'].count(16),34);margin=F(r['normalized_margin'])+saving
        rows.append(dict(ranks=r['ranks'],old_margin=r['normalized_margin'],normalized_saving=saving,
            normalized_margin=margin,absolute_margin=r['link_vertices']*margin))
    bad=[r for r in rows if r['normalized_margin']<0];zero=[r for r in rows if r['normalized_margin']==0]
    assert len(bad)==34 and not zero
    assert min(r['normalized_margin']for r in rows)==-F(973,11520)
    return dict(complete_profiles=len(rows),rows=rows,negative_profiles=bad,zero_profiles=zero,
        old_negative_count=41,new_negative_count=34,new_order17_cost='113',new_normalized_envelope='113/17',
        no_reconvolution=True,no_negative_surplus_extension=True)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    assert digest(f'results/{C.TAG}.md')==FINAL_COUNT_MD
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for path,h in d['sources'].items():assert digest(path)==h,(s,path)
    proof=local(docs);lits=literal(docs,proof['target']);ff=finite(docs[C.TAG]);ii=ideal(docs[C.TAG]['ideal'])
    assert actual_cross(('D',5),('A',12))==64*7
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_ORDER17_D5_ALTERNATE_HALF_ZERO_TRANSFER',sources={p:digest(p)for p in sorted(paths)},
        local_proof=proof,current_literal_interface=lits,finite=ff,ideal=ii,
        new_coupons=dict(D5_long_A12_endpoint_increment='1/8',D5_A12_component_increment='1/4',
            A1_D5_long_increment='1/8 only if A11 and not both A9,A2',A9_A2_branch_unchanged='1'),
        triangle_minimum='2',quadrilateral_minimum='1',completed_dimension_bound=57,new_dimension_bound=None,
        scope='Four full whole-factor exchanges, retained current universal outputs, complete changed-family bound and saved ideal update; no new dimension56 theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run()
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['new_finite_maximum'],out['ideal']['new_negative_count'],flush=True)
