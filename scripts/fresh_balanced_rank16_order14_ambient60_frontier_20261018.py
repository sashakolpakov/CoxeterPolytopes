#!/usr/bin/env python3
"""Count the compatible A/D16--17 exchanges and coupon-free affine C13 zeros.

Finite costs use eighths. Generate complete new Bellman certificates, or
check every saved upper inequality and attaining predecessor on replay.
Reconstruct all ideal convolutions and current literal interfaces.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_d16_d17_eighth_balance_20261018 as B
import fresh_a16_a17_full_payment_exchange_20261018 as A
import fresh_order14_c_distance10_terminal_projection_20261018 as I

C=B.C; Z=C.Z; N=B.N; ROOT=B.ROOT; DIM=60
TAG='fresh_balanced_rank16_order14_ambient60_frontier_20261018'
STEMS=(C.TAG,N.TAG,A.TAG,B.TAG,I.TAG,Z.TAG,Z.DEPTH,Z.Z.COMPLETE)
SHIFT={('A',16):-8,('D',16):-6,('A',17):8,('D',17):-2}


def base(t):return 2*C.base(t)+SHIFT.get(t,0)
def actual_cross(a,b):return 2*C.actual_cross(a,b)+(16 if {a,b}=={('D',5),('A',11)} else 0)
def upper_cross(a,b):return 2*C.upper_cross(a,b)
def item(r,key):return 2*C.item(r,key)-(2 if r==17 else 0)


def table(key):
    values=[0]+[-10**9]*DIM
    for r in Z.Z.ITEMS:
        w=item(r,key)
        for n in range(r,DIM+1):values[n]=max(values[n],values[n-r]+w)
    last=[0]+[next(r for r in Z.Z.ITEMS if r<=n and values[n]==values[n-r]+item(r,key)) for n in range(1,DIM+1)]
    return dict(key=list(key),values=values,last=last)


def check_table(t,key):
    assert t['key']==list(key)
    v,last=t['values'],t['last'];assert len(v)==len(last)==61 and v[0]==last[0]==0
    count=0
    for n in range(1,61):
        for r in Z.Z.ITEMS:
            if r<=n:assert v[n]>=v[n-r]+item(r,key);count+=1
        r=last[n];assert r in Z.Z.ITEMS and r<=n and v[n]==v[n-r]+item(r,key)
    return count


def finite(saved=None):
    rep=Z.Z.B.B.representative;types=[t for t in Z.Z.B.B.types() if t[1]<=DIM]
    for a in types:
        assert base(rep(a))>=base(a)
        for b in types:
            assert actual_cross(a,b)==actual_cross(b,a)
            assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)
    assert actual_cross(('D',5),('A',11))==16
    assert upper_cross(('A',5),('D',11))==76
    domain=Z.Z.large_domain();keys=sorted({key for cs,used,key in domain})
    assert len(domain)==1044 and len(keys)==792
    tables={tuple(t['key']):t for t in saved['tables']} if saved else {key:table(key) for key in keys}
    assert set(tables)==set(keys) and (saved is None or len(saved['tables'])==len(keys))
    checks=sum(check_table(tables[key],key) for key in keys)
    residual=[rep(('A' if r<=8 else 'D',r)) for r in Z.Z.ITEMS]
    assert all(upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];items=0
    for cs,used,key in domain:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for k,r in fixed)==used
        assert all(upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(Z.Z.ITEMS,residual):
            assert item(r,key)==base(t)+sum(upper_cross(t,b) for b in fixed);items+=1
        rem=DIM-used;t=tables[key];value=sum(base(x) for x in fixed)+t['values'][rem]
        forest=fixed[:]
        while rem:
            r=t['last'][rem];forest.append(rep(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==DIM
        upper=sum(base(t) for t in forest)+sum(upper_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        actual=sum(base(t) for t in forest)+sum(actual_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        assert upper==value and actual<=upper
        rows.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
            octupled_upper=value,octupled_actual=actual,deep7_deletions_included=sum(C.deep(t) for t in forest)))
    assert max(r['octupled_upper'] for r in rows)==3276
    witnesses=[]
    for ranks in ((16,16,28),(16,17,27),(17,17,26)):
        forest=[('D',r) for r in ranks]
        value=sum(base(t) for t in forest)+sum(actual_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        assert value==3276;witnesses.append(dict(forest=['D'+str(r) for r in ranks],cost='819/2'))
    assert sum(r['octupled_upper']>3240 for r in rows)==15
    return dict(exact_maximum='819/2',previous_maximum='411',curvature='405',finite_deficit='9/2',
        complete_fixed_vectors=1044,interaction_keys=792,simultaneous_type_pairs=len(types)**2,
        residual_item_checks=items,Bellman_upper_inequalities=checks,attaining_predecessors=60*len(keys),
        tables=[tables[k] for k in keys],alternatives=rows,above_curvature=15,actual_attainers=witnesses,
        exceptional_costs=dict(A16='110',D16='441/4',A17='117',D17='469/4'),
        new_D5_A11_actual_cost='2',unchanged_A5_D11_upper_allowance='19/2',finite_cost_unit='1/8',
        deep7_deletions_counted_in_upper=True,new_tables_generated=True,replay_optimizes=False)


def ideal(docs):
    previous=C.ideal(docs[Z.TAG]);assert previous==docs[C.TAG]['ideal']
    rows=[]
    for old in previous['rows']:
        row=dict(old);saving=F(row['ranks'].count(13),14)
        row['normalized_cost']=str(F(row['normalized_cost'])-saving)
        row['normalized_margin']=str(F(row['normalized_margin'])+saving)
        row['absolute_margin']=str(row['link_vertices']*F(row['normalized_margin']))
        row['order14_normalized_saving']=str(saving);rows.append(row)
    summary=Z.Z.summarize(rows)
    assert len(rows)==2113 and summary['negative_count']==17 and summary['zero_count']==0
    assert summary['normalized_minimum']=='-1117/7680'
    return dict(complete_profiles=2113,previous_negative_profiles=28,resolved_negative_profiles=11,
        rows=rows,summary=summary,low9_concentration_only=True,
        negative_terminal_absolute_margin_is_universal_bound=False)


def literal(docs):
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in C.L.STEMS}
    previous={(r['source'],C.L.N.normalize(r['case'])):r for r in docs[B.TAG]['literal_interface']['rows']}
    old_a={(r['source'],C.L.N.normalize(r['case'])):r for r in docs[A.TAG]['literal_interface']['rows']}
    target=docs[I.TAG]['target'];drops=set(map(tuple,docs[I.TAG]['deleted_pairs']))
    rows=[];changed=[];debits=credits=0
    for stem,index,g,white,cs in B.T.literal_inputs(inputs):
        key=(stem,C.L.N.normalize(index));before=previous[key];ac=old_a[key];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==before['corners'][j]['marks']==ac['corners'][j]['marks']
            ids=sorted(set(white)|set(marks));images=[];adebit=F();fdelta=F()
            sig=B.Q.inertia(B.A.sub(g,ids))
            if sig==[len(ids),0,0]:
                blocks=B.D.finite_blocks(g,ids);shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
                if shared and shared['type'] in ('A16','A17'):
                    typ=shared['type'];images=B.T.mark_pairs(g,shared['nodes'],marks,int(typ[1:]),'A')
                    change=A.ZERO if typ=='A16' else A.RESTORE
                    if set(images)&change:
                        assert set(images)<=change;fdelta=F(-1,2) if typ=='A16' else F(1,2)
                assert fdelta==F(ac['corners'][j]['delta'])
            elif sig[1]==0:
                shared=next((p for p in B.A.components(g,ids) if set(marks)<=set(p)),None)
                if shared and len(shared)==14:
                    graph=B.T.T.labels(g,shared)
                    if graph is not None:
                        maps=B.P.embeddings(graph,target)
                        images=sorted({tuple(sorted(mp[shared.index(m)] for m in marks)) for mp in maps})
                        if set(images)&drops:assert set(images)<=drops;adebit=F(1)
            value=F(before['corners'][j]['retained_payment_lower_bound'])+fdelta-adebit;assert value>=0
            delta=fdelta-adebit;debits+=delta<0;credits+=delta>0
            corners.append(dict(marks=marks,canonical_images=images,A16_A17_delta=fdelta,C14_debit=adebit,
                retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['A16_A17_delta'] or c['C14_debit'] for c in corners):changed.append(row)
    assert len(rows)==216 and debits==credits==0
    return dict(configurations=216,debited_corners=debits,credited_corners=credits,changed=changed,rows=rows,
        inherited_D16_eighth_credits=docs[B.TAG]['literal_interface']['credited_corners'])


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[N.TAG]['new_coupon']['component_pair_cost_increase']=='2'
    assert docs[A.TAG]['new_costs']==dict(A16='110',A17='117')
    assert docs[B.TAG]['new_costs']==dict(D16='441/4',D17='469/4')
    assert docs[I.TAG]['new_order14_family_costs']==dict(A=70,B=82,C=81,D=82)
    for s in (N.TAG,A.TAG,B.TAG,I.TAG):
        assert docs[s]['triangle_minimum']==2 and docs[s]['quadrilateral_minimum']==1
    assert docs[Z.DEPTH]['feature_depth']==10 and docs[Z.Z.COMPLETE]['dimension_bound']==60
    for name,doc,field in (('A16',A.TAG,'A16_weighted_pairs'),('A17',A.TAG,'A17_weighted_pairs'),
                          ('D16',B.TAG,'D16_weighted_pairs'),('D17',B.TAG,'D17_weighted_pairs')):
        assert sum(F(r['weight']) for r in docs[doc][field])==F(base((name[0],int(name[1:]))),8)
    assert F(5*sum(min(j,59-j) for j in range(60))-59*60,2)==405
    assert sum(2*j-59 for j in range(60) if 2*j>59)==900
    paths={Path(__file__),Path(B.__file__),Path(A.__file__),Path(I.__file__),Path(C.__file__),
           Path(Z.__file__),Path(Z.Z.__file__),Path(Z.Z.B.__file__),Path(Z.Z.B.B.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    return B.J.enc(dict(status='PASS_BALANCED_RANK16_ORDER14_AMBIENT60_FRONTIER',sources=sources,
        completed_unrestricted_dimension_bound=60,ambient_dimension=60,link_dimension=59,
        finite=finite(saved['finite'] if saved else None),ideal=ideal(docs),combined_literal_interface=literal(docs),
        current_ideal_order14_costs=dict(A='70',B='82',C='81',D='82'),
        current_ideal_order17_costs=docs[C.TAG]['current_ideal_order17_costs'],
        strict_global_constant='4500',constant_counted_once=True,new_dimension_bound=None,new_finite_volume_example=False,
        scope='Complete simultaneous actual/representative finite costs,1044 fixed vectors/792 fresh Bellman certificates,2113 exact ideal convolutions and216 current literal interfaces. All new finite payments and restorations counted. Universal local proofs are frozen inputs; no new geometric census or dimension59 theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    out=run(saved)
    if args.verify:assert out==saved
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'negative',out['ideal']['summary']['negative_count'],flush=True)
