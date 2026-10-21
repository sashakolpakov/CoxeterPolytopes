#!/usr/bin/env python3
"""Complete ambient60 count after the finite-rank15 and affine-order15 rules.

All three actual coupon increments and the A17 restoration are included.
Generation constructs new residual DP tables; replay checks certificates.
Every ideal convolution is reconstructed in either mode.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_ambient60_payment_frontier_20261018 as Z

ROOT=Z.ROOT; N=60
TAG='fresh_rank15_order15_ambient60_frontier_20261018'
FINITE='fresh_rank15_sparse_payment_exchange_20261018'
IDEAL='fresh_order15_bd_terminal_coupon_transfer_20261018'
DEPTH='fresh_depth10_finite7_zero_transfer_20261018'
STEMS=(FINITE,IDEAL,DEPTH,Z.TAG,Z.COMPLETE,Z.B.TAG)
NEW_BASE=dict(A15='101',B15='97',D15='102',A17='116')
EXTRA={frozenset((('A',2),('D',13))):4,
       frozenset((('A',4),('D',11))):4,
       frozenset((('A',5),('D',11))):4}


def base(t):
    name=t[0]+str(t[1])
    return int(2*F(NEW_BASE[name])) if name in NEW_BASE else Z.B.base(t)


def actual_cross(a,b):return Z.current_cross(a,b)+EXTRA.get(frozenset((a,b)),0)


def upper_cross(a,b):
    return Z.B.upper_cross(a,b)+(4 if frozenset((a,b)) in set(EXTRA)-{frozenset((('A',2),('D',13)))} else 0)


def item(r,key):return Z.B.item(r,key)+(4*key[3] if r in (4,5) else 0)


def table(key):
    values=[0]+[-10**9]*N
    for r in Z.ITEMS:
        w=item(r,key)
        for s in range(r,N+1):values[s]=max(values[s],values[s-r]+w)
    last=[0]+[next(r for r in Z.ITEMS if r<=s and values[s]==values[s-r]+item(r,key)) for s in range(1,N+1)]
    return dict(key=list(key),values=values,last=last)


def check_table(t,key):
    assert t['key']==list(key)
    v,last=t['values'],t['last'];assert len(v)==len(last)==N+1 and v[0]==last[0]==0
    count=0
    for s in range(1,N+1):
        for r in Z.ITEMS:
            if r<=s:assert v[s]>=v[s-r]+item(r,key);count+=1
        r=last[s];assert r in Z.ITEMS and r<=s and v[s]==v[s-r]+item(r,key)
    return count


def finite(old,saved):
    types=[t for t in Z.B.B.types() if t[1]<=N]
    for a in types:
        assert base(Z.B.B.representative(a))>=base(a)
        for b in types:
            assert actual_cross(a,b)==actual_cross(b,a)
            assert actual_cross(a,b)<=upper_cross(Z.B.B.representative(a),Z.B.B.representative(b))
    dom=Z.large_domain();keys=sorted({key for cs,used,key in dom})
    assert len(dom)==1044 and len(keys)==792
    assert {cs for cs,used,key in dom}=={tuple(r['counts']) for r in old['finite']['alternatives']}
    tables={tuple(t['key']):t for t in saved['tables']} if saved else {k:table(k) for k in keys}
    assert set(tables)==set(keys) and (saved is None or len(saved['tables'])==len(keys))
    inequalities=sum(check_table(tables[k],k) for k in keys)
    residual=[Z.B.B.representative(('A' if r<=8 else 'D',r)) for r in Z.ITEMS]
    assert all(upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];items=0
    for cs,used,key in dom:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(t[1] for t in fixed)==used
        assert all(upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(Z.ITEMS,residual):
            assert item(r,key)==base(t)+sum(upper_cross(t,b) for b in fixed);items+=1
        rem=N-used;t=tables[key];value=sum(base(x) for x in fixed)+t['values'][rem]
        forest=fixed[:]
        while rem:
            r=t['last'][rem];forest.append(Z.B.B.representative(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==N
        assert value==sum(base(t) for t in forest)+sum(upper_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        # The deep finite7 theorem deletes these unit pairs in classical
        # factors. It is not needed for the upper, but actual witness costs
        # below must include it.
        deep=sum(max(0,r-{'A':27,'B':28,'D':29}[k]) for k,r in forest if k in 'ABD')
        actual=sum(base(t) for t in forest)+sum(actual_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])-2*deep
        assert actual<=value
        rows.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
            doubled_upper=value,doubled_actual=actual,deleted_deep_finite7_pairs=deep))
    assert max(r['doubled_upper'] for r in rows)==826
    witness=[('D',17),('D',17),('D',26)]
    assert sum(base(t) for t in witness)+sum(actual_cross(a,b) for i,a in enumerate(witness) for b in witness[i+1:])==826
    assert all(max(0,r-29)==0 for k,r in witness)
    return dict(exact_maximum='413',curvature='405',complete_fixed_vectors=1044,
        interaction_keys=792,simultaneous_type_pairs=len(types)**2,residual_item_checks=items,
        Bellman_upper_inequalities=inequalities,attaining_predecessors=len(keys)*N,
        tables=[tables[k] for k in keys],alternatives=rows,
        above_curvature=sum(r['doubled_upper']>810 for r in rows),
        actual_attainer=['D17','D17','D26'],new_tables_generated=True,replay_optimizes=False,
        added_actual_component_pair_costs=[dict(types=sorted(k+str(r) for k,r in pair),increment=str(F(v,2))) for pair,v in EXTRA.items()],
        new_A4_D11_actual_cost=str(F(actual_cross(('A',4),('D',11)),2)),
        new_A5_D11_actual_cost=str(F(actual_cross(('A',5),('D',11)),2)),
        new_A2_D13_actual_cost=str(F(actual_cross(('A',2),('D',13)),2)))


def ideal(old):
    domain=set(Z.high_parts(59))
    for low in range(1,9):domain.update((low,)+p for p in Z.high_parts(59-low))
    before={tuple(r['ranks']):r for r in old['ideal']['rows']}
    assert set(before)==domain and len(domain)==2113
    rows=[]
    for ranks in sorted(domain):
        row=Z.scalar(ranks);assert row==before[ranks]
        saving=F(2*ranks.count(14),15)
        row['normalized_cost']=str(F(row['normalized_cost'])-saving)
        row['normalized_margin']=str(F(row['normalized_margin'])+saving)
        row['absolute_margin']=str(row['link_vertices']*F(row['normalized_margin']))
        row['order15_normalized_saving']=str(saving);rows.append(row)
    summary=Z.summarize(rows)
    assert summary['negative_count']==49 and summary['zero_count']==0
    assert summary['normalized_minimum']=='-1823/11560'
    targets={11:F(1,2),12:F(1,2),13:F(1,4),15:F(5,4),16:F(3,2)}
    conditional=Z.conditional(rows,targets)
    assert conditional['summary']['negative_count']==conditional['summary']['zero_count']==0
    assert conditional['summary']['normalized_minimum']=='73/104000'
    return dict(complete_profiles=2113,rows=rows,summary=summary,
        old_negative_profiles=77,resolved_negative_profiles=28,
        low9_concentration_only=True,high_factor_concentration_used=False,
        negative_terminal_absolute_margin_is_universal_bound=False,
        remaining_conditional_target=conditional)


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[FINITE]['new_costs']==NEW_BASE and docs[FINITE]['finite_rank60_count']['exact_maximum']=='413'
    assert docs[FINITE]['upgraded_coupon']['component_pair_cost_increase']=='2'
    assert docs[IDEAL]['new_family_costs']==dict(A=75,B=92,C=91,D=92)
    assert docs[IDEAL]['uses_finite_same_component_units'] is False
    coupons=docs[IDEAL]['new_cross_component_coupons']
    assert [(r['small_type'],r['small_role'],r['small_eligible_marks'],r['large_type'],r['large_role'],r['component_pair_cost_increment']) for r in coupons]==[
        ('whole A4','endpoint',[0,3],'whole D11','short_fork_leaf','2'),
        ('whole A5','second_from_endpoint',[1,3],'whole D11','short_fork_leaf','2')]
    assert all(docs[s]['triangle_minimum']==2 and docs[s]['quadrilateral_minimum']==1 for s in (FINITE,IDEAL,DEPTH))
    assert docs[Z.COMPLETE]['dimension_bound']==60 and docs[DEPTH]['feature_depth']==10
    paths={Path(__file__),Path(Z.__file__),Path(Z.B.__file__),Path(Z.B.B.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    assert sum(2*j-59 for j in range(60) if 2*j>59)==900
    return dict(status='PASS_RANK15_ORDER15_AMBIENT60_FRONTIER',sources=sources,
        completed_unrestricted_dimension_bound=60,ambient_dimension=60,link_dimension=59,
        current_exceptional_finite_costs=NEW_BASE,current_ideal_order15_costs=docs[IDEAL]['new_family_costs'],
        finite=finite(docs[Z.TAG],saved['finite'] if saved else None),ideal=ideal(docs[Z.TAG]),
        strict_global_constant='4500',constant_counted_once=True,new_dimension_bound=None,
        new_finite_volume_example=False,
        scope='Complete simultaneous finite upper and formal attainment with all three typed increments and A17 restoration, plus all2113 fresh exact ideal convolutions. Geometric proofs are frozen inputs; no new census, external review or ambient59 theorem. Remaining ideal targets are conditional and do not include prospective finite compensation.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    result=run(saved)
    if args.verify:assert result==saved
    else:
        assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'finite',result['finite']['exact_maximum'],'negative_ideal',result['ideal']['summary']['negative_count'],flush=True)
