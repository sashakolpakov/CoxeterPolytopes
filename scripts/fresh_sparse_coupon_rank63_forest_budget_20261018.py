#!/usr/bin/env python3
"""Rank63 finite-forest bound with all new sparse/ideal-transfer coupons.

Actual typed costs and a dominating representative allowance are distinct.
Saved verification checks all Bellman upper bounds and attaining predecessors.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sparse_coupon_rank63_forest_budget_20261018'
LOCAL='fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018'
B12='fresh_order12_b_center_terminal_transfer_20261018'
D13='fresh_order13_d_center_pair_half_transfer_20261018'
OLD='fresh_e6_b10_half_coupon_forest_budget_20261018'
CURRENT='fresh_order11_order17_dimension63_20261018'
STEMS=(LOCAL,B12,D13,OLD,CURRENT)
N=63


def base(t):
    kind,r=t
    if r==16:return 208 if kind=='B' else 228
    if r==17:return {'A':228,'B':210,'D':250}[kind]
    if r<=15:return r*(r-1)
    old=16*r-(66 if kind=='D' and r>=20 else 64 if kind=='D' else 68)
    deep=max(0,r-(22+{'A':8,'B':9,'D':10}[kind]))
    return old-2*deep


def representative(t):
    kind,r=t
    if t==('A',17):return t
    if r==6:return ('E',6)
    if r==9:return ('B',9)
    return ('A' if r<=8 else 'D',r)


def old_cross(a,b):
    out=0
    for small,large in ((a,b),(b,a)):
        if small[1]<=5:
            if large==('D',10):out+=2*small[1]
            elif large[0]=='D' and 11<=large[1]<=16:out+=3*small[1]
            elif large in (('A',13),('A',15)):out+=small[1]
            elif large==('A',17):out+=2*small[1]
        if (small,large) in ((('A',1),('D',11)),(('A',1),('D',12)),
                (('A',2),('B',10)),(('A',4),('D',13)),
                (('A',1),('A',11)),(('A',2),('A',11))):out+=2
        if (small,large) in ((('A',2),('D',11)),(('E',6),('D',11))):out+=4
        if (small,large)==(('A',1),('B',10)):out+=1
        if small in (('A',1),('A',2)) and large==('D',14):out+=small[1]
    return out


# These are doubled total costs, not the weight of an individual marked pair.
NEW_CROSS={
    (('A',4),('A',12)):8,
    (('A',4),('A',13)):8,
    (('D',5),('A',12)):4,
    (('A',2),('A',13)):4,
    (('A',2),('A',14)):4,
    (('E',6),('A',11)):8,
    (('A',2),('B',9)):2,
    (('A',2),('D',10)):2,
    (('A',3),('D',10)):1,
}


def actual_cross(a,b):
    return old_cross(a,b)+NEW_CROSS.get((a,b),0)+NEW_CROSS.get((b,a),0)


def upper_cross(a,b):
    # This additional allowance is ONLY in the representative upper bound.
    # The actual E6/D11 coupon remains total2; E6/A11's new total4 is absorbed.
    extra=4 if {a,b}=={('E',6),('D',11)} else 0
    return actual_cross(a,b)+extra


def domain():
    rows=[]
    for size in range(N//9+1):
        for rs in combinations_with_replacement(range(9,18),size):
            if sum(rs)>N:continue
            c=Counter(rs)
            key=(2*c[10]+3*sum(c[r] for r in range(11,17))+2*c[17],
                 c[11]+c[12],c[13],c[11],c[14],c[9],c[10])
            rows.append(dict(counts=[c[r] for r in range(9,18)],used_rank=sum(rs),key=list(key)))
    return sorted(rows,key=lambda r:r['counts'])


def residual_items(key):
    broad,extra,a4,d11,d14,b9,d10=key;out={}
    for r in tuple(range(1,9))+tuple(range(17,N+1)):
        t=representative(('A' if r<=8 else 'D',r))
        w=base(t)+(r*broad if r<=5 else 0)
        w+=(2*extra+d14 if r==1 else 4*d11+2*d14+2*b9+2*d10 if r==2
            else d10 if r==3 else 2*a4 if r==4 else 8*d11 if r==6 else 0)
        out[r]=w
    return out


def generate_table(key):
    items=residual_items(key);values=[0]+[None]*N
    # Item-first generation; verifier quantifies every possible last item.
    for r,w in items.items():
        for n in range(r,N+1):
            if values[n-r] is not None:
                value=values[n-r]+w
                if values[n] is None or value>values[n]:values[n]=value
    last=[0]+[min(r for r,w in items.items() if r<=n and values[n]==values[n-r]+w)
              for n in range(1,N+1)]
    return dict(key=list(key),values=values,last=last)


def verify_table(table,counts):
    items=residual_items(tuple(table['key']));values,last=table['values'],table['last']
    assert len(values)==len(last)==N+1 and values[0]==last[0]==0
    assert all(type(x) is int for x in values+last)
    for n in range(1,N+1):
        for r,w in items.items():
            if r<=n:
                assert values[n]>=values[n-r]+w;counts['bellman_upper']+=1
        r=last[n];assert r in items and r<=n and values[n]==values[n-r]+items[r]
        counts['attaining_predecessors']+=1


def forest_cost(forest,upper=False):
    cross=upper_cross if upper else actual_cross
    return sum(base(t) for t in forest)+sum(cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,doc in docs.items():
        assert doc['status'].startswith('PASS'),stem
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[CURRENT]['dimension_bound']==63
    assert docs[LOCAL]['new_finite_component_costs']==dict(A16=114,B16=104,D16=114,D17=125)
    expected_local=[
        ('A4',[0,3],'A12',[0,11],'1','4'),
        ('A4',[0,3],'A13',[1,11],'1','4'),
        ('D5',[0],'A12',[1,10],'1','2'),
        ('A2',[0,1],'A13',[0,12],'1/2','2'),
        ('A2',[0,1],'A14',[1,12],'1/2','2'),
        ('E6',None,'A11',[1,9],'1','4')]
    assert [(c['small'],c.get('small_nodes'),c['large'],c['large_nodes'],c['weight'],c['pair_cost'])
            for c in docs[LOCAL]['coupons']]==expected_local
    assert docs[LOCAL]['coupons'][-1]['small_role']=='length-two-arm endpoint'
    assert docs[LOCAL]['coupons'][-1]['small_role_count']==2
    assert docs[B12]['new_finite_coupon']==dict(small='whole A2',small_marks='either endpoint',
        large='whole B9',large_mark='doubled terminal',pair_weight='1/2',component_pair_cost='1')
    assert docs[D13]['new_coupons']==[
        dict(small='whole_A2',small_role='either_endpoint',large='whole_D10',large_role='fork_center',
             weight='1/2',component_pair_increment='1'),
        dict(small='whole_A3',small_role='unique_middle',large='whole_D10',large_role='fork_center',
             weight='1/2',component_pair_increment='1/2')]
    assert [F(c['pair_cost'])*2 for c in docs[LOCAL]['coupons']]==list(NEW_CROSS.values())[:6]
    assert F(docs[B12]['new_finite_coupon']['component_pair_cost'])*2==NEW_CROSS[(('A',2),('B',9))]
    assert [F(c['component_pair_increment'])*2 for c in docs[D13]['new_coupons']]==[
        NEW_CROSS[(('A',2),('D',10))],NEW_CROSS[(('A',3),('D',10))]]
    types=([('A',r) for r in range(1,N+1)]+[('B',r) for r in range(2,N+1)]+
           [('D',r) for r in range(4,N+1)]+[('E',6),('E',7),('E',8),('F',4),('H',3),('H',4),('I2*',2)])
    domination=[];pair_checks=[]
    for t in types:
        rep=representative(t);assert base(rep)>=base(t)
        domination.append(dict(type=list(t),representative=list(rep),base_gap=base(rep)-base(t)))
        for u in types:
            ru=representative(u);gap=upper_cross(rep,ru)-actual_cross(t,u)
            assert gap>=0,(t,u,rep,ru,gap)
            pair_checks.append(dict(types=[list(t),list(u)],representatives=[list(rep),list(ru)],
                                   actual_doubled=actual_cross(t,u),upper_doubled=upper_cross(rep,ru)))
    assert len(types)==192 and len(pair_checks)==36864
    dom=domain();keys=sorted({tuple(r['key']) for r in dom})
    assert len(dom)==1329 and len(keys)==993
    if saved is None:tables={k:generate_table(k) for k in keys}
    else:
        assert [r['key'] for r in saved['tables']]==list(map(list,keys))
        tables={tuple(r['key']):r for r in saved['tables']}
    counts=Counter()
    for table in tables.values():verify_table(table,counts)
    alternatives=[]
    for row in dom:
        large=[('B' if r==9 else 'A' if r==17 else 'D',r)
               for r,c in zip(range(9,18),row['counts']) for _ in range(c)]
        assert all(upper_cross(a,b)==0 for i,a in enumerate(large) for b in large[i+1:])
        key=tuple(row['key']);table=tables[key];items=residual_items(key)
        for r,w in items.items():
            t=representative(('A' if r<=8 else 'D',r))
            assert w==base(t)+sum(upper_cross(t,u) for u in large)
            counts['literal_interactions']+=1
        rem=N-row['used_rank'];value=sum(base(t) for t in large)+table['values'][rem]
        forest=list(large)
        while rem:
            r=table['last'][rem];forest.append(representative(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==N and forest_cost(forest,True)==value
        assert forest_cost(forest)<=value
        alternatives.append(dict(counts=row['counts'],doubled_upper=value,
                                 formal_upper_attainer=sorted(k+str(r) for k,r in forest),
                                 actual_doubled_at_that_forest=forest_cost(forest)))
        counts['global_alternatives']+=1
    maximum=max(r['doubled_upper'] for r in alternatives);assert maximum==898
    witness=[('D',17),('D',17),('D',29)]
    assert sum(r for k,r in witness)==N and forest_cost(witness)==forest_cost(witness,True)==898
    assert base(('D',29))==398
    actual_coupon_rows=[dict(types=[list(a),list(b)],added_component_pair_cost=str(F(v,2)))
                        for (a,b),v in NEW_CROSS.items()]
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in STEMS
          for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_RANK63_SPARSE_COUPON_FOREST_BUDGET',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                ambient_dimension=63,finite_maximum='449',finite_curvature='899/2',finite_slack='1/2',
                actual_attaining_forest=['D17','D17','D29'],actual_attainer_cost='449',
                actual_cross_coupon_lists=actual_coupon_rows,
                representative_allowance=dict(types=['E6','D11'],actual_total='2',upper_total='4',
                                              purpose='Dominates actual E6/A11 total4 after replacing A11 by D11; not an added actual coupon.'),
                domination_mode='Simultaneous pointwise replacement of both component types; compare upper_cross(Ra,Rb) against actual_cross(a,b).',
                type_domination=domination,type_pair_checks=pair_checks,
                complete_large_vectors=len(dom),complete_keys=len(keys),large_domain=dom,
                tables=[tables[k] for k in keys],alternatives=alternatives,checks=dict(counts),
                new_dimension_bound=None,
                scope='Complete actual rank63 finite-rule maximum449: every forest is bounded by the conservative representative envelope, whose maximum449 is attained by an actual-cost typed forest. All new coupons counted, B9/E6/A17 retained distinctly. No geometric or ambient theorem supplied by this counting certificate.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    out=json.loads(json.dumps(run(saved)))
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['finite_maximum'],out['checks'],flush=True)
