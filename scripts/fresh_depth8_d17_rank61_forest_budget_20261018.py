#!/usr/bin/env python3
"""Complete rank61 finite budget for depth8 zeros and finite D17 cost118.

Generation computes new residual tables by item-first dynamic programming.
Replay checks all last-item upper inequalities and attaining predecessors,
without rerunning the optimizing recurrence. No ambient bound follows.
"""
import argparse,ast
from collections import Counter,deque
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_depth8_d17_rank61_forest_budget_20261018'
DEPTH='fresh_depth8_finite8_zero_transfer_20261018'
D17='fresh_finite_d17_a4_d13_half_transfer_20261018'
OLD='fresh_finite_rank15_rank16_rank62_budget_20261018'
FRONT='fresh_ambient61_payment_frontier_20261018'
CURRENT='fresh_rank15_rank16_dimension61_20261018'
STEMS=(DEPTH,D17,OLD,FRONT,CURRENT)
N=61


def base(t,upper=False):
    k,r=t
    if r==15:return 208
    if r==16:return 208 if k=='B' else 226 if k=='A' or upper else 225
    if r==17:return {'A':230,'B':210,'D':236}[k]
    if r<=14:return r*(r-1)
    return 16*r-(66 if k=='D' and r>=20 else 64 if k=='D' else 68)-2*max(0,r-{'A':24,'B':25,'D':26}[k])


SPARSE={
    (('A',4),('A',12)):8,(('A',4),('A',13)):8,
    (('D',5),('A',12)):4,(('A',2),('A',13)):4,
    (('A',2),('A',14)):4,(('E',6),('A',11)):8,
    (('A',2),('B',9)):2,(('A',2),('D',10)):4,
    (('A',3),('D',10)):2,
}


def cross(a,b):
    value=0
    for sm,la in ((a,b),(b,a)):
        if sm[1]<=5:
            if la==('D',10):value+=2*sm[1]
            elif la[0]=='D' and 11<=la[1]<=16:value+=3*sm[1]
            elif la in (('A',13),('A',15)):value+=sm[1]
            elif la==('A',17):value+=2*sm[1]
        if (sm,la) in ((('A',1),('D',11)),(('A',1),('D',12)),
            (('A',2),('B',10)),(('A',4),('D',13)),
            (('A',1),('A',11)),(('A',2),('A',11))):value+=2
        if (sm,la)==(('A',2),('D',11)):value+=4
        if (sm,la)==(('E',6),('D',11)):value+=8
        if (sm,la)==(('A',1),('B',10)):value+=1
        if sm in (('A',1),('A',2)) and la==('D',14):value+=sm[1]
        value+=SPARSE.get((sm,la),0)
    return value


def representative(t):
    if t==('A',17):return t
    if t[1]==6:return ('E',6)
    if t[1]==9:return ('B',9)
    return ('A' if t[1]<=8 else 'D',t[1])


def types():
    return ([('A',r) for r in range(1,N+1)]+[('B',r) for r in range(2,N+1)]+
        [('D',r) for r in range(4,N+1)]+[('E',6),('E',7),('E',8),('F',4),('H',3),('H',4),('I2*',2)])


def cost(ts,upper=False):
    return sum(base(t,upper) for t in ts)+sum(cross(a,b) for i,a in enumerate(ts) for b in ts[i+1:])


def domain():
    rows=[]
    def rec(r,left,counts):
        if r==18:
            c=dict(zip(range(9,18),counts))
            key=(2*c[10]+3*sum(c[s] for s in range(11,17))+2*c[17],
                 c[11]+c[12],c[13],c[11],c[14],c[9],c[10])
            rows.append((tuple(counts),N-left,key));return
        for n in range(left//r+1):rec(r+1,left-r*n,counts+[n])
    rec(9,N,[]);return sorted(rows)


def item(rank,key):
    broad,extra,a4,d11,d14,b9,d10=key
    t=representative(('A' if rank<=8 else 'D',rank))
    out=base(t,True)+(rank*broad if rank<=5 else 0)
    out+=(2*extra+d14 if rank==1 else 4*d11+2*d14+2*b9+4*d10 if rank==2 else
          2*d10 if rank==3 else 2*a4 if rank==4 else 8*d11 if rank==6 else 0)
    return out


ITEMS=tuple(range(1,9))+tuple(range(17,N+1))


def generate_table(key):
    # Unbounded item-first generation; the checker below uses every possible
    # final item, independently of this loop order.
    values=[0]+[-10**9]*N
    for r in ITEMS:
        w=item(r,key)
        for s in range(r,N+1):values[s]=max(values[s],values[s-r]+w)
    last=[0]+[next(r for r in ITEMS if r<=s and values[s]==values[s-r]+item(r,key)) for s in range(1,N+1)]
    return dict(key=list(key),values=values,last=last)


def check_table(tab,key):
    assert tab['key']==list(key)
    values,last=tab['values'],tab['last']
    assert len(values)==len(last)==N+1 and values[0]==last[0]==0
    checks=0
    for s in range(1,N+1):
        for r in ITEMS:
            if r<=s:assert values[s]>=values[s-r]+item(r,key);checks+=1
        r=last[s];assert r in ITEMS and r<=s and values[s]==values[s-r]+item(r,key)
    return checks


def deep_count(kind,n,depth):
    adj=[[] for _ in range(n)]
    edges=[(i,i+1) for i in range(n-1)] if kind!='D' else [(i,i+1) for i in range(n-3)]+[(n-3,n-2),(n-3,n-1)]
    for a,b in edges:adj[a].append(b);adj[b].append(a)
    features={i for i in range(n) if len(adj[i])!=2}
    if kind=='B':features|={n-2,n-1}
    ds=[]
    for i in range(n):
        dist={i:0};q=deque([i])
        while q:
            a=q.popleft()
            for b in adj[a]:
                if b not in dist:dist[b]=dist[a]+1;q.append(b)
        assert len(dist)==n;ds.append(dist)
    out=sum(ds[i][j]==8 and min(ds[i][f] for f in features)>=depth and
            min(ds[j][f] for f in features)>=depth for i in range(n) for j in range(i+1,n))
    assert out==max(0,n-(2*depth+{'A':8,'B':9,'D':10}[kind]))
    return out


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[CURRENT]['dimension_bound']==61
    assert docs[D17]['new_cost']=='118' and docs[D17]['new_cross_component_coupons']==[]
    assert docs[D17]['all_other_pair_weights_unchanged']
    assert docs[OLD]['finite_maximum']=='434' and docs[OLD]['geometric_inputs_final']
    assert docs[OLD]['actual_costs']==dict(A15='104',B15='104',D15='104',A16='113',D16='225/2',A17='115',D17='120')
    assert docs[OLD]['new_cross_component_coupons']==[]
    old_ast=ast.parse((ROOT/f'scripts/{OLD}.py').read_text())
    old_sparse=next(ast.literal_eval(n.value) for n in old_ast.body if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='SPARSE_CROSS' for t in n.targets))
    assert old_sparse==SPARSE
    assert docs[FRONT]['finite']['exact_maximum']=='427'
    link_dimension=N-1;vertices=N
    G=sum(min(j,link_dimension-j) for j in range(link_dimension+1))
    curvature=F(5*G-link_dimension*vertices,2)
    assert (G,curvature)==(900,F(420))
    assert docs[FRONT]['finite']['curvature']==str(curvature)
    depth=docs[DEPTH]
    assert depth['retained_contact']['feature_depth']==8
    assert depth['current_baseline']==CURRENT and depth['previous_feature_depth']==11
    assert depth['feature_depth']==depth['finite_distance']==8 and depth['new_weight']=='0'
    assert depth['ideal_weights_unchanged'] and depth['all_other_finite_weights_and_coupons_unchanged']
    dom=domain();assert len(dom)==1132
    assert {x[0] for x in dom}=={tuple(r['counts']) for r in docs[FRONT]['finite']['alternatives']}
    keys=sorted({r[2] for r in dom});assert len(keys)==854
    ts=types()
    for a in ts:
        assert base(representative(a),True)>=base(a)
        for b in ts:
            assert cross(a,b)==cross(b,a)
            assert cross(representative(a),representative(b))>=cross(a,b)
    counts=[]
    for kind in 'ABD':
        for n in range(18,N+1):
            new=deep_count(kind,n,8);old=deep_count(kind,n,11)
            counts.append(dict(type=kind+str(n),old_depth11=old,new_depth8=new,additional_zeros=new-old))
    tables={tuple(t['key']):t for t in saved['tables']} if saved else {k:generate_table(k) for k in keys}
    assert set(tables)==set(keys) and (saved is None or len(saved['tables'])==len(keys))
    inequalities=sum(check_table(tables[k],k) for k in keys)
    records=[];item_checks=0
    for cs,used,key in dom:
        large=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for k,r in large)==used
        assert all(cross(a,b)==0 for i,a in enumerate(large) for b in large[i+1:])
        for r in ITEMS:
            t=representative(('A' if r<=8 else 'D',r))
            assert item(r,key)==base(t,True)+sum(cross(t,b) for b in large);item_checks+=1
        rem=N-used;tab=tables[key];total=sum(base(t,True) for t in large)+tab['values'][rem]
        forest=large[:]
        while rem:
            r=tab['last'][rem];forest.append(representative(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==N and cost(forest,True)==total
        actual=cost(forest);assert actual<=total
        records.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
            doubled_upper=total,doubled_actual=actual))
    maximum=max(r['doubled_upper'] for r in records);assert maximum==850
    witness=[('D',15)]*3+[('A',16)]
    assert sum(r for k,r in witness)==N and cost(witness)==maximum
    assert base(('D',61))==840
    sources={f'scripts/{TAG}.py':sha256(Path(__file__).read_bytes()).hexdigest()}
    for s in STEMS:
        for folder,ext in (('scripts','py'),('results','json'),('results','md')):
            p=f'{folder}/{s}.{ext}';sources[p]=sha256((ROOT/p).read_bytes()).hexdigest()
    if saved:assert saved['sources']==sources
    return dict(status='PASS_DEPTH8_D17_RANK61_COMPLETE_FINITE_BUDGET',sources=sources,
        ambient_dimension=N,finite_maximum=str(F(maximum,2)),finite_curvature=str(curvature),
        finite_slack=str(curvature-F(maximum,2)),finite_link_G=G,finite_link_vertices=vertices,
        actual_attaining_forests=[['A16','D15','D15','D15']],singleton_D61_cost='420',
        proved_changes=dict(feature_depth_old=11,feature_depth_new=8,D17_old='120',D17_new='118',new_coupons=[]),
        current_component_costs=dict(A15='104',B15='104',D15='104',A16='113',D16='225/2',A17='115',D17='118'),
        representative_D16_upper='113',simultaneous_type_pairs=len(ts)**2,deep_pair_counts=counts,
        complete_large_vectors=len(dom),interaction_keys=len(keys),tables=[tables[k] for k in keys],
        Bellman_upper_inequalities=inequalities,attaining_predecessors=len(keys)*N,
        residual_item_checks=item_checks,alternatives=records,
        representative_alternatives_above_curvature=sorted((r for r in records if r['doubled_upper']>840),
            key=lambda r:(-r['doubled_upper'],r['forest'])),
        new_residual_tables=True,generation_algorithm='item-first unbounded exact dynamic program',
        replay_algorithm='every last-item upper inequality and a saved attaining predecessor',
        inherited_old_tables_reused=False,new_dimension_bound=None,
        scope='Complete finite rank61 counting theorem for the simultaneous current61 rule with proveddepth8 zeros and D17=118 only. All current coupons and A17 restoration retained; conservative D16upper113 distinguished from actual225/2. Maximum425 is attained by actual formal3D15+A16, with no physical realization claim. It exceeds curvature420, so no ambient60 theorem follows.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';saved=json.loads(p.read_text()) if args.verify else None
    result=run(saved)
    if args.verify:assert result==saved
    else:
        assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'maximum',result['finite_maximum'],'vectors',result['complete_large_vectors'],
          'tables',result['interaction_keys'],'above420',len(result['representative_alternatives_above_curvature']),flush=True)
