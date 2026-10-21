#!/usr/bin/env python3
"""Cycle14 distance6/7 transfer using the strict spectator-energy spectrum."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import sympy as s
import fresh_finite8_cycle_replacement_20261017 as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_cycle14_ideal5_transfer_20261017'
A,B,C=R.A,R.B,R.C

def run():
    rows=[];single=[];special=[]
    ps=B.cycle_profiles()
    for P,Q in((6,8),(7,7)):
        for p in ps:
            if p['name']in('isolated','both_ends3'):continue
            side=next(i for i,x in enumerate(p['position'])if x)
            pos=p['position'][side];contact=pos if side==0 else 14-pos
            w2=2 if p['name'].startswith('end4')else 1
            edges=[(i,i+1,1)for i in range(1,13)]+[(14,contact,w2)]
            ids=list(range(1,15));adj={i:set()for i in ids}
            for i,j,_ in edges:adj[i].add(j);adj[j].add(i)
            leaf=13 if side==0 else 1;near=next(iter(adj[leaf]))
            assert len(adj[leaf])==1 and len(adj[near])==2
            dist=C.distances(adj,ids,14)[P];assert dist<=8
            single.append(dict(P=P,Q=Q,profile=p['name'],marked_distance=dist,
                fixed_leaf=leaf,degree_two_neighbor=near,positive_pivots=C.tree_pd(ids,edges)))
        keep=[p for p in ps if p['name']in('isolated','both_ends3')]
        for c,d in product(keep,repeat=2):
            vals=[]
            for ub,vb in product((0,1),repeat=2):
                val=A.value(A.Q(F(1,P)),A.Q(F(1,Q)),c,d,c['umax']if ub else A.Q(),d['umax']if vb else A.Q())
                vals.append(val)
                rows.append(dict(P=P,Q=Q,c=c['name'],d=d['name'],upper_reciprocals=[ub,vb],value=val.wire()))
            low=R.minimum(vals)
            if c['name']!=d['name']or c['name']!='isolated':assert low.sign()>0
        alpha=F(1,P)+F(1,Q);gate=1-3*alpha;required=4-1/alpha
        assert gate>0
        special.append(dict(P=P,Q=Q,one_high_height_gap=str(gate),forced_heights=['1','1'],required_D=str(required)))
    assert len(rows)==32 and len(single)==12
    assert [r['required_D']for r in special]==['4/7','1/2']
    spectrum=json.loads((ROOT/'results/fresh_small_spectator_difference_gap_20261017.json').read_text())
    assert spectrum['status']=='PASS_DISCRETE_SPECTATOR_ENERGY_AT_MOST_THREE_QUARTERS'
    assert all(s.simplify(s.Rational(4,7)-s.sympify(x))!=0 for x in spectrum['possible_energies'])
    # Reconstruct the sealed P7/Q7 equality, not the old complete cone census.
    n=17;g=[[F(2 if i==j else 0)for j in range(n)]for i in range(n)]
    for i,j in[(i,(i+1)%14)for i in range(14)]+[(0,14),(7,15),(14,16)]:g[i][j]=g[j][i]=F(-1)
    W=[i for i in range(14)if i not in(0,7)]+[16]
    assert C.inertia(g)==[15,1,1]
    assert C.inertia([[g[i][j]for j in W]for i in W])==[13,0,0]
    corners=[]
    for p,q in((0,7),(7,14),(14,15),(15,0)):
        ids=sorted(set(W)|{p,q});sig=C.inertia([[g[i][j]for j in ids]for i in ids])
        if(p,q)==(0,7):
            assert sig==[14,0,1]
            corners.append(dict(marks=[p,q],kind='ideal',coupon='0'));continue
        assert sig==[15,0,0]
        adj={i:{j for j in ids if i!=j and g[i][j]}for i in ids}
        parts=[]
        for mark in(p,q):
            ds=C.distances(adj,ids,mark);leaves=[i for i in ds if len(adj[i]&ds.keys())==1]
            assert all(len(adj[i]&ds.keys())<=2 for i in ds)
            parts.append(dict(mark=mark,nodes=sorted(ds),type=f'A{len(ds)}',endpoint_distances=sorted(ds[i]for i in leaves)))
        assert not(set(parts[0]['nodes'])&set(parts[1]['nodes']))
        sizes=sorted(len(t['nodes'])for t in parts);coupon=F()
        if sizes in([1,13],[2,13]):
            assert next(t for t in parts if len(t['nodes'])==13)['endpoint_distances']==[6,6]
            coupon=F(1,2)
        corners.append(dict(marks=[p,q],kind='finite',signature=sig,components=parts,coupon=str(coupon)))
    assert sum(F(r['coupon'])for r in corners)==1
    paths=[Path(__file__)]
    for stem in('fresh_small_spectator_difference_gap_20261017','fresh_finite8_ideal5_local_replacement_20261017',
                'fresh_ideal_distance6_mixed_replacement_20261017'):
        paths += [ROOT/f'{folder}/{stem}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    retained=[[i,j]for i in range(14)for j in range(i+1,14)if min(j-i,14-j+i)<=5]
    removed=[[i,j]for i in range(14)for j in range(i+1,14)if min(j-i,14-j+i)>5]
    assert len(retained)==70 and len(removed)==21
    return dict(status='PASS_UNIVERSAL_CYCLE14_IDEAL5_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        single_profiles=single,rectangles=rows,isolated_cases=special,
        spectrum=spectrum['possible_energies'],equality=dict(signature=[15,1,1],white=W,corners=corners),
        cycle_pair_rule=dict(retained=retained,transferred=removed,cost=70),
        counts=dict(single_profiles=12,rectangle_values=32,isolated_cases=2,finite_equality_corners=3),
        theorem='An actual ideal cycle14 pair at distance6 or7 transfers to finite distance<=8, ideal short paths, or two sealed whole-small A1/A2 x whole-A13-middle half-coupons. No simple-cusp/arithmetic assumption.',
        scope='Local transfer only; no complete high-core overlap, small-factor coverage, total budget or dimension claim. Old equality-core realization/cone census is inherited and not rerun.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
