#!/usr/bin/env python3
"""Vertex-truncated dodecahedral-prism bases for actual H9 multiwedges."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
from collections import Counter
import argparse,hashlib,json
import fresh_h9_nonpolygon_construction_20261014 as B
from fresh_h9_nonproduct_incidence_20261014_multiwedge import join
ROOT=B.ROOT;TAG='fresh_h9_nonpolygon_construction_20261014_fourbase'
Z=(F(0),F(0));O=(F(1),F(0))
def scale(x,c):return(x[0]*c,x[1]*c)
def inverse(x):
    a,b=x;n=a*a+a*b-b*b;assert n
    return((a+b)/n,-b/n)
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return[enc(v)for v in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def make(target):
    base,tri,_,_=B.base();V=[[tuple(map(F,a))for a in v]+[Z]for v in base['vertices_of_dual']]+[[Z,Z,Z,scale(O,s)]for s in[-1,1]]
    facets=sorted(tuple(sorted(t+(p,)))for t in tri for p in[12,13]);N=[]
    for f in facets:
        t=[i for i in f if i<12];p=next(i for i in f if i>=12)
        n=[B.mul(tuple(sum(V[i][d][k]for i in t)for k in range(2)),inverse((F(2),F(3))))for d in range(3)]+[scale(O,-1 if p==12 else 1)]
        vals=[B.dot(n,v)for v in V]
        assert all(v==O if i in f else B.sign(B.add(v,(-1,0)))<0 for i,v in enumerate(vals));N.append(n)
    cut=tuple(sorted(target+(12,)));assert cut in facets
    p=[scale(tuple(sum(V[i][d][k]for i in cut)for k in range(2)),F(101,400))for d in range(4)]
    vals=[B.dot(n,p)for n in N]
    assert all(B.sign(B.add(v,(-1,0)))==(1 if f==cut else-1)for f,v in zip(facets,vals,strict=True))
    # A point beyond exactly one simplicial facet stacks over that facet;
    # dualizing is an actual small vertex truncation of the simple prism.
    V.append(p);facets=[f for f in facets if f!=cut]+[tuple(sorted(t+(14,)))for t in combinations(cut,3)];facets=sorted(facets)
    assert len(facets)==43 and all(len(f)==4 for f in facets)
    assert set(Counter(t for f in facets for t in combinations(f,3)).values())=={2}
    face_sets=[set(f)for f in facets]
    def face(s):return any(set(s)<=f for f in face_sets)
    missing=[]
    for k in range(2,6):
        for s in combinations(range(15),k):
            if not face(s)and all(face(t)for t in combinations(s,k-1)):missing.append(s)
    chosen=tri[0];w=[2 if i in chosen or i in[12,13] else 1 for i in range(15)];assert sum(w)==20
    blocks=[];n=0
    for k in w:blocks.append(tuple(range(n,n+k)));n+=k
    M=sorted(tuple(sorted(i for j in m for i in blocks[j]))for m in missing)
    assert max(map(len,M))<=5
    Vert=set()
    for f in facets:
        for omit in product(*(blocks[j]for j in range(15)if j not in f)):
            Vert.add(tuple(i for i in range(n)if i not in omit))
    assert all(len(v)==9 and not any(set(m)<=set(v)for m in M)for v in Vert)
    assert set(Counter(t for v in Vert for t in combinations(v,8)).values())=={2}
    twin={}
    for i in range(n):twin.setdefault(tuple(j for j,m in enumerate(M)if i in m),[]).append(i)
    assert sorted(map(tuple,twin.values()))==sorted(blocks)
    # Minimal nonfaces involving the new facet link the former product
    # factors. Connectivity of this hypergraph rules out a product split.
    adj=[set()for _ in range(n)]
    for m in M:
        for i in m:adj[i].update(set(m)-{i})
    seen={0};todo=[0]
    for i in todo:
        for j in adj[i]-seen:seen.add(j);todo.append(j)
    assert len(seen)==n
    return {'old_weighted_triangle':chosen,'truncated_vertex_facets':cut,'base_dual_vertices':V,
        'old_supporting_normals':N,'new_point_old_support_values':vals,'base_dual_facets':facets,
        'base_minimal_nonfaces':missing,'weights':w,'blocks':blocks,'minimal_nonfaces':M,'actual_vertices':sorted(Vert)}
def run(maxstates):
    _,tri,_,group=B.base();chosen=tri[0];targets=[t for t in tri if not set(t)&set(chosen)]
    stabilizer=[p for p in group if set(p[i]for i in chosen)==set(chosen)]
    reps=sorted({min(tuple(sorted(p[i]for i in t))for p in stabilizer)for t in targets});rows=[]
    print('target vertices',len(targets),'symmetry recipes',len(reps),flush=True)
    for i,t in enumerate(reps):
        r=make(t);print('case',i,'vertices',len(r['actual_vertices']),'missing',dict(Counter(map(len,r['minimal_nonfaces']))),flush=True)
        rows.append({'case':i,**r,**join(r['minimal_nonfaces'],r['actual_vertices'],maxstates)})
    files=[f'scripts/{TAG}.py','scripts/fresh_h9_nonpolygon_construction_20261014.py','scripts/fresh_h9_nonproduct_incidence_20261014_multiwedge.py','scripts/independent_h9_14_record30_a1_labels.py']
    return enc({'status':'EXACT_VERTEX_TRUNCATED_PRISM_H9_MULTIWEDGE_CHECKPOINT','sources':{p:sha(p)for p in files},
        'max_states':maxstates,'eligible_base_vertices':targets,'symmetry_representatives':reps,'cases':rows,
        'scope':'Five wedges over the three facets meeting at one dodecahedral vertex and both interval facets, after truncating a disjoint dodecahedral vertex on one interval end. All target-vertex symmetry classes in this recipe are tested. Other bases and truncations remain outside scope.'})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--max-states',type=int,default=10000);a=p.parse_args();z=run(a.max_states);out=ROOT/f'results/{TAG}.json'
    if a.verify:assert z==json.loads(out.read_text())
    else:out.write_text(json.dumps(z,indent=2)+'\n')
    print('FINAL',Counter(r['verdict']for r in z['cases']),flush=True)
