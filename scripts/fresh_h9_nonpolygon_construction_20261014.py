#!/usr/bin/env python3
"""Actual H9 multiwedge incidences from the simple dodecahedron.

Exact mandatory Lanner4/5 joins precede any metric search. The source
polytope is explicitly certified; unknown meeting labels are relaxed only
for necessary finite-vertex tests, as in the previously audited helper.
"""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
from fractions import Fraction as F
import argparse,hashlib,json,math,time
import networkx as nx
from fresh_h9_nonproduct_incidence_20261014_multiwedge import join

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_nonpolygon_construction_20261014'
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def mul(a,b):return(a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])
def dot(a,b):
    z=(0,0)
    for x,y in zip(a,b,strict=True):z=add(z,mul(x,y))
    return z
def sign(a):
    # a0+a1*phi=(2a0+a1+a1*sqrt5)/2.
    x,y=2*a[0]+a[1],a[1]
    s=lambda z:int(z>0)-int(z<0)
    if not y:return s(x)
    if not x:return s(y)
    if s(x)==s(y):return s(x)
    return s(x)*s(x*x-5*y*y)
def neg(a):return(-a[0],-a[1])
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def base():
    V=[]
    for zero in range(3):
        for a,b in product([-1,1],repeat=2):
            v=[(0,0)]*3;v[(zero+1)%3]=(a,0);v[(zero+2)%3]=(0,b);V.append(v)
    G=[[dot(v,w)for w in V]for v in V]
    assert all(G[i][i]==(2,1)for i in range(12))
    assert all(sign(add(G[i][j],neg(G[i][i])))<0 for i in range(12)for j in range(12)if i!=j)
    edges=[(i,j)for i,j in combinations(range(12),2)if G[i][j]==(0,1)]
    facets=[t for t in combinations(range(12),3)if all(tuple(sorted(p))in edges for p in combinations(t,2))]
    assert len(edges)==30 and len(facets)==20
    for t in facets:
        n=[tuple(sum(V[i][d][k]for i in t)for k in range(2))for d in range(3)]
        values=[dot(n,v)for v in V];assert all(values[i]==(2,3)for i in t)
        assert all(sign(add(values[i],(-2,-3)))<0 for i in range(12)if i not in t)
    assert Counter(e for t in facets for e in combinations(t,2))==Counter({e:2 for e in edges})
    # Each vertex is extreme by the strict norm support above. Twenty
    # triangular supporting facets attain the 3D bound 2*12-4, hence are
    # the complete boundary. Central symmetry puts zero in the interior.
    assert all(any(w==[neg(x)for x in v]for w in V)for v in V)
    graph=nx.Graph();graph.add_nodes_from(range(12));graph.add_edges_from(edges)
    group=sorted({tuple(p[i]for i in range(12))for p in nx.algorithms.isomorphism.GraphMatcher(graph,graph).isomorphisms_iter()})
    assert len(group)==120
    FS=set(facets)
    assert all({tuple(sorted(p[i]for i in t))for t in facets}==FS for p in group)
    GS=set(group);assert tuple(range(12))in GS
    assert all(tuple(p[q[i]]for i in range(12))in GS for p in group for q in group)
    # A triangle flag determines an automorphism: across any shared edge,
    # its other incident triangle is forced. Thus 20*6=120 is also an
    # upper bound, certifying that this is the full base symmetry group.
    reached={facets[0]};todo=list(reached)
    for t in todo:
        for u in facets:
            if len(set(t)&set(u))==2 and u not in reached:reached.add(u);todo.append(u)
    assert len(reached)==20
    nonfaces=[(i,j)for i,j in combinations(range(12),2)if(i,j)not in edges]
    assert len(nonfaces)==36
    assert {t for t in combinations(range(12),3)if not any(set(m)<=set(t)for m in nonfaces)}==FS
    assert not any(all(tuple(sorted(p))in edges for p in combinations(t,2))for t in combinations(range(12),4))
    return {'base':'polar of the displayed icosahedron; simple dodecahedron','dimension':3,
        'vertices_of_dual':V,'dual_facets':facets,'minimal_nonfaces':nonfaces,'symmetry_group':group},facets,nonfaces,group

def weights(nonfaces,group):
    def compositions(total,n,prefix=()):
        if n==1:yield prefix+(total,);return
        for a in range(1,total-n+2):yield from compositions(total-a,n-1,prefix+(a,))
    vectors=list(compositions(18,12));assert len(vectors)==math.comb(17,11)==12376
    valid=[w for w in vectors if all(sum(w[i]for i in m)<=5 for m in nonfaces)]
    reps={min(tuple(w[p[i]]for i in range(12))for p in group)for w in valid}
    # Favor distributed weights and fewer large mandatory diagrams: these
    # are constructive search priorities, not completeness filters.
    out=sorted(reps,key=lambda w:(max(w),sum(sum(w[i]for i in m)==5 for m in nonfaces),sum(sum(w[i]for i in m)==4 for m in nonfaces),w))
    return out,{'all_positive_weight_compositions':len(vectors),'missing_face_order_bound_passes':len(valid),'symmetry_types':len(out)}

def incidence(w,facets,missing):
    blocks=[];n=0
    for k in w:blocks.append(tuple(range(n,n+k)));n+=k
    M=sorted(tuple(sorted(i for j in m for i in blocks[j]))for m in missing)
    V=set()
    for f in facets:
        outside=[blocks[j]for j in range(12)if j not in f]
        for omit in product(*outside):V.add(tuple(i for i in range(n)if i not in omit))
    assert n==18 and all(len(v)==9 for v in V)
    masks=[sum(1<<i for i in m)for m in M]
    assert all(not any((sum(1<<i for i in v)&m)==m for m in masks)for v in V)
    assert set(Counter(e for v in V for e in combinations(v,8)).values())=={2}
    # Recovering blocks from equal minimal-nonface membership gives the
    # icosahedral base; this cannot be an old polygon multiwedge or a product.
    membership={}
    for i in range(n):membership.setdefault(tuple(j for j,m in enumerate(M)if i in m),[]).append(i)
    assert sorted(map(tuple,membership.values()))==sorted(blocks)
    assert any(set(a)&set(b)for a,b in combinations(M,2))
    return blocks,M,sorted(V)

def run(limit,maxstates):
    B,Faces,M0,Group=base();allw,counts=weights(M0,Group);rows=[];start=time.monotonic()
    print('base types',counts,flush=True)
    for case,w in enumerate(allw[:limit]):
        blocks,M,V=incidence(w,Faces,M0)
        print('case',case,'weights',w,'vertices',len(V),'missing',dict(Counter(map(len,M))),flush=True)
        trial=join(M,V,maxstates)
        rows.append({'case':case,'weights':w,'blocks':blocks,'minimal_nonfaces':M,'actual_vertices':V,**trial})
        files=[f'scripts/{TAG}.py','scripts/fresh_h9_nonproduct_incidence_20261014_multiwedge.py','scripts/independent_h9_14_record30_a1_labels.py']
        data={'status':'EXACT_DODECAHEDRAL_MULTIWEDGE_CONSTRUCTIVE_LABEL_CHECKPOINT','sources':{p:sha(p)for p in files},
            'settings':{'limit':limit,'max_states':maxstates},'base':B,'weight_counts':counts,'available_weights':allw,
            'cases':rows,'all_weight_types_processed':len(rows)==len(allw),
            'scope':'Actual eighteen-facet H9 incidences by six ordinary wedges over a simple dodecahedron. Exact mandatory Lanner4/5 joins and necessary finite-vertex restrictions. Omitted constraints are relaxed; an empty untruncated join excludes only its displayed type. Positive states require full labels, metric rank/signature, and compactness checks.'}
    return json.loads(json.dumps(data))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=30);p.add_argument('--max-states',type=int,default=10000);p.add_argument('--verify',action='store_true');a=p.parse_args()
    assert a.limit>0 and a.max_states>0;z=run(a.limit,a.max_states);out=ROOT/f'results/{TAG}.json'
    if a.verify:assert z==json.loads(out.read_text())
    else:out.write_text(json.dumps(z,indent=2)+'\n')
    print('FINAL',len(z['cases']),z['weight_counts']['symmetry_types'],Counter(r['verdict']for r in z['cases']),flush=True)
