#!/usr/bin/env python3
"""Constructive edge-truncated-cube multiwedges, with exact label gates.

One lexicographic truncation order for each six-original-edge subset orbit
is tested. This is a bounded constructive family, not all 3D polytopes.
"""
from pathlib import Path
from itertools import combinations,product,permutations
from fractions import Fraction as F
from collections import Counter
import argparse,hashlib,json
from fresh_h9_nonproduct_incidence_20261014_multiwedge import join
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_nonpolygon_construction_20261014_edgecut'
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def sub(a,b):return tuple(x-y for x,y in zip(a,b,strict=True))
def cross(a,b):return(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def plane(V,f):
    a,b,c=[V[i]for i in f];n=cross(sub(b,a),sub(c,a));d=dot(n,a);assert d
    if d<0:n=tuple(-x for x in n);d=-d
    assert all(dot(n,v)==d if i in f else dot(n,v)<d for i,v in enumerate(V))
    return n,d
def seed():
    V=[]
    for d in range(3):
        for s in[-1,1]:V.append(tuple(F(s if j==d else 0)for j in range(3)))
    Faces=sorted(tuple(t)for t in product(*[(0,1),(2,3),(4,5)]))
    Edges=sorted({e for f in Faces for e in combinations(f,2)})
    group=[]
    for p in permutations(range(3)):
        for flips in product([0,1],repeat=3):group.append(tuple(2*p[i//2]+((i%2)^flips[i//2])for i in range(6)))
    assert len(set(group))==48
    assert all({tuple(sorted(p[i]for i in f))for f in Faces}==set(Faces)for p in group)
    reps={min(tuple(sorted(tuple(sorted(p[i]for i in e))for e in ss))for p in group)for ss in combinations(Edges,6)}
    return V,Faces,sorted(reps)
def truncated_base(order):
    V,Faces,_=seed();steps=[]
    for edge in order:
        adjacent=[f for f in Faces if set(edge)<=set(f)];assert len(adjacent)==2
        midpoint=tuple((V[edge[0]][j]+V[edge[1]][j])/2 for j in range(3))
        planes={f:plane(V,f)for f in Faces};upper=[]
        for f,(n,d)in planes.items():
            a=dot(n,midpoint)
            if f in adjacent:assert a==d
            elif a>0:upper.append(d/a)
        bound=min(upper)if upper else F(3);assert bound>1
        scale=(1+bound)/2;v=tuple(scale*x for x in midpoint);new=len(V)
        assert all(dot(n,v)>d if f in adjacent else dot(n,v)<d for f,(n,d)in planes.items())
        V.append(v);Faces=[f for f in Faces if f not in adjacent]+[tuple(sorted((i,next(k for k in f if k not in edge),new)))for f in adjacent for i in edge]
        Faces=sorted(Faces);assert len(Faces)==2*len(V)-4
        for f in Faces:plane(V,f)
        steps.append({'edge':edge,'new_vertex':v,'radial_scale':scale,'removed_faces':adjacent})
    assert len(V)==12 and len(Faces)==20
    E={p for f in Faces for p in combinations(f,2)}
    M=[p for p in combinations(range(12),2)if p not in E]
    assert len(E)==30 and len(M)==36
    assert {f for f in combinations(range(12),3)if all(p in E for p in combinations(f,2))}==set(Faces)
    assert not any(all(p in E for p in combinations(f,2))for f in combinations(range(12),4))
    return V,Faces,M,steps
def incidences(Faces,M0):
    w=[2]*6+[1]*6;blocks=[];n=0
    for k in w:blocks.append(tuple(range(n,n+k)));n+=k
    M=sorted(tuple(sorted(i for j in m for i in blocks[j]))for m in M0)
    V=set()
    for f in Faces:
        for omit in product(*(blocks[j]for j in range(12)if j not in f)):
            V.add(tuple(i for i in range(n)if i not in omit))
    assert all(len(v)==9 and not any(set(m)<=set(v)for m in M)for v in V)
    assert set(Counter(e for v in V for e in combinations(v,8)).values())=={2}
    twin={}
    for i in range(n):twin.setdefault(tuple(j for j,m in enumerate(M)if i in m),[]).append(i)
    return blocks,M,sorted(V),sorted(twin.values())
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return[enc(v)for v in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def run(maxstates):
    _,_,orders=seed();rows=[];print('six-original-edge subset orbit recipes',len(orders),flush=True)
    for ix,order in enumerate(orders):
        v,f,m,steps=truncated_base(order);blocks,M,V,twins=incidences(f,m)
        print('case',ix,'baseclasses',len(twins),'vertices',len(V),flush=True)
        if len(twins)==12:r=join(M,V,maxstates)
        else:r={'verdict':'LOWER_DIMENSIONAL_MULTIWEDGE_BASE_NOT_PURSUED','terminal_states':[]}
        rows.append({'case':ix,'old_edge_truncation_order':order,'base_vertices':v,'base_facets':f,
            'exact_truncation_steps':steps,'weights':[2]*6+[1]*6,'blocks':blocks,
            'minimal_nonfaces':M,'actual_vertices':V,'recovered_blocks':twins,**r})
    sources=[f'scripts/{TAG}.py','scripts/fresh_h9_nonproduct_incidence_20261014_multiwedge.py','scripts/independent_h9_14_record30_a1_labels.py']
    return enc({'status':'EXACT_EDGE_TRUNCATED_CUBE_MULTIWEDGE_CHECKPOINT','sources':{p:sha(p)for p in sources},
        'max_states':maxstates,'all_six_old_edge_subsets':924,'subset_orbit_recipes':len(orders),'cases':rows,
        'scope':'One lexicographic truncation order for every six-original-edge subset orbit of the octahedron. Exact rational convex bases, six wedges on original facets, and necessary Lanner/vertex gates. Other truncation orders and weight patterns are not classified.'})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--max-states',type=int,default=5000);a=p.parse_args();z=run(a.max_states);out=ROOT/f'results/{TAG}.json'
    if a.verify:assert z==json.loads(out.read_text())
    else:out.write_text(json.dumps(z,indent=2)+'\n')
    print('FINAL',Counter(r['verdict']for r in z['cases']),flush=True)
