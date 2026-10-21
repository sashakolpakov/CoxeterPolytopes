#!/usr/bin/env python3
"""Exact rational double-24-cell certificate; no second-minimum assertion."""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/fresh_ideal_ra4_double24_volume_20261018.json'

def dot(a, b):
    return -a[0]*b[0] + sum(x*y for x,y in zip(a[1:], b[1:]))

def rank(rows):
    a = [[Q(x) for x in r] for r in rows]
    if not a: return 0
    k = 0
    for j in range(len(a[0])):
        p = next((i for i in range(k,len(a)) if a[i][j]), None)
        if p is None: continue
        a[k],a[p] = a[p],a[k]
        c = a[k][j]; a[k] = [x/c for x in a[k]]
        for i in range(k+1,len(a)):
            c = a[i][j]
            if c: a[i] = [x-c*y for x,y in zip(a[i],a[k])]
        k += 1
        if k == len(a): break
    return k

def solve(a,b):
    a = [[Q(x) for x in row]+[Q(y)] for row,y in zip(a,b)]
    n = len(b)
    for j in range(n):
        p = next((i for i in range(j,n) if a[i][j]),None)
        if p is None: return None
        a[j],a[p] = a[p],a[j]
        c=a[j][j]; a[j]=[x/c for x in a[j]]
        for i in range(n):
            if i!=j:
                c=a[i][j]; a[i]=[x-c*y for x,y in zip(a[i],a[j])]
    return tuple(row[-1] for row in a)

def reflect(x,n):
    return tuple(a-2*dot(x,n)*b/dot(n,n) for a,b in zip(x,n))

def normalize(x):
    assert x[0]>0
    return tuple(a/x[0] for a in x)

def encode(x):
    if isinstance(x,Q): return str(x)
    if isinstance(x,(tuple,list)): return [encode(v) for v in x]
    if isinstance(x,dict): return {str(k):encode(v) for k,v in x.items()}
    return x

def polytope(walls, vertices):
    walls, vertices = sorted(walls), sorted(vertices)
    assert all(dot(n,n)==1 for n in walls)
    assert all(v[0]==1 and dot(v,v)==0 for v in vertices)
    assert all(dot(v,n)<=0 for v in vertices for n in walls)
    assert all(dot(a,b)<=0 for a,b in combinations(walls,2))
    facets = [frozenset(i for i,v in enumerate(vertices) if dot(v,n)==0) for n in walls]
    assert len(set(facets))==len(walls)
    full=frozenset(range(len(vertices))); faces={full}; frontier={full}
    while frontier:
        new={f&h for f in frontier for h in facets}-{frozenset()}-faces
        faces.update(new); frontier=new
    bydim={i:[] for i in range(5)}
    for f in faces:
        base=vertices[min(f)]
        d=rank([[vertices[i][j]-base[j] for j in range(1,5)] for i in f])
        bydim[d].append(f)
    for d in bydim: bydim[d].sort(key=lambda f:tuple(sorted(f)))
    assert bydim[3]==sorted(facets,key=lambda f:tuple(sorted(f)))
    assert len(bydim[4])==1 and len(bydim[0])==len(vertices)
    incidence=[]
    for vi,v in enumerate(vertices):
        fs=[i for i,f in enumerate(facets) if vi in f]
        assert len(fs)==6
        parallel=[(a,b) for a,b in combinations(fs,2) if dot(walls[a],walls[b])==-1]
        assert len(parallel)==3 and sorted(a for p in parallel for a in p)==fs
        assert all(dot(walls[a],walls[b]) in (0,-1) for a,b in combinations(fs,2))
        assert sum(vi in e for e in bydim[1])==8
        assert sum(vi in f for f in bydim[2])==12
        incidence.append({'vertex':vi,'facets':fs,'parallel_pairs':parallel})
    for e in bydim[1]:
        assert len(e)==2 and sum(e<=f for f in facets)==3
    for f in bydim[2]:
        fs=[i for i,h in enumerate(facets) if f<=h]
        assert len(fs)==2 and dot(walls[fs[0]],walls[fs[1]])==0
    facet_types=Counter(); facet_cusps=[]
    for fi,f in enumerate(facets):
        es=[e for e in bydim[1] if e<=f]
        rs=[r for r in bydim[2] if r<=f]
        hist=Counter(map(len,rs))
        facet_types[(len(f),len(es),len(rs),tuple(sorted(hist.items())))]+=1
        cusp_shapes=[]
        for vi in sorted(f):
            adj=[j for j,h in enumerate(facets) if j!=fi and vi in h and f&h in bydim[2]]
            assert len(adj)==4
            opposite=[(a,b) for a,b in combinations(adj,2) if dot(walls[a],walls[b])==-1]
            assert len(opposite)==2 and sorted(a for p in opposite for a in p)==sorted(adj)
            sums=[tuple(a+b for a,b in zip(walls[j],walls[k])) for j,k in opposite]
            # Each sum of parallel unit wall normals equals kappa*v.
            assert all(s==tuple(s[0]*x for x in vertices[vi]) for s in sums)
            ratio=sums[0][0]/sums[1][0]
            assert ratio>0
            cusp_shapes.append({'vertex':vi,'aspect_ratio':max(ratio,1/ratio)})
        facet_cusps.append({'facet':fi,'vertices':len(f),'cusps':cusp_shapes,
          'nonsquare_cusps':sum(c['aspect_ratio']!=1 for c in cusp_shapes)})
    fv=[len(bydim[d]) for d in range(4)]
    ridge_hist=Counter(map(len,bydim[2]))
    defect=sum((k-3)*v for k,v in ridge_hist.items())
    assert fv[1]==4*fv[0] and sum((-1)**i*x for i,x in enumerate(fv))==0
    assert sum(k*v for k,v in ridge_hist.items())==12*fv[0]
    assert defect==3*(fv[0]-fv[3])
    return {'walls':walls,'vertices':vertices,'faces':{d:[sorted(f) for f in bydim[d]] for d in range(4)},
            'f_vector':fv,'ridge_histogram':dict(sorted(ridge_hist.items())),
            'facet_types':[{'vertices':k[0],'edges':k[1],'ridges':k[2],'ridge_histogram':dict(k[3]),'count':v} for k,v in sorted(facet_types.items())],
            'cusp_incidence':incidence,'facet_cusp_shapes':facet_cusps,'ridge_defect':defect,
            'volume_coefficient_pi_squared':Q(fv[0]-fv[3]+4,3),
            'reflection_orbifold_euler_characteristic':Q(fv[0]-fv[3]+4,4)}

def build():
    walls=[]
    for ij in combinations(range(4),2):
        for ss in product((-1,1),repeat=2):
            x=[Q(0)]*4
            for i,s in zip(ij,ss):x[i]=Q(s)
            walls.append((Q(1),*x))
    vertices=set()
    for i in range(4):
        for s in (-1,1):
            x=[Q(0)]*4;x[i]=Q(s);vertices.add((Q(1),*x))
    vertices.update((Q(1),*(Q(s,2) for s in ss)) for ss in product((-1,1),repeat=4))
    # All vertices of the bounded Klein halfspace intersection. Boundedness
    # follows from |x_i|+|x_j|<=1; the origin is a strict interior point.
    found=set(); feasible_bases=0
    for ids in combinations(range(24),4):
        x=solve([walls[i][1:] for i in ids],[1]*4)
        if x is not None and all(sum(a*b for a,b in zip(x,n[1:]))<=1 for n in walls):
            found.add((Q(1),*x));feasible_bases+=1
    assert found==vertices
    n=(Q(1),Q(1),Q(1),Q(0),Q(0));assert n in walls
    reflected_vertices={normalize(reflect(v,n)) for v in vertices}
    doubled_vertices=vertices|reflected_vertices
    doubled_walls=set(walls)-{n}
    doubled_walls.update(reflect(m,n) for m in walls if m!=n)
    original=polytope(walls,vertices)
    double=polytope(doubled_walls,doubled_vertices)
    assert original['f_vector']==[24,96,96,24]
    assert double['f_vector']==[42,168,164,38]
    assert double['ridge_histogram']=={3:152,4:12}
    assert all(f['nonsquare_cusps']==(0 if f['vertices']==6 else 3) for f in double['facet_cusp_shapes'])
    assert double['volume_coefficient_pi_squared']==2*original['volume_coefficient_pi_squared']==Q(8,3)
    shared=sorted(vertices&reflected_vertices);assert len(shared)==6
    return encode({'status':'EXACT_DOUBLE_CONSTRUCTION_AND_DISCRETE_VOLUME_REDUCTION',
      'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
      'published_volume_formula':'https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf',
      'metric':[-1,1,1,1,1],'original_active_four_subsets':10626,'feasible_active_bases':feasible_bases,
      'gluing_wall':n,'shared_cusps':shared,'original':original,'double':double,
      'universal_identities':{'edges':'4*v','ridges':'3*v+facets','ridge_defect':'3*(v-facets)',
        'volume':'(v-facets+4)*pi^2/3','facet_defect_sum':'6*(v-facets)'},
      'possible_intermediate_volume_coefficients_pi_squared':[Q(5,3),Q(2),Q(7,3)],
      'intermediate_ridge_defects':[3,6,9],
      'second_minimum_proved':False,'equality_classification_proved':False})

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    data=build()
    if args.verify: assert json.loads(OUT.read_text())==data
    else: OUT.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'status':'PASS','verify':args.verify,'double_f_vector':data['double']['f_vector'],
      'volume_pi_squared':data['double']['volume_coefficient_pi_squared'],'second_minimum_proved':False}))

if __name__=='__main__':main()
