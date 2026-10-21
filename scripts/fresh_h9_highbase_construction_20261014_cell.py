#!/usr/bin/env python3
"""Exact 120-cell-base H9 multiwedges and complete Lanner3/4/5 joins.

The new label 7 is a category for every integer label >=7. This category
preserves the Lanner-three inequality and finite-diagram classification.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
from independent_h9_14_record30_a1_labels import lanner_codes,finite_component

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_highbase_construction_20261014_cell'
OUT=ROOT/f'results/{TAG}.json'
Z=(F(0),F(0));O=(F(1),F(0))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[0],-a[1])
def mul(a,b):return(a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])
def inv(a):
    d=a[0]*a[0]+a[0]*a[1]-a[1]*a[1];assert d
    return((a[0]+a[1])/d,-a[1]/d)
def dot(a,b):
    out=Z
    for x,y in zip(a,b):out=add(out,mul(x,y))
    return out
def sign(a):
    x,y=2*a[0]+a[1],a[1]
    sg=lambda t:int(t>0)-int(t<0)
    if not y:return sg(x)
    if not x:return sg(y)
    if sg(x)==sg(y):return sg(x)
    return sg(x)*sg(x*x-5*y*y)
def inverse(A):
    n=len(A);B=[list(row)+[O if i==j else Z for j in range(n)] for i,row in enumerate(A)]
    for i in range(n):
        j=next(j for j in range(i,n) if B[j][i]!=Z);B[i],B[j]=B[j],B[i]
        c=inv(B[i][i]);B[i]=[mul(c,x) for x in B[i]]
        for j in range(n):
            if i==j:continue
            c=B[j][i];B[j]=[add(x,neg(mul(c,y))) for x,y in zip(B[j],B[i])]
    return [row[n:] for row in B]
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(y) for y in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x
def digest(x):return sha256(json.dumps(enc(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()


def base():
    points=set()
    for i in range(4):
        for sg in [-1,1]:
            v=[Z]*4;v[i]=(F(2*sg),F(0));points.add(tuple(v))
    for ss in product([-1,1],repeat=4):points.add(tuple((F(x),F(0)) for x in ss))
    even=[p for p in permutations(range(4)) if sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2==0]
    for ss in product([-1,1],repeat=3):
        v=[Z,(F(ss[0]),F(0)),(F(0),F(ss[1])),(F(-ss[2]),F(ss[2]))]
        for p in even:points.add(tuple(v[i] for i in p))
    V=sorted(points);assert len(V)==120
    G=[[dot(v,w) for w in V] for v in V]
    assert all(G[i][i]==(4,0) for i in range(120))
    edges={(i,j) for i,j in combinations(range(120),2) if G[i][j]==(0,2)}
    adjacency=[{j for j in range(120) if tuple(sorted((i,j))) in edges} for i in range(120)]
    assert set(map(len,adjacency))=={12} and len(edges)==720
    cliques=[[(i,) for i in range(120)]]
    while cliques[-1]:
        nxt=[]
        for c in cliques[-1]:
            for j in range(c[-1]+1,120):
                if all(j in adjacency[i] for i in c):nxt.append(c+(j,))
        cliques.append(nxt)
    assert list(map(len,cliques))==[120,720,1200,600,0]
    facets=cliques[3]
    supports=[]
    for f in facets:
        n=tuple(tuple(sum(V[i][k][a] for i in f) for a in range(2)) for k in range(4))
        vals=[dot(n,v) for v in V]
        assert all(vals[i]==(4,6) for i in f)
        assert all(sign(add(vals[i],(-4,-6)))<0 for i in range(120) if i not in f)
        # The four position vectors are independent, so their supporting
        # hyperplane section is a genuine tetrahedron.
        inverse([[G[i][j] for j in f] for i in f])
        supports.append(n)
    ridge=Counter(t for f in facets for t in combinations(f,3));assert set(ridge.values())=={2}
    reached={facets[0]};todo=list(reached)
    by_ridge={t:[] for t in ridge}
    for f in facets:
        for t in combinations(f,3):by_ridge[t].append(f)
    for f in todo:
        for t in combinations(f,3):
            for g in by_ridge[t]:
                if g not in reached:reached.add(g);todo.append(g)
    assert len(reached)==600
    assert all(tuple(neg(x) for x in v) in points for v in V)
    # Certified supporting tetrahedra form a connected closed facet patch.
    # A convex-polytope facet graph is connected, so no further facet exists.
    f=facets[0];I=inverse([[G[i][j] for j in f] for i in f]);lookup={v:i for i,v in enumerate(V)}
    coordinates=[[dot(row,[G[i][j] for j in f]) for row in I] for i in range(120)]
    group=[]
    for p in permutations(range(4)):
        perm=[]
        for co in coordinates:
            y=tuple(dot(co,[V[f[p[j]]][k] for j in range(4)]) for k in range(4))
            assert y in lookup;perm.append(lookup[y])
        assert len(set(perm))==120
        assert {tuple(sorted(perm[i] for i in t)) for t in facets}==set(facets)
        group.append(tuple(perm))
    assert len(set(group))==24
    GS=set(group);assert all(tuple(p[q[i]] for i in range(120)) in GS for p in group for q in group)
    remaining=set(range(120))-set(f);orbits=[]
    while remaining:
        i=min(remaining);orb={p[i] for p in group};assert orb<=remaining
        orbits.append(sorted(orb));remaining-=orb
    missing=[p for p in combinations(range(120),2) if p not in edges]
    recipes=[]
    for orb in orbits:
        chosen=tuple(sorted((*f,orb[0])))
        l4=sum(set(m)<=set(chosen) for m in missing)
        assert 1<=l4<=4
        recipes.append(dict(doubled_vertices=chosen,fifth_vertex_orbit=orb,mandatory_lanner4=l4))
    recipes.sort(key=lambda r:(r['mandatory_lanner4'],r['doubled_vertices']))
    b=dict(name='simple 120-cell dual to the exact 600-cell',dimension=4,
           dual_vertices=V,dual_facets=facets,exact_facet_normals=supports,
           support_level=(4,6),minimal_nonfaces=missing,chosen_tetrahedron=f,
           tetrahedron_stabilizer=group,fifth_vertex_orbits=orbits)
    return b,recipes


def incidence(b,chosen):
    weights=[2 if i in chosen else 1 for i in range(120)]
    blocks=[];n=0
    for w in weights:blocks.append(tuple(range(n,n+w)));n+=w
    assert n==125
    M=[tuple(sorted(i for j in m for i in blocks[j])) for m in b['minimal_nonfaces']]
    V=set()
    for f in b['dual_facets']:
        outside=[j for j in range(120) if j not in f]
        for omit in product(*(blocks[j] for j in outside)):
            V.add(tuple(i for i in range(n) if i not in omit))
    assert all(len(v)==9 for v in V)
    assert set(Counter(t for v in V for t in combinations(v,8)).values())=={2}
    return weights,blocks,sorted(M),sorted(V)


@lru_cache(maxsize=1048576)
def finite(order,code):
    labels=dict(zip(combinations(range(order),2),code));left=set(range(order));cc=[]
    while left:
        c={min(left)};left-=c;todo=list(c)
        for i in todo:
            nn={j for j in left if labels[tuple(sorted((i,j)))]>2}
            left-=nn;c|=nn;todo+=sorted(nn)
        cc.append(c)
    return all(len(c)<=2 or finite_component(c,labels) for c in cc)


@lru_cache(None)
def codes(k):
    if k!=3:return lanner_codes(k)
    return tuple(c for c in product(range(2,8),repeat=3) if sum(F(1,m) for m in c)<1)


def join(M,V,maxstates):
    supports=[m for m in M if len(m)>=3];remaining=set(range(len(supports)));used=set();states=[()];stages=[]
    while states and remaining:
        k=max(remaining,key=lambda j:(len(set(supports[j])&used),len(supports[j]),-j))
        B=supports[k];used.update(B);pairs=list(combinations(B,2));opts=[tuple(zip(pairs,c)) for c in codes(len(B))]
        restrictions=sorted({tuple(i for i in v if i in used) for v in V},key=lambda a:(-len(a),a))
        maximal=[]
        for a in restrictions:
            if not any(set(a)<=set(c) for c in maximal):maximal.append(a)
        passing=set();examined=compatible=0;truncated=False
        for state in states:
            old=dict(state)
            for option in opts:
                examined+=1
                if any(p in old and old[p]!=m for p,m in option):continue
                compatible+=1;new=old|dict(option)
                if all(finite(len(a),tuple(new.get(p,2) for p in combinations(a,2))) for a in maximal):passing.add(tuple(sorted(new.items())))
                if len(passing)>maxstates:truncated=True;break
            if truncated:break
        stages.append(dict(missing_face=B,option_count=len(opts),input_states=len(states),
                           extensions_examined=examined,compatible_extensions=compatible,
                           vertex_restrictions=len(maximal),passing_states=len(passing),
                           state_hash=digest(sorted(passing)),truncated=truncated))
        states=sorted(passing);remaining.remove(k)
        print(' join',B,'states',len(states),'used',len(used),'truncated',truncated,flush=True)
        if truncated:break
    return dict(stages=stages,unprocessed_missing_faces=len(remaining),state_limit_reached=bool(stages[-1]['truncated']),
                terminal_states=states if states and not remaining else [],
                verdict='EXACT_LANNER_VERTEX_OBSTRUCTION' if not states else 'LABEL_SURVIVOR_OR_EXPLICIT_LIMIT')


def run(limit,maxstates):
    b,recipes=base();print('600-cell base, fifth-vertex orbits',len(recipes),flush=True)
    rows=[]
    for i,recipe in enumerate(recipes[:limit]):
        w,B,M,V=incidence(b,recipe['doubled_vertices'])
        print('case',i,recipe['mandatory_lanner4'],'L4',len(V),'vertices',flush=True)
        rows.append(dict(case=i,**recipe,weights=w,blocks=B,minimal_nonfaces=M,actual_vertices=V,
                         **join(M,V,maxstates)))
    files=[f'scripts/{TAG}.py','scripts/independent_h9_14_record30_a1_labels.py']
    return enc(dict(status='EXACT_120CELL_MULTIWEDGE_LABEL_CHECKPOINT',
                    sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                    settings=dict(limit=limit,max_states=maxstates),base=b,all_recipes=recipes,cases=rows,
                    complete_selected_recipe_family=len(rows)==len(recipes),lanner3_category_count=len(codes(3)),
                    scope='Five wedges on a tetrahedron plus a fifth vertex of the 600-cell, modulo its complete S4 stabilizer. '
                          'Labels 7 encode all integer labels >=7 in necessary Lanner3/finite tests. '
                          'Empty untruncated joins exclude only their incidence recipes. Surviving states need global metric and volume checks.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=100);p.add_argument('--max-states',type=int,default=20000);p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run(a.limit,a.max_states)))
    if a.verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],Counter(r['verdict'] for r in out['cases']))
