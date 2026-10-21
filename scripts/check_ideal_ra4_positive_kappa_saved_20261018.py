#!/usr/bin/env python3
"""Independent saved reflection/face checks; inherits the original 24-cell census."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json

R=Path(__file__).resolve().parents[1]
D=json.loads((R/'results/fresh_ideal_ra4_positive_kappa_doubling_20261018.json').read_text())
assert D['source_sha256']==sha256((R/'scripts/fresh_ideal_ra4_positive_kappa_doubling_20261018.py').read_bytes()).hexdigest()
for p,h in D['sources'].items():assert sha256((R/p).read_bytes()).hexdigest()==h
def qset(rows):return {tuple(map(F,x)) for x in rows}
def bil(a,b):return sum(a[i]*b[i]*(1 if i else -1) for i in range(5))
def ref(v,n):return tuple(v[i]-2*bil(v,n)*n[i] for i in range(5))
def norm(v):
    assert v[0]>0
    return tuple(x/v[0] for x in v)
def reverse_rank(rows):
    a=[list(map(F,row)) for row in rows];r=0
    for c in (4,3,2,1,0):
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];z=a[r][c]
        for i in range(r+1,len(a)):
            y=a[i][c];a[i]=[z*x-y*b for x,b in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r

base=json.loads((R/'results/fresh_ideal_ra4_double24_volume_20261018.json').read_text())['double']
assert len(D['stages'])==4 and len(D['steps'])==3
assert qset(D['stages'][0]['walls'])==qset(base['walls'])
assert qset(D['stages'][0]['vertices'])==qset(base['vertices'])
for number,s in enumerate(D['stages']):
    walls=[tuple(map(F,n)) for n in s['walls']];vertices=[tuple(map(F,v)) for v in s['vertices']]
    assert len(set(walls))==len(walls) and len(set(vertices))==len(vertices)
    assert all(bil(n,n)==1 for n in walls)
    assert all(v[0]==1 and bil(v,v)==0 for v in vertices)
    assert all(bil(a,b)<=0 for a,b in combinations(walls,2))
    assert all(bil(v,n)<=0 for v in vertices for n in walls)
    fs=[{i for i,v in enumerate(vertices) if bil(v,n)==0} for n in walls]
    assert [sorted(f) for f in fs]==s['facet_vertices']
    assert len(set(map(frozenset,fs)))==len(walls)
    assert all(reverse_rank([vertices[i] for i in f])==4 for f in fs)
    c=tuple(map(F,s['strict_center']))
    assert c==tuple(sum(v[i] for v in vertices) for i in range(5))
    assert c[0]>0 and bil(c,c)<0 and all(bil(c,n)<0 for n in walls)
    for v in vertices:
        incident=[n for n in walls if bil(v,n)==0];assert len(incident)==6
        assert all(bil(a,b) in (0,-1) for a,b in combinations(incident,2))
        assert all(sum(bil(a,b)==-1 for b in incident)==1 for a in incident)
    v,f=len(vertices),len(walls)
    assert s['f_vector']==[v,4*v,3*v+f,f]
    assert s['kappa']==3*v-4*f+24 and s['delta']==v-f
    assert F(s['volume_coefficient_pi_squared'])==F(v-f+4,3)==F(4*s['cells'],3)
    if number:
        old=D['stages'][number-1];step=D['steps'][number-1]
        ns=qset(old['walls']);vs=qset(old['vertices']);n=tuple(map(F,step['wall']))
        assert n in ns and n==tuple(map(F,D['path'][number-1]))
        reflected={norm(ref(v,n)) for v in vs}
        assert set(vertices)==vs|reflected
        assert set(walls)==(ns-{n})|{ref(m,n) for m in ns if m!=n}
        shared={v for v in vs if bil(v,n)==0}
        assert shared==vs&reflected==qset(step['shared_vertices'])
        neighbors={m for m in ns if m!=n and bil(m,n)==0}
        assert neighbors==qset(step['neighbor_walls'])
        b=len(shared);assert b==step['gluing_facet_vertices']
        assert len(neighbors)==step['merged_facet_count']==b+2
        assert v==2*len(vs)-b and f==2*len(ns)-b-4
        assert s['kappa']==2*old['kappa']+b-8

start=D['stages'][0];sizes=list(map(len,start['facet_vertices']));ns=[tuple(map(F,n)) for n in start['walls']]
assert set(sizes)=={6,9}
assert all(len(set(start['facet_vertices'][i])&set(start['facet_vertices'][j]))==4
    for i,j in combinations(range(38),2) if sizes[i]==sizes[j]==9 and bil(ns[i],ns[j])==0)
assert D['bounded_doubling_level_checks']==[
 {'four_cell_gluing_size':6,'four_cell_kappa':-6,'max_facet_vertices':15,'eight_cell_kappa_upper':-5},
 {'four_cell_gluing_size':9,'four_cell_kappa':-3,'max_facet_vertices':14,'eight_cell_kappa_upper':0}]
assert D['conclusions']=={'positive_kappa':6,'final_f_vector':[258,1032,972,198],
 'volume_coefficient_pi_squared':'64/3','first_positive_doubling_level':16,
 'minimal_positive_kappa_polytope_globally':False,'ideal5_existence':False,'new_dimension_bound':False}
print(json.dumps({'status':'PASS','independent_saved_stages':4,'reflection_steps':3,
 'final_vertices':258,'final_facets':198,'positive_kappa':6,'new_dimension_bound':False}))
