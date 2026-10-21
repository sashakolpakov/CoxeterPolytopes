#!/usr/bin/env python3
"""Independent exact wall, facet-interior, cusp-link and escape replay."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json

R=Path(__file__).resolve().parents[1]
D=json.loads((R/'results/fresh_ideal_ra5_positive_carrier_one_layer_20261018.json').read_text())
assert D['source_sha256']==sha256((R/'scripts/fresh_ideal_ra5_positive_carrier_one_layer_20261018.py').read_bytes()).hexdigest()
for p,h in D['sources'].items():assert sha256((R/p).read_bytes()).hexdigest()==h
def bil(x,y):return sum((1 if i else -1)*x[i]*y[i] for i in range(6))
def rr(rows):
    a=[list(x) for x in rows];r=0
    for c in range(5,-1,-1):
        p=next((i for i in range(r,len(a)) if a[i][c]),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r];z=a[r][c]
        for i in range(r+1,len(a)):
            t=a[i][c];a[i]=[z*x-t*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r

old=json.loads((R/'results/fresh_ideal_ra4_positive_kappa_doubling_20261018.json').read_text())['stages'][-1]
walls=[tuple(map(F,x)) for x in D['walls']];vs=[tuple(map(F,x)) for x in D['original_ideal_vertices']]
assert len(walls)==len(set(walls))==457 and len(vs)==len(set(vs))==258
assert vs==[tuple(map(F,v))+(F(0),) for v in old['vertices']]
assert walls[1:199]==[tuple(map(F,m))+(F(0),) for m in old['walls']]
assert walls[0]==(0,0,0,0,0,1) and F(D['lambda'])==48
assert walls[199:]==[tuple(48*x-(1 if j==5 else 0) for j,x in enumerate(v)) for v in vs]
assert all(v[0]==1 and bil(v,v)==0 for v in vs)
assert all(bil(n,n)==1 for n in walls)
assert min(-bil(v,n) for v in vs for n in walls[1:199] if bil(v,n)<0)==F(D['minimum_nonzero_vertex_wall_gap'])==F(1,25)
assert min(-bil(v,w) for v,w in combinations(vs,2))==F(D['minimum_distinct_vertex_pair_gap'])==F(1,1150)
counts={'total':0,'orthogonal':0,'parallel':0,'ultraparallel':0}
for n,m in combinations(walls,2):
    z=bil(n,m);assert z==0 or z<=-1;counts['total']+=1
    counts['orthogonal']+=z==0;counts['parallel']+=z==-1;counts['ultraparallel']+=z<-1
assert counts==D['pair_counts'] and counts['total']==104196
points=[tuple(map(F,x)) for x in D['facet_relative_interior_witnesses']]
assert len(points)==457
for i,x in enumerate(points):
    assert x[0]>0 and bil(x,x)<0 and bil(x,walls[i])==0
    assert all(bil(x,n)<0 for j,n in enumerate(walls) if j!=i)
center=tuple(map(F,D['strict_center']));u=tuple(map(F,D['strict_ideal_escape']))
assert center==(1,0,0,0,0,F(-1,2)) and u==(1,0,0,0,0,-1)
assert bil(center,center)==F(-3,4) and bil(u,u)==0
assert all(bil(center,n)<0 and bil(u,n)<0 for n in walls)
assert len(D['cusp_links'])==258
for vi,(v,c) in enumerate(zip(vs,D['cusp_links'])):
    assert c['original_vertex']==vi
    ids=[i for i,n in enumerate(walls) if bil(v,n)==0]
    assert ids==c['incident_walls'] and len(ids)==8 and 0 in ids and 199+vi in ids
    assert all(bil(v,n)<0 for i,n in enumerate(walls) if i not in ids)
    pp=[list(t) for t in combinations(ids,2) if bil(walls[t[0]],walls[t[1]])==-1]
    assert pp==c['opposite_pairs'] and len(pp)==4 and sorted(j for t in pp for j in t)==ids
    assert rr([walls[i] for i in ids])==5
    for pair,k in zip(pp,c['positive_null_scales']):
        k=F(k);assert k>0
        assert tuple(x+y for x,y in zip(walls[pair[0]],walls[pair[1]]))==tuple(k*x for x in v)
assert D['conclusions']=={'genuine_facets':457,'original_Q_is_actual_facet':True,
 'closed_original_cusp_links':258,'volume':'infinite','finite_volume_RA5_example':False,
 'local_attachment_obstruction':False,'finite_volume_completion_excluded':False,'new_dimension_bound':False}
print(json.dumps({'status':'PASS','facets':457,'genuine_facet_witnesses':457,'closed_cusp_links':258,
 'pair_checks':104196,'strict_ideal_escape':True,'volume':'infinite'}))
