#!/usr/bin/env python3
"""Exact inductive reflection certificate for a positive-kappa ideal RA4 polytope."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_ideal_ra4_positive_kappa_doubling_20261018.json'
PINS={
 'results/fresh_ideal_ra4_double24_volume_20261018.json':'386d9650b8f669e3edcb3cdee56ed53d1c9ea4a90af2b4ba5c8c1b83f9d7bfa3',
 'results/ideal_ra4_second_volume_20261018.md':'9dff11c375bddcf62f16f3f01ddabd7f5b5c137756cf02ebff275076e6869c25'}

def dot(a,b):return -a[0]*b[0]+sum(x*y for x,y in zip(a[1:],b[1:]))
def reflection(x,n):return tuple(a-2*dot(x,n)*b for a,b in zip(x,n))
def project(x):
    assert x[0]>0
    return tuple(a/x[0] for a in x)
def rank(a):
    a=[list(map(Q,r)) for r in a];k=0
    for j in range(5):
        p=next((i for i in range(k,len(a)) if a[i][j]),None)
        if p is None:continue
        a[k],a[p]=a[p],a[k];d=a[k][j];a[k]=[x/d for x in a[k]]
        for i in range(k+1,len(a)):
            d=a[i][j];a[i]=[x-d*y for x,y in zip(a[i],a[k])]
        k+=1
        if k==len(a):break
    return k
def encoded(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,(tuple,list)):return [encoded(y) for y in x]
    if isinstance(x,dict):return {str(k):encoded(v) for k,v in x.items()}
    return x

def inspect(walls,vertices,cells):
    ns=sorted(walls);vs=sorted(vertices)
    assert all(dot(n,n)==1 for n in ns)
    assert all(v[0]==1 and dot(v,v)==0 for v in vs)
    assert all(dot(a,b)<=0 for a,b in combinations(ns,2))
    products=[[dot(v,n) for n in ns] for v in vs]
    assert all(x<=0 for r in products for x in r)
    facets=[{i for i,r in enumerate(products) if r[j]==0} for j in range(len(ns))]
    assert len(set(map(frozenset,facets)))==len(ns)
    assert all(rank([vs[i] for i in f])==4 for f in facets)
    center=tuple(sum(v[j] for v in vs) for j in range(5))
    assert center[0]>0 and dot(center,center)<0 and all(dot(center,n)<0 for n in ns)
    for r in products:
        ids=[i for i,x in enumerate(r) if x==0];assert len(ids)==6
        pp=[(i,j) for i,j in combinations(ids,2) if dot(ns[i],ns[j])==-1]
        assert len(pp)==3 and sorted(i for p in pp for i in p)==ids
        assert all(dot(ns[i],ns[j]) in (0,-1) for i,j in combinations(ids,2))
    v,f=len(vs),len(ns);delta=v-f;kappa=3*v-4*f+24
    assert Q(delta+4,3)==Q(4*cells,3)
    return {'cells':cells,'walls':ns,'vertices':vs,'facet_vertices':[sorted(f) for f in facets],
      'f_vector':[v,4*v,3*v+f,f],'delta':delta,'kappa':kappa,
      'volume_coefficient_pi_squared':Q(delta+4,3),'strict_center':center}

def build():
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    old=json.loads((ROOT/next(iter(PINS))).read_text())['double']
    ns={tuple(Q(x) for x in n) for n in old['walls']}
    vs={tuple(Q(x) for x in v) for v in old['vertices']}
    path=[(Q(1),Q(0),Q(1),Q(-1),Q(0)),
          (Q(1),Q(0),Q(1),Q(0),Q(-1)),
          (Q(1),Q(0),Q(1),Q(1),Q(0))]
    assert all(n in ns for n in path)
    fsets=[{v for v in vs if dot(v,n)==0} for n in path]
    assert list(map(len,fsets))==[9,9,9]
    assert len(fsets[0]&fsets[1])==len(fsets[1]&fsets[2])==4
    assert dot(path[0],path[1])==dot(path[1],path[2])==0
    assert dot(path[0],path[2])==-1
    stages=[inspect(ns,vs,2)];steps=[]
    for number,n in enumerate(path):
        assert n in ns;oldns=set(ns);oldvs=set(vs)
        shared={v for v in vs if dot(v,n)==0};b=len(shared)
        neighbors=[m for m in ns if m!=n and dot(m,n)==0]
        assert len(neighbors)==b+2
        rvs={project(reflection(v,n)) for v in vs}
        assert oldvs&rvs==shared
        ns=(ns-{n})|{reflection(m,n) for m in ns if m!=n}
        vs=vs|rvs
        row=inspect(ns,vs,2**(number+2));prev=stages[-1]
        assert len(vs)==2*len(oldvs)-b
        assert len(ns)==2*len(oldns)-b-4
        assert row['kappa']==2*prev['kappa']+b-8
        steps.append({'wall':n,'gluing_facet_vertices':b,'merged_facet_count':len(neighbors),
          'neighbor_walls':sorted(neighbors),'shared_vertices':sorted(shared)})
        stages.append(row)
    assert [(s['cells'],s['f_vector'][0],s['f_vector'][3],s['kappa']) for s in stages]==[
      (2,42,38,-2),(4,75,63,-3),(8,136,108,0),(16,258,198,6)]
    # A bounded level statement, not minimality among all RA4 polytopes.
    # Every first double is the frozen one by regular 24-cell facet transitivity.
    sizes=[len(f) for f in stages[0]['facet_vertices']]
    assert set(sizes)=={6,9}
    double_ns=stages[0]['walls'];double_fs=stages[0]['facet_vertices']
    assert all(len(set(double_fs[i])&set(double_fs[j]))==4
      for i,j in combinations(range(38),2) if sizes[i]==sizes[j]==9 and dot(double_ns[i],double_ns[j])==0)
    # Doubling an octahedron gives kappa=-6 and maximum new facet size <=15.
    # Doubling a D9 gives kappa=-3 and maximum new facet size <=14.
    level_bounds=[{'four_cell_gluing_size':6,'four_cell_kappa':-6,'max_facet_vertices':15,'eight_cell_kappa_upper':-5},
      {'four_cell_gluing_size':9,'four_cell_kappa':-3,'max_facet_vertices':14,'eight_cell_kappa_upper':0}]
    for r in level_bounds:assert r['eight_cell_kappa_upper']==2*r['four_cell_kappa']+r['max_facet_vertices']-8
    return encoded({'status':'POSITIVE_KAPPA_IDEAL_RA4_BY_THREE_FURTHER_DOUBLINGS',
      'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'sources':PINS,
      'path':path,'steps':steps,'stages':stages,'bounded_doubling_level_checks':level_bounds,
      'conclusions':{'positive_kappa':6,'final_f_vector':[258,1032,972,198],
        'volume_coefficient_pi_squared':Q(64,3),'first_positive_doubling_level':16,
        'minimal_positive_kappa_polytope_globally':False,'ideal5_existence':False,
        'new_dimension_bound':False},
      'scope':'inductive actual halfspace-union doubling; original bounded vertex census inherited, not rerun'} )

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=build()
    if a.verify:assert json.loads(OUT.read_text())==d
    else:
        assert not OUT.exists();OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({'status':'PASS','verify':a.verify,'stages':[(s['cells'],s['kappa']) for s in d['stages']],
      'final_f_vector':d['conclusions']['final_f_vector']}))
