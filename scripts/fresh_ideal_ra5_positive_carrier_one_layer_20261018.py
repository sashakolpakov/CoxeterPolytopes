#!/usr/bin/env python3
"""Exact genuine-facet/cusp attachment; the resulting RA5 chamber has infinite volume."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_ideal_ra5_positive_carrier_one_layer_20261018.json'
PINS={
 'results/fresh_ideal_ra4_positive_kappa_doubling_20261018.json':'cf5e71eb3720d16070852c19d3274f0a62b052ed38035d3f9d51adce666deb1c',
 'results/ideal_ra4_positive_kappa_doubling_20261018.md':'3ce0bfd400920c7fc6c4fd9b5b5e2cfab3b0e33fa8f4018a34ce34a84bfed724'}

def dot(a,b):return -a[0]*b[0]+sum(x*y for x,y in zip(a[1:],b[1:]))
def rank(rows):
    a=[list(r) for r in rows];k=0
    for j in range(6):
        p=next((i for i in range(k,len(a)) if a[i][j]),None)
        if p is None:continue
        a[k],a[p]=a[p],a[k];z=a[k][j];a[k]=[x/z for x in a[k]]
        for i in range(k+1,len(a)):
            z=a[i][j];a[i]=[x-z*y for x,y in zip(a[i],a[k])]
        k+=1
        if k==len(a):break
    return k
def enc(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,(list,tuple)):return [enc(y) for y in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x

def build():
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    old=json.loads((ROOT/next(iter(PINS))).read_text())['stages'][-1]
    ms=[tuple(map(Q,m))+(Q(0),) for m in old['walls']]
    vs=[tuple(map(Q,v))+(Q(0),) for v in old['vertices']]
    assert len(ms)==198 and len(vs)==258 and all(m[0]>0 for m in ms)
    n=(Q(0),)*5+(Q(1),);p=(Q(1),)+(Q(0),)*5;L=Q(48)
    h=[tuple(L*u-z for u,z in zip(v,n)) for v in vs]
    walls=[n]+ms+h;assert len(set(walls))==457
    a=min(-dot(v,m) for v in vs for m in ms if dot(v,m)<0)
    b=min(-dot(v,w) for v,w in combinations(vs,2))
    assert a==Q(1,25) and b==Q(1,1150)
    assert L*a>1 and L*L*b>2
    assert all(dot(x,x)==1 for x in walls)
    pair_count=0;orthogonal=0;parallel=0;ultraparallel=0
    for x,y in combinations(walls,2):
        z=dot(x,y);assert z==0 or z<=-1
        pair_count+=1;orthogonal+=z==0;parallel+=z==-1;ultraparallel+=z<-1
    assert pair_count==104196
    center=p[:-1]+(Q(-1,2),);escape=p[:-1]+(Q(-1),)
    assert dot(center,center)==Q(-3,4) and all(dot(center,m)<0 for m in walls)
    assert dot(escape,escape)==0 and all(dot(escape,m)<0 for m in walls)
    witnesses=[p]
    for m in ms:
        fv=[v for v in vs if dot(v,m)==0]
        q=tuple(sum(v[j] for v in fv)/len(fv) for j in range(6))
        assert dot(q,q)<0
        eps=min(Q(1,4),-dot(q,q)/4,min(-L*dot(v,q)/4 for v in vs))
        assert eps>0
        witnesses.append(tuple(x-eps*y for x,y in zip(q,n)))
    for v in vs:
        witnesses.append(tuple(u+(x-L*z)/(L*L) for u,x,z in zip(v,p,n)))
    assert len(witnesses)==len(walls)
    for i,x in enumerate(witnesses):
        assert x[0]>0 and dot(x,x)<0 and dot(x,walls[i])==0
        assert all(dot(x,m)<0 for j,m in enumerate(walls) if i!=j)
    cusps=[]
    for vi,v in enumerate(vs):
        ids=[i for i,m in enumerate(walls) if dot(v,m)==0]
        assert len(ids)==8 and 0 in ids and 199+vi in ids
        assert all(dot(v,m)<0 for i,m in enumerate(walls) if i not in ids)
        pairs=[(i,j) for i,j in combinations(ids,2) if dot(walls[i],walls[j])==-1]
        assert len(pairs)==4 and sorted(i for pair in pairs for i in pair)==ids
        assert all(dot(walls[i],walls[j])==0 for i,j in combinations(ids,2) if (i,j) not in pairs)
        assert rank([walls[i] for i in ids])==5
        scales=[]
        for i,j in pairs:
            s=tuple(x+y for x,y in zip(walls[i],walls[j]));k=s[0]
            assert k>0 and s==tuple(k*x for x in v);scales.append(k)
        cusps.append({'original_vertex':vi,'incident_walls':ids,'opposite_pairs':pairs,'positive_null_scales':scales})
    return enc({'status':'GENUINE_RA5_ONE_LAYER_WITH_258_CUSPS_AND_OPEN_IDEAL_ESCAPE',
      'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'sources':PINS,
      'lambda':L,'minimum_nonzero_vertex_wall_gap':a,'minimum_distinct_vertex_pair_gap':b,
      'wall_order':'Q facet; 198 old ridge-neighbor walls; 258 opposite-cusp walls',
      'walls':walls,'strict_center':center,'facet_relative_interior_witnesses':witnesses,
      'original_ideal_vertices':vs,'cusp_links':cusps,'strict_ideal_escape':escape,
      'pair_counts':{'total':pair_count,'orthogonal':orthogonal,'parallel':parallel,'ultraparallel':ultraparallel},
      'conclusions':{'genuine_facets':457,'original_Q_is_actual_facet':True,'closed_original_cusp_links':258,
        'volume':'infinite','finite_volume_RA5_example':False,'local_attachment_obstruction':False,
        'finite_volume_completion_excluded':False,'new_dimension_bound':False},
      'scope':'all one-layer walls and original cusp links; no claim of complete ideal boundary or finite-volume completion'} )

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=build()
    if a.verify:assert json.loads(OUT.read_text())==d
    else:
        assert not OUT.exists();OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({'status':'PASS','verify':a.verify,'facets':457,'closed_cusps':258,
      'pairs':104196,'lambda':48,'volume':'infinite'}))
