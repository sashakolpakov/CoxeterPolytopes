#!/usr/bin/env python3
"""Merge eight certified A9 facets and verify their actual corner geometry."""
import hashlib,json
from itertools import combinations
from fresh_h9_a9_next_20261013_profiles import *
from fresh_h9_quartic_edge_walk_20261013 import finite_vertex
SOURCE_PROFILES=OUT
ORBIT=ROOT/'results/fresh_h9_a9_next_20261013_orbits.json'
OUT=ROOT/'results/fresh_h9_a9_next_20261013_geometry.json'


def robust_geometry(walls,inside):
    import mpmath as mp,numpy as np
    from scipy.optimize import linprog
    mp.mp.dps=100;t=mp.sqrt(2+mp.sqrt(2))
    def val(x):return float(sum((mp.mpf(v.numerator)/v.denominator)*t**i for i,v in enumerate(x.coeffs)))
    cov=[mv(G,r)for r in walls];M=np.array([[val(z)/max(1.,abs(val(r[0])))for z in c]for r,c in zip(walls,cov)])
    # Positive numerical row scaling is solely for LP conditioning; exact
    # cone certification below always uses the original walls.
    rng=np.random.default_rng(20261013);seen=set()
    for _ in range(300):
        lp=linprog(rng.normal(size=9),A_ub=M[:,1:],b_ub=-M[:,0],bounds=(None,None),method='highs')
        if not lp.success:continue
        ids=tuple(map(int,np.flatnonzero(abs(M[:,1:]@lp.x+M[:,0])<1e-7)))
        if len(ids)!=9 or ids in seen:continue
        seen.add(ids)
        try:v=(O,*mv(inverse(tuple(tuple(cov[i][j]for j in range(1,10))for i in ids)),tuple(-cov[i][0]for i in ids)))
        except(StopIteration,ZeroDivisionError):continue
        if any(dot(c,v)>0 for c in cov)or pair(v,G,v)<=0:continue
        eps=O
        while pair(tuple(x+eps*y for x,y in zip(v,inside)),G,tuple(x+eps*y for x,y in zip(v,inside)))<=0:eps/=2
        strict=tuple(x+eps*y for x,y in zip(v,inside));assert all(dot(c,strict)<0 for c in cov)and pair(strict,G,inside)<0
        return {'status':'EXACT_STRICT_SPACELIKE_ESCAPE','vector':enc(v),'norm':pair(v,G,v).to_json(),'active_facets':list(ids),
          'wall_pairings':enc(tuple(dot(c,v)for c in cov)),'strict_vector':enc(strict),'strict_norm':pair(strict,G,strict).to_json(),
          'strict_center_pairing':pair(strict,G,inside).to_json()}
    return {'status':'UNDECIDED','supports_tested':len(seen)}

def run():
    sd,pd,od=[json.loads(p.read_text())for p in(SOURCE,SOURCE_PROFILES,ORBIT)]
    walls=[tuple(map(K.from_json,r))for r in sd['all_walls']];inside=tuple(map(K.from_json,sd['strict_interior']));oldcount=len(walls)
    rev=lambda r:(r[0],*reversed(r[1:]));records=[]
    for case,pcase in zip(od['cases'],pd['cases'],strict=True):
        assert case['status']=='COMPLETE_FIRST_ENDPOINT_WITH_CERTIFIED_NORM_INTERFACE'
        root=tuple(map(K.from_json,case['selected']['vector']));point=tuple(map(K.from_json,case['selected']['point']));edge=[tuple(map(K.from_json,r))for r in pcase['edge_roots']]
        for swap in(False,True):
            r,p,es=(rev(root),rev(point),[rev(e)for e in edge])if swap else(root,point,edge)
            assert r not in walls;index=len(walls);walls.append(r)
            records.append({'case':case['case'],'reversed':swap,'wall_index':index,'point':enc(p),'edge':menc(es)})
    assert len(walls)==19
    for i,r in enumerate(walls):
        q=pair(r,G,r);assert q==2 and all(v.is_integral for v in r)and not all((v/T).is_integral for v in r)
        assert all((2*v/q).is_integral for v in mv(G,r))and compatible(r,q,walls[:i],G)
        assert pair(inside,G,r)<0
    proofs=[]
    for row in records:
        point=tuple(map(K.from_json,row['point']));ids=[walls.index(tuple(map(K.from_json,r)))for r in row['edge']]+[row['wall_index']]
        proof=finite_vertex(walls,G,inside,ids);v=tuple(map(K.from_json,proof['vector']))
        assert tuple(x/v[0]for x in v)==point and all(pair(point,G,r)<=0 for r in walls)
        proofs.append({**row,'certificate':proof})
    oldv=tuple(map(K.from_json,sd['geometry']['strict_vector']))
    if all(pair(oldv,G,r)<0 for r in walls):
        geo={'status':'EXACT_STRICT_SPACELIKE_ESCAPE','vector':enc(oldv),'norm':pair(oldv,G,oldv).to_json(),'wall_pairings':enc(tuple(pair(oldv,G,r)for r in walls)),
          'strict_vector':enc(oldv),'strict_norm':pair(oldv,G,oldv).to_json(),'strict_center_pairing':pair(oldv,G,inside).to_json(),'active_facets':[],'inherited_escape':True}
    else:geo=robust_geometry(walls,inside)
    out={'status':'EXACT_NINETEEN_WALL_A9_FIRST_EDGE_EXTENSION','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,SOURCE_PROFILES,ORBIT)},
      'gram':menc(G),'all_walls':menc(walls),'wall_count':len(walls),'new_wall_count':len(walls)-oldcount,'strict_interior':enc(inside),
      'actual_finite_vertices':proofs,'geometry':geo,
      'scope':'Eight coherent primitive norm-two facets from four completely certified first-edge steps and the integral reversal. Every pair and all eight actual finite vertices are checked exactly. A surviving exact spacelike escape proves this finite prefix has infinite volume, while full-host reflectivity remains unresolved unless a separate global certificate is supplied.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],geo['status'],len(proofs),flush=True)

if __name__=='__main__':run()
