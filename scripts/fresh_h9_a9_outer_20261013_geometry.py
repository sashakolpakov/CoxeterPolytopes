#!/usr/bin/env python3
"""Merge a completely certified outer vertex star and its reversal."""
import hashlib,json
from fresh_h9_a9_outer_20261013_profiles import *
from fresh_h9_quartic_edge_walk_20261013 import finite_vertex
from fresh_h9_a9_next_20261013_geometry import robust_geometry
PROFILE=OUT
ORBIT=ROOT/'results/fresh_h9_a9_outer_20261013_orbits.json'
OUT=ROOT/'results/fresh_h9_a9_outer_20261013_geometry.json'

def run():
    sd,pd,od=[json.loads(p.read_text())for p in(SOURCE,PROFILE,ORBIT)]
    assert od['status']=='COMPLETE_SEVEN_OUTER_EDGE_NORM_INTERFACE_ROUND'
    walls=[tuple(map(K.from_json,r))for r in sd['all_walls']];inside=tuple(map(K.from_json,sd['strict_interior']));count=len(walls)
    rev=lambda r:(r[0],*reversed(r[1:]));records=[]
    for case,pcase in zip(od['cases'],pd['cases'],strict=True):
        assert case['status']=='COMPLETE_FIRST_ENDPOINT_WITH_CERTIFIED_NORM_INTERFACE'
        root=tuple(map(K.from_json,case['selected']['vector']));point=tuple(map(K.from_json,case['selected']['point']));edge=[tuple(map(K.from_json,r))for r in pcase['edge_roots']]
        for swap in(False,True):
            r,p,es=(rev(root),rev(point),[rev(e)for e in edge])if swap else(root,point,edge)
            if r not in walls:walls.append(r)
            records.append({'case':case['case'],'released_wall':case['released_wall'],'reversed':swap,'wall_index':walls.index(r),'point':enc(p),'edge':menc(es)})
    for i,r in enumerate(walls):
        q=pair(r,G,r);assert q==2 and all(v.is_integral for v in r)
        assert not all((v/T).is_integral for v in r)and not all((v/5).is_integral for v in r)
        assert all((2*v/q).is_integral for v in mv(G,r))and compatible(r,q,walls[:i],G)and pair(inside,G,r)<0
    proofs=[]
    for row in records:
        point=tuple(map(K.from_json,row['point']));ids=[walls.index(tuple(map(K.from_json,r)))for r in row['edge']]+[row['wall_index']]
        proof=finite_vertex(walls,G,inside,ids);v=tuple(map(K.from_json,proof['vector']))
        assert tuple(x/v[0]for x in v)==point and all(pair(point,G,r)<=0 for r in walls)
        proofs.append({**row,'certificate':proof})
    oldv=tuple(map(K.from_json,sd['geometry']['strict_vector']))
    if all(pair(oldv,G,r)<0 for r in walls):
        geo={'status':'EXACT_STRICT_SPACELIKE_ESCAPE','strict_vector':enc(oldv),'strict_norm':pair(oldv,G,oldv).to_json(),
          'strict_center_pairing':pair(oldv,G,inside).to_json(),'wall_pairings':enc(tuple(pair(oldv,G,r)for r in walls)),'inherited_escape':True}
    else:geo=robust_geometry(walls,inside)
    out={'status':'EXACT_OUTER_A9_VERTEX_STAR_EXTENSION','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,PROFILE,ORBIT)},
      'gram':menc(G),'all_walls':menc(walls),'wall_count':len(walls),'new_wall_count':len(walls)-count,'strict_interior':enc(inside),
      'actual_finite_vertices':proofs,'geometry':geo,
      'scope':'Seven completely certified first-edge steps and their reversals give fourteen actual finite vertices. Every saved wall and cross pair is exact. The central H8 facet star at the chosen outer vertex is complete, while the H9 edge releasing that facet is not treated. The volume assertion concerns this finite prefix only.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],'walls',len(walls),'new',len(walls)-count,'vertices',len(proofs),geo['status'],flush=True)

if __name__=='__main__':run()
