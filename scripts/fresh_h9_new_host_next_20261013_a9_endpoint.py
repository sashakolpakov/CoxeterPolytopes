#!/usr/bin/env python3
"""Exact middle-edge cap on the changed A9 host; no nearest-root claim."""
import hashlib
from fresh_h9_new_host_next_20261013_a9 import *
from fresh_h9_new_host_20261013_completion import checked_geometry
from fresh_h9_quartic_edge_walk_20261013 import finite_vertex

SOURCE=ROOT/'results/fresh_h9_new_host_next_20261013_a9_h650.json'
RAW=ROOT/'results/fresh_h9_new_host_next_20261013_a9_pell_raw.txt'
OUT=ROOT/'results/fresh_h9_new_host_next_20261013_a9_endpoint.json'

def run():
    d=json.loads(SOURCE.read_text());raw=RAW.read_text();assert 'CERTIFY 1\n'in raw and'CLASS 2\n'in raw
    reps=[tuple(map(K.from_json,json.loads(line.split(' ',1)[1])))for line in raw.splitlines()if line.startswith('INTEGRAL [')]
    assert len(reps)==1;z,x=reps[0];D=10*ALPHA
    if z<0:z=-z
    if x<0:x=-x
    assert z*z-D*x*x==20 and z>0 and x>0
    k=z/10;direction=tuple(K(i)for i in(1,2,3,4,5,4,3,2,1));assert pair(direction,A,direction)==10
    r=(x,*(k*y for y in direction));assert pair(r,G,r)==2 and all(v.is_integral for v in r)
    assert not all((v/T).is_integral for v in r)and all((2*v/2).is_integral for v in mv(G,r))
    walls=[tuple(map(K.from_json,w))for w in d['all_walls']]+[r];inside=tuple(map(K.from_json,d['strict_interior']))
    assert all(compatible(w,pair(w,G,w),walls[:i],G)for i,w in enumerate(walls))and all(pair(w,G,inside)<0 for w in walls)
    c=ALPHA*x/z;point=(O,*(c*y for y in direction));assert pair(point,G,r)==0
    assert pair(point,G,point)==-20*ALPHA/(z*z)and pair(point,G,point)<0
    assert all(pair(point,G,w)<=0 for w in walls)
    ids=[i for i in range(9)if i!=4]+[len(walls)-1]
    vertex=finite_vertex(walls,G,inside,ids)
    vv=tuple(map(K.from_json,vertex['vector']));assert tuple(y/vv[0]for y in vv)==point
    frontier=[]
    cov=[mv(G,w)for w in walls]
    for released in ids:
        edge=[i for i in ids if i!=released];order=edge+[released]
        v=(Z,*mv(inverse(tuple(tuple(cov[i][j]for j in range(1,10))for i in order)),tuple([Z]*8+[-O])))
        hits=[(-dot(point,p)/dot(v,p),i)for i,p in enumerate(cov)if dot(v,p)>0]
        assert all(t>0 for t,i in hits)
        row={'released_wall':released,'edge_indices':edge,'direction':enc(v)}
        if hits:
            t=min(t for t,i in hits);p=tuple(a+t*b for a,b in zip(point,v));q=pair(p,G,p)
            row.update(parameter=t.to_json(),point=enc(p),norm=q.to_json(),existing_walls=[i for s,i in hits if s==t],
               status='FINITE_EXISTING_ENDPOINT'if q<0 else'IDEAL_EXISTING_ENDPOINT'if q==0 else'LIGHT_CONE_BEFORE_EXISTING_WALL')
        else:row['status']='NO_EXISTING_WALL_AHEAD'
        frontier.append(row)
    out={'status':'EXACT_A9_MIDDLE_EDGE_VERTEX_CONSTRUCTION','source':str(SOURCE.relative_to(ROOT)),
      'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'discovery_source':str(RAW.relative_to(ROOT)),
      'discovery_source_sha256':hashlib.sha256(RAW.read_bytes()).hexdigest(),
      'relative_discriminant':D.to_json(),'relative_target':20,'z':z.to_json(),'time':x.to_json(),
      'root':enc(r),'root_norm':K(2).to_json(),'root_height':(ALPHA*x*x/2).to_json(),
      'all_walls':menc(walls),'wall_count':len(walls),'gram':menc(G),'strict_interior':enc(inside),
      'endpoint':enc(point),'endpoint_norm':pair(point,G,point).to_json(),'actual_finite_vertex':vertex,
      'incident_type':'A4+A4+A1','geometry':checked_geometry(walls,G,inside,9),
      'outgoing_existing_wall_frontier':frontier,
      'scope':'An exact primitive norm-two reflection cap, compatible with the entire saved A9 prefix, gives an actual finite H9 vertex with positive-definite rank-nine incident Gram. PARI supplied a proposal; the displayed polynomial norm identity and every geometric claim are replayable without PARI. This does not certify firstness along the edge or the complete full-group stabilizer at the proposed endpoint. Full-host reflectivity remains unresolved.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],len(walls),out['geometry']['status'],[r['status']for r in frontier],flush=True)

if __name__=='__main__':run()
