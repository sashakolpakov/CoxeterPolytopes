#!/usr/bin/env python3
"""Exact existing-wall lookahead at the corrected finite H5 corner; no root search."""
from fresh_h9_quartic_edge_walk_20261013 import *
import hashlib

SOURCE=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_completion.json'
OUT=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_frontier.json'


def run():
    data=json.loads(SOURCE.read_text())['forced_H5_face'];g=Host(1).g
    walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    cov=[mv(g,r)for r in walls];v=tuple(map(K.from_json,data['finite_vertex']['vector']))
    v=tuple(z/v[0]for z in v);active=data['finite_vertex']['all_active_indices'];assert len(active)==5
    rows=[]
    for dropped in active:
        edge=[i for i in active if i!=dropped]
        order=edge+[dropped];matrix=tuple(tuple(cov[i][j]for j in range(1,6))for i in order)
        direction=(Z,*mv(inverse(matrix),(Z,Z,Z,Z,-O)))
        assert all(dot(direction,cov[i])==0 for i in edge)and dot(direction,cov[dropped])==-1
        crossings=[(-dot(v,c)/dot(direction,c),i)for i,c in enumerate(cov)if dot(direction,c)>0]
        assert all(t>0 for t,i in crossings)
        row={'edge_indices':edge,'released_starting_wall':dropped,'start':enc(v),'direction':enc(direction),
             'edge_gram':menc(tuple(tuple(pair(walls[i],g,walls[j])for j in edge)for i in edge))}
        if crossings:
            first=min(t for t,i in crossings);hit=[i for t,i in crossings if t==first]
            p=tuple(x+first*y for x,y in zip(v,direction));q=pair(p,g,p)
            assert all(dot(p,c)<=0 for c in cov)
            row.update(first_existing_wall_parameter=first.to_json(),first_existing_wall_indices=hit,
                       first_existing_wall_point=enc(p),first_existing_wall_point_norm=q.to_json(),
                       status='LIGHT_CONE_BEFORE_ANY_EXISTING_WALL'if q>0 else'FINITE_EXISTING_ENDPOINT'if q<0 else'IDEAL_EXISTING_ENDPOINT')
            if q>0:row['spacelike_certificate']=certificate(walls,g,tuple(map(K.from_json,data['strict_interior'])),p)
        else:row['status']='NO_EXISTING_WALL_AHEAD'
        rows.append(row)
    out={'status':'EXACT_EXISTING_WALL_FRONTIER_ONLY','source':str(SOURCE.relative_to(ROOT)),
         'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'starting_vertex_indices':active,
         'starting_vertex':enc(v),'edges':rows,
         'scope':'Every outgoing direction at this simple corrected H5 corner, with exact first intersection against the existing 76 walls. Unknown arithmetic mirrors are not enumerated here; open directions are the next residual-census tasks, not rays of the full reflection chamber.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print([(r['edge_indices'],r['status'],r.get('first_existing_wall_indices'))for r in rows],flush=True)


if __name__=='__main__':run()
