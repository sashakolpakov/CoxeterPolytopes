#!/usr/bin/env python3
"""Replace far proposals by the certified first mirrors and expose next edges."""
import hashlib,json,time
from fresh_h9_new_host_20261013_completion import symmetries,checked_geometry
from fresh_h9_new_host_20261013_edge import *
from fresh_h9_quartic_edge_walk_20261013 import finite_vertex

FIRST=ROOT/'results/fresh_h9_new_host_20261013_first.json'
OUT=ROOT/'results/fresh_h9_new_host_20261013_corrected.json'

def frontier(walls,g,proof):
    cov=[mv(g,w)for w in walls];v=tuple(map(K.from_json,proof['vector']));v=tuple(z/v[0]for z in v);active=proof['all_active_indices'];assert len(active)==5
    rows=[]
    for dropped in active:
        edge=[j for j in active if j!=dropped];order=edge+[dropped]
        d=(Z,*mv(inverse(tuple(tuple(cov[i][j]for j in range(1,6))for i in order)),(Z,Z,Z,Z,-O)))
        assert all(dot(d,cov[i])==0 for i in edge)and dot(d,cov[dropped])==-1
        crossings=[(-dot(v,c)/dot(d,c),i)for i,c in enumerate(cov)if dot(d,c)>0];assert all(t>0 for t,i in crossings)
        row={'edge_indices':edge,'released_wall':dropped,'start':enc(v),'direction':enc(d)}
        if crossings:
            t=min(t for t,i in crossings);p=tuple(x+t*y for x,y in zip(v,d));q=pair(p,g,p)
            assert all(dot(p,c)<=0 for c in cov)
            row.update(parameter=t.to_json(),first_existing_walls=[i for s,i in crossings if s==t],point=enc(p),norm=q.to_json(),
                status='LIGHT_CONE_BEFORE_EXISTING_WALL'if q>0 else'FINITE_EXISTING_ENDPOINT'if q<0 else'IDEAL_EXISTING_ENDPOINT')
        else:row['status']='NO_EXISTING_WALL_AHEAD'
        rows.append(row)
    return rows

def run():
    started=time.monotonic();data=json.loads(SOURCE.read_text());fd=json.loads(FACE.read_text());first=json.loads(FIRST.read_text())
    g,fg=Host(2).g,Host(1).g;old=[tuple(map(K.from_json,r))for r in data['all_walls']]
    r=tuple(map(K.from_json,first['selected_choice']['vector']));root=(r[0],*r[1:5],*(Z,)*4,r[5]);new=symmetries(root);assert len(new)==4
    walls=old+new;inside=tuple(map(K.from_json,data['strict_interior']))
    for i,w in enumerate(walls):assert compatible(w,pair(w,g,w),walls[:i],g)and pair(w,g,inside)<0
    fw=[tuple(map(K.from_json,r))for r in fd['all_walls']]+[(r[0],*r[1:5],r[9])for r in new if not any(r[5:9])]
    fi=tuple(map(K.from_json,fd['strict_interior']));fv=[];vertices=[]
    for i,r in enumerate(new):
        block=0 if any(r[1:5])else 1
        proofs=[]
        for released in(1+4*block,2+4*block):
            try:proof=finite_vertex(walls,g,inside,[j for j in range(9)if j!=released]+[14+i])
            except(AssertionError,StopIteration,ZeroDivisionError):continue
            proofs.append(proof)
        assert len(proofs)==1;vertices.append(proofs[0])
    for j in(8,9):
        for released in(1,2):
            try:proof=finite_vertex(fw,fg,fi,[i for i in range(5)if i!=released]+[j])
            except(AssertionError,StopIteration,ZeroDivisionError):continue
            fv.append(proof)
    assert len(fv)==2
    out={'status':'EXACT_COHERENT_FIRST_MIRROR_EXTENSION','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,FACE,FIRST)},
      'all_walls':menc(walls),'wall_count':len(walls),'new_roots':menc(new),'gram':menc(g),'strict_interior':enc(inside),
      'finite_vertices':vertices,'geometry':checked_geometry(walls,g,inside,9),
      'forced_H5_face':{'all_walls':menc(fw),'wall_count':len(fw),'gram':menc(fg),'strict_interior':enc(fi),
        'finite_vertices':fv,'geometry':checked_geometry(fw,fg,fi,5)},
      'next_existing_wall_frontier':frontier(fw,fg,fv[0]),'seconds':round(time.monotonic()-started,3),
      'scope':'Four symmetry copies of the completely certified first mirror replace the earlier far proposals. These are facets of the same full arithmetic reflection chamber. This18-wall prefix remains infinite-volume; the full host is unresolved. The next frontier is exact lookahead against existing walls only, not a census of unknown mirrors.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print({'walls':len(walls),'face_walls':len(fw),'new_vertices':len(vertices),'geometry':out['geometry']['status'],'face_geometry':out['forced_H5_face']['geometry']['status'],'frontier':[(r['edge_indices'],r['status'])for r in out['next_existing_wall_frontier']],'seconds':out['seconds']},flush=True)

if __name__=='__main__':run()
