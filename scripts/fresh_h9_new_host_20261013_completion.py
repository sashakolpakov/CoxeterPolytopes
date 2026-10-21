#!/usr/bin/env python3
"""A coherent 18-wall construction with four new finite endpoint vertices."""
import hashlib,json,time
from fresh_h9_new_host_20261013_pell import *
from fresh_h9_quartic_edge_walk_20261013 import finite_vertex

PELL=OUT
OUT=ROOT/'results/fresh_h9_new_host_20261013_completion.json'

def symmetries(r):
    return sorted({(r[0],*a,*b,r[9])for a,b in((r[1:5],r[5:9]),(r[5:9],r[1:5]))for a,b in((a,b),(tuple(reversed(a)),tuple(reversed(b))))})

def checked_geometry(walls,g,inside,n):
    # Row rescaling is solely to condition the numerical proposals. Every saved
    # vector is subsequently checked against the original exact normals.
    scaled=[tuple(x/(r[0]if r[0]>1 else O)for x in r)for r in walls]
    result=geometry(scaled,g,inside,n)
    if result['status']=='EXACT_STRICT_SPACELIKE_ESCAPE':
        v=tuple(map(K.from_json,result['vector']));s=tuple(map(K.from_json,result['strict_vector']))
        assert all(pair(v,g,w)<=0 for w in walls)and all(pair(s,g,w)<0 for w in walls)
        result['wall_pairings']=enc(tuple(pair(v,g,w)for w in walls))
    return result

def run():
    started=time.monotonic();source=json.loads(SOURCE.read_text());fd=json.loads(FACE.read_text());pd=json.loads(PELL.read_text())
    g,fg=Host(2).g,Host(1).g;old=[tuple(map(K.from_json,r))for r in source['all_walls']]
    case=min(range(len(pd['compatible_candidates'])),key=lambda j:K.from_json(pd['compatible_candidates'][j]['height']))
    record=pd['compatible_candidates'][case];root=tuple(map(K.from_json,record['vector']));new=symmetries(root);assert len(new)==4
    walls=old+new;inside=tuple(map(K.from_json,source['strict_interior']))
    for i,r in enumerate(walls):
        q=pair(r,g,r);assert compatible(r,q,walls[:i],g)
        assert all(z.is_integral for z in r)and not all((z/T).is_integral for z in r)
        assert all((2*z/q).is_integral for z in mv(g,r))and pair(r,g,inside)<0
    fw=[tuple(map(K.from_json,r))for r in fd['all_walls']]+[(r[0],*r[1:5],r[9])for r in new if not any(r[5:9])]
    fi=tuple(map(K.from_json,fd['strict_interior']));assert len(fw)==10
    vertices=[]
    for i,r in enumerate(new):
        block=0 if any(r[1:5])else 1
        # Try the two middle directions; positive-definite Gram and exact
        # inequalities identify the finite endpoint without a numerical gate.
        candidates=[]
        for released in(1+4*block,2+4*block):
            ids=[j for j in range(9)if j!=released]+[14+i]
            try:proof=finite_vertex(walls,g,inside,ids)
            except(AssertionError,StopIteration,ZeroDivisionError):continue
            candidates.append(proof)
        assert candidates;vertices.append(candidates[0])
    out={'status':'EXACT_COHERENT_DYADIC_TWIST_ENDPOINT_EXTENSION','source':str(SOURCE.relative_to(ROOT)),
         'face_source':str(FACE.relative_to(ROOT)),'pell_source':str(PELL.relative_to(ROOT)),
         'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,FACE,PELL)},
         'pell_candidate_index':case,'new_roots':menc(new),'all_walls':menc(walls),'wall_count':len(walls),
         'gram':menc(g),'strict_interior':enc(inside),'finite_vertices':vertices,
         'geometry':checked_geometry(walls,g,inside,9),
         'forced_H5_face':{'all_walls':menc(fw),'wall_count':len(fw),'gram':menc(fg),'strict_interior':enc(fi),
                           'geometry':checked_geometry(fw,fg,fi,5)},
         'scope':'Four mutually Coxeter-compatible primitive reflection mirrors added to the complete H100 prefix of the changed arithmetic host. Each gives an actual finite H9 vertex with a positive definite rank-nine incident Gram. Neither first-mirror completeness nor finite volume is asserted.',
         'seconds':round(time.monotonic()-started,3)}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print({'walls':len(walls),'face_walls':len(fw),'vertices':len(vertices),'geometry':out['geometry']['status'],'face_geometry':out['forced_H5_face']['geometry']['status'],'seconds':out['seconds']})

if __name__=='__main__':run()
