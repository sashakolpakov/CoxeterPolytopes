#!/usr/bin/env python3
"""Coherent four-image completion by the newly certified earliest pi mirror."""
from fresh_h9_quartic_edge_walk_20261013 import *
import hashlib

SOURCE=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_third.json'
OUT=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_completion.json'


def run():
    started=time.monotonic();step=json.loads(SOURCE.read_text())
    path=ROOT/step['base_source'];old=json.loads(path.read_text())['cases'][step['base_source_case']]
    g=Host(2).g;fg=Host(1).g
    walls=[tuple(map(K.from_json,r))for r in old['all_walls']]
    face=[tuple(map(K.from_json,r))for r in old['forced_H5_face']['all_walls']]
    root=tuple(map(K.from_json,step['nearest_root']));q=K.from_json(step['nearest_norm'])
    r9=(root[0],*root[1:5],*(Z,)*4,root[5]);new=sorted(orbit(r9))
    assert q==2+T and len(new)==4
    cov=[mv(g,r)for r in walls];norms=[dot(r,c)for r,c in zip(walls,cov)]
    for r in new:
        assert pair(r,g,r)==q and all(z.is_integral for z in r)
        assert not all((z/T).is_integral for z in r)
        assert all((2*z/q).is_integral for z in mv(g,r))
        assert coxeter(r,q,cov,norms)
        cov.append(mv(g,r));norms.append(q)
    combined=walls+new;fc=[(*w[:5],w[9])for w in combined if not any(w[5:9])]
    inside=tuple(map(K.from_json,old['strict_interior']));fi=(*inside[:5],inside[9])
    edge_ids=step['edge_indices'];edges=[face[i]for i in edge_ids]
    ids9=[walls.index((e[0],*e[1:5],*(Z,)*4,e[5]))for e in edges]
    ids9 += [i for i,w in enumerate(walls)if not w[0]and any(w[5:9])]
    oldv=tuple(map(K.from_json,old['geometry']['vector']))
    geo=certificate(combined,g,inside,oldv)if all(pair(w,g,oldv)<=0 for w in combined)else scaled_geometry(combined,g,inside,9)
    oldf=tuple(map(K.from_json,old['forced_H5_face']['geometry']['vector']))
    fgeo=certificate(fc,fg,fi,oldf)if all(pair(w,fg,oldf)<=0 for w in fc)else scaled_geometry(fc,fg,fi,5)
    vertex=finite_vertex(fc,fg,fi,edge_ids+[fc.index(root)])
    normalized=tuple(z/tuple(map(K.from_json,vertex['vector']))[0]for z in map(K.from_json,vertex['vector']))
    assert enc(normalized)==step['nearest_point']
    out={'status':'EXACT_NEW_FIRST_MIRROR_COHERENT_COMPLETION','source':str(SOURCE.relative_to(ROOT)),
        'base_source':step['base_source'],'base_source_case':step['base_source_case'],
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,path)},
        'norm':q.to_json(),'height':(ALPHA*root[0]*root[0]/q).to_json(),
        'nearest_parameter':step['nearest_parameter'],'constructed_H5_root':enc(root),'new_roots':menc(new),
        'all_walls':menc(combined),'wall_count':len(combined),'strict_interior':enc(inside),
        'finite_H9_vertex':finite_vertex(combined,g,inside,ids9+[combined.index(r9)]),'geometry':geo,
        'old_H5_escape_pairing':pair(root,fg,oldf).to_json(),
        'forced_H5_face':{'all_walls':menc(fc),'wall_count':len(fc),'strict_interior':enc(fi),
                          'finite_vertex':vertex,'geometry':fgeo},
        'seconds':round(time.monotonic()-started,3),
        'scope':'The new norm-(2+T) mirror is the first mirror on the completely enumerated third edge. Its four symmetry images retain the original 182 walls and pass all cross-angle gates. This is an alternative 186-wall cone to the earlier norm-one unit-grid proposal. Exact escapes decide only the displayed cones, not the full host.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print({'status':out['status'],'walls':len(combined),'H5_walls':len(fc),'height':K.from_json(out['height']).real(),
           'geometry':geo['status'],'H5_geometry':fgeo['status'],'seconds':out['seconds']},flush=True)


if __name__=='__main__':run()
