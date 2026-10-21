#!/usr/bin/env python3
"""Merge the two complete next-edge orbits and certify their actual vertices."""
from fresh_h9_quartic_edge_walk_20261013 import *
import hashlib

SOURCE=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_completion.json'
OUT=ROOT/'results/fresh_h9_quartic_next_frontier_20261013_merge.json'


def run(with_long_proposals=False):
    started=time.monotonic();old=json.loads(SOURCE.read_text());g,fg=Host(2).g,Host(1).g
    walls=[tuple(map(K.from_json,r))for r in old['all_walls']]
    face=[tuple(map(K.from_json,r))for r in old['forced_H5_face']['all_walls']]
    inside=tuple(map(K.from_json,old['strict_interior']));fi=tuple(map(K.from_json,old['forced_H5_face']['strict_interior']))
    cov=[mv(g,r)for r in walls];norms=[pair(r,g,r)for r in walls];combined=walls[:];orbits=[];sources=[SOURCE]
    for i in(2,3):
        source=ROOT/f'results/fresh_h9_quartic_next_frontier_20261013_census_edge{i}.json';sources.append(source)
        step=json.loads(source.read_text());r=tuple(map(K.from_json,step['nearest_root']));q=K.from_json(step['nearest_norm'])
        r9=(r[0],*r[1:5],*(Z,)*4,r[5]);new=sorted(orbit(r9));assert len(new)==4
        for image in new:
            assert coxeter(image,q,cov,norms)
            combined.append(image);cov.append(mv(g,image));norms.append(q)
        orbits.append({'source':str(source.relative_to(ROOT)),'edge_index':i,'norm':q.to_json(),'new_roots':menc(new),
                       'height':(ALPHA*r[0]*r[0]/q).to_json(),'constructed_H5_root':enc(r)})
    if with_long_proposals:
        path=ROOT/'results/fresh_h9_quartic_next_frontier_20261013.json';sources.append(path)
        proposals=json.loads(path.read_text())
        for dictionary_index in(2,7):
            item=proposals['compatible_dictionary'][dictionary_index]
            r9=tuple(map(K.from_json,item['vector']));r=(*r9[:5],r9[9]);q=K.from_json(item['norm'])
            new=[tuple(map(K.from_json,rr))for rr in item['new_roots']]
            for image in new:
                assert coxeter(image,q,cov,norms)
                combined.append(image);cov.append(mv(g,image));norms.append(q)
            orbits.append({'source':str(path.relative_to(ROOT)),'edge_index':item['edge_index'],
                           'proposal_dictionary_index':dictionary_index,'norm':q.to_json(),'new_roots':menc(new),
                           'height':item['height'],'constructed_H5_root':enc(r),'nearest_mirror_certified':False})
    fc=[(*w[:5],w[9])for w in combined if not any(w[5:9])]
    for row in orbits:
        step=json.loads((ROOT/row['source']).read_text());r=tuple(map(K.from_json,row['constructed_H5_root']))
        edge_ids=step['edges'][row['edge_index']]['edge_indices']if 'proposal_dictionary_index'in row else step['edge_indices']
        edges=[face[i]for i in edge_ids];r9=(r[0],*r[1:5],*(Z,)*4,r[5])
        ids9=[walls.index((e[0],*e[1:5],*(Z,)*4,e[5]))for e in edges]
        ids9 += [i for i,w in enumerate(walls)if not w[0]and any(w[5:9])]
        row['finite_H9_vertex']=finite_vertex(combined,g,inside,ids9+[combined.index(r9)])
        row['finite_H5_vertex']=finite_vertex(fc,fg,fi,edge_ids+[fc.index(r)])
        point=tuple(map(K.from_json,row['finite_H5_vertex']['vector']))
        expected=step['compatible_dictionary'][row['proposal_dictionary_index']]['point']if 'proposal_dictionary_index'in row else step['nearest_point']
        assert enc(tuple(x/point[0]for x in point))==expected
    oldv=tuple(map(K.from_json,old['geometry']['vector']))
    geo=certificate(combined,g,inside,oldv)if all(pair(w,g,oldv)<=0 for w in combined)else scaled_geometry(combined,g,inside,9)
    oldf=tuple(map(K.from_json,old['forced_H5_face']['geometry']['vector']))
    fgeo=certificate(fc,fg,fi,oldf)if all(pair(w,fg,oldf)<=0 for w in fc)else scaled_geometry(fc,fg,fi,5)
    out={'status':'EXACT_MERGED_TWO_CERTIFIED_NEXT_EDGE_ORBITS','source':str(SOURCE.relative_to(ROOT)),
         'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in sources},
         'orbits':orbits,'all_walls':menc(combined),'wall_count':len(combined),'strict_interior':enc(inside),'geometry':geo,
         'forced_H5_face':{'all_walls':menc(fc),'wall_count':len(fc),'strict_interior':enc(fi),'geometry':fgeo},
         'seconds':round(time.monotonic()-started,3),'scope':'Two complete full-H9 nearest-edge certificates supply eight mutually compatible new facets. Both endpoints are actual finite vertices of the merged cone. Volume statements concern this displayed cone only.'}
    output=OUT.with_name(OUT.stem.replace('_merge','_four_edges')+'.json')if with_long_proposals else OUT
    if with_long_proposals:
        out['status']='EXACT_FOUR_OUTGOING_EDGE_SUBGROUP_COMPLETION'
        out['scope']='Two complete nearest-edge orbits and two additional long relative-unit proposals, all mutually compatible. The four outgoing edges of the corrected starting corner now have actual finite endpoints in this selected cone. First-mirror completeness is proved only for the two B3+A1 edges; the two long A2+A1+A1 proposals are not claimed nearest. Volume statements concern only the displayed cone.'
    output.write_text(json.dumps(out,indent=2)+'\n')
    print({'walls':len(combined),'face_walls':len(fc),'geometry':geo['status'],'face_geometry':fgeo['status'],'seconds':out['seconds']},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--with-long-proposals',action='store_true');args=parser.parse_args()
    run(args.with_long_proposals)
