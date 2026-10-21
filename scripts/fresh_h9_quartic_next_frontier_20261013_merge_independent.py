#!/usr/bin/env python3
"""Separate exact root, angle, finite-vertex and escape replay for both new cones."""
from fresh_h9_quartic_edge_walk_20261013_independent import *
import hashlib

OUT=ROOT/'results/fresh_h9_quartic_next_frontier_20261013_merge_independent.json'


def replay_one(suffix):
    path=f'results/fresh_h9_quartic_next_frontier_20261013_{suffix}.json'
    data=json.loads((ROOT/path).read_text());base=json.loads((ROOT/data['source']).read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==s for p,s in data['source_sha256'].items())
    new=[];vertices=[]
    for row in data['orbits']:
        source=json.loads((ROOT/row['source']).read_text());root=vec(row['constructed_H5_root']);q=dec(row['norm'])
        if 'proposal_dictionary_index'in row:
            item=source['compatible_dictionary'][row['proposal_dictionary_index']]
            r9=vec(item['vector']);assert root==(*r9[:5],r9[9])and q==dec(item['norm'])
            expected=vec(item['point'])
        else:
            assert root==vec(source['nearest_root'])and q==dec(source['nearest_norm'])
            expected=vec(source['nearest_point']);r9=(root[0],*root[1:5],*(Z,)*4,root[5])
        assert divide(mul(ALPHA,sq(root[0])),q)==dec(row['height'])
        rs=tuple(vec(r)for r in row['new_roots']);assert set(rs)==orbit(r9)and len(rs)==4
        new.extend(rs)
        v9=check_vertex(row['finite_H9_vertex'],data,Host(2))
        v5=check_vertex(row['finite_H5_vertex'],data['forced_H5_face'],Host(1))
        point=vec(row['finite_H5_vertex']['vector']);assert tuple(divide(x,point[0])for x in point)==expected
        vertices.append({'edge_index':row['edge_index'],'H9':v9,'H5':v5,'certified_nearest':'proposal_dictionary_index'not in row})
    assert tuple(vec(r)for r in data['all_walls'])==tuple(vec(r)for r in base['all_walls'])+tuple(new)
    assert tuple(vec(r)for r in data['forced_H5_face']['all_walls'])==tuple((*r[:5],r[9])for r in map(vec,data['all_walls'])if all(z==Z for z in r[5:9]))
    h9=geometry_check(data,2);h5=geometry_check(data['forced_H5_face'],1)
    return {'source':path,'source_sha256':hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),
            'new_roots':len(new),'H9':h9,'H5':h5,'actual_vertices':vertices}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();start=time.monotonic()
    result={'status':'PASS_INDEPENDENT_TWO_NEW_FRONTIER_CONES','cases':[replay_one(s)for s in('merge','four_edges')],
            'scope':'Separate Fraction primitive-root/reflection-integrality, all pair angle gates, strict timelike interiors, exact strict escapes, and six actual finite-vertex certificates. Four-edge long proposals are not certified first mirrors.'}
    if args.verify:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print({'status':result['status'],'wall_counts':[(r['H9']['walls'],r['H5']['walls'])for r in result['cases']],
           'seconds':round(time.monotonic()-start,3)},flush=True)
