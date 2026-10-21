#!/usr/bin/env python3
"""Separate Fraction replay of the corrected third-edge orbit and cone."""
from fresh_h9_quartic_edge_walk_20261013_independent import *
import hashlib

SOURCE=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_completion.json'
OUT=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_completion_independent.json'


def replay():
    data=json.loads(SOURCE.read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==s for p,s in data['source_sha256'].items())
    step=json.loads((ROOT/data['source']).read_text())
    base=json.loads((ROOT/data['base_source']).read_text())['cases'][data['base_source_case']]
    assert data['nearest_parameter']==step['nearest_parameter']and data['norm']==step['nearest_norm']
    root=vec(data['constructed_H5_root']);assert root==vec(step['nearest_root'])
    q=add(scale(O,2),T);assert dec(data['norm'])==q
    assert divide(mul(ALPHA,sq(root[0])),q)==dec(data['height'])
    r9=(root[0],*root[1:5],*(Z,)*4,root[5]);new=tuple(vec(r)for r in data['new_roots'])
    assert set(new)==orbit(r9)and len(new)==4
    combined=tuple(vec(r)for r in data['all_walls'])
    assert combined==tuple(vec(r)for r in base['all_walls'])+new
    face=data['forced_H5_face'];assert tuple(vec(r)for r in face['all_walls'])==tuple((*r[:5],r[9])for r in combined if all(z==Z for z in r[5:9]))
    h9=geometry_check(data,2);h5=geometry_check(face,1)
    v9=check_vertex(data['finite_H9_vertex'],data,Host(2));v5=check_vertex(face['finite_vertex'],face,Host(1))
    point=vec(face['finite_vertex']['vector']);assert tuple(divide(z,point[0])for z in point)==vec(step['nearest_point'])
    frontier=json.loads((ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_frontier.json').read_text())
    assert frontier['source']==str(SOURCE.relative_to(ROOT))and frontier['source_sha256']==hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    start=vec(frontier['starting_vertex']);assert start==vec(step['nearest_point'])
    active=face['finite_vertex']['all_active_indices'];assert active==frontier['starting_vertex_indices']
    fw=tuple(vec(r)for r in face['all_walls']);fh=Host(1)
    seen=set();open_count=0
    for row in frontier['edges']:
        released=row['released_starting_wall'];seen.add(released)
        edge=row['edge_indices'];assert edge==[i for i in active if i!=released]
        d=vec(row['direction']);assert d[0]==Z
        assert all(fh.pair(d,fw[i])==Z for i in edge)and fh.pair(d,fw[released])==neg(O)
        hits=[(divide(neg(fh.pair(start,r)),fh.pair(d,r)),i)for i,r in enumerate(fw)if sign(fh.pair(d,r))>0]
        assert all(sign(t)>0 for t,i in hits)
        t=dec(row['first_existing_wall_parameter'])
        assert all(sign(sub(x,t))>=0 for x,i in hits)
        assert [i for x,i in hits if x==t]==row['first_existing_wall_indices']
        p=tuple(add(x,mul(t,y))for x,y in zip(start,d));norm=fh.pair(p,p)
        assert p==vec(row['first_existing_wall_point'])and norm==dec(row['first_existing_wall_point_norm'])
        assert all(sign(fh.pair(p,r))<=0 for r in fw)
        if sign(norm)>0:
            assert row['status']=='LIGHT_CONE_BEFORE_ANY_EXISTING_WALL';open_count+=1
            geo=row['spacelike_certificate'];strict=vec(geo['strict_vector'])
            assert sign(fh.pair(strict,strict))>0 and all(sign(fh.pair(strict,r))<0 for r in fw)
        else:assert sign(norm)<0 and row['status']=='FINITE_EXISTING_ENDPOINT'
    assert seen==set(active)and open_count==4
    forged=list(root);forged[0]=add(forged[0],O);assert Host(1).pair(tuple(forged),tuple(forged))!=q
    return {'status':'PASS_INDEPENDENT_CORRECTED_THIRD_EDGE_COMPLETION','source':str(SOURCE.relative_to(ROOT)),
            'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'new_norm':data['norm'],
            'new_mirrors':4,'nearest_parameter':data['nearest_parameter'],'H9':h9,'H5':h5,
            'H9_vertex':v9,'H5_vertex':v5,'forged_root_rejected':True,
            'existing_wall_frontier':{'outgoing_edges':5,'light_cone_before_existing_wall':4,'previous_finite_corner':1},
            'scope':'All 17,205 H9 and 2,850 H5 pairs, primitive reflection integrality, strict interiors, exact strict escapes, actual finite vertices, and matching to the complete third-edge census. This certifies the displayed cones and nearest edge step, not global completion.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=replay()
    if args.verify:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
