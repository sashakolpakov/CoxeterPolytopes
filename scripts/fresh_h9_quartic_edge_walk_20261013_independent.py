#!/usr/bin/env python3
"""Fraction-only isometry, root, finite-vertex and full-cone replay."""
from fresh_h9_quartic_mixed_20261013_target_independent import *
from fresh_h9_quartic_mixed_20261013_endpoint_independent import pd

OUT=ROOT/'results/fresh_h9_quartic_edge_walk_20261013_independent.json'


def matrix(raw):return tuple(vec(row)for row in raw)
def mm(a,b):return tuple(tuple(dot(row,tuple(b[k][j]for k in range(len(b))))for j in range(len(b[0])))for row in a)
def ident(n):return tuple(tuple(O if i==j else Z for j in range(n))for i in range(n))


def check_vertex(proof,data,host):
    walls=tuple(vec(r)for r in data['all_walls']);center=vec(data['strict_interior']);point=vec(proof['vector'])
    inds=proof['indices'];assert len(inds)==host.n and len(set(inds))==len(inds)
    gram=tuple(tuple(host.pair(walls[i],walls[j])for j in inds)for i in inds)
    assert gram==matrix(proof['gram']);pd(gram)
    values=tuple(host.pair(point,r)for r in walls)
    assert values==vec(proof['all_wall_pairings'])and all(sign(v)<=0 for v in values)
    assert all(values[i]==Z for i in inds)
    assert [i for i,v in enumerate(values)if v==Z]==proof['all_active_indices']
    norm=host.pair(point,point);cp=host.pair(point,center)
    assert norm==dec(proof['norm'])and sign(norm)<0
    assert cp==dec(proof['center_pairing'])and sign(cp)<0
    return {'rank':len(inds),'all_inequalities':len(walls),'actual_finite_vertex':True,'norm':proof['norm']}


def replay_one(source):
    data=json.loads(source.read_text());old=json.loads((ROOT/data['source']).read_text());host=Host(2);face=Host(1)
    fw=tuple(vec(r)for r in old['forced_H5_face']['all_walls']);base=tuple(vec(r)for r in old['all_walls'])
    h,hi=matrix(data['translation']),matrix(data['translation_inverse']);j,p=matrix(data['generator']),matrix(data['plane_projection'])
    assert mm(h,hi)==ident(6)and mm(mm(tuple(zip(*h)),face.g),h)==face.g
    assert mm(p,p)==p and mm(j,p)==j
    assert mm(j,j)==tuple(tuple(mul(ALPHA,z)for z in row)for row in p)
    assert all(mv(h,fw[i])==fw[i]for i in data['fixed_H5_edge_indices'])
    y,s=dec(data['relative_unit']['y']),dec(data['relative_unit']['s'])
    assert sub(sq(y),mul(ALPHA,sq(s)))==O
    y0=divide(T,sub(scale(O,2),T));s0=divide(scale(O,2),sub(scale(O,2),T))
    power=data['relative_unit']['power_of_small_unit']
    assert (y,s)==((y0,s0)if power==1 else(add(sq(y0),mul(ALPHA,sq(s0))),scale(mul(y0,s0),2)))
    assert h==tuple(tuple(add(add(O if i==k else Z,mul(sub(y,O),p[i][k])),mul(s,j[i][k]))for k in range(6))for i in range(6))
    assert all(integral(z)for row in h for z in row)==data['translation_integral']
    integral_square=matrix(data['integral_translation_square'])
    assert integral_square==(mm(h,h)if power==1 else h)
    assert all(integral(z)for row in integral_square for z in row)
    assert all(integral(z)for row in a.minverse(integral_square)for z in row)
    transported=mv(h,fw[data['transported_H5_wall']]);assert transported==vec(data['raw_transported_root'])and sign(transported[0])<0
    oriented=tuple(neg(z)for z in transported);fixed=fw[data['reflection_H5_wall']]
    multiplier=divide(scale(face.pair(oriented,fixed),2),face.pair(fixed,fixed))
    assert multiplier==dec(data['reflection_multiplier'])
    root=tuple(sub(x,mul(multiplier,z))for x,z in zip(oriented,fixed))
    assert root==vec(data['constructed_H5_root'])and face.pair(root,root)==scale(O,2)
    assert divide(mul(ALPHA,sq(root[0])),scale(O,2))==dec(data['height'])
    r9=(root[0],*root[1:5],*(Z,)*4,root[5]);new=tuple(vec(r)for r in data['new_roots'])
    assert set(new)==orbit(r9)and len(new)==4
    assert tuple(vec(r)for r in data['all_walls'])==base+new
    assert tuple(vec(r)for r in data['forced_H5_face']['all_walls'])==tuple((*r[:5],r[9])for r in base+new if all(z==Z for z in r[5:9]))
    cut=face.pair(root,vec(old['forced_H5_face']['geometry']['vector']))
    assert cut==dec(data['old_H5_escape_pairing'])and sign(cut)>0
    g9=geometry_check(data,2);g5=geometry_check(data['forced_H5_face'],1)
    v9=check_vertex(data['finite_H9_vertex'],data,host)
    v5=check_vertex(data['forced_H5_face']['finite_vertex'],data['forced_H5_face'],face)
    forged=list(root);forged[0]=add(forged[0],O);assert face.pair(tuple(forged),tuple(forged))!=scale(O,2)
    return {'source':str(source.relative_to(ROOT)),'unit_power':power,'translation_integral':data['translation_integral'],
            'integral_square_and_inverse':True,'isometry_and_four_fixed_walls':True,'root_transport_and_reflection':True,
            'H9':g9,'H5':g5,'H9_vertex':v9,'H5_vertex':v5,'forged_root_rejected':True}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();start=time.monotonic()
    cases=[replay_one(ROOT/f'results/fresh_h9_quartic_edge_walk_20261013{suffix}.json')for suffix in('','_halfstep')]
    result={'status':'PASS_INDEPENDENT_TWO_EDGE_WALK_CONSTRUCTIONS','cases':cases,'seconds':round(time.monotonic()-start,3),
            'scope':'Separate Fraction field arithmetic; no producer import. Both displayed cones retain strict spacelike escapes. Actual finite vertices, not compactness, are certified.'}
    if args.verify:
        saved=json.loads(OUT.read_text());saved.pop('seconds');result.pop('seconds');assert saved==result
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
