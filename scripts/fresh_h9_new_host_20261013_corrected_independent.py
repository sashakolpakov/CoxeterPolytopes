#!/usr/bin/env python3
"""Independent corrected-cone geometry and existing-wall edge lookahead."""
import argparse,hashlib,json,copy
from functools import cmp_to_key
import fresh_h9_new_host_20261013_independent as a

SOURCE=a.ROOT/'results/fresh_h9_new_host_20261013_corrected.json'
OUT=a.ROOT/'results/fresh_h9_new_host_20261013_corrected_independent.json'

def check_vertices(data,host):
    walls=[a.vec(r)for r in data['all_walls']];center=a.vec(data['strict_interior'])
    for row in data['finite_vertices']:
        point=a.vec(row['vector']);vals=[host.pair(point,r)for r in walls]
        assert host.pair(point,point)==a.dec(row['norm'])and a.sign(host.pair(point,point))<0
        assert host.pair(point,center)==a.dec(row['center_pairing'])and a.sign(host.pair(point,center))<0
        assert vals==list(a.vec(row['all_wall_pairings']))and all(a.sign(v)<=0 for v in vals)
        assert [i for i,v in enumerate(vals)if v==a.Z]==row['all_active_indices']
        ids=row['indices'];assert len(ids)==host.n and all(vals[i]==a.Z for i in ids)
        gram=tuple(tuple(host.pair(walls[i],walls[j])for j in ids)for i in ids)
        assert gram==tuple(a.vec(r)for r in row['gram']);a.positive_definite(gram)

def replay(verify=False):
    data=json.loads(SOURCE.read_text());host=a.Host(2);fh=a.Host(1)
    for p,h in data['sources'].items():assert hashlib.sha256((a.ROOT/p).read_bytes()).hexdigest()==h
    first=json.loads((a.ROOT/'results/fresh_h9_new_host_20261013_first.json').read_text());r=a.vec(first['selected_choice']['vector']);r=(r[0],*r[1:5],*(a.Z,)*4,r[5])
    orbit={(r[0],*u,*v,r[9])for u,v in((r[1:5],r[5:9]),(r[5:9],r[1:5]))for u,v in((u,v),(u[::-1],v[::-1]))}
    assert orbit==set(a.vec(r)for r in data['new_roots'])and len(orbit)==4
    old=json.loads((a.ROOT/'results/fresh_h9_new_host_20261013_h9_h100.json').read_text())
    assert data['all_walls'][:14]==old['all_walls']and data['all_walls'][14:]==data['new_roots']
    h9check=host.check_geometry(data);h5check=fh.check_geometry(data['forced_H5_face'])
    check_vertices(data,host);check_vertices(data['forced_H5_face'],fh)
    face=data['forced_H5_face'];walls=[a.vec(r)for r in face['all_walls']];v=a.vec(face['finite_vertices'][0]['vector']);v=tuple(a.divide(z,v[0])for z in v)
    active=face['finite_vertices'][0]['all_active_indices'];assert len(active)==5
    assert {row['released_wall']for row in data['next_existing_wall_frontier']}==set(active)
    statuses=[]
    for row in data['next_existing_wall_frontier']:
        assert a.vec(row['start'])==v and row['edge_indices']==[i for i in active if i!=row['released_wall']]
        d=a.vec(row['direction']);assert d[0]==a.Z and all(fh.pair(d,walls[i])==a.Z for i in row['edge_indices'])
        assert fh.pair(d,walls[row['released_wall']])==a.neg(a.O)
        crossings=[(a.divide(a.neg(fh.pair(v,r)),fh.pair(d,r)),i)for i,r in enumerate(walls)if a.sign(fh.pair(d,r))>0]
        if not crossings:assert row['status']=='NO_EXISTING_WALL_AHEAD'
        else:
            t=min((t for t,i in crossings),key=cmp_to_key(lambda x,y:a.sign(a.sub(x,y))))
            assert a.sign(t)>0 and t==a.dec(row['parameter'])and [i for s,i in crossings if s==t]==row['first_existing_walls']
            point=tuple(a.add(x,a.mul(t,y))for x,y in zip(v,d));assert point==a.vec(row['point'])and all(a.sign(fh.pair(point,r))<=0 for r in walls)
            norm=fh.pair(point,point);assert norm==a.dec(row['norm'])
            expected='LIGHT_CONE_BEFORE_EXISTING_WALL'if a.sign(norm)>0 else'FINITE_EXISTING_ENDPOINT'if a.sign(norm)<0 else'IDEAL_EXISTING_ENDPOINT'
            assert expected==row['status']
        statuses.append(row['status'])
    forged=copy.deepcopy(data);forged['all_walls'][14][0][0]=str(a.Q(forged['all_walls'][14][0][0])+1)
    try:host.check_geometry(forged)
    except AssertionError:pass
    else:raise AssertionError('forged normal accepted')
    forged=copy.deepcopy(data);forged['geometry']['strict_vector']=forged['strict_interior']
    try:host.check_geometry(forged)
    except AssertionError:pass
    else:raise AssertionError('timelike escape accepted')
    result={'status':'PASS_CORRECTED_FIRST_MIRROR_CONE_AND_FRONTIER','source':str(SOURCE.relative_to(a.ROOT)),
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'H9':h9check,'H5':h5check,
        'actual_H9_vertices':len(data['finite_vertices']),'actual_H5_vertices':len(face['finite_vertices']),
        'frontier':statuses,'corruption_tests':{'forged_normal':'REJECTED','timelike_escape':'REJECTED'},'scope':'Separate Fraction replay of all corrected roots, reflection integrality, all selected Coxeter pairs, strict interiors, exact escapes, six finite vertices, and all five outgoing existing-wall lookahead directions. The complete first-mirror census has a separate independent certificate.'}
    if verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],h9check['walls'],statuses,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');replay(p.parse_args().verify)
