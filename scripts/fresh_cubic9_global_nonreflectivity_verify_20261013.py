#!/usr/bin/env python3
"""Fraction-only full-chamber nonreflectivity certificate for one cubic9 host.

Repeats every dominant root through height3000 and tests both endpoints.
The matrix proposal source is untrusted; all its identities are replayed here.
"""
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

import fresh_cubic9_dihedral_independent_audit_20261013 as f
import fresh_cubic9_ellipsoid_audit_20261013 as e

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_cubic9_dihedral_h3_theta_beta_h3000_20261013.json'
SCOUT=ROOT/'results/fresh_cubic9_vertex_symmetry_theta_beta_h3000_20261013.json'
OUTPUT=ROOT/'results/fresh_cubic9_global_nonreflectivity_20261013.json'

def transpose(A):return tuple(zip(*A))
def mm(A,B):return tuple(tuple(f.dot(r,c) for c in transpose(B)) for r in A)
def encode_matrix(A):return [f.vecencode(row) for row in A]
def eye(n):return tuple(tuple(f.ONE if i==j else f.ZERO for j in range(n)) for i in range(n))

def local_check(d,matrix,point):
    alpha=f.decode(d['alpha']);assert alpha==f.mul(f.THETA,f.BETA)
    assert d['blocks']==1 and d['scalar'] and d['dimension']==3
    G=f.make_gram(alpha,1,True);H=f.decode(d['height_bound']);assert H==f.scale(f.ONE,3000)
    U=tuple(f.vecdecode(row) for row in matrix);p=f.vecdecode(point)
    assert all(f.integral(x) for row in U for x in row)
    assert mm(mm(transpose(U),G),U)==G
    invG=((f.neg(f.inverse(alpha)),f.ZERO,f.ZERO,f.ZERO),
          (f.ZERO,f.divide(f.TWO,f.BETA),f.divide(f.THETA,f.BETA),f.ZERO),
          (f.ZERO,f.divide(f.THETA,f.BETA),f.divide(f.TWO,f.BETA),f.ZERO),
          (f.ZERO,f.ZERO,f.ZERO,f.ONE))
    inverse=mm(mm(invG,transpose(U)),G)
    assert mm(U,inverse)==mm(inverse,U)==eye(4)
    assert all(f.integral(x) for row in inverse for x in row)
    tr=f.sum_field(U[i][i] for i in range(4));assert tr==f.add(f.add(f.ONE,f.scale(f.THETA,3)),f.mul(f.THETA,f.THETA))
    assert f.sign(f.sub(tr,f.scale(f.ONE,4)))>0
    image=f.mv(U,p);q=f.pair(p,G,p)
    assert f.sign(q)<0 and f.pair(image,G,image)==q
    assert f.sign(p[0])>0 and f.sign(image[0])>0
    walls=tuple(f.vecdecode(w) for w in d['all_walls'])
    for u in (p,image):assert all(f.sign(f.pair(u,G,r))<0 for r in walls)
    heights=tuple(f.sub(f.divide(f.mul(alpha,f.mul(u[0],u[0])),f.neg(q)),f.ONE) for u in (p,image))
    assert all(f.sign(f.sub(H,h))>0 for h in heights)
    # Explicitly place both endpoints in the same regular zero-height sector.
    zero=[]
    for start in ((f.ONE,f.ZERO),(f.neg(f.D),f.D)):
        v=start
        for _ in range(18):
            zero.append((f.ZERO,*v,f.ZERO));v=(f.sub(f.mul(f.THETA,v[0]),v[1]),v[0])
        assert v==start
    zero.extend(((f.ZERO,f.ZERO,f.ZERO,f.ONE),(f.ZERO,f.ZERO,f.ZERO,f.neg(f.ONE))))
    assert len(set(zero))==38
    for r in zero:
        assert f.sign(f.pair(p,G,r))==f.sign(f.pair(image,G,r))!=0
    return G,U,inverse,p,image,q,heights,tr

def run(verify=False):
    start=time.monotonic();d=json.loads(SOURCE.read_text());scout=json.loads(SCOUT.read_text())
    matrix=scout['successful_matrices'][0]['matrix']
    point=f.vecencode((f.ONE,f.scale(f.ONE,F(55,32)),f.scale(f.ONE,F(29,16)),f.scale(f.ONE,F(5,16))))
    G,U,inverse,p,image,q,heights,tr=local_check(d,matrix,point)
    bad=copy.deepcopy(matrix);bad[0][0]=f.encoded(f.add(f.decode(bad[0][0]),f.ONE))
    try:local_check(d,bad,point)
    except AssertionError:pass
    else:raise AssertionError('corrupt isometry was accepted')
    bad_point=f.vecencode((f.ONE,f.ZERO,f.ZERO,f.ZERO))
    try:local_check(d,matrix,bad_point)
    except AssertionError:pass
    else:raise AssertionError('nonregular endpoint was accepted')
    # Re-run the independent full prefix enumeration, reusing its live caches.
    audit_path=ROOT/'results/fresh_cubic9_h3_h3000_full_enumeration_audit_20261013.json'
    e.verify(SOURCE,audit_path)
    audit=json.loads(audit_path.read_text());assert audit['complete_spatial_census'] and audit['complete_time_norm_census']
    checks=[];total=0
    for row in d['levels']:
        x=f.decode(row['time']);norm=f.decode(row['norm']);roots,count=e.census.candidates(d,x,norm)
        assert len(roots)==row['candidate_roots'];root_checks=[]
        for r in sorted(roots):
            assert f.pair(r,G,r)==norm
            ps=(f.pair(p,G,r),f.pair(image,G,r))
            assert all(f.sign(t)<0 for t in ps)
            root_checks.append({'root':f.vecencode(r),'endpoint_pairings':f.vecencode(ps)})
        checks.append({'time':row['time'],'norm':row['norm'],'height':row['height'],'candidate_roots':root_checks});total+=len(roots)
    result={'status':'EXACT_FULL_HOST_NONREFLECTIVITY_CERTIFICATE',
        'host':'<-T*(4-T^2)> + [[2,-T],[-T,2]] + <1>, T^3-3T-1=0, T>1',
        'dimension':3,'source':str(SOURCE.relative_to(ROOT)),
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'complete_enumeration_audit':str(audit_path.relative_to(ROOT)),
        'matrix':matrix,'inverse_matrix':encode_matrix(inverse),'trace':f.encoded(tr),
        'point':point,'image':f.vecencode(image),'common_point_norm':f.encoded(q),
        'endpoint_height_bounds':f.vecencode(heights),'height_cutoff':3000,
        'complete_zero_oriented_roots':38,'complete_positive_controls':len(checks),
        'all_dominant_positive_candidates':total,'candidate_endpoint_checks':checks,
        'corrupt_isometry_rejected':True,'nonregular_endpoint_rejected':True,
        'scope':'Global nonreflectivity of this full integral H3 reflection group. The proof uses the complete primitive norm/zero alphabet, Cauchy separating-height bound, and finite zero-Weyl folding; every dominant root through the bound is checked at both endpoints. The forced-face H9 transfer is a separate mathematical corollary.'}
    if verify:assert json.loads(OUTPUT.read_text())==result
    else:OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'controls':len(checks),'dominant_roots':total,'seconds':round(time.monotonic()-start,3)}),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
