#!/usr/bin/env python3
"""Exact chamber-stabilizer certificate, using separate Fraction arithmetic.

The two finite segment censuses are independently replayed by the cited
checkers.  This script anchors their source hashes, joins their oriented
sectors, and constructs integral corner involutions directly from the
incident Gram.  It does not infer chamber preservation from an arbitrary
finite-prefix interior.
"""
import argparse,hashlib,json
import fresh_h9_new_host_20261013_independent as a
from fresh_h9_quartic_first_endpoint_20261013_independent import closure,line

ROOT=a.ROOT
FIRST=ROOT/'results/fresh_h9_new_host_20261013_first.json'
FIRST_AUDIT=ROOT/'results/fresh_h9_new_host_20261013_first_independent.json'
NEXT=ROOT/'results/fresh_h9_new_host_next_20261013_census_e1_c64.json'
NEXT_AUDIT=ROOT/'results/fresh_h9_new_host_next_20261013_e1_independent.json'
BASE=ROOT/'results/fresh_h9_new_host_20261013_corrected.json'
OUT=ROOT/'results/fresh_h9_new_host_next_20261013_global.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def identity(n):return tuple(tuple(a.O if i==j else a.Z for j in range(n))for i in range(n))
def mm(A,B):return tuple(tuple(a.dot(row,col)for col in zip(*B))for row in A)
def enc(M):return [[a.enc(v)for v in row]for row in M]
def lift(r):return(r[0],*r[1:5],*(a.Z,)*4,r[5])
def liftmatrix(M):
    ids=[0,1,2,3,4,9];A=[list(row)for row in identity(10)]
    for i,u in enumerate(ids):
        for j,v in enumerate(ids):A[u][v]=M[i][j]
    return tuple(map(tuple,A))
def involution(M,g):
    assert all(a.integral(z)for row in M for z in row)
    assert mm(M,M)==identity(len(M))
    assert mm(mm(tuple(zip(*M)),g),M)==g

def run(verify=False):
    f,n,b=[json.loads(p.read_text())for p in(FIRST,NEXT,BASE)]
    fa,na=[json.loads(p.read_text())for p in(FIRST_AUDIT,NEXT_AUDIT)]
    assert fa['status']=='PASS_INDEPENDENT_FULL_H9_FIRST_EDGE'and fa['source_sha256']==sha(FIRST)
    assert na['status']=='PASS_INDEPENDENT_SHIFTED_PAIRING_CENSUS'and len(na['cases'])==1
    nr=na['cases'][0];assert nr['source_sha256']==sha(NEXT)and nr['complete_H9_lift']
    assert n['status']=='COMPLETE_FIRST_EDGE_MIRROR'and nr['status']=='PASS_COMPLETE_FIRST_EDGE_MIRROR'
    assert n['source_sha256']==sha(BASE)and n['H9_lift_source_sha256']==sha(FIRST)
    for p,h in f['sources'].items():assert sha(ROOT/p)==h
    for p,h in b['sources'].items():assert sha(ROOT/p)==h
    host,h9=a.Host(1),a.Host(2);I=identity(6)
    walls=[a.vec(r)for r in b['forced_H5_face']['all_walls']]
    J=tuple(tuple(a.O if j==[0,4,3,2,1,5][i]else a.Z for j in range(6))for i in range(6))
    involution(J,host.g)
    origin=a.vec(f['start']);p=a.vec(n['start']);endpoint=a.vec(n['endpoint']['point'])
    assert origin==(a.O,*(a.Z,)*5)and a.mv(J,origin)==origin
    # J preserves the initial oriented simple roots, hence their full chamber.
    initial=walls[:5];jperm=[3,2,1,0,4]
    assert [a.mv(J,r)for r in initial]==[initial[i]for i in jperm]
    assert len(closure(initial,initial,host))==25
    assert a.mv(J,a.vec(f['end']))==p
    first_incident=[a.vec(r)for r in f['edge_roots']]+[a.vec(f['selected_choice']['vector'])]
    start_ids=n['edge_indices']+[n['released_wall']]
    assert set(a.mv(J,r)for r in first_incident)==set(walls[i]for i in start_ids)
    assert a.mv(J,a.vec(f['selected_choice']['vector']))==walls[8]
    assert nr['start_mirrors']==20 and nr['endpoint_mirrors']==33 and nr['H9_endpoint_mirrors']==57
    edge=[walls[i]for i in n['edge_indices']];new=a.vec(n['endpoint']['root']);incident=edge+[new]
    assert all(host.pair(endpoint,r)==a.Z for r in incident)
    assert a.sign(host.pair(p,new))<0
    # Exact old-wall inequalities and outward orientation on the whole segment.
    assert all(a.sign(host.pair(endpoint,r))<=0 for r in walls)
    assert all(host.compatible(r,incident[:i])for i,r in enumerate(incident))
    gram=tuple(tuple(host.pair(r,s)for s in incident)for r in incident)
    a.positive_definite(gram)
    assert gram==tuple(a.vec(row)for row in n['endpoint']['incident_gram'])
    assert len(closure(incident,incident,host))==33
    perm=[1,0,2,4,3];assert tuple(tuple(gram[perm[i]][perm[j]]for j in range(5))for i in range(5))==gram
    gi=a.minverse(gram);cov=[a.mv(host.g,r)for r in incident]
    # Fix the timelike orthogonal complement and permute the spacelike frame.
    H=tuple(tuple(a.add(I[i][j],a.dot(tuple(a.sub(incident[perm[k]][i],incident[k][i])for k in range(5)),
       a.mv(gi,tuple(cov[l][j]for l in range(5)))))for j in range(6))for i in range(6))
    involution(H,host.g)
    assert a.mv(H,endpoint)==endpoint
    assert [a.mv(H,r)for r in incident]==[incident[i]for i in perm]
    assert a.sign(host.pair(endpoint,endpoint))<0 and a.sign(endpoint[0])>0
    assert a.sign(a.mv(H,origin)[0])>0 # future-preserving; J fixes origin.
    product=mm(H,J);trace=a.Z
    for i in range(6):trace=a.add(trace,product[i][i])
    assert trace==a.dec([-56022,-30320,95640,51760])and a.sign(a.sub(trace,a.scale(a.O,6)))>0
    J9,H9=liftmatrix(J),liftmatrix(H)
    involution(J9,h9.g);involution(H9,h9.g)
    other=[tuple(a.neg(a.O)if i==j else a.Z for i in range(10))for j in range(5,9)]
    assert all(a.mv(J9,r)==r and a.mv(H9,r)==r for r in other)
    groups=[]
    for gen,expected in((initial,49),([walls[i]for i in start_ids],44),(incident,57)):
        full=[*map(lift,gen),*other];roots=closure(full,full,h9);assert len(roots)==expected
        groups.append(len(roots))
    assert a.mv(H9,lift(endpoint))==lift(endpoint)
    assert [a.mv(H9,lift(r))for r in incident]==[lift(incident[i])for i in perm]
    # A local inward direction identifies the oriented chamber sector, rather
    # than just its unoriented collection of mirrors.
    coeff=a.mv(gi,tuple(a.neg(a.O)for _ in incident))
    inward=tuple(a.dot(coeff,tuple(r[j]for r in incident))for j in range(6))
    assert all(host.pair(inward,r)==a.neg(a.O)for r in incident)
    assert a.mv(H,inward)==inward
    result={'status':'EXACT_FULL_H9_NONREFLECTIVITY_CERTIFICATE',
      'sources':{str(p.relative_to(ROOT)):sha(p)for p in(FIRST,FIRST_AUDIT,NEXT,NEXT_AUDIT,BASE)},
      'alpha':a.enc(a.ALPHA),'gram_H5':enc(host.g),'gram_H9':enc(h9.g),
      'initial_F4_reversal':enc(J),'initial_oriented_permutation':jperm,
      'endpoint':list(map(a.enc,endpoint)),'endpoint_norm':a.enc(host.pair(endpoint,endpoint)),
      'endpoint_incident_roots':enc(incident),'endpoint_incident_gram':enc(gram),
      'endpoint_oriented_permutation':perm,'endpoint_integral_involution':enc(H),
      'endpoint_inward_tangent':list(map(a.enc,inward)),
      'initial_inverse':enc(J),'endpoint_inverse':enc(H),'product':enc(product),
      'product_trace_H5':a.enc(trace),'product_trace_H9':a.enc(a.add(trace,a.scale(a.O,4))),
      'trace_H5_greater_than_6':True,'integral_both_directions':True,'future_preserving':True,
      'H9_full_corner_mirror_counts':groups,'same_sector_path_verified':True,
      'H9_identity_extension_verified':True,
      'scope':'The complete all-norm first-edge censuses identify one full arithmetic reflection chamber from the initial F4+A1 corner, through its J-reflected I2(16)+A2+A1 corner, to an I2(16)^2+A1 corner. J and H are integral lattice involutions preserving the oriented local sector at corners of this same chamber. Their product has trace>6 on H5 and trace>10 after identity extension to H9, hence infinite order. They stabilize the full reflection chamber, whose volume is therefore infinite. This excludes reflectivity and finite-volume fundamental Coxeter chambers of these two lattices; it is not a no-go theorem for other arithmetic forms or arbitrary Coxeter polytopes.'}
    if verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],'H9 corner mirrors',groups,'trace',a.enc(trace),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
