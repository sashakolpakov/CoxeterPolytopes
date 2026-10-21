#!/usr/bin/env python3
"""Certify a full-chamber symmetry at a near fixed point; transport actual walls."""
import argparse
import json
from hashlib import sha256
from pathlib import Path
import fresh_h9_cubic_corner_20261014 as P
K=P.K;Z,O=P.Z,P.O;ROOT=P.ROOT;TAG='fresh_h9_cubic_corner_20261014_extension'

def run():
    source='results/fresh_cubic7_h9_seed1_h100_20261013.json'
    mapsource='results/fresh_h9_cubic_corner_20261014.json'
    old=json.loads((ROOT/source).read_text());data=json.loads((ROOT/mapsource).read_text())
    T=[[K.decode(x)for x in r]for r in data['seeds'][1]['integral_sector_maps'][1]['matrix']]
    N=[[K.decode(x)for x in r]for r in old['all_walls']];a=K.decode(old['timelike_coefficient'])
    I=[[O if i==j else Z for j in range(10)]for i in range(10)]
    assert P.mm(T,T)==I and P.integral(T)
    u=[O]+[Z]*9;Tu=P.mv(T,u);H=[K.add(x,y)for x,y in zip(u,Tu)]
    assert P.mv(T,H)==H and K.sign(H[0])>0
    norm=K.pairing(H,H,a);assert K.sign(norm)<0
    height=K.sub(K.divide(K.mul(a,K.mul(H[0],H[0])),K.neg(norm)),O)
    assert K.sign(K.sub(K.scale(O,100),height))>0
    ips=[K.pairing(H,r,a)for r in N];assert all(K.sign(x)<=0 for x in ips)
    incident=[i for i,x in enumerate(ips)if x==Z];keys={tuple(r):i for i,r in enumerate(N)}
    perm=[keys[tuple(P.mv(T,N[i]))]for i in incident];assert sorted(perm)==incident
    # Complete Vinberg height100 plus height(H)<100 proves all full chamber
    # walls incident to H are already present. A lattice isometry fixing H
    # and positively permuting those walls preserves its local chamber sector.
    # It therefore preserves the full connected reflection chamber.
    transported=[P.mv(T,r)for r in N];allroots=N[:]
    for r in transported:
        if r not in allroots:allroots.append(r)
    norms=[K.pairing(r,r,a)for r in allroots]
    assert all(q in(O,K.scale(O,2))for q in norms)
    C=K.sub(K.add(K.mul(K.THETA,K.THETA),K.THETA),O)
    finite=[Z,O,K.scale(O,2),K.scale(O,3),K.add(K.scale(O,2),K.THETA),K.add(K.scale(O,2),C)]
    counts={'right':0,'finite_nonright':0,'dotted':0}
    for i,r in enumerate(allroots):
        for j,s in enumerate(allroots[:i]):
            p=K.pairing(r,s,a);assert K.sign(p)<=0
            c=K.divide(K.scale(K.mul(p,p),4),K.mul(norms[i],norms[j]))
            if p==Z:counts['right']+=1
            elif K.sign(K.sub(c,K.scale(O,4)))>0:counts['dotted']+=1
            else:assert c in finite;counts['finite_nonright']+=1
    h=[K.decode(x)for x in old['strict_interior']]
    assert all(K.sign(K.pairing(h,r,a))<0 for r in allroots)
    oldescape=[K.decode(x)for x in old['geometry']['vector']]
    escape=[K.pairing(oldescape,r,a)for r in allroots]
    # No volume conclusion is made if a new escape has not been certified.
    preserved=all(K.sign(x)<=0 for x in escape)
    paths=[source,mapsource,'scripts/fresh_h9_cubic_corner_20261014.py',
           'scripts/audit_fresh_cubic7_h9_no_import_20261013.py',f'scripts/{TAG}.py']
    return P.encode(dict(status='EXACT_FIXED_POINT_FULL_CHAMBER_SYMMETRY_AND_WALL_TRANSPORT',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},seed=1,
        isometry=T,fixed_point=H,fixed_point_norm=norm,fixed_point_height=height,
        incident_original_indices=incident,incident_permutation=perm,
        old_wall_count=len(N),wall_count=len(allroots),normals=allroots,pair_counts=counts,
        strict_center=h,old_escape_preserved=preserved,old_escape_pairings=escape,
        timelike_coefficient=a,
        scope='The finite-order integral symmetry is certified through its near fixed-point sector, using the inherited complete Vinberg height100 prefix. Every transported wall is a wall of the same full arithmetic chamber. The distant corner is thereby validated by transport. This alone proves neither reflectivity nor nonreflectivity.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'walls',out['wall_count'],'old escape survives',out['old_escape_preserved'])
