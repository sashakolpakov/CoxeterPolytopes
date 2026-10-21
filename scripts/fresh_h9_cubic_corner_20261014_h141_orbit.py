#!/usr/bin/env python3
"""Close the new complete height141 prefix under the certified chamber involution."""
import argparse,json
from hashlib import sha256
import fresh_h9_cubic_corner_20261014 as P
K=P.K;Z,O=P.Z,P.O;ROOT=P.ROOT;TAG='fresh_h9_cubic_corner_20261014_h141_orbit'
def run():
    source='results/fresh_h9_cubic_corner_20261014_h141.json';tsrc='results/fresh_h9_cubic_corner_20261014_extension.json'
    d=json.loads((ROOT/source).read_text());td=json.loads((ROOT/tsrc).read_text())
    T=[[K.decode(x)for x in r]for r in td['isometry']];N=[[K.decode(x)for x in r]for r in d['normals']];a=K.decode(d['timelike_coefficient'])
    for r in N[:]:
        v=P.mv(T,r)
        if v not in N:N.append(v)
    h=[K.decode(x)for x in td['strict_center']]
    assert all(K.sign(K.pairing(h,r,a))<0 for r in N)
    keys={tuple(r):i for i,r in enumerate(N)};perm=[keys[tuple(P.mv(T,r))]for r in N];assert sorted(perm)==list(range(len(N)))
    C=K.sub(K.add(K.mul(K.THETA,K.THETA),K.THETA),O)
    finite=[Z,O,K.scale(O,2),K.scale(O,3),K.add(K.scale(O,2),K.THETA),K.add(K.scale(O,2),C)]
    norms=[K.pairing(r,r,a)for r in N];counts={'right':0,'finite_nonright':0,'dotted':0}
    for i,r in enumerate(N):
        assert norms[i]in(O,K.scale(O,2))
        for j,s in enumerate(N[:i]):
            p=K.pairing(r,s,a);assert K.sign(p)<=0;c=K.divide(K.scale(K.mul(p,p),4),K.mul(norms[i],norms[j]))
            if p==Z:counts['right']+=1
            elif K.sign(K.sub(c,K.scale(O,4)))>0:counts['dotted']+=1
            else:assert c in finite;counts['finite_nonright']+=1
    paths=[source,tsrc,'scripts/fresh_h9_cubic_corner_20261014.py',f'scripts/{TAG}.py']
    return P.encode(dict(status='EXACT_HEIGHT141_PREFIX_AND_FULL_CHAMBER_SYMMETRY_ORBIT',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        normals=N,wall_count=len(N),timelike_coefficient=a,isometry=T,wall_permutation=perm,strict_center=h,pair_counts=counts,
        complete_distance_bound=141,scope='All walls are actual walls of the same full arithmetic chamber. The two complete distance neighborhoods have bound141. No finite-volume conclusion.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['wall_count'])
