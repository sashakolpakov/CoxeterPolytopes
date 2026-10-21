#!/usr/bin/env python3
"""Narrow saved-packet and transport audit; no new candidate enumeration."""
import argparse,json
from hashlib import sha256
from pathlib import Path
from fresh_h9_cubic_corner_20261014_review import C,ZERO,ONE,matrix,vector,mv
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_cubic_corner_20261014_h141_review'

def run():
    names=['results/fresh_h9_cubic_corner_20261014_h141.json',
           'results/fresh_h9_cubic_corner_20261014_h141_orbit.json',
           'results/fresh_cubic7_h9_seed1_h100_20261013.json',
           'results/fresh_h9_cubic_corner_20261014_review.json']
    new,orb,old,review=[json.loads((ROOT/p).read_text())for p in names]
    for d in (new,orb,review):
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert new['height_bound']==orb['complete_distance_bound']==141 and new['previous_height_bound']==100
    assert new['complete_new_levels']==len(new['levels'])==2
    a=C.from_json(new['timelike_coefficient']);assert a==2*C(0,1)
    Q=[[ZERO for _ in range(10)]for _ in range(10)];Q[0][0]=-a;Q[9][9]=ONE
    for i in range(1,9):Q[i][i]=C(2)
    for i,j in ((1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(3,8)):Q[i][j]=Q[j][i]=-ONE
    ip=lambda x,y:sum((u*v for u,v in zip(x,mv(Q,y))),ZERO)
    N=matrix(old['all_walls']);assert len(N)==17
    heights=[];packet_sizes=[]
    for row in new['levels']:
        q=C.from_json(row['norm']);x=C.from_json(row['control']);h=C.from_json(row['height'])
        assert q in (ONE,C(2)) and h==a*x*x/q and 100<h<=141
        if heights:assert heights[-1]<h
        heights.append(h);packet=matrix(row['new_walls']);packet_sizes.append(len(packet))
        assert row['candidates']>=len(packet)
        for r in packet:
            assert r[0]==x and ip(r,r)==q and all(t.is_integral for t in r)
            assert all((2*t/q).is_integral for t in mv(Q,r))
            assert all(ip(r,s)<=0 for s in N)
        assert all(ip(r,s)<=0 for i,r in enumerate(packet)for s in packet[:i])
        N+=packet;assert len(N)==row['wall_count']
    assert packet_sizes==[5,1] and N==matrix(new['normals']) and len(N)==23
    T=matrix(orb['isometry']);full=N[:]
    for r in N:
        im=mv(T,r)
        if im not in full:full.append(im)
    assert full==matrix(orb['normals']) and len(full)==32
    assert [full.index(mv(T,r))for r in full]==orb['wall_permutation']
    h=vector(orb['strict_center']);assert ip(h,h)<0 and all(ip(h,r)<0 for r in full)
    norms=[ip(r,r)for r in full];allowed={ZERO,ONE,C(2),C(3),C(2,1),C(1,1,1)}
    counts=dict(right=0,finite_nonright=0,dotted=0)
    for i,r in enumerate(full):
        assert norms[i]in (ONE,C(2))
        for j,s in enumerate(full[:i]):
            p=ip(r,s);assert p<=0
            k=4*p*p/(norms[i]*norms[j])
            if not p:counts['right']+=1
            elif k>4:counts['dotted']+=1
            else:assert k in allowed;counts['finite_nonright']+=1
    assert counts==orb['pair_counts']==dict(right=200,finite_nonright=36,dotted=260)
    paths=names+['scripts/fresh_h9_cubic_corner_20261014_h141.py',
                 'scripts/fresh_h9_cubic_corner_20261014_h141_orbit.py',
                 'scripts/fresh_h9_cubic_corner_20261014_review.py',f'scripts/{TAG}.py']
    return dict(status='PASS_SCOPED_HEIGHT141_SOURCE_AND_SAVED_PACKET_REVIEW',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                new_control_levels=2,level_heights=[h.to_json()for h in heights],new_packet_sizes=packet_sizes,
                raw_wall_count=23,transported_wall_count=32,pair_counts=counts,
                scope='Independent saved-packet norms/reflection-integrality/chronology and transported pair arithmetic in Cubic7. Source review of regenerated old-control matching and complete new-level enumeration. Does not independently repeat the candidate/control census or the separate finite-vertex/edge enumeration.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();dest=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['new_packet_sizes'],d['transported_wall_count'])
