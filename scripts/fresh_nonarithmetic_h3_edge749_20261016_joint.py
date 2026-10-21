#!/usr/bin/env python3
"""Finite-incidence screen in the recentered H3 Lorentz chart."""
import argparse,json
from collections import Counter
from hashlib import sha256
import numpy as np
import fresh_nonarithmetic_h3_a2_profiles_20261016_joint as J
C=J.C;L=J.L;ETA=J.ETA;ROOT=J.ROOT;BASE='fresh_nonarithmetic_h3_edge749_20261016';TAG=BASE+'_joint'

def build():
    d=json.loads((ROOT/f'results/{BASE}.json').read_text());g=d['recentered_geometry'];N=np.array(g['old_projected_normals']);pool={};skipped=0
    for pr in d['profiles']:
        co=np.array([L.product(pr['label'])if node==pr['node']else 0. for node in(5,6)])
        q=2-co@J.J2@co
        for k,r in enumerate(pr['candidates']):
            v=np.array(r['vector'])
            if abs(v@ETA@v-q)>6e-5:skipped+=1;continue
            key=(tuple(np.round(v,7)),tuple(np.round(co,10)))
            pool.setdefault(key,dict(vector=v,A2_contacts=co,residual_square=q,
                origin=dict(source='recentered',node=pr['node'],label=pr['label'],candidate=k,recipe=r['recipe'])))
    rows=list(pool.values());V=np.array([r['vector']for r in rows]);co=np.array([r['A2_contacts']for r in rows]);q=np.array([r['residual_square']for r in rows])
    rays=[dict(mask=sum(1<<i for i in r['active']),v=np.array(r['point']),kind=r['kind'])for r in g['face_rays']]
    edges,ends,sides=L.ledger(rays,N);points=np.array([e['point']for e in ends]);oldq=np.einsum('ij,jk,ik->i',N,ETA,N)
    oldfinite=q[:,None]*oldq[sides]-(V@ETA@N[sides].T)**2>2e-6;olddegree=oldfinite.sum(axis=1);adj=[];total=0
    for start in range(0,len(rows),256):
        pair=V[start:start+256]@ETA@V.T;native=pair+co[start:start+256]@J.J2@co.T
        ok=J.fits(native)&(q[start:start+256,None]*q[None,:]-pair*pair>2e-6)
        for k,line in enumerate(ok,start):
            line[k]=False;js=np.flatnonzero(line).tolist();adj.append(js);total+=len(js)
    assert total%2==0
    alive=np.ones(len(rows),dtype=bool);rounds=[]
    while True:
        degree=np.array([int(olddegree[i])+sum(alive[j]for j in adj[i])for i in range(len(rows))]);remove=alive&(degree<3)
        if not remove.any():break
        rounds.append(int(remove.sum()));alive[remove]=False
    core=np.flatnonzero(alive).tolist();stats=dict(raw_candidates=len(rows),uncertain_norm_skipped=skipped,finite_new_edges=total//2,
        old_degree_counts={str(k):int(v)for k,v in Counter(map(int,olddegree)).items()},pruning_rounds=rounds,core_count=len(core),
        core_profile_counts=dict(Counter(str((rows[i]['origin']['node'],rows[i]['origin']['label']))for i in core)))
    print('RECENTERED_CORE',stats,flush=True)
    return rows,V,co,q,N,rays,points,olddegree,adj,core,stats

def run():
    J.build=build;d=J.run();paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',f'scripts/{J.TAG}.py',f'results/{J.TAG}.json']
    d['sources']={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths}
    d['status']='BOUNDED_NUMERICAL_RECENTERED_H3_JOINT_INCIDENCE_SCREEN'
    d['scope']='Finite floating-point pool of all seven A2 profiles in the recentered chart, retaining the newly checked common center. Native mutual projection offsets and intrinsic finite-neighbor degrees are enforced. No exact geometry, universal continuous obstruction, or finite-volume claim.'
    return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved recentered joint screen')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['joint_cover']['status'],flush=True)
