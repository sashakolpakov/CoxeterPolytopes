#!/usr/bin/env python3
"""Numerical joint finite-incidence screen of zero/nonzero A2 profiles."""
import argparse,json,itertools,math
from collections import Counter
from hashlib import sha256
import numpy as np
import mpmath as mp
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import csr_matrix
import fresh_nonarithmetic_h3_a2_profiles_20261016 as A
C=A.C;B=A.B;L=A.L;ROOT=A.ROOT;ETA=A.ETA;TAG=A.TAG+'_joint'
J2=np.array([[2.,1.],[1.,2.]])/3

def fits(g):
    angle=np.pi/np.arccos(np.clip(-g/2,-1+1e-15,1-1e-15))
    return(g<-2-4e-5)|(abs(g)<4e-5)|((g<-1+4e-5)&(g>-2+4e-5)&(abs(angle-np.rint(angle))<4e-5))

def build():
    d=json.loads((ROOT/f'results/{C.TAG}.json').read_text());new=json.loads((ROOT/f'results/{A.TAG}.json').read_text());N=np.array(d['old_projected_normals'])
    blocked=set(json.loads((ROOT/f'results/{C.TAG}_affine_gate.json').read_text())['excluded_family_indices']);pool={};counts=Counter()
    def admit(v,contacts,origin):
        contacts=np.array(contacts);q=2-contacts@J2@contacts
        if abs(v@ETA@v-q)>6e-5:counts['uncertain_norm_skipped']+=1;return
        key=(tuple(np.round(v,7)),tuple(np.round(contacts,10)))
        pool.setdefault(key,dict(vector=v,A2_contacts=contacts,residual_square=q,origin=origin))
    for f in d['families']:
        if f['index']in blocked:continue
        for rec in f['intervals']:
            for s in rec['samples']:admit(np.array(s['vector']),[0,0],dict(source='two_contact',family=f['index'],parameter=s['parameter']))
    for pr in new['profiles']:
        contacts=[A.L.product(pr['label'])if node==pr['node']else 0. for node in(5,6)]
        for k,cand in enumerate(pr['candidates']):admit(np.array(cand['vector']),contacts,dict(source='nonzero_A2',node=pr['node'],label=pr['label'],candidate=k,recipe=cand['recipe']))
    # The old bank was constructed with six trial walls. Their coordinates
    # are used to reconstruct its recipes, without retaining those six walls.
    trial=json.loads((ROOT/'results/fresh_nonarithmetic_h3_closure_20261016_trial.json').read_text())
    Q,unused,H,E=C.P.geometry();basis=[mp.matrix([mp.mpf(x)for x in v])for v in d['basis']]
    def coord(v):return np.array([float(-(v.T*Q*basis[0])[0])]+[float((v.T*Q*z)[0])for z in basis[1:]])
    extra=np.array([coord(mp.matrix([mp.mpf(x)for x in v]))for v in trial['normal_intervals_or_proposals']]);allN=np.vstack([N,extra])
    bank=json.loads((ROOT/'results/fresh_nonarithmetic_h3_joint_cover_20261016.json').read_text())
    for k,r in enumerate(bank['pool']):
        anchors=[i+2 if i>=53 else i for i in r['anchors']];frame=allN[anchors];Ji=np.linalg.inv(frame@ETA@frame.T)
        c=np.array([L.product(m)for m in r['labels']]);p=c@Ji@frame;h=np.array([1.,0,0,0]);u=h-(h@ETA@frame.T)@Ji@frame
        rho=(2-p@ETA@p)/(u@ETA@u);assert rho>0;v=p+r['branch']*math.sqrt(rho)*u
        assert v[0]>0 and all(C.compatible(g,4e-5)for g in N@ETA@v)
        admit(v,[0,0],dict(source='old4329',index=k))
    rows=list(pool.values());V=np.array([r['vector']for r in rows]);co=np.array([r['A2_contacts']for r in rows]);q=np.array([r['residual_square']for r in rows])
    rays=[dict(mask=sum(1<<i for i in r['active']),v=np.array(r['point']),kind=r['kind'])for r in d['actual_face_rays']]
    edges,ends,sides=L.ledger(rays,N);points=np.array([r['point']for r in ends]);oldq=np.einsum('ij,jk,ik->i',N,ETA,N)
    oldfinite=q[:,None]*oldq[sides]-(V@ETA@N[sides].T)**2>2e-6;olddegree=oldfinite.sum(axis=1)
    adj=[];total=0
    for start in range(0,len(rows),256):
        g=V[start:start+256]@ETA@V.T;native=g+co[start:start+256]@J2@co.T
        ok=fits(native)&(q[start:start+256,None]*q[None,:]-g*g>2e-6)
        for k,line in enumerate(ok,start):
            line[k]=False;js=np.flatnonzero(line).tolist();adj.append(js);total+=len(js)
        if start%4096==0:print('FINITE_GRAPH',start,len(rows),flush=True)
    assert total%2==0
    alive=np.ones(len(rows),dtype=bool);rounds=[]
    while True:
        degree=np.array([int(olddegree[i])+sum(alive[j]for j in adj[i])for i in range(len(rows))])
        remove=alive&(degree<3)
        if not remove.any():break
        rounds.append(int(remove.sum()));alive[remove]=False
    core=np.flatnonzero(alive).tolist();print('INCIDENCE_CORE',len(rows),len(core),'edges',total//2,'rounds',rounds,flush=True)
    stats=dict(raw_candidates=len(rows),candidate_origin_counts=dict(Counter(r['origin']['source']for r in rows)),old_degree_counts={str(k):int(v)for k,v in Counter(map(int,olddegree)).items()},
        finite_new_edges=total//2,pruning_rounds=rounds,core_count=len(core),core_origin_counts=dict(Counter(rows[i]['origin']['source']for i in core)),extra_counts=dict(counts))
    return rows,V,co,q,N,rays,points,olddegree,adj,core,stats

def run():
    rows,V,co,q,N,rays,points,olddegree,adj,core,stats=build();cv=V[core];cc=co[core];position={j:i for i,j in enumerate(core)}
    matrix=[];lo=[];hi=[]
    coverage=cv@ETA@points.T>2e-7
    stats['raw_initial_edge_end_cover_counts']=[int(x)for x in (V@ETA@points.T>2e-7).sum(axis=0)]
    stats['core_initial_edge_end_cover_counts']=[int(x)for x in coverage.sum(axis=0)]
    stats['uncovered_core_initial_edge_ends']=np.flatnonzero(~coverage.any(axis=0)).tolist()
    for p in range(len(points)):
        ids=np.flatnonzero(coverage[:,p]);matrix.append({int(i):1 for i in ids});lo.append(1);hi.append(np.inf)
    for j,i in enumerate(core):
        if olddegree[i]>=3:continue
        row={j:float(3-olddegree[i])}
        for k in adj[i]:
            if k in position:row[position[k]]=-1
        matrix.append(row);lo.append(-np.inf);hi.append(0)
    for start in range(0,len(core),256):
        ok=fits(cv[start:start+256]@ETA@cv.T+cc[start:start+256]@J2@cc.T)
        for j,line in enumerate(ok,start):
            for k in np.flatnonzero(~line[j+1:])+j+1:
                matrix.append({j:1,int(k):1});lo.append(-np.inf);hi.append(1)
    print('JOINT_ROWS',len(matrix),'variables',len(core),flush=True)
    rr=[];cs=[];vv=[]
    for i,row in enumerate(matrix):
        for j,z in row.items():rr.append(i);cs.append(j);vv.append(z)
    mat=csr_matrix((vv,(rr,cs)),shape=(len(matrix),len(core)))
    ans=milp(c=np.ones(len(core)),integrality=np.ones(len(core)),bounds=Bounds(0,1),constraints=LinearConstraint(mat,lo,hi),options=dict(time_limit=60,mip_rel_gap=0))
    result=dict(status=int(ans.status),message=ans.message,constraint_rows=len(matrix),selected=None)
    if ans.x is not None:
        selected=[core[j]for j in np.flatnonzero(ans.x>.5)];result['selected']=selected
        current=rays
        for wall,i in enumerate(selected,55):current=B.clip(current,V[i],wall)
        result['face_counts']=dict(Counter(r['kind']for r in current));result['final_rays']=[dict(active=[j for j in range(r['mask'].bit_length())if r['mask']>>j&1],point=r['v'],kind=r['kind'])for r in current]
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',f'scripts/{L.TAG}.py',f'results/{C.TAG}.json',f'results/{C.TAG}_affine_gate.json',
        'results/fresh_nonarithmetic_h3_closure_20261016_trial.json','results/fresh_nonarithmetic_h3_joint_cover_20261016.json']
    return C.clean(dict(status='BOUNDED_NUMERICAL_A2_PROFILE_JOINT_INCIDENCE_SCREEN',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        counts=stats,core=[dict(index=i,**rows[i],old_finite_degree=int(olddegree[i]),finite_new_neighbors=[j for j in adj[i]if j in position],initial_edge_coverage=np.flatnonzero(V[i]@ETA@points.T>2e-7).tolist())for i in core],
        joint_cover=result,scope='Finite floating-point pool only. Includes all sampled frozen zero-profile conics, six new A2 profiles and the old4329 bank, with native offset corrections in all new-new products. Degree>=3 is necessary for every actual polygon facet of a compact H3, but insufficient. No exact or all-continuous nonexistence claim. Numerical affine exclusions are not complete; any surviving proposal still requires exact cusp screening and geometry.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved A2 joint screen')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts']['core_count'],d['joint_cover']['status'],flush=True)
