#!/usr/bin/env python3
"""Numerical five-old-edge seed and finite three-side closure diagnostic."""
import argparse,hashlib,itertools,json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from fresh_h20_two_anchor_families_20261016 import point,real_roots
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_open_edge_closure_20261016'
BASE='results/fresh_h20_two_anchor_families_20261016.json'

def run():
    D=json.loads((ROOT/BASE).read_text());G=geometry158();M=G['metric'];H=(G['H']@G['g']@G['F'].T)@M
    ids0=[1097,1993,2238,2483,1505];ff=[D['families'][i]for i in ids0]
    tt=[f['intervals'][0][0]*1.01 for f in ff]
    r0=point(ff[0],tt[0]);roots=real_roots(np.array(ff[1]['numerator'])@M@r0+np.array(ff[1]['denominator']))
    legal=[t for t in roots if t>ff[1]['intervals'][0][0]];assert len(legal)==1;tt[1]=legal[0]
    discovery=np.array([point(f,t)for f,t in zip(ff,tt)])
    # Stable promotion interface: four rational H pairings, then one exact
    # pi/3 pairing. The SVD kernel is only a numerical branch selector here.
    hf={0:-57208,2:-133898,3:-58604,4:-95777}
    support=[[151,95,142],[151,142,152],[151,145,152],[151,145,157],[151,95,157]]
    seed=discovery.copy();recipes=[]
    for i in [0,2,3,4,1]:
        U=np.r_[G['X'][support[i]],[H if i!=1 else seed[0]]]
        target=np.array([0.,0.,0.,hf[i]/1000 if i!=1 else -1.])
        B=U@M@U.T;part=target@np.linalg.inv(B)@U;ker=null_space(U@M)[:,0]
        sq=(2-part@M@part)/(ker@M@ker);assert sq>0
        branches=[part+s*np.sqrt(sq)*ker for s in [-1,1]]
        branch=min(range(2),key=lambda j:np.linalg.norm(branches[j]-discovery[i]))
        seed[i]=branches[branch]
        recipes.append(dict(wall=158+i,right_support=support[i],fourth='inherited_H'if i!=1 else'wall158',
            fourth_pairing=[hf[i],1000]if i!=1 else[-1,1],norm=[2,1],numerical_branch=branch,
            inherited_H_pairing=float(seed[i]@M@H),root=seed[i].tolist()))
    co=-seed@M@seed.T/2;np.fill_diagonal(co,0)
    assert gate_mask(co).all() and abs(co[0,1]-.5)<1e-7
    assert gate_mask(-(seed@M@G['X'].T)/np.sqrt(2*G['norms'])).all()
    assert max(abs(np.einsum('ij,jk,ik->i',seed,M,seed)-2))<1e-5
    assert max(seed@M@H)<0
    X=np.r_[G['X'],seed];q=np.r_[G['norms'],np.ones(5)*2]
    base,V,A,b,ids=face_data(G);seed_ar=-seed[:,1:]/seed[:,0,None]
    AA=np.r_[A,seed_ar@V];bb=np.r_[b,1-seed_ar@base];ii=ids+list(range(158,163))
    face=poly(base,V,AA,bb,ii)
    oldedges=[(95,142),(142,152),(145,152),(145,157),(95,157)]
    edgeproof=[]
    for j,k in oldedges:
        ends=[v['active']for v in face['vertices']if j in v['active']and k in v['active']and v['radius2']<1-1e-6]
        assert len(ends)==2;edgeproof.append(dict(edge=[j,k],finite_endpoints=ends))
    h=X@M@X[151];Xbar=X-h[:,None]*X[151]/2;qbar=np.einsum('ij,jk,ik->i',Xbar,M,Xbar)
    sides=np.where(qbar>1e-6)[0].tolist();assert sides==[95,142,145,150,152,156,157,158,159,160,161,162]
    openedges=[]
    for j,k in itertools.combinations(sides,2):
        g=Xbar[[j,k]]@M@Xbar[[j,k]].T
        if min(np.linalg.eigvalsh(g))<=1e-7:continue
        ends=[v['active']for v in face['vertices']if j in v['active']and k in v['active']and v['radius2']<1-1e-6]
        if len(ends)<2:openedges.append(dict(edge=[j,k],finite_endpoints=ends))
    labels=[2,3,4,5,6];patterns=np.array(list(itertools.product(labels,repeat=3)))
    cs=np.cos(np.pi/patterns);cs[patterns==2]=0
    stats=Counter();candidates=[];seen=set();framecounts=[]
    for triple in itertools.combinations(sides,3):
        if max(triple)<158:continue
        frame=[151,*triple];U=X[frame];B=U@M@U.T
        if min(abs(np.linalg.eigvalsh(B)))<1e-8:stats['singular_frame']+=1;continue
        inv=np.linalg.inv(B);ker=null_space(U@M)[:,0];kn=ker@M@ker
        target=-np.sqrt(2*q[frame])*np.c_[np.zeros(len(patterns)),cs]
        p=target@inv@U;rad=(2-np.einsum('ij,jk,ik->i',p,M,p))/kn
        real=np.where(rad>=-1e-8)[0];stats['patterns']+=len(patterns);start=len(candidates)
        for sign in [-1,1]:
            R=p[real]+sign*np.sqrt(np.maximum(rad[real],0))[:,None]*ker
            stats['real_branches']+=len(R)
            good=(R@M@H<0)&(R[:,0]>0)
            good &= gate_mask(-(R@M@X.T)/np.sqrt(2*q))
            for ri in np.where(good)[0]:
                r=R[ri]
                if abs(r@M@r-2)>1e-4:stats['numeric_norm_reject']+=1;continue
                key=tuple(np.round(r,7))
                if key in seen:continue
                seen.add(key)
                rb=r-(r@M@X[151])*X[151]/2;qr=rb@M@rb
                neighbors=[j for j in sides if qr*qbar[j]-(rb@M@Xbar[j])**2>1e-7]
                if len(neighbors)<3:stats['too_few_intrinsic_neighbors']+=1;continue
                vertices=[]
                for j,k in itertools.combinations(neighbors,2):
                    T=np.array([rb,Xbar[j],Xbar[k]])
                    if min(np.linalg.eigvalsh(T@M@T.T))>1e-7:vertices.append([j,k])
                ar=-r[1:]/r[0];cuts=[i for i,v in enumerate(face['vertices'])if ar@np.array(v['z'])>1+1e-7]
                if not cuts:stats['no_face_cut']+=1;continue
                candidates.append(dict(id=len(candidates),frame=frame,labels=[2,*map(int,patterns[real[ri]])],branch=sign,
                    root=r.tolist(),constraint=ar.tolist(),finite_neighbors=neighbors,new_finite_vertices=vertices,
                    old_face_vertex_cuts=cuts))
        framecounts.append(dict(frame=frame,accepted=len(candidates)-start))
    scored=[]
    chosen=sorted(candidates,key=lambda p:(-len(p['old_face_vertex_cuts']),-len(p['new_finite_vertices']),p['id']))[:64]
    for r in chosen:
        ar=np.array(r['constraint']);f=poly(base,V,np.r_[AA,[ar@V]],np.r_[bb,1-ar@base],ii+[163])
        scored.append(dict(candidate=r['id'],face=f))
    src=[f'scripts/{TAG}.py',BASE,'scripts/fresh_h20_two_anchor_families_20261016.py','scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py','scripts/fresh_h20_release144_cap_20261016_residual_scout.py']
    return dict(status='PASS_BOUNDED_OPEN_EDGE_SEED_AND_CLOSURE_SCOUT',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},
        seed=dict(discovery_families=ids0,discovery_parameters=tt,promotion_recipes=recipes,
            inherited_H=H.tolist(),inherited_H_norm=float(H@M@H),roots=seed.tolist(),
            mutual_cosines=co.tolist(),face=face,closed_old_edges=edgeproof,remaining_open_edges=openedges),
        side_domain=sides,label_domain=labels,stats=dict(stats),frames=framecounts,candidates=candidates,scored=scored,
        scope='Numerical five-cap seed closes all five old open intrinsic edges. Bounded subsequent three-side contacts, labels2..6, with at least one seed side; no arbitrary-label or continuous exclusion. Every candidate retains all163 old gates and has at least three potential finite neighbors. Exact wall promotion and full-face closure remain separate.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['seed']['remaining_open_edges']),len(d['candidates']),len(d['scored']),flush=True)
