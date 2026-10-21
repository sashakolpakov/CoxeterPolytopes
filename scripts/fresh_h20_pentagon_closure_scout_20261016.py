#!/usr/bin/env python3
"""Bounded simultaneous wall search on S16+{151,156}; numerical only."""
import argparse, hashlib, itertools, json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection
from fresh_h20_release144_cap_20261016_residual_scout import geometry156, recipe, constraint, close

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_pentagon_closure_scout_20261016'

def geometry157():
    G = geometry156()
    r = recipe(G, 3)
    G.update(N=np.r_[G['N'], [r]], norms=np.r_[G['norms'], 2],
             A=np.r_[G['A'], [constraint(G, r)]], b=np.r_[G['b'], 1],
             labels=G['labels']+[156])
    # Orthonormal Lorentz coordinates avoid solving in the near-null old chart.
    F = np.r_[[G['e0']], G['ee']]
    metric = np.diag([-1., 1., 1., 1., 1.])
    X = (G['N'] @ G['g'] @ F.T) @ metric
    assert np.max(abs(F @ G['g'] @ F.T-metric)) < 1e-7
    G.update(F=F, X=X, metric=metric)
    return G

def relation(c):
    if abs(c) < 2e-6: return 'right'
    if c > 1+2e-6: return 'dotted'
    if abs(c-1) < 2e-6: return 'near_parallel'
    if c < 0: return 'obtuse'
    m = np.pi/np.arccos(c)
    if m >= 2.9999 and abs(m-round(m)) < 2e-4: return 'finite'
    return 'nonintegral'

def face_chart(G):
    idx = [G['labels'].index(i) for i in [151, 156]]
    eq = G['A'][idx]
    base = eq.T @ np.linalg.solve(eq @ eq.T, G['b'][idx])
    V = null_space(eq)
    assert V.shape == (4, 2)
    A = G['A'] @ V
    b = G['b']-G['A'] @ base
    live = np.linalg.norm(A, axis=1) > 1e-8
    assert np.min(b[~live]) > -1e-7
    A, b = A[live], b[live]
    lengths = np.linalg.norm(A, axis=1)
    A, b = A/lengths[:,None], b/lengths
    lp = linprog([0,0,-1], A_ub=np.c_[A,np.ones(len(A))], b_ub=b,
                 bounds=[(None,None)]*2+[(0,None)], method='highs')
    assert lp.success and lp.x[2]>1e-9
    return base, V, A, b, lp.x[:2]

def polygon(base, V, A, b):
    lengths=np.linalg.norm(A,axis=1)
    assert min(lengths)>1e-10
    aa,bb=A/lengths[:,None],b/lengths
    lp=linprog([0,0,-1], A_ub=np.c_[aa,np.ones(len(A))], b_ub=bb,
               bounds=[(None,None)]*2+[(0,None)],method='highs')
    assert lp.success and lp.x[2]>1e-10
    yy=HalfspaceIntersection(np.c_[aa,-bb],lp.x[:2]).intersections
    unique=[]
    for y in yy:
        if not any(np.linalg.norm(y-t)<1e-8 for t in unique):unique.append(y)
    yy=np.array(unique);zz=base+yy@V.T
    return dict(vertices=zz.tolist(), squared_radii=np.sum(zz*zz,axis=1).tolist(),
                finite=int(np.sum(np.sum(zz*zz,axis=1)<1-1e-6)),
                near_sphere=int(np.sum(abs(np.sum(zz*zz,axis=1)-1)<=1e-6)),
                exterior=int(np.sum(np.sum(zz*zz,axis=1)>1+1e-6)))

def run(max_label):
    G=geometry157();X,M=G['X'],G['metric'];norms=G['norms']
    base,V,A,b,center=face_chart(G)
    old=polygon(base,V,A,b)
    assert (old['finite'],old['exterior'])==(2,3)
    angles=range(2,max_label+1)
    cs={m:(0. if m==2 else np.cos(np.pi/m)) for m in angles}
    support_labels=[(a,b) for a in angles for b in angles if cs[a]**2+cs[b]**2<1-1e-9]
    sides=[31,95,140,142,145]
    H=(G['H']@G['g']@G['F'].T)@M
    failures=Counter();candidates=[];trials=0;real_branches=0
    oldpoints=np.array(old['vertices'])
    for sidepair in itertools.combinations(sides,2):
        contacts=[151,156,*sidepair]
        U=X[contacts];gram=U@M@U.T
        if abs(np.linalg.det(gram))<1e-9:
            failures['singular_contact_frame']+=1;continue
        inv=np.linalg.inv(gram);ker=null_space(U@M)[:,0];kn=ker@M@ker
        assert abs(kn)>1e-9
        for support in support_labels:
            for outer in itertools.product(angles,repeat=2):
                labels=(*support,*outer);trials+=1
                target=-np.sqrt(2*norms[contacts])*[cs[m] for m in labels]
                p=(inv@target)@U;rad=(2-p@M@p)/kn
                if rad < -1e-8:failures['no_real_branch']+=1;continue
                for sign in [-1,1]:
                    real_branches+=1;r=p+sign*np.sqrt(max(rad,0))*ker
                    assert abs(r@M@r-2)<1e-6
                    if r@M@H>=0:failures['center']+=1;continue
                    co=-(X@M@r)/np.sqrt(2*norms)
                    bad=next(((j,relation(c)) for j,c in enumerate(co)
                              if relation(c) in ['obtuse','nonintegral']),None)
                    if bad:failures[str(bad)]+=1;continue
                    if any(np.linalg.norm(r-np.array(c['orthonormal_chart']))<1e-7 for c in candidates):
                        failures['duplicate']+=1;continue
                    ar=r[1:]/r[0] # Xpoint=sqrt(k)*(e0-sum z_i*e_i).
                    assert r[0]>0
                    # Hence Xpoint.r = -sqrt(k)*(r0+z.rspace), so constraint is -rspace.z/r0<=1.
                    ar=-ar
                    cut=(oldpoints@ar)>1+1e-7
                    if not any(cut):failures['no_old_exterior_cut']+=1;continue
                    normal=r@G['F'];check=constraint(G,normal)
                    assert max(abs(ar-check))<2e-6
                    fa=ar@V;fb=1-ar@base
                    pol=polygon(base,V,np.r_[A,[fa]],np.r_[b,fb])
                    candidates.append(dict(id=len(candidates),contacts=contacts,labels=list(labels),
                       branch=sign,orthonormal_chart=r.tolist(),old_chart=normal.tolist(),
                       residual_constraint=ar.tolist(),face_constraint=fa.tolist(),face_rhs=float(fb),
                       old_vertex_cuts=np.where(cut)[0].tolist(),face=pol,
                       contacts_by_kind=dict(Counter(relation(c) for c in co))))
    pairs=[]
    for a,brow in itertools.combinations(candidates,2):
        pair=float(np.array(a['orthonormal_chart'])@M@np.array(brow['orthonormal_chart']))
        kind=relation(-pair/2)
        if kind not in ['obtuse','nonintegral']:
            pairs.append(dict(ids=[a['id'],brow['id']],pair=pair,relation=kind))
    src=[f'scripts/{TAG}.py','scripts/fresh_h20_release144_cap_20261016_residual_scout.py',
         'scripts/fresh_h20_triangle_adjacent_h3_20261016.py','results/fresh_h20_release144_cap_20261016.json',
         'results/fresh_h20_release144_cap_20261016_residual_scout.json']
    return dict(status='PASS_BOUNDED_PENTAGON_CONTACT_SEARCH',max_label=max_label,
       sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in src},
       support=[151,156],side_pairs=list(map(list,itertools.combinations(sides,2))),
       support_label_pairs=list(map(list,support_labels)),trials=trials,real_branches=real_branches,
       failures=dict(failures),old_face=old,candidates=candidates,compatible_pairs=pairs,
       scope='Bounded numerical contact patterns, all157 pair tests. No exact caps, full label completeness, or finite-volume claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--max-label',type=int,default=12)
    ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run(args.max_label);path=ROOT/f'results/{TAG}.json'
    if args.verify:assert close(d,json.loads(path.read_text()))
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['trials'],len(d['candidates']),len(d['compatible_pairs']),flush=True)
