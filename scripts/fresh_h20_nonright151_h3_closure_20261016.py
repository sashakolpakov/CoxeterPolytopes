#!/usr/bin/env python3
"""Bounded nonright151 wall search for the actual H3 face S16+151."""
import argparse, hashlib, itertools, json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection
from fresh_h20_pentagon_closure_scout_20261016 import geometry157
from fresh_h20_release144_cap_20261016_residual_scout import close

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_nonright151_h3_closure_20261016'

def geometry158():
    G=geometry157()
    p=json.loads((ROOT/'results/fresh_h20_pentagon_closure_scout_20261016.json').read_text())['candidates'][3309]
    assert p['contacts']==[151,156,95,145] and p['labels']==[2,2,6,3]
    assert json.loads((ROOT/'results/fresh_h20_compact_quad_cap_20261016.json').read_text())['status']=='PASS_EXACT_158_WALL_COMPACT_QUADRILATERAL_CAP'
    r=np.array(p['old_chart']);a=np.array(p['residual_constraint']);x=np.array(p['orthonormal_chart'])
    G.update(N=np.r_[G['N'],[r]],norms=np.r_[G['norms'],2],X=np.r_[G['X'],[x]],
             A=np.r_[G['A'],[a]],b=np.r_[G['b'],1],labels=G['labels']+[157])
    return G

def face_data(G):
    idx=G['labels'].index(151);eq=G['A'][idx]
    base=eq*G['b'][idx]/(eq@eq);V=null_space(eq[None,:])
    A=G['A']@V;b=G['b']-G['A']@base
    live=np.linalg.norm(A,axis=1)>1e-8
    assert min(b[~live])>-1e-7
    return base,V,A[live],b[live],[G['labels'][i] for i in np.where(live)[0]]

def poly(base,V,A,b,ids):
    lens=np.linalg.norm(A,axis=1);aa=A/lens[:,None];bb=b/lens
    lp=linprog([0,0,0,-1],A_ub=np.c_[aa,np.ones(len(A))],b_ub=bb,
               bounds=[(None,None)]*3+[(0,None)],method='highs')
    if not lp.success or lp.x[3]<=1e-10:return dict(degenerate=True)
    yy=HalfspaceIntersection(np.c_[aa,-bb],lp.x[:3]).intersections
    points=[]
    for y in yy:
        z=base+V@y
        if not any(np.linalg.norm(z-np.array(p['z']))<1e-6 for p in points):
            points.append(dict(z=z.tolist(),radius2=float(z@z),
                active=sorted({151}|{ids[i] for i in range(len(ids)) if abs(A[i]@y-b[i])<1e-7})))
    points.sort(key=lambda p:(p['active'],p['z']))
    return dict(vertices=points,finite=sum(p['radius2']<1-1e-6 for p in points),
       near_sphere=sum(abs(p['radius2']-1)<=1e-6 for p in points),
       exterior=sum(p['radius2']>1+1e-6 for p in points),
       maximum_squared_radius=max(p['radius2'] for p in points))

def gate_mask(co):
    valid=(abs(co)<2e-6)|(co>=1-2e-6)
    inside=(co>0)&(co<1-2e-6)
    m=np.zeros_like(co);m[inside]=np.pi/np.arccos(co[inside])
    valid|=inside&(m>=2.9999)&(abs(m-np.rint(m))<2e-4)
    return np.all(valid,axis=1)

def run():
    G=geometry158();X,M=G['X'],G['metric'];base,V,A,b,ids=face_data(G)
    old=poly(base,V,A,b,ids)
    outside=np.array([p['z'] for p in old['vertices'] if p['radius2']>1+1e-6])
    assert len(outside)>=7
    sides=[31,95,132,142,145,150,152,156,157]
    labels=[2,3,4,5,6,8,10]
    patterns=np.array([(first,*rest) for first in [2,3,4,5,6] for rest in itertools.product(labels,repeat=3)])
    triples=patterns[:,1:]
    cs=np.cos(np.pi/patterns);cs[patterns==2]=0
    H=(G['H']@G['g']@G['F'].T)@M
    pool=[];seen=set();stats=Counter();frame_counts=[]
    for side3 in itertools.combinations(sides,3):
        frame=[151,*side3];U=X[frame];gram=U@M@U.T
        if abs(np.linalg.det(gram))<1e-10:
            stats['singular_frames']+=1;continue
        inv=np.linalg.inv(gram);ker=null_space(U@M)[:,0];kn=ker@M@ker
        assert abs(kn)>1e-10
        targets=-np.sqrt(2*G['norms'][frame])*cs
        particular=targets@inv@U
        rad=(2-np.einsum('ij,jk,ik->i',particular,M,particular))/kn
        real=np.where(rad>=-1e-8)[0];stats['patterns']+=len(triples)
        local=0
        for sign in [-1,1]:
            roots=particular[real]+sign*np.sqrt(np.maximum(rad[real],0))[:,None]*ker
            stats['real_branches']+=len(roots)
            mask=(roots@M@H<0)&(roots[:,0]>0)
            co=-(roots@M@X.T)/np.sqrt(2*G['norms'])
            mask&=gate_mask(co)
            for ri in np.where(mask)[0]:
                r=roots[ri];assert abs(r@M@r-2)<1e-5
                a=-r[1:]/r[0]
                cuts=(outside@a)>1+1e-7
                if not any(cuts):continue
                key=tuple(np.round(r,6))
                if key in seen:continue
                seen.add(key);local+=1
                pool.append(dict(id=len(pool),contacts=frame,labels=list(map(int,patterns[real[ri]])),
                   branch=sign,orthonormal_chart=r.tolist(),old_chart=(r@G['F']).tolist(),
                   residual_constraint=a.tolist(),old_exterior_cuts=np.where(cuts)[0].tolist(),
                   minimum_old_exterior_margin=float(min(outside@a-1))))
        frame_counts.append(dict(contacts=frame,unique_accepted=local))
        print('FRAME',frame,'new',local,'pool',len(pool),flush=True)
    # Every all-old-exterior cutter is checked for replacement rays. If none
    # exists, check the128 strongest candidates by old-ray coverage and margin.
    allcut=[p for p in pool if len(p['old_exterior_cuts'])==len(outside)]
    chosen=allcut if allcut else sorted(pool,key=lambda p:(len(p['old_exterior_cuts']),p['minimum_old_exterior_margin']),reverse=True)[:128]
    outcomes=[]
    for p in chosen:
        a=np.array(p['residual_constraint']);aa=a@V;bb=1-a@base
        result=poly(base,V,np.r_[A,[aa]],np.r_[b,bb],ids+[158])
        outcomes.append(dict(candidate=p['id'],face=result))
    sources=[f'scripts/{TAG}.py','scripts/fresh_h20_pentagon_closure_scout_20261016.py',
       'results/fresh_h20_pentagon_closure_scout_20261016.json',
       'results/fresh_h20_compact_quad_cap_20261016.json',
       'scripts/fresh_h20_release144_cap_20261016_residual_scout.py']
    return dict(status='PASS_BOUNDED_ADJACENT_H3_CONTACT_SEARCH',
       sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
       support='S16+151',new_label_to151=[2,3,4,5,6],label_domain=labels,side_domain=sides,stats=dict(stats),frames=frame_counts,
       old_face=old,candidates=pool,all_old_exterior_cutters=len(allcut),tested_full_faces=outcomes,
       scope='Numerical bounded4contact patterns with151 label2through6 and rightS16, all158old angle gates. Complete3Dhalfspace intersection for selected candidates only, tolerance1e-6. No exact newwall, cusp or finite-volume claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert close(d,json.loads(path.read_text()))
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['candidates']),d['all_old_exterior_cutters'],
          [(x['candidate'],x['face'].get('exterior')) for x in d['tested_full_faces'] if x['face'].get('exterior')==0],flush=True)
