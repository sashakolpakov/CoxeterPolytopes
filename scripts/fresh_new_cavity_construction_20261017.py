#!/usr/bin/env python3
"""Four-center changed-pentagon patch: mandatory shared endpoint discovery.

The interior cycle is A-R-B-C-A. R covers old edges8-44 and44-63;
A covers7-8, B covers63-30, C covers30-7. Either A-B or R-C
triangulates its central quadrilateral. Numerical geometry only.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
import fresh_three_cap_pentagon_changed_seed_20261017 as S

A,H,B,G,M,ROOT=S.A,S.H,S.B,S.G,S.M,S.ROOT
TAG='fresh_new_cavity_construction_20261017'

def finite_at8(r,m):
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for i,l in enumerate(G.profile_labels(r['recipe']['profile'])):g[i][4]=g[4][i]=l
    g[2][4]=g[4][2]=r['recipe']['labels'][0]
    g[3][4]=g[4][3]=m
    return A.components(tuple(v for row in g for v in row))

def norm_roots(frame,rhs,square=2.):
    f=np.array(frame);gram=f@M@f.T
    if abs(np.linalg.det(gram))<1e-10:return [],dict(status='unresolved_singular_frame')
    inv=np.linalg.inv(gram);part=np.array(rhs)@inv@f
    k=np.linalg.svd(f@M)[2][-1];kq=float(k@M@k);remainder=float(square-part@M@part)
    if abs(kq)<1e-10:return [],dict(status='unresolved_null_kernel')
    ratio=remainder/kq
    if abs(ratio)<1e-10:return [],dict(status='unresolved_tangent_norm',ratio=ratio)
    if ratio<0:return [],dict(status='no_real_norm_root',ratio=ratio)
    return [(s,part+s*np.sqrt(ratio)*k)for s in[-1,1]],dict(status='two_real_norm_branches',ratio=ratio,gram=gram.tolist(),kernel_square=kq)

def test_root(y,sign,labels,ns,h,r,oldpair,m):
    values={i:float(y@M@ns[i])for i in labels if i not in S.DELETE}
    relations={i:B.relation(v)for i,v in values.items()}
    for i in oldpair:relations[i]=2
    rv=float(y@M@np.array(r['projection']))
    bad=[dict(wall=i,product=values[i])for i in values if relations[i]is None or relations[i]=='parallel_or_unresolved']
    return dict(branch=sign,projection=y.tolist(),center_product=float(y@M@h),R_product=rv,
        assigned_R_label=m,blockers=bad,finite_old_contacts=[[i,x]for i,x in relations.items()if isinstance(x,int)])

def run():
    prior=json.loads((ROOT/'results'/(S.TAG+'.json')).read_text())
    z,labels,ns,h,profiles,primary=H.setup();actual=[i for i in z['actual_sides']if i not in S.DELETE]
    assert prior['deleted_walls']==S.DELETE and len(prior['accepted_low_caps'])==94
    assert prior['deletion_boundary']['counts']=={'finite':15,'exterior':7}and prior['deletion_boundary']['open_ends']==11
    rows=[];accepted=[];counts=Counter();unresolved=[]
    for r in prior['accepted_low_caps']:
        assert r['frame']==[8,44,63]and not r['blockers']and r['center_product']<0
        R=np.array(r['projection']);frame=[ns[7],ns[8],R]
        low=[m for m in range(2,6)if finite_at8(r,m)is not None]
        high=finite_at8(r,6)is not None;trials=[];tails=[]
        for m in low:
            roots,norm=norm_roots(frame,[0.,0.,-B.theta(m)])
            counts['low_norm_equations']+=1;counts['real_low_branches']+=len(roots)
            if norm['status'].startswith('unresolved'):unresolved.append(dict(R=r['index'],label=m,**norm))
            rec=dict(label=m,finite_corner=finite_at8(r,m),norm=norm,branches=[])
            for sign,y in roots:
                t=test_root(y,sign,labels,ns,h,r,[7,8],m);rec['branches'].append(t)
                if not t['blockers']and t['center_product']<-1e-8:
                    nd=dict(ns);nd.update({67:R,68:y});t['boundary']=S.C.detailed_cone(nd,actual+[67,68])
                    accepted.append(dict(R_index=r['index'],R_recipe=r['recipe'],A=t));counts['coherent_low_A']+=1
            trials.append(rec)
        if high:
            assert r['recipe']['profile']==0 and r['recipe']['labels'][0]==2
            counts['complete_high_label_shape_families']+=1
            for c in np.linspace(np.sqrt(3.),2.,65):
                roots,norm=norm_roots(frame,[0.,0.,-float(c)])
                if norm['status'].startswith('unresolved'):unresolved.append(dict(R=r['index'],cosine=float(c),**norm))
                branches=[]
                for sign,y in roots:
                    t=test_root(y,sign,labels,ns,h,r,[7,8],'continuous_high_tail')
                    branches.append(t)
                    if not t['blockers']and t['center_product']<-1e-8:counts['coherent_high_samples']+=1
                tails.append(dict(cosine=float(c),norm=norm,branches=branches));counts['continuous_high_samples']+=1
        rows.append(dict(R_index=r['index'],R_recipe=r['recipe'],R_projection=r['projection'],
            complete_low_labels=low,unbounded_isolated_R_label=high,low_trials=trials,high_continuous_samples=tails))
        if len(rows)%20==0:print('FOUR_CENTER_ENDPOINT',len(rows),'of94','accepted',len(accepted),'counts',dict(counts),flush=True)
    paths=[Path(__file__),Path(S.__file__),Path(A.__file__),Path(H.__file__),ROOT/'results'/(S.TAG+'.json'),H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_NUMERICAL_FOUR_CENTER_CHANGED_PENTAGON_ENDPOINT_DOMAIN',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=S.DELETE,boundary_cycle=S.CYCLE,deletion_cost=prior['exact_old_deletion_cost'],deletion_boundary=prior['deletion_boundary'],
        patch=dict(old_cycle=S.CYCLE,cap_boundary_arcs={'A':[7,8],'R':[8,44,63],'B':[63,30],'C':[30,7]},
            interior_cycle=['A','R','B','C'],possible_diagonals=[['A','B'],['R','C']],
            mandatory_new_finite_corner=[8,'A','R'],A_profile=0,A_old_right_contacts=[7,8]),
        R_bank_size=94,endpoint_trials=rows,accepted_low_pairs=accepted,unresolved_norm_events=unresolved,counts=dict(counts),
        equation='For each literal R, A lies in the three-contact frame7,8,R with products0,0,-2cos(pi/m). Its norm2 equation has at most two real branches. The required finite cornerS6+8+A+R gives the complete labels2..5 or an explicitly isolated I2(m),m>=6. Both branches are tested against every retained native old wall and the old common center.',
        scope='Fresh four-cap topology, with its necessary A-R endpoint tested over all94 coherent low R seeds. Numerical roots/native gates and65-point continuous high-tail scouts only. No all-parameter exclusion, exact new wall, ambient finite census, or finite-volume claim. Further cap equations are attempted only if the mandatory A-R finite corner survives.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()));dest=ROOT/'results'/(TAG+'.json')
    if a.verify:assert out==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],'accepted',len(out['accepted_low_pairs']),flush=True)
