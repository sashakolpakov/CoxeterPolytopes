#!/usr/bin/env python3
"""Numerical adjacent-H3 geometry and bounded finite-partner target.
All old data frozen; no exact cone/angle proof is asserted by this scout.
"""
import argparse,hashlib,itertools,json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection,ConvexHull
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_triangle_adjacent_h3_20261016'
def load(tag):return json.loads((ROOT/'results'/f'{tag}.json').read_text())
def mid(x):return sum(float(F(t))for t in x)/2

def geometry():
    D=load('fresh_h20_residual_h4_corner_20261016');M=load('fresh_h20_mixed_isolated_cap_20261016');P=load('fresh_h20_h4_parallel_partner_20261016');T=load('fresh_h20_second_cusp_completion_20261016')
    R=np.array([[float(F(x))for x in row]for row in D['definition']['four_recipes']]);q=np.array([float(F(x))for x in D['definition']['four_norms']]);C=np.array([float(F(x))for x in D['new_finite_corner']['chart']]);k=-float(F(D['new_finite_corner']['norm']));g=np.diag([-1,1,1,.5,.5])
    cp=np.array([mid(x)for x in D['new_finite_corner']['all145_pair_bounds']]);pp=np.zeros((145,4))
    for row in M['pairs']:pp[row['wall'],0]=mid(row['pair_bounds'])
    for row in D['pairs']:
     a,b=row['walls'];pp[a,b-141]=mid(row['pair_bounds'])
    for i in range(4):pp[141+i,i]=q[i]
    N=np.array([sum((pp[j,i]/q[i])*R[i]for i in range(4))-cp[j]*C/k for j in range(145)])
    assert np.max(np.abs(N[141:]-R))<1e-8
    nu=C/np.sqrt(k)-R[3]/np.sqrt(q[3]);z=mid(P['definition']['z_bounds']);r=-np.sqrt(2/q[0])*R[0]+np.sqrt(2)*z*nu
    N=np.r_[N,[r],[3*nu],[3*nu],[3*nu],[-R[1]+5*nu],[-R[2]+6*nu]]
    norms=np.r_[np.full(140,2.),1.46,q,[2,2,2,2,q[1],q[2]]];assert len(N)==len(norms)==151
    E5=np.array([1,-.01,-.98,-.3,-.04]);E4=E5.copy();E4[-1]=-.05
    S16=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,61,87]

    U=N[[95,142,145]];JU=np.linalg.inv(U@g@U.T)
    proj=lambda w:(w@g@U.T)@JU@U
    v=N[150]-proj(N[150]);vn=v@g@v
    raw=np.array([0,0,0,0,1.]);kv=raw-proj(raw);kv-=((kv@g@v)/vn)*v;kn=kv@g@kv
    target=-np.sqrt(2*q[2])*np.cos(np.pi/8)
    rr=(target/vn)*v-np.sqrt((2-target*target/vn)/kn)*kv
    NN=np.r_[N,[rr]]
    hp=np.array([mid(x)for x in D['common_center']['all145_pair_bounds']]);hc=mid(D['new_finite_corner']['center_pair_bounds'])
    H=sum(hp[141+i]/q[i]*R[i]for i in range(4))-hc*C/k
    V=NN[[142,145,150,151]];Jc=np.linalg.inv(V@g@V.T)
    C8=H-(H@g@V.T)@Jc@V;k8=-C8@g@C8
    rr=-N[142]/np.sqrt(q[1])-np.sqrt(3)*N[145]/2+C8/np.sqrt(2*k8)
    N=np.r_[NN,[rr]];norms=np.r_[norms,2.,2.]
    B=-(R@g@N.T).T*np.sqrt(k/q);bb=-N@g@C
    ids=[i for i in range(153)if i not in S16 and i not in range(141,145)]
    A=np.r_[B[ids]/bb[ids,None],-np.eye(4)];b=np.r_[np.ones(len(ids)),np.zeros(4)]
    labels=ids+list(range(141,145))
    return dict(N=N,norms=norms,g=g,R=R,q=q,C=C,k=k,H=H,C8=C8,S16=S16,A=A,b=b,labels=labels)

def face_census(G,face):
    A,b=G['A'],G['b'];idx=G['labels'].index(face);aa=A[idx];bb=b[idx]
    base=aa*bb/(aa@aa);V=null_space(aa[None,:]);A3=A@V;b3=b-A@base
    nn=np.linalg.norm(A3,axis=1);active=np.where(nn>1e-8)[0]
    assert min(b3[nn<=1e-8],default=0)>-1e-7
    aa3=A3[active]/nn[active,None];bb3=b3[active]/nn[active]
    lp=linprog([0,0,0,-1],A_ub=np.c_[aa3,np.ones(len(aa3))],b_ub=bb3,bounds=[(None,None)]*3+[(0,None)],method='highs')
    assert lp.success and lp.x[3]>1e-9
    vv=[]
    for y in HalfspaceIntersection(np.c_[aa3,-bb3],lp.x[:3]).intersections:
        z=base+V@y
        if any(np.linalg.norm(z-w['z'])<1e-7 for w in vv):continue
        act=[G['labels'][i]for i in range(len(b))if abs(A[i]@z-b[i])<1e-7]
        vv.append(dict(z=z.tolist(),squared_radius=float(z@z),active=sorted(set(act))))
    vv.sort(key=lambda v:(v['active'],v['z']))
    counts=dict(total=len(vv),finite=sum(v['squared_radius']<1-1e-6 for v in vv),near_sphere=sum(abs(v['squared_radius']-1)<=1e-6 for v in vv),exterior=sum(v['squared_radius']>1+1e-6 for v in vv))
    # Reproducible volume estimate: uniformly distributed Sobol points in the affine-coordinate bounding box.
    from scipy.stats import qmc
    yy=np.array([(np.array(v['z'])-base)@V for v in vv]);lo=yy.min(0);hi=yy.max(0)
    sample=lo+qmc.Sobol(3,scramble=False).random_base2(18)*(hi-lo)
    inpoly=np.max(sample@aa3.T-bb3,axis=1)<=1e-10
    zz=base+sample@V.T;outside=np.sum(zz*zz,axis=1)>1
    box=float(np.prod(hi-lo));volume=ConvexHull(yy).volume
    return dict(support=G['S16']+[face],vertices=vv,counts=counts,maximum_squared_radius=max(v['squared_radius']for v in vv),relative_inradius=float(lp.x[3]),affine_euclidean_volume=float(volume),estimated_exterior_euclidean_volume=box*int(np.sum(inpoly&outside))/len(sample),volume_method='Deterministic 2^18 Sobol box samples in an orthonormal face chart; numerical only')

def projected_escape_test(G,face,allvertices):
    N,g=G['N'],G['g'];E=np.array([1,-.01,-.98,-.28,-.05]);n=N[face];q=G['norms'][face]
    Ep=E-(E@g@n/q)*n
    finite=[v for v in allvertices if v['squared_radius']<1-1e-6]
    points=np.array([G['C']-(np.array(v['z'])*np.sqrt(G['k']/G['q']))@G['R']for v in finite])
    # Permit an unrestricted multiple of the released normal; right-to-face walls see the same target.
    mat=np.c_[points.T,n,-n]
    lp=linprog(np.r_[np.ones(len(points)),0.,0.],A_eq=mat,b_eq=Ep,bounds=[(0,None)]*mat.shape[1],method='highs')
    ans=dict(projection_norm=float(Ep@g@Ep),projection_chart=Ep.tolist(),known_finite_generators=len(points),cone_feasible=bool(lp.success))
    if lp.success:
        ans.update(max_abs_residual=float(np.max(np.abs(mat@lp.x-Ep))),normal_coefficient=float(lp.x[-2]-lp.x[-1]),positive_terms=[dict(active=finite[i]['active'],coefficient=float(a))for i,a in enumerate(lp.x[:-2])if a>1e-9])
    else:ans['lp_status']=lp.message
    return ans

def run():
    G=geometry();old=load('fresh_h20_pi8_triangle_face_20261016_residual_scout')['new153']['vertices'];faces={}
    for f in (150,151):
        ans=face_census(G,f);ans['E6_projection_cone_test']=projected_escape_test(G,f,old);faces[str(f)]=ans
    src=[f'scripts/{TAG}.py']+[f'results/{x}.json'for x in ['fresh_h20_residual_h4_corner_20261016','fresh_h20_mixed_isolated_cap_20261016','fresh_h20_h4_parallel_partner_20261016','fresh_h20_second_cusp_completion_20261016','fresh_h20_second_cusp_pi8_cap_20261016','fresh_h20_pi8_triangle_face_20261016','fresh_h20_pi8_triangle_face_20261016_residual_scout']]
    return dict(status='PASS_BOUNDED_NUMERICAL_ADJACENT_H3_CENSUS',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},faces=faces,scope='Numerical relative-face vertex and escape-cone diagnostic. No complete exact cone or new wall certified; finite vertex cone tests alone use the saved numerical residual representatives.')
def close(a,b):
    if isinstance(a,float):return isinstance(b,(int,float))and abs(a-b)<=1e-7*(1+abs(a))
    if isinstance(a,dict):return a.keys()==b.keys()and all(close(v,b[k])for k,v in a.items())
    if isinstance(a,list):return len(a)==len(b)and all(close(x,y)for x,y in zip(a,b))
    return a==b
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],{f:(v['counts'],v['maximum_squared_radius'],v['E6_projection_cone_test'])for f,v in d['faces'].items()},flush=True)
