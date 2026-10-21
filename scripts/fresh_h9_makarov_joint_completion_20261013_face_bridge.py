#!/usr/bin/env python3
"""Exact first PD7 cap and numerical coupled face-bridge proposals.

The first normal's seven right relations, twenty strict dotted relations,
strict center, and common-escape cut are exact. Further caps are explicitly
numerical proposals and are never claimed as an exact chamber.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from scipy.linalg import null_space
import fresh_h9_makarov_joint_20261013_exact_independent as I
from fresh_h9_makarov_joint_completion_20261013_profiles import setup,ROOT,P,CORE
from fresh_h9_makarov_joint_completion_20261013_continuous import intervals,intersect
OUT=ROOT/f'results/{P}_face_bridge.json'
SUPPORT=[0,2,3,5,7,8,9]
Z=I.Z
def solve(A,b):
    A=[row[:]+[x]for row,x in zip(A,b)]
    for j in range(len(A)):
        i=next(i for i in range(j,len(A))if I.sign(A[i][j]));A[i],A[j]=A[j],A[i]
        q=I.reciprocal(A[j][j]);A[j]=[I.mul(x,q)for x in A[j]]
        for i in range(len(A)):
            if i!=j:
                q=A[i][j];A[i]=[I.sub(x,I.mul(q,y))for x,y in zip(A[i],A[j])]
    return [row[-1]for row in A]
def dot(v,w,Q):return sumk(I.mul(I.mul(x,Q[i][j]),y)for i,x in enumerate(v)for j,y in enumerate(w))
def sumk(vals):
    x=Z
    for y in vals:x=I.add(x,y)
    return x
def enc(x):return list(map(str,x))
def build():
    src,Qf,Wf,cf,Ef=setup();Q=[[I.k(*a)for a in row]for row in src['form']]
    normals=[[(I.k(*a),Z)for a in row]for row in src['seed_normals']];ds=[Z]*11
    for row in src['caps']:
        normals.append([(I.k(*a),I.k(*b))for a,b in row['normal']]);ds.append(I.k(*row['radicand']))
    c=[I.k(*a)for a in src['strict_center']]
    gate=ROOT/f'results/{P}_old_bank_gate.json';E=[I.k(a)for a in json.loads(gate.read_text())['escape']['point']]
    # Coordinates 0,2,3,5,6,7,8 are fixed by the seven exact right equations.
    fixed=[0,2,3,5,6,7,8];free=[1,4,9]
    # This rational free-coordinate proposal is checked exactly below.
    numerical=[F('51.86498164'),F('9.93349942'),F('-32.18449178')]
    accepted=None
    for den in [1,2,4,8,16,32,64,128,256,512,1000,10000]:
        v=[Z]*10
        for i,x in zip(free,numerical):v[i]=I.k(F(round(x*den),den))
        rhs=[I.neg(sumk(I.mul(Q[i][j],v[j])for j in free))for i in fixed]
        for i,x in zip(fixed,solve([[Q[i][j]for j in fixed]for i in fixed],rhs)):v[i]=x
        q=dot(v,v,Q)
        if I.sign(q)<=0 or I.sign(dot(v,c,Q))>=0 or I.sign(dot(v,E,Q))<=0:continue
        d=I.mul(I.k(2),I.reciprocal(q));r=[(Z,x)for x in v];pairs=[]
        for j,w in enumerate(normals):
            p=I.pair(r,w,Q);sg=I.psign(p,d,ds[j])
            if j in SUPPORT:
                if sg:break
                tag='m2'
            else:
                oldq=I.k(F(7,11),F(-1,11))if j==6 else I.k(2)
                sq=I.pmul(p,p,d,ds[j]);gap=I.psign(I.psub(sq,I.scalar(I.mul(I.k(2),oldq))),d,ds[j])
                if sg>=0 or gap<=0:break
                tag='dotted'
            pairs.append({'old_wall':j,'type':tag,'pairing':I.encode(p)})
        else:accepted=(v,q,d,r,pairs,den);break
    assert accepted is not None
    v,q,d,r,pairs,den=accepted
    f=lambda x:float(x[0])+np.sqrt(5)*float(x[1]);rf=np.array([f(x)for x in v])*np.sqrt(f(d))
    assert max(Wf@Qf@rf)<1e-7 and rf@Qf@Ef>0
    # On the intrinsic H2 face, prescribe a finite angle to the first cap.
    T=null_space(Wf[SUPPORT]@Qf);p=cf-(Wf[SUPPORT]@Qf@cf)@np.linalg.inv(Wf[SUPPORT]@Qf@Wf[SUPPORT].T)@Wf[SUPPORT]
    p-= (p@Qf@rf)/2*rf;p/=np.sqrt(-p@Qf@p)
    b=T[:,0]+(T[:,0]@Qf@p)*p-(T[:,0]@Qf@rf)/2*rf
    if b@Qf@b<1e-8:b=T[:,1]+(T[:,1]@Qf@p)*p-(T[:,1]@Qf@rf)/2*rf
    b/=np.sqrt(b@Qf@b);assert abs(p@Qf@b)<1e-6
    outside=sorted(set(range(27))-set(SUPPORT));regions=[];caps=[]
    for m in[2,3,4,5,6,8,10]:
        k=-2*np.cos(np.pi/m)if m>2 else 0;u=k/2*rf;sq=np.sqrt(2-k*k/2)
        A=p@Qf@Wf.T;B=b@Qf@Wf.T;D=u@Qf@Wf.T
        for branch in[-1,1]:
            allowed=[(0.,np.inf)]
            for j in outside:
                allowed=intersect(allowed,intervals(sq*(A[j]+branch*B[j])/2,D[j]+2,sq*(-A[j]+branch*B[j])/2))
                if not allowed:break
            for lo,hi in allowed:
                idx=len(regions);regions.append({'label':m,'branch':branch,'interval':[lo,None if not np.isfinite(hi)else hi]})
                xs=[lo+max(1e-5,abs(lo)*1e-5)]
                xs +=[(lo+hi)/2,hi-max(1e-5,abs(hi)*1e-5)]if np.isfinite(hi)else[max(lo*1.2,lo+1),max(lo*2,lo+2)]
                for x in xs:
                    if not lo<x<hi:continue
                    rr=u+sq*((x-1/x)/2*p+branch*(x+1/x)/2*b)
                    if rr@Qf@cf>=0 or max((Wf@Qf@rr)[outside])>=-2+1e-7:continue
                    caps.append({'id':len(caps),'region':idx,'boost_exp':float(x),'normal':rr.tolist(),'pair_with_first':float(rr@Qf@rf),'target_pairing':float(rr@Qf@Ef)})
    out={'status':'EXACT_PD7_CUTTER_WITH_NUMERICAL_COUPLED_SCOUT','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[CORE,gate,Path(I.__file__),Path(__file__)]},'support':SUPPORT,'free_coordinate_denominator':den,'projective_normal':[enc(x)for x in v],'projective_norm':enc(q),'radicand':enc(d),'normal':[[enc(a),enc(b)]for a,b in r],'all27_pairs':pairs,'center_pairing':enc(dot(v,c,Q)),'target_projective_pairing':enc(dot(v,E,Q)),'numerical_coupled_regions':regions,'numerical_coupled_caps':caps,'scope':'Exactly one norm2 cap with seven right and twenty strict dotted old relations cuts the exact common-bank escape. Coupled cap regions/proposals use numerical predicates only. No global volume, compactness, or complete continuous-family claim.'}
    print('exact first cap denominator',den,'coupled regions',len(regions),'caps',len(caps),flush=True)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=build()
    if a.verify:assert json.loads(OUT.read_text())==out;print('PASS saved exact first-cap replay')
    else:OUT.write_text(json.dumps(out,indent=2)+'\n');print(OUT)
