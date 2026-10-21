#!/usr/bin/env python3
"""Two new exact right-coupled caps in the PD7 residual H2 face.

This is a finite constructive checkpoint, not a completion: the file also
replays a strict spacelike escape of all29 walls. Proposal generation uses
floating point only to choose rational free coordinates and the escape.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib
import numpy as np
from scipy.optimize import linprog
import fresh_h9_makarov_joint_20261013_exact_independent as I
import fresh_h9_makarov_joint_completion_20261013_face_bridge as B
import fresh_h9_makarov_joint_completion_20261013_profiles as A
from fresh_h9_makarov_joint_completion_20261013_face_bridge import solve,dot,sumk,enc,SUPPORT
from fresh_h9_makarov_joint_completion_20261013_profiles import setup,ROOT,P,CORE
FIRST=ROOT/f'results/{P}_face_bridge.json'
GATE=ROOT/f'results/{P}_old_bank_gate.json'
OUT=ROOT/f'results/{P}_coupled_exact.json'
Z=I.Z
def encodev(v):return [enc(x)for x in v]
def cov(v,Q):return [sumk(I.mul(v[i],Q[i][j])for i in range(10))for j in range(10)]
def pivot_columns(M):
    A=[r[:]for r in M];piv=[];i=0
    for j in range(len(A[0])):
        choices=[k for k in range(i,len(A))if I.sign(A[k][j])]
        if not choices:continue
        k=choices[0];A[i],A[k]=A[k],A[i];u=I.reciprocal(A[i][j]);A[i]=[I.mul(x,u)for x in A[i]]
        for k in range(len(A)):
            if k!=i:
                u=A[k][j];A[k]=[I.sub(x,I.mul(u,y))for x,y in zip(A[k],A[i])]
        piv.append(j);i+=1
        if i==len(A):break
    return piv
def load():
    src=json.loads(CORE.read_text());first=json.loads(FIRST.read_text());Q=[[I.k(*x)for x in row]for row in src['form']]
    W=[[(I.k(*x),Z)for x in row]for row in src['seed_normals']];D=[Z]*11
    for row in src['caps']:W.append([(I.k(*a),I.k(*b))for a,b in row['normal']]);D.append(I.k(*row['radicand']))
    return src,first,Q,W,D,[I.k(*x)for x in src['strict_center']]
def newpair(v,d,w,e,Q,expected='dotted',oldnorm=None):
    p=I.pair([(Z,x)for x in v],w,Q);sg=I.psign(p,d,e)
    if expected=='m2':assert sg==0
    else:
        assert sg<0
        if oldnorm is None:oldnorm=I.k(2)
        gap=I.psign(I.psub(I.pmul(p,p,d,e),I.scalar(I.mul(I.k(2),oldnorm))),d,e)
        assert gap>0
    return I.encode(p)
def scout_second(first,Q,W,D,c):
    v=[I.k(*x)for x in first['projective_normal']];d=I.k(*first['radicand'])
    M=[cov([a for a,b in W[i]],Q)for i in SUPPORT]+[cov(v,Q)]
    piv=pivot_columns(M);assert len(piv)==8;free=[i for i in range(10)if i not in piv]
    proposal=first['numerical_coupled_caps'][1]['normal']
    for den in[1,2,4,8,16,32,64,128,256,512,1000,10000]:
        u=[Z]*10
        for i in free:u[i]=I.k(F(round(proposal[i]*den),den))
        rhs=[I.neg(sumk(I.mul(row[j],u[j])for j in free))for row in M]
        for i,x in zip(piv,solve([[row[j]for j in piv]for row in M],rhs)):u[i]=x
        q=dot(u,u,Q)
        if I.sign(q)<=0 or I.sign(dot(u,c,Q))>=0:continue
        e=I.mul(I.k(2),I.reciprocal(q))
        try:
            for j,w in enumerate(W):newpair(u,e,w,D[j],Q,'m2'if j in SUPPORT else'dotted',I.k(F(7,11),F(-1,11))if j==6 else I.k(2))
            newpair(u,e,[(Z,x)for x in v],d,Q,'m2')
        except AssertionError:continue
        return {'projective_normal':encodev(u),'projective_norm':enc(q),'radicand':enc(e),'normal':[[enc(Z),enc(x)]for x in u],'free_coordinates':free,'free_coordinate_denominator':den}
    raise AssertionError('no exact second proposal')
def scout_escape(Q,W,D,c):
    f=lambda x:float(x[0])+np.sqrt(5)*float(x[1]);qf=np.array([[f(x)for x in row]for row in Q]);cf=np.array([f(x)for x in c])
    wf=np.array([[f(a)+np.sqrt(f(d))*f(b)for a,b in w]for w,d in zip(W,D)]);A=wf@qf
    rng=np.random.default_rng(2709027);best=(-1e100,None)
    for i in range(120):
        direction=qf@best[1]if i%4==0 and best[1]is not None else rng.normal(size=10)
        sol=linprog(-direction,A_ub=A,b_ub=np.zeros(len(W)),A_eq=[cf@qf],b_eq=[-1],bounds=[(None,None)]*10,method='highs')
        if sol.success:
            q=sol.x@qf@sol.x
            if q>best[0]:best=(q,sol.x)
    assert best[0]>0
    for eps in[.001,.0001,.00001]:
        x=(1-eps)*best[1]+eps*cf/(-cf@qf@cf)
        for den in[10**6,10**8,10**10]:
            y=[I.k(int(round(a*den)))for a in x];point=[(a,Z)for a in y]
            if I.sign(dot(y,y,Q))<=0 or I.sign(dot(y,c,Q))>=0:continue
            if all(I.psign(I.pair(point,w,Q),Z,d)<0 for w,d in zip(W,D)):
                return [int(a[0])for a in y]
    raise AssertionError('strict escape rounding failed')
def verify(data):
    src,first,Q,W,D,c=load();assert I.inertia(Q)==[9,1]
    expected={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[CORE,FIRST,GATE,Path(I.__file__),Path(B.__file__),Path(A.__file__),Path(__file__)]}
    assert data['sources']==expected
    assert data['new_caps'][0]['projective_normal']==first['projective_normal']
    caps=[];pairrows=[]
    for k,row in enumerate(data['new_caps']):
        v=[I.k(*x)for x in row['projective_normal']];q=dot(v,v,Q);d=I.k(*row['radicand'])
        assert I.sign(q)>0 and row['projective_norm']==enc(q)and I.mul(q,d)==I.k(2)
        assert row['normal']==[[enc(Z),enc(x)]for x in v]
        assert I.sign(dot(v,c,Q))<0
        for j,w in enumerate(W):
            tag='m2'if j in SUPPORT or j>=27 else'dotted'
            pairrows.append({'walls':[j,27+k],'type':tag,'pairing':newpair(v,d,w,D[j],Q,tag,I.k(F(7,11),F(-1,11))if j==6 else I.k(2))})
        W.append([(Z,x)for x in v]);D.append(d);caps.append((v,q))
    assert len(pairrows)==55
    target=[I.k(a)for a in json.loads(GATE.read_text())['escape']['point']]
    targetpairs=[dot(v,target,Q)for v,q in caps];assert all(I.sign(x)>0 for x in targetpairs)
    # Rational-over-K projection to the PD7+2A1 corner; both new roots are
    # orthogonal to the support and to one another, so scales cancel.
    base=[[a for a,b in W[i]]for i in SUPPORT];gram=[[dot(v,w,Q)for w in base]for v in base]
    assert I.inertia(gram)==[7,0]
    coeff=solve(gram,[dot(v,c,Q)for v in base]);p=[I.sub(c[i],sumk(I.mul(a,v[i])for a,v in zip(coeff,base)))for i in range(10)]
    for v,q in caps:
        a=I.mul(dot(v,c,Q),I.reciprocal(q));p=[I.sub(x,I.mul(a,y))for x,y in zip(p,v)]
    assert I.sign(dot(p,p,Q))<0 and I.sign(dot(p,c,Q))<0
    signs=[I.psign(I.pair([(a,Z)for a in p],w,Q),Z,d)for w,d in zip(W,D)]
    assert max(signs)==0 and [i for i,x in enumerate(signs)if not x]==SUPPORT+[27,28]
    x=data['escape']['point'];assert len(x)==10 and all(type(a)is int for a in x)
    y=[I.k(a)for a in x];assert I.sign(dot(y,y,Q))>0 and I.sign(dot(y,c,Q))<0
    es=[I.psign(I.pair([(a,Z)for a in y],w,Q),Z,d)for w,d in zip(W,D)];assert es==[-1]*29
    out={'new_pair_counts':{'m2':15,'dotted':40},'new_pairs':pairrows,'target_projective_pairings':[enc(x)for x in targetpairs],
         'finite_vertex':{'point':encodev(p),'norm':enc(dot(p,p,Q)),'center_pairing':enc(dot(p,c,Q)),'incident':SUPPORT+[27,28],'all29_wall_signs':signs,'positive_gram_proof':'PD7 old support orthogonal to both positive norm new normals, which are mutually right'},
         'escape':{'point':x,'norm':enc(dot(y,y,Q)),'center_pairing':enc(dot(y,c,Q)),'all29_wall_signs':es}}
    return out
def run(replay=False):
    src,first,Q,W,D,c=load()
    if replay:data=json.loads(OUT.read_text())
    else:
        cap0={k:first[k]for k in['projective_normal','projective_norm','radicand','normal']};cap1=scout_second(first,Q,W,D,c)
        new=[cap0,cap1]
        for row in new:W.append([(I.k(*a),I.k(*b))for a,b in row['normal']]);D.append(I.k(*row['radicand']))
        data={'status':'EXACT_29_WALL_PD7_COUPLED_CORNER_WITH_STRICT_ESCAPE','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[CORE,FIRST,GATE,Path(I.__file__),Path(B.__file__),Path(A.__file__),Path(__file__)]},'new_caps':new,'escape':{'point':scout_escape(Q,W,D,c)},'scope':'Retains all27 frozen facets, adds two exact norm2 caps with one actual PD7+2A1 finite vertex and cuts the old-bank common escape. All29 retain an exact strict spacelike escape, so this chamber has infinite volume. No full finite vertex or cusp census, or completion claim.'}
    exact=verify(data)
    if replay:
        for k,v in exact.items():assert data[k]==v,k
        print('PASS exact29 coupled replay')
    else:data.update(exact);OUT.write_text(json.dumps(data,indent=2)+'\n');print('PASS written',OUT)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
