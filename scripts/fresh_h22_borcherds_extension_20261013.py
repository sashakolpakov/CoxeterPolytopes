#!/usr/bin/env python3
"""Numerical H22 extension scout, with exact reconstruction of old vertices.

Keep Borcherds's 210 H21 normals at height zero. Add a symmetric pair over
the center and caps over sampled original cusps/finite vertices. This is a
proposal search, not a complete mirror or vertex enumeration.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd,lcm
from pathlib import Path
import numpy as np
import scipy.linalg as la
from scipy.optimize import linprog
import sympy as sp
import fresh_h20_global_fraction_verify_20261013 as base

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/enumerations/fresh_h22_borcherds_extension_20261013_scout.json'

def prepare():
    Q,W,G,h=base.source();Q=np.array(Q,dtype=np.int64);W=np.array(W,dtype=np.int64);h=np.array(h,dtype=np.int64)
    Z=la.null_space((h@Q)[None,:]);d,V=la.eigh(Z.T@Q@Z);assert min(d)>0
    X=Z@V@np.diag(1/np.sqrt(d));assert np.max(abs(X.T@Q@X-np.eye(21)))<1e-7
    A=W@Q@X;b=-W@Q@h/np.sqrt(14);Qn=la.block_diag(Q,1)
    R=np.pad(W,((0,0),(0,1))).astype(float)
    return Q,W,np.array(G),h,X,A,b,Qn,R

def gate(R,Q,tol=2e-6,maxlabel=120):
    H=R@Q@R.T;n=np.diag(H);assert min(n)>0
    C=-H/np.sqrt(np.outer(n,n));allowed=np.cos(np.pi/np.arange(2,maxlabel+1))
    D=np.min(abs(C[:,:,None]-allowed),axis=2)
    ok=(C>=1-tol)|(D<tol);np.fill_diagonal(ok,True)
    return ok,C

def reconstruct(v,active,Q,W,h):
    r=[F(float(x)).limit_denominator(2000) for x in v]
    if any(sum(F(int(a))*b for a,b in zip(row,r))>0 for row in W@Q) or sum(F(int(a))*b for a,b in zip(h@Q,r))!=-1:
        M=sp.Matrix((W[list(active)]@Q).tolist());ids=M.T.rref()[1];M=M[list(ids),:]
        if M.rows!=21:return None
        r=list(M.col_join(sp.Matrix([list(map(int,h@Q))])).inv()*sp.Matrix([0]*21+[-1]))
        r=[F(int(x.p),int(x.q)) for x in r]
    ps=[sum(F(int(a))*b for a,b in zip(row,r)) for row in W@Q]
    if max(ps)>0:return None
    q=sum(r[i]*int(Q[i,j])*r[j] for i in range(22) for j in range(22));assert q<=0
    den=lcm(*(x.denominator for x in r));v=[int(x*den) for x in r];g=gcd(*v);v=[x//g for x in v]
    return {'vector':v,'norm':str(q*den*den/g**2),'incident_walls':[i for i,x in enumerate(ps) if not x]}

def run(count=1000,seed=22013):
    Q,W,G,h,X,A,b,Qn,R=prepare();rng=np.random.default_rng(seed);vertices={};counts=Counter()
    for k in range(count):
        c=rng.normal(size=21);lp=linprog(c,A_ub=A,b_ub=b,bounds=[(None,None)]*21,method='highs')
        if not lp.success:counts['lp_failure']+=1;continue
        active=tuple(np.flatnonzero(abs(A@lp.x-b)<1e-6));v=(h/np.sqrt(14)+X@lp.x)/np.sqrt(14)
        row=reconstruct(v,active,Q,W,h)
        if row is None:counts['reconstruction_failure']+=1;continue
        vertices[tuple(row['vector'])]=row;counts['samples']+=1
        if k%100==0:print('old vertices',k,len(vertices),flush=True)
    rows=list(vertices.values());caps=[]
    # Cusp caps: integer old direction, height chosen by a short-wall label.
    # Finite projections: fix norm via one or two independently prescribed
    # squared cosines. Retain only all210-pair numerical gate survivors.
    labels=[3,4,5,6,7,8,10,12,16,24,40,80,120]
    for vi,row in enumerate(rows):
        v=np.array(row['vector'],dtype=float);q=float(F(row['norm']));p=v@Q@W.T;oldn=G.diagonal()
        possible=set()
        for j in np.flatnonzero(p<-1e-8):
            for m in labels:
                newq=p[j]**2/(oldn[j]*np.cos(np.pi/m)**2)
                if newq-q>0:possible.add(round(float(np.sqrt(newq-q)),10))
        # One all-dotted proposal, keeping the original projected direction.
        positive=(-p/np.sqrt(oldn))[p<-1e-8]
        if len(positive):
            newq=min(positive)**2
            if newq-q>0:possible.add(round(float(np.sqrt(newq-q)),10))
        for t in sorted(possible):
            newq=q+t*t
            if newq<=1e-8:continue
            cs=-p/np.sqrt(newq*oldn);allowed=np.cos(np.pi/np.arange(2,121))
            if np.any((cs<1-2e-6)&(np.min(abs(cs[:,None]-allowed),axis=1)>2e-6)):continue
            for sign in(1,-1):caps.append({'old_vertex':vi,'height':sign*t,'normal':list(v)+[sign*t],'norm':newq})
    print('old directions',len(rows),'caps',len(caps),'types',Counter(F(x['norm'])==0 for x in rows),flush=True)
    bank=np.array([r['normal'] for r in caps]);bankn=np.array([r['norm'] for r in caps]);B=bank[:,:22]@Q@X if len(bank) else np.empty((0,21))
    bankb=-bank[:,:22]@Q@h/np.sqrt(14) if len(bank) else np.empty(0)
    runs=[]
    for m in [3,4,5,6,8,12,24]:
        norm=1/(2*np.cos(np.pi/m)**2);T=np.sqrt(14+norm);center=np.array([list(h)+[T],list(h)+[-T]])
        centerpairs=bank@Qn@center.T if len(bank) else np.empty((0,2));cs=-centerpairs/np.sqrt(bankn[:,None]*norm) if len(bank) else np.empty((0,2))
        allowed=np.cos(np.pi/np.arange(2,121));ok=np.all((cs>=1-2e-6)|(np.min(abs(cs[:,:,None]-allowed),axis=2)<2e-6),axis=1)
        ids=np.flatnonzero(ok);allA=np.vstack([np.column_stack([A,np.zeros(210)]),np.array([[0]*21+[T],[0]*21+[-T]]),np.column_stack([B[ids],bank[ids,-1]])]);allb=np.concatenate([b,[np.sqrt(14)]*2,bankb[ids]])
        escapes=[]
        for k in range(100):
            lp=linprog(rng.normal(size=22),A_ub=allA,b_ub=allb,bounds=[(None,None)]*22,method='highs')
            if lp.success and lp.x@lp.x>1+1e-7:
                v=np.r_[h/np.sqrt(14)+X@lp.x[:21],lp.x[-1]]
                escapes.append({'norm':float(lp.x@lp.x-1),'point':v.tolist(),'active_old':np.flatnonzero(abs(A@lp.x[:21]-b)<1e-6).tolist()});break
        row={'center_short_label':m,'center_norm':norm,'center_height':T,'compatible_bank_count':len(ids),'escape':escapes[0] if escapes else None,'note':'Whole-bank relaxation; bank members may conflict with one another.'}
        runs.append(row);print('m',m,'bank',len(ids),'escape',None if not escapes else escapes[0]['norm'],flush=True)
    result={'status':'BOUNDED_NUMERICAL_H22_EXTENSION_SCOUT','source':str(base.SOURCE.relative_to(ROOT)),'samples':dict(counts),'old_vertices':rows,'caps':caps,'runs':runs,'scope':'Bounded LP sampling and finite-angle profile generation, max numerical label120. Exact old vertex reconstruction only; cap gates and escaped H22 points are numerical. Whole banks need not be compatible.'}
    OUT.write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--samples',type=int,default=1000);a=p.parse_args();run(a.samples)
