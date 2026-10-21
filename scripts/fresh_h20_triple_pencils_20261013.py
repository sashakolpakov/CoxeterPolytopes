#!/usr/bin/env python3
"""Exact actual-ridge constraints on central three-wall pencils.

Each unrooted triple has three sign patterns. Both orders of the positive
terms on 0<t<=1 cover every positive parameter. Constraints are obtained from
actual ridges with polynomial timelike points, allowing redundant zero walls
only with exact positive-cone coefficients. No raw pair is called a ridge.
"""
from __future__ import annotations
import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import numpy as np
import sympy as sp
import fresh_h20_triple_central_20261013 as source
import fresh_h20_global_fraction_verify_20261013 as arithmetic

ROOT=source.ROOT
OUTPUT=ROOT/'data/enumerations/fresh_h20_triple_pencils_20261013.json'
t=sp.Symbol('t')

def coeff(p):return [str(x) for x in p.all_coeffs()]
def domain(p):
    a=arithmetic.readp(coeff(p));return arithmetic.transform_unit(a)
def positive(p,strict=True):
    v=domain(p)
    return all(x>=0 for x in v) and (not strict or (any(v) and p.eval(1)>0))
def pp(x):return sp.Poly(x,t,domain=sp.QQ)

def setup(Q,W,G,h,i,j,k):
    cp=[int(x@Q@h) for x in W]
    U=np.array([int(-8//cp[l])*W[l] for l in range(210)],dtype=object)
    e0=U[i]-U[k];e1=U[j]-U[k]
    qe=pp(int(e0@Q@e0)+2*t*int(e0@Q@e1)+t*t*int(e1@Q@e1))
    v=[pp(int(w@Q@e0)+t*int(w@Q@e1)) for w in W]
    def gr(a,b):return qe*int(G[a,b])-v[a]*v[b]
    return cp,e0,e1,qe,v,gr

def candidates(Q,W,G,h,sd,value=1.0):
    cp,e0,e1,qe,v,gr=sd
    Qf=np.array(Q,dtype=float);Wf=np.array(W,dtype=float);hf=np.array(h,dtype=float)
    e=np.array(e0+value*e1,dtype=float);q=e@Qf@e;R=Wf-np.outer(Wf@Qf@e,e)/q
    H=R@Qf@R.T;c=R@Qf@hf
    candidates=[]
    for a,b in combinations(range(210),2):
        aa,bb,ab=H[a,a],H[b,b],H[a,b];D=aa*bb-ab*ab
        if aa<1e-8 or bb<1e-8 or D<1e-8:continue
        ratio=ab*ab/(aa*bb)
        if ab<=0 and (ratio>=.7499 or min(abs(ratio-z) for z in (0,.25,.5,(3+5**.5)/8,(3-5**.5)/8))<1e-8):continue
        A=bb*c[a]-ab*c[b];B=aa*c[b]-ab*c[a]
        values=D*c-A*H[:,a]-B*H[:,b];scale=max(1,abs(values).max())
        if values.max()>1e-8*scale:continue
        candidates.append((a,b))
    return candidates

def constraint(Q,W,G,h,sd,a,b):
    cp,e0,e1,qe,v,gr=sd
    aa,bb,ab=gr(a,a),gr(b,b),gr(a,b);D=aa*bb-ab*ab
    if not all(positive(p) for p in (qe,aa,bb,D)):return None
    obtuse=positive(ab)
    if not obtuse and not positive(-ab,strict=False):return None
    N=ab*ab;den=aa*bb
    if not obtuse and not positive(3*den-4*N):return None
    A=bb*cp[a]-ab*cp[b];B=aa*cp[b]-ab*cp[a]
    vals=[D*cp[l]-A*gr(l,a)-B*gr(l,b) for l in range(210)]
    zeros=[l for l,p in enumerate(vals) if p.is_zero]
    if a not in zeros or b not in zeros:return None
    for l,p in enumerate(vals):
        if l not in zeros and not positive(-p):return None
    cone=[]
    for l in zeros:
        if l in (a,b):continue
        alpha=bb*gr(l,a)-ab*gr(l,b);beta=aa*gr(l,b)-ab*gr(l,a)
        if not positive(alpha,strict=False) or not positive(beta,strict=False):return None
        residual=D*gr(l,l)-alpha*gr(l,a)-beta*gr(l,b)
        if not residual.is_zero:return None
        cone.append({'ambient_wall_one_based':l+1,'alpha_numerator':coeff(alpha),'beta_numerator':coeff(beta),'denominator':coeff(D)})
    if N.is_zero:return None
    P=pp(1) if obtuse else N*(4*N-den)*(2*N-den)*(4*N-3*den)*(16*N*N-12*N*den+den*den)
    # D*h - A*(qe*r_a-v_a*e) - B*(qe*r_b-v_b*e).
    point=[D*int(h[l])-A*(qe*int(W[a,l])-v[a]*pp(int(e0[l])+t*int(e1[l])))-B*(qe*int(W[b,l])-v[b]*pp(int(e0[l])+t*int(e1[l]))) for l in range(22)]
    hp=D*int(h@Q@h)-qe*(A*cp[a]+B*cp[b]);square=D*hp
    if not positive(-hp) or not positive(-square):return None
    return {'ridge_walls_one_based':[a+1,b+1],'projected_gram':[coeff(aa),coeff(ab),coeff(bb)],
            'point_coordinate_polynomials':[coeff(p) for p in point],
            'point_square':coeff(square),'point_center_pairing':coeff(hp),
            'all_ambient_pairings':[coeff(p) for p in vals],'zero_walls_one_based':[l+1 for l in zeros],
            'redundant_zero_wall_cones':cone,'angle_rejection':'obtuse_everywhere' if obtuse else 'necessary_Coxeter_polynomial','necessary_angle_polynomial':coeff(P)}

_FACE_CACHE={}
def ambient_face(Q,W,G,h,ids):
    key=tuple(sorted(set(ids)))
    if len(key)!=3:return None
    if key in _FACE_CACHE:return _FACE_CACHE[key]
    M=sp.Matrix(G[list(key),:][:,list(key)])
    if M[:2,:2].det()<=0 or M.det()<=0:return None
    weights=M.inv()*sp.Matrix([int(w@Q@h) for w in W[list(key)]])
    z=sp.Matrix(h)-sp.Matrix(W[list(key)].T)*weights
    den=sp.ilcm(*[x.q for x in z]);z=source.base.primitive(np.array(z,dtype=object).reshape(-1)*int(den))
    vals=W@Q@z
    if [l for l,v in enumerate(vals) if v==0]!=list(key) or any(v>0 for v in vals):return None
    assert z@Q@z<0 and z@Q@h<0
    _FACE_CACHE[key]=z
    return z

def interpolated_constraint(Q,W,G,h,sd,a,b,i,j,k):
    if i not in (a,b) and j not in (a,b):return None
    remaining=j if i in (a,b) else i
    z1=ambient_face(Q,W,G,h,(a,b,remaining));z2=ambient_face(Q,W,G,h,(a,b,k))
    if z1 is None or z2 is None:return None
    cp,e0,e1,qe,v,gr=sd
    ep1=pp(int(e0@Q@z1)+t*int(e1@Q@z1));ep2=pp(int(e0@Q@z2)+t*int(e1@Q@z2))
    if not positive(ep1) or not positive(-ep2):return None
    aa,bb,ab=gr(a,a),gr(b,b),gr(a,b);D=aa*bb-ab*ab
    if not all(positive(p) for p in (qe,aa,bb,D)):return None
    obtuse=positive(ab)
    N=ab*ab;den=aa*bb
    if not obtuse and (not positive(-ab,False) or not positive(3*den-4*N)):return None
    if N.is_zero:return None
    point=[ep1*int(z2[l])-ep2*int(z1[l]) for l in range(22)]
    vals=[ep1*int(w@Q@z2)-ep2*int(w@Q@z1) for w in W]
    zeros=[l for l,p in enumerate(vals) if p.is_zero]
    if zeros!=[a,b] or any(not positive(-p) for l,p in enumerate(vals) if l not in zeros):return None
    square=ep1*ep1*int(z2@Q@z2)-2*ep1*ep2*int(z1@Q@z2)+ep2*ep2*int(z1@Q@z1)
    hp=ep1*int(z2@Q@h)-ep2*int(z1@Q@h)
    if not positive(-square) or not positive(-hp):return None
    P=pp(1) if obtuse else N*(4*N-den)*(2*N-den)*(4*N-3*den)*(16*N*N-12*N*den+den*den)
    return {'ridge_walls_one_based':[a+1,b+1],'projected_gram':[coeff(aa),coeff(ab),coeff(bb)],
            'point_coordinate_polynomials':[coeff(p) for p in point],'point_square':coeff(square),'point_center_pairing':coeff(hp),
            'all_ambient_pairings':[coeff(p) for p in vals],'zero_walls_one_based':[l+1 for l in zeros],
            'redundant_zero_wall_cones':[],'angle_rejection':'obtuse_everywhere' if obtuse else 'necessary_Coxeter_polynomial',
            'necessary_angle_polynomial':coeff(P),'point_construction':'interpolated_ambient_faces',
            'ambient_face_points':[list(map(int,z1)),list(map(int,z2))]}

def build():
    Q,W,G,h=source.base.load();raw=json.loads(source.OUTPUT.read_text());rows=[]
    cached=json.loads(OUTPUT.read_text())['half_pencils'] if OUTPUT.exists() else []
    previous={(tuple(r['positive_walls_one_based']),r['negative_wall_one_based']):r for r in cached}
    cached={key:r for key,r in previous.items() if r['status']=='EXCLUDED_CONTINUOUS_HALF_PENCIL'}
    for number,old in enumerate(raw['sections']):
        k=old['negative_wall_one_based']-1;i,j=[x-1 for x in old['positive_walls_one_based']]
        for swap in (False,True):
            if swap:i,j=j,i
            key=((i+1,j+1),k+1)
            if key in cached:
                rows.append(cached[key]);continue
            sd=setup(Q,W,G,h,i,j,k);found=[];gcd=pp(0);tried=0
            proposed=[]
            old=previous.get(key)
            if old is not None and old['possible_parameters_0_1']:
                for factor in old['possible_parameter_factors']:
                    poly=sp.Poly.from_list([sp.Rational(c) for c in factor['factor']],t)
                    for (lo,hi),multiplicity in poly.intervals(eps=sp.Rational(1,10**15)):
                        if hi>0 and lo<=1:proposed.extend(candidates(Q,W,G,h,sd,float((lo+hi)/2)))
            proposed.extend(candidates(Q,W,G,h,sd));proposed=list(dict.fromkeys(proposed))
            for a,b in proposed:
                tried+=1
                r=constraint(Q,W,G,h,sd,a,b)
                if r is None:r=interpolated_constraint(Q,W,G,h,sd,a,b,i,j,k)
                if r is not None:
                    found.append(r);gcd=sp.gcd(gcd,sp.Poly.from_list([sp.Rational(x) for x in r['necessary_angle_polynomial']],t,domain=sp.QQ))
                    nonzero=gcd
                    while not nonzero.is_zero and nonzero.eval(0)==0:nonzero=nonzero.exquo(pp(t))
                    if not nonzero.is_zero and nonzero.count_roots(0,1)==0:break
                if tried>=48 or len(found)>=6:break
            roots=[];remaining=None
            if not gcd.is_zero:
                nonzero=gcd
                while nonzero.eval(0)==0:nonzero=nonzero.exquo(pp(t))
                remaining=int(nonzero.count_roots(0,1))
                if remaining:
                    roots=[{'factor':[str(c) for c in factor.monic().all_coeffs()],'roots_0_1':int(factor.count_roots(0,1))} for factor,power in sp.factor_list(nonzero)[1] if factor.count_roots(0,1)]
            row={'positive_walls_one_based':[i+1,j+1],'negative_wall_one_based':k+1,'domain':'0<t<=1',
                 'candidate_pairs_tried':tried,'constraints':found,'gcd':coeff(gcd.monic()) if not gcd.is_zero else ['0'],
                 'possible_parameters_0_1':remaining,'possible_parameter_factors':roots,
                 'status':'EXCLUDED_CONTINUOUS_HALF_PENCIL' if remaining==0 else 'UNRESOLVED'}
            rows.append(row)
            print('pencil',number+1,swap,[i+1,j+1,k+1],row['status'],'constraints',len(found),'tried',tried,'roots',remaining,flush=True)
            if len(rows)%10==0:OUTPUT.write_text(json.dumps({'schema':'central-triple-pencils-v1','complete':False,'half_pencils':rows},indent=2)+'\n')
    out={'schema':'central-triple-pencils-v1','complete':True,'central_source':str(source.OUTPUT.relative_to(ROOT)),
         'normal_definition':'e(t)=u_i+t*u_j-(1+t)*u_k, u_l=8*r_l/(-(h,r_l)); t>0',
         'half_pencil_count':len(rows),'half_pencils':rows,
         'all_half_pencils_excluded':all(r['status']=='EXCLUDED_CONTINUOUS_HALF_PENCIL' for r in rows)}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(rows),sum(r['status']=='UNRESOLVED' for r in rows),flush=True)

if __name__=='__main__':build()
