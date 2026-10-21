#!/usr/bin/env python3
"""Noncentral algebraic sections chosen to repair a persistent triple ridge.

For each elliptic triple and each metrically symmetric positive pair, solve
exactly for e=r_i+r_j-t*r_k (t>0) so its persistent ridge has label2,3,4,5,6.
Every retained algebraic parameter has a strict timelike interior certificate.
Global failure needs an exact actual forbidden ridge, allowing certified
redundant zero inequalities. This is a bounded angle alphabet and pencil set.
"""
from __future__ import annotations
from itertools import combinations
import json
import numpy as np
import sympy as sp
import fresh_h20_triple_central_20261013 as source
from fresh_h20_global_section_parameters_20261013 import Field,coeffs

ROOT=source.ROOT
OUTPUT=ROOT/'data/enumerations/fresh_h20_triple_angle_targets_20261013.json'
t=sp.Symbol('t')
P=lambda x:sp.Poly(x,t,domain=sp.QQ)


def data(Q,W,G,h,i,j,k):
    e0=W[i]+W[j];e1=-W[k]
    e=[P(int(e0[l])+t*int(e1[l])) for l in range(22)]
    qe=P(int(e0@Q@e0)+2*t*int(e0@Q@e1)+t*t*int(e1@Q@e1))
    v=[P(int(w@Q@e0)+t*int(w@Q@e1)) for w in W]
    def gr(a,b):return qe*int(G[a,b])-v[a]*v[b]
    return e,qe,v,gr


def pdot(Q,x,y,F):
    out=P(0)
    if all(p.degree()<=0 for p in y):
        raw=[p.eval(0) for p in y]
        dual=[sum(int(Q[a,b])*raw[b] for b in range(22)) for a in range(22)]
        return F.rem(sum((dual[a]*x[a] for a in range(22)),P(0)))
    for a in range(22):
        for b in range(22):
            if Q[a,b]:out+=int(Q[a,b])*x[a]*y[b]
    return F.rem(out)


def interior(Q,W,G,h,i,j,k,F,sd):
    e,qe,v,gr=sd
    he=sum((int(h[l])*(sum((int(Q[l,m])*e[m] for m in range(22)),P(0))) for l in range(22)),P(0))
    C=[F.rem(qe*int(h[l])-he*e[l]) for l in range(22)]
    def valid(C):
        if F.sign(pdot(Q,C,C,F))>=0 or F.sign(pdot(Q,C,[P(int(x)) for x in h],F))>=0:return False
        return all(F.sign(pdot(Q,C,[P(int(x)) for x in w],F))<0 for w in W)
    if valid(C):return C,'projected_original_center'
    ids=[i,j,k];M=sp.Matrix(G[ids,:][:,ids]);weights=M.inv()*sp.Matrix([int(w@Q@h) for w in W[ids]])
    zz=sp.Matrix(h)-sp.Matrix(W[ids].T)*weights
    weights=M.inv()*sp.Matrix([-t,-t,-2]);yy=sp.Matrix(W[ids].T)*weights
    for exponent in range(60):
        C=[F.rem(P(zz[l]+yy[l]/2**exponent)) for l in range(22)]
        if valid(C):return C,'persistent_face_perturbation_2^(-%d)'%exponent
    raise AssertionError('constructive strict center did not settle')


def candidates(Q,W,h,sd,C,F):
    value=float((F.lo+F.hi)/2);Qf=np.array(Q,dtype=float);Wf=np.array(W,dtype=float)
    e=np.array([float(x.eval(value)) for x in sd[0]]);q=e@Qf@e
    R=Wf-np.outer(Wf@Qf@e,e)/q
    center=np.array([float(x.eval(value)) for x in C]);H=R@Qf@R.T;cp=R@Qf@center
    for a,b in combinations(range(210),2):
        aa,bb,ab=H[a,a],H[b,b],H[a,b];D=aa*bb-ab*ab
        if aa<1e-8 or bb<1e-8 or D<1e-8:continue
        ratio=ab*ab/(aa*bb)
        if ab<=0 and (ratio>=.74999 or min(abs(ratio-v) for v in (0,.25,.5,(3+5**.5)/8,(3-5**.5)/8))<1e-8):continue
        A=bb*cp[a]-ab*cp[b];B=aa*cp[b]-ab*cp[a]
        vals=D*cp-A*H[:,a]-B*H[:,b];scale=max(1,abs(vals).max())
        if vals.max()<1e-8*scale:yield a,b


def witness(Q,W,h,sd,C,F,a,b,cp=None):
    e,qe,v,gr=sd;aa,bb,ab=map(F.rem,(gr(a,a),gr(b,b),gr(a,b)));D=F.rem(aa*bb-ab*ab)
    if min(F.sign(aa),F.sign(bb),F.sign(D))<=0:return None
    N=F.rem(ab*ab);den=F.rem(aa*bb)
    if F.sign(ab)>0:reason='obtuse'
    elif F.sign(3*den-4*N)>0 and F.sign(N)!=0 and all(F.sign(p)!=0 for p in (4*N-den,2*N-den,16*N*N-12*N*den+den*den)):
        reason='acute_cosine_squared_below_3/4_but_not_Coxeter'
    else:return None
    if cp is None:cp=[pdot(Q,C,[P(int(x)) for x in w],F) for w in W]
    A=F.rem(bb*cp[a]-ab*cp[b]);B=F.rem(aa*cp[b]-ab*cp[a])
    vals=[F.rem(D*cp[l]-A*gr(l,a)-B*gr(l,b)) for l in range(210)];signs=[F.sign(p) for p in vals]
    if max(signs)>0:return None
    zeros=[l for l,s in enumerate(signs) if s==0]
    if a not in zeros or b not in zeros:return None
    cone=[]
    for l in zeros:
        if l in (a,b):continue
        alpha=F.rem(bb*gr(l,a)-ab*gr(l,b));beta=F.rem(aa*gr(l,b)-ab*gr(l,a))
        if min(F.sign(alpha),F.sign(beta))<0:return None
        if not F.rem(D*gr(l,l)-alpha*gr(l,a)-beta*gr(l,b)).is_zero:return None
        cone.append({'wall_one_based':l+1,'alpha_numerator':coeffs(alpha),'beta_numerator':coeffs(beta),'denominator':coeffs(D)})
    point=[F.rem(D*C[l]-A*(qe*int(W[a,l])-v[a]*e[l])-B*(qe*int(W[b,l])-v[b]*e[l])) for l in range(22)]
    square=pdot(Q,point,point,F);hp=pdot(Q,point,[P(int(x)) for x in h],F)
    assert F.sign(square)<0 and F.sign(hp)<0 and pdot(Q,point,e,F).is_zero
    for l,w in enumerate(W):assert pdot(Q,point,[P(int(x)) for x in w],F)==vals[l]
    return {'ridge_walls_one_based':[a+1,b+1],'projected_gram':[coeffs(aa),coeffs(ab),coeffs(bb)],
            'point_coordinate_polynomials':[coeffs(p) for p in point],'point_square':coeffs(square),'point_h_pairing':coeffs(hp),
            'zero_walls_one_based':[l+1 for l in zeros],'redundant_zero_wall_cones':cone,'angle_rejection':reason}


def pencil(Q,W,G,h,i,j,k):
    sd=data(Q,W,G,h,i,j,k);e,qe,v,gr=sd
    num=gr(i,j)**2;den=gr(i,i)*gr(j,j)
    equations={2:num,3:4*num-den,4:2*num-den,5:16*num*num-12*num*den+den*den,6:4*num-3*den}
    factors={}
    for m,poly in equations.items():
        assert not poly.is_zero
        for f,power in sp.factor_list(poly)[1]:f=f.monic();factors[tuple(coeffs(f))]=f
    rows=[];conjugates=0
    for key,polynomial in sorted(factors.items()):
        for (lo,hi),multiplicity in polynomial.intervals(eps=sp.Rational(1,10**18)):
            if hi<=0:continue
            if lo==hi:lo-=sp.Rational(1,10**20);hi+=sp.Rational(1,10**20)
            assert lo>0
            F=Field(polynomial,lo,hi);labels=[]
            if min(F.sign(gr(i,i)),F.sign(gr(j,j)))<=0:continue
            if F.sign(gr(i,j))>0:continue
            for m,poly in equations.items():
                if F.sign(poly)!=0:continue
                if m==5 and F.sign(2*num-den)<=0:continue
                labels.append(m)
            if not labels:conjugates+=1;continue
            C,method=interior(Q,W,G,h,i,j,k,F,sd);bad=None
            cp=[pdot(Q,C,[P(int(x)) for x in w],F) for w in W]
            for a,b in candidates(Q,W,h,sd,C,F):
                bad=witness(Q,W,h,sd,C,F,a,b,cp)
                if bad is not None:break
            row={'minimal_polynomial':list(key),'isolating_interval':[str(F.lo),str(F.hi)],
                 'approximate_parameter':float((F.lo+F.hi)/2),'repaired_ridge_labels':labels,
                 'interior_point':[coeffs(p) for p in C],'interior_method':method,
                 'actual_forbidden_ridge':bad,'status':'NON_COXETER_ACTUAL_RIDGE' if bad else 'UNRESOLVED'}
            rows.append(row);print([i+1,j+1,k+1],labels,row['approximate_parameter'],row['status'],flush=True)
    return {'positive_walls_one_based':[i+1,j+1],'negative_wall_one_based':k+1,'target_labels':[2,3,4,5,6],
            'target_polynomials':{str(m):coeffs(p) for m,p in equations.items()},'discarded_conjugate_roots':conjugates,'candidates':rows}


def build():
    Q,W,G,h=source.base.load();orbits=json.loads(source.ORBIT_OUTPUT.read_text())['orbits'];rows=[]
    for orbit in orbits:
        ids=[x-1 for x in orbit['representative']];M=sp.Matrix(G[ids,:][:,ids])
        if M[:2,:2].det()<=0 or M.det()<=0:continue
        for k in ids:
            i,j=[l for l in ids if l!=k]
            if G[i,i]!=G[j,j] or G[i,k]!=G[j,k]:continue
            rows.append(pencil(Q,W,G,h,i,j,k))
            OUTPUT.write_text(json.dumps({'schema':'symmetric-triple-angle-targets-v1','complete':False,'pencils':rows},indent=2)+'\n')
    out={'schema':'symmetric-triple-angle-targets-v1','complete':True,'source':str(source.base.SOURCE.relative_to(ROOT)),
         'normal_definition':'r_i+r_j-t*r_k, t>0','scope':'All17metrically symmetric positive-pair choices among15elliptic triple orbits, target labels2,3,4,5,6. Not all triples, coefficient ratios, or labels.',
         'pencils':rows,'candidate_count':sum(len(r['candidates']) for r in rows)}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(rows),out['candidate_count'],flush=True)

if __name__=='__main__':build()
