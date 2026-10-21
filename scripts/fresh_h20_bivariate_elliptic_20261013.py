#!/usr/bin/env python3
"""Exact two-parameter obstruction for every elliptic three-wall section.

The section normal is ri+s*rj-t*rk with independent s,t>0.  Actual ridges
are certified by interpolating strict ambient triple-face points.  This is
not a statement about affine/Lorentzian triple spans or higher support.
"""
from __future__ import annotations
import json
import numpy as np
import sympy as sp
import fresh_h20_triple_central_20261013 as source
from fresh_h20_triple_pencils_20261013 import ambient_face
from fresh_h20_global_section_parameters_20261013 import Field, coeffs
import fresh_h20_triple_angle_targets_20261013 as old

ROOT=source.ROOT
OUTPUT=ROOT/'data/enumerations/fresh_h20_bivariate_elliptic_20261013.json'
s,t=sp.symbols('s t')

def enc(p):
    return [[int(a),int(b),str(c)] for (a,b),c in sp.Poly(p,s,t).terms() if c]

def positive(p):
    p=sp.Poly(p,s,t)
    return not p.is_zero and all(c>=0 for c in p.coeffs())

def data(G,i,j,k):
    ids=(i,j,k);c=(1,s,-t)
    q=sp.expand(sum(a*b*int(G[x,y]) for a,x in zip(c,ids) for b,y in zip(c,ids)))
    v=[sp.expand(sum(a*int(G[l,x]) for a,x in zip(c,ids))) for l in range(210)]
    def gr(a,b):return sp.expand(q*int(G[a,b])-v[a]*v[b])
    return q,v,gr

def ridge(Q,W,G,h,i,j,k,a,l,sd):
    assert a in (i,j)
    other=j if a==i else i
    z1=ambient_face(Q,W,G,h,(a,other,l))
    z2=ambient_face(Q,W,G,h,(a,k,l))
    if z1 is None or z2 is None:return None
    q,v,gr=sd
    def ep(z):return sp.expand(int(W[i]@Q@z)+s*int(W[j]@Q@z)-t*int(W[k]@Q@z))
    u=ep(z1);v2=ep(z2)
    assert positive(u) and positive(-v2)
    point=[sp.expand(u*int(z2[d])-v2*int(z1[d])) for d in range(22)]
    square=sp.expand(u*u*int(z2@Q@z2)-2*u*v2*int(z1@Q@z2)+v2*v2*int(z1@Q@z1))
    assert positive(-square)
    return {'ridge_walls_one_based':[a+1,l+1],
            'ambient_face_points':[list(map(int,z1)),list(map(int,z2))],
            'face_normal_pairings':[enc(u),enc(v2)],
            'point_coordinate_polynomials':[enc(x) for x in point],
            'point_square':enc(square),
            'projected_gram':[enc(gr(a,a)),enc(gr(a,l)),enc(gr(l,l))]}

def global_single(Q,W,G,h,i,j,k,sd):
    q,v,gr=sd;seen=set()
    for a,other in ((i,j),(j,i)):
        for l in range(210):
            if l in (i,j,k):continue
            key=(a,*(int(G[l,x]) for x in (l,i,j,k)))
            if key in seen:continue
            seen.add(key)
            aa,bb,ab=gr(a,a),gr(l,l),gr(a,l)
            N=sp.expand(ab*ab);den=sp.expand(aa*bb)
            reason=None
            if positive(ab):reason='obtuse_on_positive_quadrant'
            elif positive(-ab) and positive(den-4*N):reason='0<cosine_squared<1/4'
            elif positive(-ab) and positive(4*N-den) and positive(den-2*N):reason='1/4<cosine_squared<1/2'
            if reason is None or not all(positive(p) for p in (aa,bb,den-N)):continue
            proof=ridge(Q,W,G,h,i,j,k,a,l,sd)
            if proof is None:continue
            proof['reason']=reason
            return proof

def leaves(Q,W,G,h,i,j,k,sd):
    found=[]
    for a,b in ((i,j),(j,i)):
        proof=None
        for l in range(210):
            if l in (i,j,k) or G[l,a]!=0 or G[l,k]!=0 or G[l,b]>=0:continue
            proof=ridge(Q,W,G,h,i,j,k,a,l,sd)
            if proof is not None:break
        if proof is None:return None
        proof['conditional_reason']='obtuse_when_own_wall_normal_pairing_positive'
        found.append(proof)
    return found

def candidate(Q,W,G,h):
    # Independent coefficients, determined simultaneously by a right-angle
    # equation and a pi/8 equation.  Here old.t is the algebraic parameter s.
    i,j,k=0,47,19;u=old.t
    polynomial=old.P(23*u**4-52*u**3+44*u**2-16*u+2).monic()
    intervals=[(lo,hi) for (lo,hi),_ in polynomial.intervals(eps=sp.Rational(1,10**20))
               if lo>0 and hi<sp.Rational(1,2)]
    assert len(intervals)==1
    F=Field(polynomial,*intervals[0]);S=old.P(u);T=old.P(1-2*u)
    e=[old.P(int(W[i,d])+S.as_expr()*int(W[j,d])-T.as_expr()*int(W[k,d])) for d in range(22)]
    q=old.pdot(Q,e,e,F)
    v=[old.pdot(Q,e,[old.P(int(x)) for x in w],F) for w in W]
    def gr(a,b):return F.rem(q*int(G[a,b])-v[a]*v[b])
    sd=(e,q,v,gr)
    he=old.pdot(Q,e,[old.P(int(x)) for x in h],F)
    assert F.sign(he)!=0 and F.sign(S)>0 and F.sign(T)>0
    z=ambient_face(Q,W,G,h,(i,j,k));assert z is not None
    weights=sp.Matrix(G[[i,j,k],:][:,[i,j,k]]).inv()*sp.Matrix([-T.as_expr(),-T.as_expr(),-(1+S.as_expr())])
    y=sp.Matrix(W[[i,j,k]].T)*weights
    C=None
    for exponent in range(80):
        point=[F.rem(old.P(int(z[d])+y[d]/2**exponent)) for d in range(22)]
        if F.sign(old.pdot(Q,point,point,F))>=0:continue
        if any(F.sign(old.pdot(Q,point,[old.P(int(x)) for x in w],F))>=0 for w in W):continue
        C=point;break
    assert C is not None and old.pdot(Q,C,e,F).is_zero
    # Genuine repaired triple ridge and independent interpolated right ridge.
    N=F.rem(gr(i,j)**2);den=F.rem(gr(i,i)*gr(j,j))
    assert F.rem(16*N*N-16*N*den+2*den*den).is_zero
    assert F.sign(4*N-3*den)>0 and F.sign(gr(i,j))<0
    assert gr(j,21).is_zero
    bad=old.witness(Q,W,h,sd,C,F,i,18)
    assert bad is not None
    return {'status':'TWO_ACTUAL_ANGLES_REPAIRED_BUT_NON_COXETER',
            'positive_walls_one_based':[i+1,j+1],'negative_wall_one_based':k+1,
            'minimal_polynomial':coeffs(polynomial),'isolating_interval':[str(F.lo),str(F.hi)],
            's_polynomial':coeffs(S),'t_polynomial':coeffs(T),
            'approximate_s':float((F.lo+F.hi)/2),'approximate_t':float(1-F.lo-F.hi),
            'section_normal_coordinate_polynomials':[coeffs(x) for x in e],
            'section_normal_center_pairing':coeffs(he),
            'strict_interior':[coeffs(x) for x in C],
            'strict_interior_perturbation_exponent':exponent,
            'repaired_ridges':[{'walls_one_based':[i+1,j+1],'label':8,'point':list(map(int,z)),
                               'redundant_zero_wall_one_based':k+1},
                              {'walls_one_based':[j+1,22],'label':2,
                               'point_source':'exceptional_positive_branch_interpolated_point'}],
            'actual_forbidden_ridge':bad}

def build():
    Q,W,G,h=source.base.load();orbits=json.loads(source.ORBIT_OUTPUT.read_text())['orbits']
    rows=[];elliptic=[]
    for orbit in orbits:
        ids=[x-1 for x in orbit['representative']];M=sp.Matrix(G[ids,:][:,ids])
        if M[:2,:2].det()<=0 or M.det()<=0:continue
        elliptic.append(orbit)
        for k in ids:
            i,j=[x for x in ids if x!=k];sd=data(G,i,j,k)
            entry={'positive_walls_one_based':[i+1,j+1],'negative_wall_one_based':k+1,
                   'triple_one_based':orbit['representative'],'triple_orbit_size':orbit['size'],
                   'triple_face_point':list(map(int,ambient_face(Q,W,G,h,(i,j,k))))}
            proof=global_single(Q,W,G,h,i,j,k,sd)
            if proof is not None:entry.update(kind='one_global_forbidden_ridge',ridges=[proof])
            else:
                proofs=leaves(Q,W,G,h,i,j,k,sd)
                if proofs is not None:entry.update(kind='two_positive_covector_branches',ridges=proofs)
                else:
                    assert (i,j,k)==(0,47,19)
                    proofs=[ridge(Q,W,G,h,i,j,k,*p,sd) for p in ((j,21),(j,18),(i,18))]
                    assert all(p is not None for p in proofs)
                    entry.update(kind='covector_sign_and_zero_line',ridges=proofs,
                                 branch_function=enc(2*s+t-1),zero_line='t=1-2*s, 0<s<1/2',
                                 zero_line_forbidden_cosine_squared='(2*s-1)^2/(2*(5*s^2-6*s+2))')
            entry['status']='ALL_POSITIVE_S_T_EXCLUDED';rows.append(entry)
            print([i+1,j+1,k+1],entry['kind'],flush=True)
    assert len(elliptic)==15 and len(rows)==45
    out={'schema':'h20-bivariate-elliptic-v1','complete':True,
         'source':str(source.base.SOURCE.relative_to(ROOT)),
         'orbit_source':str(source.ORBIT_OUTPUT.relative_to(ROOT)),
         'normal_definition':'e=ri+s*rj-t*rk, independently s,t>0',
         'elliptic_triple_orbits':elliptic,'elliptic_triple_count':sum(x['size'] for x in elliptic),
         'families':rows,'two_angle_candidate':candidate(Q,W,G,h),
         'scope':'All interior H20 sections whose normal lies in an elliptic span of at most three original Borcherds walls. Two-wall boundary cases use the prior separate theorem. No claim for affine/Lorentzian three-wall spans or larger support.'}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    print('DONE',out['elliptic_triple_count'],len(rows),flush=True)

if __name__=='__main__':build()
