#!/usr/bin/env python3
"""Exclude every positive wall-pair parameter by exact actual-ridge proofs.

All 21,945 pencils e=r_i-t*r_j, t>0, reduce to nine symmetry orbits. Six
elliptic orbits have persistent triple-face angle equations with no common
positive root. Three other orbits have explicit polynomial ridge points and
forbidden angle intervals; one uses reciprocal-parameter symmetry.

The archived 28 angle-targeted algebraic candidates are exploratory history,
superseded by this universal wall-pair-pencil obstruction. No conclusion is
claimed for normals involving three or more walls, or general H20 polytopes.
"""
from __future__ import annotations
import argparse
from functools import reduce
from itertools import combinations
import json
from pathlib import Path
import numpy as np
import sympy as sp
import fresh_h20_global_interior_sections_20261013 as central

ROOT=central.ROOT
OUTPUT=ROOT/'data/enumerations/fresh_h20_global_section_parameters_20261013.json'
t=sp.Symbol('t');z=sp.Symbol('z')

def coeffs(p):
    return [str(x) for x in sp.Poly(p,t,domain=sp.QQ).all_coeffs()]

def polynomial(c):return sp.Poly.from_list([sp.Rational(x) for x in c],t,domain=sp.QQ)

def persistent_test(Q,W,G,h,i,j):
    a,b,c=map(int,(G[i,i],G[j,j],G[i,j]));delta=a*b-c*c
    assert delta>0
    rows=[];seen=set();gcd=sp.Poly(0,t,domain=sp.QQ)
    for k in range(210):
        if k in (i,j):continue
        d,u,v=map(int,(G[k,k],G[i,k],G[j,k]));M=sp.Matrix([[a,c,u],[c,b,v],[u,v,d]])
        if M.det()<=0 or (d,u,v) in seen:continue
        seen.add((d,u,v))
        n=sp.Poly(((a*v-c*u)+t*(b*u-c*v))**2,t,domain=sp.QQ)
        den=sp.Poly(delta*(d*(a-2*t*c+t*t*b)-(u-t*v)**2),t,domain=sp.QQ)
        bound=3*den-4*n
        if bound.is_zero or not all(x>=0 for x in bound.all_coeffs()):continue
        if n.is_zero:continue
        P=n*(4*n-den)*(2*n-den)*(4*n-3*den)*(16*n*n-12*n*den+den*den)
        gcd=sp.gcd(gcd,P)
        ids=[i,j,k];cp=sp.Matrix([int(W[x]@Q@h) for x in ids]);weights=M.inv()*cp
        point=sp.Matrix(h)-sp.Matrix(W[ids].T)*weights
        denominator=sp.ilcm(*[x.q for x in point]);point=central.primitive(np.array(point,dtype=object).reshape(-1)*int(denominator))
        vals=W@Q@point
        assert [l for l,x in enumerate(vals) if x==0]==sorted(ids) and all(x<=0 for x in vals)
        assert point@Q@point<0 and point@Q@h<0
        rows.append({'third_wall_one_based':k+1,'triple_gram':[[int(x) for x in row] for row in M.tolist()],
                     'cosine_squared_numerator':coeffs(n),'cosine_squared_denominator':coeffs(den),
                     'three_denominator_minus_four_numerator':coeffs(bound),
                     'necessary_angle_polynomial':coeffs(P),'persistent_face_point':list(map(int,point)),
                     'persistent_face_point_square':int(point@Q@point)})
    assert rows and not gcd.is_zero
    nonzero_gcd=gcd
    while nonzero_gcd.eval(0)==0:nonzero_gcd=nonzero_gcd.exquo(sp.Poly(t,t,domain=sp.QQ))
    positive_roots=nonzero_gcd.count_roots(0,sp.oo)
    assert positive_roots==0
    return {'pair_one_based':[i+1,j+1],'triple_constraints':rows,'gcd':coeffs(gcd.monic()),
            'gcd_positive_real_roots':int(positive_roots),'status':'ALL_POSITIVE_PARAMETERS_EXCLUDED'}

class Field:
    def __init__(self,p,lo,hi):
        self.p=sp.Poly(p,t,domain=sp.QQ).monic();self.lo=sp.Rational(lo);self.hi=sp.Rational(hi)
        assert self.lo<self.hi
        assert self.p.count_roots(self.lo,self.hi)==1
    def rem(self,v):return sp.Poly(v,t,domain=sp.QQ).rem(self.p)
    def interval(self,p):
        lo=hi=sp.Rational(0)
        for c in p.all_coeffs():
            v=[lo*self.lo,lo*self.hi,hi*self.lo,hi*self.hi]
            lo=min(v)+c;hi=max(v)+c
        return lo,hi
    def sign(self,v):
        p=self.rem(v)
        if p.is_zero:return 0
        for _ in range(2000):
            lo,hi=self.interval(p)
            if lo>0:return 1
            if hi<0:return -1
            self.lo,self.hi=sp.polys.polytools.refine_root(self.p,self.lo,self.hi,steps=4)
        raise RuntimeError('root sign refinement did not settle')

def section_data(G,h,W,Q,i,j):
    assert int(G[i,i])==int(G[j,j])==4
    a,b,c=map(int,(G[i,i],G[j,j],G[i,j]));qe=sp.Poly(a-2*c*t+b*t*t,t,domain=sp.QQ)
    values=[sp.Poly(int(G[l,i])-t*int(G[l,j]),t,domain=sp.QQ) for l in range(210)]
    cp=[sp.Poly((1+t)*int(W[l]@Q@h)+2*int(G[l,i])+2*t*int(G[l,j]),t,domain=sp.QQ) for l in range(210)]
    def gram(a,b):return qe*int(G[a,b])-values[a]*values[b]
    return qe,values,cp,gram

def numeric_bad_pair(Q,W,G,h,i,j,value):
    q=np.array(Q,dtype=float);w=np.array(W,dtype=float);hh=np.array(h,dtype=float)
    e=w[i]-value*w[j];qe=e@q@e;R=w-np.outer(w@q@e,e)/qe
    center=(1+value)*hh+2*w[i]+2*value*w[j]
    H=R@q@R.T;cp=R@q@center
    for a,b in combinations(range(210),2):
        aa,bb,ab=H[a,a],H[b,b],H[a,b];D=aa*bb-ab*ab
        if aa<1e-9 or bb<1e-9 or D<1e-9:continue
        ratio=ab*ab/(aa*bb)
        if ab<0 and (ratio>=.749999 or min(abs(ratio-v) for v in (0,.25,.5,(3+5**.5)/8,(3-5**.5)/8))<1e-8):continue
        C=bb*cp[a]-ab*cp[b];E=aa*cp[b]-ab*cp[a]
        vals=D*cp-C*H[:,a]-E*H[:,b]
        scale=max(abs(vals).max(),1.)
        others=np.delete(vals,[a,b])
        if others.max() < -1e-8*scale:yield a,b

def exact_witness(Q,W,G,h,i,j,F,sd,a,b):
    qe,v,cp,gram=sd
    aa,bb,ab=map(F.rem,(gram(a,a),gram(b,b),gram(a,b)));D=F.rem(aa*bb-ab*ab)
    if min(F.sign(aa),F.sign(bb),F.sign(D))<=0:return None
    num=F.rem(ab*ab);den=F.rem(aa*bb)
    if F.sign(ab)>0:reason='obtuse'
    elif F.sign(3*den-4*num)>0 and F.sign(num)!=0 and all(F.sign(p)!=0 for p in (4*num-den,2*num-den,16*num*num-12*num*den+den*den)):
        reason='cosine_squared_below_3/4_and_not_any_Coxeter_value_in_that_range'
    else:return None
    A=F.rem(bb*cp[a]-ab*cp[b]);B=F.rem(aa*cp[b]-ab*cp[a])
    vals=[F.rem(D*cp[l]-A*gram(l,a)-B*gram(l,b)) for l in range(210)]
    signs=[F.sign(x) for x in vals]
    if [l for l,s in enumerate(signs) if s==0]!=[a,b] or max(signs)>0:return None
    e=[sp.Poly(int(W[i,k])-t*int(W[j,k]),t,domain=sp.QQ) for k in range(22)]
    center=[sp.Poly((1+t)*int(h[k])+2*int(W[i,k])+2*t*int(W[j,k]),t,domain=sp.QQ) for k in range(22)]
    point=[F.rem(D*center[k]-A*(qe*int(W[a,k])-v[a]*e[k])-B*(qe*int(W[b,k])-v[b]*e[k])) for k in range(22)]
    def dot(x,y):
        out=sp.Poly(0,t,domain=sp.QQ)
        for k in range(22):
            for l in range(22):
                if Q[k,l]:out+=int(Q[k,l])*x[k]*y[l]
        return F.rem(out)
    q=dot(point,point);hp=dot(point,[sp.Poly(int(x),t,domain=sp.QQ) for x in h])
    assert F.sign(q)<0 and F.sign(hp)<0 and dot(point,e).is_zero
    # Independent reconstruction of all ambient wall pairings from coordinates.
    for l in range(210):
        r=[sp.Poly(int(x),t,domain=sp.QQ) for x in W[l]]
        assert dot(point,r)==vals[l]
    return {'ambient_ridge_walls_one_based':[a+1,b+1],'point_coordinate_polynomials':[coeffs(x) for x in point],
            'point_square':coeffs(q),'center_pairing':coeffs(hp),'projected_gram_polynomials':[coeffs(x) for x in (aa,ab,bb)],
            'angle_rejection':reason,'all_other_208_ambient_wall_pairings_strictly_negative':True}

def targeted_test(Q,W,G,h,row):
    i,j=[x-1 for x in row['pair_one_based']]
    a,b=[x[0]-1 for x in row['actual_forbidden_ridge']['normal_ambient_wall_labels_one_based']]
    sd=section_data(G,h,W,Q,i,j);qe,v,cp,gram=sd
    num=gram(a,b)**2;den=gram(a,a)*gram(b,b)
    factors={};target_rows=[]
    for m in range(2,13):
        minpoly=sp.Poly(sp.minimal_polynomial(sp.cos(sp.pi/m)**2,z),z,domain=sp.QQ)
        degree=minpoly.degree();P=sp.Poly(0,t,domain=sp.QQ)
        for k,c in enumerate(reversed(minpoly.all_coeffs())):P+=c*num**k*den**(degree-k)
        assert not P.is_zero
        rowf=[]
        for p,power in sp.factor_list(P)[1]:
            p=p.monic();key=tuple(coeffs(p));factors[key]=p;rowf.append(list(key))
        target_rows.append({'m':m,'cosine_minimal_polynomial':str(minpoly.as_expr()),'parameter_polynomial':coeffs(P),'irreducible_factors':rowf})
    candidates=[]
    for key,p in sorted(factors.items()):
        intervals=p.intervals(eps=sp.Rational(1,10**18))
        for (lo,hi),multiplicity in intervals:
            if hi<=0:continue
            if lo==hi:
                # Rational roots can use an irreducible linear field with a
                # nondegenerate isolating interval as well.
                radius=sp.Rational(1,10**20);lo-=radius;hi+=radius
            if lo<=0:lo,hi=sp.polys.polytools.refine_root(p,lo,hi,eps=sp.Rational(1,10**30))
            assert lo>0
            F=Field(p,lo,hi);value=float((lo+hi)/2)
            witness=None
            for aa,bb in numeric_bad_pair(Q,W,G,h,i,j,value):
                witness=exact_witness(Q,W,G,h,i,j,F,sd,aa,bb)
                if witness is not None:break
            candidates.append({'minimal_polynomial':list(key),'isolating_interval':[str(F.lo),str(F.hi)],
                               'approximate_parameter':value,'actual_forbidden_ridge':witness,
                               'status':'NON_COXETER_ACTUAL_RIDGE' if witness else 'UNRESOLVED'})
            print('parameter',i+1,j+1,value,candidates[-1]['status'],flush=True)
    return {'pair_one_based':[i+1,j+1],'targeted_central_defect_walls_one_based':[a+1,b+1],
            'target_labels':list(range(2,13)),'target_equations':target_rows,'algebraic_candidates':candidates,
            'scope':'Superset: all positive parameters where the selected pair has squared cosine cos²(pi/m), m=2,...,12; conjugate squared-cosine roots and degenerate pairs may also be included.'}

def domain_poly(p,unit_interval):
    p=sp.Poly(p,t,domain=sp.QQ)
    if not unit_interval:return p
    return sp.Poly(sp.cancel((1+t)**p.degree()*p.as_expr().subs(t,t/(1+t))),t,domain=sp.QQ)

def positive_on_domain(p,unit_interval):
    q=domain_poly(p,unit_interval)
    assert not q.is_zero and all(x>=0 for x in q.all_coeffs())
    if unit_interval:assert sp.Poly(p,t).eval(1)>0

def direct_nonelliptic_test(Q,W,G,h,row,gens):
    i,j=[x-1 for x in row['pair_one_based']]
    a,b=[x[0]-1 for x in row['actual_forbidden_ridge']['normal_ambient_wall_labels_one_based']]
    unit_interval=int(G[i,j])==-8
    sd=section_data(G,h,W,Q,i,j);qe,v,cp,gram=sd
    aa,bb,ab=gram(a,a),gram(b,b),gram(a,b);D=aa*bb-ab*ab
    for p in (qe,aa,bb,D,-ab):positive_on_domain(p,unit_interval)
    A=bb*cp[a]-ab*cp[b];B=aa*cp[b]-ab*cp[a]
    vals=[D*cp[l]-A*gram(l,a)-B*gram(l,b) for l in range(210)]
    assert [l for l,p in enumerate(vals) if p.is_zero]==[a,b]
    for l,p in enumerate(vals):
        if l not in (a,b):positive_on_domain(-p,unit_interval)
    num=ab*ab;den=aa*bb
    bounds=(num,den-4*num) if unit_interval else (4*num-den,den-2*num)
    for p in bounds:positive_on_domain(p,unit_interval)
    e=[sp.Poly(int(W[i,k])-t*int(W[j,k]),t,domain=sp.QQ) for k in range(22)]
    center=[sp.Poly((1+t)*int(h[k])+2*int(W[i,k])+2*t*int(W[j,k]),t,domain=sp.QQ) for k in range(22)]
    point=[D*center[k]-A*(qe*int(W[a,k])-v[a]*e[k])-B*(qe*int(W[b,k])-v[b]*e[k]) for k in range(22)]
    def dot(x,y):return sum((int(Q[k,l])*x[k]*y[l] for k in range(22) for l in range(22) if Q[k,l]),sp.Poly(0,t,domain=sp.QQ))
    square=dot(point,point);hp=dot(point,[sp.Poly(int(x),t,domain=sp.QQ) for x in h])
    assert dot(point,e).is_zero
    for p in (-square,-hp):positive_on_domain(p,unit_interval)
    word=[]
    if unit_interval:
        seen={(i,j):[]};queue=[(i,j)]
        while (j,i) not in seen:
            assert queue
            old=queue.pop(0)
            for k,g in enumerate(gens):
                perm=g['wall_permutation'];nxt=(perm[old[0]],perm[old[1]])
                if nxt not in seen:seen[nxt]=seen[old]+[k];queue.append(nxt)
        word=seen[(j,i)]
    return {'pair_one_based':[i+1,j+1],'direct_domain':'0<t<=1' if unit_interval else 't>0',
            'swapping_generator_word':word,'ridge_ambient_walls_one_based':[a+1,b+1],
            'projected_gram_polynomials':[coeffs(p) for p in (aa,ab,bb)],
            'ridge_point_coordinate_polynomials':[coeffs(p) for p in point],
            'ridge_point_square':coeffs(square),'ridge_point_center_pairing':coeffs(hp),
            'all_ambient_wall_pairing_polynomials':[coeffs(p) for p in vals],
            'squared_cosine_simplified':str(sp.factor(num.as_expr()/den.as_expr())),
            'forbidden_open_interval':'(0,1/4)' if unit_interval else '(1/4,1/2)',
            'positive_angle_bound_polynomials':[coeffs(p) for p in bounds],
            'status':'ALL_POSITIVE_PARAMETERS_EXCLUDED'}


def build():
    Q,W,G,h=central.load();raw=json.loads(central.OUTPUT.read_text());persistent=[];targeted=[]
    for row in raw['sections']:
        i,j=[x-1 for x in row['pair_one_based']]
        if G[i,i]*G[j,j]-G[i,j]**2>0:
            persistent.append(persistent_test(Q,W,G,h,i,j));print('all-positive exclusion',i+1,j+1,flush=True)
        else:targeted.append(targeted_test(Q,W,G,h,row))
    out={'schema':'interior-H20-one-parameter-sections-v1','central_source':str(central.OUTPUT.relative_to(ROOT)),
         'parameter_convention':'e = r_i - t r_j, t>0; all sections meet the strict interior and have finite volume',
         'all_parameter_exclusions':persistent,'bounded_angle_targeted_tests':targeted,
         'nonelliptic_all_parameter_exclusions':[direct_nonelliptic_test(Q,W,G,h,r,raw['ambient_symmetry_generators']) for r in raw['sections'] if G[r['pair_one_based'][0]-1,r['pair_one_based'][0]-1]*G[r['pair_one_based'][1]-1,r['pair_one_based'][1]-1]-G[r['pair_one_based'][0]-1,r['pair_one_based'][1]-1]**2<=0]}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    print('ARCHIVED',OUTPUT,flush=True)

def verify():
    saved=json.loads(OUTPUT.read_text());Q,W,G,h=central.load();raw=json.loads(central.OUTPUT.read_text())
    lookup={tuple(r['pair_one_based']):r for r in raw['sections']}
    for row in saved['all_parameter_exclusions']:
        i,j=[x-1 for x in row['pair_one_based']]
        assert persistent_test(Q,W,G,h,i,j)==row
    for row in saved['bounded_angle_targeted_tests']:
        i,j=[x-1 for x in row['pair_one_based']];sd=section_data(G,h,W,Q,i,j)
        # Recompute factor/root census as well as every field certificate.
        rebuilt=targeted_test(Q,W,G,h,lookup[(i+1,j+1)])
        assert rebuilt==row
    for row in saved['nonelliptic_all_parameter_exclusions']:
        assert direct_nonelliptic_test(Q,W,G,h,lookup[tuple(row['pair_one_based'])],raw['ambient_symmetry_generators'])==row
    print('PASS: all-positive triple obstructions, complete bounded algebraic root census, and actual ridge certificates',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    if args.verify:verify()
    else:build()
