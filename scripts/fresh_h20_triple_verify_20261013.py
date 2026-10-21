#!/usr/bin/env python3
"""Independent stdlib/Fraction replay of central three-wall section proofs.

No numpy, sympy, graph package, or search module is imported. Polynomial root
counts use independently implemented positive pseudo-remainder Sturm chains.
"""
from __future__ import annotations
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd,lcm
from pathlib import Path
import fresh_h20_global_fraction_verify_20261013 as A

ROOT=A.ROOT
ORB=ROOT/'data/enumerations/fresh_h20_triple_orbits_20261013.json'
CENTRAL=ROOT/'data/enumerations/fresh_h20_triple_central_20261013.json'
PENCILS=ROOT/'data/enumerations/fresh_h20_triple_pencils_20261013.json'
TARGETS=ROOT/'data/enumerations/fresh_h20_triple_angle_targets_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_triple_verify_20261013.json'
Z=A.Z;T=A.T;ONE=A.ONE
add,sub,mul,scale,power,readp=A.add,A.sub,A.mul,A.scale,A.power,A.readp


def iprim(p):
    p=A.trim(p);d=lcm(*(x.denominator for x in p));v=[int(x*d) for x in p];g=gcd(*v)
    return tuple(x//g for x in v) if g else (0,)
def itrim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return tuple(p)
def prempos(a,b):
    a=iprim(a);b=iprim(b);assert b!=(0,)
    while len(a)>=len(b) and a!=(0,):
        shift=len(a)-len(b);c=a[-1]*(1 if b[-1]>0 else -1);v=[abs(b[-1])*x for x in a]
        for k,x in enumerate(b):v[k+shift]-=c*x
        a=iprim(itrim(v))
    return a
def gcdpoly(a,b):
    a,b=iprim(a),iprim(b)
    while b!=(0,):a,b=b,prempos(a,b)
    return scale(a,F(1,a[-1])) if a!=(0,) else Z

def evalp(p,x):
    y=F(0)
    for a in reversed(p):y=y*x+a
    return y
def sign(x):return (x>0)-(x<0)
def sturm(p):
    p=iprim(p);assert p!=(0,)
    if len(p)==1:return [p]
    result=[p,iprim([k*p[k] for k in range(1,len(p))])]
    while True:
        r=prempos(result[-2],result[-1])
        if r==(0,):break
        result.append(tuple(-x for x in r))
    return result

def variation(chain,x):
    signs=[sign(evalp(p,x)) for p in chain];signs=[s for s in signs if s]
    return sum(a!=b for a,b in zip(signs,signs[1:]))
def rootcount(p,lo,hi):
    assert lo<hi and evalp(p,lo)!=0 and evalp(p,hi)!=0
    chain=sturm(p);return variation(chain,lo)-variation(chain,hi)
def no_roots_unit(p):
    while len(p)>1 and p[0]==0:p=p[1:]
    assert evalp(p,0)!=0 and evalp(p,1)!=0 and rootcount(p,F(0),F(1))==0


def sturm_regressions():
    cases=[(A.trim((-F(1,2),1)),1),
           (A.trim((-F(1,2),0,1)),1),
           (mul(A.trim((-F(1,3),1)),A.trim((-F(2,3),1))),2),
           (power(A.trim((-F(1,2),1)),2),1),
           (mul(A.trim((1,0,1)),A.trim((1,1))),0),
           (scale(A.trim((-F(1,2),1)),-7),1)]
    for p,n in cases:assert rootcount(p,F(0),F(1))==n
    no_roots_unit(power(T,2))
    return len(cases)+1

def triples(perms):
    enc=lambda a,b,c:44100*a+210*b+c
    def dec(v):a,v=divmod(v,44100);b,c=divmod(v,210);return a,b,c
    remaining={enc(*x) for x in combinations(range(210),3)};rows=[]
    while remaining:
        first=min(remaining);seen={first};queue=[first]
        while queue:
            a,b,c=dec(queue.pop())
            for g in perms:
                nxt=enc(*sorted((g[a],g[b],g[c])))
                if nxt not in seen:seen.add(nxt);queue.append(nxt)
        assert seen<=remaining;remaining-=seen
        rows.append({'representative':[x+1 for x in dec(first)],'size':len(seen)})
    assert len(rows)==49 and sum(r['size'] for r in rows)==1521520
    return rows


def project_data(Q,W,G,h,i,j,k):
    U=[[(-8//A.dot(w,h,Q))*x for x in w] for w in W]
    e0=[x-y for x,y in zip(U[i],U[k])];e1=[x-y for x,y in zip(U[j],U[k])]
    e=[A.trim((a,b)) for a,b in zip(e0,e1)]
    qe=A.trim((A.dot(e0,e0,Q),2*A.dot(e0,e1,Q),A.dot(e1,e1,Q)))
    v=[A.trim((A.dot(w,e0,Q),A.dot(w,e1,Q))) for w in W]
    def gram(a,b):return sub(scale(qe,G[a][b]),mul(v[a],v[b]))
    return e,qe,v,gram


def check_central(data,Q,W,G,h,orbits):
    assert data['triple_count']==1521520 and data['triple_orbit_count']==49 and len(data['sections'])==147
    expected={(tuple(r['representative']),k) for r in orbits for k in r['representative']}
    actual={(tuple(r['triple_one_based']),r['negative_wall_one_based']) for r in data['sections']}
    assert actual==expected
    cp=[A.dot(w,h,Q) for w in W];U=[[(-8//cp[l])*x for x in w] for l,w in enumerate(W)]
    for row in data['sections']:
        i,j=[x-1 for x in row['positive_walls_one_based']];k=row['negative_wall_one_based']-1
        assert sorted((i+1,j+1,k+1))==row['triple_one_based']
        e=A.primitive([a+b-2*c for a,b,c in zip(U[i],U[j],U[k])]);q=A.dot(e,e,Q)
        assert e==row['section_normal'] and q==row['section_normal_square']>0 and A.dot(e,h,Q)==0
        proof=row['actual_forbidden_ridge'];labels=proof['normal_ambient_wall_labels_one_based']
        normals=[]
        for ids in labels:
            l=ids[0]-1;v=A.primitive([q*x-A.dot(W[l],e,Q)*y for x,y in zip(W[l],e)])
            normals.append(v)
            for index in ids:
                l=index-1;assert v==A.primitive([q*x-A.dot(W[l],e,Q)*y for x,y in zip(W[l],e)])
        r,s=normals;aa,bb,ab=A.dot(r,r,Q),A.dot(s,s,Q),A.dot(r,s,Q)
        assert proof['gram']==[[aa,ab],[ab,bb]] and aa>0 and bb>0 and aa*bb>ab*ab
        assert F(proof['cosine_squared'])==F(ab*ab,aa*bb)
        assert ab>0 or F(ab*ab,aa*bb) not in {F(0),F(1,4),F(1,2),F(3,4)}
        p=proof['point'];vals=[A.dot(w,p,Q) for w in W]
        assert A.dot(p,e,Q)==0 and A.dot(p,p,Q)==proof['point_square']<0 and A.dot(p,h,Q)<0
        assert [l+1 for l,x in enumerate(vals) if x==0]==sorted(sum(labels,[]))
        assert all(x<=0 for x in vals) and row['status']=='NON_COXETER_ACTUAL_RIDGE'


def point_products(point,Q,W):
    dual=[sum_polys([scale(point[j],Q[i][j]) for j in range(22) if Q[i][j]]) for i in range(22)]
    vals=[sum_polys([scale(dual[j],w[j]) for j in range(22) if w[j]]) for w in W]
    square=sum_polys([mul(point[j],dual[j]) for j in range(22)])
    return dual,vals,square

def sum_polys(polys):
    out=Z
    for p in polys:out=add(out,p)
    return out


def check_constraint(proof,Q,W,G,h,sd):
    e,qe,v,gram=sd;a,b=[x-1 for x in proof['ridge_walls_one_based']]
    aa,bb,ab=gram(a,a),gram(b,b),gram(a,b);D=sub(mul(aa,bb),mul(ab,ab))
    assert [aa,ab,bb]==[readp(x) for x in proof['projected_gram']]
    for p in (qe,aa,bb,D):A.positive(p,True)
    point=[readp(p) for p in proof['point_coordinate_polynomials']]
    dual,vals,square=point_products(point,Q,W)
    assert vals==[readp(p) for p in proof['all_ambient_pairings']]
    assert square==readp(proof['point_square']);A.positive(A.neg(square),True)
    hp=sum_polys([scale(dual[j],h[j]) for j in range(22)])
    assert hp==readp(proof['point_center_pairing']);A.positive(A.neg(hp),True)
    assert sum_polys([mul(dual[j],e[j]) for j in range(22)])==Z
    zeros=[l for l,p in enumerate(vals) if p==Z]
    assert [l+1 for l in zeros]==proof['zero_walls_one_based'] and a in zeros and b in zeros
    for l,p in enumerate(vals):
        if l not in zeros:A.positive(A.neg(p),True)
    cones=proof['redundant_zero_wall_cones']
    assert sorted(r['ambient_wall_one_based']-1 for r in cones)==sorted(set(zeros)-{a,b})
    for row in cones:
        l=row['ambient_wall_one_based']-1
        alpha=sub(mul(bb,gram(l,a)),mul(ab,gram(l,b)));beta=sub(mul(aa,gram(l,b)),mul(ab,gram(l,a)))
        assert alpha==readp(row['alpha_numerator']) and beta==readp(row['beta_numerator']) and D==readp(row['denominator'])
        for p in (alpha,beta):assert all(c>=0 for c in A.transform_unit(p)) and sum(p)>=0
        assert sub(sub(mul(D,gram(l,l)),mul(alpha,gram(l,a))),mul(beta,gram(l,b)))==Z
    N=mul(ab,ab);den=mul(aa,bb)
    if proof.get('angle_rejection')=='obtuse_everywhere':
        A.positive(ab,True);poly=ONE
    else:
        assert all(c>=0 for c in A.transform_unit(A.neg(ab))) and sum(ab)<=0
        A.positive(sub(scale(den,3),scale(N,4)),True)
        poly=N
        for p in (sub(scale(N,4),den),sub(scale(N,2),den),sub(scale(N,4),scale(den,3)),add(sub(scale(mul(N,N),16),scale(mul(N,den),12)),mul(den,den))):poly=mul(poly,p)
    assert poly==readp(proof['necessary_angle_polynomial'])
    return poly


def check_pencils(data,Q,W,G,h,orbits):
    assert data['complete'] and data['all_half_pencils_excluded'] and data['half_pencil_count']==294
    expected=set()
    for r in orbits:
        for k in r['representative']:
            i,j=[x for x in r['representative'] if x!=k];expected.add(((i,j),k));expected.add(((j,i),k))
    rows=data['half_pencils'];assert len(rows)==294
    assert {(tuple(r['positive_walls_one_based']),r['negative_wall_one_based']) for r in rows}==expected
    for number,row in enumerate(rows):
        i,j=[x-1 for x in row['positive_walls_one_based']];k=row['negative_wall_one_based']-1
        sd=project_data(Q,W,G,h,i,j,k);gcd=Z
        for proof in row['constraints']:gcd=gcdpoly(gcd,check_constraint(proof,Q,W,G,h,sd))
        assert gcd==readp(row['gcd']);no_roots_unit(gcd)
        assert row['status']=='EXCLUDED_CONTINUOUS_HALF_PENCIL' and row['possible_parameters_0_1']==0
        if number%30==0:print('verified continuous half-pencils',number+1,flush=True)


def main():
    sturm_checks=sturm_regressions()
    Q,W,G,h=A.source();old=json.loads(A.CENTRAL.read_text());perms,_=A.verify_symmetries(old,Q,W,h)
    saved_orbits=json.loads(ORB.read_text());orbits=triples(perms);assert orbits==saved_orbits['orbits']
    central=json.loads(CENTRAL.read_text());check_central(central,Q,W,G,h,orbits)
    pencils=json.loads(PENCILS.read_text());check_pencils(pencils,Q,W,G,h,orbits)
    forged=deepcopy(pencils['half_pencils'][0]['constraints'][0]);forged['point_coordinate_polynomials'][0][-1]=str(F(forged['point_coordinate_polynomials'][0][-1])+1)
    row=pencils['half_pencils'][0];i,j=[x-1 for x in row['positive_walls_one_based']];k=row['negative_wall_one_based']-1
    try:check_constraint(forged,Q,W,G,h,project_data(Q,W,G,h,i,j,k))
    except AssertionError:mutation_rejected=True
    else:raise AssertionError('forged polynomial point was accepted')
    out={'status':'PASS','standard_library_only':True,'ambient_triple_count':1521520,'triple_orbits':49,
         'balanced_sections':147,'continuous_half_pencils':294,'all_continuous_half_pencils_excluded':True,
         'forged_parametric_point_rejected':mutation_rejected,'sturm_calibration_checks':sturm_checks,
         'scope':'Every hyperplane through the fixed h whose normal has three nonzero ambient-wall coefficients; two-wall boundary cases have the earlier separate certificate.'}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
