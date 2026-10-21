#!/usr/bin/env python3
"""Independent stdlib/Fraction replay of all 21,945 positive wall-pair pencils.

Imports no search module, numpy, sympy, graph or geometric package. The known
Borcherds chamber source is the input. Exact matrices, orbit coverage, actual
ridges and polynomial coefficient/gcd certificates are checked from scratch.
The auxiliary bounded algebraic sample is superseded by the universal proof
and is deliberately not a dependency of this checker.
"""
from __future__ import annotations
from collections import deque
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd,comb
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
CENTRAL=ROOT/'data/enumerations/fresh_h20_global_interior_sections_20261013.json'
PARAMETRIC=ROOT/'data/enumerations/fresh_h20_global_section_parameters_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_global_fraction_verify_20261013.json'

# Ascending coefficient order, normalized to a unique zero representation.
def trim(a):
    a=list(map(F,a))
    while len(a)>1 and a[-1]==0:a.pop()
    return tuple(a) if a else (F(0),)
def readp(a):return trim(reversed(a))
def add(a,b):return trim((a[k] if k<len(a) else 0)+(b[k] if k<len(b) else 0) for k in range(max(len(a),len(b))))
def neg(a):return trim(-x for x in a)
def sub(a,b):return add(a,neg(b))
def scale(a,c):return trim(c*x for x in a)
def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return trim(out)
def power(a,n):
    out=(F(1),)
    for _ in range(n):out=mul(out,a)
    return out
def divmodp(a,b):
    assert b!=(F(0),)
    r=list(a);q=[F(0)]*max(1,len(a)-len(b)+1)
    while len(r)>=len(b) and trim(r)!=(F(0),):
        degree=len(r)-len(b);c=r[-1]/b[-1];q[degree]=c
        for k,x in enumerate(b):r[k+degree]-=c*x
        r=list(trim(r))
    return trim(q),trim(r)
def gcdp(a,b):
    while b!=(F(0),):a,b=b,divmodp(a,b)[1]
    return scale(a,1/a[-1]) if a!=(F(0),) else a
Z=(F(0),);ONE=(F(1),);T=(F(0),F(1))
def transform_unit(a):
    n=len(a)-1;out=Z
    for k,c in enumerate(a):out=add(out,(F(0),)*k+tuple(c*comb(n-k,l) for l in range(n-k+1)))
    return out
def positive(a,unit=False):
    b=transform_unit(a) if unit else a
    assert b!=Z and all(x>=0 for x in b)
    if unit:assert sum(a)>0

def matmul(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def transpose(A):return list(map(list,zip(*A)))
def matvec(A,v):return [sum(a*b for a,b in zip(row,v)) for row in A]
def dot(x,y,Q):return sum(a*b for a,b in zip(x,matvec(Q,y)))
def pdot(x,y,Q):
    out=Z
    for i,row in enumerate(Q):
        for j,c in enumerate(row):
            if c:out=add(out,scale(mul(x[i],y[j]),c))
    return out

def primitive(v):
    d=gcd(*v);assert d>0
    return [x//d for x in v]
def source():
    raw=json.loads(SOURCE.read_text())['coxeter_input']
    Q=[[-int(x) for x in row] for row in raw['basis_gram']]
    W=[[int(x) for x in row] for row in raw['wall_coordinates_as_columns_transposed']]
    G=matmul(matmul(W,Q),transpose(W));h=primitive([sum(col) for col in zip(*W)])
    assert dot(h,h,Q)==-14 and set(dot(w,h,Q) for w in W)=={-1,-8}
    return Q,W,G,h

def verify_symmetries(data,Q,W,h):
    perms=[]
    for row in data['ambient_symmetry_generators']:
        A=row['matrix'];g=row['wall_permutation']
        assert sorted(g)==list(range(210))
        assert matmul(matmul(transpose(A),Q),A)==Q
        assert matmul(W,transpose(A))==[W[k] for k in g]
        assert matvec(A,h)==h
        perms.append(g)
    remaining=set(combinations(range(210),2));orbits=[]
    while remaining:
        first=min(remaining);seen={first};queue=[first]
        while queue:
            i,j=queue.pop()
            for g in perms:
                v=tuple(sorted((g[i],g[j])))
                if v not in seen:seen.add(v);queue.append(v)
        assert seen<=remaining;remaining-=seen;orbits.append((first,len(seen)))
    assert len(orbits)==data['pair_orbit_count']==9
    assert sum(n for pair,n in orbits)==data['pair_count']==21945
    assert len(data['sections'])==9
    return perms,orbits

def central_ridges(data,Q,W,G,h,orbits):
    assert data['center']==h and data['center_square']==-14
    cp=[dot(w,h,Q) for w in W]
    for row,((i,j),size) in zip(data['sections'],orbits):
        assert row['pair_one_based']==[i+1,j+1] and row['orbit_size']==size
        e=primitive([cp[j]*x-cp[i]*y for x,y in zip(W[i],W[j])]);q=dot(e,e,Q)
        assert e==row['section_normal'] and q==row['section_normal_square']>0 and dot(e,h,Q)==0
        normals=[];labels=[];lookup={};removed=0
        for l,w in enumerate(W):
            ew=dot(e,w,Q);v=[q*x-ew*y for x,y in zip(w,e)]
            assert dot(v,h,Q)<0
            if dot(v,v,Q)<=0:removed+=1;continue
            v=tuple(primitive(v))
            if v in lookup:labels[lookup[v]].append(l+1)
            else:lookup[v]=len(normals);normals.append(v);labels.append([l+1])
        assert len(normals)==row['positive_projected_normal_rays']
        assert row['discarded'].get('nonpositive_norm',0)==removed
        proof=row['actual_forbidden_ridge'];a,b=proof['normal_indices'];r,s=normals[a],normals[b]
        assert proof['normal_ambient_wall_labels_one_based']==[labels[a],labels[b]]
        aa,bb,ab=dot(r,r,Q),dot(s,s,Q),dot(r,s,Q)
        assert aa>0 and bb>0 and aa*bb>ab*ab
        assert proof['gram']==[[aa,ab],[ab,bb]] and F(proof['cosine_squared'])==F(ab*ab,aa*bb)
        assert ab>0 or F(ab*ab,aa*bb) not in {F(0),F(1,4),F(1,2),F(3,4)}
        p=proof['point'];vals=[dot(v,p,Q) for v in normals]
        assert dot(p,p,Q)==proof['point_square']<0 and dot(p,h,Q)==proof['center_pairing']<0
        assert dot(p,e,Q)==0 and [k for k,v in enumerate(vals) if v==0]==[a,b]
        assert max(v for v in vals if v!=0)==proof['maximum_other_wall_pairing']<0

def verify_elliptic(row,Q,W,G,h):
    i,j=[x-1 for x in row['pair_one_based']];a,b,c=G[i][i],G[j][j],G[i][j];delta=a*b-c*c
    assert delta>0;gcdpoly=Z
    for proof in row['triple_constraints']:
        k=proof['third_wall_one_based']-1;d,u,v=G[k][k],G[i][k],G[j][k]
        M=[[a,c,u],[c,b,v],[u,v,d]]
        assert M==proof['triple_gram']
        determinant=a*b*d+2*c*u*v-a*v*v-b*u*u-d*c*c
        assert a>0 and delta>0 and determinant>0
        point=proof['persistent_face_point'];vals=[dot(w,point,Q) for w in W]
        assert dot(point,point,Q)==proof['persistent_face_point_square']<0 and dot(point,h,Q)<0
        assert [l for l,x in enumerate(vals) if x==0]==sorted((i,j,k))
        assert all(x<0 for l,x in enumerate(vals) if l not in (i,j,k))
        lin=trim((a*v-c*u,b*u-c*v));assert all(x<=0 for x in lin) and lin!=Z
        N=mul(lin,lin);D=scale(sub(scale(trim((a,-2*c,b)),d),power(trim((u,-v)),2)),delta)
        assert N==readp(proof['cosine_squared_numerator']) and D==readp(proof['cosine_squared_denominator'])
        positive(D);bound=sub(scale(D,3),scale(N,4));positive(bound)
        assert bound==readp(proof['three_denominator_minus_four_numerator'])
        P=N
        for p in (sub(scale(N,4),D),sub(scale(N,2),D),sub(scale(N,4),scale(D,3)),add(sub(scale(mul(N,N),16),scale(mul(N,D),12)),mul(D,D))):P=mul(P,p)
        assert P==readp(proof['necessary_angle_polynomial']);gcdpoly=gcdp(gcdpoly,P)
    assert gcdpoly==readp(row['gcd'])
    assert gcdpoly[-1]==1 and all(x==0 for x in gcdpoly[:-1])
    assert row['gcd_positive_real_roots']==0 and row['status']=='ALL_POSITIVE_PARAMETERS_EXCLUDED'

def verify_nonelliptic(row,Q,W,G,h,perms):
    i,j=[x-1 for x in row['pair_one_based']];a,b=[x-1 for x in row['ridge_ambient_walls_one_based']]
    assert G[i][i]==G[j][j]==4 and G[i][j]<=-4
    unit=row['direct_domain']=='0<t<=1'
    assert row['direct_domain'] in ('0<t<=1','t>0')
    e=[trim((W[i][k],-W[j][k])) for k in range(22)]
    qe=trim((4,-2*G[i][j],4));assert pdot(e,e,Q)==qe
    center=[trim((h[k]+2*W[i][k],h[k]+2*W[j][k])) for k in range(22)]
    assert pdot(center,e,Q)==Z
    positive(neg(pdot(center,center,Q)),unit)
    v=[trim((G[l][i],-G[l][j])) for l in range(210)]
    def gram(x,y):return sub(scale(qe,G[x][y]),mul(v[x],v[y]))
    aa,bb,ab=gram(a,a),gram(b,b),gram(a,b)
    assert [aa,ab,bb]==[readp(p) for p in row['projected_gram_polynomials']]
    D=sub(mul(aa,bb),mul(ab,ab))
    for p in (qe,aa,bb,D,neg(ab)):positive(p,unit)
    point=[readp(p) for p in row['ridge_point_coordinate_polynomials']]
    square=pdot(point,point,Q);hp=pdot(point,[(F(x),) for x in h],Q)
    assert square==readp(row['ridge_point_square']) and hp==readp(row['ridge_point_center_pairing'])
    assert pdot(point,e,Q)==Z
    positive(neg(square),unit);positive(neg(hp),unit)
    zeros=[]
    for l,w in enumerate(W):
        val=pdot(point,[(F(x),) for x in w],Q)
        assert val==readp(row['all_ambient_wall_pairing_polynomials'][l])
        if val==Z:zeros.append(l)
        else:positive(neg(val),unit)
    assert zeros==[a,b]
    N=mul(ab,ab);den=mul(aa,bb)
    bounds=(N,sub(den,scale(N,4))) if unit else (sub(scale(N,4),den),sub(den,scale(N,2)))
    assert list(bounds)==[readp(p) for p in row['positive_angle_bound_polynomials']]
    for p in bounds:positive(p,unit)
    assert row['forbidden_open_interval']==('(0,1/4)' if unit else '(1/4,1/2)')
    if unit:
        pair=(i,j)
        for k in row['swapping_generator_word']:pair=(perms[k][pair[0]],perms[k][pair[1]])
        assert pair==(j,i)
    else:assert row['swapping_generator_word']==[]
    assert row['status']=='ALL_POSITIVE_PARAMETERS_EXCLUDED'

def verify_all(central,param):
    Q,W,G,h=source();perms,orbits=verify_symmetries(central,Q,W,h)
    central_ridges(central,Q,W,G,h,orbits)
    rows=param['all_parameter_exclusions'];others=param['nonelliptic_all_parameter_exclusions']
    assert len(rows)==6 and len(others)==3
    assert sorted(tuple(r['pair_one_based']) for r in rows+others)==sorted((i+1,j+1) for (i,j),n in orbits)
    for r in rows:verify_elliptic(r,Q,W,G,h)
    for r in others:verify_nonelliptic(r,Q,W,G,h,perms)
    return Q,W,G,h,perms

def main():
    central=json.loads(CENTRAL.read_text());param=json.loads(PARAMETRIC.read_text())
    Q,W,G,h,perms=verify_all(central,param)
    rejected=[]
    forged=deepcopy(central['sections'][0]['actual_forbidden_ridge']['point']);forged[0]+=1
    modified=deepcopy(central);modified['sections'][0]['actual_forbidden_ridge']['point']=forged
    try:central_ridges(modified,Q,W,G,h,verify_symmetries(central,Q,W,h)[1])
    except AssertionError:rejected.append('forged_central_ridge_point')
    bad=deepcopy(param['nonelliptic_all_parameter_exclusions'][0]);bad['ridge_point_coordinate_polynomials'][0][-1]=str(F(bad['ridge_point_coordinate_polynomials'][0][-1])+1)
    try:verify_nonelliptic(bad,Q,W,G,h,perms)
    except AssertionError:rejected.append('forged_parametric_ridge_point')
    assert len(rejected)==2
    result={'status':'PASS','standard_library_only':True,'pair_orbits':9,'all_wall_pairs':21945,
            'all_positive_parameter_families_excluded':9,'central_actual_ridges':9,
            'elliptic_persistent_triple_families':6,'nonelliptic_direct_ridge_families':3,
            'forged_certificates_rejected':rejected,
            'scope':'Every e=r_i-t*r_j, t>0, for two distinct ambient walls; no claim about hyperplanes whose normals involve three or more walls.'}
    OUTPUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
