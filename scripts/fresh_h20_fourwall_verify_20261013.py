#!/usr/bin/env python3
"""Independent standard-library replay of four-wall section certificates.

Rebuilds every normal, the two orthogonality constraints, the complete real
root census on every declared pencil, and actual ridge/volume certificates.
Numerical proposal selection is deliberately outside the validation scope.
"""
from collections import Counter,defaultdict
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd,isqrt,lcm
import fresh_h20_triple_verify_20261013 as V
import fresh_h20_triple_angle_verify_20261013 as U
import fresh_h20_bivariate_nonelliptic_verify_20261013 as B

A=V.A;ROOT=V.ROOT
SOURCE=ROOT/'data/enumerations/fresh_h20_fourwall_sections_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_fourwall_verify_20261013.json'


class Field(U.Field):
    def __init__(self,p,lo,hi):
        self.p=A.trim(p);self.lo=F(lo);self.hi=F(hi);self.rational=None
        assert self.p[-1]==1 and self.lo<self.hi
        assert V.evalp(self.p,self.lo)!=0 and V.evalp(self.p,self.hi)!=0
        assert V.rootcount(self.p,self.lo,self.hi)==1
        self.chain=V.sturm(self.p)


def rref(rows):
    rows=[list(map(F,r)) for r in rows];rank=0;pivots=[]
    if not rows:return rows,pivots
    for col in range(len(rows[0])):
        pivot=next((i for i in range(rank,len(rows)) if rows[i][col]),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank];q=rows[rank][col]
        rows[rank]=[x/q for x in rows[rank]]
        for i in range(len(rows)):
            if i!=rank:
                q=rows[i][col];rows[i]=[a-q*b for a,b in zip(rows[i],rows[rank])]
        pivots.append(col);rank+=1
        if rank==len(rows):break
    return rows,pivots


def positive_definite(M):
    M=[list(map(F,r)) for r in M]
    while M:
        if M[0][0]<=0:return False
        M=[[M[i][j]-M[i][0]*M[0][j]/M[0][0] for j in range(1,len(M))] for i in range(1,len(M))]
    return True


def no_rational_root(p):
    denominator=lcm(*(x.denominator for x in p));coefficients=[int(x*denominator) for x in p]
    divisor=gcd(*coefficients);coefficients=[x//divisor for x in coefficients]
    if coefficients[0]==0:return False
    def divisors(n):
        n=abs(n);out=set()
        for d in range(1,isqrt(n)+1):
            if n%d==0:out.update((d,n//d))
        return out
    return all(V.evalp(p,sign*F(a,b))!=0 for a in divisors(coefficients[0])
               for b in divisors(coefficients[-1]) for sign in (-1,1))


def plane_small_support(E,W):
    """Exact test: can any <=3 original wall vectors span the rational plane?"""
    reduced,pivots=rref(E);assert len(pivots)==2
    zeros=[];groups=defaultdict(list)
    for i,w in enumerate(W):
        v=list(map(F,w))
        for row,p in zip(reduced,pivots):
            c=v[p];v=[a-c*b for a,b in zip(v,row)]
        pivot=next((x for x in v if x),None)
        if pivot is None:zeros.append(i)
        else:groups[tuple(x/pivot for x in v)].append(i)
    cases=[zeros]+[ids+zeros for ids in groups.values()]
    for ids in cases:
        for size in (2,3):
            for js in combinations(ids,size):
                rows=[W[i] for i in js];rank=len(rref(rows)[1])
                if len(rref(rows+E)[1])==rank:return [i+1 for i in js]
    return None


def section(row,Q,W,G,h):
    ids=[i-1 for i in row['support_walls_one_based']]
    assert len(ids)==4 and ids==sorted(set(ids)) and all(0<=i<210 for i in ids)
    assert len(rref([W[i] for i in ids])[1])==4
    assert row['support_gram']==[[G[i][j] for j in ids] for i in ids]
    assert not positive_definite(row['support_gram'])
    assert len(rref(row['support_gram'])[1])==4
    basis=[list(map(F,v)) for v in row['normal_coefficient_basis']]
    assert len(basis)==2 and all(len(v)==4 for v in basis) and len(rref(basis)[1])==2
    assert all(basis[0][i] or basis[1][i] for i in range(4))
    assert any(c==0 for c in basis[1])  # t=infinity uses at most three listed walls.
    aux=[i-1 for i in row['auxiliary_orthogonal_walls_one_based']]
    assert len(aux)==2 and len(set(aux))==2 and not set(aux)&set(ids)
    assert len(rref([[G[l][i] for i in ids] for l in aux])[1])==2
    E=[[sum(v[k]*W[i][d] for k,i in enumerate(ids)) for d in range(22)] for v in basis]
    e=[A.trim((E[0][d],E[1][d])) for d in range(22)]
    q=A.trim((A.dot(E[0],E[0],Q),2*A.dot(E[0],E[1],Q),A.dot(E[1],E[1],Q)))
    vals=[A.trim((A.dot(w,E[0],Q),A.dot(w,E[1],Q))) for w in W]
    he=A.trim((A.dot(h,E[0],Q),A.dot(h,E[1],Q)))
    assert all(vals[l]==A.Z for l in aux) and he!=A.Z
    assert he==A.readp(row['center_pairing_polynomial'])
    assert sum(v==A.Z for v in vals)==row['identically_orthogonal_wall_count']
    def gr(a,b):return A.sub(A.scale(q,G[a][b]),A.mul(vals[a],vals[b]))
    a,b=[i-1 for i in row['target_pair_one_based']]
    N=A.mul(gr(a,b),gr(a,b));den=A.mul(gr(a,a),gr(b,b))
    equations={2:N,3:A.sub(A.scale(N,4),den),4:A.sub(A.scale(N,2),den),
               5:A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den)),
               6:A.sub(A.scale(N,4),A.scale(den,3))}
    assert row['target_labels']==[2,3,4,5,6]
    assert {int(m):A.readp(p) for m,p in row['target_equations'].items()}==equations
    return (e,q,vals,gr,N,den,equations),he,basis,aux,plane_small_support(E,W)


def root_census(row,sd):
    product=A.ONE
    for eq in sd[-1].values():assert eq!=A.Z;product=A.mul(product,eq)
    derivative=A.trim(k*product[k] for k in range(1,len(product)))
    squarefree=A.divmodp(product,V.gcdpoly(product,derivative))[0]
    maximum=max((abs(c/squarefree[-1]) for c in squarefree[:-1]),default=F(0))
    bound=F(maximum.numerator//maximum.denominator+2)
    expected=V.rootcount(squarefree,-bound,bound)
    intervals=[]
    for candidate in row['candidates']:
        field=Field(A.readp(candidate['minimal_polynomial']),*candidate['isolating_interval'])
        assert field.rem(product)==A.Z
        intervals.append((field.lo,field.hi))
    intervals.sort()
    assert all(l[1]<r[0] for l,r in zip(intervals,intervals[1:]))
    assert len(intervals)==expected
    return expected


def check(row,pencil,Q,W,G,h,sd,he,basis,aux):
    e,q,vals,gr,N,den,equations=sd
    field=Field(A.readp(row['minimal_polynomial']),*row['isolating_interval'])
    coefficients=[A.trim((basis[0][i],basis[1][i])) for i in range(4)]
    status=row['status']
    if status=='LOWER_SUPPORT_PARAMETER':
        assert any(field.sign(p)==0 for p in coefficients);return False
    assert all(field.sign(p)!=0 for p in coefficients)
    if status=='CENTRAL_PARAMETER':assert field.sign(he)==0;return False
    assert field.sign(he)!=0
    if status=='NO_SPACELIKE_NORMAL':assert field.sign(q)<=0;return False
    assert field.sign(q)>0
    a,b=[i-1 for i in pencil['target_pair_one_based']]
    if status=='TARGET_NOT_ACUTE_ELLIPTIC_PAIR':
        assert min(field.sign(gr(a,a)),field.sign(gr(b,b)))<=0 or field.sign(gr(a,b))>0;return False
    assert min(field.sign(gr(a,a)),field.sign(gr(b,b)))>0 and field.sign(gr(a,b))<=0
    labels=[m for m,p in equations.items() if field.rem(p)==A.Z and (m!=5 or field.sign(A.sub(A.scale(N,2),den))>0)]
    if status=='CONJUGATE_ANGLE_BRANCH':assert not labels;return False
    assert labels==row['target_labels'] and labels
    C=[field.rem(A.sub(A.scale(q,h[d]),A.mul(he,e[d]))) for d in range(22)]
    dual,cp,square=U.field_products(C,Q,W,field)
    if status=='UNRESOLVED_STRICT_CENTER_SCOUT':assert any(field.sign(p)>=0 for p in cp);return False
    if 'strict_interior' in row:
        if row['interior_method']=='projected_original_center':assert C==[A.readp(p) for p in row['strict_interior']]
        else:assert row['interior_method']=='rational_LP_point_projected_exactly'
        C=[A.readp(p) for p in row['strict_interior']];dual,cp,square=U.field_products(C,Q,W,field)
    assert field.sign(square)<0 and all(field.sign(p)<0 for p in cp)
    assert field.rem(V.sum_polys(A.mul(dual[d],e[d]) for d in range(22)))==A.Z
    assert field.sign(V.sum_polys(A.scale(dual[d],h[d]) for d in range(22)))<0
    if status in ('UNRESOLVED_TARGET_RIDGE_SCOUT','UNRESOLVED_RIGHT_RIDGE_SCOUT'):return False
    assert C==[A.readp(p) for p in row['strict_interior']]
    repaired=row['repaired_target_ridge'] is not None and len(row['repaired_right_ridges'])==2
    if row['repaired_target_ridge'] is not None:
        assert row['repaired_target_ridge']['ridge_walls_one_based']==pencil['target_pair_one_based']
        B.actual(row['repaired_target_ridge'],Q,W,G,h,sd,field)
    rights=row['repaired_right_ridges'];assert len(rights)<=2
    assert len({tuple(sorted(r['ridge_walls_one_based'])) for r in rights})==len(rights)
    if repaired:
        assert len({tuple(proof['zero_walls_one_based']) for proof in [row['repaired_target_ridge'],*rights]})==3
    for l,proof in zip(aux,rights):
        assert l+1 in proof['ridge_walls_one_based']
        aa,ab,bb=B.actual(proof,Q,W,G,h,sd,field);assert ab==A.Z
    assert status==('NON_COXETER_ACTUAL_RIDGE' if repaired else 'NON_COXETER_WITH_INCOMPLETE_REPAIR_CERTIFICATES')
    proof=row['actual_forbidden_ridge'];aa,ab,bb=B.actual(proof,Q,W,G,h,sd,field)
    N=field.rem(A.mul(ab,ab));den=field.rem(A.mul(aa,bb))
    if proof['angle_rejection']=='obtuse':assert field.sign(ab)>0
    else:
        assert field.sign(ab)<=0 and field.sign(N)>0 and field.sign(A.sub(A.scale(den,3),A.scale(N,4)))>0
        for p in (A.sub(A.scale(N,4),den),A.sub(A.scale(N,2),den),
                  A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den))):assert field.sign(p)!=0
    return True


def main():
    Q,W,G,h=A.source();data=json.loads(SOURCE.read_text());assert data['complete']
    assert len(data['pencils'])==data['pencil_count']==data['scouting']['selected']
    statuses=Counter();accepted=total=0;degrees=Counter();cases=[];plane_reductions=[];genuine_four=0
    for pencil in data['pencils']:
        sd,he,basis,aux,reduction=section(pencil,Q,W,G,h);total+=root_census(pencil,sd)
        assert reduction is None
        plane_reductions.append({'support':pencil['support_walls_one_based'],'three_wall_spanning_set':reduction})
        for row in pencil['candidates']:
            yes=check(row,pencil,Q,W,G,h,sd,he,basis,aux);accepted+=yes;statuses[row['status']]+=1
            if yes:
                degrees[len(row['minimal_polynomial'])-1]+=1;cases.append((row,pencil,sd,he,basis,aux))
                if len(row['minimal_polynomial'])>2:
                    assert no_rational_root(A.readp(row['minimal_polynomial']));genuine_four+=1
    assert dict(statuses)==data['statuses'] and cases
    row,pencil,sd,he,basis,aux=cases[0];forged=deepcopy(row)
    forged['repaired_target_ridge']['point_coordinate_polynomials'][0][-1]=str(F(forged['repaired_target_ridge']['point_coordinate_polynomials'][0][-1])+1)
    try:check(forged,pencil,Q,W,G,h,sd,he,basis,aux)
    except AssertionError:mutation=True
    else:raise AssertionError('forged repaired ridge point accepted')
    forged=deepcopy(pencil);forged['normal_coefficient_basis'][0][0]=str(F(forged['normal_coefficient_basis'][0][0])+1)
    try:section(forged,Q,W,G,h)
    except AssertionError:normal_mutation=True
    else:raise AssertionError('forged auxiliary orthogonality accepted')
    out={'status':'PASS','standard_library_only':True,'pencils':len(data['pencils']),
         'complete_real_angle_equation_root_census':total,'statuses':dict(statuses),
         'finite_volume_sections_with_three_actual_angle_repairs':statuses['NON_COXETER_ACTUAL_RIDGE'],
         'finite_volume_sections':accepted,'actual_forbidden_ridge_rejections':accepted,
         'accepted_parameter_degrees':dict(degrees),'normal_plane_three_wall_reductions':plane_reductions,
         'irrational_candidates_proved_not_in_any_three_original_wall_span':genuine_four,
         'forged_repaired_ridge_point_rejected':mutation,'forged_normal_coefficient_rejected':normal_mutation,
         'scope':'Exact validation of the literal declared pencils, root censuses and archived geometry. No replay/completeness claim for numerical pencil selection. Normal-plane test excludes all <=3-wall representations for irrational parameters only when no spanning set is found; rational specializations require separate checking.'}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':main()
