#!/usr/bin/env python3
"""Independent Fraction replay of bounded nonelliptic two-angle proposals.

Reconstructs the declared lines and angle equations, counts their positive
roots, rejects the wrong angle conjugates, and checks all actual ridge points
of accepted algebraic candidates. No producer, SymPy or NumPy is imported.
"""
from copy import deepcopy
from fractions import Fraction as F
import json
import fresh_h20_triple_verify_20261013 as V
import fresh_h20_triple_angle_verify_20261013 as U
from fresh_h20_bivariate_verify_20261013 import elliptic

A=V.A;ROOT=V.ROOT
SOURCE=ROOT/'data/enumerations/fresh_h20_bivariate_nonelliptic_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_bivariate_nonelliptic_verify_20261013.json'

def section(row,Q,W,G,h):
    i,j=[x-1 for x in row['positive_walls_one_based']];k=row['negative_wall_one_based']-1
    l=row['auxiliary_orthogonal_wall_one_based']-1
    assert i<j and sorted((i+1,j+1,k+1))==row['triple_one_based']
    assert not elliptic(G,(i,j,k))
    assert row['triple_gram']==[[G[x][y] for y in (i,j,k)] for x in (i,j,k)]
    a,b=map(F,row['line_t_coefficients']);assert a>=0 and b>=0 and a+b>0
    assert G[l][k]<0 and a==F(G[l][i],G[l][k]) and b==F(G[l][j],G[l][k])
    e=[A.trim((W[i][d]-a*W[k][d],W[j][d]-b*W[k][d])) for d in range(22)]
    dual=[[A.scale(e[jj],Q[ii][jj]) for jj in range(22) if Q[ii][jj]] for ii in range(22)]
    dual=[V.sum_polys(p) for p in dual]
    q=V.sum_polys(A.mul(e[d],dual[d]) for d in range(22))
    vals=[V.sum_polys(A.scale(dual[d],w[d]) for d in range(22)) for w in W]
    he=V.sum_polys(A.scale(dual[d],h[d]) for d in range(22))
    assert vals[l]==A.Z and he!=A.Z and he==A.readp(row['section_center_pairing'])
    def gr(a,b):return A.sub(A.scale(q,G[a][b]),A.mul(vals[a],vals[b]))
    aa,bb=[x-1 for x in row['target_pair_one_based']]
    N=A.mul(gr(aa,bb),gr(aa,bb));den=A.mul(gr(aa,aa),gr(bb,bb))
    equations={3:A.sub(A.scale(N,4),den),4:A.sub(A.scale(N,2),den),
               5:A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den)),
               6:A.sub(A.scale(N,4),A.scale(den,3))}
    assert row['target_labels']==[3,4,5,6]
    assert {int(m):A.readp(p) for m,p in row['target_equations'].items()}==equations
    return (e,q,vals,gr,N,den,equations),he,l

def actual(proof,Q,W,G,h,sd,field):
    e,q,vals,gr,*_=sd;a,b=[x-1 for x in proof['ridge_walls_one_based']]
    aa,ab,bb=map(field.rem,(gr(a,a),gr(a,b),gr(b,b)));D=field.rem(A.sub(A.mul(aa,bb),A.mul(ab,ab)))
    assert [aa,ab,bb]==[A.readp(p) for p in proof['projected_gram']]
    assert min(field.sign(aa),field.sign(bb),field.sign(D))>0
    p=[A.readp(x) for x in proof['point_coordinate_polynomials']]
    dual,pairings,square=U.field_products(p,Q,W,field)
    assert square==A.readp(proof['point_square']) and field.sign(square)<0
    hp=field.rem(V.sum_polys(A.scale(dual[d],h[d]) for d in range(22)))
    assert hp==A.readp(proof['point_h_pairing']) and field.sign(hp)<0
    assert field.rem(V.sum_polys(A.mul(dual[d],e[d]) for d in range(22)))==A.Z
    signs=[field.sign(x) for x in pairings];zeros=[d for d,v in enumerate(signs) if v==0]
    assert [x+1 for x in zeros]==proof['zero_walls_one_based'] and a in zeros and b in zeros
    assert max(signs)<=0
    cones=proof['redundant_zero_wall_cones']
    assert sorted(r['wall_one_based']-1 for r in cones)==sorted(set(zeros)-{a,b})
    for row in cones:
        l=row['wall_one_based']-1
        al=field.rem(A.sub(A.mul(bb,gr(l,a)),A.mul(ab,gr(l,b))))
        be=field.rem(A.sub(A.mul(aa,gr(l,b)),A.mul(ab,gr(l,a))))
        assert al==A.readp(row['alpha_numerator']) and be==A.readp(row['beta_numerator']) and D==A.readp(row['denominator'])
        assert min(field.sign(al),field.sign(be))>=0
        residual=field.rem(A.sub(A.sub(A.mul(D,gr(l,l)),A.mul(al,gr(l,a))),A.mul(be,gr(l,b))))
        assert residual==A.Z
        # The residual normal lies in p^perp and has zero square; p is
        # timelike, so positive definiteness makes the residual vector zero.
    return aa,ab,bb

def census(row,sd):
    equations=sd[-1];product=A.ONE
    for p in equations.values():product=A.mul(product,p)
    squarefree=U.squarefree_positive(product)
    maximum=max(abs(c/squarefree[-1]) for c in squarefree[:-1]) if len(squarefree)>1 else F(0)
    bound=F(maximum.numerator//maximum.denominator+2)
    expected=V.rootcount(squarefree,F(0),bound) if len(squarefree)>1 else 0
    intervals=[]
    for candidate in row['candidates']:
        field=U.Field(A.readp(candidate['minimal_polynomial']),*candidate['isolating_interval'])
        assert field.rem(product)==A.Z
        intervals.append((field.lo,field.hi))
    intervals.sort()
    assert all(left[1]<right[0] for left,right in zip(intervals,intervals[1:]))
    assert len(intervals)==expected
    return expected

def check_candidate(row,pencil,Q,W,G,h,sd,he,orthogonal_wall):
    e,q,vals,gr,N,den,equations=sd
    field=U.Field(A.readp(row['minimal_polynomial']),*row['isolating_interval'])
    aa,bb=[x-1 for x in pencil['target_pair_one_based']]
    assert field.sign(q)>0 and min(field.sign(gr(aa,aa)),field.sign(gr(bb,bb)))>0
    assert field.sign(gr(aa,bb))<=0
    labels=[m for m,p in equations.items() if field.rem(p)==A.Z and (m!=5 or field.sign(A.sub(A.scale(N,2),den))>0)]
    if row['status']=='CONJUGATE_ANGLE_BRANCH':
        assert not labels and field.rem(equations[5])==A.Z and field.sign(A.sub(A.scale(N,2),den))<0
        return False
    assert row['status']=='NON_COXETER_ACTUAL_RIDGE' and row['target_labels']==labels and labels
    assert field.sign(he)!=0
    converted={**row,'repaired_ridge_labels':labels,'interior_point':row['strict_interior']}
    U.verify_candidate(converted,Q,W,G,h,sd)
    assert row['repaired_target_ridge']['ridge_walls_one_based']==pencil['target_pair_one_based']
    actual(row['repaired_target_ridge'],Q,W,G,h,sd,field)
    a,b=[x-1 for x in row['repaired_right_ridge']['ridge_walls_one_based']]
    assert orthogonal_wall in (a,b)
    aa,ab,bb=actual(row['repaired_right_ridge'],Q,W,G,h,sd,field)
    assert ab==A.Z
    return True

def main():
    Q,W,G,h=A.source();data=json.loads(SOURCE.read_text());assert data['complete']
    assert data['pencil_count']==len(data['pencils'])==data['selection_limit']==12
    total=accepted=conjugates=0;case=None
    for row in data['pencils']:
        sd,he,l=section(row,Q,W,G,h);total+=census(row,sd)
        for candidate in row['candidates']:
            yes=check_candidate(candidate,row,Q,W,G,h,sd,he,l)
            accepted+=yes;conjugates+=not yes
            if yes:case=(candidate,row,sd,he,l)
    assert total==12 and accepted==2 and conjugates==10 and case is not None
    candidate,row,sd,he,l=case;forged=deepcopy(candidate)
    forged['repaired_target_ridge']['point_coordinate_polynomials'][0][-1]=str(F(forged['repaired_target_ridge']['point_coordinate_polynomials'][0][-1])+1)
    try:check_candidate(forged,row,Q,W,G,h,sd,he,l)
    except AssertionError:mutation=True
    else:raise AssertionError('forged actual target ridge accepted')
    out={'status':'PASS','standard_library_only':True,'declared_noncentral_nonelliptic_pencils':12,
         'complete_positive_angle_equation_root_census':total,'wrong_conjugate_angle_roots':conjugates,
         'two_actual_angles_repaired':accepted,'exact_actual_forbidden_ridge_rejections':accepted,
         'forged_repaired_ridge_point_rejected':mutation,
         'scope':'Verifies all declared pencil root denominators and accepted algebraic geometry. Does not replay or claim completeness of the numerical/deterministic pencil-selection scout; no global nonelliptic exclusion.'}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
