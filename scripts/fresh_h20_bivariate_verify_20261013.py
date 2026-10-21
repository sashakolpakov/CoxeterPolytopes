#!/usr/bin/env python3
"""Independent stdlib/Fraction replay of all elliptic three-wall parameters.

No producer, sympy, numpy, numerical root solver, or project geometry import.
Sparse bivariate polynomials reconstruct every ridge identity.  The isolated
quartic candidate uses the separately implemented Fraction field verifier.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import comb
import json
import fresh_h20_triple_verify_20261013 as V
import fresh_h20_triple_angle_verify_20261013 as U

A=V.A;ROOT=V.ROOT
SOURCE=ROOT/'data/enumerations/fresh_h20_bivariate_elliptic_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_bivariate_verify_20261013.json'
ZERO={};ONE={(0,0):F(1)};S={(1,0):F(1)};T={(0,1):F(1)}

def clean(p):return {k:F(v) for k,v in p.items() if v}
def const(c):return clean({(0,0):c})
def scale(p,c):return clean({k:v*c for k,v in p.items()})
def add(p,q):
    r=dict(p)
    for k,c in q.items():r[k]=r.get(k,F(0))+c
    return clean(r)
def sub(p,q):return add(p,scale(q,-1))
def mul(p,q):
    r={}
    for (a,b),c in p.items():
        for (d,e),f in q.items():r[a+d,b+e]=r.get((a+d,b+e),F(0))+c*f
    return clean(r)
def power(p,n):
    out=ONE
    for _ in range(n):out=mul(out,p)
    return out
def summ(ps):
    out=ZERO
    for p in ps:out=add(out,p)
    return out
def read(p):
    assert len({(a,b) for a,b,c in p})==len(p)
    return clean({(int(a),int(b)):F(c) for a,b,c in p})
def positive(p):assert p and all(c>=0 for c in p.values())
def pdot(x,y,Q):
    return summ(scale(mul(x[i],y[j]),Q[i][j]) for i in range(22) for j in range(22) if Q[i][j])
def cv(x):return [const(c) for c in x]

def det3(M):
    a,b,c=M[0];d,e,f=M[1];g,h,i=M[2]
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
def elliptic(G,ids):
    M=[[G[i][j] for j in ids] for i in ids]
    return M[0][0]>0 and M[0][0]*M[1][1]>M[0][1]**2 and det3(M)>0

def section(Q,W,G,i,j,k):
    e=[sub(add(const(W[i][d]),scale(S,W[j][d])),scale(T,W[k][d])) for d in range(22)]
    q=pdot(e,e,Q)
    v=[pdot(e,cv(w),Q) for w in W]
    def gram(a,b):return sub(scale(q,G[a][b]),mul(v[a],v[b]))
    return e,q,v,gram

def check_face(z,ids,Q,W,h):
    assert len(z)==22 and A.dot(z,z,Q)<0 and A.dot(z,h,Q)<0
    vals=[A.dot(w,z,Q) for w in W]
    assert [l for l,p in enumerate(vals) if p==0]==sorted(ids)
    assert all(p<=0 for p in vals)

def check_ridge(proof,Q,W,G,h,i,j,k,sd):
    e,q,v,gram=sd;a,l=[x-1 for x in proof['ridge_walls_one_based']]
    assert a in (i,j) and l not in (i,j,k)
    other=j if a==i else i
    z1,z2=proof['ambient_face_points']
    check_face(z1,(a,other,l),Q,W,h);check_face(z2,(a,k,l),Q,W,h)
    assert A.dot(z1,z2,Q)<0
    u=pdot(e,cv(z1),Q);v2=pdot(e,cv(z2),Q)
    assert [u,v2]==[read(p) for p in proof['face_normal_pairings']]
    positive(u);positive(scale(v2,-1))
    p=[sub(scale(u,z2[d]),scale(v2,z1[d])) for d in range(22)]
    assert p==[read(x) for x in proof['point_coordinate_polynomials']]
    assert pdot(p,e,Q)==ZERO
    square=pdot(p,p,Q);assert square==read(proof['point_square']);positive(scale(square,-1))
    positive(scale(pdot(p,cv(h),Q),-1))
    for d,w in enumerate(W):
        value=pdot(p,cv(w),Q)
        if d in (a,l):assert value==ZERO
        else:positive(scale(value,-1))
    aa,ab,bb=gram(a,a),gram(a,l),gram(l,l)
    assert [aa,ab,bb]==[read(x) for x in proof['projected_gram']]
    # All three e,ra,rl lie in the positive-definite p^perp. They are
    # independent: ra,rl are independent (elliptic ambient face), and e
    # pairs strictly positively with z1 whereas both ra,rl vanish there.
    # Thus their projected pair is positive definite and a true ridge.
    assert G[a][a]*G[l][l]>G[a][l]**2
    return a,l,aa,ab,bb

def line(p):
    """Substitute s=u,t=1-2u, returning low-first univariate coefficients."""
    out=A.Z
    for (a,b),c in p.items():
        term=A.mul(A.power(A.T,a),A.power(A.trim((1,-2)),b))
        out=A.add(out,A.scale(term,c))
    return out

def positive_half(p):
    """Exact positivity on 0<u<1/2 via u=x/[2(1+x)]."""
    p=A.trim(p);d=len(p)-1;out=[F(0)]*(d+1)
    for a,c in enumerate(p):
        for b in range(d-a+1):out[a+b]+=c*F(comb(d-a,b),2**a)
    assert any(out) and all(c>=0 for c in out)

def check_family(row,Q,W,G,h):
    i,j=[x-1 for x in row['positive_walls_one_based']];k=row['negative_wall_one_based']-1
    assert i<j and sorted((i+1,j+1,k+1))==row['triple_one_based']
    assert elliptic(G,(i,j,k));check_face(row['triple_face_point'],(i,j,k),Q,W,h)
    sd=section(Q,W,G,i,j,k);e,q,v,gram=sd
    ridges=[check_ridge(p,Q,W,G,h,i,j,k,sd) for p in row['ridges']]
    if row['kind']=='one_global_forbidden_ridge':
        assert len(ridges)==1
        a,l,aa,ab,bb=ridges[0];N=mul(ab,ab);den=mul(aa,bb)
        reason=row['ridges'][0]['reason']
        if reason=='obtuse_on_positive_quadrant':positive(ab)
        elif reason=='0<cosine_squared<1/4':
            positive(scale(ab,-1));positive(sub(den,scale(N,4)))
        elif reason=='1/4<cosine_squared<1/2':
            positive(scale(ab,-1));positive(sub(scale(N,4),den));positive(sub(den,scale(N,2)))
        else:raise AssertionError(reason)
    elif row['kind']=='two_positive_covector_branches':
        assert len(ridges)==2
        for (a,b),r in zip(((i,j),(j,i)),ridges):
            aa,l,_,ab,_=r;assert aa==a and G[l][a]==G[l][k]==0 and G[l][b]<0
            positive(scale(v[l],-1));assert ab==scale(mul(v[a],v[l]),-1)
        # At least one of v_i,v_j is positive. The displayed sum is a
        # strictly positive completed square plus nonnegative t-terms.
        Gii,Gij,Gjj=G[i][i],G[i][j],G[j][j]
        assert Gii*Gjj-Gij*Gij>0 and G[i][k]<=0 and G[j][k]<=0
        completed=add(const(F(Gii*Gjj-Gij*Gij,Gjj)),scale(power(add(S,const(F(Gij,Gjj))),2),Gjj))
        completed=sub(sub(completed,scale(T,G[i][k])),scale(mul(S,T),G[j][k]))
        assert add(v[i],mul(S,v[j]))==completed
    elif row['kind']=='covector_sign_and_zero_line':
        assert (i,j,k)==(0,47,19) and len(ridges)==3
        B=sub(add(scale(S,2),T),ONE);assert read(row['branch_function'])==B
        assert ridges[0][0:2]==(47,21) and ridges[0][3]==scale(B,2)
        assert ridges[1][0:2]==(47,18) and ridges[1][3]==scale(mul(T,B),-4)
        # B>0 and B<0 each give an obtuse actual ridge. B=0 gives
        # t=1-2s, with exactly 0<s<1/2, handled without endpoint omissions.
        a,l,aa,ab,bb=ridges[2];assert (a,l)==(0,18)
        aa,ab,bb=map(line,(aa,ab,bb));N=A.mul(ab,ab);den=A.mul(aa,bb)
        positive_half(A.neg(ab));positive_half(A.sub(den,A.scale(N,4)))
        wantedN=A.power(A.trim((-1,2)),2);wantedD=A.trim((4,-12,10))
        assert A.mul(N,wantedD)==A.mul(den,wantedN)
    else:raise AssertionError(row['kind'])
    assert row['status']=='ALL_POSITIVE_S_T_EXCLUDED'
    return sd

def check_candidate(row,Q,W,G,h,exceptional):
    field=U.Field(A.readp(row['minimal_polynomial']),*row['isolating_interval'])
    i,j=0,47;k=19
    S1=A.readp(row['s_polynomial']);T1=A.readp(row['t_polynomial'])
    assert S1==A.T and T1==A.trim((1,-2))
    assert field.sign(S1)>0 and field.sign(T1)>0
    assert field.p==A.trim((F(2,23),F(-16,23),F(44,23),F(-52,23),1))
    # The primitive quartic reduces to 2*(1+x+x^2+x^3+x^4) modulo3.
    # No monic linear or quadratic divides it, proving irreducibility.
    for degree in (1,2):
        for code in range(3**degree):
            divisor=[]
            for _ in range(degree):divisor.append(code%3);code//=3
            divisor.append(1);remainder=[1]*5
            while len(remainder)>=len(divisor):
                c=remainder[-1];shift=len(remainder)-len(divisor)
                for d,v in enumerate(divisor):remainder[shift+d]=(remainder[shift+d]-c*v)%3
                while remainder and remainder[-1]==0:remainder.pop()
            assert remainder
    e=[A.sub(A.add(A.trim((W[i][d],)),A.scale(S1,W[j][d])),A.scale(T1,W[k][d])) for d in range(22)]
    assert e==[A.readp(p) for p in row['section_normal_coordinate_polynomials']]
    dual,vals,q=U.field_products(e,Q,W,field);assert field.sign(q)>0
    hp=field.rem(V.sum_polys(A.scale(dual[d],h[d]) for d in range(22)))
    assert hp==A.readp(row['section_normal_center_pairing']) and field.sign(hp)!=0
    def gr(a,b):return field.rem(A.sub(A.scale(q,G[a][b]),A.mul(vals[a],vals[b])))
    N=field.rem(A.mul(gr(i,j),gr(i,j)));den=field.rem(A.mul(gr(i,i),gr(j,j)))
    equation=field.rem(A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),16)),A.scale(A.mul(den,den),2)))
    assert equation==A.Z and field.sign(gr(i,j))<0 and field.sign(A.sub(A.scale(N,4),A.scale(den,3)))>0
    repaired=row['repaired_ridges'];assert [r['label'] for r in repaired]==[8,2]
    check_face(repaired[0]['point'],(i,j,k),Q,W,h)
    # The third zero inequality is a nonnegative combination, because
    # t*r'_k=r'_i+s*r'_j and s,t are both strictly positive.
    projected=lambda a:[field.rem(A.sub(A.scale(q,W[a][d]),A.mul(vals[a],e[d]))) for d in range(22)]
    ri,rj,rk=map(projected,(i,j,k))
    assert all(field.rem(A.sub(A.sub(A.mul(T1,rk[d]),ri[d]),A.mul(S1,rj[d])))==A.Z for d in range(22))
    assert gr(j,21)==A.Z
    point=[field.rem(line(read(p))) for p in exceptional['ridges'][0]['point_coordinate_polynomials']]
    _,pairings,square=U.field_products(point,Q,W,field)
    assert field.sign(square)<0 and [d for d,p in enumerate(pairings) if field.sign(p)==0]==[21,47]
    assert all(field.sign(p)<=0 for p in pairings)
    # Reuse only the independent Fraction geometry checker, supplying this
    # independently reconstructed section and exact pi/8 equation.
    converted={**row,'repaired_ridge_labels':[8], 'interior_point':row['strict_interior'],
               'status':'NON_COXETER_ACTUAL_RIDGE'}
    U.verify_candidate(converted,Q,W,G,h,(e,q,vals,gr,N,den,{8:equation}))
    assert row['status']=='TWO_ACTUAL_ANGLES_REPAIRED_BUT_NON_COXETER'

def main():
    Q,W,G,h=A.source();prior=json.loads(A.CENTRAL.read_text())
    perms,_=A.verify_symmetries(prior,Q,W,h)
    orbits=V.triples(perms);assert orbits==json.loads(V.ORB.read_text())['orbits']
    ell=[r for r in orbits if elliptic(G,[x-1 for x in r['representative']])]
    data=json.loads(SOURCE.read_text());assert data['complete'] and data['elliptic_triple_orbits']==ell
    assert len(ell)==15 and data['elliptic_triple_count']==sum(r['size'] for r in ell)==266000
    rows=data['families'];expected={(tuple(r['representative']),k) for r in ell for k in r['representative']}
    assert len(rows)==45 and {(tuple(r['triple_one_based']),r['negative_wall_one_based']) for r in rows}==expected
    counts={}
    for row in rows:
        check_family(row,Q,W,G,h);counts[row['kind']]=counts.get(row['kind'],0)+1
    exceptional=next(r for r in rows if r['kind']=='covector_sign_and_zero_line')
    check_candidate(data['two_angle_candidate'],Q,W,G,h,exceptional)
    forged=deepcopy(rows[0]);forged['ridges'][0]['ambient_face_points'][0][0]+=1
    try:check_family(forged,Q,W,G,h)
    except AssertionError:face_mutation=True
    else:raise AssertionError('forged ambient face point accepted')
    forged=deepcopy(data['two_angle_candidate']);forged['s_polynomial'][-1]='1'
    try:check_candidate(forged,Q,W,G,h,exceptional)
    except AssertionError:candidate_mutation=True
    else:raise AssertionError('forged coefficient parameter accepted')
    out={'status':'PASS','standard_library_only':True,'elliptic_triple_orbits':15,
         'elliptic_triple_count':266000,'independent_positive_parameter_families':45,
         'proof_kinds':counts,'actual_ridge_certificates':sum(len(r['ridges']) for r in rows),
         'unequal_coefficient_two_angle_quartic_candidate':'PASS; exact actual forbidden ridge',
         'forged_face_point_rejected':face_mutation,'forged_parameter_rejected':candidate_mutation,
         'scope':data['scope']}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
