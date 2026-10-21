#!/usr/bin/env python3
"""Independent Fraction algebraic-field replay of 64 noncentral angle repairs.

Checks the complete positive parameter census for the declared 17 pencils
and labels2,3,4,5,6, every strict timelike interior, the repaired local angle,
and every other actual forbidden ridge. No search module or sympy is used.
"""
from __future__ import annotations
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
import json
import fresh_h20_triple_verify_20261013 as V
A=V.A
add,sub,mul,scale,power,readp=V.add,V.sub,V.mul,V.scale,V.power,V.readp
Z=V.Z;ONE=V.ONE;T=V.T
OUTPUT=V.ROOT/'data/enumerations/fresh_h20_triple_angle_verify_20261013.json'

class Field:
    def __init__(self,p,lo,hi):
        self.p=A.trim(p);self.lo=F(lo);self.hi=F(hi);self.rational=None
        assert self.p[-1]==1 and self.lo>0 and self.lo<self.hi
        assert V.rootcount(self.p,self.lo,self.hi)==1
        self.chain=V.sturm(self.p)
    def rem(self,p):return A.divmodp(p,self.p)[1]
    def interval(self,p):
        lo=hi=F(0)
        for c in reversed(p):
            values=(lo*self.lo,lo*self.hi,hi*self.lo,hi*self.hi)
            lo=min(values)+c;hi=max(values)+c
        return lo,hi
    def refine(self):
        mid=(self.lo+self.hi)/2
        if V.evalp(self.p,mid)==0:self.rational=mid;return
        if V.variation(self.chain,self.lo)-V.variation(self.chain,mid)==1:self.hi=mid
        else:self.lo=mid
    def sign(self,p):
        p=self.rem(p)
        if p==Z:return 0
        for _ in range(4000):
            if self.rational is not None:return V.sign(V.evalp(p,self.rational))
            lo,hi=self.interval(p)
            if lo>0:return 1
            if hi<0:return -1
            self.refine()
        raise AssertionError('algebraic sign did not settle')


def squarefree_positive(p):
    p=A.trim(p)
    while len(p)>1 and p[0]==0:p=p[1:]
    assert p!=Z
    derivative=A.trim(k*p[k] for k in range(1,len(p)))
    if derivative!=Z:p=A.divmodp(p,V.gcdpoly(p,derivative))[0]
    return p

def positive_roots(p):
    p=squarefree_positive(p)
    if len(p)==1:return []
    ratio=max(abs(F(c,p[-1])) for c in p[:-1]);bound=F(ratio.numerator//ratio.denominator+2)
    assert V.evalp(p,bound)!=0
    count=V.rootcount(p,F(0),bound);queue=[(F(0),bound,count)];out=[]
    while queue:
        lo,hi,count=queue.pop()
        if count==0:continue
        if count==1:
            # Field requires a strictly positive lower endpoint.
            while lo==0:
                mid=(lo+hi)/2
                if V.evalp(p,mid)==0:mid=(lo+2*hi)/3
                if V.evalp(p,mid)==0:mid=(lo+3*hi)/4
                assert V.evalp(p,mid)!=0
                if V.rootcount(p,mid,hi)==1:lo=mid
                else:hi=mid
            out.append((scale(p,1/p[-1]),lo,hi));continue
        mid=(lo+hi)/2
        for denominator in range(3,20):
            if V.evalp(p,mid)!=0:break
            mid=(lo+(denominator-1)*hi)/denominator
        assert V.evalp(p,mid)!=0
        left=V.rootcount(p,lo,mid)
        queue.extend(((lo,mid,left),(mid,hi,count-left)))
    return out


def section(Q,W,G,i,j,k):
    e0=[x+y for x,y in zip(W[i],W[j])];e1=[-x for x in W[k]]
    e=[A.trim((x,y)) for x,y in zip(e0,e1)]
    qe=A.trim((A.dot(e0,e0,Q),2*A.dot(e0,e1,Q),A.dot(e1,e1,Q)))
    v=[A.trim((A.dot(w,e0,Q),A.dot(w,e1,Q))) for w in W]
    def gram(a,b):return sub(scale(qe,G[a][b]),mul(v[a],v[b]))
    N=mul(gram(i,j),gram(i,j));den=mul(gram(i,i),gram(j,j))
    equations={2:N,3:sub(scale(N,4),den),4:sub(scale(N,2),den),5:add(sub(scale(mul(N,N),16),scale(mul(N,den),12)),mul(den,den)),6:sub(scale(N,4),scale(den,3))}
    return e,qe,v,gram,N,den,equations


def field_products(point,Q,W,field):
    dual=[field.rem(V.sum_polys([scale(point[j],Q[i][j]) for j in range(22) if Q[i][j]])) for i in range(22)]
    vals=[field.rem(V.sum_polys([scale(dual[j],w[j]) for j in range(22) if w[j]])) for w in W]
    square=field.rem(V.sum_polys([mul(point[j],dual[j]) for j in range(22)]))
    return dual,vals,square


def solve(M,b):
    rows=[[F(x) for x in row]+[F(y)] for row,y in zip(M,b)];n=len(M)
    for col in range(n):
        pivot=next(i for i in range(col,n) if rows[i][col])
        rows[col],rows[pivot]=rows[pivot],rows[col];q=rows[col][col];rows[col]=[x/q for x in rows[col]]
        for i in range(n):
            if i!=col:
                q=rows[i][col];rows[i]=[a-q*b for a,b in zip(rows[i],rows[col])]
    return [r[-1] for r in rows]


def local_face(Q,W,G,h,i,j,k):
    ids=[i,j,k];M=[[G[a][b] for b in ids] for a in ids]
    weights=solve(M,[A.dot(W[l],h,Q) for l in ids]);z=[F(h[a])-sum(weights[b]*W[ids[b]][a] for b in range(3)) for a in range(22)]
    vals=[A.dot(w,z,Q) for w in W]
    assert A.dot(z,z,Q)<0 and A.dot(z,h,Q)<0
    assert [l for l,v in enumerate(vals) if v==0]==sorted(ids) and all(v<=0 for v in vals)


def verify_candidate(row,Q,W,G,h,sd):
    e,qe,v,gram,N,den,equations=sd
    field=Field(readp(row['minimal_polynomial']),*row['isolating_interval'])
    assert field.sign(qe)>0
    for m in row['repaired_ridge_labels']:
        assert m in equations and field.rem(equations[m])==Z
        if m==5:assert field.sign(sub(scale(N,2),den))>0
    assert row['repaired_ridge_labels']
    C=[readp(p) for p in row['interior_point']];dual,vals,square=field_products(C,Q,W,field)
    assert field.sign(square)<0 and all(field.sign(p)<0 for p in vals)
    assert field.rem(V.sum_polys([mul(dual[j],e[j]) for j in range(22)]))==Z
    assert field.sign(V.sum_polys([scale(dual[j],h[j]) for j in range(22)]))<0
    proof=row['actual_forbidden_ridge'];assert proof is not None
    a,b=[x-1 for x in proof['ridge_walls_one_based']]
    aa,bb,ab=map(field.rem,(gram(a,a),gram(b,b),gram(a,b)));D=field.rem(sub(mul(aa,bb),mul(ab,ab)))
    assert [aa,ab,bb]==[readp(p) for p in proof['projected_gram']]
    assert min(field.sign(aa),field.sign(bb),field.sign(D))>0
    point=[readp(p) for p in proof['point_coordinate_polynomials']];dual,vals,square=field_products(point,Q,W,field)
    assert square==readp(proof['point_square']) and field.sign(square)<0
    hp=field.rem(V.sum_polys([scale(dual[j],h[j]) for j in range(22)]))
    assert hp==readp(proof['point_h_pairing']) and field.sign(hp)<0
    assert field.rem(V.sum_polys([mul(dual[j],e[j]) for j in range(22)]))==Z
    signs=[field.sign(p) for p in vals];zeros=[l for l,s in enumerate(signs) if s==0]
    assert [l+1 for l in zeros]==proof['zero_walls_one_based'] and max(signs)<=0 and a in zeros and b in zeros
    cones=proof['redundant_zero_wall_cones'];assert sorted(r['wall_one_based']-1 for r in cones)==sorted(set(zeros)-{a,b})
    for r in cones:
        l=r['wall_one_based']-1;alpha=field.rem(sub(mul(bb,gram(l,a)),mul(ab,gram(l,b))));beta=field.rem(sub(mul(aa,gram(l,b)),mul(ab,gram(l,a))))
        assert alpha==readp(r['alpha_numerator']) and beta==readp(r['beta_numerator']) and D==readp(r['denominator'])
        assert min(field.sign(alpha),field.sign(beta))>=0
        assert field.rem(sub(sub(mul(D,gram(l,l)),mul(alpha,gram(l,a))),mul(beta,gram(l,b))))==Z
    N=field.rem(mul(ab,ab));den=field.rem(mul(aa,bb))
    if proof['angle_rejection']=='obtuse':assert field.sign(ab)>0
    else:
        assert field.sign(ab)<=0 and field.sign(sub(scale(den,3),scale(N,4)))>0 and field.sign(N)>0
        for p in (sub(scale(N,4),den),sub(scale(N,2),den),add(sub(scale(mul(N,N),16),scale(mul(N,den),12)),mul(den,den))):assert field.sign(p)!=0
    assert row['status']=='NON_COXETER_ACTUAL_RIDGE'
    return field


def main():
    V.sturm_regressions()
    Q,W,G,h=A.source();raw=json.loads(V.TARGETS.read_text());orbits=json.loads(V.ORB.read_text())['orbits']
    expected=set()
    for r in orbits:
        ids=[x-1 for x in r['representative']];a,b,c=ids
        det=G[a][a]*G[b][b]*G[c][c]+2*G[a][b]*G[a][c]*G[b][c]-G[a][a]*G[b][c]**2-G[b][b]*G[a][c]**2-G[c][c]*G[a][b]**2
        if G[a][a]*G[b][b]-G[a][b]**2<=0 or det<=0:continue
        for k in ids:
            i,j=[x for x in ids if x!=k]
            if G[i][i]==G[j][j] and G[i][k]==G[j][k]:expected.add(((i+1,j+1),k+1))
    assert len(expected)==17 and raw['complete'] and len(raw['pencils'])==17
    assert {(tuple(r['positive_walls_one_based']),r['negative_wall_one_based']) for r in raw['pencils']}==expected
    total=0;degrees=set()
    for index,pencil in enumerate(raw['pencils']):
        i,j=[x-1 for x in pencil['positive_walls_one_based']];k=pencil['negative_wall_one_based']-1
        local_face(Q,W,G,h,i,j,k);sd=section(Q,W,G,i,j,k);e,qe,v,gram,N,den,equations=sd
        assert {int(m):readp(p) for m,p in pencil['target_polynomials'].items()}==equations
        fields=[]
        for row in pencil['candidates']:
            field=verify_candidate(row,Q,W,G,h,sd);fields.append(field);degrees.add(len(field.p)-1);total+=1
            assert field.sign(gram(i,i))>0 and field.sign(gram(j,j))>0 and field.sign(gram(i,j))<=0
        for a,b in combinations(range(len(fields)),2):
            fa,fb=fields[a],fields[b];common=V.gcdpoly(fa.p,fb.p)
            if len(common)==1:continue
            lo=max(fa.lo,fb.lo);hi=min(fa.hi,fb.hi)
            if lo<hi:assert V.rootcount(common,lo,hi)==0
        for m,equation in equations.items():
            intervals=positive_roots(equation)
            if m==5:
                intervals=[r for r in intervals if Field(*r).sign(sub(scale(N,2),den))>0]
            actual=sum(m in row['repaired_ridge_labels'] for row in pencil['candidates'])
            assert actual==len(intervals),(i,j,k,m,actual,len(intervals))
        print('verified noncentral pencil',index+1,'candidates',len(fields),flush=True)
    assert total==raw['candidate_count']==64
    pencil=raw['pencils'][0];bad=deepcopy(pencil['candidates'][0]);bad['actual_forbidden_ridge']['point_coordinate_polynomials'][0][-1]=str(F(bad['actual_forbidden_ridge']['point_coordinate_polynomials'][0][-1])+1)
    i,j=[x-1 for x in pencil['positive_walls_one_based']];k=pencil['negative_wall_one_based']-1
    try:verify_candidate(bad,Q,W,G,h,section(Q,W,G,i,j,k))
    except AssertionError:mutation=True
    else:raise AssertionError('forged algebraic ridge point accepted')
    out={'status':'PASS','standard_library_only':True,'pencils':17,'target_labels':[2,3,4,5,6],
         'complete_positive_angle_repair_parameter_census':True,'algebraic_candidates':64,'parameter_degrees':sorted(degrees),
         'strict_timelike_interiors':64,'actual_forbidden_ridges':64,'forged_algebraic_point_rejected':mutation,
         'scope':'Declared17noncentral symmetric elliptic-triple pencils and labels2..6 only; no all-parameter noncentral theorem.'}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
