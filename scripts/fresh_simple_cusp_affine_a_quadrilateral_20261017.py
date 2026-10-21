#!/usr/bin/env python3
"""All-rank affine-A quadrilateral payment; no new ambient dimension bound."""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import permutations, product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_simple_cusp_affine_a_quadrilateral_20261017'
INTERFACE='fresh_simple_cusp_vinberg_weight_frontier_20261017'
AFFINE='fresh_simple_cusp_affine_exclusion_20261017'

class Q:
    """The ordered physical field Q(sqrt(2)), with exact rational signs."""
    def __init__(self,a=0,b=0):self.a,self.b=F(a),F(b)
    def __add__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return Q(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,o):return self+-coerce(o)
    def __rsub__(self,o):return coerce(o)+-self
    def __mul__(self,o):
        o=coerce(o);return Q(self.a*o.a+2*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=coerce(o);d=o.a*o.a-2*o.b*o.b
        assert d
        return self*Q(o.a/d,-o.b/d)
    def __pow__(self,k):
        assert k>=0
        r=Q(1)
        for _ in range(k):r=r*self
        return r
    def __eq__(self,o):
        o=coerce(o);return self.a==o.a and self.b==o.b
    def sign(self):
        a,b=self.a,self.b
        if not b:return (a>0)-(a<0)
        if not a:return (b>0)-(b<0)
        if a*b>0:return (a>0)-(a<0)
        d=a*a-2*b*b
        assert d
        return ((a>0)-(a<0))*((d>0)-(d<0))
    def wire(self):return [str(self.a),str(self.b)]

def coerce(x):return x if isinstance(x,Q) else Q(x)
Z=Q();ONE=Q(1);SQ=Q(0,1)

def positive_definite(M):
    A=[r[:] for r in M]
    for k in range(len(A)):
        d=A[k][k]
        if d.sign()<=0:return False
        for i in range(k+1,len(A)):
            for j in range(k+1,len(A)):
                A[i][j]=A[i][j]-A[i][k]*A[k][j]/d
    return True

def affine_witness(name,n,edges,kernel):
    M=[[Q(2) if i==j else Z for j in range(n)] for i in range(n)]
    for a,b,w in edges:M[a][b]=M[b][a]=-w
    assert all(x.sign()>0 for x in kernel)
    assert all(sum((M[i][j]*kernel[j] for j in range(n)),Z)==0 for i in range(n))
    assert all(positive_definite([[M[i][j] for j in range(n) if j!=k]
                                 for i in range(n) if i!=k]) for k in range(n))
    return dict(name=name,order=n,edges=[[a,b,w.wire()] for a,b,w in edges],
                positive_kernel=[x.wire() for x in kernel],all_one_deletions_PD=True)

# Small exact multivariate polynomials, for the independent formal determinant expansion.
NV=8;ZERO_MON=(0,)*NV
def pc(c):return {} if not c else {ZERO_MON:F(c)}
def pv(i):
    e=[0]*NV;e[i]=1;return {tuple(e):F(1)}
def pa(*args):
    out={}
    for p in args:
        for e,c in p.items():out[e]=out.get(e,F(0))+c
    return {e:c for e,c in out.items() if c}
def pn(p):return {e:-c for e,c in p.items()}
def pm(p,q):
    out={}
    for e,c in p.items():
        for f,d in q.items():
            k=tuple(x+y for x,y in zip(e,f));out[k]=out.get(k,F(0))+c*d
    return {e:c for e,c in out.items() if c}
def pdet(M):
    out={};n=len(M)
    for p in permutations(range(n)):
        inv=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        t=pc(-1 if inv%2 else 1)
        for i,j in enumerate(p):t=pm(t,M[i][j])
        out=pa(out,t)
    return out

def determinant_identity():
    alpha,A,B,C,D,E,h,k=[pv(i) for i in range(NV)]
    M=[[alpha,pn(alpha),pa(B,pn(h)),pn(A)],
       [pn(alpha),alpha,pn(B),pa(A,pn(k))],
       [pa(B,pn(h)),pn(B),C,pn(E)],
       [pn(A),pa(A,pn(k)),pn(E),D]]
    l=pa(pm(h,k),pn(pm(k,B)),pn(pm(h,A)))
    rhs=pa(pm(l,l),pn(pm(alpha,pa(pm(C,pm(k,k)),pm(D,pm(h,h)),
                              pm(pc(2),pm(E,pm(h,k)))))))
    actual=pdet(M);assert actual==rhs
    return dict(variables=['alpha','A','B','C','D','E','h','k'],
                terms=[dict(exponents=list(e),coefficient=str(c)) for e,c in sorted(actual.items())],
                identity='det S=(hk-kB-hA)^2-alpha*(C*k^2+D*h^2+2*E*h*k)')

def profiles():
    # moment=(magnitude times distance from the replaced old black on each arc).
    u=SQ-1
    return [dict(name='isolated',c0=Q(2),moment=[Z,Z],umax=SQ/2,minimum_height=SQ),
        *[dict(name=f'end3_{j}',c0=ONE,moment=[ONE if i==j else Z for i in range(2)],
               umax=u,minimum_height=1+SQ) for j in range(2)],
        *[dict(name=f'end4_{j}',c0=Z,moment=[SQ if i==j else Z for i in range(2)],
               umax=u,minimum_height=1+SQ) for j in range(2)],
        *[dict(name=f'next3_{j}',c0=Z,moment=[Q(2) if i==j else Z for i in range(2)],
               umax=u,minimum_height=1+SQ) for j in range(2)],
        dict(name='both_ends3',c0=Z,moment=[ONE,ONE],umax=1-SQ/2,minimum_height=2+SQ)]

def value(a,b,c,d,u,v):
    l=u*c['moment'][0]+v*d['moment'][0]
    m=u*c['moment'][1]+v*d['moment'][1]
    return 1-2*(a*l+b*m)-(a+b)*(c['c0']*u*u+d['c0']*v*v)-a*b*(l-m)**2

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    interface=json.loads((ROOT/f'results/{INTERFACE}.json').read_text())
    assert interface['status']=='PASS_AFFINE_3N_WEIGHT_INTERFACE_QUADRILATERAL_PAYMENT_OPEN'
    for p,h in interface['sources'].items():assert sha(p)==h,p
    witnesses=[affine_witness('A2',3,[(0,1,ONE),(1,2,ONE),(2,0,ONE)],[ONE]*3),
        affine_witness('A3',4,[(i,(i+1)%4,ONE) for i in range(4)],[ONE]*4),
        affine_witness('B3',4,[(0,1,SQ),(0,2,ONE),(0,3,ONE)],[Q(2),SQ,ONE,ONE])]
    e8edges=[];e8kernel=[Q(6)];nextnode=1
    for length in (1,2,5):
        previous=0
        for j in range(1,length+1):
            e8edges.append((previous,nextnode,ONE))
            e8kernel.append(Q(F(6*(length+1-j),length+1)))
            previous=nextnode;nextnode+=1
    witnesses.append(affine_witness('E8',9,e8edges,e8kernel))
    records=[];zeros=[];P=profiles()
    for c,d in product(P,repeat=2):
        for ub,vb in product((0,1),repeat=2):
            u=c['umax'] if ub else Z;v=d['umax'] if vb else Z
            f=value(Q(F(1,4)),Q(F(1,4)),c,d,u,v)
            assert f.sign()>=0,(c['name'],d['name'],ub,vb,f.wire())
            row=dict(c=c['name'],d=d['name'],upper_u=ub,upper_v=vb,F=f.wire())
            records.append(row)
            if f==0:zeros.append(row)
    assert len(records)==256 and len(zeros)==3
    # Along either reciprocal arc parameter, F is strictly decreasing for u,v>0.
    # Nonnegative coefficients are literal in the displayed formula; the only
    # zero-moment profile is isolated and has c0=2, so both derivatives are strict.
    assert all(p['c0'].sign()>=0 and all(x.sign()>=0 for x in p['moment']) for p in P)
    assert all(any(x.sign()>0 for x in p['moment']) or p['c0'].sign()>0 for p in P)
    for p in P:assert p['minimum_height']*p['umax']==1
    paths=[f'scripts/{TAG}.py',f'scripts/{INTERFACE}.py',f'results/{INTERFACE}.json',
           f'results/{INTERFACE}.md',f'scripts/{AFFINE}.py',f'results/{AFFINE}.json',
           f'results/{AFFINE}.md',f'results/{AFFINE}_review.md']
    wireprofiles=[{k:([x.wire() for x in v] if isinstance(v,list) else v.wire()
                     if isinstance(v,Q) else v) for k,v in p.items()} for p in P]
    return dict(status='PASS_ALL_RANK_AFFINE_A_QUADRILATERAL_PAYMENT',
        sources={p:sha(p) for p in paths},minimum_dimension=14,
        profiles=wireprofiles,affine_exclusion_witnesses=witnesses,
        formal_determinant=determinant_identity(),rectangle_checks=records,zero_checks=zeros,
        counts=dict(profiles=8,ordered_profile_pairs=64,rectangle_corners=256,
                    zero_rectangle_values=3,small_affine_witnesses=4),
        all_rank_argument=dict(arc_lengths='P,Q>=4; P+Q=n>=14',
            reciprocal_relaxation='0<a=1/P<=1/4,0<b=1/Q<=1/4; at least one strict',
            height_rectangle='0<u<=profile.umax,0<v<=profile.umax',
            joint_concavity='negative diagonal quadratic plus negative square',
            strict_monotonicity='dF/da<0 and dF/db<0 for u,v>0',
            conclusion='native c.d=0 would force detS/(h^2*k^2)=F=0, but F>0'),
        theorem='At an actual quadrilateral in an all-simple finite-volume Coxeter n-polytope, n>=14, an affine A_(n-1) corner is paid by distance<=3, or its opposite corner has adjacent marked roots and is paid.',
        unresolved='Quadrilaterals with other affine corner types; no new ambient dimension bound.',
        scope='Complete all-rank algebraic A-cycle branch. Finite profile and rectangle checks certify analytic arguments, not a tail cutoff. Inputs include the inherited global exclusion of affine subdiagrams of order below n.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],result['counts'])
