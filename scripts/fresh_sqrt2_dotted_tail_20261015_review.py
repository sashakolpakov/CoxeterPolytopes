#!/usr/bin/env python3
"""Independent Fraction/cycle-cover review of the C5 pendant-tail classification.

No producer imports or symbolic library. Double-tail integrality is proved by
the sqrt2 coefficient, independently of the producer's square-root membership.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt2_dotted_tail_20261015_review'
PRIMARY='fresh_sqrt2_dotted_tail_20261015'


class Q:
    def __init__(self,a=0,b=0):
        if isinstance(a,Q):self.a,self.b=a.a,a.b
        else:self.a,self.b=F(a),F(b)
    def __add__(self,y):
        y=Q(y);return Q(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,y):return self+-Q(y)
    def __rsub__(self,y):return Q(y)+-self
    def __mul__(self,y):
        y=Q(y);return Q(self.a*y.a+2*self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def __truediv__(self,y):
        y=Q(y);d=y.a*y.a-2*y.b*y.b;assert d
        z=self*Q(y.a,-y.b);return Q(z.a/d,z.b/d)
    def __eq__(self,y):
        y=Q(y);return (self.a,self.b)==(y.a,y.b)
    def __bool__(self):return bool(self.a or self.b)
    def sign(self):
        a,b=(self.a>0)-(self.a<0),(self.b>0)-(self.b<0)
        if not a:return b
        if not b or a==b:return a
        n=self.a*self.a-2*self.b*self.b
        return a*((n>0)-(n<0))
    def conj(self):return Q(self.a,-self.b)
    def integral(self):return self.a.denominator==self.b.denominator==1
    def enc(self):return [str(self.a),str(self.b)]


def sub(M,I):return [[M[i][j]for j in I]for i in I]


def inertia(M):
    # Congruence using only one-dimensional pivots. A zero-diagonal block
    # with entry(i,j)!=0 is first changed by e_i -> e_i+e_j.
    A=[[Q(z)for z in row]for row in M];counts=[0,0,0]
    while A:
        n=len(A);k=next((i for i in range(n)if A[i][i]),None)
        if k is None:
            edge=next(((i,j)for i,j in combinations(range(n),2)if A[i][j]),None)
            if edge is None:counts[2]+=n;break
            i,j=edge
            for h in range(n):A[i][h]=A[i][h]+A[j][h]
            for h in range(n):A[h][i]=A[h][i]+A[h][j]
            k=i
        d=A[k][k];assert d
        counts[0 if d.sign()>0 else 1]+=1
        I=[i for i in range(n)if i!=k]
        A=[[A[i][j]-A[i][k]*A[k][j]/d for j in I]for i in I]
    return counts


def inverse(M):
    n=len(M);A=[[Q(x)for x in row]+[Q(i==j)for j in range(n)]for i,row in enumerate(M)]
    for k in range(n):
        i=next(i for i in range(k,n)if A[i][k]);A[i],A[k]=A[k],A[i]
        d=A[k][k];A[k]=[z/d for z in A[k]]
        for j in range(n):
            if j!=k:
                d=A[j][k];A[j]=[a-d*b for a,b in zip(A[j],A[k])]
    return [r[n:]for r in A]


def mv(M,x):return [sum((a*b for a,b in zip(row,x)),Q())for row in M]
def dot(x,y,M):return sum((a*b for a,b in zip(x,mv(M,y))),Q())


def symbolic_determinant(double):
    # Each matrix entry is a monomial in p,q,r, with coefficient in Qsqrt2.
    z=(0,0,0);M=[[None]*7 for _ in range(7)]
    for i in range(7):M[i][i]=(Q(2),(1,0,0)if i==0 else(0,1,0)if i==1 else z)
    for i,j,v in [(2,3,Q(0,-1)),(2,4,Q(-1)),(4,5,Q(-1)),(5,6,Q(-1)),(6,3,Q(-1)),(0,5,Q(-1))]+([(1,5,Q(-1))]if double else[]):
        M[i][j]=M[j][i]=(v,z)
    M[0][1]=M[1][0]=(Q(-1),(0,0,1))
    result={}
    for perm in permutations(range(7)):
        if any(M[i][j]is None for i,j in enumerate(perm)):continue
        sign=(-1)**sum(perm[i]>perm[j]for i,j in combinations(range(7),2))
        value=Q(sign);powers=[0,0,0]
        for i,j in enumerate(perm):
            c,e=M[i][j];value=value*c;powers=[a+b for a,b in zip(powers,e)]
        key=tuple(powers);result[key]=result.get(key,Q())+value
    result={k:v for k,v in result.items()if v}
    factor=2*Q(1,1)
    expected={(0,0,2):factor,(1,1,0):-4*factor}
    if double:
        for key in [(0,0,1),(1,0,0),(0,1,0)]:expected[key]=factor*Q(1,-1)
    else:expected[(0,1,0)]=factor*Q(1,-1)
    assert result==expected
    return [dict(powers=list(k),coefficient=v.enc())for k,v in sorted(result.items())]


def literal(ell,double,contact):
    n=7+ell;M=[[Q(2*(i==j))for j in range(n)]for i in range(n)]
    E=[(2,3,Q(0,-1)),(2,4,Q(-1)),(4,5,Q(-1)),(5,6,Q(-1)),(6,3,Q(-1)),(0,5,Q(-1))]
    if double:E.extend([(1,5,Q(-1)),(0,1,-contact)])
    else:
        M[1][1]=2*contact;E.append((0,1,-contact))
    for k in range(ell):E.append((0 if k==0 else 6+k,7+k,Q(-1)))
    for i,j,v in E:M[i][j]=M[j][i]=v
    return M


def local_geometry(ell,double,contact,saved):
    M=literal(ell,double,contact);n=len(M);W=[4,5,6]+list(range(7,n));B=[0,2,1,3]
    assert [[z.enc()for z in row]for row in M]==saved['gram']
    assert inertia(M)==saved['physical_inertia']==[n-2,1,1]
    assert inertia([[z.conj()for z in row]for row in M])==saved['conjugate_inertia']==[n-1,0,1]
    tests=0
    for size in range(1,n+1):
        for I in combinations(range(n),size):
            iv=inertia(sub(M,I));assert iv[1]or not iv[2];tests+=1
    assert tests==saved['no_affine_subset_tests']
    J=inverse(sub(M,W));R=[[M[a][b]-dot([M[a][w]for w in W],[M[b][w]for w in W],J)for b in B]for a in B]
    G=sub(R,[0,1,2]);assert inertia(G)==[2,1,0]
    rays=[];pairings=[];weights=[]
    for i,j in [(0,1),(1,2),(2,3),(3,0)]:
        I=W+[B[i],B[j]]
        assert inertia(sub(M,I))==[len(W)+2,0,0]
        adjacency={u:[v for v in I if v!=u and M[u][v]]for u in I}
        reached={B[i]};stack=[B[i]]
        while stack:
            u=stack.pop()
            for v in adjacency[u]:
                if v not in reached:reached.add(v);stack.append(v)
        weight=F(0);component=None
        if B[j]in reached:
            center=[u for u in reached if len(adjacency[u])==3]
            assert len(center)==1
            paths={center[0]:[]};stack=[center[0]]
            while stack:
                u=stack.pop()
                for v in adjacency[u]:
                    if v not in paths:paths[v]=paths[u]+[v];stack.append(v)
            arms=sorted(len(paths[u])for u in reached if len(adjacency[u])==1)
            component='D5'if arms==[1,1,2]else'E6'if arms==[1,2,2]else None
            assert component is not None
            marked_paths=[paths[B[i]],paths[B[j]]]
            assert sorted(map(len,marked_paths))==[1,2]
            assert marked_paths[0][0]!=marked_paths[1][0]
            if component=='D5':assert all(len(adjacency[u])==1 for u in[B[i],B[j]])
            else:assert sorted(len(adjacency[u])for u in[B[i],B[j]])==[1,2]
            weight=F(1,2)
        weights.append(dict(marked=[B[i],B[j]],entire_component=component,halfweight=str(weight)))
        unused=next(k for k in range(3)if k not in[i,j])
        rows=[[R[k][h]for k in range(3)]for h in[i,j,unused]]
        v=mv(inverse(rows),[Q(),Q(),Q(-1)])
        ps=[sum((v[k]*R[k][h]for k in range(3)),Q())for h in range(4)]
        assert all(not ps[k]if k in[i,j]else ps[k].sign()<0 for k in range(4))
        assert dot(v,v,G).sign()<0
        rays.append(v);pairings.append(ps)
    assert all(dot(v,w,G).sign()<0 for v in rays for w in rays)
    assert sum(F(row['halfweight'])for row in weights)>=1
    for pair in [(0,2),(1,3)]:assert inertia(sub(M,W+[B[k]for k in pair]))[1]==1
    return dict(ell=ell,double=double,contact=contact.enc(),physical_inertia=inertia(M),
                conjugate_inertia=inertia([[z.conj()for z in row]for row in M]),
                principal_no_affine_tests=tests,residual_Gram=[[z.enc()for z in row]for row in R],
                four_future_vertex_rays=[[z.enc()for z in v]for v in rays],
                strict_other_side_pairings=[[z.enc()for z in v]for v in pairings],
                literal_entire_component_halfweights=weights)


def run():
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    symbols=[symbolic_determinant(d)for d in(False,True)]
    ps=[F(1),F(3,4),F(2,3),F(5,8)]
    tails=[]
    for ell in range(5):
        if ell:
            A=[[Q(2*(i==j)-int(abs(i-j)==1))for j in range(ell)]for i in range(ell)]
            p=1-inverse(A)[0][0]/2
        else:p=Q(1)
        assert p==F(ell+2,2*(ell+1))
        assert p==F(primary['tail_diagonals'][ell]['reduced_unit_diagonal'])
        M=literal(ell,False,Q(3,1));I=[4,5,6]+list(range(7,7+ell))+[0,2]
        iv=inertia(sub(M,I));assert iv==([ell+5,0,0]if ell<4 else[8,0,1])
        assert iv==primary['corner_tails'][ell]['inertia']
        tails.append(dict(length=ell,reduced_unit_diagonal=p.enc(),corner_inertia=iv))
    # Integral double-cycle contact forces p+q integral, without extracting
    # any quadratic radicals: coefficient sqrt2 is2AB+B-A-(p+q).
    double_cases=[]
    for ell,m in product(range(4),repeat=2):
        p,q=ps[ell],ps[m];integral_sum=(p+q).denominator==1
        assert integral_sum==(ell==m==0)
        row=primary['double_family_cases'][len(double_cases)]
        assert row['tails']==[ell,m] and row['integral']==integral_sum
        if row['contact'] is not None:
            r=Q(*row['contact'])
            assert r*r+Q(1,-1)*r-4*p*q+Q(1,-1)*(p+q)==0
            assert r.sign()>0 and r.integral()==integral_sum
        double_cases.append(dict(tails=[ell,m],sum_reduced_diagonals=str(p+q),necessary_integral_sum=integral_sum))
    r=Q(1,1);assert r*r+Q(1,-1)*r-4+2*Q(1,-1)==0
    assert Q(-2)*Q(-2)+Q(1,-1)*Q(-2)-4+2*Q(1,-1)==0
    assert r.sign()>0
    single=[]
    for ell in range(4):
        R=F(ell+3,ell+1);D=Q(R,1);special=Q(F(1,2),F(-1,4))*D
        assert not special.integral()
        assert D.integral()==(ell<=1)
        saved=primary['single_family_cases'][ell]
        assert D==Q(*saved['contact_square']) and special==Q(*saved['I2_8_y_attachment_contact_square'])
        assert saved['integral']==D.integral()
        single.append(dict(length=ell,at_q_one=D.enc(),at_I2_8_q=special.enc(),q_one_integral=D.integral()))
    # Exact affine fork obstruction for pendants at either4-bond endpoint.
    fork=[[Q(i==j)for j in range(4)]for i in range(4)]
    for j,c in [(1,Q(0,F(-1,2))),(2,Q(F(-1,2))),(3,Q(F(-1,2)))]:fork[0][j]=fork[j][0]=c
    null=[Q(1),Q(0,F(1,2)),Q(F(1,2)),Q(F(1,2))]
    assert all(not z for z in mv(fork,null)) and all(z.sign()>0 for z in null)
    assert inertia(fork)==[3,0,1]
    local=[local_geometry(ell,d,c,saved)for(ell,d,c),saved in zip([(0,False,Q(3,1)),(1,False,Q(2,1)),(0,True,Q(1,1))],primary['allowed_local_stars'])]
    assert sum(r['principal_no_affine_tests']for r in local)==509
    for p,h in primary['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    paths=[f'scripts/{TAG}.py',f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json']
    return dict(status='PASS_INDEPENDENT_CYCLE_COVER_FRACTION_TAIL_CLASSIFICATION_AND_HOSTILE_SCOPE_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                symbolic_determinants=symbols,tail_diagonals_and_corner_inertia=tails,
                single_family=single,double_integral_sum_obstruction=double_cases,
                double_positive_integral_contact=r.enc(),affine_fork_null=[z.enc()for z in null],
                local_stars=local,no_affine_principal_tests=509,
                independent_methods=['permutation cycle-cover polynomial expansion; no SymPy',
                                     'direct integral sqrt2-coefficient obstruction for15 of16 double-tail cases',
                                     'one-dimensional Fraction congruence pivots and exact residual vertex signs'],
                scope='Complete only for the specified C5 single/double middle attachment and pendant white trees attached to one core node each. No multiple-core external component, extra core contact, different core labels or global dimension result is certified.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text());print('PASS saved replay')
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
