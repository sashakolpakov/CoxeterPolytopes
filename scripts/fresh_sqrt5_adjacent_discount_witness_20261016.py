#!/usr/bin/env python3
"""Exact arithmetic local quadrilateral and a candidate-payment failure.

Geometry uses a standalone Fraction Q(sqrt5) field. Payment compilation alone
uses the frozen proposed-rule graph helper and preserved decimal coefficients.
"""
import argparse
import json
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_adjacent_discount_witness_20261016'


@dataclass(frozen=True)
class K:
    a: F = F(0)
    b: F = F(0)
    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a)); object.__setattr__(self, 'b', F(self.b))
    def __add__(self, y):
        y = cast(y); return K(self.a+y.a, self.b+y.b)
    __radd__ = __add__
    def __neg__(self): return K(-self.a, -self.b)
    def __sub__(self, y): return self+-cast(y)
    def __rsub__(self, y): return cast(y)+-self
    def __mul__(self, y):
        y = cast(y); return K(self.a*y.a+5*self.b*y.b, self.a*y.b+self.b*y.a)
    __rmul__ = __mul__
    def __truediv__(self, y):
        y = cast(y); d = y.a*y.a-5*y.b*y.b; assert d
        return self*K(y.a/d, -y.b/d)
    def __bool__(self): return bool(self.a or self.b)
    def conj(self): return K(self.a, -self.b)
    def sign(self):
        if not self: return 0
        if self.a >= 0 and self.b >= 0: return 1
        if self.a <= 0 and self.b <= 0: return -1
        d = self.a*self.a-5*self.b*self.b
        assert d
        return (1 if self.a > 0 else -1)*(1 if d > 0 else -1)
    def enc(self): return [str(self.a), str(self.b)]
    def integral(self): return (self.a-self.b).denominator == 1 and (2*self.b).denominator == 1


def cast(x): return x if isinstance(x, K) else K(x)
Z, O = K(), K(1)


def dot(x, y): return sum((a*b for a,b in zip(x,y)), Z)
def mv(M, v): return [dot(r,v) for r in M]
def pair(M, x, y): return dot(x,mv(M,y))
def sub(M, ids): return [[M[i][j] for j in ids] for i in ids]
def enc(M): return [[x.enc() for x in r] for r in M]


def rref(M):
    A = [list(map(cast,r)) for r in M]; pivots=[]; k=0
    for j in range(len(A[0])):
        p = next((i for i in range(k,len(A)) if A[i][j]), None)
        if p is None: continue
        A[k],A[p] = A[p],A[k]; d=A[k][j]; A[k]=[x/d for x in A[k]]
        for i in range(len(A)):
            if i != k:
                d=A[i][j]; A[i]=[x-d*y for x,y in zip(A[i],A[k])]
        pivots.append(j); k+=1
        if k == len(A): break
    return A,pivots


def inverse(M):
    n=len(M); A,p=rref([r+[K(i==j) for j in range(n)] for i,r in enumerate(M)])
    assert p==list(range(n))
    return [r[n:] for r in A]


def kernel(M):
    A,p=rref(M); n=len(M[0]); out=[]
    for j in range(n):
        if j in p: continue
        v=[Z]*n;v[j]=O
        for i,k in enumerate(p):v[k]=-A[i][j]
        assert all(not x for x in mv(M,v));out.append(v)
    return out


def inertia(M):
    A=[r[:] for r in M]; ans=[0,0,0]
    while A:
        k=next((i for i in range(len(A)) if A[i][i]),None)
        if k is not None:
            d=A[k][k];ans[0 if d.sign()>0 else 1]+=1;ids=[i for i in range(len(A)) if i!=k]
            A=[[A[i][j]-A[i][k]*A[k][j]/d for j in ids] for i in ids]
        else:
            p=next(((i,j) for i,j in combinations(range(len(A)),2) if A[i][j]),None)
            if p is None:ans[2]+=len(A);break
            i,j=p;d=A[i][j];ans[0]+=1;ans[1]+=1;ids=[k for k in range(len(A)) if k not in p]
            A=[[A[k][l]-(A[k][i]*A[j][l]+A[k][j]*A[i][l])/d for l in ids] for k in ids]
    return ans


def geometry():
    phi=K(F(1,2),F(1,2));G=[[K(2*(i==j)) for j in range(8)] for i in range(8)];G[6][6]=K(4)
    bonds=[(0,5,5),(1,3,5),(2,3,3),(2,5,3),(3,7,3),(4,5,3),(6,7,4)]
    for i,j,m in bonds:G[i][j]=G[j][i]=-{3:O,4:K(2),5:phi}[m]
    assert inertia(G)==[6,1,1] and inertia([[x.conj() for x in r] for r in G])==[7,0,1]
    ids=list(range(1,8));Q=sub(G,ids);Qi=inverse(Q)
    N=[mv(Qi,[G[i][0] for i in ids])]+[[K(i==j) for j in range(7)] for i in range(7)]
    assert [[pair(Q,x,y) for y in N] for x in N]==G
    relation=[K(F(-3,4),F(-1,4)),K(F(1,4),F(1,4)),K(F(1,4),F(-1,4)),O,
              K(F(-1,4),F(-1,4)), -phi,K(F(1,2)),O]
    assert all(not x for x in mv(G,relation))
    assert all(not dot(relation,[r[j] for r in N]) for j in range(7))
    H=mv(Qi,[-O]*7); hp=[pair(Q,H,n) for n in N]; hq=pair(Q,H,H)
    assert hp[0]==K(F(-21,2),F(9,2)) and hp[1:]==[-O]*7
    assert hq==K(F(225,2),F(-113,2)) and hq.sign()<0 and max(x.sign() for x in hp)<0
    assert inertia(Q)==[6,1,0] and inertia([[x.conj() for x in r] for r in Q])==[7,0,0]
    white=[4,5,6,7];assert inertia(sub(G,white))==[4,0,0]
    vertices=[]
    for a,b in product((0,2),(1,3)):
        S=white+[a,b];SG=sub(G,S);assert inertia(SG)==[6,0,0]
        coeff=mv(inverse(SG),[hp[i] for i in S])
        v=[H[j]-sum((c*N[i][j] for c,i in zip(coeff,S)),Z) for j in range(7)]
        products=[pair(Q,v,n) for n in N]
        assert all(not products[i] for i in S) and all(products[i].sign()<0 for i in range(8) if i not in S)
        assert pair(Q,v,v).sign()<0 and pair(Q,v,H).sign()<0
        vertices.append(dict(black=[a,b],active=sorted(S),coordinates=[x.enc() for x in v],
                             norm=pair(Q,v,v).enc(),products=[x.enc() for x in products]))
    # Complete three-dimensional face cone: every extreme ray has two
    # independent black equalities in addition to the four white equalities.
    cone=[];accepted=[]
    for black in combinations(range(4),2):
        S=white+list(black);ker=kernel([mv(Q,N[i]) for i in S]);assert len(ker)==1
        orientations=[]
        for sign in (-1,1):
            v=[sign*x for x in ker[0]];p=[pair(Q,v,n) for n in N]
            if max(x.sign() for x in p)<=0:
                assert pair(Q,v,v).sign()<0 and pair(Q,v,H).sign()<0
                orientations.append(sign);accepted.append(list(black))
        cone.append(dict(black=list(black),feasible_orientations=orientations))
    assert sorted(accepted)==sorted([sorted(r['black']) for r in vertices]) and len(accepted)==4
    # Explicit spectator forest I2(5)^2 + A11, of rank15, realizes the exact
    # common context (0,2,0,0,0,0,0) at ambient rank21.
    S=[[K(2*(i==j)) for j in range(15)] for i in range(15)]
    for i,j in [(0,1),(2,3)]:S[i][j]=S[j][i]=-phi
    for i in range(4,14):S[i][i+1]=S[i+1][i]=-O
    assert inertia(S)==[15,0,0] and inertia([[x.conj() for x in r] for r in S])==[15,0,0]
    hs=[-x/K(1024) for x in mv(inverse(S),[O]*15)]
    assert mv(S,hs)==[-K(F(1,1024))]*15
    full_hq=hq+pair(S,hs,hs);assert full_hq.sign()<0
    big=[[Z]*23 for _ in range(23)]
    for i in range(8):
        for j in range(8):big[i][j]=G[i][j]
    for i in range(15):
        for j in range(15):big[8+i][8+j]=S[i][j]
    assert inertia(big)==[21,1,1] and inertia([[x.conj() for x in r] for r in big])==[22,0,1]
    # All reflections preserve the O_K-span of the23 generating roots.
    cartan=[[2*big[i][j]/big[j][j] for j in range(23)] for i in range(23)]
    assert all(x.integral() for r in cartan for x in r)
    for vertex in vertices:
        assert len(vertex['active'])+15==21 and 4 in vertex['active']
    return dict(field='a+b*sqrt(5)',scaled_gram=enc(G),basis_indices=ids,basis_gram=enc(Q),
                root_coordinates=enc(N),kernel=[x.enc() for x in relation],
                physical_inertia=[6,1,1],conjugate_inertia=[7,0,1],
                strict_center=[x.enc() for x in H],center_products=[x.enc() for x in hp],
                center_norm=hq.enc(),white=white,facet_mark=4,vertices=vertices,
                complete_face_cone=cone,spectator_forest=['I2(5)','I2(5)','A11'],
                spectator_gram=enc(S),full_center_norm=full_hq.enc(),
                rank21_physical_inertia=[21,1,1],rank21_conjugate_inertia=[22,0,1],
                integral_cartan_entries=23*23),bonds


def payment(bonds):
    import fresh_sqrt5_marked_h3_center_cross_20261016 as M
    path='data/enumerations/fresh_sqrt5_marked_h3_center_secondary_20261016_checkpoint.json'
    checkpoint=json.loads((ROOT/path).read_text());x=[F(str(z)) for z in checkpoint['latest_numeric_weights'][:M.N]]
    edges=bonds+[(8,9,5),(10,11,5)]+[(i,i+1,3) for i in range(12,22)]
    white=list(range(4,23));compiled=[M.compile_corner(white+list(p),edges,p,[4]) for p in product((0,2),(1,3))]
    common,tables=M.add_bonus(compiled,(0,)*7);profile=common
    for table in tables:profile=M.plus(profile,table[4])
    expected={-1:F(1),1:F(-1),14612:F(1),24603:F(1)}
    assert profile==expected and x[0]==0 and x[14612]==x[24603]==0
    paid=M.value(profile,x);assert paid==F(1820357345943253,5000000000000000)<1
    return dict(common_spectators=['I2(5)','I2(5)','A11'],external_counts=[0,2,0,0,0,0,0],
                facet_root=4,profile=[[j,str(z)] for j,z in sorted(profile.items())],
                exact_decimal_coefficients={str(j):str(x[j]) for j in (0,1,14612,24603)},
                payment=str(paid),required='1',necessary_symbolic_inequality='w14612+w24603 >= delta_adjacent3',
                coefficient_scope='Exact rational interpretation of a preserved numerical vector; not an exact feasible global primal.')


def generate():
    g,bonds=geometry();p=payment(bonds)
    paths=[f'scripts/{TAG}.py','results/fresh_sqrt5_adjacent_discount_frontier_20261016.json',
           'scripts/fresh_sqrt5_marked_h3_center_cross_20261016.py',
           'data/enumerations/fresh_sqrt5_marked_h3_center_secondary_20261016_checkpoint.json']
    return dict(status='PASS_EXACT_ARITHMETIC_LOCAL_QUAD_AND_PAYMENT_FAILURE',
                sources={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
                minimal_class=4202,geometry=g,payment=p,
                scope='Actual compact quadrilateral face of a local Coxeter chamber inH6 and its orthogonal-spectator extension inH21. All reflections preserve an admissible Q(sqrt5) integral root module. This does not establish a compact or finite-volume ambient polytope, global realizability inside one, or impossibility of repairing the payment model.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=generate();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text());print('PASS saved exact witness replay')
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['payment']['payment'])
