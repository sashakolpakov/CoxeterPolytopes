#!/usr/bin/env python3
"""Separate recursive-Fraction E7 face replay; no producer arithmetic imports."""
import argparse
import json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_face_20261014_review'
FACE = 'results/fresh_h20_e7_face_20261014.json'
BASE = 'results/fresh_h20_e8_interior_caps_20261014_exact.json'
RADICALS = (2, 79, 11)
ZERO = (F(0),)*8


def sha(p):
    return sha256((ROOT/p).read_bytes()).hexdigest()


def load(p):
    d = json.loads((ROOT/p).read_text())
    for f, h in d.get('sources', {}).items():
        assert sha(f) == h, f
    return d


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def scale(a,s):
    return tuple(x*s for x in a)


def mul(a,b):
    if len(a) == 1:
        return (a[0]*b[0],)
    n = len(a)//2
    d = RADICALS[n.bit_length()-1]
    aa, bb = mul(a[:n],b[:n]), mul(a[n:],b[n:])
    ab, ba = mul(a[:n],b[n:]), mul(a[n:],b[:n])
    return add(aa,scale(bb,d)) + add(ab,ba)


def sign(a):
    if len(a) == 1:
        return (a[0]>0)-(a[0]<0)
    n = len(a)//2
    x, y = a[:n], a[n:]
    sx, sy = sign(x), sign(y)
    if not sx:
        return sy
    if not sy or sx == sy:
        return sx
    d = RADICALS[n.bit_length()-1]
    return sx*sign(add(mul(x,x),scale(mul(y,y),-d)))


def c(x):
    return (F(x),)+(F(0),)*7


def vec(v):
    return tuple(tuple(map(F,x)) for x in v)


def enc(x):
    if isinstance(x,F):
        return str(x)
    if isinstance(x,dict):
        return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [enc(v) for v in x]
    return x


def run():
    d,f = load(BASE),load(FACE)
    Q = [[F(x) for x in row] for row in d['form']]
    N = list(map(vec,d['normals']))
    H = vec(d['strict_center'])
    assert len(N) == 87 and len(Q) == 21

    @lru_cache(None)
    def metric(y):
        return tuple(tuple(sum(Q[i][j]*y[j][k] for j in range(21))
                           for k in range(8)) for i in range(21))

    @lru_cache(None)
    def dot(x,y):
        z=ZERO
        for a,b in zip(x,metric(y)):
            if a != ZERO and b != ZERO:
                z=add(z,mul(a,b))
        return z

    def linear(rows,co):
        v=[ZERO]*21
        for a,r in zip(co,rows):
            for i,x in enumerate(r):
                if a != ZERO and x != ZERO:
                    v[i]=add(v[i],mul(a,x))
        return tuple(v)

    S = f['E7_support']
    A = f['face_walls']
    assert S == [12,13,14,15,19,33,67]
    assert A == list(range(11))+[18,83]
    B = sp.Matrix([[dot(N[i],N[j])[0] for j in S] for i in S])
    assert all(dot(N[i],N[j]) == c(B[a,b]) for a,i in enumerate(S) for b,j in enumerate(S))
    assert B.tolist() == [[F(x) for x in row] for row in f['gram']]
    _,D = B.LDLdecomposition(hermitian=False)
    assert all(D[i,i]>0 for i in range(7)) and B.det() == 2
    I = [[F(str(x)) for x in row] for row in B.inv().tolist()]
    assert I == [[F(x) for x in row] for row in f['inverse']]
    assert all(x>0 for row in I for x in row)
    diag=[I[i][i] for i in range(7)]
    assert diag == list(map(F,['7/2','12','6','2','15/2','3/2','4']))
    adj={i:[j for j in range(7) if i!=j and B[i,j]] for i in range(7)}
    fork=next(i for i in adj if len(adj[i])==3)
    arms=[]
    for j in adj[fork]:
        prev,cur,k=fork,j,1
        while len(adj[cur])==2:
            nxt=next(x for x in adj[cur] if x!=prev)
            prev,cur,k=cur,nxt,k+1
        arms.append(k)
    assert sorted(arms)==[1,2,3]

    def project(v):
        p=[dot(v,N[i]) for i in S]
        co=[tuple(sum(I[a][j]*p[j][k] for j in range(7)) for k in range(8)) for a in range(7)]
        shift=linear([N[i] for i in S],co)
        return tuple(add(x,scale(y,-1)) for x,y in zip(v,shift))

    P=project(H)
    assert P==vec(f['face_strict_center'])
    assert sign(dot(P,P))<0 and sign(dot(P,H))<0
    assert all(dot(P,N[i])==ZERO if i in S else sign(dot(P,N[i]))<0 for i in range(87))
    causal={r['wall']:r for r in f['causal_nonfacets']}
    kinds={'zero':[],'positive':[],'negative':[]}
    for i,r in enumerate(N):
        p=project(r)
        assert p==vec(f['projected_walls'][i])
        assert all(dot(p,N[j])==ZERO for j in S)
        q=dot(p,p);sg=sign(q)
        if sg==0:
            assert p==(ZERO,)*21
            kinds['zero'].append(i)
        elif sg>0:
            assert p==r
            kinds['positive'].append(i)
        else:
            assert i in causal and q==tuple(map(F,causal[i]['projection_norm']))
            assert dot(p,P)==tuple(map(F,causal[i]['orientation'])) and sign(dot(p,P))<0
            kinds['negative'].append(i)
    assert kinds['zero']==S and kinds['positive']==A
    assert len(kinds['negative'])==67 and set(causal)==set(kinds['negative'])
    C=sp.Matrix([[2 if i==j else -1 if ((i<11 and j<11) or (i>=11 and j>=11)) and abs(i-j)==1 else 0
                  for j in range(13)] for i in range(13)])
    assert C.det()==36
    assert all(dot(N[i],N[j])==c(C[a,b]) for a,i in enumerate(A) for b,j in enumerate(A))
    assert C.tolist()==[[F(x) for x in row] for row in f['face_gram']]
    _,DD=C.LDLdecomposition(hermitian=False)
    assert all(DD[i,i]>0 for i in range(13))
    V=vec(f['finite_apex'])
    assert dot(V,V)==tuple(map(F,f['apex_norm'])) and sign(dot(V,V))<0 and sign(dot(V,H))<0
    assert all(dot(V,N[i])==ZERO if i in S+A else sign(dot(V,N[i]))<0 for i in range(87))
    CI=[[F(str(x)) for x in row] for row in C.inv().tolist()]
    Y=linear([N[i] for i in A],[c(-sum(row)) for row in CI])
    assert Y==vec(f['strict_outgoing_spatial_direction'])
    assert dot(Y,Y)==c(145) and dot(Y,V)==ZERO
    assert all(dot(Y,N[i])==c(-1) for i in A) and all(dot(Y,N[i])==ZERO for i in S)
    normalization=f['intrinsic_normalization']
    s=tuple(map(F,normalization['apex_scale']))
    assert s==tuple(map(F,['-267/2','0','-51/8','0','5537/24','0','89/6','0']))
    assert sign(s)>0 and dot(V,V)==scale(mul(s,s),-2)
    T=vec(normalization['timelike_basis'])
    assert V==tuple(mul(s,x) for x in T)
    assert dot(T,T)==c(-2) and sign(dot(T,H))<0
    assert all(dot(T,N[i])==ZERO for i in S+A)
    intrinsic=sp.diag(C,sp.Matrix([[-2]]))
    assert intrinsic.det()==-72
    assert intrinsic.tolist()==[[F(x) for x in row] for row in normalization['form']]
    assert F(normalization['determinant'])==-72 and normalization['spatial_basis_walls']==A
    basis=[N[i] for i in A]+[T]
    assert normalization['embedding_field']=='Q(sqrt79,sqrt11)'
    assert normalization['ambient_coefficient_basis']==[1,2,79,158,11,22,869,1738]
    assert all(x[j]==0 for b in basis for x in b for j in (1,3,5,7))
    assert all(dot(basis[i],basis[j])==c(intrinsic[i,j]) for i in range(14) for j in range(14))
    assert len(normalization['projected_wall_coordinates'])==87
    for i,r in enumerate(N):
        p=project(r)
        pairs=[dot(p,N[j]) for j in A]
        co=[tuple(sum(CI[a][b]*pairs[b][k] for b in range(13)) for k in range(8)) for a in range(13)]
        co.append(scale(dot(p,T),F(-1,2)))
        assert tuple(co)==vec(normalization['projected_wall_coordinates'][i])
        assert linear(basis,co)==p
    # The all-label Schur gate uses only these rational sharp lower bounds.
    single_minimum_cost=[x/2 for x in diag]
    assert [S[i] for i,x in enumerate(single_minimum_cost) if x<1]==[33]
    assert I[5][5]==F(3,2) and I[5][5]*F(1,2)==F(3,4)
    assert I[5][5]*1>1  # all labels m>=4 have 2*cos(pi/m)^2 >= 1
    assert all((I[i][i]+I[j][j]+2*I[i][j])/2>1 for i in range(7) for j in range(i))
    omega=linear([N[i] for i in S],[c(x) for x in I[5]])
    profile=f['unique_nonright_meeting_profile']
    assert omega==vec(profile['fundamental_weight']) and dot(omega,omega)==c(F(3,2))
    assert all(dot(omega,N[i])==c(int(i==33)) for i in S)
    assert profile['wall']==33 and profile['label']==3 and profile['other_six_right']
    assert F(profile['projected_norm_fraction'])==F(1,4)
    # For an intrinsic normal of square z>0, add -sqrt(2z)*omega:
    # total square z+2z*(3/2)=4z, and its normalized wall33 cosine is 1/2.
    assert 1+2*F(3,2)==4 and F(2,4*2)==F(1,4)
    paths=[BASE,FACE,f'scripts/{TAG}.py','scripts/fresh_h20_e7_face_20261014.py']
    return enc(dict(status='PASS_SEPARATE_FRACTION_E7_FACE_AND_ALL_LABEL_SCHUR_REVIEW',
                    sources={p:sha(p) for p in paths},
                    arithmetic='Recursive exact Fraction tower Q(sqrt2,sqrt79,sqrt11); no producer imports or numerical sign tests. SymPy independently inverts the rational Gram matrices.',
                    E7_determinant=2,E7_arm_lengths=sorted(arms),inverse_diagonal=diag,
                    all_strict_positive_inverse_entries=49,projections_reconstructed=87,
                    future_timelike_nonfacets=67,zero_projections=7,actual_face_facets=13,
                    intrinsic_spatial_type='A11+A2',intrinsic_spatial_determinant=36,
                    intrinsic_form_determinant=-72,intrinsic_embedding_field='Q(sqrt79,sqrt11)',
                    reduced_projection_coordinates_reconstructed=87,
                    apex_inequalities_checked=87,direction_norm=145,
                    minimum_single_contact_costs=single_minimum_cost,
                    unique_positive_residual_nonright_profile=dict(wall=33,label=3,residual_fraction=F(1,4)),
                    scope='Narrow separate arithmetic and proof review of the E7 face, all 87 projections, one actual apex and the all-label face-facet profile gate. It does not repeat the 87-wall pair or complete ambient vertex/edge census and does not prove ambient noncompletion.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert out==json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
