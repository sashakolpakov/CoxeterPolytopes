#!/usr/bin/env python3
"""Exact changed-local-Witt H9 host, norm pool and complete zero chamber.

No Vinberg height search. Standard-library Fraction arithmetic only.
The local-field deductions and all-parameter splitting lemma are in the MD.
"""
import argparse
import json
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_witt_twist_20261016'

class K:
    def __init__(self,a=0,b=0):self.a,self.b=Q(a),Q(b)
    @staticmethod
    def of(x):return x if isinstance(x,K) else K(x)
    def __add__(self,x):x=K.of(x);return K(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self):return K(-self.a,-self.b)
    def __sub__(self,x):return self+-K.of(x)
    def __rsub__(self,x):return K.of(x)+-self
    def __mul__(self,x):
        x=K.of(x);return K(self.a*x.a+self.b*x.b,self.a*x.b+self.b*x.a+self.b*x.b)
    __rmul__=__mul__
    def conjugate(self):return K(self.a+self.b,-self.b)
    def norm(self):return self.a*self.a+self.a*self.b-self.b*self.b
    def inverse(self):
        n=self.norm();assert n;z=self.conjugate();return K(z.a/n,z.b/n)
    def __truediv__(self,x):return self*K.of(x).inverse()
    def __eq__(self,x):x=K.of(x);return self.a==x.a and self.b==x.b
    def sign(self):
        a,b=2*self.a+self.b,self.b
        s=lambda t:(t>0)-(t<0)
        if not a:return s(b)
        if not b or s(a)==s(b):return s(a)
        return s(a)*s(a*a-5*b*b)
    def integral(self):return self.a.denominator==self.b.denominator==1
    def data(self):return dict(a=str(self.a),b=str(self.b))

PHI=K(0,1);BETA=4-PHI

def matrix_data(M):return [[K.of(x).data() for x in row] for row in M]
def inverse(M):
    n=len(M);A=[[Q(x) for x in row]+[Q(i==j) for j in range(n)] for i,row in enumerate(M)]
    for j in range(n):
        k=next(i for i in range(j,n) if A[i][j]);A[j],A[k]=A[k],A[j]
        v=A[j][j];A[j]=[x/v for x in A[j]]
        for i in range(n):
            if i!=j:
                v=A[i][j];A[i]=[x-v*y for x,y in zip(A[i],A[j])]
    return [row[n:] for row in A]
def det(M):
    A=[list(map(Q,row)) for row in M];d=Q(1);n=len(A)
    for j in range(n):
        k=next(i for i in range(j,n) if A[i][j])
        if k!=j:A[j],A[k]=A[k],A[j];d=-d
        p=A[j][j];d*=p
        for i in range(j+1,n):
            v=A[i][j]/p
            for l in range(j+1,n):A[i][l]-=v*A[j][l]
    return d
def dot(x,y):return sum((a*b for a,b in zip(x,y)),Q())
def unpack(x):return K(x['a'],x['b'])

def e8_and_family():
    roots=[[Q(1,2),Q(-1,2),Q(-1,2),Q(-1,2),Q(-1,2),Q(-1,2),Q(-1,2),Q(1,2)],
           [1,1,0,0,0,0,0,0],[-1,1,0,0,0,0,0,0]]
    for i in range(1,6):
        v=[0]*8;v[i]=-1;v[i+1]=1;roots.append(v)
    C=[[dot(x,y) for y in roots] for x in roots]
    edges={(0,2),(1,3),(2,3),(3,4),(4,5),(5,6),(6,7)}
    assert C==[[Q(2) if i==j else Q(-1) if tuple(sorted((i,j))) in edges else Q(0) for j in range(8)] for i in range(8)]
    assert abs(det(roots))==1 and det(C)==1
    J=inverse(C);assert all(x.denominator==1 for row in J for x in row)
    inputs=['results/fresh_h9_beta11_vinberg_20261014.json',
            'results/fresh_h9_beta55_pilot_20261014.json',
            'results/fresh_h9_beta131_pilot_20261016.json']
    rows=[]
    for m,path in enumerate(inputs,1):
        old=json.loads((ROOT/path).read_text());G=[[unpack(x) for x in row] for row in old['form']]
        expected=[[K(0) for _ in range(10)] for _ in range(10)]
        for i in range(8):
            for j in range(8):expected[i+1][j+1]=K(C[i][j])
        expected[0][0]=-2*PHI;expected[0][9]=expected[9][0]=K(-1);expected[9][9]=K(2*m)
        assert G==expected
        delta=1+4*m*PHI
        assert delta.norm()==1+4*m-16*m*m
        rows.append(dict(m=m,source=path,delta=delta.data(),determinant=(-delta).data(),determinant_norm=str(delta.norm())))
    assert [r['determinant_norm'] for r in rows]==['-11','-55','-131']
    silver_path='results/fresh_frontier_construction_20261015.json'
    silver=json.loads((ROOT/silver_path).read_text())
    # Identical coefficient-pair layout, but this file uses a+b*sqrt2.
    assert silver['form']==json.loads((ROOT/inputs[0]).read_text())['form']
    assert 1-2*4*4==-31
    inputs.append(silver_path)
    # Modulo2, phi has order3 and all three nonzero residues are units.
    residue=[]
    units=[K(1),PHI,PHI*PHI]
    for a,b in ((1,0),(0,1),(1,1)):
        c=K(a,b)
        u=next(u for u in units if (u*c+1).a%2==0 and (u*c+1).b%2==0)
        k=(u*c+1)/(2*PHI)
        assert k.integral() and u*c-2*PHI*k==-1 and abs(u.norm())==1
        residue.append(dict(c=c.data(),unit=u.data(),shift=k.data()))
    return dict(E8_euclidean_simple_roots=[[str(x) for x in r] for r in roots],
                E8_gram=[[str(x) for x in r] for r in C],
                E8_integral_inverse=[[str(x) for x in r] for r in J],
                old_family=rows,
                silver_analogue=dict(source=silver_path,field='Q(sqrt2)',binary=[['-2sqrt2','-1'],['-1','2']],determinant_norm=-31),
                nonzero_mod2_gluing_classes=residue,
                all_parameter_split='E8 unimodular: [[C,U],[U^T,B]] is integrally congruent to C orthogonal (B-U^T C^-1 U).',
                relative_fixed_source_normal_forms='E8 orthogonal [[-2phi,0],[0,2d]] or [[-2phi,-1],[-1,2d]], d in O; three tested integer d values do not exhaust this family.'),inputs

def local_obstruction():
    assert PHI*PHI==PHI+1 and BETA.norm()==11
    assert BETA.sign()>0 and BETA.conjugate().sign()>0
    squares=sorted({a*a%11 for a in range(11)})
    assert 10 not in squares and 4 in squares
    zeros=[(a,b) for a,b in product(range(11),repeat=2) if (a*a+b*b)%11==0]
    assert zeros==[(0,0)]
    assert (4*4-4-1)%11==0 and (2*4-1)%11!=0
    assert (2*2-4)%11==0 and (2*2)%11!=0
    iso=[1,1,3,0];assert sum(x*x for x in iso)%11==0 and 2*iso[0]%11!=0
    return dict(prime_generator=BETA.data(),prime_norm=11,residue_phi=4,
                residue_squares=squares,sum_two_squares_zeroes=zeros,
                phi_square_Hensel_seed=2,I4_isotropic_Hensel_seed=iso,
                candidate_local_Witt_index=3,E8_plus_binary_local_Witt_index_at_least=4,
                anisotropic_quaternary='I2 orthogonal beta I2; mod-beta descent forces every coordinate divisible beta.',
                maximal_integral='Only beta can divide a proper overlattice index. On beta^-1 I2, integrality requires a²+b²=0 mod beta, hence the vector was already integral.')

DIAG=[-PHI]+[K(1)]*7+[BETA]*2
def pair(x,y):return sum((g*a*b for g,a,b in zip(DIAG,x,y)),K())
def ordinary_roots(n):
    out=[]
    for i in range(n):
        for s in (-1,1):
            v=[0]*n;v[i]=s;out.append(v)
    for i,j in combinations(range(n),2):
        for s,t in product((-1,1),repeat=2):
            v=[0]*n;v[i],v[j]=s,t;out.append(v)
    assert len(out)==2*n*n
    return out

def sum_squares_coefficients(n):
    # A=sum(a_i²+b_i²)<=2 bounds each pair. Track the phi coefficient
    # B=sum(2a_i b_i+b_i²); the complete norm1/2 condition is B=0.
    alphabet=[(a,b) for a,b in product(range(-1,2),repeat=2) if a*a+b*b<=2]
    states={(0,0):1}
    for _ in range(n):
        nxt={}
        for (A,B),count in states.items():
            for a,b in alphabet:
                aa=A+a*a+b*b;bb=B+2*a*b+b*b
                if aa<=2:nxt[aa,bb]=nxt.get((aa,bb),0)+count
        states=nxt
    assert states[1,0]==2*n and states[2,0]==2*n*(n-1)
    return dict(rank=n,norm1=states[1,0],norm2=states[2,0])

def norm_and_zero():
    norms=[K(1),K(2),BETA,2*BETA]
    dual_exponent_entries=[BETA/g for g in DIAG]
    assert all(v.integral() for v in dual_exponent_entries)
    assert (Q(4,11)<1)
    roots=[]
    for v in ordinary_roots(7):roots.append([K(0)]+list(map(K,v))+[K(0)]*2)
    for v in ordinary_roots(2):roots.append([K(0)]*8+list(map(K,v)))
    assert len(roots)==106
    counts={str(x.data()):0 for x in norms}
    for r in roots:
        q=pair(r,r);assert any(q==v for v in norms)
        counts[str(q.data())]+=1
        assert all((2*g*x/q).integral() for g,x in zip(DIAG,r))
    simple=[]
    for i in range(1,7):
        r=[K(0)]*10;r[i],r[i+1]=K(-1),K(1);simple.append(r)
    r=[K(0)]*10;r[7]=K(-1);simple.append(r)
    r=[K(0)]*10;r[8],r[9]=K(-1),K(1);simple.append(r)
    r=[K(0)]*10;r[9]=K(-1);simple.append(r)
    G=[[pair(x,y) for y in simple] for x in simple]
    relations=[]
    for i,j in combinations(range(9),2):
        cyc=4*G[i][j]*G[i][j]/(G[i][i]*G[j][j])
        assert G[i][j].sign()<=0
        label=2 if cyc==0 else 3 if cyc==1 else 4 if cyc==2 else None
        assert label is not None
        relations.append(dict(pair=[i,j],label=label))
    assert [r['pair'] for r in relations if r['label']==4]==[[5,6],[7,8]]
    return dict(primitive_norm_classes=[x.data() for x in norms],
                beta_times_inverse_diagonal=[v.data() for v in dual_exponent_entries],
                complete_norm_argument='beta*G^-1 is integral; reflection integrality gives2beta*r/q integral. Primitive-coordinate Bezout therefore gives q|2beta. Positive units are squares, leaving1,2,beta,2beta.',
                zero_mixing_bound='Norm(nonzero coordinate)^2<=4/11<1 after the beta residue condition, impossible.',
                ordinary_coefficient_counts=[sum_squares_coefficients(n) for n in (7,2)],
                root_count_convention='Unit-normalized to norms1,2,beta,2beta; opposite directions counted separately.',
                zero_root_count=106,zero_mirror_count=53,zero_type='B7+B2',
                zero_roots=matrix_data(roots),zero_norm_counts=counts,
                zero_simple_walls=matrix_data(simple),zero_gram=matrix_data(G),zero_relations=relations,
                chamber='y1>=y2>=...>=y7>=0; z1>=z2>=0')

def run():
    family,inputs=e8_and_family()
    determinant=-PHI*BETA*BETA
    assert determinant.norm()==-121
    assert (determinant/(-PHI))==BETA*BETA
    assert DIAG[0].sign()<0 and all(x.sign()>0 for x in DIAG[1:])
    assert all(x.conjugate().sign()>0 for x in DIAG)
    src=f'scripts/{TAG}.py'
    refs=[src]+inputs+['results/fresh_h9_h8_cap_height_20261014_norm_pool_review.json']
    return dict(status='PASS_DISTINCT_LOCAL_WITT_H9_HOST_COMPLETE_NORMS_AND_ZERO_SECTOR',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in refs},
                field='Q(phi), phi²=phi+1, physical phi>1',coordinates=['x']+[f'y{i}' for i in range(1,8)]+['z1','z2'],
                form=matrix_data([[DIAG[i] if i==j else K(0) for j in range(10)] for i in range(10)]),
                basepoint=[K(int(i==0)).data() for i in range(10)],
                determinant=determinant.data(),determinant_norm='-121',
                signatures=[[9,1],[10,0]],family_comparison=family,
                local_obstruction=local_obstruction(),norm_and_zero=norm_and_zero(),
                no_reflectivity_claim=True,
                scope='New exact admissible maximal-integral H9 lattice, not similar to ANY E8-plus-binary form over Qsqrt5. Complete primitive norm alphabet and zero sector. No positive-height search, finite-volume example, or universal nonreflectivity claim for the earlier one-parameter family.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text());print('PASS fresh saved local-Witt host replay')
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['norm_and_zero']['zero_type'],out['norm_and_zero']['zero_root_count'])
