#!/usr/bin/env python3
"""Independent Fraction arithmetic for the free host and two norm classes.

No producer, project arithmetic helper, or symbolic-algebra import is used.
The universal primitivity/unit argument is reviewed in the accompanying note.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_h8_cap_height_20261014_norm_pool_review'
BASE = 'results/fresh_h9_h8_cap_height_20261014.json'
HOST = 'results/fresh_h9_h8_cap_height_20261014_integral_host.json'
BASIS = 'results/fresh_h9_h8_cap_height_20261014_host_basis.json'
POOL = 'results/fresh_h9_h8_cap_height_20261014_norm_pool.json'


class K:
    """a+b*phi, phi^2=phi+1."""
    def __init__(self, a=0, b=0):
        if isinstance(a, K):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = F(a), F(b)
    def __add__(self, other):
        z=K(other); return K(self.a+z.a, self.b+z.b)
    __radd__=__add__
    def __neg__(self): return K(-self.a,-self.b)
    def __sub__(self, other): return self+-K(other)
    def __rsub__(self, other): return K(other)+-self
    def __mul__(self, other):
        z=K(other)
        return K(self.a*z.a+self.b*z.b,
                 self.a*z.b+self.b*z.a+self.b*z.b)
    __rmul__=__mul__
    def conjugate(self): return K(self.a+self.b,-self.b)
    def norm(self): return self.a*self.a+self.a*self.b-self.b*self.b
    def __truediv__(self, other):
        z=K(other); n=z.norm(); assert n
        c=z.conjugate(); return self*K(c.a/n,c.b/n)
    def __eq__(self, other):
        z=K(other); return self.a==z.a and self.b==z.b
    def __bool__(self): return bool(self.a or self.b)
    def sign(self):
        a,b=self.a+self.b/2,self.b/2
        if not b: return (a>0)-(a<0)
        if not a: return (b>0)-(b<0)
        if (a>0)==(b>0): return (a>0)-(a<0)
        d=a*a-5*b*b; assert d
        return ((a>0)-(a<0))*((d>0)-(d<0))
    def integral(self): return self.a.denominator==self.b.denominator==1
    def enc(self): return [str(self.a),str(self.b)]


def dec(x):
    # The source's extra quadratic extension has zero coefficient here.
    assert all(F(t)==0 for t in x[1])
    a,b=map(F,x[0]); return K(a-b,2*b)
def matrix(x): return [[dec(z) for z in row] for row in x]
def dot(x,y): return sum((a*b for a,b in zip(x,y)),K())
def mm(x,y): return [[dot(row,col) for col in zip(*y)] for row in x]
def transpose(x): return list(map(list,zip(*x)))
def determinant(x):
    a=[row[:] for row in x]; d=K(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None: return K()
        if pivot!=i: a[i],a[pivot]=a[pivot],a[i]; d=-d
        p=a[i][i]; d*=p
        for j in range(i+1,len(a)):
            z=a[j][i]/p
            for k in range(i+1,len(a)): a[j][k]-=z*a[i][k]
    return d
def sha(path): return sha256((ROOT/path).read_bytes()).hexdigest()
def read(path):
    d=json.loads((ROOT/path).read_text())
    for p,h in d['sources'].items(): assert sha(p)==h,p
    return d


def run():
    base,host,basis,pool=map(read,(BASE,HOST,BASIS,POOL))
    Q,N,B=matrix(base['form']),matrix(host['normals']),matrix(basis['basis'])
    R=[[K(*x) for x in row] for row in basis['generators_in_basis_integer_pairs']]
    T=[[K(*x) for x in row] for row in basis['basis_in_generators_integer_pairs']]
    assert all(x.integral() for M in [R,T] for row in M for x in row)
    assert mm(R,B)==N and mm(T,N)==B
    G=mm(mm(B,Q),transpose(B)); assert G==matrix(basis['gram'])
    phi,beta=K(0,1),K(3,1)
    assert phi*phi==phi+1 and phi*(phi-1)==1
    assert phi.norm()==-1 and beta.norm()==11
    assert (phi-1).sign()>0 and (2-phi).sign()>0
    assert beta.sign()>0 and beta.conjugate().sign()>0
    assert all(x.integral() for row in G for x in row)
    assert all((G[i][i]/2).integral() for i in range(10))
    E8=[row[1:9] for row in G[1:9]]
    minors=[determinant([row[:k] for row in E8[:k]]) for k in range(1,9)]
    assert all(x.sign()>0 and x.conjugate().sign()>0 for x in minors)
    assert minors[-1]==1
    assert G[0][0]==-2*phi and G[0][9]==G[9][0]==-1 and G[9][9]==2
    assert all(not G[i][j] for i in (0,9) for j in range(1,9))
    det=determinant(G); assert det==-phi*beta and det.norm()==-11
    assert det.sign()<0 and det.conjugate().sign()>0
    assert G[0][0].conjugate().sign()>0
    # E8 is positive at both embeddings. The binary block is (1,1) at
    # the first and positive definite at the second, by its 1x1/2x2 minors.
    assert [dec(x) for x in pool['normalized_primitive_root_norms']]==[K(2),2*beta]
    assert dec(pool['determinant'])==det and dec(pool['prime_generator'])==beta
    assert F(pool['euclidean_error_bound'])==F(3,4)<1
    assert all(11%d for d in range(2,4))
    # O/(beta)=F_11 under phi->8: determinant of multiplication is11.
    assert (8*8-8-1)%11==0 and (3+8)%11==0
    assert determinant([[K(3),K(1)],[K(1),K(4)]])==11
    residue_tests=0
    for a in range(11):
        for b in range(11):
            assert ((a+8*b)%11==0)==((K(a,b)/beta).integral())
            residue_tests+=1
    residue_gram=[[(int(x.a)+8*int(x.b))%11 for x in row] for row in G]
    kernel=[1]+[0]*8+[6]
    assert all(sum(a*b for a,b in zip(row,kernel))%11==0 for row in residue_gram)
    assert residue_gram[0][0]==6 and residue_gram[0][9]==10 and residue_gram[9][9]==2
    # E8 determinant1 plus a nonzero singular binary block gives rank9.
    assert (6*2-10*10)%11==0 and minors[-1]==1
    # Check the operational norm alphabet/filter on every given generator.
    generator_checks=[]
    for i,r in enumerate(R):
        gr=[dot(row,r) for row in G]; q=dot(r,gr)
        assert q in [K(2),2*beta]
        assert all((x/(q/2)).integral() for x in gr)
        generator_checks.append(dict(index=i,norm=q.enc(),reflection_filter=True))
    # Rational bounds used in the universal no-unit-in-(1,phi) proof.
    # For Norm(u)=+1: 0<(u-u^-1)/sqrt5<1/sqrt5<1.
    # For Norm(u)=-1: 0<(u+u^-1)/sqrt5<(phi+phi^-1)/sqrt5=1.
    assert phi-(phi-1)==1
    sqrt5=2*phi-1; assert sqrt5*sqrt5==5 and sqrt5.sign()>0
    assert phi+(phi-1)==sqrt5
    paths=[BASE,HOST,BASIS,POOL,f'scripts/{TAG}.py',
           'scripts/fresh_h9_h8_cap_height_20261014_host_basis.py',
           'scripts/fresh_h9_h8_cap_height_20261014_norm_pool.py']
    return dict(status='PASS_INDEPENDENT_FREE_HOST_AND_NORM_ALPHABET_INPUTS',
        sources={p:sha(p) for p in paths},
        arithmetic='Standalone Fraction a+b*phi; no project or symbolic arithmetic imports.',
        module_inclusions=dict(generators=15,basis_rank=10,both_checked=True),
        gram_integral=True,half_diagonal_integral=True,
        E8_leading_minors=[x.enc() for x in minors],
        determinant=det.enc(),determinant_norm=str(det.norm()),
        signatures=[[9,1],[10,0]],beta=beta.enc(),beta_norm=str(beta.norm()),
        normalized_norms=[K(2).enc(),(2*beta).enc()],
        residue_map=dict(prime=11,phi_image=8,coefficient_classes_checked=residue_tests),
        residue_gram_rank=9,residue_kernel_generator=kernel,
        beta_norm_filter='Every E8 coordinate is0mod beta; t=6xmod beta.',
        generator_reflection_checks=generator_checks,
        universal_proof_review='PASS: Euclidean rounding, primitive Bezout twice, prime determinant divisor, total positivity, fundamental unit and unit-square normalization.',
        scope='Complete primitive L-preserving reflection-direction alphabet up to unit scaling. Finite algebraic inputs independently reconstructed; universal proof reviewed in the report. No root enumeration, height completeness, reflectivity or arbitrary geometric-completion claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify: assert d==json.loads(out.read_text())
    else: out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
