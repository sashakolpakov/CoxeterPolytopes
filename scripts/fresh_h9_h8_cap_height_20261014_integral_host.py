#!/usr/bin/env python3
"""An explicit integral reflection module for the maximal-height chamber."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as s
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_cap_height_20261014_integral_host'
BASE='results/fresh_h9_h8_cap_height_20261014.json'
REVIEW='results/fresh_h9_h8_cap_height_20261014_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def integer_pair(x):
    assert not x.b
    # a+b sqrt5=(a-b)+(2b)phi, O_K=Z[phi].
    u,v=x.a.a-x.a.b,2*x.a.b
    assert u.denominator==v.denominator==1
    return[int(u),int(v)]
def run():
    d,r=[json.loads((ROOT/p).read_text())for p in(BASE,REVIEW)]
    for obj in(d,r):
        for p,h in obj['sources'].items():assert sha(p)==h,p
    Q,N=A.dm(d['form']),A.dm(d['normals'])
    N[11]=[2*x for x in N[11]]
    G=[[A.pair(n,m,Q)for m in N]for n in N];norms=[G[i][i]for i in range(15)]
    assert norms[:11]==[A.E(2)]*11 and norms[12:]==[A.E(2)]*3
    assert norms[11]==A.E(A.Q(7,1))
    gram_integer=[[integer_pair(x)for x in row]for row in G]
    cartan=[[2*G[i][j]/norms[j]for j in range(15)]for i in range(15)]
    cartan_integer=[[integer_pair(x)for x in row]for row in cartan]
    def expr(x):
        assert not x.b
        return s.Rational(x.a.a.numerator,x.a.a.denominator)+s.Rational(x.a.b.numerator,x.a.b.denominator)*s.sqrt(5)
    K=s.QQ.algebraic_field(s.sqrt(5));M=s.polys.matrices.DomainMatrix.from_Matrix(s.Matrix([[expr(x)for x in row]for row in N])).convert_to(K)
    assert M.rank()==10
    # Independently verify every ring entry from exact symbolic products.
    Qs=s.Matrix([[expr(x)for x in row]for row in Q]);Ns=s.Matrix([[expr(x)for x in row]for row in N])
    Gs=Ns*Qs*Ns.T
    phi=(1+s.sqrt(5))/2
    for i in range(15):
        for j in range(15):
            u,v=gram_integer[i][j];assert s.simplify(Gs[i,j]-u-v*phi)==0
            u,v=cartan_integer[i][j];assert s.simplify(2*Gs[i,j]/Gs[j,j]-u-v*phi)==0
    assert A.inertia(Q)==[9,1]
    sigma=[[A.E(A.Q(x.a.a,-x.a.b))for x in row]for row in Q]
    assert A.inertia(sigma)==[10,0]
    assert cartan[11][12]==A.E(-1)*A.E(A.Q(F(7,2),F(1,2)))
    assert cartan[12][11]==-1
    paths=[BASE,REVIEW,f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
    return A.enc(dict(status='PASS_MAXHEIGHT_INTEGRAL_REFLECTION_MODULE',sources={p:sha(p)for p in paths},
        ring='Z[phi],phi=(1+sqrt5)/2; integer pairs encode u+v*phi',
        scaled_bottom_index=11,normals=N,gram_integer_pairs=gram_integer,cartan_integer_pairs=cartan_integer,
        generator_count=15,field_rank=10,signature=[9,1],conjugate_signature=[10,0],
        proof='Let L be the Z[phi]-span of the15scaled normals. It has full field rank10 and integral Gram. Every reflection sends generator r_i to r_i−Cartan[i,j]*r_j, so every generator reflection preserves L. The bottom scaling changes no hyperplane. Separate symbolic products verify all225Gram and225Cartan ring entries.',
        scope='An explicit crystallographic reflection module in an admissible quadratic form, suitable for a new arithmetic root search. No claim that the15walls form the full reflection chamber, that the lattice is reflective, or that its reflection subgroup has finite index. The chamber still has the certified strict escape. No unimodular classification is applied.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['field_rank'])
