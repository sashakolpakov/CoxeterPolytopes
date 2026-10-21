#!/usr/bin/env python3
"""Finite exact inputs to the complete primitive-reflector norm alphabet."""
import argparse,json
from hashlib import sha256
from pathlib import Path
import sympy as s
import fresh_h9_mixed_d4_coupling_20261014_review as A
from fresh_h9_h8_cap_height_20261014_integral_host import integer_pair
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_cap_height_20261014_norm_pool'
BASIS='results/fresh_h9_h8_cap_height_20261014_host_basis.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    d=json.loads((ROOT/BASIS).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    G=A.dm(d['gram']);phi=A.PHI;beta=3+phi
    even_diagonal=[integer_pair(G[i][i]/2)for i in range(10)]
    integral_gram=[[integer_pair(x)for x in row]for row in G]
    det=A.E(A.Q(*d['determinant']))
    assert det==-phi*beta and d['determinant_field_norm']==-11
    assert beta.sign()>0 and A.Q(beta.a.a,-beta.a.b).sign()>0
    assert beta.a.a*beta.a.a-5*beta.a.b*beta.a.b==11 and s.isprime(11)
    assert phi*phi==phi+1 and phi*(phi-1)==1
    assert phi.a.a*phi.a.a-5*phi.a.b*phi.a.b==-1
    assert (phi-1).sign()>0 and (2-phi).sign()>0
    # Euclidean error bound: |a²+ab−b²| <= 1/4+1/4+1/4 < 1.
    euclidean_bound=s.Rational(3,4);assert euclidean_bound<1
    paths=[BASIS,f'scripts/{TAG}.py','scripts/fresh_h9_h8_cap_height_20261014_integral_host.py',
           'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
    return A.enc(dict(status='PASS_EXACT_INPUTS_TO_TWO_PRIMITIVE_ROOT_NORM_CLASSES',sources={p:sha(p)for p in paths},
        half_diagonal_integer_pairs=even_diagonal,integral_gram=integral_gram,
        determinant=det,prime_generator=beta,prime_generator_norm=11,
        normalized_primitive_root_norms=[A.E(2),2*beta],euclidean_error_bound='3/4',
        proof='O is normEuclidean by rounding its1,phi coefficients, so every reflection line has a primitive vector. Evenness givesq=2b. Reflectionintegrality givesGr in bO^10; adjugate and primitiveBezout imply b dividesdetG=−phi*beta. Normbeta11 isprime. Both q andbeta aretotallypositive. Units are±phi^n: no positive unit lies strictlybetween1 andphi, using the integer coefficient(u−u_sigma)/sqrt5 and Normu=±1. Totallypositive units arephi^(2k). Unit rescaling therefore givesq=2 or2beta.',
        scope='Complete norm alphabet for primitive L-preserving reflection directions, up to unit rescaling. The finite algebraic inputs are checked here; the adjugate, Euclidean and unit arguments are explicit mathematical deductions in the report. No height bound, root enumeration, wall completion or reflectivity claim.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
