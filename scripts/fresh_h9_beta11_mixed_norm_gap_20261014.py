#!/usr/bin/env python3
"""Exact inputs to orthogonality/disjointness of the two root norm classes."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_beta11_mixed_norm_gap_20261014'
POOL='results/fresh_h9_h8_cap_height_20261014_norm_pool_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    d=json.loads((ROOT/POOL).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    beta=A.Q(F(7,2),F(1,2));sigma=A.Q(F(7,2),F(-1,2))
    assert beta*sigma==11 and beta+sigma==7 and (beta-4).sign()>0 and sigma.sign()>0
    traces=[t for t in range(1,8)if t*t>=44];assert traces==[7]
    assert 7*7-44==5
    norm_choices=[11*k*k for k in range(1,3)if 11*k*k<16];assert norm_choices==[11]
    paths=[POOL,f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
    return dict(status='PASS_MIXED_NORM_ROOT_HYPERPLANES_ORTHOGONAL_OR_DISJOINT',
        sources={p:sha(p)for p in paths},beta=beta.enc(),conjugate_beta=sigma.enc(),
        hypothetical_norm_choices=norm_choices,hypothetical_trace_choices=traces,
        proof='For norm2 r andnorm2beta s, longroot reflectionintegrality givesg=(r,s)=beta*z withz inO. Ifnonzero, a=g²/beta=beta*z² istotallypositive, Norm(a)=11Norm(z)². Theconjugateform ispositive definite, so0<a_sigma<4. Intersectionorparallelism atthereal hyperbolicembedding wouldgive0<a<=4. Thus0<Norm(a)<16 forcesNorm(a)=11 andintegerTrace(a)<8. Nonnegative quadraticdiscriminant forcestrace7, hencea=(7±sqrt5)/2. Oneembeddingthenexceeds4, contradiction. Normbeta11nonsquare excludescoincidentrootlines.',
        scope='All integral reflection directions in the two normalized norm classes, not only searched roots or bounded heights. Mixed classes have no nonright finite or parallel contacts. The finite exact inputs are checked here; reflection integrality and conjugate positive-definiteness are inherited from the independent lattice review. No reflectivity or finitevolume conclusion.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
