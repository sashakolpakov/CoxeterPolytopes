#!/usr/bin/env python3
"""Separate Fraction and 128-bit review of the low-polar attachment gate."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
from fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate import I, ei, SQ5, BITS

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_h8_corner_anchors_20261014_polar_review'
PROOF = 'results/fresh_h9_h8_corner_anchors_20261014_polar_gate.json'
PAIR = 'results/fresh_h9_h8_corner_anchors_20261014_pair_review.json'
BASE = 'results/fresh_h9_mixed_d4_coupling_20261014.json'

def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def p_add(*terms):
    out = {}
    for p in terms:
        for k, v in p.items(): out[k] = out.get(k, F(0))+v
    return {k: v for k, v in out.items() if v}
def p_mul(a, b):
    out = {}
    for k, v in a.items():
        for l, w in b.items():
            m = tuple(x+y for x, y in zip(k, l)); out[m] = out.get(m, F(0))+v*w
    return {k: v for k, v in out.items() if v}
def p_scale(a, x): return {k: v*x for k, v in a.items() if v*x}
def p_diff(a, i):
    return {tuple(n-(j == i) for j, n in enumerate(k)): v*k[i] for k, v in a.items() if k[i]}

def run():
    pending, seen = [PROOF, PAIR], set()
    while pending:
        p = pending.pop()
        if p in seen: continue
        seen.add(p); d = json.loads((ROOT/p).read_text())
        for path, h in d.get('sources', {}).items():
            assert sha(path) == h, path
            if path.endswith('.json'): pending.append(path)
    proof = json.loads((ROOT/PROOF).read_text())
    pair_review = json.loads((ROOT/PAIR).read_text())
    base = json.loads((ROOT/BASE).read_text())
    Q, N = A.dm(base['form']), A.dm(base['final_normals'])
    norm = [A.pair(n, n, Q) for n in N]
    assert all(q == 2 for q in norm[:11]) and norm[11].sign() > 0
    zeros = [(0, 2), (0, 4), (2, 4), (0, 11), (2, 11)]
    assert all(not A.pair(N[i], N[j], Q) for i, j in zeros)
    rights = [i for i in range(11) if not A.pair(N[11], N[i], Q)]
    assert rights == [0, 1, 2, 3, 5, 6, 7, 9, 10]
    assert all(A.pair(N[11], N[i], Q).sign() < 0 for i in (4, 8))
    C = -2*ei(A.pair(N[4], N[11], Q))/ei(norm[4]*norm[11]).sqrt()
    a, B = (SQ5-1)/4, SQ5-1
    L2, DA = 2*(B+2), 9+4*SQ5
    margin = dict(C_lower=C-F(44, 5), C_upper=9-C,
                  high_contact=L2-6, low_contact=F(169, 25)-L2)
    gram = A.dm(pair_review['source_exceptional_gram'])
    ae = A.E(A.Q(F(-1, 4), F(1, 4)))
    dae = A.E(A.Q(9, 4)); dbe = A.E(A.Q(17, 8))
    assert all(gram[i][i] == -ae for i in range(3))
    assert -gram[0][1]/ae == dae == -gram[0][2]/ae
    assert -gram[1][2]/ae == dbe and (dbe-dae).sign() > 0
    def reach(x, y, c):
        d = (((x+2)*(y+2)).sqrt()+c)/(x*y).sqrt()
        assert (d-1).lo > 0
        return d*(1+2/y).sqrt()+(d.square()-1).sqrt()*(2/y).sqrt()
    comparisons = dict(right_pair=reach(a, B, I(0)), own_ball=(1+2/a).sqrt(),
                       high_ball_diameter=1+4/B, wall4_if_contact_at_least4=reach(a, I(6), I(9)))
    margin.update({key: DA-value for key, value in comparisons.items()})
    assert all(value.lo > 0 for value in margin.values())
    # Independent polynomial arithmetic verifies the Lorentz-block inverse
    # numerator and its derivative. Variable order is A,D,C.
    pa, pd, pc = [{tuple(int(i == j) for i in range(3)): F(1)} for j in range(3)]
    a2, d2, c2 = [p_mul(p, p) for p in (pa, pd, pc)]
    ad = p_mul(pa, pd)
    numerator = p_scale(p_add(a2, d2, p_mul(ad, pc)), 2)
    denominator = p_add(c2, {(0, 0, 0): F(-4)})
    derivative_numerator = p_add(p_mul(p_diff(numerator, 2), denominator),
                                 p_scale(p_mul(numerator, p_diff(denominator, 2)), -1))
    expected = p_add(p_scale(p_mul(ad, p_add(c2, {(0, 0, 0): F(4)})), -2),
                     p_scale(p_mul(pc, p_add(a2, d2)), -4))
    assert derivative_numerator == expected
    assert all(v < 0 for v in expected.values())
    assert all(v > 0 for v in p_diff(numerator, 0).values())
    assert all(v > 0 for v in p_diff(numerator, 1).values())
    # Multiplication by [[2,C],[C,2]] gives (4-C²) times the identity.
    diag = p_add({(0, 0, 0): F(4)}, p_scale(c2, -1))
    assert p_add(diag, denominator) == {}
    assert p_add(p_scale(pc, 2), p_scale(pc, -2)) == {}
    upper_A, upper_D, lower_C = F(13, 5), F(4), F(44, 5)
    penalty = 2*(upper_A**2+upper_D**2+lower_C*upper_A*upper_D)/(lower_C**2-4)
    lower = 6-penalty
    assert penalty == F(2857, 918) and lower == F(2651, 918) > 2
    assert str(lower) == proof['rational_total_norm_lower']
    paths = [PROOF, PAIR, BASE, f'scripts/{TAG}.py',
             'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
             'scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py']
    return dict(status='PASS_SEPARATE_LOW_POLAR_PROJECTION_REVIEW',
                sources={p: sha(p) for p in paths}, interval_bits=BITS,
                strict_margins={k: v.enc() for k, v in margin.items()},
                exact_zero_pairs=[list(p) for p in zeros], rational_norm_lower=str(lower),
                scope='No polar-gate producer, SymPy or its192-bit evaluator imported. Separate128-bit interval bounds, exact literal pairings and sparse-polynomial identities reconstruct the four-root projection contradiction. Reviewed source-gap and at-most-one-low-center results remain dependencies. Only raw twelve-old-ball coverage is excluded; arbitrary additional nonperpendicular walls remain outside scope.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--verify', action='store_true'); args = parser.parse_args()
    result = run(); out = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['rational_norm_lower'], '> 2')
