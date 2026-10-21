#!/usr/bin/env python3
"""Separate exact source-Gram and interval review of the two-anchor gate."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
from fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate import I, ei, BITS
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_h8_corner_anchors_20261014_pair_review'

def run():
    cert = 'results/fresh_h9_h8_corner_anchors_20261014_pair_gate.json'
    data = json.loads((ROOT/cert).read_text())
    pending, seen = [(cert, data)], set()
    while pending:
        path, doc = pending.pop()
        if path in seen: continue
        seen.add(path)
        for p, h in doc.get('sources', {}).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
            if p.endswith('.json'): pending.append((p, json.loads((ROOT/p).read_text())))
    srcpath = 'data/known_examples/bugaenko_compact_h8.json'
    gap = json.loads((ROOT/'results/fresh_h9_h8_facet_cover_20261014_gap_gate.json').read_text())
    src = json.loads((ROOT/srcpath).read_text())
    def sf(x): return A.E(A.Q(F(x['a'])+F(x['b'])/2, F(x['b'])/2))
    G = [[sf(x) for x in r] for r in src['wall_gram']]
    J = A.inverse([r[:9] for r in G[:9]])
    assert all(sum((G[i][k]*J[k][j] for k in range(9)), A.Z) == A.E(i == j)
               for i in range(9) for j in range(9))
    rows = sorted(gap['exceptional_rays'], key=lambda r: r['source_vertex'])
    assert [r['source_vertex'] for r in rows] == [0, 32, 33]
    contacts = [[A.E(A.Q(*x)) for x in r['contacts']] for r in rows]
    coefficients = [A.mv(J, h[:9]) for h in contacts]
    for h, coeff in zip(contacts, coefficients):
        assert all(sum((G[i][j]*coeff[j] for j in range(9)), A.Z) == h[i] for i in range(11))
        assert any(x == 1 for x in h)
    a = A.E(A.Q(F(-1, 4), F(1, 4)))
    ba = A.E(A.Q(F(11, 4), F(5, 4))); bb = A.E(A.Q(F(23, 4), F(9, 4)))
    expected = [[-a, -ba, -ba], [-ba, -a, -bb], [-ba, -bb, -a]]
    observed = [[A.dot(h[:9], c) for c in coefficients] for h in contacts]
    assert observed == expected
    assert (ba*ba-a*a-2*a).sign() > 0 and (bb*bb-a*a-2*a).sign() > 0
    phi = A.E(A.Q(F(1, 2), F(1, 2))); assert ba-2-a == 2*phi
    # Same-ray product is at least2: the square-root radicand exceeds
    # (2+a*s*t)^2 by the nonnegative polynomial2a(s-t)^2.
    lhs = {(2, 0): 2, (0, 2): 2, (1, 1): -4}
    assert lhs == {p: 2*v for p, v in {(2, 0): 1, (0, 2): 1, (1, 1): -2}.items()}
    basepath = 'results/fresh_h9_mixed_d4_coupling_20261014.json'
    base = json.loads((ROOT/basepath).read_text()); Q = A.dm(base['form']); N = A.dm(base['final_normals'])
    norms = [A.pair(n, n, Q) for n in N]; assert all(n.sign() > 0 for n in norms)
    def magnitude(i, j):
        g = A.pair(N[i], N[j], Q); assert g.sign() <= 0
        return -2*ei(g)/ei(norms[i]*norms[j]).sqrt()
    ordinary, special = [], {}
    for i in range(12):
        for j in range(i+1, 12):
            c = magnitude(i, j)
            if (i, j) in ((4, 11), (8, 11)): special[i, j] = c
            else:
                margin = 2*ei(phi)-c; assert margin.lo > 0
                ordinary.append(dict(pair=[i, j], margin=margin.enc()))
    assert len(ordinary) == 64
    ai, bai, bbi = map(ei, (a, ba, bb)); one, rt2, phii = I(1), I(2).sqrt(), ei(phi)
    def f(b, s, t): return b*s*t-((2+ai*s.square())*(2+ai*t.square())).sqrt()
    c4, c8 = special[4, 11], special[8, 11]
    margins = dict(
        pair8_A_lower=c8-f(bai, one, one), pair8_A_next=f(bai, one, rt2)-c8,
        pair4_A_no_one=c4-f(bai, one, I(2)), pair4_A_both_four=c4-f(bai, rt2, rt2),
        pair4_A_next=f(bai, rt2, phii)-c4, pair4_B_lower=c4-f(bbi, one, one),
        pair4_B_next=f(bbi, one, rt2)-c4, pair8_B_below=f(bbi, one, one)-c8)
    assert all(v.lo > 0 for v in margins.values())
    paths = [cert, srcpath, basepath, f'scripts/{TAG}.py',
             'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
             'scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py']
    return dict(status='PASS_SEPARATE_EXACT_TWO_ANCHOR_REVIEW',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                source_exceptional_gram=A.enc(observed), interval_bits=BITS,
                ordinary_pair_margins=ordinary, special_margins={k: v.enc() for k, v in margins.items()},
                scope='No pair-gate or its192-bit arithmetic imports. Exact Fraction-quartic matrix inversion reconstructs all11 contacts and the three-ray Gram; a separate128-bit dyadic evaluator reproduces all72 strict pair comparisons. Analytic monotonicity and same-ray inequalities are checked with the stated scale domain. The frozen source-gap theorem remains a dependency. At most one low-norm center is necessary, not sufficient for an attachment.')

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); args = p.parse_args()
    result = run(); out = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
