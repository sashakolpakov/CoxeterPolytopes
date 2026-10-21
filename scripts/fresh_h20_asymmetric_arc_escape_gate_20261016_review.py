#!/usr/bin/env python3
"""Independent symbolic review of the abelian-base upgrade and weak escape."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_asymmetric_arc_escape_gate_20261016_review'
PRIMARY = 'results/fresh_h20_asymmetric_arc_escape_gate_20261016.json'
CHAIN = 'results/fresh_h20_s61_multicap_20261014_exact.json'
CAP = 'results/fresh_h20_dotted_arc_cap_20261016.json'
CAP_REVIEW = 'results/fresh_h20_dotted_arc_cap_20261016_review.json'
ABEL_REVIEW = 'results/fresh_h20_chain_abelian_gate_20261016_review.json'


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def read(p):
    d = json.loads((ROOT / p).read_text())
    for q, expected in d.get('sources', {}).items():
        assert sha(q) == expected, (p, q)
    return d


def rat(v):
    return S.Rational(str(v))


def radical(v):
    return rat(v[0]) + rat(v[1]) * S.sqrt(2)


def positive(v):
    return S.simplify(v).is_positive is True


def run():
    primary, chain, cap, caprev, abelrev = map(
        read, (PRIMARY, CHAIN, CAP, CAP_REVIEW, ABEL_REVIEW))
    assert caprev['sources'][CAP] == sha(CAP)
    factors = {r['prime']: r for r in abelrev['modular_factorizations']}
    assert factors[2017]['factor_degrees'] == [50]
    assert factors[359]['factor_degrees'] == [13, 37]
    assert abelrev['independently_reconstructed_endpoint_polynomial_degree'] == 50
    for r in factors.values():
        assert all(r[k] for k in ('simple_tower_roots', 'leading_unit', 'squarefree'))
    # A separately labelled Frobenius representative: first a37-cycle, then13.
    sigma = list(range(1, 37)) + [0] + list(range(38, 50)) + [37]
    tau = []
    for i in range(50):
        j = i
        for _ in range(13):
            j = sigma[j]
        tau.append(j)
    fixed = sum(i == j for i, j in enumerate(tau))
    orbit, j = [], 0
    while j not in orbit:
        orbit.append(j)
        j = tau[j]
    assert fixed == 13 and len(orbit) == 37 and j == 0
    proper_blocks = [(b, 50 // b) for b in range(2, 50) if 50 % b == 0]
    assert proper_blocks == [(2, 25), (5, 10), (10, 5), (25, 2)]
    assert all(max(b, count) < 37 for b, count in proper_blocks)
    # The normal-subgroup and base-change conclusions are proved in the report.
    upgrade = primary['Galois_upgrade']
    assert upgrade['degree_over_K_Qab'] == 50
    assert upgrade['no_nonzero_polynomial_over_Qab_of_degree_below'] == 50

    c, A, B, D = S.symbols('c A B D', positive=True)
    triple = S.Matrix([[2, 0, -4*c-2], [0, 2, 0], [-4*c-2, 0, 2]])
    contacts = S.Matrix([-B, -A, -D])
    q = S.factor((contacts.T * triple.inv() * contacts)[0])
    formula = A*A/2 - ((B+D)**2 + 4*c*B*D)/(8*c*(c+1))
    assert S.simplify(q-formula) == 0
    assert S.factor(triple.det()) == -32*c*(c+1)
    Q = S.symbols('Q')
    relation = S.Poly(S.cancel(8*c*(c+1)*(q-2+Q)), c)
    assert relation.degree() == 2 and S.simplify(relation.TC()+(B+D)**2) == 0
    # At B=D=0 the projection is -A/2 times the middle root.
    assert triple.inv()*S.Matrix([0, -A, 0]) == S.Matrix([0, -A/2, 0])
    assert S.simplify((-A/2)*(-4*c-2)-A*(2*c+1)) == 0
    expected_next = [[j, j+2 if j <= 130 else j-2] for j in range(108,132)]
    assert primary['triple']['zero_outer_next_walls'] == expected_next
    outer = S.Matrix([1,0,1])
    assert (outer.T*triple*outer)[0] == -8*c
    # Replacing A by0 and B,D by2 maximizes u=-q on the nonnegative box.
    u = -q
    assert S.simplify(S.diff(u,B) - (B+D+2*c*D)/(4*c*(c+1))) == 0
    assert S.simplify(S.diff(u,D) - (B+D+2*c*B)/(4*c*(c+1))) == 0
    assert S.simplify(u.subs({A:0,B:2,D:2})) == 2/c

    support = chain['support']
    G = S.Matrix([[rat(v) for v in row] for row in chain['support_gram']])
    assert G.is_positive_definite is True
    assert all(G[i,j] <= 0 for i in range(18) for j in range(18) if i != j)
    J = G.inv()
    assert J == S.Matrix([[rat(v) for v in row] for row in chain['support_inverse']])
    assert all(v >= 0 for v in J)
    assert min(J[i,i] for i in range(18)) == S.Rational(1,2)
    inside = [support.index(i) for i in range(1,11)]
    outside = sorted(set(range(18))-set(inside))
    assert J.extract(inside,inside) == S.Matrix([
        [S.Rational(min(i,j)*(11-max(i,j)),11) for j in range(1,11)]
        for i in range(1,11)])
    assert all(J[i,j] == 0 for i in inside for j in outside)
    lo, hi = map(rat, chain['parameter']['isolating_interval'])
    assert lo > 0 and lo*lo > 280
    total, small = S.Rational(281,70), S.Rational(141,70)
    assert 4+4/(lo*lo) < total and small == total-2
    assert total-S.Rational(40,11) < S.Rational(1,2)
    assert 4*S.Rational(18,11) > total
    assert S.Rational(24,11) > small
    assert S.Rational(10,11)*S.Rational(5,2) > small
    assert S.Rational(18,11)*2 > small

    saved = cap['strict_escape']
    assert saved['support'] == support
    assert list(map(rat,saved['frame_coefficients'])) == [1,-S.Rational(1,100),0]
    eps = S.Rational(1,100000)
    eb = S.Matrix(list(map(rat,saved['dual_coefficients'])))
    assert list(eb) == [-1 if wall == 5 else -eps for wall in support]
    e0 = S.Matrix([-1 if wall == 5 else 0 for wall in support])
    assert all(v >= 0 for v in e0-eb)
    # Their normal coefficients are also nonnegative: E0-Eeps=Σpositive*n_i.
    normal_difference = J*(e0-eb)
    assert all(v >= 0 for v in normal_difference) and any(v > 0 for v in normal_difference)
    assert len(saved['all134_pair_bounds']) == 134
    assert all(rat(pair[1]) < 0 for pair in saved['all134_pair_bounds'])
    assert rat(saved['strict_center_pair_bounds'][1]) < 0
    K = S.Rational(9999,10000)
    norm = (e0.T*J*e0)[0]-K
    assert norm == S.Rational(190011,110000)
    assert norm == rat(primary['weak_escape']['norm'])
    z = -J*e0
    assert all(z[j] == 0 for j in outside)
    endpoints = []
    for wall in (10,1):
        j = support.index(wall)
        lead = J[j,j]*K-z[j]**2
        lower = 4*J[j,j]-2
        margin = lower*lead-2*z[j]**2
        assert lead > 0 and margin > 0
        row = next(r for r in primary['dotted_A10_endpoint_cases'] if r['wall'] == wall)
        assert rat(row['linear_coefficient']) == lead
        assert rat(row['strict_margin']) == margin
        endpoints.append(dict(wall=wall,linear_coefficient=str(lead),margin=str(margin)))

    expected = {tuple(r['labels']):r for r in primary['A10_profiles']}
    rejected = {tuple(r['labels']):r for r in primary['A10_excluded_profiles']}
    kept, excluded = [], 0
    slots = [10,9,2,1]
    choices = [(0,1,S.sqrt(2)),(0,1),(0,1),(0,1,S.sqrt(2))]
    for values in product(*choices):
        d = S.zeros(18,1)
        for wall,value in zip(slots,values):
            d[support.index(wall)] = value
        cost = S.expand((d.T*J*d)[0])
        pay = S.expand((d.T*z)[0])
        labels = tuple(reversed([0 if v == 0 else 3 if v == 1 else 4 for v in values]))
        if cost == small or positive(cost-small):
            excluded += 1
            assert S.simplify(cost-radical(rejected[labels]['cost'])) == 0
            continue
        assert positive(small-cost)
        margin = S.expand(K*cost-pay*pay)
        assert cost == 0 and pay == 0 or positive(margin)
        row = expected[labels]
        assert cost == rat(row['A10_cost'])
        assert S.simplify(pay-radical(row['escape_A10_pair_magnitude'])) == 0
        assert S.simplify(margin-radical(row['strict_squared_margin'])) == 0
        kept.append(dict(labels=list(labels),cost=str(cost),pairing=str(pay),margin=str(margin)))
    assert len(kept) == 8 and excluded == 28
    assert {tuple(r['labels']) for r in kept} == set(expected)
    paths = [f'scripts/{TAG}.py',PRIMARY,CHAIN,CAP,CAP_REVIEW,ABEL_REVIEW]
    return dict(
        status='PASS_INDEPENDENT_ABELIAN_BASE_DEGREE50_AND_ASYMMETRIC_COMMON_ESCAPE',
        sources={p:sha(p) for p in paths},
        inherited_modular_factorizations=abelrev['modular_factorizations'],
        proper_block_size_and_count=[list(row) for row in proper_blocks],
        powered_Frobenius=dict(fixed_points=fixed,nontrivial_orbit_length=len(orbit)),
        degree_over_K_Qab=50,
        symbolic_triple_norm=str(q),triple_determinant=str(S.factor(triple.det())),
        zero_outer_exception_next_walls=expected_next,
        upper_bound_for_u='2/C',support_budget=str(total),
        weak_escape_norm=str(norm),
        weak_minus_strict_normal_coefficients=list(map(str,normal_difference)),
        continuous_endpoint_checks=endpoints,
        A10_profiles=sorted(kept,key=lambda r:r['labels']),excluded_profiles=excluded,
        scope='Independent symbolic Gram inversion and box bound, permutation/block arithmetic and reviewed normal-subgroup/base-change proof, rational support inversion, weak escape inheritance, complete36-profile census and continuous endpoint margins. Inherits independently checked modular factorization and physical chamber. Applies to any finite coherent collection of caps with consecutive nondotted triples; no new dimension bound or finite-volume construction.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    out = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
