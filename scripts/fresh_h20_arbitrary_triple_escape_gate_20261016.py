#!/usr/bin/env python3
"""Exact arbitrary-three-chain-contact obstruction and common weak escape.

Standard-library polynomial/RREF/Fraction arithmetic. The abelian-base
degree-50 theorem and physical retained chamber are explicitly inherited.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_arbitrary_triple_escape_gate_20261016'
BASE = 'results/fresh_h20_s61_multicap_20261014_exact.json'
UPGRADE = 'results/fresh_h20_asymmetric_arc_escape_gate_20261016.json'
UPGRADE_SOURCE = 'scripts/fresh_h20_asymmetric_arc_escape_gate_20261016.py'
UPGRADE_REVIEW = 'results/fresh_h20_asymmetric_arc_escape_gate_20261016_review.json'
ORDER = ['q_y', 'f1_squared', 'f2_squared', 'f3_squared', 'f1_f2', 'f1_f3', 'f2_f3']


def digest(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(a, b):
    c = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        c[i] += x
    for i, x in enumerate(b):
        c[i] += x
    return trim(c)


def scale(a, z):
    return trim([z * x for x in a])


def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return trim(c)


def interval_add(a, b):
    return a[0] + b[0], a[1] + b[1]


def interval_mul(a, b):
    vals = [x * y for x in a for y in b]
    return min(vals), max(vals)


def evaluate(p, x):
    out = (F(0), F(0))
    for c in reversed(p):
        out = interval_add(interval_mul(out, x), (F(c), F(c)))
    return out


def encoded_interval(x):
    return [str(v) for v in x]


def rref_kernel(rows):
    a = [[F(v) for v in row] for row in rows]
    pivots = []
    k = 0
    for col in range(7):
        candidates = [i for i in range(k, len(a)) if a[i][col]]
        if not candidates:
            continue
        j = candidates[0]
        a[k], a[j] = a[j], a[k]
        den = a[k][col]
        a[k] = [v / den for v in a[k]]
        for i in range(len(a)):
            if i != k and a[i][col]:
                c = a[i][col]
                a[i] = [u - c * v for u, v in zip(a[i], a[k])]
        pivots.append(col)
        k += 1
    basis = []
    for col in range(7):
        if col in pivots:
            continue
        v = [F(0)] * 7
        v[col] = 1
        for i, p in enumerate(pivots):
            v[p] = -a[i][col]
        assert all(sum(x * y for x, y in zip(row, v)) == 0 for row in rows)
        basis.append(v)
    assert len(basis) + len(pivots) == 7
    return pivots, [[str(x) for x in row] for row in basis]


EXPECTED = {
    'consecutive': [['1/2', '0', '1', '0', '0', '0', '0'],
                    ['0', '-1', '0', '1', '0', '0', '0'],
                    ['0', '0', '0', '0', '1', '0', '0'],
                    ['0', '0', '0', '0', '0', '0', '1']],
    'left_1_right_2': [['-1/2', '-1', '-1', '1', '0', '0', '0'],
                       ['0', '2', '0', '0', '0', '1', '0']],
    'left_2_right_1': [['1/2', '-1', '1', '1', '0', '0', '0'],
                       ['-1', '2', '-2', '0', '0', '1', '0']],
    'left_1': [['-1/2', '-1', '-1', '1', '0', '0', '0']],
    'right_1': [['1/2', '-1', '1', '1', '0', '0', '0']],
    'equal_gaps': [['0', '-1', '0', '1', '0', '0', '0'],
                   ['0', '0', '0', '0', '-1', '0', '1']],
    'full_rank': [],
}


def category(p, q):
    if (p, q) == (1, 1):
        return 'consecutive'
    if (p, q) == (1, 2):
        return 'left_1_right_2'
    if (p, q) == (2, 1):
        return 'left_2_right_1'
    if p == 1:
        return 'left_1'
    if q == 1:
        return 'right_1'
    if p == q:
        return 'equal_gaps'
    return 'full_rank'


def qadd(x, y):
    return x[0] + y[0], x[1] + y[1]


def qmul(x, y):
    return x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def qscale(x, z):
    return x[0] * z, x[1] * z


def qsign(x):
    a, b = x
    if not a:
        return (b > 0) - (b < 0)
    if not b or a * b > 0:
        return (a > 0) - (a < 0)
    v = a * a - 2 * b * b
    assert v
    return ((a > 0) - (a < 0)) * ((v > 0) - (v < 0))


def manifest(paths):
    seen = set()
    manifests = references = 0

    def visit(p):
        nonlocal manifests, references
        if p in seen:
            return
        seen.add(p)
        if not p.endswith('.json'):
            return
        d = json.loads((ROOT / p).read_text())
        if 'sources' in d:
            manifests += 1
        for q, h in d.get('sources', {}).items():
            references += 1
            assert digest(q) == h, (p, q)
            visit(q)

    for p in paths:
        visit(p)
    return dict(manifests=manifests, references=references, files=len(seen))


def budget(base):
    S = base['support']
    J = [[F(x) for x in row] for row in base['support_inverse']]
    assert len(S) == 18 and all(x >= 0 for row in J for x in row)
    assert min(J[i][i] for i in range(18)) == F(1, 2)
    A = [S.index(i) for i in range(1, 11)]
    B = [i for i in range(18) if i not in A]
    assert all(J[i][j] == F(min(i + 1, j + 1) * (10 - max(i, j)), 11)
               for i in range(10) for j in range(10))
    assert A == list(range(10)) and all(J[i][j] == 0 for i in A for j in B)
    total = F(281, 70)
    residual = total - 2
    alo = F(base['parameter']['isolating_interval'][0])
    assert 4 + 4 / (alo * alo) < total
    assert total - F(40, 11) < F(1, 2)
    assert 4 * F(18, 11) > total
    assert F(24, 11) > residual
    assert F(10, 11) * F(5, 2) > residual
    assert 2 * F(18, 11) > residual
    K = F(9999, 10000)
    z = [J[i][S.index(5)] for i in range(18)]
    assert all(z[i] == 0 for i in B)
    endpoints = []
    for wall in (1, 10):
        j = S.index(wall)
        lead = K * J[j][j] - z[j] ** 2
        margin = (4 * J[j][j] - 2) * lead - 2 * z[j] ** 2
        assert lead > 0 and margin > 0
        endpoints.append(dict(wall=wall, lead=str(lead), margin=str(margin)))
    domains = [[0, 3, 4], [0, 3], [0, 3], [0, 3, 4]]
    kept, rejected = [], []
    zero = (F(0), F(0))
    angle = {0: zero, 3: (F(1), F(0)), 4: (F(0), F(1))}
    for labels in product(*domains):
        d = [zero] * 18
        for wall, label in zip((1, 2, 9, 10), labels):
            d[S.index(wall)] = angle[label]
        cost = pay = zero
        for i in A:
            pay = qadd(pay, qscale(d[i], z[i]))
            for j in A:
                cost = qadd(cost, qscale(qmul(d[i], d[j]), J[i][j]))
        if qsign(qadd(cost, (-residual, F(0)))) >= 0:
            rejected.append(dict(labels=list(labels), cost=list(map(str, cost))))
            continue
        assert cost[1] == 0
        margin = qadd(qscale(cost, K), qscale(qmul(pay, pay), -1))
        if cost[0]:
            assert qsign(margin) > 0
        else:
            assert pay == margin == zero
        kept.append(dict(labels=list(labels), cost=str(cost[0]),
                         pay=list(map(str, pay)), margin=list(map(str, margin))))
    assert len(kept) == 8 and len(rejected) == 28
    return dict(total_cost_upper=str(total), A10_cost_upper_if_dotted_outside=str(residual),
                endpoint_cases=endpoints, A10_profiles=kept, excluded_profiles=rejected,
                E0_norm=str(J[S.index(5)][S.index(5)] - K), frame_negative_norm=str(K))


def run():
    paths = [BASE, UPGRADE, UPGRADE_SOURCE, UPGRADE_REVIEW]
    base = json.loads((ROOT / BASE).read_text())
    upgrade = json.loads((ROOT / UPGRADE).read_text())
    review = json.loads((ROOT / UPGRADE_REVIEW).read_text())
    assert upgrade['status'] == 'PASS_DEGREE50_AB_BASE_AND_ALL_CONSECUTIVE_ARCS_COMMON_WEAK_ESCAPE'
    assert upgrade['Galois_upgrade']['no_nonzero_polynomial_over_Qab_of_degree_below'] == 50
    assert upgrade['Galois_upgrade']['degree_over_K_Qab'] == 50
    assert upgrade['weak_escape']['off_support_old134_strict']
    assert upgrade['weak_escape']['strict_center_pair_negative']
    assert upgrade['weak_escape']['norm'] == '190011/110000'
    assert review['status'] == 'PASS_INDEPENDENT_ABELIAN_BASE_DEGREE50_AND_ASYMMETRIC_COMMON_ESCAPE'
    assert review['sources'][UPGRADE] == digest(UPGRADE)
    assert review['degree_over_K_Qab'] == 50
    audit = manifest(paths)
    alo, ahi = map(F, base['parameter']['isolating_interval'])
    clo, chi = alo * alo / (alo * alo + 2), ahi * ahi / (ahi * ahi + 2)
    denom = 10 ** 30
    CI = (F(clo * denom // 1, denom), F(chi * denom // 1 + 1, denom))
    assert 0 < CI[0] < clo < chi < CI[1] < 1
    T = [[1], [0, 1]]
    for k in range(2, 26):
        T.append(add(mul([0, 2], T[-1]), scale(T[-2], -1)))
    g = {}
    for k in range(1, 26):
        numerator = scale(add(T[k], [0, -1]), 2)
        # Divide by (1-X): the quotient coefficients are cumulative sums.
        quotient, running = [], 0
        for c in numerator:
            running += c
            quotient.append(running)
        assert quotient[-1] == 0
        g[k] = trim(quotient[:-1] or [0])
        assert mul(g[k], [1, -1]) == numerator
    assert g[1] == [0] and g[2] == [-2, -4]
    gap_bounds = []
    for k in range(1, 26):
        bounds = evaluate(g[k], CI)
        if k >= 2:
            assert bounds[1] < -2
        diff = evaluate(add(g[k], scale(g[k - 1], -1)), CI) if k > 1 else None
        if diff:
            assert diff[1] < 0
        gap_bounds.append(dict(gap=k, polynomial=g[k], bounds=encoded_interval(bounds),
                               difference_bounds=encoded_interval(diff) if diff else None))
    cases = []
    counts = Counter()
    for p in range(1, 25):
        for q in range(1, 26 - p):
            a, b, c = g[p], g[q], g[p + q]
            det = add(add([8], scale(mul(mul(a, b), c), 2)),
                      scale(add(add(mul(a, a), mul(b, b)), mul(c, c)), -2))
            columns = [det, add(mul(b, b), [-4]), add(mul(c, c), [-4]),
                       add(mul(a, a), [-4]),
                       add(scale(mul(b, c), -2), scale(a, 4)),
                       add(scale(mul(a, b), -2), scale(c, 4)),
                       add(scale(mul(a, c), -2), scale(b, 4))]
            degree = max(len(z) - 1 for z in columns)
            assert degree <= 2 * (p + q) - 2 <= 48
            rows = [[col[i] if i < len(col) else 0 for col in columns]
                    for i in range(degree + 1)]
            pivots, kernel = rref_kernel(rows)
            kind = category(p, q)
            assert kernel == EXPECTED[kind], (p, q, kernel)
            bounds = evaluate(det, CI)
            assert bounds[1] < 0
            counts[kind] += 1
            cases.append(dict(gaps=[p, q], translations=26 - p - q, degree=degree,
                              coefficient_columns=columns, pivot_columns=pivots,
                              kernel_basis=kernel, classification=kind,
                              determinant_bounds=encoded_interval(bounds)))
    assert len(cases) == 300 and sum(row['translations'] for row in cases) == 2600
    assert counts == dict(consecutive=1, left_1_right_2=1, left_2_right_1=1,
                          left_1=22, right_1=22, equal_gaps=11, full_rank=242)
    exceptions = []
    for mid in range(1, 25):
        witness = mid + 2 if mid + 2 < 26 else mid - 2
        assert 0 <= witness < 26 and abs(witness - mid) == 2
        exceptions.append(dict(middle_chain_index=mid, positive_pair_chain_index=witness,
                               positive_pair_factor=[1, 2]))
    exact_budget = budget(base)
    assert exact_budget['E0_norm'] == upgrade['weak_escape']['norm']
    paths.append(f'scripts/{TAG}.py')
    return dict(status='PASS_ARBITRARY_THREE_CHAIN_CONTACTS_COMMON_WEAK_ESCAPE',
                sources={p: digest(p) for p in paths}, recursive_manifest_audit=audit,
                inherited_degree_over_K_Qab=50,
                coefficient_variable_order=ORDER, polynomial_coefficient_order='ascending powers of C',
                C_interval=encoded_interval(CI), chain_pair_polynomials=gap_bounds,
                gap_case_count=300, actual_triple_count=2600, maximum_degree=48,
                classification_counts=dict(sorted(counts.items())), cases=cases,
                exceptional_consecutive_acuteness_witnesses=exceptions,
                kernel_rank_one_reductions={
                    'full_rank': 'The coefficient vector is zero.',
                    'left_1': 'f2^2+f3^2=0 forces both zero, then all coordinates zero.',
                    'right_1': 'f1^2+f2^2=0 forces both zero, then all coordinates zero.',
                    'equal_gaps': 'f1^2+f3^2=0 and f2^2=0; actual cross products then vanish.',
                    'left_1_right_2': 'f2^2+f3^2=0; then f1*f3=0 kills the second basis coefficient.',
                    'left_2_right_1': 'f1^2+f2^2=0; then f1*f3=0 kills the second basis coefficient.',
                    'consecutive': 'f1^2+f3^2=0; actual products vanish, leaving q_y=f2^2/2.'},
                coherent_zero_exclusion='Three independent chain normals imply y=0. Then r=-sum((Jd)_i*n_i), with Jd>=0 and r nonzero, pairs positively with every common strict interior point.',
                consecutive_nonzero_exclusion='y=-f2*n_middle/2; the saved separation2 wall has positive pairing f2*(1+2C).',
                metric_bound=dict(assumptions='Any two nondotted chain contacts with index separation n>=2, through25.',
                    projection='u=-y^2 <= (2*f1^2+2*f2^2-2*g_n*f1*f2)/(g_n^2-4)',
                    box_bound='For 0<=f1,f2<=2, u<=8/(-g_n-2)<=2/C.',
                    mechanism='The pair spans signature(1,1); its orthogonal complement in Lorentz3 is positive definite.',
                    future_sheet='Once a dotted support contact forces u>0, y pairs nonpositively with the future timelike outer-root sum, and hence is future timelike.'),
                weak_escape=upgrade['weak_escape'], independently_recomputed_budget=exact_budget,
                conclusion='Every coherent Coxeter cap retaining the old chamber and having any three nondotted contacts among its26 chain walls pairs strictly negatively with E0. Thus any cap cutting E0 has at most2 nondotted chain contacts and at least24 dotted chain contacts.',
                collection_conclusion='Every finite mutually coherent collection of caps each having at least3 nondotted chain contacts preserves a common strict spacelike escape after perturbation toward a shared interior point; the intersection remains infinite-volume.',
                scope='All300 positive gap pairs and all2600 chain triples; all integer finite labels and parallel contacts; arbitrary18 support Coxeter contacts and arbitrary real cap coefficients. This is a necessary condition for cutting this particular E0 in the fixed retained chamber. It does not imply every cap has24 dotted contacts, does not cover a cutter with at most2 nondotted chain contacts, and gives no dimension bound or new finite-volume polytope. Degree50 over K Qab and physical old134 geometry are inherited dependencies.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['gap_case_count'], result['actual_triple_count'], flush=True)
