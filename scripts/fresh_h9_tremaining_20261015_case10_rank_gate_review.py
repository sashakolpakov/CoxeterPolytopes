#!/usr/bin/env python3
"""Independent literal-incidence and exact-field review of the T10 rank gate.

No project arithmetic imports. SymPy's primitive-element number field differs
from the producer's four-coordinate Fraction implementation. Determinants use
elimination/direct cofactor formulas, and Chebyshev values use a closed sum.
"""
import argparse
import json
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_tremaining_20261015_case10_rank_gate_review'
PRIMARY = 'fresh_h9_tremaining_20261015_case10_rank_gate'
HIGH = 'fresh_h9_tremaining_20261015_case10_high_bound'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'
STATE = 'results/fresh_h9_tremaining_20261015_search_case_10_batch_001.json'
DOM = 'results/fresh_h9_tremaining_20261015_case10_domains.json'
K = sp.QQ.algebraic_field(sp.sqrt(2), sp.sqrt(3))
zero, one = K.zero, K.one
half = K.convert(sp.Rational(1, 2))
rt2, rt3 = map(K.from_sympy, (sp.sqrt(2), sp.sqrt(3)))
COS = {2: zero, 3: half, 4: rt2 * half, 6: rt3 * half}


def decode(xs):
    a, b, c, d = [K.convert(sp.Rational(x)) for x in xs]
    return a + b * rt2 + c * rt3 + d * rt2 * rt3


def literal(x):
    return str(K.to_sympy(x))


def sign(x):
    if x == zero:
        return 0
    y = K.to_sympy(x)
    if y.is_positive is True:
        return 1
    assert y.is_negative is True, y
    return -1


def det(A):
    A = [row[:] for row in A]
    ans = one
    for j in range(len(A)):
        pivot = next((i for i in range(j, len(A)) if A[i][j] != zero), None)
        if pivot is None:
            return zero
        if pivot != j:
            A[j], A[pivot] = A[pivot], A[j]
            ans = -ans
        p = A[j][j]
        ans *= p
        for i in range(j + 1, len(A)):
            f = A[i][j] / p
            for k in range(j + 1, len(A)):
                A[i][k] -= f * A[j][k]
            A[i][j] = zero
    return ans


def chebyshev_closed(n, x):
    # T_n(x) = n/2 sum_k (-1)^k (n-k-1)!/[k!(n-2k)!] (2x)^(n-2k).
    return sum((K.convert(sp.Rational(n * (-1)**k * comb(n-k, k), 2*(n-k)))
                * (2*x)**(n-2*k) for k in range(n//2+1)), zero)


def run():
    checked = set()
    def sha(p):
        return sha256((ROOT / p).read_bytes()).hexdigest()
    def read(p):
        data = json.loads((ROOT / p).read_text())
        if p not in checked:
            checked.add(p)
            for source, digest in data.get('sources', {}).items():
                assert sha(source) == digest, source
                if source.endswith('.json'):
                    read(source)
        return data
    primary = read(f'results/{PRIMARY}.json')
    high = read(f'results/{HIGH}.json')
    high_review = read(f'results/{HIGH}_review.json')
    geo, state, domains = map(read, (GEOM, STATE, DOM))
    assert primary['status'] == 'PASS_T10_FIRST_MANDATORY_FAMILY_ALL_HIGH_LABELS_EXCLUDED'
    assert high['status'] == 'PASS_T10_SAVED_FAMILY_HIGH_LABEL_AT_MOST_ELEVEN'
    assert high_review['status'] == 'PASS_INDEPENDENT_T10_HIGH_LABEL_DETERMINANT_REVIEW'
    assert high_review['remaining_high_labels'] == [7, 8, 9, 10, 11]
    assert state['status'] == 'FULL_MANDATORY_BRIDGE_SURVIVOR'
    assert state['case'] == 10 and len(state['pending_frames']) == 10
    fixed = {(a, b): m for a, b, m in state['solution']['pair_labels']}
    assert state['solution']['pair_labels'] == primary['fixed_mandatory_family']
    assert state['solution']['pair_labels'] == high['fixed_mandatory_family']
    geometry = geo['representatives'][10]['geometry']
    assert geometry['dimension'] == 9
    vertices = [set(v['incident']) for v in geometry['vertices']]
    blocks = [(3, 4), (6, 7), (9, 10), (12, 13)]
    S = sum((list(b) for b in blocks), [])
    assert [fixed[b] for b in blocks] == [4, 4, 7, 4]
    assert primary['positive_support'] == S
    required_right = set(combinations(S, 2)) - set(blocks)
    required_right.update([(9, 11), (10, 11), (3, 16), (4, 16),
                           (6, 16), (7, 16), (11, 16)])
    assert len(required_right) == 31
    assert {tuple(row['pair']) for row in high_review['forced_right_guards']} == required_right
    for a, b in [(2, 3), (2, 4), (11, 12), (11, 13), (12, 16), (13, 16), (10, 16)]:
        assert fixed[a, b] == 3
    assert fixed[9, 16] == 2
    core = S + [11, 16, 2]
    assert len(set(core)) == 11
    assert all(any(set(pair) <= vertex for vertex in vertices) for pair in combinations(core, 2))

    # Independent leading principal minors for the two new upper-contact guards.
    domain_rows = {tuple(row['pair']): row for row in domains['all_domains']}
    upper_guards = []
    for pair in [(2, 11), (2, 16)]:
        row = domain_rows[pair]
        assert row['domain'] == [2, 3]
        witness = next(w for w in row['failures'] if w['label'] == 4)
        support = witness['support']
        assert len(support) == 9 and set(support) == vertices[witness['actual_vertex']]
        A = []
        for a in support:
            line = []
            for b in support:
                p = tuple(sorted((a, b)))
                m = 4 if p == pair else fixed.get(p, 2)
                if m == 7:
                    m = 6
                line.append(one if a == b else -COS[m])
            A.append(line)
        minors = [det([row[:j] for row in A[:j]]) for j in range(1, 10)]
        assert all(sign(x) > 0 for x in minors[:-1]) and sign(minors[-1]) <= 0
        stored = next(w for w in primary['new_upper_guards'] if w['pair'] == list(pair))
        assert stored['support'] == support
        running = one
        for x, pivot in zip(minors, stored['pivots']):
            running *= decode(pivot)
            assert x == running
        upper_guards.append(dict(pair=list(pair), actual_vertex=witness['actual_vertex'],
                                 support=support, leading_minors=list(map(literal, minors))))

    # The retained seven-node positive blocks in each budget have nonnegative inverses.
    # Reconstruct their low-B2 inverse from a literal matrix instead of importing J.
    B2 = sp.Matrix([[1, -sp.sqrt(2)/2], [-sp.sqrt(2)/2, 1]])
    B2inv = B2.inv()
    assert B2inv == sp.Matrix([[2, sp.sqrt(2)], [sp.sqrt(2), 2]])
    assert all(x >= 0 for x in B2inv)
    assert sp.simplify(1/(1-sp.cos(sp.pi/6)**2)) == 4
    budgets = []
    for leaf, own, foreign in [(11, (12, 13), (3, 4, 6, 7)),
                                (2, (3, 4), (6, 7, 12, 13))]:
        for removed in own:
            support = sorted((set(S) - {removed}) | {leaf})
            stored = next(w for w in primary['budget_supports']
                          if w['leaf'] == leaf and w['removed_own_root'] == removed)
            assert stored['support'] == support
            assert set(support) <= vertices[stored['actual_vertex']]
            retained = next(x for x in own if x != removed)
            assert fixed[tuple(sorted((leaf, retained)))] == 3
            assert all(tuple(sorted((a, b))) in required_right
                       for a, b in combinations(set(S)-{removed}, 2)
                       if tuple(sorted((a, b))) not in blocks)
            budgets.append(dict(leaf=leaf, removed=removed, support=support,
                                actual_vertex=stored['actual_vertex']))
    cost_bounds = dict(retained_own=sp.Rational(1,4), low_single3=sp.Rational(1,2),
                       two_low3=1, low_single4=1, high_single3=1)
    assert all(cost_bounds['retained_own'] + cost_bounds[k] > 1
               for k in ['two_low3', 'low_single4', 'high_single3'])

    # Build low-block projected scalar products directly, with the high contribution
    # separated as d. Vector slots use the reversed support order, independently of
    # the producer's eight-coordinate J storage.
    order = list(reversed(S))
    low_blocks = [(3, 4), (6, 7), (12, 13)]
    def profile(own, hit):
        return {i: half if i in own or i == hit else zero for i in order}
    def cost(p, q):
        return sum((2*(p[a]*q[a]+p[b]*q[b]) + rt2*(p[a]*q[b]+p[b]*q[a])
                    for a, b in low_blocks), zero)
    own16 = profile((12, 13), None)
    u = one + rt2*half
    assert cost(own16, own16) == u
    saved_cases = {(r['anchor11_foreign'], r['third2_foreign'], r['pair_11_2'], r['pair_16_2']): r
                   for r in primary['cases']}
    assert len(saved_cases) == 100
    outcomes = []
    survivors = []
    for hit2, hit11, label2, label1 in product([None, 6, 7, 12, 13],
                                              [None, 3, 4, 6, 7], [2, 3], [2, 3]):
        p, t = profile((12, 13), hit11), profile((3, 4), hit2)
        aa, tt = one-cost(p,p), one-cost(t,t)
        av, at = -cost(p,own16), -COS[label1]-cost(p,t)
        vt = -COS[label2]-cost(own16,t)
        slope = at*at-aa*tt
        intercept = aa*((one-u)*tt-vt*vt)-av*(av*tt-at*vt)+at*(av*vt-at*(one-u))
        assert slope != zero
        d = -intercept/slope
        assert d != zero
        z = one-one/(2*d)
        key = (hit11, hit2, label1, label2)
        saved = saved_cases[key]
        for name, value in [('coefficient_A', slope), ('coefficient_B', intercept),
                            ('unique_d_root', d), ('cosine2theta', z)]:
            assert decode(saved[name]) == value, (key, name)
        accepted = []
        for m in range(7, 12):
            value = chebyshev_closed(m, z)
            assert value == decode(saved['chebyshev_values'][str(m)]), (key, m)
            if value == one:
                accepted.append(m)
        assert accepted == saved['necessary_integer_labels']
        outcome = dict(profile=list(key), necessary_labels=accepted)
        if accepted:
            assert accepted == [8] and d == u
            minor = aa*tt-at*at
            assert sign(aa) < 0 and sign(tt) < 0 and sign(minor) > 0
            assert det([[aa,av,at],[av,one-u-d,vt],[at,vt,tt]]) == zero
            rec = saved['negative_plane_11_2']
            assert [decode(rec[k]) for k in ['norm11','norm2','inner','determinant']] == [aa,tt,at,minor]
            outcome['negative_plane'] = list(map(literal, (aa,tt,at,minor)))
            survivors.append(list(key))
        outcomes.append(outcome)
    assert sorted(survivors, key=str) == sorted([[6,None,2,2],[6,7,2,3],[7,None,2,2],[7,6,2,3]], key=str)
    # Supplementary source-binder audit: literal twin permutations, without
    # importing or running the producer binder and without depending on its JSON.
    vertex_sets = {frozenset(v) for v in vertices}
    missing = {frozenset(m) for m in geometry['minimal_nonfaces']}
    transported = set()
    for mask in range(16):
        perm = list(range(geometry['facet_count']))
        for j, (a,b) in enumerate(blocks):
            if mask & (1 << j):
                perm[a], perm[b] = b,a
        assert {frozenset(perm[i] for i in v) for v in vertex_sets} == vertex_sets
        assert {frozenset(perm[i] for i in m) for m in missing} == missing
        family = tuple(sorted((*sorted((perm[a],perm[b])),m) for (a,b),m in fixed.items()))
        assert sum(m == 7 for _,_,m in family) == 1
        transported.add(family)
    assert len(vertex_sets) == 444 and len(transported) == 2
    paths = [f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json',GEOM,STATE,DOM,
             f'scripts/{HIGH}.py',f'results/{HIGH}.json',f'scripts/{HIGH}_review.py',
             f'results/{HIGH}_review.json',
             'scripts/fresh_h9_tremaining_20261015_case10_metric_gates.py',f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_T10_RANK_GATE_REVIEW',
                sources={p:sha(p) for p in paths},
                arithmetic='SymPy primitive-element QQ(sqrt2,sqrt3); no project arithmetic imports; '
                           'independent elimination and closed-sum Chebyshev evaluation',
                literal_case=10,core_meeting_pairs=55,inherited_reviewed_right_guards=31,
                new_upper_guards=upper_guards,budget_faces=budgets,
                cost_lower_bounds={k:str(v) for k,v in cost_bounds.items()},
                residual_cases=outcomes,
                counts=dict(residual_cases=100,exact_chebyshev_tests=500,angle_survivors=4,
                            negative_planes=4,new_upper_guards=2,budget_faces=4),
                supplementary_binder_check=dict(literal_twin_maps=16,complete_vertices=444,
                                                distinct_full_family_images=2,
                                                symbolic_category7='all finite integers >=7',
                                                producer_binder_json_used=False),
                source_review=dict(all_labels='m>=12 inherited from frozen independent high-label review; '
                                              'm7..11 exhaust the remainder',
                                   rank='PD support rank8 in Lorentz space dimension10 leaves residual rank at most2',
                                   monotonicity='Actual nonpositive off-diagonal magnitudes dominate each relaxed '
                                                'guard, so Perron monotonicity preserves non-PD',
                                   profiles='Actual face budgets forbid high contacts, two low contacts, '
                                            'or one low label>=4; untested feasibility only enlarges the case set',
                                   incidence_scope='Only the first saved T10 mandatory family is excluded; '
                                                   'ten pending search frames remain; no traversal rerun'),
                recursive_manifests_checked=len(checked))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    path=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
