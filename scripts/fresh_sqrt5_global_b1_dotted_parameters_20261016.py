#!/usr/bin/env python3
"""Complete necessary arithmetic distance domains for all minimal B1 dots.

The one-dot Schur bound and the separate two-dot product bound are exact.
These are parameter seeds, not an attachment closure or a dimension theorem.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path

import fresh_sqrt5_dotted_star_closure_20261015 as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_global_b1_dotted_parameters_20261016'
INPUT = 'results/fresh_sqrt5_adjacent_discount_frontier_20261016.json'
D, Q, K = S.D, S.Q, S.K
PHI = Q.PHI


def floor_k(x):
    """Exact floor, using rational bounds and bisection, never a float."""
    lo = (x.a - 3 * abs(x.b)).__floor__() - 1
    hi = (x.a + 3 * abs(x.b)).__floor__() + 2
    assert (x - lo).sign() > 0 and (hi - x).sign() > 0
    while hi - lo > 1:
        mid = (hi + lo) // 2
        if (x - mid).sign() >= 0:
            lo = mid
        else:
            hi = mid
    assert (x - lo).sign() >= 0 and (lo + 1 - x).sign() > 0
    return lo


def ceil_k(x):
    return -floor_k(-x)


def embedding_box(bound):
    """All t in Z[phi] with 4<t<=bound and 0<=t'<=4."""
    if (bound - 4).sign() <= 0:
        return []
    answer = []
    # t-t'=B*sqrt5, so B is positive and B*sqrt5<=bound.
    for b in range(1, floor_k(bound / K(0, 1)) + 1):
        low, high = ceil_k(-b * PHI.conj()), floor_k(4 - b * PHI.conj())
        for a in range(low, high + 1):
            t = K(a) + b * PHI
            if (t - 4).sign() <= 0 or (bound - t).sign() < 0:
                continue
            assert t.conj().sign() >= 0 and (4 - t.conj()).sign() >= 0
            assert D.integral(t)
            answer.append(dict(A=a, B=b, value=t.enc()))
    assert len(answer) == len({(r['A'], r['B']) for r in answer})
    return answer


def matrix_record(M):
    n = len(M)
    p = Q.inertia(M)
    c = Q.inertia([[x.conj() for x in row] for row in M])
    out = dict(physical=p, conjugate=c)
    if p[1] > 1:
        return dict(**out, status='REJECT_TWO_NEGATIVE')
    if c[1]:
        return dict(**out, status='REJECT_CONJUGATE')
    assert p[1] == 1 and p[2] in (0, 1)
    assert c == [n - p[2], 0, p[2]]
    cycles = D.cycles(M)
    bad = next((r for r in cycles if not r['integral']), None)
    if bad:
        return dict(**out, status='REJECT_CYCLIC_INTEGRALITY', witness=bad)
    out.update(gram=D.enc(M), cycles=cycles)
    if p[2]:
        vectors = Q.kernel(M)
        assert len(vectors) == 1
        signs = [x.sign() for x in vectors[0]]
        out.update(kernel=[x.enc() for x in vectors[0]], kernel_signs=signs)
        if all(x >= 0 for x in signs) or all(x <= 0 for x in signs):
            return dict(**out, status='REJECT_NONNEGATIVE_ROOT_DEPENDENCE')
        # A complete quadrilateral star must have this sign pattern.
        # Record rather than silently discard a different pattern.
        out['black_kernel_alternates'] = (signs[0] == signs[2] and
                                          signs[1] == signs[3] and
                                          signs[0] * signs[1] == -1)
    return dict(**out, status='ACCEPT_SINGULAR' if p[2] else 'ACCEPT_NONSINGULAR')


def one_dot(row):
    preliminary = D.evaluate(row)
    if preliminary['status'].startswith('REJECT_'):
        return dict(status=preliminary['status'], preliminary=preliminary,
                    trials=[], accepted=[])
    assert 'schur' in preliminary, preliminary['status']
    schur = preliminary['schur']
    a, b, c = (K(*schur[k]) for k in ('a', 'b', 'c'))
    M = Q.tree_scaled(row['nodes'], {p: 2 if m == 0 else m
                                   for p, m in Q.edges(row).items()})
    bound = 8 * (a * b + c * c) / (M[0][0] * M[2][2])
    # If the two finite components are disconnected, c=0 and x²<=ab.
    if not preliminary['finite_path']:
        assert not c
        bound /= 2
    trials, accepted = [], []
    for r in embedding_box(bound):
        t = K(*r['value'])
        G, error = S.matrix(row, t / 4)
        z = dict(status=error) if error else matrix_record(G)
        trial = dict(**r, **z)
        if trial['status'].startswith('ACCEPT_'):
            accepted.append(dict(trial_index=len(trials), ratio=(t / 4).enc(),
                                 singular=trial['status'] == 'ACCEPT_SINGULAR'))
        trials.append(trial)
    return dict(status='COMPLETE_ONE_DOT_NECESSARY_PARAMETER_DOMAIN',
                preliminary=preliminary, upper_bound=bound.enc(),
                integer_upper_bound=ceil_k(bound), trials=trials, accepted=accepted)


def double_matrix(t, u):
    """Field-valued gauge for the path 0--dot--2--3--3--dot--1."""
    M = [[K() for _ in range(4)] for _ in range(4)]
    M[0][0], M[1][1], M[2][2], M[3][3] = t / 2, u / 2, K(2), K(2)
    M[0][2] = M[2][0] = -t / 2
    M[1][3] = M[3][1] = -u / 2
    M[2][3] = M[3][2] = K(-1)
    return M


def two_dot(row):
    assert row['nodes'] == 4 and row['labels'] == [2, 0, 2, 2, 0, 3]
    minimum = 3 + PHI
    # B>=3 forces t>=3sqrt5>minimum; B=1,2 are checked explicitly.
    assert (3 * K(0, 1) - minimum).sign() > 0
    low_box = embedding_box(minimum)
    assert low_box == [dict(A=3, B=1, value=minimum.enc())]
    bound = 4 + K(4) / (minimum - 4)
    assert bound == 4 + 4 * PHI
    box = embedding_box(bound)
    trials, accepted = [], []
    for tr, ur in product(box, repeat=2):
        t, u = K(*tr['value']), K(*ur['value'])
        physical_product = (t - 4) * (u - 4)
        conjugate_product = (4 - t.conj()) * (4 - u.conj())
        z = matrix_record(double_matrix(t, u))
        # Independent scalar Schur predicates for the explicit four-node gauge.
        assert (z['physical'][1] <= 1) == ((4 - physical_product).sign() >= 0)
        assert (z['conjugate'][1] == 0) == ((conjugate_product - 4).sign() >= 0)
        trial = dict(t=tr, u=ur, physical_product=physical_product.enc(),
                     conjugate_product=conjugate_product.enc(), **z)
        if trial['status'].startswith('ACCEPT_'):
            accepted.append(dict(trial_index=len(trials),
                                 dotted_ratios={'0,2': (t / 4).enc(), '1,3': (u / 4).enc()},
                                 singular=trial['status'] == 'ACCEPT_SINGULAR'))
        trials.append(trial)
    return dict(status='COMPLETE_TWO_DOT_NECESSARY_PARAMETER_DOMAIN',
                minimum_integral_parameter=minimum.enc(),
                lower_bound_finite_check=low_box, upper_bound=bound.enc(),
                parameter_box=box, trials=trials, accepted=accepted)


def run():
    src = json.loads((ROOT / INPUT).read_text())
    rows = [r for r in src['records'] if r['analysis']['baseline'] == 1
            and r['inherited_status'] == 'DOTTED_STRUCTURAL_ONE_WHITE_DOMAIN']
    assert len(rows) == 82
    records, counts, trials = [], Counter(), Counter()
    for i, r in enumerate(rows):
        n_dots = r['row']['labels'].count(0)
        assert n_dots in (1, 2)
        z = one_dot(r['row']) if n_dots == 1 else two_dot(r['row'])
        rec = dict(index=i, class_index=r['class_index'], union_index=r['union_index'],
                   row=r['row'], dotted_pairs=n_dots, result=z)
        records.append(rec)
        counts[z['status']] += 1
        trials.update(t['status'] for t in z['trials'])
        print('B1_DOTTED_DOMAIN', i, r['class_index'], n_dots,
              z['status'], len(z['trials']), len(z['accepted']), flush=True)
    paths = [f'scripts/{TAG}.py', INPUT, f'scripts/{S.TAG}.py',
             f'scripts/{D.TAG}.py', 'scripts/fresh_sqrt5_quad_contacts_20261015_review.py',
             'scripts/fresh_sqrt5_singular_extension_frontier_20261015.py']
    return dict(status='PASS_COMPLETE_MINIMAL_B1_DOTTED_PARAMETER_DOMAINS',
                sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                minimal_classes=len(rows), counts=dict(counts), trial_counts=dict(trials),
                accepted_parameter_seeds=sum(len(r['result']['accepted']) for r in records),
                records=records,
                scope='Complete necessary integral distance domains for all82 baseline1 minimal dotted cores in the inherited complete frontier:81 one-dot and one two-dot. Physical/conjugate signatures and every cyclic product are exact. Principal rejections and singular seed relations are explicit. All distances compatible with these necessary tests are retained; arbitrary white attachments, face payment, and global compact realization remain separate. No new dimension bound.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    data = json.loads(json.dumps(run()))
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], data['accepted_parameter_seeds'], data['trial_counts'], flush=True)
