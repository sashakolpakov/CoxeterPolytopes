#!/usr/bin/env python3
"""Separate tuple-Fraction review of all82 new minimal dotted domains.

No primary parameter producer or closure imports. Reuses frozen independent
minimal-gauge arithmetic and elementary cycle enumeration.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import fresh_sqrt5_dotted_parameters_20261015_review as P
import fresh_sqrt5_nonsingular_rank29_exhaustion_20261015_review as I

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_sqrt5_global_b1_dotted_parameters_20261016'
TAG = BASE + '_review'
R = P.R
Z, q, dec = R.Z, P.q, P.dec


def integer_sign(a, b):
    """Sign of a+b*sqrt5 with integer coefficients."""
    sign = lambda x: (x > 0) - (x < 0)
    if not a:
        return sign(b)
    if not b or sign(a) == sign(b):
        return sign(a)
    return sign(a) * sign(a * a - 5 * b * b)


def box(bound, ceiling):
    assert R.sign(R.sub(q(ceiling), bound)) >= 0
    assert R.sign(R.sub(q(ceiling - 1), bound)) < 0
    out = {}
    # A deliberately larger symmetric rectangle, unlike the producer's
    # exact-floor strip. Integer signs make its3million tests inexpensive.
    for a in range(-ceiling, ceiling + 1):
        for b in range(-ceiling, ceiling + 1):
            p = 2 * a + b
            if integer_sign(p - 8, b) <= 0:
                continue
            if integer_sign(p, -b) < 0 or integer_sign(8 - p, b) < 0:
                continue
            t = (F(p, 2), F(b, 2))
            if R.sign(R.sub(bound, t)) >= 0:
                out[a, b] = t
    return out


def review_matrix(G, saved):
    p = R.inertia(G)
    c = R.inertia([[R.conj(x) for x in row] for row in G])
    assert p == saved['physical'] and c == saved['conjugate']
    if p[1] > 1:
        return 'REJECT_TWO_NEGATIVE', 0
    if c[1]:
        return 'REJECT_CONJUGATE', 0
    cycles = P.cycles(G)
    if any(not ok for value, ok in cycles.values()):
        return 'REJECT_CYCLIC_INTEGRALITY', len(cycles)
    archived = {P.canonical_cycle(r['cycle']): (dec(r['value']), r['integral'])
                for r in saved['cycles']}
    assert cycles == archived
    if p[2]:
        kernel = I.kernel(G)
        assert len(kernel) == 1
        signs = [R.sign(x) for x in kernel[0]]
        if all(x >= 0 for x in signs) or all(x <= 0 for x in signs):
            return 'REJECT_NONNEGATIVE_ROOT_DEPENDENCE', len(cycles)
        alternating = signs[0] == signs[2] and signs[1] == signs[3] and signs[0] * signs[1] == -1
        assert alternating == saved['black_kernel_alternates']
    return ('ACCEPT_SINGULAR' if p[2] else 'ACCEPT_NONSINGULAR'), len(cycles)


def one_dot(row, saved):
    finite = P.finite_row(row)
    colors = R.coloring(finite)
    assert colors is not None
    M = R.matrix(finite, colors)
    n = row['nodes']
    principal = []
    for omitted in (0, 2):
        ids = [i for i in range(n) if i != omitted]
        B = [[M[i][j] for j in ids] for i in ids]
        p, c = R.inertia(B), R.inertia([[R.conj(x) for x in r] for r in B])
        principal.append(dict(omitted=omitted, physical=p, conjugate=c))
        status = ('REJECT_PHYSICAL_PRINCIPAL' if p[1] > 1 else
                  'REJECT_CONJUGATE_PRINCIPAL' if c[1] else None)
        if status:
            assert saved['status'] == status and saved['trials'] == saved['accepted'] == []
            assert principal == saved['preliminary']['cliques']
            return dict(status=status, principal=principal, trials=0, accepted=0), Counter(), 0
    ids = [i for i in range(n) if i not in (0, 2)]
    B = [[M[i][j] for j in ids] for i in ids]
    assert R.inertia(B) == [n - 3, 1, 0]
    inv = R.inverse(B)
    u, v = [M[0][i] for i in ids], [M[2][i] for i in ids]
    a = R.sub(M[0][0], P.dot(u, [P.dot(r, u) for r in inv]))
    b = R.sub(M[2][2], P.dot(v, [P.dot(r, v) for r in inv]))
    c = P.dot(u, [P.dot(r, v) for r in inv])
    disconnected = 2 not in P.component(row, 0)
    if disconnected:
        assert c == Z
    bound = R.div(R.mul(q(4 if disconnected else 8), R.add(R.mul(a, b), R.mul(c, c))),
                  R.mul(M[0][0], M[2][2]))
    assert bound == dec(saved['upper_bound'])
    candidates = box(bound, saved['integer_upper_bound'])
    table = {(r['A'], r['B']): r for r in saved['trials']}
    assert len(table) == len(saved['trials']) and set(table) == set(candidates)
    counts, accepted, cycle_count = Counter(), set(), 0
    for key, t in candidates.items():
        s = table[key]
        assert t == dec(s['value'])
        G, error = P.matrix(row, R.div(t, q(4)))
        if error:
            status = error
        else:
            status, cy = review_matrix(G, s)
            cycle_count += cy
        assert status == s['status'], (key, status, s['status'])
        counts[status] += 1
        if status.startswith('ACCEPT_'):
            accepted.add(key)
    archived = { (saved['trials'][r['trial_index']]['A'], saved['trials'][r['trial_index']]['B'])
                 for r in saved['accepted'] }
    assert accepted == archived and len(archived) == len(saved['accepted'])
    return dict(status='PASS_ONE_DOT_DOMAIN', upper_bound=R.enc(bound),
                trials=len(table), accepted=len(accepted)), counts, cycle_count


def two_dot(row, saved):
    assert row['nodes'] == 4 and row['labels'] == [2, 0, 2, 2, 0, 3]
    minimum = R.add(q(3), P.PHI)
    assert R.sign(R.sub((F(0), F(3)), minimum)) > 0
    assert box(minimum, 5) == {(3, 1): minimum}
    bound = R.add(q(4), R.div(q(4), R.sub(minimum, q(4))))
    assert bound == dec(saved['upper_bound'])
    candidates = box(bound, 11)
    assert set(candidates) == {(r['A'], r['B']) for r in saved['parameter_box']}
    assert len(candidates) == 12 and len(saved['trials']) == 144
    table = {(r['t']['A'], r['t']['B'], r['u']['A'], r['u']['B']): r for r in saved['trials']}
    counts, accepted, cycle_count = Counter(), set(), 0
    for kt, t in candidates.items():
        for ku, u in candidates.items():
            key = kt + ku
            s = table[key]
            # Different diagonal gauge: root0/root1 are twice as long.
            G = [[Z for _ in range(4)] for _ in range(4)]
            G[0][0], G[1][1], G[2][2], G[3][3] = R.mul(q(2), t), R.mul(q(2), u), q(2), q(2)
            G[0][2] = G[2][0] = R.neg(t)
            G[1][3] = G[3][1] = R.neg(u)
            G[2][3] = G[3][2] = q(-1)
            status, cy = review_matrix(G, s)
            cycle_count += cy
            assert status == s['status']
            counts[status] += 1
            if status.startswith('ACCEPT_'):
                accepted.add(key)
    archived = set()
    for r in saved['accepted']:
        t = saved['trials'][r['trial_index']]
        archived.add((t['t']['A'], t['t']['B'], t['u']['A'], t['u']['B']))
    assert accepted == archived and len(accepted) == 17
    return dict(status='PASS_TWO_DOT_DOMAIN', trials=144, accepted=17), counts, cycle_count


def run():
    primary = P.read(f'results/{BASE}.json')
    origin = P.read('results/fresh_sqrt5_adjacent_discount_frontier_20261016.json')
    source = [r for r in origin['records'] if r['analysis']['baseline'] == 1 and
              r['inherited_status'] == 'DOTTED_STRUCTURAL_ONE_WHITE_DOMAIN']
    assert len(source) == len(primary['records']) == 82
    records, counts, cycles, accepted = [], Counter(), 0, 0
    for origin, rec in zip(source, primary['records']):
        assert origin['row'] == rec['row'] and origin['class_index'] == rec['class_index']
        result, local, ncy = (one_dot(rec['row'], rec['result']) if rec['dotted_pairs'] == 1
                              else two_dot(rec['row'], rec['result']))
        records.append(dict(class_index=rec['class_index'], **result))
        counts.update(local)
        cycles += ncy
        accepted += result['accepted']
        print('INDEPENDENT_B1_DOTTED', rec['index'], rec['class_index'], result['status'], flush=True)
    assert dict(counts) == primary['trial_counts'] and accepted == 377
    paths = [f'scripts/{TAG}.py', f'results/{BASE}.json', f'scripts/{BASE}.py',
             'results/fresh_sqrt5_adjacent_discount_frontier_20261016.json',
             f'scripts/{P.TAG}.py', f'scripts/{I.TAG}.py',
             'scripts/fresh_sqrt5_quad_filters_20261015_review.py']
    return dict(status='PASS_INDEPENDENT_COMPLETE_MINIMAL_B1_DOTTED_PARAMETER_DOMAINS',
                sources=P.manifest(paths), classes=82, accepted_parameter_seeds=accepted,
                counts=dict(counts), checked_integral_cycle_products=cycles, records=records,
                scope='Separate tuple-Fraction field, minimal radical gauges, full principal inverses and inertias, larger symmetric integer boxes and elementary cycle enumeration verify all82 minimal domains and377seeds. No primary parameter/closure imports. The new two-dot bound is separately reconstructed. No arbitrary attachment, payment, compact realization or dimension claim.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    data = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], data['accepted_parameter_seeds'], data['checked_integral_cycle_products'])
