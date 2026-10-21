#!/usr/bin/env python3
"""Exact depth-six / distance-eight corridor diagnostic.

This certifies the retained-profile and adaptive Schur calculations only.
Current short-seven/cycle payment composition remains a separate obligation.
"""
import argparse
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import sympy as S
import fresh_depth8_finite7_profile_transfer_20261018 as P
import fresh_depth6_finite7_adaptive_transfer_20261018 as O

ROOT = P.ROOT
TAG = 'fresh_depth6_finite8_adaptive_pilot_20261018'
CAT = O.CAT
CURRENT = 'fresh_a15_d15_dimension58_20261018'
N, A, B = 21, 6, 14
SQ = S.sqrt(2)


def enc(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {str(k): enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [enc(v) for v in x]
    return x


def profiles():
    white = [i for i in range(N) if i not in (A, B)]
    allowed = [i for i in white if not
               (i > B and i - B <= 6 or A < i < B and B - i <= 6)]
    assert allowed == [0, 1, 2, 3, 4, 5, 7]
    rows, kept = [], []
    for digits in product((0, 1, 2), repeat=len(allowed)):
        if sum(bool(v) for v in digits) > 2:
            continue
        x = {i: S.sqrt(v) for i, v in zip(allowed, digits) if v}
        edge = P.projection(x, white)
        corner = P.projection(x, [i for i in range(N) if i != A])
        accepted = P.field(2-edge)['sign'] > 0 and P.field(2-corner)['sign'] >= 0
        q = P.projection(x, list(range(N)))
        rows.append(dict(contacts={i: str(v) for i, v in x.items()},
                         edge=P.field(edge), adjacent=P.field(corner),
                         selected_projection=P.field(q), accepted=accepted))
        if accepted:
            kept.append(x)
    assert len(rows) == 99
    assert kept == [{}, {7: 1}, {7: SQ}, {5: 1}, {5: 1, 7: 1},
                    {5: SQ}, {4: 1}, {3: 1}, {2: 1}, {1: 1},
                    {0: 1}, {0: 1, 7: 1}, {0: SQ}]
    low = [i for i, x in enumerate(kept)
           if P.field(P.projection(x, list(range(N)))-2)['sign'] <= 0]
    assert low == [0, 9, 10, 12]
    assert min(P.inv(k, i, i) for k in (6, 7) for i in range(k)) == S.Rational(6, 7)
    assert P.inv(N, A, A) == S.Rational(105, 22) > 2
    return kept, dict(candidates=rows, candidate_count=len(rows), accepted_count=len(kept),
                      profiles=[{i: str(v) for i, v in x.items()} for x in kept],
                      allowed=allowed, low_projection_indices=low,
                      three_contacts_projection_lower='18/7',
                      label_at_least5_projection_strict_lower='15/7')


def small_cores(cat):
    rows, survivors = [], []
    for record in cat['diagrams']:
        n = record['order']
        graph = {i: {} for i in range(n)}
        for i, j, label in record['edges']:
            graph[i][j] = graph[j][i] = label
        assert set(record['finite_deletions']) == {
            i for i, deleted in enumerate(record['deletions'])
            if all(t['kind'] == 'finite' for t in deleted)}
        for a, c in product(record['finite_deletions'], repeat=2):
            if a == c:
                continue
            distances, todo = {a: 0}, deque([a])
            while todo:
                i = todo.popleft()
                for j in graph[i]:
                    if j != c and j not in distances:
                        distances[j] = distances[i]+1
                        todo.append(j)
            contacts = [i for i in distances if c in graph[i]]
            assert contacts
            distance = min(distances[i]+1 for i in contacts)
            features = [i for i in distances if
                        sum(j != c for j in graph[i]) >= 3 or
                        any(label != 3 and j != c for j, label in graph[i].items())]
            depth = min((distances[i] for i in features), default=99)
            row = dict(diagram=record['index'], order=n, a=a, c=c,
                       component_distances=distances, component_contact_nodes=contacts,
                       distance=distance, branch_or_double_depth=depth)
            rows.append(row)
            if depth >= 6 and distance >= 6:
                survivors.append(row)
    assert len(cat['diagrams']) == 26 and len(rows) == 768
    assert [(r['diagram'], r['a'], r['c']) for r in survivors] == [(25, 8, 1), (25, 8, 2)]
    core = cat['diagrams'][25]
    assert core['order'] == 10 and core['index'] == 25
    gram = S.diag(*([2]*10))
    for i, j, label in core['edges']:
        gram[i, j] = gram[j, i] = -1 if label == 3 else -SQ
    patches = []
    for kind, cnode, mapping in (
        ('extend', 1, {0: 1, 1: 21, 2: 0, 3: 23, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7}),
        ('isolated', 2, {0: 1, 1: 0, 2: 21, 3: 23, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7})):
        g = S.diag(P.path(N), S.diag(2, 2, 2))
        g[1, 21] = g[21, 1] = -1
        j = 0 if kind == 'extend' else 21
        g[j, 23] = g[23, j] = -1
        ids = [mapping[i] for i in range(10)]
        assert g.extract(ids, ids) == gram
        patches.append(dict(kind=kind, c_node=cnode, map_to_retained_path=mapping,
                            extra_white=23, extra_white_contact=j,
                            contact_magnitude='1', old_path_contact_of_c=1))
    assert P.inv(6, 1, 1) == S.Rational(10, 7)
    assert 1/(2-P.inv(6, 5, 5)) == S.Rational(7, 8)
    assert S.Rational(10, 7)+S.Rational(7, 8) == S.Rational(129, 56) > 2
    return dict(all_ordered_finite_deletion_pairs=rows, survivors=survivors, patches=patches,
                catalogue_diagrams=26, isolated_patch_forbidden_far_endpoint_projection='129/56',
                geometry_input='Minimal indefinite order at most10, finite deletion at both marks; branch/nonsimple depth at least6. The inherited all-label classification is not regenerated.')


def schur(ps):
    t, u, T, U = S.symbols('t u T U', nonnegative=True)
    rows, excluded, short, bad, constants = [], [], [], [], []

    def inner(x, y, n, shift):
        return S.expand(sum(v*w*P.inv(n, i+shift, j+shift)
                            for i, v in x.items() for j, w in y.items()))

    for ci, c in enumerate(ps):
        for di, original in enumerate(ps):
            d = {N-1-i: v for i, v in original.items()}
            paths = [(i, j) for i in c for j in d if abs(i-j)+2 <= 7
                     and not any(min(i, j) < k < max(i, j) for k in (A, B))]
            if paths:
                short.append(dict(c=ci, d=di, witness=paths[0]))
                continue
            x, y = {**c, A: t}, {**d, B: u}
            original_q, original_r = inner(x, x, N, 0), inner(y, y, N, 0)
            for tc, ud in product((0, 1), repeat=2):
                clow = P.field(original_q.subs(t, tc)-2)['sign'] <= 0
                dlow = P.field(original_r.subs(u, ud)-2)['sign'] <= 0
                if (clow and (tc or c != {1: 1})) or (dlow and (ud or d != {N-2: 1})):
                    excluded.append(dict(c=ci, d=di, t=tc, u=ud, c_low=clow, d_low=dlow))
                    continue
                for cp, dp in product(('extend', 'isolated') if clow else ('none',),
                                      ('extend', 'isolated') if dlow else ('none',)):
                    shift = int(cp == 'extend')
                    n = N+shift+int(dp == 'extend')
                    q = inner(x, x, n, shift)+S.Rational(int(cp == 'isolated'), 2)
                    r = inner(y, y, n, shift)+S.Rational(int(dp == 'isolated'), 2)
                    cross = inner(x, y, n, shift)
                    q0, r0 = P.field(q.subs(t, tc)-2), P.field(r.subs(u, ud)-2)
                    assert q0['sign'] > 0 and r0['sign'] > 0
                    polynomial = S.Poly(S.expand(((q-2)*(r-2)-cross**2).subs(
                        {t: 1+T if tc else 0, u: 1+U if ud else 0})), T, U)
                    coefficients = [dict(powers=list(m), value=P.field(v)) for m, v in polynomial.terms()]
                    constant = polynomial.eval({T: 0, U: 0})
                    row = dict(c=ci, d=di, t=tc, u=ud, c_patch=cp, d_patch=dp,
                               q_minus2_lower=q0, r_minus2_lower=r0,
                               coefficients=coefficients, constant=P.field(constant))
                    rows.append(row)
                    if P.field(constant)['sign'] <= 0 or any(z['value']['sign'] < 0 for z in coefficients):
                        bad.append(row)
                    else:
                        constants.append(constant)
    assert not bad, bad
    least = constants[0]
    for value in constants[1:]:
        if P.field(value-least)['sign'] < 0:
            least = value
    assert P.inv(N, 7, 7) == S.Rational(56, 11)
    assert P.inv(N, 7, 13) == S.Rational(32, 11)
    assert (P.inv(N, 7, 7)-2)**2-P.inv(N, 7, 13)**2 == S.Rational(12, 11)
    return dict(branch_count=len(rows), positive_branches=len(constants), bad_branches=bad,
                minimum_positive_constant=P.field(least), branches=rows,
                small_core_excluded=excluded, short_profile_pairs=short,
                central_singleton_determinant='12/11',
                arbitrary_opposite_products='Zero or1+arbitrary nonnegative parameter; all coefficients signed exactly in Q(sqrt2).')


def counts():
    rows = []
    for kind in 'ABD':
        for n in range(1 if kind == 'A' else 2 if kind == 'B' else 4, 81):
            edges = ([(i, i+1, 3) for i in range(n-1)] if kind != 'D' else
                     [(i, i+1, 3) for i in range(n-3)]+[(n-3, n-2, 3), (n-3, n-1, 3)])
            if kind == 'B':
                edges[-1] = (n-2, n-1, 4)
            graph = {i: [] for i in range(n)}
            for i, j, label in edges:
                graph[i].append(j)
                graph[j].append(i)
            features = {i for i in graph if len(graph[i]) != 2} | {
                i for a, b, label in edges if label != 3 for i in (a, b)}
            distances = []
            for i in graph:
                dd, todo = {i: 0}, deque([i])
                while todo:
                    v = todo.popleft()
                    for z in graph[v]:
                        if z not in dd:
                            dd[z] = dd[v]+1
                            todo.append(z)
                distances.append(dd)
            depths = [min(distances[i][f] for f in features) for i in graph]
            pairs = [(i, j) for i, j in combinations(range(n), 2)
                     if distances[i][j] == 8 and min(depths[i], depths[j]) >= 6]
            old = [(i, j) for i, j in pairs if min(depths[i], depths[j]) >= 8]
            assert len(pairs) == max(0, n-{'A': 20, 'B': 21, 'D': 22}[kind])
            assert len(old) == max(0, n-{'A': 24, 'B': 25, 'D': 26}[kind])
            rows.append(dict(kind=kind, rank=n, depth6_pairs=pairs, depth8_pairs=old,
                             hypothetical_additional_deletions=len(pairs)-len(old)))
    return dict(depth6_formulas=dict(A='max(0,n-20)', B='max(0,n-21)', D='max(0,n-22)'),
                depth8_formulas=dict(A='max(0,n-24)', B='max(0,n-25)', D='max(0,n-26)'), rows=rows)


def run():
    stems = (P.TAG, O.TAG, CAT, CURRENT, 'fresh_depth8_finite8_zero_transfer_20261018')
    docs = {s: json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for stem, doc in docs.items():
        assert doc['status'].startswith('PASS'), stem
        for path, digest in doc['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, (stem, path)
    ps, profile_record = profiles()
    print('PROFILES', len(ps), flush=True)
    core_record = small_cores(docs[CAT])
    schur_record = schur(ps)
    print('SCHUR', schur_record['branch_count'], 'minimum', schur_record['minimum_positive_constant'], flush=True)
    paths = {Path(__file__), Path(P.__file__), Path(O.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in stems
                 for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md')))
    return enc(dict(status='PASS_DEPTH6_FINITE8_ADAPTIVE_PILOT',
                    sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
                    profiles=profile_record, minimal_core_reduction=core_record,
                    schur=schur_record, counts=counts(), current_dimension_bound=58,
                    universal_payment_theorem=False, new_dimension_bound=None,
                    short7_and_cycle_composition_complete=False,
                    scope='Exact retained-profile/adaptive-Schur diagnostic conditional on inherited actual-star and minimal-core geometry. All branches are strict. Short-seven/cycle replacement and current-rule composition remain open; no local payment theorem or dimension improvement is asserted.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = ROOT/f'results/{TAG}.json'
    result = run()
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        assert not path.exists()
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], flush=True)
