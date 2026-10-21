#!/usr/bin/env python3
"""Bounded numerical compatibility/corner/global-gap diagnostics of a saved bank.

No new event search and no exact acceptance claims. The bank is read-only.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import networkx as nx
import fresh_h20_nonright_s61_20261014_scout as SOURCE

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_nonright_s61_20261014_selection'
BANK = 'results/fresh_h20_nonright_s61_20261014_scout.json'
BASE = 'results/fresh_h20_s61_multicap_20261014_exact.json'
TOL = 1e-10


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def relation(g):
    if g < -2 - TOL:
        return 'dotted'
    if abs(g + 2) <= TOL:
        return 'parallel'
    if abs(g) <= TOL:
        return 'right'
    if g > 0:
        return 'obtuse'
    m = int(round(np.pi / np.arccos(-g / 2)))
    if m >= 3 and abs(g + 2 * np.cos(np.pi / m)) <= TOL:
        return f'pi/{m}'
    return 'forbidden_numeric'


def run():
    bank = json.loads((ROOT / BANK).read_text())
    base = json.loads((ROOT / BASE).read_text())
    for doc in [bank, base]:
        for p, h in doc['sources'].items():
            assert sha(p) == h, p
    proposals = bank['proposals']
    Jexact = [[F(x) for x in row] for row in base['support_inverse']]
    J = np.array(Jexact, float)
    C = np.array([r['contacts'] for r in proposals])
    Y = np.array([r['residual_frame'] for r in proposals])
    newG = C @ J @ C.T + (Y * [-1, 1, 1]) @ Y.T
    norm_error = float(np.max(np.abs(np.diag(newG) - 2)))
    assert norm_error < 1e-10
    graph = nx.Graph()
    graph.add_nodes_from(range(len(proposals)))
    counts = Counter()
    finite = []
    for i, j in combinations(range(len(proposals)), 2):
        kind = relation(float(newG[i, j]))
        counts[kind] += 1
        if kind not in ['obtuse', 'forbidden_numeric']:
            graph.add_edge(i, j)
        if kind not in ['obtuse', 'forbidden_numeric', 'dotted', 'parallel']:
            finite.append(dict(indices=[i, j], relation=kind, pairing=float(newG[i, j])))
    assert [r['indices'] for r in finite] == [[103, 811], [104, 812], [105, 813]]
    assert all(r['relation'] == 'right' for r in finite)
    assert all(set(graph.neighbors(i)) == {j} for i, j in [r['indices'] for r in finite])
    assert all(graph.degree(i) == 0 for i in [44, 45, 63])
    cliques = list(nx.find_cliques(graph))
    clique_sizes = Counter(map(len, cliques))
    best = max(cliques, key=lambda x: (len(x), sum(proposals[i]['escape_pair'] for i in x)))
    print('PAIR_GRAPH', dict(counts), 'maximal_cliques', dict(clique_sizes), flush=True)

    # Exact rational cancellation for the common-residual right-pair template.
    ca, cb = [list(map(lambda x: F(str(x)), proposals[i]['contacts'])) for i in [103, 811]]
    def jc(a, b):
        return sum(a[i] * Jexact[i][j] * b[j] for i in range(18) for j in range(18))
    qa, qb, cross = jc(ca, ca), jc(cb, cb), jc(ca, cb)
    assert qa == qb == F(46, 11) and cross == F(24, 11)
    assert cross + 2 - qa == 0
    assert all(proposals[i]['anchors'] == proposals[j]['anchors'] == [65, 95]
               and proposals[i]['anchor_labels'] == proposals[j]['anchor_labels']
               and min(proposals[i]['anchor_labels']) >= 5 for i, j in [r['indices'] for r in finite])

    # Reconstruct the old metric from the high-precision physical loader,
    # then use a well-conditioned orthonormal Euclidean/Lorentz chart.
    d = SOURCE.load()
    P = np.array(d['P'].tolist(), float)
    W = np.array(d['W'].tolist(), float)
    Gs = np.linalg.inv(J)
    L = np.linalg.cholesky(Gs)
    old_spatial = (L.T @ W).T
    oldA = np.column_stack([old_spatial, P[:, 1:]])
    oldb = -P[:, 0]
    newA = np.column_stack([-np.linalg.solve(L, C.T).T, Y[:, 1:]])
    newb = Y[:, 0]
    assert np.max(np.abs(np.einsum('ij,ij->i', newA, newA) - newb**2 - 2)) < 1e-9
    finite_corners = []
    for record in finite:
        ia, ib = record['indices']
        common = [j for j in range(133) if proposals[ia]['relations'][j] not in ['dotted', 'parallel']
                  and proposals[ib]['relations'][j] not in ['dotted', 'parallel']]
        # Evaluate only the small common-neighbor Gram, at the loader's100 digits.
        G = np.array([[float(d['dot'](d['N'][i], d['N'][j])) for j in common]
                      for i in common])
        checked = 0
        hits = []
        for ix in combinations(range(len(common)), 18):
            H0 = G[np.ix_(ix, ix)]
            if np.any(H0[np.triu_indices(18, 1)] <= -2 + 1e-9):
                continue
            checked += 1
            ids = [common[j] for j in ix]
            H = np.zeros((20, 20))
            H[:18, :18] = H0
            for j, cap in enumerate([ia, ib], 18):
                vals = [proposals[cap]['all_pairs'][i] for i in ids]
                H[j, :18] = H[:18, j] = vals
                H[j, j] = 2
            H[18, 19] = H[19, 18] = newG[ia, ib]
            eig = np.linalg.eigvalsh(H)
            if eig[0] > 1e-9:
                hits.append(dict(old_support=ids, minimum_eigenvalue=float(eig[0])))
        finite_corners.append(dict(indices=[ia, ib], common_old_finite_neighbors=common,
                                   nondotted_18_subsets_tested=checked, positive_definite_supports=hits))
    print('FINITE_CORNER_TESTS', [(x['indices'], len(x['positive_definite_supports'])) for x in finite_corners], flush=True)

    # Ignore mutual cap compatibility: if even the whole intersection has a
    # strict spacelike point, every subset has that same numerical global gap.
    A = np.vstack([oldA, newA])
    b = np.r_[oldb, newb]
    row_scale = np.linalg.norm(A, axis=1)
    An, bn = A / row_scale[:, None], b / row_scale
    center = linprog(np.r_[np.zeros(20), -1], A_ub=np.column_stack([An, np.ones(len(A))]),
                     b_ub=bn, bounds=[(None, None)] * 20 + [(0, None)], method='highs')
    assert center.success and center.x[-1] > 1e-8
    h = center.x[:20]
    interior_norm = float(h @ h - 1)
    rng = np.random.default_rng(20261014)
    directions = list(np.eye(20)) + list(-np.eye(20))
    directions += [x / np.linalg.norm(x) for x in rng.normal(size=(160, 20))]
    directions += [newA[i] / np.linalg.norm(newA[i]) for i in [44, 45, 63, 103, 104, 105, 811, 812, 813]]
    best_gap = None
    statuses = Counter()
    for u in directions:
        result = linprog(-u, A_ub=An, b_ub=bn, bounds=[(None, None)] * 20, method='highs')
        statuses[result.status] += 1
        if result.status == 3:
            result = linprog(-u, A_ub=An, b_ub=bn, bounds=[(-100, 100)] * 20, method='highs')
        if not result.success:
            continue
        z = (1 - 1e-5) * result.x + 1e-5 * h
        q = float(z @ z - 1)
        margin = float(np.min(bn - An @ z))
        if q > 1e-6 and margin > 1e-9 and (best_gap is None or q > best_gap['norm']):
            best_gap = dict(spatial_chart=z.tolist(), time_coordinate=1,
                            norm=q, minimum_scaled_margin=margin,
                            maximum_unscaled_pairing=float(np.max(A @ z - b)))
    if best_gap is not None:
        # Physical reconstruction is only a numerical proposal for later exactification.
        z = np.array(best_gap['spatial_chart'])
        support_coeff = np.linalg.solve(L.T, z[:18])
        physical = sum((float(a) * np.array(list(d['N'][i]), float) for a, i in zip(support_coeff, d['S'])), np.zeros(21))
        physical += np.array(list(d['frame'][0]), float)
        physical += z[18] * np.array(list(d['frame'][1]), float) + z[19] * np.array(list(d['frame'][2]), float)
        best_gap['physical_point_double'] = physical.tolist()
    print('WHOLE_BANK_GAP', None if best_gap is None else (best_gap['norm'], best_gap['minimum_scaled_margin']), flush=True)
    paths = [BANK, BASE, f'scripts/{TAG}.py', 'scripts/fresh_h20_nonright_s61_20261014_scout.py',
             'scripts/fresh_h20_s61_multicap_20261014_scout.py', 'scripts/fresh_h20_coupled_roles_20261014.py']
    return dict(status='BOUNDED_NUMERICAL_847_BANK_SELECTION_AND_GLOBAL_GAP_DIAGNOSTIC',
                sources={p: sha(p) for p in paths}, proposal_count=len(proposals),
                pair_tolerance=TOL, maximum_norm_error=norm_error, pair_counts=dict(counts),
                robust_finite_pairs=finite, finite_pair_no_third_member=True,
                isolated_priority_candidates=[44, 45, 63], maximal_clique_size_counts=dict(clique_sizes),
                largest_compatible_clique=best, finite_corner_screens=finite_corners,
                rational_common_residual_identity=dict(first_cost=str(qa), second_cost=str(qb), cross_cost=str(cross),
                    consequence='For exactly equal residuals of square2−46/11, the two cap normals have exactly zero mutual pairing. Equality of the chosen physical anchor recipes requires a separate exact check.'),
                conditional_no_PSD20='If the cap formulas and complete common finite-neighbor lists are exactly certified, no PSD20 support can contain either right pair. Including65or95 gives a3-node star with mutually right caps and both equal anchor labels>=5, hence Schur value2−4cos(pi/m)^2<0. Omitting both anchors forces all18S roots, but S+either cap has residual−24/11. This is conditional algebra, not exact certification of the847-proposal bank.',
                whole_bank_relaxation=dict(LP_queries=len(directions), LP_statuses=dict(statuses),
                    common_interior_scaled_margin=float(center.x[-1]), common_interior_chart=h.tolist(),
                    common_interior_norm=interior_norm, strict_spacelike_proposal=best_gap),
                scope='Numerical diagnostic only. The original847 normals remain proposals. Gram tolerances, LP witnesses, finite support tests and all-label nearest-angle fits are not exact acceptance or exclusion certificates. The rational J-cost identity is exact conditional algebra. No event bank expansion or primary-data edits.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    data = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        # Only stable discrete outcomes and input manifests are replayed; LP
        # floating point coordinates are diagnostic, not bitwise certificates.
        saved = json.loads(path.read_text())
        for key in ['sources', 'pair_counts', 'finite_corner_screens',
                    'rational_common_residual_identity']:
            assert data[key] == saved[key], key
        assert [(r['indices'], r['relation']) for r in data['robust_finite_pairs']] == [
            (r['indices'], r['relation']) for r in saved['robust_finite_pairs']]
        gap = data['whole_bank_relaxation']['strict_spacelike_proposal']
        assert gap is not None and gap['norm'] > 1 and gap['minimum_scaled_margin'] > 1e-8
        assert data['whole_bank_relaxation']['common_interior_norm'] < 0
    else:
        path.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], flush=True)
