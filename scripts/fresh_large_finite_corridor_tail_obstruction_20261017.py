#!/usr/bin/env python3
"""Deep corridors force short7 alternatives on all 64 one-block open corners.

This is not a pair-payment deletion. The complete 251-context remaining
frontier is inventoried; only its nonsingular, one-block simple-path corners
with marked distance eight are treated by the exact tail argument.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import fresh_distance7_high_supporting_frontier_20261017 as S
import fresh_finite_a17_endpoint_payment_support_20261017 as A
import fresh_d16_d17_d14_coupon_forest_budget_20261017 as BUDGET

ROOT = Path(__file__).resolve().parents[1]
H = S.H
TAG = 'fresh_large_finite_corridor_tail_obstruction_20261017'
JOINT = 'fresh_distance7_joint_contact_domains_20261017'
DIRECT = 'fresh_distance7_direct_contact_features_20261017'
DEPTH, TAIL_END = 20, 17


def load(stem):
    d = json.loads((ROOT / f'results/{stem}.json').read_text())
    for p, h in d['sources'].items():
        assert sha256((ROOT / p).read_bytes()).hexdigest() == h, (stem, p)
    return d


def quadratic(j, x, y):
    return sum((j[a][b] for a in x for b in y), F())


def graph(n, kind):
    edges = [(i, i + 1, 3) for i in range(n - 1)]
    if kind == 'B':
        edges[0] = (0, 1, 4)
    if kind == 'D':
        edges = [(i, i + 1, 3) for i in range(n - 3)]
        edges += [(n - 3, n - 2, 3), (n - 3, n - 1, 3)]
    adj = {i: set() for i in range(n)}
    features = set()
    for a, b, m in edges:
        adj[a].add(b)
        adj[b].add(a)
        if m != 3:
            features.update((a, b))
    features.update(i for i in adj if len(adj[i]) != 2)
    return adj, features


def distances(adj, starts):
    d = {i: 0 for i in starts}
    queue = list(starts)
    for i in queue:
        for k in sorted(adj[i] - d.keys()):
            d[k] = d[i] + 1
            queue.append(k)
    return d


def candidate_counts():
    """Only conditional costs; this checker does not authorize the change."""
    rows = []
    for n in range(18, 76):
        for kind, offset in [('A', 48), ('B', 49), ('D', 50)]:
            adj, features = graph(n, kind)
            fd = distances(adj, features)
            pairs = [(a, b) for a in range(n) for b in range(a + 1, n)
                     if fd[a] >= DEPTH and fd[b] >= DEPTH
                     and distances(adj, [a])[b] == 8]
            assert len(pairs) == max(0, n - offset)
            base = 8 * n - (32 if kind == 'D' and n < 20 else 33 if kind == 'D' else 34)
            assert BUDGET.base((kind, n)) == 2 * base
            rows.append(dict(rank=n, kind=kind, candidate_pairs=pairs,
                current_same_component_cost=base,
                hypothetical_half_weight_cost=str(F(base) - F(len(pairs), 2))))
    return rows


def run():
    parent, joint, direct = map(load, [S.TAG, JOINT, DIRECT])
    load(BUDGET.TAG)  # Bind the current cost formula, without replaying its DP.
    old_pending = {r['parent_index']: r for r in joint['cases']
                   if r['status'] == 'unresolved_supporting_extensions'}
    remaining = {r['parent_index']: r for r in direct['cases']
                 if r['status'] == 'unresolved_supporting_extensions'}
    assert len(old_pending) == 548 and len(remaining) == 251
    assert set(remaining) <= old_pending.keys()

    # Physical norm-two A20 inverse. Its displayed nodes are all supporting
    # whites, at distances 1,...,20 outward from one marked black.
    a20 = [[F(2 if i == k else -1 if abs(i-k) == 1 else 0)
            for k in range(DEPTH)] for i in range(DEPTH)]
    j20 = H.inverse(a20)
    assert all(x > 0 for row in j20 for x in row)
    diagonals = [F(d * (DEPTH + 1 - d), DEPTH + 1) for d in range(1, DEPTH + 1)]
    assert [j20[i][i] for i in range(DEPTH)] == diagonals
    assert all(diagonals[d-1] > 2 for d in range(3, 19))
    assert diagonals[1] > 1 and 2 * diagonals[1] > 2

    frames = {}
    inventory = Counter()
    excluded = Counter()
    records, contexts = [], []
    for index, pending in sorted(remaining.items()):
        old = parent['cases'][index]
        core = parent['cores'][old['core_key']]
        g = A.matrix(core['gram'])
        context = dict(parent_index=index, core_key=old['core_key'],
                       base_context=old['base_context'], corners=[])
        for c in old['all_corners']:
            ns = c['marked_component']
            path = A.simple_path(g, ns) if ns else None
            if core['radical']:
                why = 'radical_core_outside_scope'
            elif ns is None:
                why = 'marks_disconnected_in_displayed_corner'
            elif len(c['blocks']) != 1:
                why = 'old_spectator_blocks_outside_scope'
            elif path is None:
                why = 'displayed_marked_component_not_simple_path'
            elif abs(path.index(c['preferred_corner'][0]) - path.index(c['preferred_corner'][1])) != 8:
                why = 'displayed_mark_distance_not_eight'
            else:
                why = 'eligible_one_block_simple_path'
            inventory[why] += 1
            context['corners'].append(dict(marks=c['preferred_corner'], classification=why))
            if why != 'eligible_one_block_simple_path':
                continue

            key = old['core_key']
            if key not in frames:
                native = A.matrix(core['native_graph'])
                gg, sc = H.native_to_rational(native)
                assert gg == g and sc == core['gauge'] and set(sc) == {1}
                assert H.inertia(g) == core['signature'] == [len(g)-1, 1, 0]
                iv = H.inverse(g)
                assert H.encode(iv) == core['inverse']
                frames[key] = native, iv
            native, iv = frames[key]
            black = old['left_marks'] + old['right_marks']
            a, b = c['preferred_corner']
            assert H.inspect_corner(g, core['gauge'], black, a, b) == c
            omitted = c['omitted_black']
            assert len(omitted) == 2
            assert len(set(omitted) & set(old['left_marks'])) == 1
            assert len(set(omitted) & set(old['right_marks'])) == 1
            # This is the essential no-old-spectator guard, not just connected
            # marks. Every old corner root belongs to this single path.
            assert len(c['blocks']) == 1
            assert sorted(path) == sorted(set(range(len(g))) - set(omitted))
            mark_positions = sorted(path.index(x) for x in (a, b))
            assert mark_positions[1] - mark_positions[0] == 8
            near = [mark_positions[0], len(path)-1-mark_positions[1]]
            assert all(1 <= d <= 7 for d in near)
            lengths = [TAIL_END-d for d in near]
            assert all(m >= 10 for m in lengths)
            du = H.graph_dist(native, old['left_marks'])
            dv = H.graph_dist(native, old['right_marks'])
            assert min(du[x] for x in old['right_marks']) > S.THRESHOLD == 7

            sides = []
            for endpoint, d, m in zip((path[0], path[-1]), near, lengths):
                supports = [[endpoint]]
                # At d=1 the first new root has outward distance2. Any old
                # omitted-black contact there is simple; both omitted blacks
                # at once would give the opposite-black path of length2.
                if d == 1:
                    supports += [sorted([endpoint, x]) for x in omitted]
                side = []
                for ss in supports:
                    shortcut = min((min(du[i]+dv[k], dv[i]+du[k])+2
                                    for i, k in product(ss, ss)), default=999)
                    q = quadratic(iv, ss, ss)
                    diag = F(m+1, m) - q
                    side.append(dict(support=ss, core_projection=str(q),
                        tail_schur_diagonal=str(diag), shortcut_bound=shortcut,
                        short7_excluded=shortcut <= 7))
                sides.append(side)

            profiles = []
            for left, right in product(*sides):
                cross = quadratic(iv, left['support'], right['support'])
                dl, dr = F(left['tail_schur_diagonal']), F(right['tail_schur_diagonal'])
                det = dl * dr - cross * cross
                if left['short7_excluded'] or right['short7_excluded']:
                    reason = 'single_entry_short7'
                elif dl < 0 or dr < 0:
                    reason = 'negative_tail_schur_diagonal'
                elif det < 0:
                    reason = 'negative_tail_schur_determinant'
                else:
                    raise AssertionError(('unresolved_corridor_profile', index, (a, b), left, right))
                excluded[reason] += 1
                profiles.append(dict(supports=[left['support'], right['support']],
                    diagonals=[str(dl), str(dr)], cross_projection=str(cross),
                    determinant=str(det), exclusion=reason))
            records.append(dict(parent_index=index, core_key=key, marks=[a, b],
                all_displayed_corner_blocks=c['blocks'], omitted_black=omitted,
                path=path, marked_path_positions=mark_positions,
                old_endpoint_distances=near, forced_tail_lengths=lengths,
                side_profiles=sides, all_joint_profiles=profiles,
                conclusion='deep_corridor_implies_an_opposite_black_short7_path_in_the_full_supporting_star'))
        contexts.append(context)
    assert sum(inventory.values()) == 4 * 251
    assert len(records) == inventory['eligible_one_block_simple_path'] == 64
    assert len({r['parent_index'] for r in records}) == 64
    paths = [Path(__file__), Path(A.__file__), Path(H.__file__),
             ROOT / 'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    paths += [Path(BUDGET.Q.__file__), Path(BUDGET.B.__file__)]
    paths += [ROOT/f'{folder}/{stem}.{ext}' for stem in (S.TAG, JOINT, DIRECT, BUDGET.TAG)
              for folder, ext in [('scripts', 'py'), ('results', 'json'), ('results', 'md')]]
    return dict(status='PASS_CONDITIONAL_LARGE_FINITE_CORRIDOR_TAIL_OBSTRUCTION',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        current_pending_contexts=251, inventoried_corners=1004,
        corner_inventory=dict(inventory), selected_corners=len(records),
        selected_contexts=len(records), reconstructed_nonsingular_cores=len(frames),
        native_feature_distance=DEPTH, forced_white_path_order=DEPTH,
        white_path_inverse_diagonals=list(map(str, diagonals)),
        omitted_black_contact_forbidden_distances=list(range(3, 19)),
        first_entry_distance2_maximum_label=3, last_used_tail_distance=TAIL_END,
        no_short7_branch_means='no opposite-black path of length<=7 anywhere in the full actual supporting star, not merely in the saved old core',
        profile_exclusions=dict(excluded), cases=records, complete_pending_inventory=contexts,
        conditional_all_rank_pair_counts=candidate_counts(),
        scope='All1004 corners of the251 still-open contexts are inventoried. At64 nonsingular corners whose entire displayed corner is one simple path with marked distance8, deep classical finite A/B/D corridors with both feature distances>=20 force an opposite-black path<=7 in the full actual supporting star. Equivalently such corridors are excluded only in the full-star no-short7 branch. Forced A20 white arms give all-label omitted-black restrictions; exact Lorentz-complement tail Schur tests exclude every profile not already supplying that shortcut. Other251-context corners, previously certified corner payment maps, low/disjoint and later payment interfaces remain outside this theorem. Candidate half-costs are conditional graph counts, not a new pair rule, finite budget, global bound or construction.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = json.loads(json.dumps(run()))
    target = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(out, separators=(',', ':')) + '\n')
    print(out['status'], out['corner_inventory'], out['profile_exclusions'], flush=True)
