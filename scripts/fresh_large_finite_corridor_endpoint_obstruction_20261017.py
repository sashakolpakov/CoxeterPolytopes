#!/usr/bin/env python3
"""Endpoint-mark deep corridors force full-star short7 alternatives.

All188 endpoint corners are inventoried. The116 physical-norm2 nonsingular
Lorentz principal cores admit complete seven/three entry domains and exact
at-most-four-by-four Q(sqrt(2)) tail Schur tests. Other cases stay open.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import fresh_large_finite_corridor_principal_obstruction_20261017 as P
import fresh_quadratic_sqrt2_exact_20261017 as Q

T, H, A, ROOT = P.T, P.H, P.A, P.ROOT
TAG = 'fresh_large_finite_corridor_endpoint_obstruction_20261017'


def endpoint_profiles(endpoint, distance, omitted):
    """Complete profiles in the full-star no-short7 branch.

    Two nonzero omitted contacts on one endpoint arm either make a short2/3
    path or violate strict edge PD on its two consecutive white roots.
    """
    m = 17-distance
    assert 0 <= distance <= 7
    if distance == 0:
        choices = [dict(kind='no_omitted_contact', vectors=[{endpoint: Q.q(1)}, {}])]
        for b in omitted:
            for label, weight in [(3, Q.q(1)), (4, Q.q(0, 1))]:
                choices.append(dict(kind='first_white_contact', omitted_black=b, label=label,
                    vectors=[{endpoint: Q.q(1), b: weight}, {}]))
            choices.append(dict(kind='second_white_contact', omitted_black=b, label=3,
                vectors=[{endpoint: Q.q(1)}, {b: Q.q(1)}]))
        base = [[Q.q(2), Q.q(-1)], [Q.q(-1), Q.q(F(m, m-1))]]
        assert len(choices) == 7
    else:
        choices = [dict(kind='no_omitted_contact', vectors=[{endpoint: Q.q(1)}])]
        if distance == 1:
            for b in omitted:
                choices.append(dict(kind='first_white_contact', omitted_black=b, label=3,
                    vectors=[{endpoint: Q.q(1), b: Q.q(1)}]))
        base = [[Q.q(F(m+1, m))]]
        assert len(choices) == (3 if distance == 1 else 1)
    return choices, base


def projected_product(inverse, positions, x, y):
    out = Q.q()
    for i, a in x.items():
        for j, b in y.items():
            out = Q.add(out, Q.mul(Q.mul(a, b), Q.q(inverse[positions[i]][positions[j]])))
    return out


def augmented_shortest_path(g, old, vectors, left_order):
    """An explicit old-core/retained-tail walk; no tail-root identity assumed."""
    n = len(g)
    adj = {i: {j for j in range(n) if i != j and g[i][j]} for i in range(n)}
    for k, v in enumerate(vectors):
        adj[n+k] = set(v)
        for i in v:
            adj[i].add(n+k)
    for start, size in [(0, left_order), (left_order, len(vectors)-left_order)]:
        for k in range(start, start+size-1):
            adj[n+k].add(n+k+1)
            adj[n+k+1].add(n+k)
    queue = list(old['left_marks'])
    previous = {i: None for i in queue}
    for i in queue:
        if i in old['right_marks']:
            path = []
            while i is not None:
                path.append(i)
                i = previous[i]
            return path[::-1]
        for k in sorted(adj[i]-previous.keys()):
            previous[k] = i
            queue.append(k)
    raise AssertionError('opposite marks disconnected in augmented graph')


def encoded_choice(choice):
    return {**{k: v for k, v in choice.items() if k != 'vectors'},
            'vectors': [{str(i): Q.encode(x) for i, x in v.items()} for v in choice['vectors']]}


def run():
    previous, parent = map(T.load, [P.TAG, T.S.TAG])
    selected = [(context, c) for context in previous['complete_inventory']
                for c in context['corners'] if c['classification'] == 'endpoint_position_outside_supported_range']
    assert len(selected) == 188
    assert all(0 in c['old_endpoint_distances'] for _, c in selected)
    helper_checks = Q.self_check()
    # All-label gates in physical norm2 gauge. At distance1, label5 has
    # squared magnitude >5/2; label4 is retained. At distance2 it is excluded.
    diag = [F(d*(21-d), 21) for d in range(1, 21)]
    assert diag[0]*F(5, 2) > 2 and diag[0]*2 < 2
    assert diag[1]*2 > 2 and all(diag[d-1] > 2 for d in range(3, 19))
    a2inverse = H.inverse([[F(2), F(-1)], [F(-1), F(2)]])
    assert all(x > 0 for row in a2inverse for x in row)
    assert sum(sum(row) for row in a2inverse) == 2
    frames, old_frames = {}, {}
    inventory, outcomes, sizes, negative_sizes = Counter(), Counter(), Counter(), Counter()
    rows = []
    for context, info in selected:
        index = context['parent_index']
        old = parent['cases'][index]
        core = parent['cores'][old['core_key']]
        if old['core_key'] not in old_frames:
            native = A.matrix(core['native_graph'])
            g, sc = H.native_to_rational(native)
            assert H.encode(g) == core['gram'] and sc == core['gauge']
            old_frames[old['core_key']] = native, g, sc
        native, g, sc = old_frames[old['core_key']]
        corner = H.inspect_corner(g, sc, old['left_marks']+old['right_marks'], *info['marks'])
        assert corner in old['all_corners']
        path = A.simple_path(g, corner['marked_component'])
        assert path == info['path']
        mark_positions = sorted(path.index(x) for x in info['marks'])
        assert mark_positions[1]-mark_positions[0] == 8
        near = [mark_positions[0], len(path)-1-mark_positions[1]]
        assert near == info['old_endpoint_distances'] and 0 in near
        assert all(0 <= d <= 7 for d in near)
        omitted = corner['omitted_black']
        ids = sorted(path+omitted)
        assert ids == info['principal_nodes'] and len(set(ids)) == len(ids)
        assert len(set(omitted) & set(old['left_marks'])) == 1
        assert len(set(omitted) & set(old['right_marks'])) == 1
        rec = dict(parent_index=index, core_key=old['core_key'], marks=info['marks'],
            path=path, old_endpoint_distances=near, omitted_black=omitted,
            principal_nodes=ids, principal_gauges=[sc[i] for i in ids],
            old_corner_blocks=corner['blocks'])
        if any(sc[i] != 1 for i in ids):
            rec['classification'] = 'nonunit_rational_gauge_outside_current_norm2_filter'
            inventory[rec['classification']] += 1
            rows.append(rec)
            continue
        frame_key = old['core_key'], tuple(ids)
        if frame_key not in frames:
            gg = A.sub(g, ids)
            sig = H.inertia(gg)
            inverse = H.inverse(gg) if sig == [len(ids)-1, 1, 0] else None
            frames[frame_key] = gg, sig, inverse
        gg, sig, inverse = frames[frame_key]
        rec.update(principal_signature=sig, principal_gram=H.encode(gg))
        if inverse is None:
            rec['classification'] = 'principal_core_not_nonsingular_lorentz'
            inventory['signature_'+str(sig)] += 1
            rows.append(rec)
            continue
        rec['classification'] = 'eligible_nonsingular_lorentz_norm2_core'
        inventory[rec['classification']] += 1
        rec['principal_inverse'] = H.encode(inverse)
        positions = {i: k for k, i in enumerate(ids)}
        left, lmat = endpoint_profiles(path[0], near[0], omitted)
        right, rmat = endpoint_profiles(path[-1], near[1], omitted)
        rec.update(side_profiles=[[encoded_choice(c) for c in side] for side in (left, right)],
            forced_tail_lengths=[17-d for d in near])
        profiles = []
        for li, lc in enumerate(left):
            for ri, rc in enumerate(right):
                vectors = lc['vectors']+rc['vectors']
                nl, n = len(lc['vectors']), len(vectors)
                assert 2 <= n <= 4
                sizes[n] += 1
                route = augmented_shortest_path(g, old, vectors, nl)
                profile = dict(side_choices=[li, ri], retained_tail_root_ids=list(range(len(g), len(g)+n)),
                    shortest_opposite_black_path=route)
                if len(route)-1 <= 7:
                    profile['outcome'] = 'full_star_short7_alternative'
                else:
                    schur = [[Q.q() for _ in range(n)] for _ in range(n)]
                    for offset, base in [(0, lmat), (nl, rmat)]:
                        for i, row in enumerate(base):
                            for j, v in enumerate(row):
                                schur[offset+i][offset+j] = v
                    for i in range(n):
                        for j in range(n):
                            schur[i][j] = Q.sub(schur[i][j],
                                projected_product(inverse, positions, vectors[i], vectors[j]))
                    witness, checked = Q.negative_principal_minor(schur)
                    assert witness is not None, ('PSD survivor', index, info['marks'], li, ri)
                    negative_sizes[len(witness['indices'])] += 1
                    profile.update(outcome='strict_negative_tail_schur_principal_minor',
                        schur_matrix=[[Q.encode(x) for x in row] for row in schur],
                        negative_minor=witness, principal_minors_checked=checked)
                outcomes[profile['outcome']] += 1
                profiles.append(profile)
        rec.update(joint_profiles=profiles,
            conclusion='deep_finite_corridor_forces_an_opposite_black_path_at_most7_in_the_full_actual_supporting_star')
        rows.append(rec)
    assert len(rows) == 188 and inventory['eligible_nonsingular_lorentz_norm2_core'] == 116
    assert outcomes == {'full_star_short7_alternative': 788, 'strict_negative_tail_schur_principal_minor': 640}
    paths = [Path(__file__), Path(Q.__file__), Path(A.__file__), Path(H.__file__),
             ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py',
             ROOT/'results/fresh_large_finite_corridor_endpoint_next_step_20261017.md']
    paths += [ROOT/f'{folder}/{stem}.{ext}' for stem in (P.TAG, T.TAG, T.S.TAG)
              for folder, ext in [('scripts', 'py'), ('results', 'json'), ('results', 'md')]]
    return dict(status='PASS_ENDPOINT_CORRIDOR_SHORT7_ALTERNATIVE_ON_COMPLETE_NORM2_LORENTZ_SUBDOMAIN',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        inherited_pending_contexts=251, endpoint_corners=188,
        endpoint_parent_contexts=len({r['parent_index'] for r in rows}),
        inventory=dict(inventory), eligible_corners=116,
        eligible_parent_contexts=len({r['parent_index'] for r in rows if 'joint_profiles' in r}),
        excluded_from_current_method_corners=72, joint_profiles=sum(outcomes.values()),
        profile_outcomes=dict(outcomes), retained_matrix_orders=dict(sizes),
        negative_minor_orders=dict(negative_sizes), exact_field='Q(sqrt(2)) in Fraction coefficient pairs',
        helper_independent_determinant_comparisons=helper_checks,
        all_label_gates=dict(A20_inverse_diagonals=list(map(str, diag)),
            distance1_labels=[3, 4], distance2_labels=[3],
            no_omitted_contacts_at_distances=list(range(3, 19)),
            A2_two_unit_contacts_projection_square='2'),
        cases=rows,
        scope='All188 endpoint-mark corners from the frozen164-corner principal inventory are reconstructed. Exactly116 pass the norm2-gauge and nonsingular-Lorentz principal-core filters; all1428 seven/three-profile combinations yield a full-star short7 alternative or a strict negative exact tail Schur principal minor. The72 gauge/signature cases remain open to different methods; they are not geometrically excluded. Deep-corridor impossibility is asserted only under an explicit full-actual-star no-short7 hypothesis. No pair payment change, global finite budget, dimension improvement or full parent-context resolution is asserted.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = json.loads(json.dumps(run()))
    target = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(out, separators=(',', ':'))+'\n')
    print(out['status'], out['inventory'], out['profile_outcomes'], flush=True)
