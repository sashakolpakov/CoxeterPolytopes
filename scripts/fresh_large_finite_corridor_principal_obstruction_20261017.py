#!/usr/bin/env python3
"""Deep-corridor short7 alternatives using only a Lorentz principal core.

Old spectator blocks may join the actual corridor. Their roots are omitted
from the tested principal core; contacts to them need not be controlled.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import fresh_large_finite_corridor_tail_obstruction_20261017 as T

ROOT, H, A = T.ROOT, T.H, T.A
TAG = 'fresh_large_finite_corridor_principal_obstruction_20261017'


def run():
    parent, direct, previous = map(T.load, [T.S.TAG, T.DIRECT, T.TAG])
    pending = [r for r in direct['cases'] if r['status'] == 'unresolved_supporting_extensions']
    assert len(pending) == 251
    assert previous['selected_corners'] == 64
    previous_keys = {(r['parent_index'], tuple(r['marks'])) for r in previous['cases']}
    old_frames, principal_frames = {}, {}
    inventory, reasons = Counter(), Counter()
    contexts, rows = [], []
    for rem in pending:
        index = rem['parent_index']
        old = parent['cases'][index]
        core = parent['cores'][old['core_key']]
        if old['core_key'] not in old_frames:
            native = A.matrix(core['native_graph'])
            g, sc = H.native_to_rational(native)
            assert H.encode(g) == core['gram'] and sc == core['gauge']
            old_frames[old['core_key']] = native, g, sc
        native, g, sc = old_frames[old['core_key']]
        du = H.graph_dist(native, old['left_marks'])
        dv = H.graph_dist(native, old['right_marks'])
        assert min(du[x] for x in old['right_marks']) > 7
        context = dict(parent_index=index, core_key=old['core_key'], corners=[])
        for c in old['all_corners']:
            ns, marks = c['marked_component'], c['preferred_corner']
            path = A.simple_path(g, ns) if ns else None
            info = dict(marks=marks)
            if ns is None:
                why = 'marks_disconnected_in_displayed_corner'
            elif path is None:
                why = 'marked_component_not_simple_path'
            elif abs(path.index(marks[0])-path.index(marks[1])) != 8:
                why = 'displayed_mark_distance_not_eight'
            else:
                positions = sorted(path.index(x) for x in marks)
                near = [positions[0], len(path)-1-positions[1]]
                ids = sorted(ns + c['omitted_black'])
                info.update(path=path, old_endpoint_distances=near, principal_nodes=ids)
                if not all(1 <= d <= 7 for d in near):
                    why = 'endpoint_position_outside_supported_range'
                elif any(sc[i] != 1 for i in ids):
                    why = 'principal_core_not_physical_norm2_gauge'
                else:
                    frame_key = old['core_key'], tuple(ids)
                    if frame_key not in principal_frames:
                        gg = A.sub(g, ids)
                        sig = H.inertia(gg)
                        iv = H.inverse(gg) if sig == [len(ids)-1, 1, 0] else None
                        principal_frames[frame_key] = gg, sig, iv
                    gg, sig, iv = principal_frames[frame_key]
                    info['principal_signature'] = sig
                    why = ('eligible_lorentz_principal_core' if iv is not None
                           else 'principal_core_not_nonsingular_lorentz')
            inventory[why] += 1
            info['classification'] = why
            context['corners'].append(info)
            if why != 'eligible_lorentz_principal_core':
                continue

            assert H.inspect_corner(g, sc, old['left_marks']+old['right_marks'], *marks) == c
            omitted = c['omitted_black']
            assert len(omitted) == 2 and not set(omitted) & set(path)
            assert len(set(omitted) & set(old['left_marks'])) == 1
            assert len(set(omitted) & set(old['right_marks'])) == 1
            assert set(ids) == set(path) | set(omitted)
            pos = {x: k for k, x in enumerate(ids)}
            lengths = [T.TAIL_END-d for d in near]
            sides = []
            for endpoint, d, m in zip((path[0], path[-1]), near, lengths):
                supports = [[endpoint]] + ([sorted([endpoint, x]) for x in omitted] if d == 1 else [])
                side = []
                for ss in supports:
                    shortcut = min(min(du[i]+dv[k], dv[i]+du[k])+2 for i, k in product(ss, ss))
                    q = T.quadratic(iv, [pos[x] for x in ss], [pos[x] for x in ss])
                    side.append(dict(support=ss, projection=str(q),
                        diagonal=str(F(m+1, m)-q), shortcut_bound=shortcut))
                sides.append(side)
            profiles = []
            for left, right in product(*sides):
                cross = T.quadratic(iv, [pos[x] for x in left['support']], [pos[x] for x in right['support']])
                dl, dr = F(left['diagonal']), F(right['diagonal'])
                det = dl*dr-cross*cross
                if left['shortcut_bound'] <= 7 or right['shortcut_bound'] <= 7:
                    reason = 'short7_alternative_in_full_supporting_star'
                elif dl < 0 or dr < 0:
                    reason = 'negative_tail_schur_diagonal'
                elif det < 0:
                    reason = 'negative_tail_schur_determinant'
                else:
                    raise AssertionError(('unresolved', index, marks, left, right, det))
                reasons[reason] += 1
                profiles.append(dict(supports=[left['support'], right['support']],
                    diagonals=[str(dl), str(dr)], cross_projection=str(cross),
                    determinant=str(det), reason=reason))
            rows.append(dict(parent_index=index, core_key=old['core_key'], marks=marks,
                path=path, old_corner_blocks=c['blocks'], omitted_black=omitted,
                dropped_old_corner_nodes=sorted(set(range(len(g)))-set(ids)),
                principal_nodes=ids, principal_gram=H.encode(gg),
                principal_inverse=H.encode(iv), principal_signature=sig,
                full_old_core_had_radical=bool(core['radical']),
                old_endpoint_distances=near, forced_tail_lengths=lengths,
                side_profiles=sides, joint_profiles=profiles,
                inherited_one_block_case=(index, tuple(marks)) in previous_keys,
                conclusion='deep_corridor_forces_an_opposite_black_path_at_most7_in_the_full_actual_supporting_star'))
        contexts.append(context)
    keys = {(r['parent_index'], tuple(r['marks'])) for r in rows}
    assert previous_keys <= keys
    assert sum(inventory.values()) == 1004 and len(rows) == 164
    assert len({r['parent_index'] for r in rows}) == 134
    assert sum(r['inherited_one_block_case'] for r in rows) == 64
    assert sum(r['full_old_core_had_radical'] for r in rows) == 4
    paths = [Path(__file__), Path(A.__file__), Path(H.__file__),
             ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    paths += [ROOT/f'{folder}/{stem}.{ext}' for stem in (T.S.TAG, T.DIRECT, T.TAG)
              for folder, ext in [('scripts', 'py'), ('results', 'json'), ('results', 'md')]]
    return dict(status='PASS_LARGE_FINITE_CORRIDOR_PRINCIPAL_CORE_SHORT7_ALTERNATIVE',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        pending_contexts=251, inventoried_corners=1004, corner_inventory=dict(inventory),
        eligible_corners=len(rows), eligible_contexts=len({r['parent_index'] for r in rows}),
        inherited_one_block_corners=64, additional_corners=100,
        old_corner_block_counts=dict(Counter(str(len(r['old_corner_blocks'])) for r in rows)),
        full_radical_old_cores_with_lorentz_principal_subcore=4,
        selected_principal_frames=len({(r['core_key'], tuple(r['principal_nodes'])) for r in rows}),
        profile_outcomes=dict(reasons), cases=rows, complete_inventory=contexts,
        scope='On164 corners of134 residual contexts, a deep finite A/B/D corridor with both marks at feature-distance>=20 forces a full-star opposite-black short7 alternative. The64 one-block cases are included. In100 additional cases all old spectator components are dropped from the tested Gram; its nonsingular Lorentz principal core is just the marked path and omitted blacks. Forced tail roots may be old spectator whites and may contact dropped roots freely. Only their products against the principal core and one another are tested. No actual whole-component isolation, pair deletion, finite budget, global dimension improvement or full251-context resolution is asserted.')


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
    print(out['status'], out['corner_inventory'], out['profile_outcomes'], flush=True)
