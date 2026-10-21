#!/usr/bin/env python3
"""Full A39 completion excludes deep corridors in the no-short7 high domain.

All 26,812 corners of the 6,703 frozen high-overlap contexts are covered.
The conclusion is a full-star short7 alternative, not a payment change.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import fresh_large_finite_corridor_tail_obstruction_20261017 as T
import fresh_quadratic_sqrt2_exact_20261017 as Q

ROOT, H, A = T.ROOT, T.H, T.A
TAG = 'fresh_full_a39_corridor_high_domain_20261017'
N, LEFT, RIGHT = 39, 15, 23
DEPTH, RETAINED_OUTWARD = 17, 15


def physical(g, gauge, i, j):
    """Undo the frozen positive diagonal congruence to physical norm two."""
    value, square = g[i][j], gauge[i]*gauge[j]
    assert square in (1, 2, 4)
    return Q.q(value) if square == 1 else Q.q(value/2) if square == 4 else Q.q(0, value/2)


def inverse_entry(n, i, j):
    return F((min(i, j)+1)*(n-max(i, j)), n+1)


@lru_cache(None)
def product_value(x, y, white=False):
    answer = Q.q()
    for i, a in x:
        if white and i in (LEFT, RIGHT):
            continue
        for j, b in y:
            if white and j in (LEFT, RIGHT):
                continue
            if white:
                block = next(((lo, hi) for lo, hi in ((0, 14), (16, 22), (24, 38))
                              if lo <= i <= hi and lo <= j <= hi), None)
                if block is None:
                    continue
                lo, hi = block
                value = inverse_entry(hi-lo+1, i-lo, j-lo)
            else:
                value = inverse_entry(N, i, j)
            answer = Q.add(answer, Q.mul(Q.mul(a, b), Q.q(value)))
    return answer


def nonsimple_witness(g, nodes, mark):
    adj = {i: {j for j in nodes if i != j and g[i][j]} for i in nodes}
    distances = T.distances(adj, [mark])
    assert len(distances) == len(nodes) <= 17
    assert max(distances.values()) <= 16 < DEPTH
    for i in nodes:
        if len(adj[i]) >= 3:
            return dict(reason='persistent_branch', node=i, marked_distance=distances[i])
        for j in sorted(adj[i]):
            if 4*g[i][j]**2 != g[i][i]*g[j][j]:
                return dict(reason='persistent_nonsimple_bond', bond=[i, j],
                            marked_distance=min(distances[i], distances[j]))
    # A connected simple-edge graph of maximum degree two is a path or cycle.
    assert all(len(adj[i]) == 2 for i in nodes)
    return dict(reason='cycle_cannot_embed_in_finite_classical_tree', nodes=nodes)


def embeddings(paths, marks, connected):
    if connected:
        path = paths[0]
        assert len(path) <= 17
        distance = abs(path.index(marks[0])-path.index(marks[1]))
        if distance != 8:
            return [], dict(reason='wrong_fixed_mark_distance', distance=distance)
        assert len(path)-9 <= 8 < RETAINED_OUTWARD
        if path.index(marks[0]) > path.index(marks[1]):
            path = path[::-1]
        return [{i: LEFT-path.index(marks[0])+j for j, i in enumerate(path)}], None
    assert sum(map(len, paths)) <= 17 and all(1 <= len(p) <= 16 for p in paths)
    assert all(len(p)-1 <= RETAINED_OUTWARD for p in paths)
    choices = [[p, p[::-1]] if len(p) > 1 else [p] for p in paths]
    out, trials = [], []
    for p, r in product(*choices):
        a, b = p.index(marks[0]), r.index(marks[1])
        inward = [len(p)-1-a, len(r)-1-b]
        gap = 7-sum(inward)
        trials.append(dict(paths=[p, r], inward_distances=inward, inserted_whites=gap))
        if gap >= 1:
            out.append(dict([(i, LEFT-a+j) for j, i in enumerate(p)] +
                            [(i, RIGHT+b-j) for j, i in enumerate(r)]))
    return out, None if out else dict(reason='no_distance8_embedding', orientation_trials=trials)


def constraints(g, gauge, mapping, omitted, is_left):
    mark = LEFT if is_left else RIGHT
    slots = [mark, mark-1, mark-2, mark+1] if is_left else [mark, mark+1, mark+2, mark-1]
    fixed = {position: Q.neg(physical(g, gauge, omitted, node)) for node, position in mapping.items()}
    assert mark in fixed
    for position, value in sorted(fixed.items()):
        assert Q.sign(value) >= 0
        if value == Q.q() or position in slots:
            continue
        reflected = position if is_left else N-1-position
        if reflected <= LEFT-3 or reflected >= RIGHT+3:
            distance = LEFT-reflected if reflected <= LEFT-3 else reflected-RIGHT
            assert 3 <= distance <= RETAINED_OUTWARD
            diagonal = inverse_entry(DEPTH, distance-1, distance-1)
            assert diagonal > 2
            reason = dict(reason='A17_strict_edge_contact_exclusion', outward_distance=distance,
                          inverse_diagonal=str(diagonal))
        else:
            distance = 1+abs(reflected-RIGHT)
            assert distance <= 7
            reason = dict(reason='full_star_short7_from_fixed_contact', route_length=distance)
        return None, dict(omitted_black=omitted, path_position=position,
                          physical_contact_magnitude=Q.encode(value), **reason)
    for k, position in enumerate(slots[1:]):
        if position not in fixed:
            continue
        value = fixed[position]
        domain = [Q.q(), Q.q(1)] + ([] if k == 1 else [Q.q(0, 1)])
        if value not in domain:
            diagonal = [F(17, 18), F(16, 9), F(7, 8)][k]
            projection = Q.mul(Q.mul(value, value), Q.q(diagonal))
            assert Q.sign(Q.sub(projection, Q.q(2))) >= 0
            return None, dict(reason='strict_actual_edge_label_bound', omitted_black=omitted,
                path_position=position, physical_contact_magnitude=Q.encode(value),
                inverse_diagonal=str(diagonal), projection_square=Q.encode(projection))
    return tuple(fixed.get(i) for i in slots), None


def side_profiles(fixed, is_left):
    mark = LEFT if is_left else RIGHT
    slots = [mark-1, mark-2, mark+1] if is_left else [mark+1, mark+2, mark-1]
    choices = []
    for k, value in enumerate(fixed[1:]):
        domain = [Q.q(), Q.q(1)] + ([] if k == 1 else [Q.q(0, 1)])
        assert value is None or value in domain
        choices.append(domain if value is None else [value])
    rows, rejected_consecutive_contacts = [], 0
    for values in product(*choices):
        if values[0] != Q.q() and values[1] != Q.q():
            rejected_consecutive_contacts += 1
            continue  # A2 positive inverse gives projection >=2; edge must be strict PD.
        vector = {mark: fixed[0]}
        vector.update(zip(slots, values))
        vector = tuple(sorted((i, value) for i, value in vector.items() if value != Q.q()))
        edge_slack = Q.sub(Q.q(2), product_value(vector, vector, True))
        rows.append((vector, edge_slack))
    return rows, rejected_consecutive_contacts


def analyse_template(fixed_pair):
    left, left_rejected = side_profiles(fixed_pair[0], True)
    right, right_rejected = side_profiles(fixed_pair[1], False)
    assert left and right
    records, outcomes = [], Counter()
    for li, (x, xe) in enumerate(left):
        for ri, (y, ye) in enumerate(right):
            rec = dict(side_choices=[li, ri], strict_edge_slacks=[Q.encode(xe), Q.encode(ye)])
            if Q.sign(xe) <= 0 or Q.sign(ye) <= 0:
                reason = 'strict_actual_edge_failure'
            else:
                a = Q.sub(Q.q(2), product_value(x, x))
                b = Q.sub(Q.q(2), product_value(y, y))
                cross = Q.neg(product_value(x, y))  # The two omitted blacks are orthogonal.
                determinant = Q.sub(Q.mul(a, b), Q.mul(cross, cross))
                assert Q.sign(a) < 0 and Q.sign(determinant) > 0, ('survivor', fixed_pair, x, y)
                reason = 'two_negative_directions'
                rec.update(schur_diagonals=[Q.encode(a), Q.encode(b)],
                           schur_cross=Q.encode(cross), determinant=Q.encode(determinant))
            rec['reason'] = reason
            outcomes[reason] += 1
            records.append(rec)
    return dict(fixed_slot_constraints=[[None if v is None else Q.encode(v) for v in row] for row in fixed_pair],
        side_profiles=[[[[i, Q.encode(v)] for i, v in vector] for vector, _ in side] for side in (left, right)],
        consecutive_outer_contact_exclusions=[left_rejected, right_rejected],
        profile_outcomes=dict(outcomes), joint_profiles=records)


def rank65_diagnostic():
    """Counts only: neither this theorem nor these numbers change a payment."""
    rows=[]
    for kind,base in [('A',486),('B',486),('D',487)]:
        adj,features=T.graph(65,kind);fd=T.distances(adj,features)
        eligible=[(i,j)for i in range(65)for j in range(i+1,65)
                  if fd[i]>=DEPTH and fd[j]>=DEPTH and T.distances(adj,[i])[j]==8]
        assert len(eligible)=={'A':23,'B':22,'D':21}[kind]
        rows.append(dict(kind=kind,rank=65,feature_nodes=sorted(features),eligible_pairs=eligible,
            eligible_count=len(eligible),inherited_same_component_cost=base,
            hypothetical_half_weight_cost=str(F(base)-F(len(eligible),2))))
    return rows


def run():
    parent = T.load(T.S.TAG)
    T.load('fresh_full_a43_corridor_high_domain_20261017')
    assert len(parent['cases']) == 6703 and len(parent['cores']) == 342
    frames = {}
    for key, core in parent['cores'].items():
        native = A.matrix(core['native_graph'])
        g, gauge = H.native_to_rational(native)
        assert H.encode(g) == core['gram'] and gauge == core['gauge']
        assert len(g) <= 19 and set(gauge) <= {1, 2}
        frames[key] = native, g, gauge
    # Independently verify the closed inverse formula used throughout.
    for n in (2, 7, 15, 17, 39):
        for i in range(n):
            for j in range(n):
                value = 2*inverse_entry(n, i, j)
                if i:
                    value -= inverse_entry(n, i-1, j)
                if i+1 < n:
                    value -= inverse_entry(n, i+1, j)
                assert value == int(i == j)
    assert inverse_entry(17, 0, 0)*F(5, 2) > 2
    assert inverse_entry(7, 0, 0)*F(5, 2) > 2
    assert inverse_entry(17, 1, 1)*2 > 2
    assert all(inverse_entry(17, d-1, d-1) > 2 for d in range(3, 16))
    assert sum(inverse_entry(2, i, j) for i in range(2) for j in range(2)) == 2
    templates, template_ids, uses = {}, {}, Counter()
    counts, profiles, rows = Counter(), Counter(), []
    orientation_count = 0
    for index, old in enumerate(parent['cases']):
        native, g, gauge = frames[old['core_key']]
        black = old['left_marks']+old['right_marks']
        assert len(set(black)) == 4
        dist = H.graph_dist(native, old['left_marks'])
        assert min(dist[i] for i in old['right_marks']) > 7
        assert {tuple(c['preferred_corner']) for c in old['all_corners']} == set(product(old['left_marks'], old['right_marks']))
        for c in old['all_corners']:
            marks = c['preferred_corner']
            omitted = [next(b for b in side if b != m)
                       for side, m in zip((old['left_marks'], old['right_marks']), marks)]
            assert set(omitted) == set(c['omitted_black']) and physical(g, gauge, *omitted) == Q.q()
            parts = A.components(g, sorted(set(range(len(g)))-set(omitted)))
            assert parts == [b['nodes'] for b in c['blocks']]
            assert sum(map(len,parts)) == len(g)-2 <= 17
            blocks = [next(b for b in parts if m in b) for m in marks]
            paths = [A.simple_path(g, b) for b in blocks]
            connected = blocks[0] == blocks[1]
            assert (c['marked_component'] is not None) == connected
            rec = dict(parent_index=index, core_key=old['core_key'], marks=marks,
                       omitted_black=omitted, marked_blocks=blocks, connected=connected, orientations=[])
            if any(path is None for path in paths):
                k = next(k for k, path in enumerate(paths) if path is None)
                rec.update(status='nonsimple_block_incompatible_with_deep_corridor',
                           obstruction=nonsimple_witness(g, blocks[k], marks[k]))
            else:
                maps, obstruction = embeddings(paths, marks, connected)
                if obstruction is not None:
                    rec.update(status=obstruction['reason'], obstruction=obstruction)
                else:
                    rec['status'] = 'covered_full_corridor_profiles'
                    for mapping in maps:
                        orientation_count += 1
                        assert len(set(mapping.values())) == len(mapping)
                        assert 0 <= min(mapping.values()) <= max(mapping.values()) < N
                        assert mapping[marks[0]] == LEFT and mapping[marks[1]] == RIGHT
                        assert all(physical(g, gauge, i, j) == Q.q(2 if i == j else -1 if abs(mapping[i]-mapping[j]) == 1 else 0)
                                   for i in mapping for j in mapping)
                        fixed, bad = [], []
                        for b, is_left in zip(omitted, (True, False)):
                            constraint, witness = constraints(g, gauge, mapping, b, is_left)
                            fixed.append(constraint)
                            if witness is not None:
                                bad.append(witness)
                        entry = dict(old_root_to_corridor_position=sorted(mapping.items()))
                        if bad:
                            entry.update(status='fixed_contact_exclusion', witnesses=bad)
                            profiles['no_contact_profile'] += 1
                        else:
                            key = tuple(fixed)
                            if key not in template_ids:
                                name = str(len(template_ids))
                                template_ids[key] = name
                                templates[name] = analyse_template(key)
                            name = template_ids[key]
                            uses[name] += 1
                            entry.update(status='complete_contact_template', template=name)
                            profiles.update(templates[name]['profile_outcomes'])
                        rec['orientations'].append(entry)
            counts[rec['status']] += 1
            rows.append(rec)
    assert counts == dict(nonsimple_block_incompatible_with_deep_corridor=16627,
                         wrong_fixed_mark_distance=2187, no_distance8_embedding=14,
                         covered_full_corridor_profiles=7984)
    assert profiles == dict(no_contact_profile=8337, strict_actual_edge_failure=57582, two_negative_directions=50137)
    assert len(rows) == 26812 and orientation_count == 13644 and len(templates) == 379
    sources = [Path(__file__), Path(Q.__file__), Path(T.__file__), Path(A.__file__), Path(H.__file__),
               ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    sources += [ROOT/f'{folder}/{T.S.TAG}.{ext}' for folder, ext in [('scripts', 'py'), ('results', 'json'), ('results', 'md')]]
    for folder, ext in [('scripts','py'),('results','json'),('results','md')]:
        sources.append(ROOT/f'{folder}/fresh_full_a43_corridor_high_domain_20261017.{ext}')
    conditional = rank65_diagnostic()
    return dict(status='PASS_FULL_A39_CORRIDOR_ALTERNATIVE_ON_COMPLETE_6703_HIGH_DOMAIN',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in sources},
        contexts=6703, cores=342, corners=26812, maximum_old_core_order=19, maximum_displayed_corner_order=17,
        corridor_order=N, marked_positions=[LEFT, RIGHT], outward_depth_assumption=DEPTH, retained_outward_whites=RETAINED_OUTWARD,
        strict_arm_order=DEPTH, unused_outer_whites=2,
        conditional_rank65_diagnostic=conditional,
        corner_outcomes=dict(counts), complete_path_orientations=orientation_count,
        weighted_profile_outcomes=dict(profiles), distinct_contact_templates=len(templates),
        distinct_joint_profile_checks=sum(len(t['joint_profiles']) for t in templates.values()),
        template_multiplicities=dict(uses), templates=templates, cases=rows,
        scope='For every corner of all6703 frozen high-overlap contexts, an actual whole finite A/B/D component with marked distance8 and both marks at feature distance>=17 forces an opposite-black path<=7 in the full actual supporting star. In the branch without that full-star path, the deep corridor is excluded. Both connected and disconnected displayed marked components, all gauges1/2, singular/non-Lorentz old principal cores and arbitrary spectator contacts are allowed. No paid-corner map, pair weight change, low/disjoint coverage, finite budget, dimension improvement or new polytope is asserted.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print(result['status'], result['corner_outcomes'], result['weighted_profile_outcomes'],
          'distinct_profiles', result['distinct_joint_profile_checks'], flush=True)
