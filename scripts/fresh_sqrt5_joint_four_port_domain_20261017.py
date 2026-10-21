#!/usr/bin/env python3
"""Exact four-port necessary domain and a new simplex/prism facet.

No optimizer. Geometric validity is the accompanying proof; the finite
enumeration is a necessary combinatorial domain, not a realizability claim.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_joint_four_port_domain_20261017'
QUADS = 'fresh_quad_product_aggregation_20261017'
PAIRS = tuple(combinations(range(4), 2))
COEFFICIENTS = (1, 2, 3, 0, 1, 1)


def dump(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'))


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def rank(rows):
    if not rows:
        return 0
    a = [[F(x) for x in row] for row in rows]
    r = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][col]
        a[r] = [x / scale for x in a[r]]
        for i in range(r + 1, len(a)):
            scale = a[i][col]
            if scale:
                a[i] = [x - scale * y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def affine_rank(points):
    return rank([[x - y for x, y in zip(p, points[0])] for p in points[1:]])


def partitions(prefix=()):
    if len(prefix) == 4:
        yield prefix
        return
    for value in range(max(prefix, default=-1) + 2):
        yield from partitions(prefix + (value,))


def close_square_groups(partition, edges):
    """Triangle plus two squares forces the two square edges into one Q."""
    parent = list(range(len(edges)))

    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i

    equalities = []
    for i, j in combinations(range(len(edges)), 2):
        common = set(edges[i]) & set(edges[j])
        if not common:
            continue
        hub = next(iter(common))
        a = next(x for x in edges[i] if x != hub)
        b = next(x for x in edges[j] if x != hub)
        if partition[a] == partition[b]:
            parent[find(j)] = find(i)
            equalities.append(dict(square_edges=[i, j], hub=hub, triangle_pair=sorted((a, b))))
    groups = defaultdict(list)
    for i, edge in enumerate(edges):
        groups[find(i)].append(edge)
    groups = sorted(sorted(group) for group in groups.values())
    checked = []
    for group in groups:
        hubs = set(group[0])
        for edge in group:
            hubs &= set(edge)
        valid_hubs = []
        for hub in sorted(hubs):
            leaves = sorted({x for edge in group for x in edge if x != hub})
            if len({partition[x] for x in leaves}) == 1:
                valid_hubs.append(dict(hub=hub, leaves=leaves))
        checked.append(dict(edges=group, valid_hubs=valid_hubs))
    return equalities, checked


def necessary_domain():
    records, points = [], set()
    counts = Counter()
    rejected_shapes = Counter()
    for part in partitions():
        blocks = [[i for i, p in enumerate(part) if p == j] for j in range(max(part) + 1)]
        triangles = {p for p in PAIRS if part[p[0]] == part[p[1]]}
        remaining = [p for p in PAIRS if p not in triangles]
        for mask in range(1 << len(remaining)):
            edges = [p for j, p in enumerate(remaining) if mask >> j & 1]
            equalities, groups = close_square_groups(part, edges)
            accepted = all(g['valid_hubs'] for g in groups)
            counts['patterns'] += 1
            counts['accepted' if accepted else 'rejected'] += 1
            rec = dict(index=len(records), triangle_partition=blocks,
                selected_square_pairs=edges, forced_same_Q=equalities,
                forced_Q_groups=groups, accepted=accepted)
            if accepted:
                ss = Counter(map(len, blocks))
                qq = Counter(len(g['edges']) for g in groups)
                point = (ss[2], ss[3], ss[4], qq[1], qq[2], qq[3])
                assert sum(a * b for a, b in zip(point, COEFFICIENTS)) <= 3
                rec['invariant_vector'] = point
                points.add(point)
            else:
                assert sorted(map(len, blocks)) == [2, 2] and len(edges) in (3, 4)
                rejected_shapes['three_edge_path' if len(edges) == 3 else 'four_cycle'] += 1
                rec['reason'] = 'FORCED_SAME_Q_CLASS_IS_NOT_A_STAR_WITH_LEAVES_IN_ONE_SIMPLEX_BLOCK'
            records.append(rec)
    points = sorted(points)
    tight = [p for p in points if sum(a * b for a, b in zip(p, COEFFICIENTS)) == 3]
    assert len(list(partitions())) == 15
    assert counts == dict(patterns=337, accepted=322, rejected=15)
    assert len(points) == 25 and affine_rank(points) == 6
    assert len(tight) == 6 and affine_rank(tight) == 5
    assert rejected_shapes == dict(three_edge_path=12, four_cycle=3)
    return dict(counts=dict(counts, triangle_partitions=15, invariant_points=len(points)),
        rejected_shapes=dict(rejected_shapes),
        coordinate_order=['s2', 's3', 's4', 'q1', 'q2', 'q3'],
        facet_coefficients=list(COEFFICIENTS), facet_capacity=3,
        affine_dimension=6, tight_affine_dimension=5,
        invariant_points=points, tight_points=tight, records=records,
        scope='All retained patterns are an enlarged necessary domain. They are not asserted geometrically realizable. For an actual selected dotted Q, its nonempty restriction is a star; if two of its edges occur, their leaves are in one triangle block, so the closure groups recover that whole restricted Q, including singleton groups.')


def template_scope():
    obj = json.loads((ROOT / f'results/{QUADS}.json').read_text())
    for p, h in obj['sources'].items():
        assert sha(p) == h, p
    records, corners, square_pairs = [], 0, 0
    for rec in obj['dotted_prism_census']['records']:
        dotted = set(rec['dotted_base'])
        assert len(dotted) == 2
        factors = list(map(set, rec['maximal_prism_factors']))
        assert dotted in factors
        leaves = next(s for s in factors if s != dotted)
        assert dotted.isdisjoint(leaves)
        footprints = defaultdict(set)
        for left, right in rec['quadrilateral_black_partitions']:
            assert len(left) == len(right) == 2
            for a in left:
                for b in right:
                    omitted = tuple(sorted((a, b)))
                    other_a = next(x for x in left if x != a)
                    other_b = next(x for x in right if x != b)
                    footprints[omitted].add(tuple(sorted((other_a, other_b))))
        assert set(footprints) == {tuple(sorted((a, b))) for a in dotted for b in leaves}
        assert len(footprints) == rec['prism_vertex_count'] == 2 * rec['prism_dimension']
        pd = {tuple(x['omitted']) for x in rec['positive_definite_two_deletions']}
        assert set(footprints) <= pd
        checks = []
        for omitted, edges in sorted(footprints.items()):
            hub = next(iter(dotted - set(omitted)))
            base_ports = sorted(leaves - set(omitted))
            expected = {tuple(sorted((hub, x))) for x in base_ports}
            assert edges == expected and len(edges) == rec['prism_dimension'] - 1
            assert len(edges) >= 1
            checks.append(dict(vertex_omitted=list(omitted), hub=hub,
                simplex_base_ports=base_ports, occupied_square_pairs=sorted(edges)))
            corners += 1
            square_pairs += len(edges)
        records.append(dict(index=rec['index'], source_index=rec['source_index'],
            prism_dimension=rec['prism_dimension'], corners=checks))
    assert len(records) == 1399
    return dict(records=records, counts=dict(dotted_records=len(records),
        actual_prism_vertex_occurrences=corners, occupied_square_pair_occurrences=square_pairs),
        scope='Literal full-star K2,2 corner union only, not every positive-definite two-deletion. Each leaf set spans the actual simplex base of its prism and hence belongs to one containing maximal simplex. This is a combinatorial recount of frozen templates, not a fresh arithmetic or geometric census.')


def fractional_witness():
    simplex = [dict(ports=[0, 1], mass='1'), dict(ports=[2, 3], mass='1')]
    quads = []
    for hub in range(4):
        leaves = [2, 3] if hub < 2 else [0, 1]
        quads.append(dict(hub=hub, leaves=leaves,
            occupied_pairs=[sorted((hub, x)) for x in leaves], mass='1/3'))
    pair_loads = []
    for pair in PAIRS:
        s = sum((F(c['mass']) for c in simplex if set(pair) <= set(c['ports'])), F())
        q = sum((F(c['mass']) for c in quads if list(pair) in c['occupied_pairs']), F())
        assert s + q <= 1
        pair_loads.append(dict(pair=pair, simplex_load=str(s), square_load=str(q)))
    triple_checks = []
    for u in combinations(range(4), 3):
        ss = sum((F(c['mass']) for c in simplex if len(set(c['ports']) & set(u)) >= 2), F())
        st = sum((F(c['mass']) for c in simplex if set(u) <= set(c['ports'])), F())
        q = sum((F(c['mass']) for c in quads if any(set(p) <= set(u) for p in c['occupied_pairs'])), F())
        wedge = sum((F(c['mass']) for c in quads if sum(set(p) <= set(u) for p in c['occupied_pairs']) >= 2), F())
        assert 2 * ss + st + q == 3
        assert wedge <= 1
        triple_checks.append(dict(ports=u, simplex_indicator=str(ss), all_three_simplex=str(st),
            quad_indicator=str(q), strengthened_joint_load=str(2 * ss + st + q),
            old_mixed_wedge_load=str(wedge)))
    r = sum((F(c['mass']) * (len(c['ports']) - 1) for c in simplex), F())
    q2 = sum((F(c['mass']) for c in quads), F())
    assert r == 2 and q2 == F(4, 3) and r + q2 > 3
    return dict(simplex_blocks=simplex, selected_prism_footprints=quads,
        pair_loads=pair_loads, triple_checks=triple_checks,
        full_simplex_packing_cover=dict(mass='1', blocks=[[0, 1], [2, 3]]),
        invariant_vector=['2', '0', '0', '0', '4/3', '0'],
        new_four_port_load=str(r + q2), capacity='3', violation='1/3',
        scope='Exact fractional separation from pair capacities, the full simplex-only disjoint-packing cone, all four strengthened three-port inequalities2s+t+q<=3, and old mixed-wedge capacities on these ports. No assertion of global flag, arithmetic-form, or geometric realization.')


def run():
    domain = necessary_domain()
    templates = template_scope()
    witness = fractional_witness()
    paths = [f'scripts/{TAG}.py', f'results/{QUADS}.json', f'results/{QUADS}.md',
        'results/fresh_prism_simplex_shared_capacity_20261017.md',
        'results/fresh_simplex_disjoint_ports_20261017.md']
    return dict(status='PASS_EXACT_SELECTED_DOTTED_PRISM_FOUR_PORT_FACET',
        sources={p:sha(p) for p in paths}, optimizer_invocations=0,
        theorem='For four ports U at a finite vertex of a compact simple Coxeter polytope, let r(U)=sum_S max(|C_S intersect U|-1,0) over maximal simplex port blocks, including singleton maximal edges. Let q2(U) count selected canonical dotted-prism aggregates occupying at least two square pairs inside U. Then r(U)+q2(U)<=3. Each selected Q footprint is a star whose leaves lie in one maximal simplex block; triangle plus two squares forces those squares into the same canonical aggregate.',
        necessary_domain=domain, frozen_template_scope_check=templates,
        fractional_relaxation_witness=witness,
        charge_recipe=dict(vertex_cost=3,
            simplex_corner_coefficient='max(|simplex_port_block intersect U|-1,0)',
            selected_quad_corner_coefficient='1 iff at least two occupied square pairs lie in U; otherwise0',
            averaging='Average the whole literal inequality over the finite forest automorphism orbit of U; retain one common actual-provenance mass.'),
        scope='This is a proved local inequality and exact necessary-domain facet. It has not been scanned on the saved rank13 dual or used in an optimizer. General non-dotted product/full-star aggregates, including Delta2 x Delta2, are outside the selected Q scope. No ambient dimension bound is claimed.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    result = json.loads(dump(run()))
    out = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(dump(result) + '\n')
    print(result['status'], result['necessary_domain']['counts'],
          result['frozen_template_scope_check']['counts'], flush=True)
