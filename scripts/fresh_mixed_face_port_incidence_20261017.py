#!/usr/bin/env python3
"""Exact mixed triangle/quadrilateral port features and a local packing cut.

This introduces no optimizer and makes no completeness claim for the inherited
face catalogue. It reconstructs literal three-port indicators for every saved
simplex and marked dotted quadrilateral aggregate.
"""
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_mixed_face_port_incidence_20261017'
TRI = 'fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017'
QUAD = 'fresh_quad_product_aggregation_20261017'


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def read(tag):
    doc = json.loads((ROOT / f'results/{tag}.json').read_text())
    for path, expected in doc.get('sources', {}).items():
        assert sha(path) == expected, (tag, path)
    return doc


def pair(a, b):
    assert a != b
    return tuple(sorted((a, b)))


def triple_indicators(edges):
    edges = set(map(tuple, edges))
    ports = sorted({v for e in edges for v in e})
    return [list(U) for U in combinations(ports, 3)
            if sum(e in edges for e in combinations(U, 2)) >= 2]


def small_packing():
    pairs = list(combinations(range(3), 2))
    masks = list(range(1, 8))
    chosen = [m for m in masks if m.bit_count() >= 2]
    # Exhaust all families of distinct nonempty footprints on a three-port set.
    # Repeated selected footprints also share a pair and are impossible.
    feasible = []
    for flags in range(1 << len(masks)):
        family = [m for i, m in enumerate(masks) if flags >> i & 1]
        if any(a & b for a, b in combinations(family, 2)):
            continue
        value = sum(m in chosen for m in family)
        assert value <= 1
        feasible.append(dict(footprints=family,indicator_sum=value))
    wedges = [m for m in masks if m.bit_count() == 2]
    mass = F(1, 2)
    loads = [sum((mass for m in wedges if m >> i & 1), F()) for i in range(3)]
    assert loads == [F(1)] * 3 and len(wedges) * mass == F(3, 2)
    return dict(port_pairs=[list(p) for p in pairs],selected_footprints=chosen,
                exhaustive_disjoint_families=feasible,
                integer_maximum=1,
                fractional_witness=[dict(footprint=m,mass=str(mass)) for m in wedges],
                pair_loads=list(map(str,loads)),fractional_indicator_sum='3/2',
                scope='A combinatorial relaxation witness, not an asserted arithmetic realization of all three wedges at one forest.')


def quad_features(record):
    at_vertex = defaultdict(set)
    quadrilaterals = record['quadrilateral_black_partitions']
    for pos, neg in quadrilaterals:
        assert len(pos) == len(neg) == 2
        for i, j in product(pos, neg):
            a = next(x for x in pos if x != i)
            b = next(x for x in neg if x != j)
            v = pair(i, j)
            e = pair(a, b)
            assert e not in at_vertex[v], ('duplicate actual quad corner', v, e)
            at_vertex[v].add(e)
    assert sum(map(len,at_vertex.values())) == 4 * len(quadrilaterals)
    corners = []
    for v, edges in sorted(at_vertex.items()):
        assert all(not(set(v) & set(e)) for e in edges)
        selected = triple_indicators(edges)
        # Every such graph is bipartite by the exact kernel-sign theorem.
        # Hence each selected triple has exactly two occupied port pairs.
        assert all(sum(e in edges for e in combinations(U,2)) == 2 for U in selected)
        corners.append(dict(vertex_omitted=list(v),occupied_pairs=[list(e) for e in sorted(edges)],
                            selected_triples=selected))
    return corners


def run():
    tri = read(TRI)
    quad = read(QUAD)
    counts = Counter()
    trace = hashlib.sha256()
    simplex_dimensions = Counter()
    for star in tri['stars']:
        S = {r['omitted'] for r in star['spherical_deletions']}
        d = star['simplex_dimension']
        assert len(S) == d + 1
        hits = 0
        for omitted in sorted(S):
            ports = sorted(S - {omitted})
            edges = list(combinations(ports,2))
            triples = triple_indicators(edges)
            assert triples == [list(U) for U in combinations(ports,3)]
            hits += len(triples)
            trace.update((json.dumps([star['index'],omitted,triples],separators=(',',':'))+'\n').encode())
        assert hits == (d + 1) * comb(d,3)
        counts['simplex_stars'] += 1
        counts['simplex_corners'] += d + 1
        counts['simplex_triple_indicators'] += hits
        simplex_dimensions[d] += 1
    records = []
    prism_dimensions = Counter()
    for source in quad['dotted_prism_census']['records']:
        corners = quad_features(source)
        k = source['prism_dimension']
        assert len(corners) == 2 * k
        assert len(source['quadrilateral_black_partitions']) == comb(k,2)
        expected_edges = k - 1
        for corner in corners:
            assert len(corner['occupied_pairs']) == expected_edges
            assert len(corner['selected_triples']) == comb(k - 1,2)
        total = sum(len(c['selected_triples']) for c in corners)
        assert total == 2 * k * comb(k - 1,2)
        records.append(dict(source_record=source['index'],prism_dimension=k,
                            aggregate_demand=source['quadrilateral_count'],
                            corners=corners,total_triple_indicators=total))
        counts['marked_dotted_stars'] += 1
        counts['dotted_aggregate_corners'] += len(corners)
        counts['dotted_quad_corner_occurrences'] += sum(len(c['occupied_pairs']) for c in corners)
        counts['dotted_triple_indicators'] += total
        prism_dimensions[k] += 1
    overlap = quad['exact_finite_overlap']
    overlap_record = dict(quadrilateral_black_partitions=overlap['distinct_quadrilaterals'])
    overlap_corners = quad_features(overlap_record)
    assert sum(len(c['occupied_pairs']) for c in overlap_corners) == 20
    paths = [f'scripts/{TAG}.py',f'results/{TRI}.json',f'results/{QUAD}.json',
             'results/fresh_sqrt5_rank13_triangle_simplex_flags_20261017.md',
             'results/fresh_simplex_port_incidence_20261017.md',f'results/{QUAD}.md']
    return dict(status='PASS_EXACT_MIXED_THREE_PORT_INCIDENCE_FEATURES',
        sources={p:sha(p) for p in paths},counts=dict(counts),
        simplex_dimension_counts=dict(sorted(simplex_dimensions.items())),
        dotted_prism_dimension_counts=dict(sorted(prism_dimensions.items())),
        simplex_feature_trace_sha256=trace.hexdigest(),
        exact_three_port_packing=small_packing(),
        dotted_corner_features=records,
        finite_overlap_corner_features=overlap_corners,
        incidence_theorem='Partition triangular two-faces by their unique maximal simplex and quadrilateral two-faces by their entire incident-wall star. At a simple vertex, each aggregate occupies the two-port pairs of its incident two-faces. Different aggregates have disjoint pair sets. For any triple U of ports, at most one aggregate occupies at least two of its three pairs. This inequality can be averaged under the labelled finite forest automorphism group and used as one shared vertex-budget charge for simplex and quadrilateral aggregate rows.',
        scope='Exact incidence theorem and literal features on saved stars only. No optimizer, complete face coverage, arithmetic coexistence claim for the fractional witness, ambient compact polytope or dimension exclusion. Marked dotted records are not deduplicated across source isomorphisms. Single quadrilateral corners and simple triangular-face corners have zero three-port charge.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    args = ap.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],flush=True)
