#!/usr/bin/env python3
"""Exact five-port necessary domain and fifteen valid nonnegative facets.

The saved verifier uses only the standard library. Floating convex hull and
small local LP discovery proposed the inequalities and witness; neither is
trusted or called here. No global model is solved.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_joint_five_port_domain_20261017'
FOUR = 'fresh_sqrt5_joint_four_port_domain_20261017'
PAIRS = tuple(combinations(range(5), 2))
# s2,s3,s4,s5,q1,q2,q3,q4,capacity. No completeness claim for this list.
FACETS = (
    (0,1,1,1,0,0,0,0,1),
    (1,0,1,2,0,0,1,1,2),
    (1,1,2,2,0,0,0,0,2),
    (0,1,2,3,0,1,1,1,3),
    (1,3,4,5,0,1,1,1,5),
    (1,7,9,10,1,1,1,1,10),
    (3,5,9,10,1,1,1,1,10),
    (3,5,8,10,1,1,2,2,10),
    (3,6,10,12,1,2,2,2,12),
    (4,7,11,14,1,2,3,3,14),
    (3,8,12,15,1,3,3,3,15),
    (5,10,17,20,2,3,3,3,20),
    (5,10,16,20,2,3,4,4,20),
    (6,12,19,24,2,4,5,5,24),
    (7,22,33,40,4,7,7,7,40),
)


def dump(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'))


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def partitions(prefix=()):
    if len(prefix) == 5:
        yield prefix
        return
    for value in range(max(prefix, default=-1) + 2):
        yield from partitions(prefix + (value,))


def affine_rank(points):
    rows = [[F(a-b) for a,b in zip(p,points[0])] for p in points[1:]]
    rank = 0
    for col in range(len(points[0])):
        pivot = next((i for i in range(rank,len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank],rows[pivot] = rows[pivot],rows[rank]
        scale = rows[rank][col]
        rows[rank] = [a/scale for a in rows[rank]]
        for i in range(rank+1,len(rows)):
            scale = rows[i][col]
            rows[i] = [a-scale*b for a,b in zip(rows[i],rows[rank])]
        rank += 1
    return rank


def square_components(part, edges):
    """Independent block-pair connected components, without union-find."""
    by_pair = defaultdict(list)
    for a,b in edges:
        assert part[a] != part[b]
        by_pair[tuple(sorted((part[a],part[b])))].append((a,b))
    groups = []
    for pair, es in sorted(by_pair.items()):
        unseen = set(es)
        while unseen:
            first = min(unseen)
            unseen.remove(first)
            component = [first]
            vertices = set(first)
            while True:
                added = sorted(e for e in unseen if vertices.intersection(e))
                if not added:
                    break
                for e in added:
                    unseen.remove(e)
                    component.append(e)
                    vertices.update(e)
            hubs = set(component[0])
            for e in component[1:]:
                hubs.intersection_update(e)
            groups.append(dict(block_pair=pair, edges=sorted(component), hubs=sorted(hubs)))
    return sorted(groups, key=lambda g:g['edges'])


@lru_cache(None)
def recursive_pair_profiles(a,b):
    """Separate unlabelled star-forest count recursion for one block pair."""
    if not a or not b:
        return ((0,0,0,0),)
    # A distinguished A vertex is isolated, a center, or a leaf of a B center.
    out = set(recursive_pair_profiles(a-1,b))
    for degree in range(1,b+1):
        for p in recursive_pair_profiles(a-1,b-degree):
            q = list(p); q[degree-1] += 1; out.add(tuple(q))
    for degree in range(2,a+1):
        for p in recursive_pair_profiles(a-degree,b-1):
            q = list(p); q[degree-1] += 1; out.add(tuple(q))
    return tuple(sorted(out))


def necessary_domain():
    points = Counter()
    records = []
    accepted = rejected = 0
    for part in partitions():
        sizes = Counter(Counter(part).values())
        simplex = tuple(sizes[k] for k in range(2,6))
        cross = [p for p in PAIRS if part[p[0]] != part[p[1]]]
        for mask in range(1 << len(cross)):
            edges = [e for j,e in enumerate(cross) if (mask >> j)&1]
            groups = square_components(part, edges)
            valid = all(g['hubs'] for g in groups)
            record = dict(partition=part, square_mask=mask, accepted=valid)
            if valid:
                counts = Counter(len(g['edges']) for g in groups)
                vector = simplex + tuple(counts[k] for k in range(1,5))
                points[vector] += 1
                record['vector'] = vector
                accepted += 1
            else:
                record['nonstar_components'] = [g for g in groups if not g['hubs']]
                rejected += 1
            records.append(record)
    # Reconstruct the invariant domain independently via component recursion.
    second = set()
    shapes = set(tuple(sorted(Counter(p).values())) for p in partitions())
    for shape in sorted(shapes):
        qs = {(0,0,0,0)}
        for a,b in combinations(shape,2):
            qs = {tuple(x+y for x,y in zip(q,r)) for q in qs
                  for r in recursive_pair_profiles(a,b)}
        simplex = tuple(shape.count(k) for k in range(2,6))
        second.update(simplex + q for q in qs)
    assert set(points) == second and len(points) == 74
    ordered = sorted(points)
    assert affine_rank(ordered) == 8
    facets = []
    for j,co in enumerate(FACETS):
        assert min(co) >= 0
        slacks = [co[-1]-sum(a*b for a,b in zip(co,p)) for p in ordered]
        assert min(slacks) == 0
        tight = [p for p,s in zip(ordered,slacks) if s == 0]
        assert affine_rank(tight) == 7
        facets.append(dict(index=j,coefficients=co[:-1],capacity=co[-1],
            tight_points=tight,tight_affine_dimension=7))
    return dict(counts=dict(labelled_partitions=52,partition_shapes=len(shapes),
        square_patterns=len(records),accepted=accepted,rejected=rejected,
        invariant_points=len(points)), coordinate_order=['s2','s3','s4','s5','q1','q2','q3','q4'],
        affine_dimension=8,invariant_points=ordered,
        point_multiplicities=[points[p] for p in ordered],facets=facets,
        records=records,independent_recursive_domain_match=True)


def separating_witness():
    # Five pair-disjoint 2-stars partition K5; no positive simplex mass.
    footprints = [((0,1),(1,2)),((0,2),(2,4)),((0,3),(3,4)),
                  ((1,3),(2,3)),((0,4),(1,4))]
    assert sorted(e for es in footprints for e in es) == list(PAIRS)
    checks = []
    for k,cap in ((2,1),(3,3),(4,3)):
        for us in combinations(range(5),k):
            u = set(us)
            counts = [sum(set(e)<=u for e in es) for es in footprints]
            load = sum(counts) if k == 2 else sum(x >= (2 if k == 4 else 1) for x in counts)
            assert load <= cap
            if k == 3:
                assert sum(x>=2 for x in counts) <= 1
            checks.append(dict(ports=us,family={2:'pair',3:'strong_three',4:'four'}[k],
                load=load,capacity=cap))
    assert sum(len(es)>=2 for es in footprints) == 5 > FACETS[3][-1]
    return dict(simplex_blocks=[],simplex_packing_cover=[],
        quad_stars=[dict(edges=es,mass='1') for es in footprints],old_constraints=checks,
        five_port_facet=3,new_load=5,new_capacity=3,
        scope='Exact local relaxation witness only. It passes all pair, strengthened three-port, four-port and old three-port wedge capacities, and the complete simplex-only packing cone. It is not geometrically realizable and is not claimed to satisfy coface flags or arithmetic-form constraints.')


def run():
    parent = json.loads((ROOT / f'results/{FOUR}.json').read_text())
    for p,h in parent['sources'].items():
        assert sha(p) == h,p
    domain = necessary_domain()
    witness = separating_witness()
    paths = [f'scripts/{TAG}.py',f'scripts/{FOUR}.py',f'results/{FOUR}.json',
        f'results/{FOUR}.md','results/fresh_quad_product_aggregation_20261017.md',
        'results/fresh_prism_simplex_shared_capacity_20261017.md',
        'results/fresh_simplex_disjoint_ports_20261017.md']
    return dict(status='PASS_EXACT_FIVE_PORT_SELECTED_PRISM_FACETS',
        sources={p:sha(p) for p in paths},necessary_domain=domain,
        separating_relaxation_witness=witness,
        verification_optimizer_invocations=0,global_optimizer_invocations=0,
        discovery='Floating local ConvexHull proposals and fifteen small local LP witness searches; every retained facet and the selected witness are rechecked exactly here. No completeness claim for the list of facets.',
        scope='Compact simple Coxeter vertices; maximal simplex port blocks and selected canonical dotted full-star prism aggregates only. Their square graphs between any two triangle blocks are disjoint unions of stars. All retained patterns form a necessary domain, not a realizability classification. No saved global dual scan or dimension theorem in this package.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args = parser.parse_args()
    result = json.loads(dump(run()))
    output = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(output.read_text())
    else:
        output.write_text(dump(result)+'\n')
    print(result['status'],result['necessary_domain']['counts'],flush=True)
