#!/usr/bin/env python3
"""Independent complete T8 traversal, exact forest Schur arithmetic.

No project/producer imports. Incidence and reviewed analytic theorems are
inherited JSON data. Assignment tuples are immutable. Q(sqrt5) Schur
elimination replaces the producer's Dynkin-diagram lookup.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_t8_sparse_search_20261014_independent_dfs'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'
BASE = 'results/fresh_h9_t8_sparse_search_20261014_all_blocks.json'
BASE_REVIEW = 'results/fresh_h9_t8_all_blocks_20261014_review.json'
HIGH = 'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks.json'
HIGH_REVIEW = 'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.json'
MAPS = 'results/fresh_h9_t8_sparse_search_20261014_bridge_variants.json'
MAP_REVIEW = 'results/fresh_h9_t8_two_high_transfers_20261014_review.json'
PRIMARY = 'results/fresh_h9_t8_sparse_search_20261014_outside_case_08.json'

def digest(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()

ZERO = (F(0), F(0))
ONE = (F(1), F(0))
SQUARES = {3: (F(1, 4), F(0)), 4: (F(1, 2), F(0)),
           5: (F(3, 8), F(1, 8))}

def sign(x):
    a, b = x
    if b == 0:
        return (a > 0) - (a < 0)
    if a == 0:
        return (b > 0) - (b < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    q = a*a - 5*b*b
    assert q != 0
    return ((q > 0) - (q < 0)) * (1 if a > 0 else -1)

def subtract(a, b):
    return (a[0]-b[0], a[1]-b[1])

def divide(a, b):
    den = b[0]*b[0]-5*b[1]*b[1]
    assert den != 0
    return ((a[0]*b[0]-5*a[1]*b[1])/den,
            (a[1]*b[0]-a[0]*b[1])/den)

@lru_cache(None)
def schur(tree):
    """Root pivot after exact leaf elimination, or None if not PD."""
    pivot = ONE
    for label, child in tree:
        value = schur(child)
        if value is None:
            return None
        pivot = subtract(pivot, divide(SQUARES[label], value))
    return pivot if sign(pivot) > 0 else None

def run():
    manifest = {}
    def read(p):
        obj = json.loads((ROOT / p).read_text())
        if p in manifest:
            return obj
        manifest[p] = digest(p)
        for q, h in obj.get('sources', {}).items():
            assert digest(q) == h, q
            if q.endswith('.json'):
                read(q)
        return obj
    geom = read(GEOM)['representatives'][8]['geometry']
    base, high, maps = [read(p) for p in (BASE, HIGH, MAPS)]
    for p in (BASE_REVIEW, HIGH_REVIEW, MAP_REVIEW):
        assert 'PASS' in read(p)['status']
    primary = read(PRIMARY)
    missing = [frozenset(m) for m in geom['minimal_nonfaces']]
    assert Counter(map(len, missing)) == {2: 11, 3: 10}
    vertices = [v for v in combinations(range(18), 9)
                if not any(m <= frozenset(v) for m in missing)]
    assert set(vertices) == {tuple(v['incident']) for v in geom['vertices']}
    assert len(vertices) == 462
    pairs = list(combinations(range(18), 2))
    number = {p: i for i, p in enumerate(pairs)}
    vertex_edges = [tuple(number[e] for e in combinations(v, 2)) for v in vertices]
    touching = [tuple(j for j, es in enumerate(vertex_edges) if i in es)
                for i in range(len(pairs))]
    triangles = [tuple(sorted(m)) for m in missing if len(m) == 3]
    triangle_edges = [tuple(number[e] for e in combinations(t, 2)) for t in triangles]
    categories = [t for t in product(range(2, 8), repeat=3)
                  if sum((F(1, m) for m in t), F(0)) < 1]
    assert len(categories) == 175
    categories.sort(key=lambda t: (max(t), sum(t), t.count(2), t))
    twins = [(1, 2), (4, 5), (7, 8), (12, 13)]
    all_twins = set(x for p in twins for x in p)
    for choices in product((False, True), repeat=4):
        perm = list(range(18))
        for (a, b), take in zip(twins, choices):
            if take:
                perm[a], perm[b] = b, a
        assert {tuple(sorted(perm[x] for x in v)) for v in vertices} == set(vertices)
        assert {frozenset(perm[x] for x in m) for m in missing} == set(missing)
    domains = [tuple(categories) for t in triangles]
    reductions = []
    for a, b in twins:
        j = next(j for j, t in enumerate(triangles) if a in t and b in t)
        leaf = next(x for x in triangles[j] if x not in (a, b))
        assert leaf not in all_twins
        edges = list(combinations(triangles[j], 2))
        u = edges.index(tuple(sorted((leaf, a))))
        v = edges.index(tuple(sorted((leaf, b))))
        domains[j] = tuple(t for t in domains[j] if t[u] <= t[v])
        reductions.append(dict(triangle=list(triangles[j]), twin=[a, b], leaf=leaf))

    @lru_cache(maxsize=None)
    def elliptic(edges):
        adjacency = {}
        for index, label in edges:
            a, b = pairs[index]
            assert 3 <= label <= 7
            adjacency.setdefault(a, []).append((b, label))
            adjacency.setdefault(b, []).append((a, label))
        unseen = set(adjacency)
        while unseen:
            root = min(unseen)
            component = {root}
            queue = [root]
            unseen.remove(root)
            for v in queue:
                for w, label in adjacency[v]:
                    if w in unseen:
                        unseen.remove(w)
                        component.add(w)
                        queue.append(w)
            if len(component) <= 2:
                continue
            if sum(len(adjacency[v]) for v in component) != 2*(len(component)-1):
                # A cycle has the all-ones quadratic witness <=0.
                return False
            if any(label >= 6 for v in component for w, label in adjacency[v]):
                # A connected 3-node extension of a >=6 edge has det<=0.
                return False
            def shape(v, parent):
                return tuple(sorted((label, shape(w, v))
                                    for w, label in adjacency[v] if w != parent))
            if schur(shape(root, None)) is None:
                return False
        return True

    clauses = []
    for a, b in combinations(missing, 2):
        if a & b or any(m <= a | b and m != a and m != b for m in missing):
            continue
        crossing = tuple(sorted(number[tuple(sorted((x, y)))] for x in a for y in b))
        assert all(touching[k] for k in crossing)
        clauses.append(crossing)
    assert len(clauses) == 60
    families = [(base['fixed_labels'], base['finite_label_lower_bounds'])]
    swap = list(range(18))
    swap[1], swap[2] = 2, 1
    families.append(([[*sorted((swap[a], swap[b])), m] for a, b, m in base['fixed_labels']],
                     [[*sorted((swap[a], swap[b])), m] for a, b, m in base['finite_label_lower_bounds']]))
    for bounds in maps['distinct_condition_sets']:
        assert sorted(m for a, b, m in bounds) == [4, 4, 7, 7]
        families.append(([], bounds))
    predicates = [[(number[a, b], True, m) for a, b, m in fixed] +
                  [(number[a, b], False, m) for a, b, m in bounds]
                  for fixed, bounds in families]
    pair_occurrences = {k: [(j, ts.index(k)) for j, ts in enumerate(triangle_edges) if k in ts]
                        for condition in predicates for k, eq, m in condition}
    assert all(pair_occurrences.values())
    counts = Counter()
    nodes = 0
    survivor = None
    leaves = 0

    def theorem_hit(assignment, active_domains):
        for family, condition in enumerate(predicates):
            for k, eq, threshold in condition:
                if assignment[k]:
                    possible = {assignment[k]}
                else:
                    projections = [{r[pos] for r in active_domains[j]}
                                   for j, pos in pair_occurrences[k]]
                    possible = set.intersection(*projections)
                if not possible or any((m != threshold if eq else m < threshold) for m in possible):
                    break
            else:
                return family
        return None

    def visit(assignment, active_domains, remaining):
        nonlocal nodes, survivor, leaves
        nodes += 1
        if nodes % 10000 == 0:
            print('independent nodes', nodes, dict(counts),
                  'forest masks', elliptic.cache_info().currsize, flush=True)
        hit = theorem_hit(assignment, active_domains)
        if hit is not None:
            counts[f'theorem_{hit}'] += 1
            return False
        finite_masks = [tuple((k, assignment[k]) for k in es if assignment[k] > 2)
                        for es in vertex_edges]
        possible_nonright = {}
        def can_be_nonright(k):
            if assignment[k]:
                return assignment[k] > 2
            if k not in possible_nonright:
                possible_nonright[k] = all(elliptic(tuple(sorted(finite_masks[j] + ((k, 3),))))
                                          for j in touching[k])
            return possible_nonright[k]
        if any(not any(can_be_nonright(k) for k in clause) for clause in clauses):
            counts['two_negative_blocks'] += 1
            return False
        if not remaining:
            leaves += 1
            survivor = [list(p) + [m] for p, m in zip(pairs, assignment) if m]
            return True
        chosen = min(remaining, key=lambda j: (
            len(active_domains[j]),
            -sum(assignment[k] != 0 for k in triangle_edges[j]),
            -sum(len(set(triangle_edges[j]) & set(triangle_edges[k])) for k in remaining if k != j), j))
        later = remaining - {chosen}
        for row in active_domains[chosen]:
            assert all(assignment[k] in (0, label) for k, label in zip(triangle_edges[chosen], row))
            child = list(assignment)
            changed = []
            for k, label in zip(triangle_edges[chosen], row):
                if not child[k]:
                    child[k] = label
                    if label > 2:
                        changed.append(k)
            child = tuple(child)
            affected = {j for k in changed for j in touching[k]}
            if not all(elliptic(tuple((k, child[k]) for k in vertex_edges[j] if child[k] > 2))
                       for j in affected):
                counts['nonpositive_vertex'] += 1
                continue
            propagated = list(active_domains)
            for j in later:
                propagated[j] = tuple(r for r in active_domains[j]
                                      if all(child[k] in (0, label) for k, label in zip(triangle_edges[j], r)))
                if not propagated[j]:
                    counts['inconsistent_triangle'] += 1
                    break
            else:
                if visit(child, tuple(propagated), later):
                    return True
        return False

    visit((0,)*len(pairs), tuple(domains), frozenset(range(len(triangles))))
    assert survivor is None and leaves == 0
    assert primary['status'] == 'EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY'
    # This comparison is a result check, not a traversal bound or stop rule.
    assert nodes == primary['nodes']
    expected = dict(nonpositive_vertex=primary['prunes']['spherical_vertex'],
                    two_negative_blocks=primary['prunes']['orthogonal_negative_blocks'],
                    theorem_0=primary['prunes']['certified_all_blocks_sparse_partial'],
                    theorem_1=primary['prunes']['certified_all_blocks_sparse_twin_partial'],
                    theorem_2=primary['prunes']['certified_two_high_placement_0_partial'],
                    theorem_3=primary['prunes']['certified_two_high_placement_1_partial'])
    assert dict(counts) == expected, (counts, expected)
    sources = [GEOM, BASE, BASE_REVIEW, HIGH, HIGH_REVIEW, MAPS, MAP_REVIEW, PRIMARY,
               f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_COMPLETE_T8_DFS_EXCLUSION',
        sources={p: digest(p) for p in sources},
        recursively_checked_manifests=len(manifest), vertex_count=len(vertices),
        categories=len(categories), triangle_count=len(triangles), bridge_clauses=len(clauses),
        simultaneous_symmetries=16, reductions=reductions, nodes=nodes, prunes=dict(counts),
        finite_masks_checked=elliptic.cache_info().currsize,
        rooted_tree_pivots_checked=schur.cache_info().currsize,
        complete=True, cutoff=None, survivor=None, unexplored_frontier=None,
        arithmetic='Separate exact Fraction Q(sqrt5) forest Schur elimination; cycles and connected >=6 edges have direct nonpositive witnesses; isolated I2 accepts every finite label.',
        scope='Independent recursion, immutable assignments, propagation, triangle categories and finite-diagram arithmetic; no producer/project imports. Literal incidence realization and reviewed analytic theorem conditions are inherited. Exhaustive empty mandatory-label search excludes compact H9 realizations of this one T8 incidence, not other incidences or noncompact chambers.')

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    d = run()
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert d == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(d, indent=2) + '\n')
    print(d['status'], d['nodes'], d['prunes'], flush=True)
