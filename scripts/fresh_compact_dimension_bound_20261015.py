#!/usr/bin/env python3
"""Exact dimension-29 bottleneck in Vinberg's planar-angle weighting.

No chamber nonexistence is claimed. The classification of spherical Coxeter
diagrams and Vinberg's global face inequality are explicit theorem inputs.
Only elementary graph, partition and h-vector arithmetic is implemented here.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import comb, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact_dimension_bound_20261015'
PRIMARY = 'https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt'


def connected(vertices, edges):
    vertices = set(vertices)
    if not vertices:
        return False
    reached = {min(vertices)}
    while True:
        enlarged = reached | {b for a, b in edges if a in reached and b in vertices}
        enlarged |= {a for a, b in edges if b in reached and a in vertices}
        if enlarged == reached:
            return reached == vertices
        reached = enlarged


def classical_end_segments(kind, n):
    if kind in ('A', 'B'):
        edges = [(i, i + 1) for i in range(n - 1)]
        ends = [frozenset(range(k)) for k in range(1, 7)]
        ends += [frozenset(range(n - k, n)) for k in range(1, 7)]
    else:
        assert kind == 'D'
        edges = [(i, i + 1) for i in range(n - 3)]
        edges += [(n - 3, n - 2), (n - 3, n - 1)]
        # An invariant end segment contains either the long end, or both
        # fork leaves. Single fork leaves are excluded by whisker invariance.
        ends = [frozenset(range(k)) for k in range(1, min(6, n - 2) + 1)]
        ends += [frozenset(range(n - k, n)) for k in range(3, 7)]
    assert all(connected(s, edges) and len(s) <= 6 for s in ends)
    # A separate subset definition checks the explicit list at the first
    # five ranks, including every possible overlap of the two end segments.
    if n <= 12:
        degree = Counter(x for e in edges for x in e)
        leaves = {x for x in range(n) if degree[x] == 1}
        brute = set()
        for k in range(1, 7):
            for row in combinations(range(n), k):
                s = frozenset(row)
                if not s & leaves or not connected(s, edges):
                    continue
                if kind == 'D' and ((n - 2 in s) != (n - 1 in s)):
                    continue
                brute.add(s)
        assert set(ends) == brute, (kind, n)
    weighted = {tuple(sorted(e)) for e in edges}
    for s in ends:
        weighted.update(combinations(sorted(s), 2))
    return dict(type=f'{kind}{n}', rank=n, weight=len(weighted),
                edges=edges, end_segments=[sorted(s) for s in ends],
                weighted_pairs=sorted(weighted))


def weight(n):
    if n <= 4:
        return n * (n - 1)
    if n <= 7:
        return comb(n, 2)
    if n == 8:
        return 24
    if n == 9:
        return 27
    return n + 19


def partitions(n, least=1):
    if n == 0:
        yield ()
    for k in range(least, n + 1):
        for tail in partitions(n - k, k):
            yield (k,) + tail


TYPES = {
    4: ('A4', 'B4', 'D4', 'F4', 'H4'),
    7: ('A7', 'B7', 'D7', 'E7'),
    8: ('A8', 'B8', 'D8', 'E8'),
    9: ('A9', 'B9', 'D9'),
    10: ('A10', 'B10', 'D10'),
}


def typed_forests(p, omit_h4=False):
    groups = []
    for n, multiplicity in sorted(Counter(p).items()):
        types = [t for t in TYPES[n] if not (omit_h4 and t == 'H4')]
        groups.append(list(combinations_with_replacement(types, multiplicity)))
    return [tuple(t for group in row for t in group) for row in product(*groups)]


def run():
    classical = [classical_end_segments(k, n)
                 for k in ('A', 'B', 'D') for n in range(8, 30)]
    assert len(classical) == 66
    assert all(r['weight'] == weight(r['rank']) for r in classical)
    # E8 is the only nonclassical component beyond rank seven. Its four
    # omitted pairs are a literal input from Vinberg, printed pp.80--81.
    assert comb(8, 2) - 4 == weight(8)
    deficits = [dict(rank=n, weight=weight(n), deficit=3*n-weight(n))
                for n in range(1, 30)]
    assert all(r['deficit'] >= 0 for r in deficits)
    assert [r['rank'] for r in deficits if r['deficit'] <= 2] == [4, 7, 8, 9, 10]
    all_partitions = list(partitions(29))
    assert len(all_partitions) == 4565
    good = [p for p in all_partitions if sum(weight(n) for n in p) > 84]
    assert len(good) == 9
    rows = []
    all_typed = []
    no_h4 = []
    for p in good:
        full = typed_forests(p)
        reduced = typed_forests(p, omit_h4=True)
        assert len(full) == prod(comb(len(TYPES[k]) + v - 1, v)
                                 for k, v in Counter(p).items())
        assert all(sum(int(t[1:]) for t in row) == 29 for row in full)
        assert len(set(full)) == len(full)
        all_typed.extend(full)
        no_h4.extend(reduced)
        rows.append(dict(ranks=p, weight=sum(weight(n) for n in p),
                         deficit=sum(3*n-weight(n) for n in p),
                         typed_count=len(full), without_H4_count=len(reduced),
                         typed_forests=full))
    assert len(all_typed) == len(set(all_typed)) == 2126
    assert len(no_h4) == len(set(no_h4)) == 1354
    assert set(no_h4) == {x for x in all_typed if 'H4' not in x}
    # Rank-two I2(m), including all unbounded finite m, has deficit four.
    # It therefore cannot occur in a rank-29 forest of weight >84.
    assert 3*2-weight(2) == 4
    # Exact odd-dimensional coefficient comparison, with h_i=h_{29-i}.
    coefficients = []
    for i in range(15):
        q = comb(29-i, 2) + comb(i, 2) - 2*98
        assert q == (14-i)*(15-i)
        assert 5*(comb(29-i, 2)+comb(i, 2)) - 2*comb(29, 2) - 2*84 == 5*q
        coefficients.append(dict(i=i, coefficient=q, five_times=5*q))
    assert Fraction(comb(29, 2), 1) / Fraction(4*29, 28) == 98
    assert 5*98-comb(29, 2) == 84
    assert coefficients[0]['five_times'] == 1050
    assert sum(r['five_times'] for r in coefficients) == 5600
    # 4cos²(pi/5)=(3+sqrt5)/2 is a two-cycle in the cyclic field.
    # Hence H4 forces sqrt5 into that field. An odd cubic field has no
    # quadratic subfield; Qsqrt2 and Qsqrt5 are distinct quadratic fields.
    return dict(
        status='PASS_EXACT_DIMENSION29_EXTREMAL_VERTEX_REDUCTION_NOT_EXCLUSION',
        sources={f'scripts/{TAG}.py': sha256((ROOT/f'scripts/{TAG}.py').read_bytes()).hexdigest()},
        theorem_inputs=[
            dict(source=PRIMARY, location='printed pp.70--72, Proposition1/2',
                 input='Simple polytope h-vector positivity/symmetry and weighted face inequality.'),
            dict(source=PRIMARY, location='printed pp.80--81, section6.1--2',
                 input='Weights: adjacency/end segments<=6/component<=7/E8 except4 pairs; double in component<=4.'),
            dict(source=PRIMARY, location='printed pp.81--91, remainder section6',
                 input='Every k-gonal two-face has weight sum >=5-k.'),
            dict(source=PRIMARY, location='printed p.91, section6.12',
                 input='E7 prevents uniform c<3 in the stated nonnegative diagram-local weighting architecture with zero weight across components. No separate rank29-only impossibility is inferred.'),
            dict(source=PRIMARY, location='printed p.69 Theorem2 and section7',
                 input='Compact arithmetic dimension>=22 fields Qsqrt2,Qsqrt5,Qcos(2pi/7).'),
            dict(input='Complete classification of irreducible finite Coxeter diagrams; B and C denote the same unoriented Coxeter type.')],
        dimension=29, classical_graph_checks=classical, classical_graph_count=len(classical),
        E8_excluded_pair_count_input=4, component_table=deficits,
        all_rank_partitions_checked=len(all_partitions), extremal_rank_partitions=rows,
        extremal_typed_count=len(all_typed), arithmetic_without_sqrt5_typed_count=len(no_h4),
        extremal_contains_E7=sum('E7' in x for x in all_typed),
        extremal_contains_E8=sum('E8' in x for x in all_typed),
        extremal_contains_H4=sum('H4' in x for x in all_typed),
        h_vector_slack_coefficients=coefficients,
        exact_balance='sum_vertices(sigma-84)=5*sum_i((14-i)*(15-i)*h_i)+sum_faces(sigma(F)-(5-k_F))',
        unconditional_balance_lower_bound=1050,
        necessary_extremal_vertex_count_lower_bound=350,
        all_label_scope='I2(m) for every finite m has rank2 deficit4, so no angle cutoff enters the extremal list.',
        open_target='Exclude all2126 corner types, or prove a global charge inequality contradicting the exact balance.',
        exclusions_proved=0,
        scope='Necessary local vertex types in any hypothetical compact H29 Coxeter polytope, with no arithmetic, facet-count, or nerve hypothesis. Not a classification of full polytopes.'
    )


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = json.loads(json.dumps(run(), sort_keys=True))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == result
    else:
        path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(result['status'])
    print('4565 rank partitions;9 extremal partitions;2126 typed forests;1354 without H4')
