#!/usr/bin/env python3
"""Compose complete Qsqrt5 face coverage with exact rank25/26/27 rules.

No queue or optimizer arithmetic is imported. The separately reviewed
coverage maps and literal graph implementation are reused read-only.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_dimension27_coverage_20261015 as C
import fresh_sqrt5_complete_face_weights_20261015_review as W

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_dimension24_coverage_20261015'
DIMENSIONS = (25, 26, 27)


class Binder(C.Binder):
    def __init__(self):
        super().__init__()
        self.paths.update([f'scripts/{TAG}.py', f'scripts/{W.TAG}.py',
                           f'scripts/{W.V.TAG}.py'])
        self.rules = {}
        self.weight_reviews = {}
        self.dimension_skips = Counter()
        for n in DIMENSIONS:
            primary = self.load(W.PRIMARY + '_n' + str(n))
            review = self.review(W.TAG + '_n' + str(n))
            assert primary['status'] == 'PASS_EXACT_FACE_WEIGHTS_WITH_COMPLETE_FINITE_QUEUE'
            assert review['status'] == 'PASS_INDEPENDENT_EXACT_WEIGHTS_ON_PRIOR_CATALOGUE_AND_ALL_FINITE_SINGULAR_STOPS'
            assert primary['dimension'] == review['dimension'] == n
            assert review['finite_queue_states'] == 5172 and review['all_finite_singular_stops'] == 2070
            assert review['sources'][f'results/{W.PRIMARY}_n{n}.json'] == sha256(
                (ROOT / f'results/{W.PRIMARY}_n{n}.json').read_bytes()).hexdigest()
            self.rules[n] = W.V.Rule(primary)
            self.weight_reviews[n] = review

    def payment(self, row, origin):
        cases = []
        for n, rule in self.rules.items():
            if row['nodes'] - 2 > n:
                self.dimension_skips[n] += 1
                continue
            for rec in rule.payment(row):
                cases.append(dict(dimension=n, external_A1=rec['external_A1'],
                                  corners=rec['corner_weights'], total=rec['payment']))
        # An over-rank row is excluded only at the indicated target dimension.
        # It remains in the arithmetic queue, in the preceding28/29 binder,
        # and in this explicit list with an empty case list if all three skip.
        self.payments.append(dict(origin=origin, nodes=row['nodes'], cases=cases))


def triangle_and_vertex_bounds(B):
    triangle = B.review('fresh_sqrt5_triangle_weights_20261015_review')
    assert len(triangle['markings']) == 76
    triangle_totals = []
    for rec in triangle['markings']:
        edges = [(i, j, m) for (i, j), m in
                 zip(combinations(range(rec['order']), 2), rec['labels'])]
        total = F(0)
        for omit in rec['black']:
            marked = [i for i in rec['black'] if i != omit]
            parts = W.R.components(W.R.graph([i for i in range(rec['order']) if i != omit], edges))
            g = next(h for h in parts if marked[0] in h)
            total += W.R.baseline(g, *marked) if marked[1] in g else 0
        assert total >= 2
        triangle_totals.append(str(total))
    records = []
    for n, rule in B.rules.items():
        models, totals, fronts, _, _ = W.V.component_bounds(rule)
        review = B.weight_reviews[n]
        assert fronts == review['forest_maxima']
        maximum = max(F(r['maximum']) for r in fronts)
        assert maximum == F(review['vertex_maximum'])
        paired = [F((n-i)*(n-i-1) + i*(i-1), 4) for i in range(n + 1)]
        minimum = min(paired)
        assert paired[0] > minimum
        threshold = 5 * minimum - F(n*(n-1), 2)
        assert threshold == {25:F(60), 26:F(65), 27:F(143,2)}[n]
        assert maximum <= threshold == F(review['strict_face_threshold'])
        records.append(dict(dimension=n, component_types=len(models), exact_A1_contexts=len(fronts),
                            maximum_vertex_weight=str(maximum), forest_maxima=fronts,
                            minimum_paired_f2_coefficient=str(minimum),
                            strict_face_weight_threshold=str(threshold),
                            margin=str(threshold-maximum),
                            strictness='h_0=h_n=1 and the extreme paired coefficient is larger'))
    return dict(triangle_markings=76, baseline_totals=triangle_totals, bounds=records)


def run(stage=False):
    if not stage:
        prior = json.loads((ROOT / f'results/{C.TAG}.json').read_text())
        assert prior['status'] == 'PASS_COMPLETE_QSQRT5_FACE_COVERAGE_EXCLUDES_COMPACT_DIMENSIONS28_AND29'
        assert prior['complete_coverage'] and prior['arithmetic_ground_field'] == 'Q(sqrt5)'
        assert prior['compact_dimension_upper_bound'] == 27 and prior['excluded_dimensions'] == [28,29]
    B = Binder()
    if not stage:
        B.load(C.TAG)
    bounds = triangle_and_vertex_bounds(B)
    old, roots, dotted_ids = C.old_branches(B)
    dotted = C.dotted_branch(B, dotted_ids)
    finite = None if stage else C.finite_branch(B, roots)
    sources = B.bind_manifests()
    case_counts = Counter(r['dimension'] for row in B.payments for r in row['cases'])
    row_counts = Counter(n for row in B.payments for n in DIMENSIONS
                         if any(r['dimension'] == n for r in row['cases']))
    return dict(status='STAGED_PAYMENTS_ONLY_NO_DIMENSION_CLAIM' if stage else
                'PASS_COMPLETE_QSQRT5_FACE_COVERAGE_EXCLUDES_COMPACT_DIMENSIONS25_THROUGH29',
                arithmetic_ground_field='Q(sqrt5)', compact_dimension_upper_bound=None if stage else 24,
                excluded_dimensions=[] if stage else [25,26,27,28,29], complete_coverage=not stage,
                sources=sources, recursive_manifests=len(B.manifests), source_references=B.references,
                baseline_and_vertex_bounds=bounds, old_reductions=old, dotted=dotted, finite=finite,
                literal_map_checks=B.maps, zero_touch_rows=len(B.payments),
                dimension_feasible_rows=dict(sorted(row_counts.items())),
                dimension_feasible_cases=dict(sorted(case_counts.items())),
                dimension_over_rank_skips=dict(sorted(B.dimension_skips.items())), payments=B.payments,
                scope='Complete-mode proof composition for compact arithmetic hyperbolic Coxeter polytopes with ground field exactly Q(sqrt5). The preceding complete binder excludes28/29. All current-rule old finite, dotted and complete finite-queue singular stops are recomputed for25/26/27; no previous payment is assumed to transfer. Arithmetic queue completeness and field gates are independently reviewed dependencies, not recomputed here. Stage mode has no dimension conclusion.')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--stage',action='store_true');a=p.parse_args()
    assert not(a.verify and a.stage)
    out=json.loads(json.dumps(run(a.stage)))
    path=ROOT/f'results/{TAG}.json'
    if a.stage:
        print(out['status'],out['dimension_feasible_cases'])
    elif a.verify:
        assert out==json.loads(path.read_text());print('PASS saved complete25-through29 coverage replay')
    else:
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
