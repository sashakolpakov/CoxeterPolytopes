#!/usr/bin/env python3
"""Complete Qsqrt5 rank24 face coverage using normalized I2(5) cross weights."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_dimension27_coverage_20261015 as C
import fresh_sqrt5_dimension24_i25_cross_20261015_review as H

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_dimension23_coverage_20261015'
PREVIOUS = 'fresh_sqrt5_dimension24_coverage_20261015'
BOUND = 24


class Binder(C.Binder):
    def __init__(self):
        super().__init__()
        self.paths.update([f'scripts/{TAG}.py', f'scripts/{H.TAG}.py', f'scripts/{H.V.TAG}.py'])
        primary = self.load(H.PRIMARY + '_hub_context')
        self.hub_review = self.review(H.TAG)
        assert self.hub_review['status'] == 'PASS_INDEPENDENT_EXACT_DIMENSION24_HUB_CONTEXT_FACE_AND_FOREST_REVIEW'
        assert self.hub_review['dimension'] == BOUND
        assert self.hub_review['sources'][f'results/{H.PRIMARY}_hub_context.json'] == sha256(
            (ROOT / f'results/{H.PRIMARY}_hub_context.json').read_bytes()).hexdigest()
        self.rule = H.Rule(primary)
        self.over_rank = 0

    def payment(self, row, origin):
        if row['nodes'] - 2 > BOUND:
            self.over_rank += 1
            cases = []
        else:
            cases = [dict(dimension=BOUND, external_A1=r['external_A1'],
                          external_I25=r['external_I25'], corners=r['corner_weights'], total=r['payment'])
                     for r in self.rule.payment(row)]
        self.payments.append(dict(origin=origin, nodes=row['nodes'], cases=cases))


def bounds(B):
    triangle = B.review('fresh_sqrt5_triangle_weights_20261015_review')
    assert len(triangle['markings']) == 76
    triangle_totals = []
    for rec in triangle['markings']:
        edges = [(i,j,m) for (i,j),m in zip(combinations(range(rec['order']),2),rec['labels'])]
        total = F(0)
        for omit in rec['black']:
            marked = [i for i in rec['black'] if i != omit]
            parts = H.R.components(H.R.graph([i for i in range(rec['order']) if i != omit],edges))
            g = next(h for h in parts if marked[0] in h)
            total += H.R.baseline(g,*marked) if marked[1] in g else 0
        assert total >= 2
        triangle_totals.append(str(total))
    models, totals, fronts, multiplicities = H.component_bounds(B.rule)
    assert fronts == B.hub_review['forest_maxima']
    maximum = max(F(r['maximum']) for r in fronts)
    assert maximum == F(B.hub_review['vertex_maximum']) == F(267527,5000)
    paired = [F((BOUND-i)*(BOUND-i-1)+i*(i-1),4) for i in range(BOUND+1)]
    minimum = min(paired)
    assert paired[0] > minimum
    threshold = 5*minimum-F(BOUND*(BOUND-1),2)
    assert threshold == 54 and maximum < threshold
    return dict(triangle_markings=76, baseline_totals=triangle_totals,
                dimension=BOUND, component_types=len(models), exact_hub_count_contexts=len(fronts),
                forest_maxima=fronts, maximum_vertex_weight=str(maximum),
                minimum_paired_f2_coefficient=str(minimum), strict_face_weight_threshold=str(threshold),
                margin=str(threshold-maximum), strictness='h_0=h_n=1 and the extreme paired coefficient is larger')


def run(stage=False):
    if not stage:
        prior=json.loads((ROOT/f'results/{PREVIOUS}.json').read_text())
        assert prior['status']=='PASS_COMPLETE_QSQRT5_FACE_COVERAGE_EXCLUDES_COMPACT_DIMENSIONS25_THROUGH29'
        assert prior['complete_coverage'] and prior['arithmetic_ground_field']=='Q(sqrt5)'
        assert prior['compact_dimension_upper_bound']==24 and prior['excluded_dimensions']==[25,26,27,28,29]
    B=Binder()
    if not stage:B.load(PREVIOUS)
    vertex=bounds(B)
    old,roots,dotted_ids=C.old_branches(B)
    dotted=C.dotted_branch(B,dotted_ids)
    finite=None if stage else C.finite_branch(B,roots)
    sources=B.bind_manifests()
    return dict(status='STAGED_HUB_PAYMENTS_ONLY_NO_DIMENSION_CLAIM' if stage else
                'PASS_COMPLETE_QSQRT5_FACE_COVERAGE_EXCLUDES_COMPACT_DIMENSIONS24_THROUGH29',
                arithmetic_ground_field='Q(sqrt5)',compact_dimension_upper_bound=None if stage else 23,
                excluded_dimensions=[] if stage else [24,25,26,27,28,29],complete_coverage=not stage,
                sources=sources,recursive_manifests=len(B.manifests),source_references=B.references,
                baseline_and_vertex_bounds=vertex,old_reductions=old,dotted=dotted,finite=finite,
                literal_map_checks=B.maps,zero_touch_rows=len(B.payments),
                dimension_feasible_rows=sum(bool(r['cases']) for r in B.payments),
                dimension_feasible_cases=sum(len(r['cases']) for r in B.payments),
                over_rank_rows=B.over_rank,payments=B.payments,
                scope='Complete-mode composition for compact arithmetic hyperbolic Coxeter polytopes with ground field exactly Q(sqrt5). The preceding complete binder excludes25..29. Every older and complete-queue zero-touch row fitting rank24 is reevaluated with common exact external A1/I2(5) counts. All-count hub DP and independently reviewed arithmetic queues are explicit dependencies. No noncompact, other-field, nonarithmetic or construction claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--stage',action='store_true');a=p.parse_args()
    assert not(a.verify and a.stage)
    out=json.loads(json.dumps(run(a.stage)))
    path=ROOT/f'results/{TAG}.json'
    if a.stage:print(out['status'],out['dimension_feasible_cases'])
    elif a.verify:assert out==json.loads(path.read_text());print('PASS saved complete24-through29 coverage replay')
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
