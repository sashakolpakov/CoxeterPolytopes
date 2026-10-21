#!/usr/bin/env python3
"""Complete Qsqrt5 dimension22 exclusion with three normalized hub families."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_dimension27_coverage_20261015 as C
import fresh_sqrt5_dimension22_h3_cross_20261016_review as H

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_h3_hub_coverage_20261016'
PREVIOUS = 'fresh_sqrt5_hub_count_coverage_20261016_bound22'
BOUND = 22
CAP = 2


class Binder(C.Binder):
    def __init__(self):
        super().__init__()
        self.paths.update([f'scripts/{TAG}.py', f'scripts/{H.TAG}.py'])
        primary = self.load(H.PRIMARY + ('_cap' + str(CAP) if CAP else '_exact'))
        self.hub_review = self.review(H.TAG)
        assert self.hub_review['status'] == 'PASS_INDEPENDENT_EXACT_THREE_HUB_FACE_AND_FOREST_REVIEW'
        assert self.hub_review['dimension'] == BOUND
        assert primary['I25_count_context_cap'] == 2 and primary['H3_count_context_cap'] == 1
        assert self.hub_review['sources'][f'results/{H.PRIMARY}' + ('_cap' + str(CAP) if CAP else '_exact') + '.json'] == sha256(
            (ROOT / (f'results/{H.PRIMARY}' + ('_cap' + str(CAP) if CAP else '_exact') + '.json')).read_bytes()).hexdigest()
        self.rule = H.Rule(primary)
        self.over_rank = 0

    def payment(self, row, origin):
        if row['nodes'] - 2 > BOUND:
            self.over_rank += 1
            cases = []
        else:
            cases = [dict(dimension=BOUND, external_A1=r['external_A1'],
                          external_I25=r['external_I25'], external_H3=r['external_H3'], corners=r['corner_weights'], total=r['payment'])
                     for r in self.rule.payment(row)]
            remaining = BOUND-row['nodes']+2
            wanted = {(a,b,c) for a in range(remaining+1)
                      for b in range((remaining-a)//2+1)
                      for c in range((remaining-a-2*b)//3+1)
                      if remaining-a-2*b-3*c != 1}
            actual = {(r['external_A1'],r['external_I25'],r['external_H3']) for r in cases}
            assert cases and len(actual)==len(cases) and actual==wanted
            for r in cases:
                values=list(map(F,r['corners']))
                assert len(values)==4 and min(values)>=0 and sum(values)==F(r['total'])>=1
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
    assert maximum == F(B.hub_review['vertex_maximum'])
    paired = [F((BOUND-i)*(BOUND-i-1)+i*(i-1),4) for i in range(BOUND+1)]
    minimum = min(paired)
    assert paired[0] > minimum
    threshold = 5*minimum-F(BOUND*(BOUND-1),2)
    assert threshold == H.TARGET and maximum <= threshold
    return dict(triangle_markings=76, baseline_totals=triangle_totals,
                dimension=BOUND, component_types=len(models), exact_hub_count_contexts=len(fronts),
                forest_maxima=fronts, maximum_vertex_weight=str(maximum),
                minimum_paired_f2_coefficient=str(minimum), strict_face_weight_threshold=str(threshold),
                margin=str(threshold-maximum), strictness='h_0=h_n=1 and the extreme paired coefficient is larger')


def run(stage=False):
    if not stage:
        prior=json.loads((ROOT/f'results/{PREVIOUS}.json').read_text())
        assert prior['status'].startswith('PASS_COMPLETE_QSQRT5_FACE_COVERAGE_EXCLUDES_COMPACT_DIMENSIONS')
        assert prior['complete_coverage'] and prior['arithmetic_ground_field']=='Q(sqrt5)'
        assert prior['compact_dimension_upper_bound']==BOUND and prior['excluded_dimensions']==list(range(BOUND+1,30))
    B=Binder()
    if not stage:B.load(PREVIOUS)
    vertex=bounds(B)
    old,roots,dotted_ids=C.old_branches(B)
    dotted=C.dotted_branch(B,dotted_ids)
    finite=None if stage else C.finite_branch(B,roots)
    sources=B.bind_manifests()
    return dict(status='STAGED_HUB_PAYMENTS_ONLY_NO_DIMENSION_CLAIM' if stage else
                f'PASS_COMPLETE_QSQRT5_FACE_COVERAGE_EXCLUDES_COMPACT_DIMENSIONS{BOUND}_THROUGH29',
                arithmetic_ground_field='Q(sqrt5)',compact_dimension_upper_bound=None if stage else BOUND-1,
                excluded_dimensions=[] if stage else list(range(BOUND,30)),complete_coverage=not stage,
                dimension_newly_excluded=BOUND,I25_count_context_cap=CAP,H3_count_context_cap=1,sources=sources,recursive_manifests=len(B.manifests),source_references=B.references,
                baseline_and_vertex_bounds=vertex,old_reductions=old,dotted=dotted,finite=finite,
                literal_map_checks=B.maps,zero_touch_rows=len(B.payments),
                dimension_feasible_rows=sum(bool(r['cases']) for r in B.payments),
                dimension_feasible_cases=sum(len(r['cases']) for r in B.payments),
                over_rank_rows=B.over_rank,payments=B.payments,
                scope='Complete-mode composition for compact arithmetic hyperbolic Coxeter polytopes with ground field exactly Q(sqrt5). The preceding complete binder excludes all higher dimensions through29. Every older and complete-queue zero-touch row fitting the target rank is reevaluated with common exact external A1/I2(5)/H3 counts. All-count hub DP and independently reviewed arithmetic queues are explicit dependencies. No noncompact, other-field, nonarithmetic or construction claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--stage',action='store_true');a=p.parse_args()
    assert not(a.verify and a.stage)
    assert H.BOUND==BOUND and H.TARGET==44
    out=json.loads(json.dumps(run(a.stage)))
    path=ROOT/f'results/{TAG}_bound{BOUND-1}.json'
    if a.stage:print(out['status'],out['dimension_feasible_cases'])
    elif a.verify:assert out==json.loads(path.read_text());print('PASS saved complete lower-dimension coverage replay')
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
