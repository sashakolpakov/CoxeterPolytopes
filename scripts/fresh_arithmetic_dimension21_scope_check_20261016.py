#!/usr/bin/env python3
"""Independent standard-library check of final quantifiers and field composition.

This does not reprove Gram arithmetic, orbit matching, queue exhaustiveness,
or the finite-forest dynamic programs. It checks the final saved composition
against its literal obligation identities, all spectator triples, rational
payment sums, strict h-vector threshold, and exact field/dimension scope.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_arithmetic_dimension21_scope_check_20261016'
COVER = 'fresh_sqrt5_h3_hub_coverage_20261016_bound21'
COMPOSE = 'fresh_arithmetic_dimension21_composition_20261016_bound21'
PRIOR = 'fresh_sqrt5_hub_count_coverage_20261016_bound22'
REVIEW = 'fresh_sqrt5_dimension22_h3_cross_20261016_review'
N = 22


class Inputs:
    def __init__(self):
        self.paths = {f'scripts/{TAG}.py'}

    def read(self, tag):
        path = f'results/{tag}.json'
        self.paths.add(path)
        return json.loads((ROOT/path).read_text())

    def manifests(self):
        checked, references, hashes = set(), 0, {}
        def digest(path):
            if path not in hashes:
                hashes[path] = sha256((ROOT/path).read_bytes()).hexdigest()
            return hashes[path]
        def visit(path):
            nonlocal references
            if path in checked or not path.endswith('.json'):
                return
            checked.add(path)
            doc = json.loads((ROOT/path).read_text())
            for key in ('sources', 'source_hashes'):
                for source, expected in doc.get(key, {}).items():
                    assert isinstance(expected, str)
                    assert digest(source) == expected, (path, source)
                    references += 1
                    visit(source)
        for path in sorted(self.paths):
            visit(path)
        return ({p:digest(p) for p in sorted(self.paths)},
                dict(manifests=len(checked), references=references,
                     unique_hashed_dependencies=len(hashes)))


def triples(rank):
    """Derive possible exact counts independently as a weighted-rank equation."""
    assert rank >= 0
    result = set()
    for c in range(rank//3+1):
        for b in range((rank-3*c)//2+1):
            for a in range(rank-3*c-2*b+1):
                residual = rank-a-2*b-3*c
                if residual == 0 or residual >= 2:
                    result.add((a,b,c))
    return result


def origin_key(origin):
    return json.dumps(origin, sort_keys=True, separators=(',', ':'))


def obligations(inputs):
    expected = {}
    def add(origin, nodes):
        key = origin_key(origin)
        assert key not in expected
        expected[key] = nodes
    base = inputs.read('fresh_sqrt5_quad_exhaustion_20261015')['labelled_union_supports']
    contacts = inputs.read('fresh_sqrt5_quad_contacts_20261015')
    for rec in contacts['records']:
        if (rec['status']=='UNPAID_FINITE_CONTACT_DOMAIN' and
                not rec['accepted'] and rec['singular']):
            add(['old_sealed',rec['class_index']],base[rec['union_index']]['nodes'])
    first = inputs.read('fresh_sqrt5_sealed_subcore_closure_20261015')
    for rec in first['records']:
        if rec['singular']:
            add(['singular_first_extension',rec['source_class_index'],rec['trial_index']],
                base[rec['source_union_index']]['nodes']+1)
    for tag, name in [
            ('fresh_sqrt5_singular_extension_frontier_20261015','singular51'),
            ('fresh_sqrt5_singular_attachment_closure_20261015','singular_closure17')]:
        for rec in inputs.read(tag)['records']:
            add([name,rec['index']],rec['row']['nodes'])
    dotted = inputs.read('fresh_sqrt5_dotted_rank29_exhaustion_20261015')
    assert dotted['processed']==dotted['marked_states']==len(dotted['records'])==438
    dotted_ids = []
    for rec in dotted['records']:
        if rec['result']['singular']:
            dotted_ids.append(rec['index'])
            add(['dotted',rec['index'],rec['ratio']],rec['row']['nodes'])
    assert dotted_ids==dotted['weight_obligations'] and len(dotted_ids)==264
    finite = inputs.read('fresh_sqrt5_nonsingular_rank29_fast_20261015')
    assert finite['complete_queue'] and not finite['counts'].get('PENDING',0)
    assert finite['state']['next_core']==finite['marked_states']==len(finite['records'])==5172
    finite_ids = []
    for rec in finite['records']:
        if rec['result']['singular']:
            finite_ids.append(rec['index'])
            add(['finite_queue',rec['index']],rec['row']['nodes'])
    assert finite_ids==finite['weight_obligations'] and len(finite_ids)==2070
    assert len(expected)==2661
    return expected


def run():
    inputs = Inputs()
    coverage, composition = inputs.read(COVER), inputs.read(COMPOSE)
    previous, review = inputs.read(PRIOR), inputs.read(REVIEW)
    assert coverage['status']=='PASS_COMPLETE_QSQRT5_FACE_COVERAGE_EXCLUDES_COMPACT_DIMENSIONS22_THROUGH29'
    assert coverage['complete_coverage'] and coverage['arithmetic_ground_field']=='Q(sqrt5)'
    assert coverage['compact_dimension_upper_bound']==21
    assert coverage['excluded_dimensions']==list(range(22,30))
    assert coverage['dimension_newly_excluded']==N
    assert previous['complete_coverage'] and previous['arithmetic_ground_field']=='Q(sqrt5)'
    assert previous['compact_dimension_upper_bound']==22
    assert previous['excluded_dimensions']==list(range(23,30))
    for tag in (PRIOR, REVIEW):
        path=f'results/{tag}.json'
        assert coverage['sources'][path]==sha256((ROOT/path).read_bytes()).hexdigest()
    assert review['status']=='PASS_INDEPENDENT_EXACT_THREE_HUB_FACE_AND_FOREST_REVIEW'
    assert review['dimension']==22 and review['face_cases']==89255

    expected = obligations(inputs)
    actual, origin_counts, by_rank = {}, Counter(), Counter()
    cases = feasible = excluded = 0
    minimum_payment = None
    for row in coverage['payments']:
        key=origin_key(row['origin'])
        assert key not in actual
        actual[key]=row['nodes']
        origin_counts[row['origin'][0]]+=1
        rank=row['nodes']-2
        if rank>N:
            assert not row['cases']
            excluded+=1
            continue
        feasible+=1
        want=triples(N-rank)
        got=set()
        assert row['cases']
        for case in row['cases']:
            assert case['dimension']==N
            triple=tuple(case[k] for k in ('external_A1','external_I25','external_H3'))
            assert all(type(x) is int and x>=0 for x in triple)
            assert triple not in got
            got.add(triple)
            corners=list(map(F,case['corners']))
            assert len(corners)==4 and min(corners)>=0
            total=sum(corners,F(0))
            assert total==F(case['total']) and total>=1
            minimum_payment=total if minimum_payment is None else min(minimum_payment,total)
            cases+=1
        assert got==want,(row['origin'],N-rank)
        by_rank[rank]+=len(got)
    assert actual==expected
    assert dict(origin_counts)==dict(old_sealed=112,singular_first_extension=147,
        singular51=51,singular_closure17=17,dotted=264,finite_queue=2070)
    assert len(actual)==coverage['zero_touch_rows']==2661
    assert feasible==coverage['dimension_feasible_rows']==1634
    assert excluded==coverage['over_rank_rows']==1027
    assert cases==coverage['dimension_feasible_cases']==91726
    assert coverage['literal_map_checks']==471

    bounds=coverage['baseline_and_vertex_bounds']
    assert bounds['dimension']==N and bounds['triangle_markings']==76
    assert len(bounds['baseline_totals'])==76
    assert min(map(F,bounds['baseline_totals']))>=2
    maxima={tuple(x['context']):F(x['maximum']) for x in bounds['forest_maxima']}
    assert len(maxima)==len(bounds['forest_maxima'])==410
    assert set(maxima)==triples(N)
    assert bounds['forest_maxima']==review['forest_maxima']
    maximum=max(maxima.values())
    assert maximum==F(bounds['maximum_vertex_weight'])==F(review['vertex_maximum'])==F(21977,500)
    paired=[F(comb(N-i,2)+comb(i,2),2) for i in range(N+1)]
    minimum=min(paired)
    assert minimum==55 and paired[0]>minimum and paired[-1]>minimum
    threshold=5*minimum-comb(N,2)
    assert threshold==44 and maximum<threshold
    assert F(bounds['minimum_paired_f2_coefficient'])==minimum
    assert F(bounds['strict_face_weight_threshold'])==threshold
    assert F(bounds['margin'])==threshold-maximum==F(23,500)

    assert composition['status']=='PASS_ARITHMETIC_REFLECTION_DIMENSION_COMPOSITION'
    assert composition['compact_arithmetic_dimension_upper_bound']==21
    assert composition['finite_volume_arithmetic_reflection_dimension_upper_bound']==21
    expected_branches={(n,field,upper) for n in range(22,30) for field,upper in
        [('Q(sqrt2)',20),('Q(sqrt5)',21),('Q(cos(2pi/7))',20)]}
    actual_branches=[(r['dimension'],r['field'],r['excluded_by_upper_bound'])
                     for r in composition['compact_high_dimension_field_branches']]
    assert len(actual_branches)==len(set(actual_branches))==24
    assert set(actual_branches)==expected_branches
    assert all(n>=22 and n>upper for n,field,upper in actual_branches)
    assert composition['noncompact_arithmetic']==dict(upper_bound=21,excluded_dimension=20,dimension21_attained=True)
    assert len(composition['published_inputs'])==3
    assert composition['sources'][f'results/{COVER}.json']==sha256((ROOT/f'results/{COVER}.json').read_bytes()).hexdigest()
    sources, manifest_summary = inputs.manifests()
    return dict(status='PASS_INDEPENDENT_FINAL_ARITHMETIC_DIMENSION21_SCOPE_AND_QUANTIFIER_CHECK',
        sources=sources,manifest_summary=manifest_summary,
        obligation_rows=len(actual),origin_counts=dict(sorted(origin_counts.items())),
        rank_feasible_rows=feasible,over_rank_rows=excluded,common_external_triples=cases,
        cases_by_core_corner_rank={str(k):v for k,v in sorted(by_rank.items())},
        minimum_saved_quadrilateral_payment=str(minimum_payment),
        complete_vertex_count_contexts=len(maxima),maximum_saved_vertex_weight=str(maximum),
        hvector=dict(minimum_paired_coefficient=str(minimum),threshold=str(threshold),
                     margin=str(threshold-maximum),strictness='h_0=h_n=1 and both extreme coefficients exceed the minimum'),
        field_dimension_branches=actual_branches,finite_volume_arithmetic_dimension_upper_bound=21,
        scope='Independent saved-composition scope check using only Python standard-library exact rational arithmetic. Every obligation identity/node count, common external triple, saved payment sum, h-vector threshold and field/dimension branch is checked. Gram arithmetic, queue exhaustiveness, graph weight evaluation, forest dynamic programming and published theorems remain explicit imported dependencies; this check does not reprove them. No external validation, formal proof-assistant verification or priority claim.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=json.loads(json.dumps(run()))
    output=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result==json.loads(output.read_text())
    else:
        output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['common_external_triples'],len(result['field_dimension_branches']))
