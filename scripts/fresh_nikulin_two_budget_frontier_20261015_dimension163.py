#!/usr/bin/env python3
"""Bind reviewed ordered endpoint budgets to the complete dimension163 split."""
import argparse,json
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_two_budget_frontier_20261015_dimension163'
BASE='fresh_nikulin_defect_frontier_20261015_dimension164'
PRIMARY='fresh_nikulin_two_budget_frontier_20261015'
REVIEW=PRIMARY+'_review'
FACET='fresh_nikulin_narrow_20261015_compact_facet'

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    tags=[BASE,PRIMARY,REVIEW,FACET]
    data={t:json.loads((ROOT/f'results/{t}.json').read_text())for t in tags}
    assert data[BASE]['status']=='PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_164_CANDIDATE'
    assert data[PRIMARY]['status']=='PASS_EXACT_ORDERED_RESIDUAL_BUDGET_CANDIDATE163_PENDING_INDEPENDENT_REVIEW'
    assert data[REVIEW]['status']=='PASS_INDEPENDENT_PHYSICAL_GRAMS_CONTINUOUS_ENVELOPE_ORDERED_PACKING_AND_163_SCOPE'
    for p in(f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json'):assert data[REVIEW]['sources'][p]==digest(p)
    p=data[PRIMARY];r=data[REVIEW];B=Q(21,20)
    assert Q(p['height_ratio_limit'])==B and Q(p['residual_endpoint_cap'])==Q(3,2)
    assert p['counts']==dict(models=91,retained_support_pairs=1116,continuous_domains=3630,positive_loss_domains=1395,pareto_records=28)
    assert r['finite_continuous_envelope']['counts']==dict(models=91,support_pairs=1116,domains=3630,positive=1395)
    budgets=[([Q(10,9),Q(4,3)],Q(20,9)),([Q(1),Q(3,2)],Q(22,9)),([Q(3,2),Q(1)],Q(22,9))]
    for (caps,loss),a,b in zip(budgets,p['ordered_budget_packing'],r['ordered_packing']):
        assert list(map(Q,a['endpoint_budgets']))==list(map(Q,b['budgets']))==caps
        assert Q(a['maximum_loss'])==Q(b['maximum'])==loss
    assert len(p['same_row_E7_exclusions'])==len(r['actual_failing_factor_split']['same_row_exclusions'])==4
    assert len(p['remaining_failure_floors'])==len(r['actual_failing_factor_split']['remaining_floors'])==13
    assert Q(16,9)+Q(20,9)==Q(14,9)+Q(22,9)==4
    assert Q(p['aggregate']['total_loss'])==4
    # N=n-1. One prescribed omitted facet costs2/9 in the target sum.
    near=(Q(2*163,9)-Q(2,9)-4)/2-2
    far=B*Q(163-43,9)
    assert near==far==14
    assert Q(1,9)>0 and B/9==Q(7,60)>0
    assert p['aggregate']['near_pairing']=='M*(N-37)/9'
    assert p['aggregate']['far_pairing']=='M*(N-43)/9'
    assert p['aggregate']['ideal_branch_dimension_bound']==r['global_threshold']['ideal_dimension']==163
    assert p['aggregate']['compact_branch_inherited_bound']==r['global_threshold']['compact_inherited_dimension']==data[FACET]['ambient_dimension_at_most']==143
    assert p['aggregate']['combined_candidate_bound']==163
    assert r['global_threshold']['near']=='t/9 + 14'and r['global_threshold']['far']=='7*t/60 + 14'
    assert (Q(2*162,9)-Q(2,9)-4)/2-2<14 and B*Q(162-43,9)<14
    paths=[f'scripts/{TAG}.py',f'results/{TAG}.md']
    for tag in tags:paths +=[f'scripts/{tag}.py',f'results/{tag}.json',f'results/{tag}.md']
    seen,checked=set(),set()
    def recurse(path):
        if path in seen:return
        seen.add(path)
        for source,h in json.loads((ROOT/path).read_text()).get('sources',{}).items():
            assert digest(source)==h,(path,source);checked.add(source)
            if source.endswith('.json'):recurse(source)
    for tag in tags:recurse(f'results/{tag}.json')
    return dict(status='PASS_REVIEWED_DEPENDENCIES_AND_CASE_SPLIT_FOR_DIMENSION_163_CANDIDATE',sources={p:digest(p)for p in paths},
        recursive_manifests_checked=len(seen),recursive_source_files_checked=len(checked),
        candidate_statement='Every finite-sided finite-volume real hyperbolic Coxeter polytope has dimension at most163.',
        arithmetic_assumption=False,compactness_assumption=False,simple_cusp_assumption=False,
        chosen_narrow_facet_cases=[dict(case='contains an ideal vertex',ambient_dimension_at_most=163,primary=PRIMARY,independent_review=REVIEW),
            dict(case='only finite vertices, hence compact',ambient_dimension_at_most=143,primary=FACET,independent_review='fresh_dimension225_hostile_20261015')],
        common_height_ratio_split='21/20',ordered_endpoint_cost_budgets=['2','2'],total_energy_loss_at_most='4',
        residual_packing=[dict(budgets=list(map(str,caps)),loss=str(loss))for caps,loss in budgets],avoid_node_deficit='2/9',
        near_pairing_lower='M*((n-1)-37)/9',far_pairing_lower='M*((n-1)-43)/9',
        first_excluded_dimension=164,combined_candidate_dimension_at_most=163,norm_two_narrow_threshold=14,
        threshold=dict(N=163,near='14',far='14',near_slope='1/9',far_slope='7/60'),
        frozen_baseline=BASE,primary_sources=data[BASE]['primary_sources'],
        scope='Reviewed local-proof dependency and complete case binder. Geometric narrowness/framing/common-height/Schur-budget and compact-facet facts remain inherited mathematical dependencies. Internal research proof candidate; external validation, novelty and current-best status remain unestablished. Frozen primary pending-review wording is historical and superseded here. No dimension162 claim.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status']);print(result['recursive_manifests_checked'],result['recursive_source_files_checked'])
