#!/usr/bin/env python3
"""Integrity and saved-headline audit only; no mathematical computation replay."""
import argparse,ast,json
from fractions import Fraction as F
from pathlib import Path
from urllib.parse import unquote
import check_dimension59_distance6_checkpoint_20261018 as P
ROOT=P.ROOT;digest=P.digest;read=P.read;links=P.links
TAG='dimension59_positive_cusps_and_a14_frontier_20261018'
CHECKPOINT=f'results/{TAG}.md';AUDIT=f'results/{TAG}_audit.json'
STATUS='PASS_DIMENSION59_POSITIVE_CUSPS_CHECKPOINT_METADATA'
I='fresh_order16_cd_sparse_sixteenth_transfer_20261018'
A='fresh_finite_a14_near_endpoint_zero_transfer_20261018'
C='fresh_a14_cd16_ambient59_frontier_20261018'
T='fresh_ambient59_exact_typed_obstruction_frontier_20261018'
R='fresh_raw112_complete_marker_family_20261018'
N='fresh_oblique_left_boundary_connector_20261018'
PACKAGES=[
 (I,'PASS_ORDER16_C101_D411_OVER4_SPARSE_SIXTEENTH_TRANSFER','23602/adec0b','11344/e9d3ab'),
 (A,'PASS_FINITE_A14_NEAR_ENDPOINT_ZERO_EXISTING_A15_UNIT','68801/6f039a','90485/336604'),
 (C,'PASS_A14_CD16_AMBIENT59_POSITIVE_IDEAL_FINITE787_OVER2','57727/f8bf79','56445/7e95fe'),
 (T,'PASS_AMBIENT59_EIGHT_ACTUAL_TYPED_OBSTRUCTIONS','56e15f','2c4774'),
 (R,'PASS_RAW112_COMPLETE_EXTRA_VERTEX_MARKER_FAMILY','83347/cc1a25','69187/b92335'),
 (N,'BOUNDED_OBLIQUE_LEFT_BOUNDARY_CONNECTOR','53001/842f40','88749/e70498')]
HELPERS=['scripts/check_raw112_complete_marker_saved_20261018.py',
 'scripts/check_oblique_left_connector_saved_20261018.py',
 'results/fresh_oblique_left_boundary_connector_20261018_pilot.json']
PINS={
 'scripts/check_dimension59_distance6_checkpoint_20261018.py':'699187b70669ab1d116c0dbdcca1fd5415656468533a267c5633784b84576bcf',
 'results/dimension59_distance6_and_contextual_ideal_frontier_20261018.md':'dcd1fa5cda4519f5c829ce59884a4a588e07f78059931704cc3e9125e0a82a9e',
 'results/dimension59_distance6_and_contextual_ideal_frontier_20261018_audit.json':'aae6a3703dd40dc5f1d6b78a2467926d93499d9f3ff5b40226c0c355bf86c781'}
MUTABLE=['README.md','RESEARCH_PLAN.md','GOAL.md']


def headlines(d):
    i=d[I];assert i['new_family_costs']==dict(A='80',B='102',C='101',D='411/4')
    assert i['new_order16_maximum']=='411/4'and i['new_coupon']['large_marks']==[11]
    assert i['new_coupon']['component_pair_increment']=='1/16'
    assert i['finite_accounting']['representative_A1_D15_allowance']=='25/16'
    assert i['current_literal_interface']['changes']==dict(ideal_debits=4,finite_credits=2)
    a=d[A];assert a['new_costs']==dict(A14='89',B14='87',D14='181/2')
    assert a['rank14_representative_upper']=='181/2'and not a['new_cross_component_coupons']
    s=a['complete_selected_support'];assert(s['complete_patterns'],s['complete_maps'],s['affected_maps'])==(600,552,8)
    assert [p['index']for p in a['selected_proofs']]==[1157,1617,1640]
    assert a['current_literal_interface']['changed_corners']==0
    c=d[C];f=c['finite'];ideal=c['ideal'];ss=ideal['summary']
    assert c['ambient_dimension']==c['completed_unrestricted_bound']==59 and c['new_dimension_bound']is None
    assert f['exact_maximum']=='787/2'and f['curvature']=='783/2'
    assert f['actual_attainer']==['D14','A15','A15','A15']
    assert f['fixed_large_vectors']==len(f['alternatives'])==961
    assert f['fixed_large_and_A5_vectors']==len(f['expanded_alternatives'])==f['interaction_keys']==2574
    assert len(f['tables'])==5148 and f['Bellman_upper_inequalities']==6723288
    assert f['attaining_predecessors']==303732 and f['residual_item_checks']==128700
    assert f['simultaneous_type_pairs']==32400 and f['context_options']==9191
    assert len(f['above_curvature'])==3 and sorted(r['sixteenths_upper']for r in f['above_curvature'])==[6272,6276,6296]
    assert ideal['complete_profiles']==len(ideal['rows'])==len({tuple(r['ranks'])for r in ideal['rows']})==1890
    assert ss['negative_count']==ss['zero_count']==0 and all(F(r['normalized_margin'])>0 for r in ideal['rows'])
    assert ss['normalized_minimum']=='1/340'and ss['absolute_minimum']==ideal['uniform_actual_cusp_surplus']=='91/2'
    assert ideal['positive_terminal_extension_to_actual_cusps']and ideal['actual_charge_bounded_by_envelope']
    assert c['strict_global_constant']=='4350'and c['constant_counted_once']
    assert c['necessary_ambient59_inequality']=='2*v_f >= 4350 + (91/2)*v_i'
    t=d[T];assert t['complete_saved_expanded_vectors']==2574 and t['actual_type_forests_checked']==1164
    assert t['actual_obstruction_count']==len(t['actual_obstructions'])==8 and t['actual_rank_patterns']==2
    assert t['removed_pattern_maximum']=='781/2'and t['exact_finite_maximum']=='787/2'
    assert t['all_other_actual_forests_at_most_curvature']and not t['geometric_realizations_claimed']
    assert [(r['actual_type_forests'],r['actual_maximum'])for r in t['typed_refinements']]==[(1080,'781/2'),(30,'787/2'),(54,'1569/4')]
    r=d[R];z=r['counts'];e=r['exported'];repair=e['repair']
    assert [z[k]for k in('classes','pointed_marker_domains','target_forests','markers','relative_cases','retained_relative_cases')]==[196,392,383,11012,162,14]
    assert [z[k]for k in('exported_columns','strictly_violated_columns','raised_allowances')]==[529,3,3]
    assert len(e['preserved_columns'])==3299 and len(e['preserved_multi_owner_costs'])==14
    assert F(repair['bound'])==F(14671331978990040,1224897727063181)<12
    assert r['optimizer_invocations']==0 and r['old_outcomes_changed']is False
    n=d[N];assert n['summary']==dict(frozen_pointed_configurations=376,distinct_recipe_prefixes=376,tested_prefixes=376,
        recipes=8968,branches=17936,native_N=0,clean_N=0,completed_disks=0,fewer_than_ten_ends=0,unresolved=0)
    assert not n['accepted']and not n['unresolved']
    assert n['wall_count_breakdown']==dict(projected_originals=52,fixed_original_whites=6,original_facets=58,
        new_pre_N_walls=6,explicit_product_gates=58,total_pre_N_facets=64)
    assert n['exact_deletion_debit']['lost_finite_count']==19
    assert n['exact_deletion_debit']['compact_faces_by_dimension']=={'2':12,'3':1}
    for s in(I,A,C):assert d[s]['triangle_minimum']==2 and d[s]['quadrilateral_minimum']==1
    return dict(completed_unrestricted_bound=59,new_dimension_bound=None,
        completed_proof_has_no_identified_gap=True,finite_maximum='787/2',finite_curvature='783/2',
        positive_ideal_profiles=1890,uniform_actual_cusp_surplus='91/2',
        actual_finite_obstructions=8,actual_rank_patterns=2,removed_representative_maximum='781/2',
        arithmetic_fixed_mass_dual_ratio=repair['bound'],arithmetic_primal_or_exclusion_claim=False,
        forward_N_native=0,forward_N_branches=17936,new_finite_volume_example=False,
        external_refereeing=False,proof_assistant_formalization=False)


def collect():
    for p,h in PINS.items():assert digest(p)==h,('predecessor',p)
    # This is itself an integrity/saved-headline replay, never mathematics.
    P.verify(json.loads((ROOT/P.AUDIT).read_text()))
    docs={};records=[];bound={};total=0
    for stem,status,generation,replay in PACKAGES:
        d=read(stem);assert d['status']==status;docs[stem]=d
        records.append(dict(stem=stem,status=status,generation=generation,fresh_replay=replay,
            scope='Recorded terminal PASS evidence; no execution by this auditor.'))
        for p,h in d['sources'].items():
            if p not in bound:bound[p]=digest(p)
            assert bound[p]==h,(stem,p)
        total+=len(d['sources'])
    return docs,records,bound,total


def artifacts():
    return list(dict.fromkeys([p for s,*_ in PACKAGES for p in
        (f'scripts/{s}.py',f'results/{s}.json',f'results/{s}.md')]+HELPERS+list(PINS)+
        [CHECKPOINT,str(Path(__file__).relative_to(ROOT))]))


def check_link(r):
    p=((ROOT/r['document']).parent/unquote(r['target'].split('#',1)[0])).resolve()
    assert p.is_file()or p==(ROOT/AUDIT).resolve(),r


def build():
    docs,records,bound,total=collect();values=headlines(docs);paths=artifacts()
    asts=[p for p in paths if p.endswith('.py')]
    for p in asts:ast.parse((ROOT/p).read_text(),filename=p)
    reports=[p for p in paths if p.endswith('.md')]+MUTABLE;ls=[r for p in reports for r in links(p)]
    for r in ls:check_link(r)
    snapshots={}
    for p in MUTABLE:
        opening=(ROOT/p).read_text()[:6000];assert CHECKPOINT in opening and 'finite-volume'in opening and '59'in opening
        snapshots[p]=dict(sha256=digest(p),opening=opening)
    return dict(status=STATUS,checkpoint=CHECKPOINT,packages=records,saved_headlines=values,
        counts=dict(new_packages=len(PACKAGES),direct_bindings=total,distinct_bound_files=len(bound),python_asts=len(asts),documents=len(reports),local_links=len(ls)),
        artifact_sha256={p:digest(p)for p in paths},bound_file_sha256=bound,python_ast_files=asts,documents=reports,
        local_links=ls,mutable_snapshots=snapshots,
        scope='Integrity and saved headline checks only, including previous metadata auditor. No mathematical producer, helper, geometry, census, Bellman proof, convolution, cone or optimizer replay. Mutable openings are historical snapshots.')


def verify(saved):
    docs,records,bound,total=collect();paths=artifacts()
    assert saved['status']==STATUS and saved['checkpoint']==CHECKPOINT and saved['packages']==records and saved['saved_headlines']==headlines(docs)
    assert set(saved['artifact_sha256'])==set(paths)and saved['bound_file_sha256']==bound
    for p,h in saved['artifact_sha256'].items():assert digest(p)==h,p
    asts=[p for p in paths if p.endswith('.py')];reports=[p for p in paths if p.endswith('.md')]+MUTABLE
    assert saved['python_ast_files']==asts and saved['documents']==reports
    assert saved['counts']==dict(new_packages=len(PACKAGES),direct_bindings=total,distinct_bound_files=len(bound),python_asts=len(asts),documents=len(reports),local_links=len(saved['local_links']))
    assert [r for r in saved['local_links']if r['document']not in MUTABLE]==[r for p in reports if p not in MUTABLE for r in links(p)]
    for p in asts:ast.parse((ROOT/p).read_text(),filename=p)
    for r in saved['local_links']:check_link(r)
    assert set(saved['mutable_snapshots'])==set(MUTABLE)
    for p,r in saved['mutable_snapshots'].items():assert CHECKPOINT in r['opening']and len(r['sha256'])==64


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();p=ROOT/AUDIT
    if args.verify:out=json.loads(p.read_text());verify(out)
    else:assert not p.exists();out=build();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
