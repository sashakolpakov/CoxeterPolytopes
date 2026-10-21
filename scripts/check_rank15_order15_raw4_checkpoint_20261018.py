#!/usr/bin/env python3
"""Checkpoint integrity and saved status only; no proof or numerical replay.

No producer/helper is imported or executed. Prior execution identifiers are
recorded evidence. Mutable README/plan hashes are historical snapshots.
"""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/rank15_order15_and_raw4_pair_barrier_20261018.md'
AUDIT=CHECKPOINT.replace('.md','_audit.json')
STATUS='PASS_RANK15_ORDER15_RAW4_CHECKPOINT_INTEGRITY_SAVED_STATUS_SYNTAX_LINKS'
PREDECESSOR='fresh_order11_order16_dimension60_20261018'
PACKAGES=[
 ('fresh_rank15_sparse_payment_exchange_20261018','PASS_RANK15_SPARSE_PAYMENT_EXCHANGE','0ca025','0a64d9',
  'Two peer full source/report mathematical reads including inherited end-projection, D15 proper/whole and A16-half interfaces. No reviewer map, inertia, census, literal, Bellman or scalar replay.'),
 ('fresh_order15_bd_terminal_coupon_transfer_20261018','PASS_ORDER15_BD92_TYPED_UNIT_TRANSFER','2941ed','9431c5',
  'Root and peer full source/report/helper mathematical reads: arbitrary contacts before zeros, exterior order14, terminal congruence, sealed typed D coupons and complete inherited interfaces. No reviewer map, inertia, census or scalar replay.'),
 ('fresh_rank15_order15_ambient60_frontier_20261018','PASS_RANK15_ORDER15_AMBIENT60_FRONTIER','71fb95','589cdf',
  'Peer full source/report mathematical composition read: simultaneous local rules, all charges, complete finite certificates, low9 domain and conditional scope. No reviewer DP, convolution, scalar or geometric census replay.'),
 ('fresh_compact10_edge_subdivision_20261018','BOUNDED_FINITE_COMPACT10_EDGE_SUBDIVISION','cdc8dc','f79a4d',
  'Root full producer/report/helper scope read including inherited bank and permanent-cusp interfaces. Helper231c6b is separate saved-data evidence; no reviewer native arithmetic, normal, cone/classifier or saved-data replay.'),
 ('fresh_new5121_raw4_all89_cover_barrier_20261018','PASS_ALL89_RAW4_CORNER_MASS_AND_PAIR_ONLY_CARDINALITY_BARRIER','09ab3f','b08afb',
  'Root full producer/report/helper mathematical read of provenance, singleton feasibility, cover/slack transport and cardinality upper/attainment scope. Helper3d871c is separate evidence; no peer review or reviewer projection, embedding, arithmetic or cardinality replay claimed.'),
 ('fresh_raw4_all89_retained_pair_barrier_20261018','PASS_ALL89_CURRENT_PAIR_GATES_RETAIN_ABOVE12_WITNESS','6da7c1','9a3eac',
  'Root full producer/report/helper mathematical read. Helper58bd31 is separate saved-data evidence. One owned unknown is retained by the stated gates, not asserted physically realizable; no reviewer native, graph, LP or arithmetic replay.'),
 ('fresh_joint10_30_native_replacement_20261018','BOUNDED_JOINT10_30_NATIVE_REPLACEMENT','29384c','8d426e',
  'Root full producer/report/helper internal mathematical and scope read including native/profile and scoring interfaces. Helperb49fcd is separate saved-data evidence; no reviewer arithmetic, norm solves, classification or cone reconstruction.'),
]
HELPERS=[
 'scripts/check_compact10_edge_subdivision_saved_20261018.py',
 'scripts/check_new5121_raw4_all89_barrier_saved_20261018.py',
 'scripts/check_raw4_all89_retained_pair_saved_20261018.py',
 'scripts/check_joint10_30_native_saved_20261018.py',
]
FROZEN_PREDECESSORS={
 'scripts/fresh_order11_order16_dimension60_20261018.py':'89d84a0f82d21fdd7c73e57fceadcc5d9859417d85be562b57617e6ef659c671',
 'results/fresh_order11_order16_dimension60_20261018.json':'019c8d9bca5dad4de0e4035acfae154f05b6f96836745370e660eaecc9add6a9',
 'results/fresh_order11_order16_dimension60_20261018.md':'9434e6c00eeb5bb52023aaf5c6411cffecd1d711bad1d99044f1b1a08c8ddeeb',
 'scripts/check_deep_finite7_raw4_checkpoint_20261018.py':'a72094b0ab94ac04736d24243a6600cbd443c3e51cd7c4a840015ee339ce4c79',
 'results/deep_finite7_and_raw4_cover_frontier_20261018.md':'519ad25a609a39340f79e10ec4ac914447fd407e182f0e58a027f8489737bd93',
 'results/deep_finite7_and_raw4_cover_frontier_20261018_audit.json':'45e31b209bd44bb8c793d998dbee7525f8dd68b6257b81802cb7600a6885eb34',
}
MUTABLE=['README.md','RESEARCH_PLAN.md']


def digest(name):return sha256((ROOT/name).read_bytes()).hexdigest()
def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())


def links_exist(document,target):
    path=(ROOT/document).parent/unquote(target.split('#',1)[0])
    assert path.exists() or path.resolve()==(ROOT/AUDIT).resolve(),(document,target)


def headlines(docs):
    p=docs[PREDECESSOR]
    assert p['status']=='PASS_COMPLETE_ORDER11_ORDER16_DIMENSION60_COMPOSITION'
    assert (p['dimension_bound'],p['excluded_dimension'],p['inherited_dimension_bound'])==(60,61,61)
    assert p['finite_maximum']==p['finite_curvature']=='420'
    assert p['actual_ideal_absolute_surplus']=='64' and p['strict_global_constant']=='4650'
    assert p['contradiction']=='0 >= 0*v_f + 64*v_i + 4650 > 0'
    assert not any(p[k] for k in ('arithmeticity_required','compactness_required','simple_cusps_required',
        'facet_count_bound_required','external_refereeing','formal_proof_assistant_verified','new_finite_volume_example'))
    a,b,c,e,u,r,j=[docs[x[0]] for x in PACKAGES]
    costs=dict(A15='101',B15='97',D15='102',A17='116')
    assert a['new_costs']==costs and a['A17_restored_unit_pair']==[1,15]
    assert a['upgraded_coupon']['component_pair_cost_increase']=='2'
    assert a['upgraded_coupon']['new_total_component_pair_cost']=='5'
    assert a['literal_configurations']==216 and a['literal_debited_corners']=={}
    assert a['finite_rank60_count']['exact_maximum']=='413' and a['ideal_weights_unchanged']
    assert b['new_family_costs']==dict(A=75,B=92,C=91,D=92)
    assert (b['B_affected_contexts'],b['B_affected_maps'],len(b['D_whole_variants']))==(53,108,20)
    assert b['B_outcomes']==dict(forced_entry_projection_exclusion=68,forced_completed_index_two=40)
    assert b['D_role_counts']==dict(A4=12,A5=8) and not b['uses_finite_same_component_units']
    assert [(x['small_type'],x['small_role'],x['component_pair_cost_increment']) for x in b['new_cross_component_coupons']]==[
        ('whole A4','endpoint','2'),('whole A5','second_from_endpoint','2')]
    assert not b['finite_forest_bound_proved'] and b['finite_budget_requires_new_coupon_increments']
    assert c['current_exceptional_finite_costs']==costs and c['current_ideal_order15_costs']==b['new_family_costs']
    f=c['finite'];s=c['ideal']['summary']
    assert (f['exact_maximum'],f['curvature'],f['complete_fixed_vectors'],f['interaction_keys'])==('413','405',1044,792)
    assert (f['new_A2_D13_actual_cost'],f['new_A4_D11_actual_cost'],f['new_A5_D11_actual_cost'])==('5','8','19/2')
    assert f['actual_attainer']==['D17','D17','D26'] and f['above_curvature']==15
    assert f['new_tables_generated'] and not f['replay_optimizes']
    assert (c['ideal']['complete_profiles'],s['negative_count'],s['zero_count'])==(2113,49,0)
    assert s['normalized_minimum']=='-1823/11560' and s['absolute_minimum']=='-123905'
    assert not c['ideal']['negative_terminal_absolute_margin_is_universal_bound']
    assert c['ideal']['remaining_conditional_target']['summary']['normalized_minimum']=='73/104000'
    assert c['strict_global_constant']=='4500' and c['constant_counted_once']
    for x in (a,b):assert x['triangle_minimum']==2 and x['quadrilateral_minimum']==1
    for x in (a,b,c):assert x['new_dimension_bound'] is None
    assert e['deleted_walls']==[30] and e['retained_old_walls']==62
    assert e['exact_deletion_debit']['lost_finite_count']==13
    assert e['exact_deletion_debit']['compact_faces_by_dimension']=={'2':10,'3':1}
    assert e['counts']==dict(norm_equations=106,norm_branches=212)
    assert [x['complete_low_recipe_count'] for x in e['banks']]==[64,42]
    assert [len(x['caps']) for x in e['banks']]==[1,1]
    assert e['selected_R']==[] and e['complete_pair_count']==0 and e['completed_connector_pairs']==[]
    assert e['topology']['removed_edge']==[10,67] and len(e['topology']['triangles'])==4
    assert len(e['boundary_ledgers'])==3
    assert (u['actual_represented_classes'],u['actual_represented_forests'])==(324,600)
    assert u['counts']==dict(actual_corners=534,additional_forests=55,classes=89,
        enumerated_component_subsets=318,existing_forests=183,independent_face_components=64,
        largest_face_component=7,provenance_kinds=dict(explicit=18,unique_saved_context=71),target_forests=238)
    assert u['minimum_cardinality_permitting_ratio_below12']==42
    assert u['all89']['candidate_pair_tests_upper']==332054
    assert u['best_minimum_batch']['candidate_pair_tests_upper']==220350
    assert u['all89']['ratio_lower_bound']=='8222765954018112/685239135430273'
    assert u['new_physical_pair_tests']==u['covering_LPs']==u['global_model_optimizations']==u['new_clique_enumerations']==0
    assert (r['authoritative_classes'],r['authoritative_forests'],r['pair_tests'])==(324,600,1)
    assert r['selected_face']==14368 and r['selected_vertex']==495
    assert r['pair_outcome']['new_exclusions']==0 and r['pair_outcome']['new_unknown_pairs']==1
    assert r['pair_outcome']['newly_owned_pair_domain']==[[306,315,None]]
    assert r['pair_cover_mass']==r['saving_after_multiplier_slack']=='1034543/2500000000'
    assert r['ratio_lower_bound']=='8222765954018112/684963831123457'
    assert r['strict_numerator_gap_over12']=='266665044719/55440000000000'
    assert r['old_graph_edges_sha256']==r['baseline']['graph_edges_sha256']
    assert r['covering_LPs']==r['global_model_optimizations']==r['new_clique_enumerations']==0
    assert u['new_dimension_bound'] is r['new_dimension_bound'] is None
    assert j['deleted_walls']==[10,30] and j['retained_old_wall_count']==61
    assert j['exact_deletion_debit']['lost_finite_count']==75 and len(j['lost_old8_finite_supports'])==44
    assert j['exact_deletion_debit']['compact_faces_by_dimension']=={'2':313,'3':325,'4':244,'5':129,'6':46,'7':10,'8':1}
    assert j['counts']==dict(new10_norm_equations=288,new10_norm_branches=576)
    assert [x['complete_low_recipe_count'] for x in j['new10_banks']]==[96,96,96]
    assert [len(x['caps']) for x in j['new10_banks']]==[0,0,1]
    assert j['complete_native_pair_count']==len(j['coherent_pairs'])==1
    control=j['coherent_pairs'][0]
    assert control['old30_restored'] and control['old10_restored'] and control['open_ends']==10
    assert control['finite_rays']==26 and control['all_old_compact_S6_faces_rebuilt']
    assert not control['changed44_quadrilateral_complete']
    assert len(j['unresolved'])==4
    return dict(completed_unrestricted_dimension_bound=60,new_dimension_bound=None,
        current_finite_rank60_maximum='413',finite_curvature='405',ideal_terminal_profiles=2113,negative_ideal_profiles=49,
        current_ideal_order15_costs=b['new_family_costs'],current_arithmetic_classes=324,current_arithmetic_forests=600,
        all89_singleton_ratio_lower_bound=u['all89']['ratio_lower_bound'],
        retained_pair_ratio_lower_bound=r['ratio_lower_bound'],retained_pair_is_physical_realization_claim=False,
        compact10_completed_connector_pairs=0,joint10_30_only_old_control=True,
        new_finite_volume_example=False,external_refereeing=False,formal_proof_assistant_verified=False)


def collect():
    assert all('PENDING' not in x for row in PACKAGES for x in row),'Await final packet run/review records.'
    for p,h in FROZEN_PREDECESSORS.items():assert digest(p)==h,('predecessor',p)
    docs={PREDECESSOR:read(PREDECESSOR)};bound={};packages=[];bindings=0
    for stem,status,generation,replay,review in PACKAGES:
        d=read(stem);assert d['status']==status,(stem,d['status']);docs[stem]=d
        packages.append(dict(package=stem,status=status,direct_bindings=len(d['sources']),
            generation=dict(identifier=generation,exit_code=0,status='PASS'),
            fresh_replay=dict(identifier=replay,exit_code=0,status='PASS'),review_scope=review))
    for stem,d in docs.items():
        assert d['sources']
        for p,h in d['sources'].items():
            assert digest(p)==h,('binding',stem,p)
            assert p not in bound or bound[p]==h;bound[p]=h
        bindings+=len(d['sources'])
    return docs,bound,packages,bindings


def artifacts():
    return list(dict.fromkeys([p for stem,*_ in PACKAGES for p in
        (f'scripts/{stem}.py',f'results/{stem}.json',f'results/{stem}.md')]+
        HELPERS+list(FROZEN_PREDECESSORS)+[CHECKPOINT,str(Path(__file__).relative_to(ROOT))]))


def build():
    docs,bound,packages,bindings=collect();values=headlines(docs)
    paths=artifacts();asts=[p for p in paths if p.endswith('.py')]
    for p in asts:ast.parse((ROOT/p).read_text(),filename=p)
    reports=[f'results/{s}.md' for s,*_ in PACKAGES]+[p for p in FROZEN_PREDECESSORS if p.endswith('.md')]+[CHECKPOINT]+MUTABLE
    links=[]
    for p in reports:
        for raw in re.findall(r'\]\(([^\n]+?)\)',(ROOT/p).read_text()):
            target=raw.split(' "',1)[0].strip('<>')
            if target.startswith(('http:','https:','mailto:','#')):continue
            links_exist(p,target);links.append(dict(document=p,target=target))
    return dict(status=STATUS,date='2026-10-18',checkpoint=CHECKPOINT,saved_headlines=values,
        packages=packages,predecessor=PREDECESSOR,
        counts=dict(completed_packages=len(packages),direct_bindings=bindings,
            predecessor_direct_bindings=len(docs[PREDECESSOR]['sources']),distinct_bound_files=len(bound),
            python_asts=len(asts),documents=len(reports),local_links=len(links)),
        artifact_sha256={p:digest(p) for p in paths},bound_file_sha256=bound,
        mutable_document_snapshots={p:digest(p) for p in MUTABLE},python_ast_files=asts,local_links=links,
        scope='Integrity, exact saved headline/status comparisons, syntax and links only. Prior run identifiers are recorded evidence, not reruns by this auditor. No mathematical, physical, numerical, optimization, classification or proof-assistant replay. Mutable document hashes are historical snapshots.')


def verify(data):
    assert data['status']==STATUS and data['checkpoint']==CHECKPOINT
    assert set(data['artifact_sha256'])==set(artifacts())
    assert data['python_ast_files']==[p for p in artifacts() if p.endswith('.py')]
    for group in ('artifact_sha256','bound_file_sha256'):
        for p,h in data[group].items():assert digest(p)==h,(group,p)
    docs,bound,packages,bindings=collect()
    assert bound==data['bound_file_sha256'] and packages==data['packages'] and headlines(docs)==data['saved_headlines']
    assert bindings==data['counts']['direct_bindings']
    assert len(packages)==data['counts']['completed_packages']==7
    assert len(bound)==data['counts']['distinct_bound_files']
    assert len(data['local_links'])==data['counts']['local_links']
    for p in data['python_ast_files']:ast.parse((ROOT/p).read_text(),filename=p)
    for r in data['local_links']:links_exist(r['document'],r['target'])


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/AUDIT
    if args.verify:result=json.loads(path.read_text());verify(result)
    else:
        assert not path.exists(),'Frozen audit exists; use --verify.'
        result=build();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
