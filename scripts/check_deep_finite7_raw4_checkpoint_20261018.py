#!/usr/bin/env python3
"""Integrity, saved headline/status, syntax and local-link audit only.

No producer/helper is imported or executed. Recorded run identifiers are prior
evidence. Mutable document hashes are snapshots, not future equality checks.
"""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/deep_finite7_and_raw4_cover_frontier_20261018.md'
AUDIT=CHECKPOINT.replace('.md','_audit.json')
STATUS='PASS_DEEP_FINITE7_RAW4_CHECKPOINT_INTEGRITY_SAVED_STATUS_SYNTAX_LINKS'
PREDECESSOR='fresh_order11_order16_dimension60_20261018'
PACKAGES=[
 ('fresh_depth11_finite7_zero_transfer_20261018','PASS_DEPTH11_FINITE7_ZERO_LOCAL_TRANSFER','2d1eb8','a2aecd',
  'Two internal full mathematical source/report reviews, including inherited minimal-core, current short6 and triangle interfaces. No reviewer matrix, triple, DP or classification replay.'),
 ('fresh_depth10_finite7_zero_transfer_20261018','PASS_DEPTH10_FINITE7_ZERO_LOCAL_TRANSFER','2c69dd','039503',
  'Two internal full source/report mathematical delta reviews: equality radical, genuine ninth white, complete contact choices and retained interfaces. No reviewer matrix, scalar, census or DP replay.'),
 ('fresh_ambient60_payment_frontier_20261018','PASS_CURRENT_AMBIENT60_PAYMENT_FRONTIER','5b0d55','bc9b8f',
  'Root and separate full source/report mathematical reviews of complete finite/ideal domains, typed domination, Bellman certificates, exact tails and conditional scope. No reviewer DP, convolution, scalar or census replay.'),
 ('fresh_new5121_raw4_two_face_local_covers_20261018','PASS_NEW5121_RAW4_TWO_FACE_EXACT_LOCAL_COVERS','643674','eed6be',
  'Root and separate full producer/report/helper mathematical reviews of provenance, native-orbit transport, fixed old graph, exact gates, complete cones and monotone repair. Saved helper0677a8 is separate arithmetic evidence; no source-review physical, graph or LP replay.'),
 ('fresh_new5121_raw4_next_four_domain_20261018','PASS_NEW5121_RAW4_NEXT_FOUR_EXACT_DOMAIN_AND_BENEFIT_CEILING','ee22d9','a8d0e8',
  'Peer full producer/report/helper and inherited-interface mathematical review of domain transport and singleton-cover ceiling. Saved helper58d280 is separate accounting evidence. Frozen report claims no parent review; no reviewer physical, graph, LP or arithmetic replay.'),
 ('fresh_joint30_33_compact10_first_20261018','BOUNDED_JOINT30_33_COMPACT10_FIRST','c5bbf0','72d297',
  'Peer full producer/report/helper scope review: actual-side profiles, complete bounded recipes/joins, duplicateold8, provisional versus final cones and exact old deletion debit. Saved helper73ba06 is separate evidence; no reviewer normal, classifier or numerical cone replay.'),
]
HELPERS=[
 'scripts/check_new5121_raw4_local_covers_saved_20261018.py',
 'scripts/check_new5121_raw4_next_four_saved_20261018.py',
 'scripts/check_joint30_33_compact10_saved_20261018.py',
]
FROZEN_PREDECESSORS = {'scripts/fresh_order11_order16_dimension60_20261018.py': '89d84a0f82d21fdd7c73e57fceadcc5d9859417d85be562b57617e6ef659c671',
 'results/fresh_order11_order16_dimension60_20261018.json': '019c8d9bca5dad4de0e4035acfae154f05b6f96836745370e660eaecc9add6a9',
 'results/fresh_order11_order16_dimension60_20261018.md': '9434e6c00eeb5bb52023aaf5c6411cffecd1d711bad1d99044f1b1a08c8ddeeb',
 'scripts/check_dimension60_checkpoint_20261018.py': '22f6af587401c701b3f57941efbc343e0e0f9bfa7e6b19be400281665e733c38',
 'results/dimension60_complete_and_current_search_20261018.md': '6ae088d6d2d225b58856e2eefe035fd65b1029eb5b59cafce272480f47ab3131',
 'results/dimension60_complete_and_current_search_20261018_audit.json': '626764609a60cdaeb18ec1b4b765ce5f3153af14cff62bbac8beda280a1ed1f9'}
MUTABLE=['README.md','RESEARCH_PLAN.md']
CHECKPOINT_REVIEW='A separate agent completed a full consistency read of the checkpoint, current README/RESEARCH_PLAN sections, all six packet reports and completed dimension60 report: PASS. No computation or mathematical certificate replay.'


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
    assert not any(p[k] for k in ('external_refereeing','formal_proof_assistant_verified','new_finite_volume_example'))
    a,b,f,c,n,j=[docs[x[0]] for x in PACKAGES]
    for d,depth,order in ((a,11,26),(b,10,24)):
        assert (d['finite_distance'],d['feature_depth'],d['new_weight'])==(7,depth,'0')
        assert d['corridor']['path_order']==order and d['triangle_minimum']==2 and d['quadrilateral_minimum']==1
        assert d['ideal_weights_unchanged'] and d['all_other_finite_weights_and_coupons_unchanged']
        assert d['new_dimension_bound'] is None and not d['new_finite_volume_example']
    assert a['triangle_interface']['complete_high_core_marked_triples']==1025
    assert len(a['triangle_interface']['retained_distance7_elementary_payments'])==b['inherited_triangle_protected_pairs']==14
    assert b['corridor']['equality_principal_inertia']==[24,1,1]
    assert b['corridor']['required_contact_difference']=='1/5' and len(b['corridor']['complete_contact_pairings'])==4
    assert b['rank60_finite_maximum_by_monotonicity_and_attainment']=='413'
    assert b['unchanged_rank60_attainer']==['D17','D17','D26']
    assert (f['inherited_completed_dimension_bound'],f['ambient_dimension'],f['link_dimension'])==(60,60,59)
    assert f['finite']['exact_maximum']=='413' and f['finite_curvature']=='405'
    assert (f['finite']['complete_large_vectors'],f['finite']['interaction_keys'])==(1044,792)
    assert f['finite']['final_A1_B10_unit_counted'] and not f['finite']['replay_optimizing_recurrence']
    assert len(f['finite']['representative_alternatives_above_curvature'])==len(f['finite']['actual_witnesses_above_curvature'])==17
    s=f['ideal']['summary'];assert f['ideal']['complete_profiles']==s['count']==2113
    assert (s['negative_count'],s['zero_count'],s['normalized_minimum'],s['absolute_minimum'])==(77,0,'-1823/11560','-192250')
    assert not f['ideal']['negative_terminal_absolute_margin_is_universal_bound'] and f['new_dimension_bound'] is None
    assert f['strict_global_constant']=='4500'
    assert (c['represented_classes'],c['positive_represented_classes'],c['represented_forests'])==(324,261,600)
    assert c['counts']==dict(new_actual_corners=12,new_exclusions=1332,new_literal_stars=28,new_pair_tests=1364,new_unknown_pairs=32)
    assert (c['new_charge_columns'],c['inherited_forest_count'],len(c['preserved_multi_owner_columns']))==(4,594,14)
    assert (c['complete_count_vectors'],c['complete_maximal_cliques'])==(89,1258)
    assert c['bound']=='205569148850452800/16961245013123269' and c['global_optimizer_invocations']==c['replay_local_LPs']==0
    assert n['counts']==dict(active_embeddings=204,actual_corner_occurrences=24,distinct_native_stars=102,faces=4,literal_orientation_occurrences=152,vertex_types=12)
    assert (n['baseline_represented_classes'],n['baseline_represented_forests'])==(324,600)
    v=n['cover_inventory'];assert (v['affected_existing_forests'],v['potential_new_forests'],v['candidate_new_pair_positions'])==(11,1,18147)
    assert v['slack_aware_ratio_lower_bound']=='205569148850452800/17073378168395869'
    assert n['pair_compatibility_tests']==0 and not any(n[k] for k in ('new_global_optimization','new_local_cover','current_compatibility_graphs_changed'))
    assert n['new_dimension_bound'] is None
    assert j['counts']==dict(norm_equations=85,norm_branches=170)
    assert j['exact_deletion_debit']['lost_finite_count']==15 and j['exact_deletion_debit']['compact_faces_by_dimension']=={'2':10,'3':1}
    assert [x['complete_low_recipe_count'] for x in j['banks']]==[64,21]
    assert j['complete_pair_count']==len(j['pairs'])==12 and j['pair_histogram']=={'None':9,'dotted':2,'2':1}
    assert len(j['coherent_pairs'])==3 and len(j['unresolved'])==13
    restored=[x for x in j['coherent_pairs'] if x['old30_restored'] and x['old33_restored']]
    assert len(restored)==1 and restored[0]['open_ends']==10
    assert all(len(x['whole_face_records'])==4 for x in j['coherent_pairs'])
    assert all(not x['compact10_preserved'] and x['open_ends']==11 for x in j['coherent_pairs'] if not x['old33_restored'])
    return dict(completed_unrestricted_dimension_bound=60,new_dimension_bound=None,
        current_finite_rank60_maximum='413',finite_curvature='405',ideal_terminal_profiles=2113,negative_ideal_profiles=77,
        current_arithmetic_classes=324,current_arithmetic_forests=600,current_arithmetic_ratio=c['bound'],
        next_four_pair_only_ratio_lower_bound=v['slack_aware_ratio_lower_bound'],next_four_pairs_tested=0,
        new_finite_volume_example=False,external_refereeing=False,formal_proof_assistant_verified=False)


def collect():
    for p,h in FROZEN_PREDECESSORS.items():assert digest(p)==h,('predecessor',p)
    assert len(FROZEN_PREDECESSORS)==6
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
        packages=packages,predecessor=PREDECESSOR,checkpoint_consistency_review=CHECKPOINT_REVIEW,
        counts=dict(completed_packages=len(packages),direct_bindings=bindings,
            predecessor_direct_bindings=len(docs[PREDECESSOR]['sources']),distinct_bound_files=len(bound),
            python_asts=len(asts),documents=len(reports),local_links=len(links)),
        artifact_sha256={p:digest(p) for p in paths},bound_file_sha256=bound,
        mutable_document_snapshots={p:digest(p) for p in MUTABLE},python_ast_files=asts,local_links=links,
        scope='Integrity, exact saved headline/status comparisons, syntax and links only. Prior run identifiers are recorded evidence, not reruns by this auditor. No mathematical, physical, numerical, optimization, classification or proof-assistant replay. Mutable document hashes are historical snapshots.')


def verify(data):
    assert data['status']==STATUS and data['checkpoint']==CHECKPOINT
    assert set(data['artifact_sha256'])==set(artifacts())
    assert data['checkpoint_consistency_review']==CHECKPOINT_REVIEW
    assert data['python_ast_files']==[p for p in artifacts() if p.endswith('.py')]
    for group in ('artifact_sha256','bound_file_sha256'):
        for p,h in data[group].items():assert digest(p)==h,(group,p)
    docs,bound,packages,bindings=collect()
    assert bound==data['bound_file_sha256'] and packages==data['packages'] and headlines(docs)==data['saved_headlines']
    assert bindings==data['counts']['direct_bindings']
    assert len(packages)==data['counts']['completed_packages']==6
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
