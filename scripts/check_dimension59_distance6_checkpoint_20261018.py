#!/usr/bin/env python3
"""Saved-headline and integrity checks only; no mathematical replay.

Mutable status files are recorded as historical snapshots. All theorem
packets and this checkpoint are immutable after authoritative generation.
"""
import argparse,ast,json,re
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[1]
TAG='dimension59_distance6_and_contextual_ideal_frontier_20261018'
CHECKPOINT=f'results/{TAG}.md';AUDIT=f'results/{TAG}_audit.json'
STATUS='PASS_DIMENSION59_DISTANCE6_CHECKPOINT_METADATA'
BASE='fresh_depth6_dimension59_20261018'
FINITE='fresh_depth6_finite6_zero_transfer_20261018'
I13='fresh_order13_bc_inner_terminal_transfer_20261018'
I17='fresh_order17_center_contextual_zero_transfer_20261018'
COMBINED='fresh_depth6_finite6_order13_order17_ambient59_frontier_20261018'
ARITH='fresh_raw108_complete_marker_family_20261018'
CON='fresh_oblique52_boundary_connector_20261018'
REVERSE='fresh_oblique52_reverse_r_closure_20261018'
PACKAGES=[
 (FINITE,'PASS_DEPTH6_FINITE6_ZERO_LOCAL_TRANSFER','22306/9f27a2','29495/ea8eac'),
 (I13,'PASS_ORDER13_B143_OVER2_C71_INNER_TERMINAL_TRANSFER','20307/0648e5','30080/b4bb29'),
 (I17,'PASS_ORDER17_CENTER_CONTEXTUAL_ZERO_TRANSFER','6908/47e0ca','73245/2c938b'),
 (COMBINED,'PASS_DEEP_FINITE6_ORDER13_ORDER17_AMBIENT59_FRONTIER','80994/1f585d','21070/56c5d1'),
 (ARITH,'PASS_RAW108_COMPLETE_EXTRA_VERTEX_MARKER_FAMILY','56256/dba596','56955/bca986'),
 (CON,'BOUNDED_OBLIQUE52_BOUNDARY_CONNECTOR','5003/99095d','64410/fb9504'),
 (REVERSE,'BOUNDED_OBLIQUE52_REVERSE_R','19327/49248a','41091/38daaa')]
HELPERS=['scripts/check_raw108_complete_marker_saved_20261018.py',
 'scripts/check_oblique52_boundary_connector_saved_20261018.py',
 'scripts/check_oblique52_reverse_r_saved_20261018.py',
 'results/fresh_oblique52_reverse_r_closure_20261018_pilot.json']
PINS={
 'scripts/fresh_depth6_dimension59_20261018.py':'c7f67913ef0d64ff6e23409701df0d1b7c64ce4e8098027f0acd4027e806d46c',
 'results/fresh_depth6_dimension59_20261018.json':'2b09dac21f742a216bd10229fe25fba3f1e4114b8f43fae1206dcb1840169bab',
 'results/fresh_depth6_dimension59_20261018.md':'9bc6ca854bc7bceb8e8f5a50915f391c0cec345f901be09e251d09917069cdd8',
 'scripts/check_dimension59_three_corner_checkpoint_20261018.py':'1d240b5e430677e8dab2b31f2182e05253af1b39f0bb337e7af2618f75569380',
 'results/dimension59_three_corner_payments_and_oblique_closure_20261018.md':'1f04c0c351930f7e14316703e59baefd59889e220bc3c53af7c7590f594e0634',
 'results/dimension59_three_corner_payments_and_oblique_closure_20261018_audit.json':'73da59059e1a740ca8e6ac284231c72b5a14ed78d78513c3cc11223e277bf327'}
MUTABLE=['README.md','RESEARCH_PLAN.md','GOAL.md']


def digest(name):
    h=sha256()
    with (ROOT/name).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())


def headlines(d):
    b=d[BASE]
    assert b['status']=='PASS_COMPLETE_DEPTH6_DIMENSION59_COMPOSITION'
    assert (b['dimension_bound'],b['excluded_dimension'])==(59,60)
    assert (b['finite_maximum'],b['finite_curvature'],b['actual_ideal_absolute_surplus'],b['strict_global_constant'])==('404','405','54','4500')
    assert b['contradiction']=='0 >= 1*v_f + 54*v_i + 4500 > 0'
    assert not any(b[k] for k in ('arithmeticity_required','compactness_required','simple_cusps_required',
        'facet_count_bound_required','new_finite_volume_example','external_refereeing','formal_proof_assistant_verified'))
    f=d[FINITE]
    assert f['feature_depth']==f['finite_distance']==6 and f['new_weight']=='0'
    assert (f['profiles']['candidate_count'],f['profiles']['accepted_count'])==(163,19)
    assert (f['adaptive_schur']['branch_count'],f['adaptive_schur']['strict_branches'])==(829,822)
    assert len(f['exceptional_payment_proofs']['branch_reduction'])==7 and len(f['exceptional_payment_proofs']['proofs'])==4
    assert f['new_coupon']['component_pair_increment']=='2' and f['new_coupon']['unchanged_representative_allowance']=='5/2'
    assert f['current_literal_interface']['changed']==[] and f['simultaneous_actual_type_pairs']==32400
    assert f['all_ideal_weights_unchanged'] and f['all_representative_cross_costs_unchanged']
    h=d[I13]
    assert h['new_order13_maximum']=='143/2' and h['new_coupon']['component_pair_increment']=='1/2'
    assert h['current_literal_interface']['changed']==[]
    z=d[I17]
    assert z['local_proof']['new_cost']=='227/2' and z['new_coupon']['component_pair_increment']=='1'
    assert z['new_coupon']['required_whole_unmarked_components']==['A1','A9'] and z['new_coupon']['witnesses_not_multiplied']
    assert z['current_literal_interface']['changes']==dict(ideal_debits=2,context_credits=2)
    assert z['finite']['exact_maximum']=='793/2' and z['pending_finite6_rule_assumed'] is False
    c=d[COMBINED];cf=c['finite'];ci=c['ideal'];lit=c['combined_literal_interface']
    assert c['ambient_dimension']==c['completed_unrestricted_bound']==59 and c['new_dimension_bound'] is None
    assert c['geometric_inputs_final'] and c['strict_global_constant']=='4350' and c['constant_counted_once']
    assert cf['exact_maximum']=='394' and cf['curvature']=='783/2' and cf['actual_attainer']==['A14','A15','A15','A15']
    assert cf['fixed_large_vectors']==len(cf['alternatives'])==961
    assert cf['fixed_large_and_A5_vectors']==len(cf['expanded_alternatives'])==2574
    assert cf['interaction_keys']==2572 and len(cf['tables'])==5144
    assert cf['Bellman_upper_inequalities']==6718064 and cf['attaining_predecessors']==303496 and cf['residual_item_checks']==128700
    assert cf['simultaneous_type_pairs']==32400 and cf['new_residual_tables'] and not cf['geometric_attainment_claimed']
    assert len(cf['above_curvature'])==3 and sorted(r['sixteenths_upper'] for r in cf['above_curvature'])==[6272,6292,6304]
    assert max(r['sixteenths_upper'] for r in cf['alternatives'])==6304
    assert ci['complete_profiles']==len(ci['rows'])==len({tuple(r['ranks']) for r in ci['rows']})==1890
    assert ci['summary']['negative_count']==sum(F(r['normalized_margin'])<0 for r in ci['rows'])==1
    assert ci['summary']['zero_count']==sum(F(r['normalized_margin'])==0 for r in ci['rows'])==0
    assert ci['summary']['normalized_minimum']=='-69/5440' and ci['summary']['negative_rows'][0]['ranks']==[9,9,9,15,16]
    assert ci['next_order16_threshold']['strict_total_cost_saving_threshold']=='69/340'
    assert ci['no_new_convolution'] and ci['negative_absolute_margins_not_extended']
    assert lit['configurations']==len(lit['rows'])==216 and lit['minimum_total']=='1'
    a=d[ARITH];ac=a['counts'];e=a['exported'];repair=e['repair']
    assert (ac['classes'],ac['pointed_marker_domains'],ac['target_forests'],ac['markers'],ac['relative_cases'])==(196,392,383,1676,32)
    assert (ac['exported_columns'],ac['strictly_violated_columns'],ac['raised_allowances'])==(415,8,8)
    assert len(e['preserved_columns'])==2884 and len(e['preserved_multi_owner_costs'])==14
    assert F(repair['bound'])==F(14671331978990040,1221031586432741)>12
    assert a['optimizer_invocations']==0 and not a['old_outcomes_changed']
    con=d[CON];ss=con['summary']
    assert ss==dict(frozen_pointed_configurations=52,distinct_recipe_prefixes=52,tested_prefixes=52,
        recipes=1271,branches=2542,native_S=449,clean_S=449,completed_disks=0,fewer_than_ten_ends=0,unresolved=0)
    assert len(con['accepted'])==449 and all(r['finite_rays']==23 and r['open_ends']==13 for r in con['accepted'])
    assert con['exact_deletion_debit']['lost_finite_count']==19 and con['exact_deletion_debit']['compact_faces_by_dimension']=={'2':12,'3':1}
    rev=d[REVERSE]
    assert rev['summary']==dict(frozen_pointed_configurations=449,distinct_recipe_prefixes=449,tested_prefixes=449,
        recipes=16894,branches=33788,native_R=376,clean_R=376,completed_disks=0,fewer_than_ten_ends=0,unresolved=0)
    assert len(rev['accepted'])==376 and all(r['finite_rays']==22 and r['open_ends']==14 for r in rev['accepted'])
    assert rev['exact_deletion_debit']==con['exact_deletion_debit'] and rev['topology']==con['topology']
    for s in (FINITE,I13,I17,COMBINED):
        assert d[s]['new_dimension_bound'] is None and d[s]['triangle_minimum']==2 and d[s]['quadrilateral_minimum']==1
    return dict(completed_bound=59,no_identified_gap_in_completed_proof=True,new_dimension_bound=None,
        finite_rank59_maximum='394',finite_curvature='783/2',finite_representative_obstructions=3,
        ideal_terminal_profiles=1890,ideal_negative_profiles=1,ideal_negative_margin='-69/5440',
        arithmetic_fixed_mass_ratio=repair['bound'],arithmetic_optimizer_invocations=0,
        numerical_S_connectors=449,all_S_open_ends=13,numerical_reverse_R_walls=376,all_reverse_R_open_ends=14,new_finite_volume_example=False,
        external_refereeing=False,proof_assistant_formalization=False)


def collect():
    assert all('PENDING' not in str(p) for p in PACKAGES),'Wait for final frozen packets.'
    for p,h in PINS.items():assert digest(p)==h,('predecessor',p)
    docs={BASE:read(BASE)};records=[];bound={};total=0
    for stem,status,generation,replay in PACKAGES:
        d=read(stem);assert d['status']==status;docs[stem]=d
        records.append(dict(stem=stem,status=status,generation=generation,fresh_replay=replay,
            scope='Recorded terminal PASS evidence; no execution by this auditor.'))
    for stem,d in docs.items():
        for p,h in d['sources'].items():
            if p not in bound:bound[p]=digest(p)
            assert bound[p]==h,(stem,p)
        total+=len(d['sources'])
    return docs,records,bound,total


def artifacts():
    return list(dict.fromkeys([p for s,*_ in PACKAGES for p in
        (f'scripts/{s}.py',f'results/{s}.json',f'results/{s}.md')]+HELPERS+list(PINS)+
        [CHECKPOINT,str(Path(__file__).relative_to(ROOT))]))


def links(document):
    out=[]
    for raw in re.findall(r'\]\(([^\n]+?)\)',(ROOT/document).read_text()):
        target=raw.split(' "',1)[0].strip('<>')
        if target.startswith(('http:','https:','mailto:','#')):continue
        out.append(dict(document=document,target=target))
    return out


def check_link(r):
    p=((ROOT/r['document']).parent/unquote(r['target'].split('#',1)[0])).resolve()
    assert p.is_file() or p==(ROOT/AUDIT).resolve(),r


def build():
    docs,records,bound,total=collect();values=headlines(docs);paths=artifacts()
    asts=[p for p in paths if p.endswith('.py')]
    for p in asts:ast.parse((ROOT/p).read_text(),filename=p)
    reports=[p for p in paths if p.endswith('.md')]+MUTABLE;ls=[r for p in reports for r in links(p)]
    for r in ls:check_link(r)
    snapshots={}
    for p in MUTABLE:
        opening=(ROOT/p).read_text()[:6000];assert CHECKPOINT in opening and 'finite-volume' in opening and '59' in opening
        snapshots[p]=dict(sha256=digest(p),opening=opening)
    return dict(status=STATUS,checkpoint=CHECKPOINT,packages=records,saved_headlines=values,
        counts=dict(new_packages=len(PACKAGES),direct_bindings=total,distinct_bound_files=len(bound),python_asts=len(asts),documents=len(reports),local_links=len(ls)),
        artifact_sha256={p:digest(p) for p in paths},bound_file_sha256=bound,python_ast_files=asts,documents=reports,
        local_links=ls,mutable_snapshots=snapshots,
        scope='Integrity and saved headline/cardinality checks only. No mathematical producer, helper, geometry, census, Bellman proof, convolution, cone or optimizer replay. Mutable openings are historical snapshots.')


def verify(saved):
    docs,records,bound,total=collect();paths=artifacts()
    assert saved['status']==STATUS and saved['checkpoint']==CHECKPOINT and saved['packages']==records and saved['saved_headlines']==headlines(docs)
    assert set(saved['artifact_sha256'])==set(paths) and saved['bound_file_sha256']==bound
    for p,h in saved['artifact_sha256'].items():assert digest(p)==h,p
    asts=[p for p in paths if p.endswith('.py')];reports=[p for p in paths if p.endswith('.md')]+MUTABLE
    assert saved['python_ast_files']==asts and saved['documents']==reports
    assert saved['counts']==dict(new_packages=len(PACKAGES),direct_bindings=total,distinct_bound_files=len(bound),python_asts=len(asts),documents=len(reports),local_links=len(saved['local_links']))
    assert [r for r in saved['local_links'] if r['document'] not in MUTABLE]==[r for p in reports if p not in MUTABLE for r in links(p)]
    for p in asts:ast.parse((ROOT/p).read_text(),filename=p)
    for r in saved['local_links']:check_link(r)
    assert set(saved['mutable_snapshots'])==set(MUTABLE)
    for p,r in saved['mutable_snapshots'].items():assert CHECKPOINT in r['opening'] and len(r['sha256'])==64


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/AUDIT
    if a.verify:out=json.loads(p.read_text());verify(out)
    else:assert not p.exists();out=build();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
