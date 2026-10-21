#!/usr/bin/env python3
"""Integrity and saved headlines only; no mathematical producer/helper replay.

Shared status documents are mutable snapshots, not future freeze guards.
"""
import argparse,ast,json,re
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[1]
TAG='dimension59_three_corner_payments_and_oblique_closure_20261018'
CHECKPOINT=f'results/{TAG}.md';AUDIT=f'results/{TAG}_audit.json'
STATUS='PASS_DIMENSION59_THREE_CORNER_CHECKPOINT_METADATA'
FINAL_READY=True
BASE='fresh_depth6_dimension59_20261018'
D5='fresh_d5_d12_alternate_coupon_rebalance_20261018'
NODE='fresh_d17_existing_three_corner_transfer_20261018'
LEAF='fresh_finite_d17_leaf_a16_half_transfer_20261018'
D16='fresh_finite_d16_alternate_half_zero_transfer_20261018'
A17='fresh_a17_five_pair_alternate_transfer_20261018'
IDEAL='fresh_order15_bd_basic_fractional_transfer_20261018'
COMBINED='fresh_order15_d16_d17_ambient59_frontier_20261018'
JOINT='fresh_raw56_raw57_raw76_joint_marker_family_20261018'
PILOT='fresh_raw14_raw112_raw77_raw108_core_incidence_pilot_20261018'
CON='fresh_moving64_oblique_r_closure_20261018'
PACKAGES=[
 (D5,'PASS_D5_D12_ALTERNATE_COUPON_REBALANCE','4937/3a01df','81754/2b1408'),
 (NODE,'PASS_D17_EXISTING_THREE_CORNER_TRANSFER','58007/1144ae','88973/f6514e'),
 (LEAF,'PASS_FINITE_D17_LEAF_A16_HALF_TRANSFER','49744/a5b164','35615/c0b025'),
 (D16,'PASS_FINITE_D16_ALTERNATE_HALF_ZERO_TRANSFER','99536/544528','94172/4b3425'),
 (A17,'PASS_A17_FIVE_PAIR_ALTERNATE_TRANSFER','10506/877466','44888/3a2934'),
 (IDEAL,'PASS_ORDER15_BD_BASIC_FRACTIONAL_TRANSFER','41241/0c3a92','29179/e9816b'),
 (COMBINED,'PASS_ORDER15_D16_D17_AMBIENT59_FRONTIER','98230/63ceeb','84439/b967dd'),
 (JOINT,'PASS_RAW56_RAW57_RAW76_JOINT_COMPLETE_FAMILY','70276/80687b','98160/68975f'),
 (PILOT,'PASS_BOUNDED_FOUR_RAW_CORE_INCIDENCE_PILOT','55998/97ceed','91924/957b47'),
 (CON,'BOUNDED_MOVING64_OBLIQUE_R_CLOSURE','13660/2ede67','29934/585488')]
HELPERS=['scripts/check_raw56_raw57_raw76_joint_saved_20261018.py',
 'scripts/check_raw_four_core_incidence_saved_20261018.py',
 'scripts/check_moving64_oblique_r_saved_20261018.py']
PINS={
 'scripts/fresh_depth6_dimension59_20261018.py':'c7f67913ef0d64ff6e23409701df0d1b7c64ce4e8098027f0acd4027e806d46c',
 'results/fresh_depth6_dimension59_20261018.json':'2b09dac21f742a216bd10229fe25fba3f1e4114b8f43fae1206dcb1840169bab',
 'results/fresh_depth6_dimension59_20261018.md':'9bc6ca854bc7bceb8e8f5a50915f391c0cec345f901be09e251d09917069cdd8',
 'scripts/check_dimension59_contextual_checkpoint_20261018.py':'d9a6ea4610984d28e15e235f177b687e1f349cee37067dfe123b6e3f79f66332',
 'results/dimension59_contextual_payments_and_moving64_20261018.md':'35204948afe5820bba89c0717794cc18adcefcc441711e0c31099601fcb6c816',
 'results/dimension59_contextual_payments_and_moving64_20261018_audit.json':'85ed131422dd2ba465e0238d687d6c4d7efb36a0d6a44355301b3b6d45b3cf94'}
MUTABLE=['README.md','RESEARCH_PLAN.md','GOAL.md']


def digest(name):
    h=sha256()
    with (ROOT/name).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())


def headlines(d):
    old=d[BASE];assert old['status']=='PASS_COMPLETE_DEPTH6_DIMENSION59_COMPOSITION'
    assert (old['dimension_bound'],old['excluded_dimension'])==(59,60)
    assert old['finite_maximum']=='404' and old['finite_curvature']=='405'
    assert old['actual_ideal_absolute_surplus']=='54' and old['strict_global_constant']=='4500'
    assert old['contradiction']=='0 >= 1*v_f + 54*v_i + 4500 > 0'
    assert not any(old[k] for k in ('arithmeticity_required','compactness_required','simple_cusps_required',
        'facet_count_bound_required','new_finite_volume_example','external_refereeing','formal_proof_assistant_verified'))
    a=d[D5];assert a['finite']['exact_maximum']=='795/2' and len(a['finite']['above_curvature'])==12
    assert a['all_old_broad_halves_preserved'] and a['ideal_weights_unchanged']
    n=d[NODE];assert (n['old_cost'],n['new_cost'],n['component_saving'])==('467/4','927/8','7/8')
    assert len(n['selected_proofs'])==4 and n['all_cross_component_costs_unchanged']
    l=d[LEAF];assert l['zero_pairs']==[[3,15],[3,16]] and l['new_cost']=='463/4'
    assert l['retained_whole_A17_unit_pairs']==[[0,13],[3,16]]
    assert l['new_cross_component_coupons'][0]['component_pair_increment']=='1'
    x=d[D16];assert x['changed_pairs']==[[2,14],[2,15]] and x['new_cost']=='439/4'
    assert x['typed_costs']['new_rank16_representative_upper']=='441/4'
    assert x['combined_literal_interface']['counts']==dict(debits=0,credits=0)
    assert len(x['selected_proofs']['proofs'])==8 and x['selected_proofs']['affected_maps']==16
    z=d[A17];assert (z['old_cost'],z['new_cost'],z['component_saving'])==('467/4','451/4','4')
    assert len(z['selected_proofs'])==5 and [r['actual_component_increment'] for r in z['new_coupons']]==['4','2']
    i=d[IDEAL];assert i['new_family_costs']==dict(A='75',B='91',C='91',D='367/4')
    assert i['new_order15_envelope']=='367/60' and i['new_coupon']['component_pair_increment']=='1/16'
    for s in (D5,NODE,LEAF,D16,A17,IDEAL):
        assert d[s]['new_dimension_bound'] is None
        assert d[s]['triangle_minimum']==2 and d[s]['quadrilateral_minimum']==1
    c=d[COMBINED];f=c['finite'];q=c['ideal'];lit=c['combined_literal_interface']
    assert c['completed_unrestricted_bound']==c['ambient_dimension']==59 and c['new_dimension_bound'] is None
    assert f['exact_maximum']=='793/2' and f['curvature']=='783/2'
    assert f['actual_attainer']==['A16','A16','D27']
    assert (f['actual_D16'],f['representative_rank16'],f['actual_A16'],f['actual_D17'],f['fixed_A17'])==('439/4','441/4','441/4','919/8','451/4')
    assert f['fixed_large_vectors']==len(f['alternatives'])==961
    assert f['fixed_large_and_A5_vectors']==len(f['expanded_alternatives'])==2574
    assert f['interaction_keys']==2572 and len(f['tables'])==5144
    assert f['Bellman_upper_inequalities']==6718064 and f['attaining_predecessors']==303496
    assert len(f['above_curvature'])==9 and all(r['sixteenths_upper']>6264 for r in f['above_curvature'])
    assert max(r['sixteenths_upper'] for r in f['alternatives'])==6344
    assert f['fixed_rank17']=='A17' and f['residual_rank17']=='D17' and not f['geometric_attainment_claimed']
    assert q['complete_profiles']==len(q['rows'])==1890 and len({tuple(r['ranks']) for r in q['rows']})==1890
    assert q['summary']['negative_count']==sum(F(r['normalized_margin'])<0 for r in q['rows'])==6
    assert q['summary']['zero_count']==sum(F(r['normalized_margin'])==0 for r in q['rows'])==0
    assert q['summary']['normalized_minimum']=='-5241/57460'
    assert F(next(r['normalized_margin'] for r in q['rows'] if r['ranks']==[8,9,9,9,9,14]))==F(8,1125)
    assert not q['fresh_convolution_performed'] and q['negative_absolute_margin_not_extended_to_arbitrary_profiles']
    assert lit['configurations']==len(lit['rows'])==216 and lit['minimum_total']=='1'
    assert lit['changes']==dict(A1_A16_credits=8,D17_debits=16)
    assert lit['bound_local_corner_delta_comparisons']==4160 and lit['local_keysets_and_marks_identical']
    assert c['strict_global_constant']=='4350' and c['constant_counted_once']
    ar=d[JOINT];ac=ar['counts'];e=ar['exported'];r=e['repair']
    assert ac==dict(cross_marker_positions=30336,exported_columns=1430,family_classes=1643,marker_orientations=24788,
        overlap_targets=53,raised_allowances=2,retained_cross_marker_pairs=9888,strictly_violated_columns=3,targets=1394)
    assert len(ar['forest_records'])==1394 and len(e['columns'])==1430
    assert len(e['preserved_multi_owner_costs'])==14 and not ar['current_graphs_changed'] and ar['optimizer_invocations']==0
    assert F(r['bound'])==F(14671331978990040,1220681367066101)>12
    p=d[PILOT];assert p['counts']==dict(core_support_tests=92,exact_compact_simplex_faces=6,exact_native_pilot_maps=40,
        extra_marker_domains=784,extra_role_orientations=12688,family_classes=1331)
    assert [(r['raw_family'],r['actual_extra_vertices']) for r in p['summaries']]==[(14,0),(112,2),(77,0),(108,2)]
    assert p['new_pair_tests']==p['optimizer_invocations']==0 and p['exported_columns']==[] and not p['old_outcomes_changed']
    assert p['current_joint_bound']==r['bound']
    w=p['future_pair_workload'];assert (w['new_marker_pairs'],w['new_old_marker_pairs'],w['new_ordinary_pairs'],w['target_forests'])==(949752,486316,1207540,673)
    j=d[CON];s=j['summary'];assert s['reconstructed_M_prefixes']==63 and s['R_recipes']==647 and s['R_branches']==1294
    assert s['native_R']==s['completed_disks']==s['fewer_than_ten_ends']==s['unresolved']==0
    assert (s['reverse_prefixes'],s['reverse_M_recipes'],s['reverse_M_branches'],s['reverse_native_M'],s['reverse_clean_M'])==(11,74,148,52,52)
    assert s['reverse_completed_disks']==s['reverse_fewer_than_ten_ends']==0 and s['old8_duplicates']==63
    assert len(j['reverse_accepted'])==52 and all(v['clean'] and v['whole_J_closed'] and v['finite_rays']==21 and v['open_ends']==13 for v in j['reverse_accepted'])
    assert j['exact_deletion_debit']['lost_finite_count']==19 and j['exact_deletion_debit']['compact_faces_by_dimension']=={'2':12,'3':1}
    return dict(completed_bound=59,no_identified_gap_in_completed_proof=True,new_dimension_bound=None,
        finite_rank59_maximum='793/2',finite_curvature='783/2',finite_representative_obstructions=9,
        ideal_terminal_profiles=1890,ideal_negative_profiles=6,arithmetic_fixed_mass_ratio=r['bound'],
        arithmetic_new_optimizer_invocations=0,new_finite_volume_example=False,external_refereeing=False,
        proof_assistant_formalization=False,clean_numerical_reverse_caps=52,all_reverse_open_ends=13)


def collect():
    assert FINAL_READY and all('PENDING' not in str(r) for r in PACKAGES),'Wait for final frozen inputs.'
    for p,h in PINS.items():assert digest(p)==h,('predecessor',p)
    docs={BASE:read(BASE)};records=[];bound={};total=0
    for stem,status,generation,replay in PACKAGES:
        d=read(stem);assert d['status']==status;docs[stem]=d
        records.append(dict(stem=stem,status=status,generation=generation,fresh_replay=replay,
            run_scope='Previously recorded terminal0/PASS evidence; not executed by this auditor.'))
    for stem,d in docs.items():
        assert d['sources']
        for p,h in d['sources'].items():
            if p not in bound:bound[p]=digest(p)
            assert bound[p]==h,(stem,p)
        total+=len(d['sources'])
    return docs,records,bound,total


def artifacts():
    return list(dict.fromkeys([p for s,*_ in PACKAGES for p in
        (f'scripts/{s}.py',f'results/{s}.json',f'results/{s}.md')]+HELPERS+list(PINS)+
        [CHECKPOINT,str(Path(__file__).relative_to(ROOT))]))


def document_links(document):
    out=[]
    for raw in re.findall(r'\]\(([^\n]+?)\)',(ROOT/document).read_text()):
        target=raw.split(' "',1)[0].strip('<>')
        if target.startswith(('http:','https:','mailto:','#')):continue
        out.append(dict(document=document,target=target))
    return out


def check_link(row):
    p=((ROOT/row['document']).parent/unquote(row['target'].split('#',1)[0])).resolve()
    assert p.is_file() or p==(ROOT/AUDIT).resolve(),row


def build():
    docs,records,bound,total=collect();values=headlines(docs);paths=artifacts()
    asts=[p for p in paths if p.endswith('.py')]
    for p in asts:ast.parse((ROOT/p).read_text(),filename=p)
    reports=[p for p in paths if p.endswith('.md')]+MUTABLE
    links=[r for p in reports for r in document_links(p)]
    for r in links:check_link(r)
    snapshots={}
    for p in MUTABLE:
        opening=(ROOT/p).read_text()[:6000];assert CHECKPOINT in opening,p
        assert '59' in opening and 'finite-volume' in opening
        snapshots[p]=dict(sha256=digest(p),opening=opening)
    return dict(status=STATUS,checkpoint=CHECKPOINT,packages=records,saved_headlines=values,
        counts=dict(new_packages=len(PACKAGES),direct_bindings=total,distinct_bound_files=len(bound),python_asts=len(asts),
            documents=len(reports),local_links=len(links)),artifact_sha256={p:digest(p) for p in paths},
        bound_file_sha256=bound,python_ast_files=asts,documents=reports,local_links=links,mutable_snapshots=snapshots,
        scope='Metadata and saved headline/cardinality checks only. No producer/helper execution, geometry, census, native arithmetic, finite Bellman proof, convolution, cone or optimizer replay. Mutable document hashes/openings are historical snapshots.')


def verify(saved):
    assert saved['status']==STATUS and saved['checkpoint']==CHECKPOINT
    docs,records,bound,total=collect();paths=artifacts()
    assert saved['packages']==records and saved['saved_headlines']==headlines(docs)
    assert set(saved['artifact_sha256'])==set(paths) and saved['bound_file_sha256']==bound
    for p,h in saved['artifact_sha256'].items():assert digest(p)==h,p
    asts=[p for p in paths if p.endswith('.py')];reports=[p for p in paths if p.endswith('.md')]+MUTABLE
    assert saved['python_ast_files']==asts and saved['documents']==reports
    assert saved['counts']==dict(new_packages=len(PACKAGES),direct_bindings=total,distinct_bound_files=len(bound),
        python_asts=len(asts),documents=len(reports),local_links=len(saved['local_links']))
    assert set(saved['mutable_snapshots'])==set(MUTABLE)
    # Replay frozen-document link completeness. Mutable openings and links
    # remain historical snapshots, so later authorized status edits are allowed.
    frozen_links=[r for p in reports if p not in MUTABLE for r in document_links(p)]
    assert [r for r in saved['local_links'] if r['document'] not in MUTABLE]==frozen_links
    for p in asts:ast.parse((ROOT/p).read_text(),filename=p)
    for r in saved['local_links']:check_link(r)
    for p,v in saved['mutable_snapshots'].items():assert CHECKPOINT in v['opening'] and len(v['sha256'])==64


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();p=ROOT/AUDIT
    if args.verify:out=json.loads(p.read_text());verify(out)
    else:assert not p.exists();out=build();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
