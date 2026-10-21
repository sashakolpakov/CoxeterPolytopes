#!/usr/bin/env python3
"""Checkpoint integrity and saved status only; no producer/helper execution.

Run identifiers and reviews are recorded prior evidence, not replayed here.
Mutable README/plan hashes are historical snapshots, not future freeze guards.
"""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/finite411_and_fractional_d17_frontier_20261018.md'
AUDIT = CHECKPOINT.replace('.md', '_audit.json')
STATUS = 'PASS_FINITE411_FRACTIONAL_D17_CHECKPOINT_INTEGRITY_SAVED_STATUS_SYNTAX_LINKS'
PREDECESSOR = 'fresh_order11_order16_dimension60_20261018'
FINITE = 'fresh_finite_d17_d5_d12_unit_transfer_20261018'
IDEAL = 'fresh_order17_d5_d12_unit_coupon_zero_transfer_20261018'
QUARTER = 'fresh_order17_d5_d12_center_quarter_transfer_20261018'
COMBINED = 'fresh_d17_fractional_ambient60_frontier_20261018'
MODEL = 'fresh_current324_full_model_20261018'
CONSTRUCTION = 'fresh_joint10_30_65_compact_pair_20261018'
PACKAGES = [
    (FINITE, 'PASS_FINITE_D17_D5_D12_UNIT_TRANSFER', '895600', '926a28',
     'Two peer full source/report mathematical reads. No reviewer selected-map, inertia, literal, census or finite-budget arithmetic replay.'),
    (IDEAL, 'PASS_ORDER17_D114_WHOLE_D5_D12_UNIT_TRANSFER', '76b709', 'c6fd44',
     'Root and peer full source/report mathematical reads, including arbitrary-contact sealing and current interfaces. No reviewer selected-map, inertia, census or scalar replay.'),
    (QUARTER, 'PASS_ORDER17_D227_OVER2_D5_D12_CENTER_QUARTER_TRANSFER', 'fcf696', '926d0e',
     'Root and peer full source/report mathematical reads of the inherited two whole-affine seals, disjoint center coupon and complete current interfaces. No reviewer map, inertia, census or arithmetic replay.'),
    (COMBINED, 'PASS_FINITE411_FRACTIONAL_IDEAL17_AMBIENT60_FRONTIER', '3a59d5', 'be16ca',
     'Two peer full source/report mathematical composition reads: disjoint actual coupon roles, complete finite upper/attainment certificates, current literal interfaces, low9 exact-tail profile domain and conditional scope. No reviewer DP, scalar, literal, native or census replay.'),
    (MODEL, 'PASS_EXACT_CURRENT324_MODEL_CERTIFICATES', '5ee521', 'cd25ad',
     'Root and peer full producer/report/helper and inherited matrix/certificate-interface mathematical reads; root final report/helper delta bb5831. Helper431c88 is separate saved-data evidence. No reviewer matrix, solver, native, clique or arithmetic replay.'),
    (CONSTRUCTION, 'BOUNDED_JOINT10_30_65_COMPACT_PAIR', '500ad0', '3886e8',
     'Root full producer/report/helper and inherited gate/score-interface mathematical and scope read. Second agent full helper/domain/ledger and associated report-scope read. Helperad6146 is separate saved-data evidence. No reviewer execution, normal reconstruction, physical interval or cone replay.'),
]
HELPERS = [
    'scripts/check_current324_model_saved_20261018.py',
    'scripts/check_joint10_30_65_compact_saved_20261018.py',
]
FROZEN_PREDECESSORS = {
    'scripts/fresh_order11_order16_dimension60_20261018.py': '89d84a0f82d21fdd7c73e57fceadcc5d9859417d85be562b57617e6ef659c671',
    'results/fresh_order11_order16_dimension60_20261018.json': '019c8d9bca5dad4de0e4035acfae154f05b6f96836745370e660eaecc9add6a9',
    'results/fresh_order11_order16_dimension60_20261018.md': '9434e6c00eeb5bb52023aaf5c6411cffecd1d711bad1d99044f1b1a08c8ddeeb',
    'scripts/check_rank15_order15_raw4_checkpoint_20261018.py': '8b3de64ee1aab3e44100356bc3a4a0f0d574129992d8cf9610c822b3dc0a832f',
    'results/rank15_order15_and_raw4_pair_barrier_20261018.md': 'c84c146b5c2f8b4c5682a7fa88e7c639e1b1ca705c72caa57445f9a7af9f578c',
    'results/rank15_order15_and_raw4_pair_barrier_20261018_audit.json': 'ebf2d8bc467d94633736e3b00914ef6315c6e7607d0b8352ed572807b55bed80',
}
MUTABLE = ['README.md', 'RESEARCH_PLAN.md']


def digest(name):
    return sha256((ROOT / name).read_bytes()).hexdigest()


def read(stem):
    return json.loads((ROOT / f'results/{stem}.json').read_text())


def sources(stem, data):
    return data['prepared']['sources'] if stem == MODEL else data['sources']


def links_exist(document, target):
    path = (ROOT / document).parent / unquote(target.split('#', 1)[0])
    assert path.exists() or path.resolve() == (ROOT / AUDIT).resolve(), (document, target)


def headlines(docs):
    p = docs[PREDECESSOR]
    assert p['status'] == 'PASS_COMPLETE_ORDER11_ORDER16_DIMENSION60_COMPOSITION'
    assert (p['dimension_bound'], p['excluded_dimension'], p['inherited_dimension_bound']) == (60, 61, 61)
    assert p['finite_maximum'] == p['finite_curvature'] == '420'
    assert p['actual_ideal_absolute_surplus'] == '64' and p['strict_global_constant'] == '4650'
    assert p['contradiction'] == '0 >= 0*v_f + 64*v_i + 4650 > 0'
    assert not any(p[k] for k in ('arithmeticity_required', 'compactness_required', 'simple_cusps_required',
        'facet_count_bound_required', 'external_refereeing', 'formal_proof_assistant_verified', 'new_finite_volume_example'))
    f, d, q, b, m, j = [docs[s] for s, *_ in PACKAGES]
    assert (f['old_cost'], f['new_cost'], f['saving']) == ('119', '235/2', '3/2')
    assert f['zero_pairs'] == [[4, 14], [4, 15], [4, 16]]
    assert f['retained_D16_alternative_pairs'] == [[1, 15], [1, 16]]
    support = f['complete_selected_support']
    assert (support['complete_patterns'], support['validated_saved_maps']) == (140, 168)
    assert len(support['affected']) == len(f['selected_proofs']) == 6
    assert f['upgraded_coupon']['component_pair_cost_increase'] == '3'
    assert f['upgraded_coupon']['small_marks'] == [3, 4] and f['upgraded_coupon']['large_marks'] == [9, 10, 11]
    assert f['upgraded_coupon']['D5_long_arm_endpoint_unchanged']
    assert f['literal_configurations'] == 216 and f['literal_changed_corners'] == {}
    assert f['all_ideal_weights_unchanged'] and f['compatible_finite7_feature_depth'] == 10
    assert f['current_order15_maximum'] == '92' and not f['finite_forest_bound_proved']
    assert (d['old_cost'], d['new_cost'], d['new_order17_envelope']) == ('116', '114', '114/17')
    assert d['new_order17_costs'] == dict(A=85, B=111, C=108, D=114)
    assert d['changed_pairs'] == [[2, 15], [2, 16], [10, 13], [10, 14]]
    assert (d['old_weight'], d['new_weight']) == ('1/2', '0')
    assert (d['complete_eligible_patterns'], d['affected_maps']) == (140, 32)
    assert d['affected_contexts'] == [1978, 2371, 2535, 2539]
    assert d['upgraded_coupon']['component_pair_cost_increase'] == '1'
    assert d['upgraded_coupon']['small_marks'] == [0] and d['upgraded_coupon']['large_marks'] == [10, 11]
    assert not d['uses_same_component_finite_payments'] and not d['finite_budget_proved']
    assert d['unchanged_order15_costs'] == dict(A=75, B=92, C=91, D=92)
    for x in (f, d):
        assert x['triangle_minimum'] == 2 and x['quadrilateral_minimum'] == 1
        assert x['new_dimension_bound'] is None
    assert (q['old_cost'], q['new_cost'], q['saving'], q['new_order17_envelope']) == ('114', '227/2', '1/2', '227/34')
    assert (q['old_weight'], q['new_weight']) == ('1/2', '1/4') and q['changed_pairs'] == [[0, 10], [2, 12]]
    assert q['affected_contexts'] == [1314, 2529] and q['affected_maps'] == 16 and q['complete_eligible_patterns'] == 140
    assert q['upgraded_coupon']['component_pair_cost_increase'] == '1/4'
    assert q['upgraded_coupon']['small_marks'] == [0] and q['upgraded_coupon']['large_marks'] == [9]
    assert q['upgraded_coupon']['simultaneous_finite_transfer_total_component_pair_cost'] == '47/4'
    assert q['compatible_finite_D17_cost'] == '235/2' and not q['uses_same_component_finite_payments']
    assert q['triangle_minimum'] == 2 and q['quadrilateral_minimum'] == 1
    assert not q['finite_budget_proved'] and q['new_dimension_bound'] is None
    assert (b['completed_unrestricted_dimension_bound'], b['ambient_dimension'], b['link_dimension']) == (60, 60, 59)
    assert b['current_ideal_order17_costs'] == q['new_order17_costs'] == dict(A='85', B='111', C='108', D='227/2')
    assert (b['combined_D5_D12_increment'], b['combined_D5_D12_cost']) == ('17/4', '47/4')
    assert b['disjoint_coupon_marked_pairs'] == dict(finite=[[3, 9], [3, 10], [3, 11], [4, 9], [4, 10], [4, 11]],
        ideal_leaf=[[0, 10], [0, 11]], ideal_center=[[0, 9]])
    finite = b['finite']; ideal = b['ideal']; summary = ideal['summary']
    assert (finite['exact_maximum'], finite['previous_maximum'], finite['curvature'], finite['finite_deficit']) == ('411', '413', '405', '6')
    assert (finite['complete_fixed_vectors'], finite['interaction_keys']) == (1044, 792)
    assert (finite['simultaneous_type_pairs'], finite['residual_item_checks'], finite['Bellman_upper_inequalities'],
        finite['attaining_predecessors'], finite['above_curvature']) == (33489, 54288, 1142064, 47520, 16)
    assert finite['actual_attainer'] == ['D16', 'D16', 'D28']
    assert (finite['D17_cost'], finite['D5_D12_actual_cost'], finite['finite_cost_unit']) == ('235/2', '47/4', '1/4')
    assert finite['deep7_deletions_counted_in_upper'] and finite['new_tables_generated'] and not finite['replay_optimizes']
    assert (ideal['complete_profiles'], ideal['previous_negative_profiles'], ideal['resolved_negative_profiles']) == (2113, 49, 21)
    assert (summary['negative_count'], summary['zero_count'], summary['normalized_minimum']) == (28, 0, '-1117/7680')
    assert len(summary['negative_rows']) == 28 and all(16 not in r['ranks'] for r in summary['negative_rows'])
    assert ideal['low9_concentration_only'] and not ideal['negative_terminal_absolute_margin_is_universal_bound']
    conditional = ideal['remaining_conditional_target']['summary']
    assert conditional['negative_count'] == conditional['zero_count'] == 0
    assert (conditional['normalized_minimum'], conditional['absolute_minimum']) == ('73/104000', '54')
    assert b['combined_literal_interface']['configurations'] == len(b['combined_literal_interface']['rows']) == 216
    assert b['combined_literal_interface']['changed_configurations'] == len(b['combined_literal_interface']['changed']) == 6
    assert b['strict_global_constant'] == '4500' and b['constant_counted_once']
    assert b['new_dimension_bound'] is None and not b['new_finite_volume_example']
    assert m['solver_invocations'] == 1 and m['solver']['success'] and m['solver']['status'] == 0
    prepared = m['prepared']; counts = prepared['counts']; audit = prepared['local_audit']
    assert counts == dict(active_variables=173996, appended_local_columns=203, face_matrix_nonzeros=829446,
        faces=82806, inherited_variables=439619, inherited_vertex_cost_signed_variables=52,
        old_dual_violated_new_columns=47, preserved_multi_owner_columns=14, reused_local_columns=444,
        variables=439822, vertex_forests=6558, zero_budget_signed_variables=96176)
    assert (audit['saved_charge_certificates'], audit['saved_charge_inequalities'], audit['zero_charge_records'],
        audit['distinct_local_charges'], audit['old_dual_violated_exports']) == (738, 3273, 59, 647, 47)
    assert prepared['dimension'] == 13 and prepared['ground_field'] == 'Q(sqrt5)' and prepared['target'] == 12
    primal, dual = m['certificates']
    assert primal['kind'] == 'EXACT_PARTIAL_MODEL_PRIMAL' and dual['kind'] == 'EXACT_PARTIAL_MODEL_DUAL'
    assert primal['maximum'] == '3088979902461/250000000000' and not primal['meets_target12']
    assert dual['bound'] == '222628992286411620/18018021785910581' and dual['excludes_target12']
    assert len(dual['face_multipliers']) == 5073 and len(dual['vertex_multipliers']) == 1095
    assert j['deleted_walls'] == [10, 30, 65] and j['retained_old_walls'] == 60
    debit = j['exact_deletion_debit']
    assert debit['lost_finite_count'] == 80 and len(j['additional_lost_supports']) == 5
    assert debit['compact_faces_by_dimension'] == {'2': 314, '3': 325, '4': 244, '5': 129, '6': 46, '7': 10, '8': 1}
    assert debit['compact_S6_sides_touched'] == [8, 10, 11, 44, 55, 63, 64, 65]
    assert tuple(j['topology'][k] for k in ('vertices', 'edges', 'faces')) == (7, 13, 7)
    assert len(j['selected_prefixes']) == len(j['records']) == 6
    assert j['counts'] == dict(J_equations=260, J_branches=520)
    assert j['complete_native_triple_count'] == 216 and j['coherent_triples'] == []
    caps = [c for record in j['records'] for c in record['J_bank']['caps']]
    assert len(caps) == 216
    assert all(c['raw_whole_N'] and not c['raw_whole_J'] and c['score'] is None for c in caps)
    assert all(c['failed_disk_corners'] == [[8, 52, 69], [52, 64, 69], [63, 64, 69]] for c in caps)
    assert len(j['unresolved']) == 223 and len(j['boundary_ledgers']) == 2
    assert all(x['counts'] == dict(finite=26, exterior=6) and x['open_ends'] == 10 and x['warnings'] == []
        for x in j['boundary_ledgers'])
    assert j['control_score'] == 0 and j['deletion_score'] is None
    return dict(completed_unrestricted_dimension_bound=60, new_dimension_bound=None,
        current_finite_rank60_maximum='411', finite_curvature='405', current_ideal_order17_maximum='227/2',
        ideal_terminal_profiles=2113, negative_ideal_profiles=28, additional_target_is_conditional=True,
        current_arithmetic_classes=324, current_arithmetic_forests=600,
        changed_global_model_variables=counts['variables'], new_global_columns=203,
        exact_global_model_lower=dual['bound'], exact_global_model_upper=primal['maximum'],
        target12_feasibility_proved=False, new_positive_support_geometry_complete=False,
        complete_native_construction_triples=216, coherent_new_construction_triples=0,
        new_finite_volume_example=False, external_refereeing=False, formal_proof_assistant_verified=False)


def collect():
    assert all('PENDING' not in x for row in PACKAGES for x in row), 'Await final packet execution/review records.'
    for p, h in FROZEN_PREDECESSORS.items():
        assert digest(p) == h, ('predecessor', p)
    docs = {PREDECESSOR: read(PREDECESSOR)}; bound = {}; packages = []; bindings = 0
    for stem, status, generation, replay, review in PACKAGES:
        d = read(stem); assert d['status'] == status, (stem, d['status']); docs[stem] = d
        packages.append(dict(package=stem, status=status, direct_bindings=len(sources(stem, d)),
            generation=dict(identifier=generation, exit_code=0, status='PASS'),
            fresh_replay=dict(identifier=replay, exit_code=0, status='PASS'), review_scope=review))
    for stem, d in docs.items():
        src = sources(stem, d); assert src
        for p, h in src.items():
            assert digest(p) == h, ('binding', stem, p)
            assert p not in bound or bound[p] == h; bound[p] = h
        bindings += len(src)
    return docs, bound, packages, bindings


def artifacts():
    return list(dict.fromkeys([p for stem, *_ in PACKAGES for p in
        (f'scripts/{stem}.py', f'results/{stem}.json', f'results/{stem}.md')] + HELPERS +
        list(FROZEN_PREDECESSORS) + [CHECKPOINT, str(Path(__file__).relative_to(ROOT))]))


def build():
    docs, bound, packages, bindings = collect(); values = headlines(docs)
    paths = artifacts(); asts = [p for p in paths if p.endswith('.py')]
    for p in asts: ast.parse((ROOT / p).read_text(), filename=p)
    reports = [f'results/{s}.md' for s, *_ in PACKAGES] + [p for p in FROZEN_PREDECESSORS if p.endswith('.md')] + [CHECKPOINT] + MUTABLE
    links = []
    for p in reports:
        for raw in re.findall(r'\]\(([^\n]+?)\)', (ROOT / p).read_text()):
            target = raw.split(' "', 1)[0].strip('<>')
            if target.startswith(('http:', 'https:', 'mailto:', '#')): continue
            links_exist(p, target); links.append(dict(document=p, target=target))
    return dict(status=STATUS, date='2026-10-18', checkpoint=CHECKPOINT, saved_headlines=values,
        packages=packages, predecessor=PREDECESSOR,
        counts=dict(completed_packages=len(packages), direct_bindings=bindings,
            predecessor_direct_bindings=len(sources(PREDECESSOR, docs[PREDECESSOR])), distinct_bound_files=len(bound),
            python_asts=len(asts), documents=len(reports), local_links=len(links)),
        artifact_sha256={p: digest(p) for p in paths}, bound_file_sha256=bound,
        mutable_document_snapshots={p: digest(p) for p in MUTABLE}, python_ast_files=asts, local_links=links,
        scope='Integrity, exact saved headline/status comparisons, syntax and local links only. Prior run identifiers and completed reviews are recorded evidence. No mathematical, physical, numerical, optimization, classification, certificate or proof-assistant replay. Mutable document hashes are historical snapshots.')


def verify(data):
    assert data['status'] == STATUS and data['checkpoint'] == CHECKPOINT
    assert set(data['artifact_sha256']) == set(artifacts())
    assert data['python_ast_files'] == [p for p in artifacts() if p.endswith('.py')]
    for group in ('artifact_sha256', 'bound_file_sha256'):
        for p, h in data[group].items(): assert digest(p) == h, (group, p)
    docs, bound, packages, bindings = collect()
    assert bound == data['bound_file_sha256'] and packages == data['packages'] and headlines(docs) == data['saved_headlines']
    assert data['counts'] == dict(completed_packages=6, direct_bindings=bindings,
        predecessor_direct_bindings=len(sources(PREDECESSOR, docs[PREDECESSOR])), distinct_bound_files=len(bound),
        python_asts=len(data['python_ast_files']), documents=11, local_links=len(data['local_links']))
    assert set(data['mutable_document_snapshots']) == set(MUTABLE)
    for p in data['python_ast_files']: ast.parse((ROOT / p).read_text(), filename=p)
    for r in data['local_links']: links_exist(r['document'], r['target'])


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    path = ROOT / AUDIT
    if args.verify:
        result = json.loads(path.read_text()); verify(result)
    else:
        assert not path.exists(), 'Frozen audit exists; use --verify.'
        result = build(); path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['counts'])
