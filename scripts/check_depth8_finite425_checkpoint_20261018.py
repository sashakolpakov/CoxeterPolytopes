#!/usr/bin/env python3
"""Checkpoint bindings, saved headlines, syntax and links; no proof replay."""
import argparse
import ast
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/depth8_finite425_and_joint_motion_frontier_20261018.md'
AUDIT = CHECKPOINT.replace('.md', '_audit.json')
STATUS = 'PASS_DEPTH8_FINITE425_LOCAL_FRONTIER_CHECKPOINT'
PACKAGES = [
    ('fresh_ambient61_payment_frontier_20261018', '811515', 'dffc45',
     'Complete saved finite/ideal frontier with scoped mathematical review; conditional changes are not local theorems. This is the preceding ideal21-negative baseline.'),
    ('fresh_depth8_finite8_zero_transfer_20261018', 'a91fd9', '9f838f',
     'Author and two full source/report mathematical reads: retained-contact cases including disconnected deletion, norm-subtracted Schur estimate, distinct Q9/Q10/tail domains and current priority interfaces. Reviewers did not repeat scalar arithmetic or inherited classification.'),
    ('fresh_finite_d17_a4_d13_half_transfer_20261018', '5bc07e', '27374f',
     'Root and construction peer full source/report mathematical reads: full radical seals, existing A4/D13 halves, complete selected maps and latest216 literal debits. No reviewer map, inertia, literal or historical census replay.'),
    ('fresh_depth8_d17_rank61_forest_budget_20261018', '8cf16d', 'dad836',
     'Root and two peers full mathematical source/report reads: actual/conservative D16 distinction, all typed coupons, simultaneous domination, complete1132/854 domain, every-last-item upper/attainment proof, deep8 counts and actual425 attainer. No reviewer DP or scalar replay.'),
    ('fresh_finite_a15_endpoint_coupon_zero_20261018', 'fe1988', 'ada911',
     'Author and two full source/report mathematical reads: two radical whole-component seals, explicit A2/A13 endpoint coupon increment2, complete saved536-pattern/482-map domain and216 current literal interfaces. No reviewer embedding, inertia or literal arithmetic replay.'),
    ('fresh_cross_vertex_support_migration_20261018', '71fc96', 'c4058a',
     'Saved support/provenance diagnostic at5121 masses; no new geometry, local optimization or clique enumeration.'),
    ('fresh_cross_vertex_current290_local_covers_20261018', 'f12278', 'c98ffe',
     'Root and peer full source/report reads; peer reread exact local-cover helper. Complete current domains at new demands, exact certificates and both-owner costs checked at recorded scopes. No reviewer LP or physical census replay.'),
    ('fresh_cross_vertex_historical32_local_covers_20261018', '3782b2', '55642f',
     'Root and peer full source/report mathematical reads; peer reread relation, current-domain, preparation and exact-cover helpers. Complete historical ownership, conservative unknowns, all affected clique domains, exact covers and preserved two-owner costs. No reviewer clique, LP or physical replay.'),
    ('fresh_wall59_remaining_low_seed_joins_20261018', 'a648c4', '17642d',
     'Root full producer/report/saved-helper read. Helper f49da5 independently checks saved normal products,60 pair hashes and exact deletion debit; only12 closest pairs have full cone scoring. No root normal or cone replay.'),
    ('fresh_wall66_ideal_junction_pilot_20261018', 'ea4921', 'd37d3b',
     'Root full producer/report read; peer full producer/report/shared-helper read. Shared helper de3a7e checks saved domains. Affine-incidence diagnostics remain unresolved and raw cone counts are not accepted boundary counts.'),
    ('fresh_wall66_boundary_bridge_flip_20261018', 'b665dc', '6c4b73',
     'Root full producer/report read; peer full producer/report/shared-helper and inherited low-bank read. Shared helper de3a7e checks saved records. Changed J frame and mandatory Y/J join, without fictitious66/67 adjacency. No reviewer physical-normal or cone replay.'),
    ('fresh_wall59_nonright55_prefix_20261018', '44412c', '8a8180',
     'Root full producer/report/shared-helper read. Helper9fc2d3 independently reconstructs exact199-support deletion debit and saved numerical ledgers. New normals and cones remain numerical.'),
    ('fresh_wall59_nonright55_companions_20261018', '47e5e3', 'd0163d',
     'Root full producer/report/shared-helper read. Helper9fc2d3 independently checks all2968 saved normal products and15 digests; all190 single caps but only3 closest pairs have complete cone scoring. No reviewer normal or cone replay.'),
    ('fresh_order16_d_distance9_half_transfer_20261018', '0be232', '635003',
     'Root and construction peer full source/report plus inherited core_constraints/complete reads; root reread residual/no-bridge proof and universal tree/cycle interfaces. All394 target domains,24 contexts/128 maps, arbitrary first-entry contacts, sealed wholeD5/D11 half and completed61/D17=120 literal baseline. No reviewer arithmetic, map, inertia or historical census replay.'),
    ('fresh_finite_b15_d15_existing_payment_transfer_20261018', '0dda0b', 'b335e4',
     'Root and peer full source/report plus inherited contact-domain reads. Four forced7/2 new-white projections, four wholeD15 seals, existing A2/D13 short-leaf halves, complete saved maps and216 simultaneous literal interfaces. No reviewer arithmetic, map or historical census replay.'),
    ('fresh_finite_a16_endpoint_sparse_half_transfer_20261018', '36d81b', '7f59e3',
     'Root and peer full source/report mathematical reads: four19-root radical seals, distinct whole-component coupons with additional costs1 and2, complete394/238 saved domain and216 unchanged literal interfaces. No reviewer arithmetic or map replay.'),
    ('fresh_joint66_67_native_prefixes_20261018', '9eab8e', '33e8ce',
     'Root full producer/report/standalone-helper read. Helper72d323 independently checks exact9/6/1 deletion debit, all350 saved normal products,299 complete saved boundary records and stored cone/PD equality. No root normal, classifier, arithmetic or cone replay.'),
    ('fresh_order16_d105_ambient61_frontier_20261018', 'cc70c7', '530d73',
     'Root full bounded source/report read. Complete independently reconstructed2365sorted-part domain, exact rank15 saving, four remaining negatives, strict order11 threshold and explicit current family-cost provenance. No root scalar, convolution or geometric census replay.'),
]
HELPERS = [
    'scripts/check_wall59_remaining_low_seed_saved_20261018.py',
    'scripts/check_wall66_ideal_and_flip_saved_20261018.py',
    'scripts/check_wall59_nonright55_saved_20261018.py',
    'scripts/check_joint66_67_saved_20261018.py',
]
PREDECESSORS = {
    'results/ambient62_ideal_closure_and_current_search_20261018.md':
        '62d6512871adb45dc0de0e717e3aefad556d9d5cd0927308e35ad1b653c00d90',
    'results/ambient62_ideal_closure_and_current_search_20261018_audit.json':
        '556836df8a47275192994a1058d7481c2488621dcc1d998be65b880cad9cb447',
    'results/fresh_rank15_rank16_dimension61_20261018.md':
        'ce1a00a2eccd0c913aa14bd70bd60f80bc89ec0f6eeeba0e406c69a21a4b6587',
}
MUTABLE = ['README.md', 'RESEARCH_PLAN.md']


def digest(name):
    return sha256((ROOT/name).read_bytes()).hexdigest()


def check_link(document, target):
    path = (ROOT/document).parent/unquote(target.split('#', 1)[0])
    assert path.exists() or path.resolve() == (ROOT/AUDIT).resolve(), (document, target)


def headlines(docs):
    d = docs['fresh_ambient61_payment_frontier_20261018']
    assert (d['inherited_completed_dimension_bound'], d['ambient_dimension'],
            d['finite_curvature'], d['strict_global_constant']) == (61, 61, '420', '4650')
    assert d['ideal']['complete_profiles'] == 2365 and d['new_dimension_bound'] is None
    d = docs['fresh_depth8_finite8_zero_transfer_20261018']
    assert (d['feature_depth'], d['cycle_ordered_pairs'], d['cycle_rectangle_values']) == (8, 77, 308)
    assert d['corridor']['strict_determinant_lower_bound'] == '1875/11264'
    d = docs['fresh_finite_d17_a4_d13_half_transfer_20261018']
    assert (d['new_cost'], d['selected_contexts'], d['affected_maps']) == ('118', 8, 16)
    assert d['new_cross_component_coupons'] == []
    d = docs['fresh_depth8_d17_rank61_forest_budget_20261018']
    assert (d['finite_maximum'], d['finite_curvature'], d['singleton_D61_cost']) == ('425', '420', '420')
    assert (d['complete_large_vectors'], d['interaction_keys'], d['Bellman_upper_inequalities']) == (1132, 854, 1276730)
    d = docs['fresh_finite_a15_endpoint_coupon_zero_20261018']
    assert d['new_cost'] == '103' and d['upgraded_coupon']['component_pair_cost_increase'] == '2'
    d = docs['fresh_cross_vertex_support_migration_20261018']
    assert len(d['face_multipliers']) == 5121 and len(d['newly_positive_historical_outside290']) == 32
    d = docs['fresh_cross_vertex_current290_local_covers_20261018']
    assert len(d['forest_records']) == 571 and len(d['columns']) == 16
    d = docs['fresh_cross_vertex_historical32_local_covers_20261018']
    assert (d['represented_classes'], d['represented_forests'], d['complete_maximal_cliques']) == (322, 600, 1597290)
    assert Fraction(d['bound']) == Fraction(205569148850452800, 16644039406300549) > 12
    assert len(d['preserved_multi_owner_columns']) == 14
    d = docs['fresh_wall59_remaining_low_seed_joins_20261018']
    assert d['complete_pair_tests'] == 18684 and d['completed_disks'] == [] and not d['unresolved']
    d = docs['fresh_wall66_ideal_junction_pilot_20261018']
    assert len(d['tests']) == 208 and len(d['native_candidates']) == 924
    d = docs['fresh_wall66_boundary_bridge_flip_20261018']
    assert len(d['YJ_tests']) == 233 and d['completions'] == []
    d = docs['fresh_wall59_nonright55_prefix_20261018']
    assert d['raw_low_domain_size'] == 336 and len(d['accepted']) == 8
    d = docs['fresh_wall59_nonright55_companions_20261018']
    assert d['complete_pair_tests'] == 2968 and d['completed_disks'] == []
    d = docs['fresh_order16_d_distance9_half_transfer_20261018']
    assert d['new_order16_maximum'] == '105' and d['new_family_totals'] == dict(A='80', B='105', C='103', D='105')
    assert (d['complete_target_maps'], d['affected_maps'], len(d['complete_target_domains'])) == (928, 128, 394)
    assert d['new_finite_coupons'] == []
    d = docs['fresh_finite_b15_d15_existing_payment_transfer_20261018']
    assert d['new_costs'] == dict(B15='103', D15='103') and d['new_cross_component_coupons'] == []
    d = docs['fresh_finite_a16_endpoint_sparse_half_transfer_20261018']
    assert d['new_cost'] == '111' and [c['component_pair_cost'] for c in d['new_cross_component_coupons']] == ['1', '2']
    d = docs['fresh_joint66_67_native_prefixes_20261018']
    assert d['complete_pair_count'] == 350 and len(d['accepted']) == 299 and not d['unresolved']
    assert d['exact_deletion_debit']['lost_finite_count'] == 9
    assert d['exact_deletion_debit']['compact_faces_by_dimension'] == {'2': 6, '3': 1}
    assert sum(not x['old66_restored'] for x in d['accepted']) == 251
    d = docs['fresh_order16_d105_ambient61_frontier_20261018']
    assert (d['complete_profiles'], d['negative_count'], d['zero_count']) == (2365, 4, 0)
    assert d['normalized_minimum'] == '-1159/42592'
    assert d['strict_order11_total_saving_threshold'] == '1159/11616'
    return dict(completed_unrestricted_dimension_bound=61, latest_complete_finite_count='425',
        finite_curvature='420', later_local_costs=dict(A15='103', B15='103', D15='103', A16='111', ideal_order16='105'),
        ideal_terminal_profiles=2365, ideal_negative_profiles_after_D105=4,
        later_local_costs_included_in_425=False, compact_arithmetic_dimension_bound=21,
        arbitrary_compact_dimension_bound=29, new_finite_volume_example=False,
        new_dimension_bound=None)


def build():
    for name, expected in PREDECESSORS.items():
        assert digest(name) == expected, ('changed predecessor', name)
    bound, packages, artifacts, docs = {}, [], [], {}
    bindings = 0
    for stem, generation, replay, review in PACKAGES:
        names = [f'scripts/{stem}.py', f'results/{stem}.json', f'results/{stem}.md']
        data = json.loads((ROOT/names[1]).read_text())
        assert data['status'].startswith('BOUNDED' if stem.startswith(('fresh_wall', 'fresh_joint')) else 'PASS'), (stem, data['status'])
        sources = data['sources']
        assert sources
        for name, expected in sources.items():
            assert digest(name) == expected, ('binding', stem, name)
            assert name not in bound or bound[name] == expected
            bound[name] = expected
        bindings += len(sources)
        artifacts.extend(names)
        docs[stem] = data
        packages.append(dict(package=stem, status=data['status'], direct_bindings=len(sources),
            generation=dict(identifier=generation, exit_code=0),
            fresh_replay=dict(identifier=replay, exit_code=0), review_scope=review))
    values = headlines(docs)
    artifacts += HELPERS + list(PREDECESSORS) + [CHECKPOINT, str(Path(__file__).relative_to(ROOT))]
    artifacts = list(dict.fromkeys(artifacts))
    asts = [name for name in artifacts if name.endswith('.py')]
    for name in asts:
        ast.parse((ROOT/name).read_text(), filename=name)
    reports = [f'results/{stem}.md' for stem, *_ in PACKAGES] + [CHECKPOINT] + MUTABLE
    links = []
    for name in reports:
        for raw in re.findall(r'\]\(([^\n]+?)\)', (ROOT/name).read_text()):
            target = raw.split(' "', 1)[0].strip('<>')
            if target.startswith(('http:', 'https:', 'mailto:', '#')):
                continue
            check_link(name, target)
            links.append(dict(document=name, target=target))
    return dict(status=STATUS, date='2026-10-18', checkpoint=CHECKPOINT,
        saved_headlines=values, packages=packages,
        counts=dict(completed_packages=len(packages), direct_bindings=bindings,
            distinct_bound_files=len(bound), python_asts=len(asts), documents=len(reports), local_links=len(links)),
        artifact_sha256={name: digest(name) for name in artifacts}, bound_file_sha256=bound,
        mutable_document_snapshots={name: digest(name) for name in MUTABLE},
        python_ast_files=asts, local_links=links,
        scope='Integrity, saved headline assertions, syntax and local links. Execution identifiers record prior evidence, not execution by this auditor. No new mathematical, geometric, numerical, optimization or proof-assistant verification. Mutable document hashes are historical snapshots.')


def verify(data):
    assert data['status'] == STATUS
    for group in ('artifact_sha256', 'bound_file_sha256'):
        for name, expected in data[group].items():
            assert digest(name) == expected, (group, name)
    for name in data['python_ast_files']:
        ast.parse((ROOT/name).read_text(), filename=name)
    for row in data['local_links']:
        check_link(row['document'], row['target'])
    docs = {stem: json.loads((ROOT/f'results/{stem}.json').read_text()) for stem, *_ in PACKAGES}
    assert headlines(docs) == data['saved_headlines']


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    path = ROOT/AUDIT
    if args.verify:
        result = json.loads(path.read_text())
        verify(result)
    else:
        assert not path.exists(), 'Frozen audit exists; use --verify.'
        result = build()
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['counts'])
