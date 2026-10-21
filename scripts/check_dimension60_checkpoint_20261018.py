#!/usr/bin/env python3
"""Check checkpoint integrity and saved conclusions; no mathematical replay."""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/dimension60_complete_and_current_search_20261018.md'
AUDIT = CHECKPOINT.replace('.md', '_audit.json')
STATUS = 'PASS_DIMENSION60_COMPLETE_AND_CURRENT_SEARCH_CHECKPOINT'
PACKAGES = [
    ('fresh_finite_d16_fork_exchange_20261018', 'ba04f0', 'f4242e',
     'Author and two full source/report mathematical reads: arbitrary contacts, radical seals, D17 restoration, explicit A4/A12 charge, complete saved maps and current216 literal interfaces. No reviewer map, inertia or historical census replay.'),
    ('fresh_rank61_complete_finite420_budget_20261018', '37b845', '01904f',
     'Author and two full source/report mathematical reads: all actual/conservative charges, simultaneous two-end domination, complete1132/854 domain, new item costs, last-item upper/attainment induction and actual420 witnesses. No reviewer DP, scalar or census replay.'),
    ('fresh_order11_bc_zero_transfer_20261018', '287a9e', 'b625d5',
     'Author and peers full source/report and inherited profile reads: complete C contact domain, strict gap, B extreme/equality cases, whole-component sealing, distinct retained outputs, coupon increment and noncircular corner priorities. No reviewer profile, determinant or literal arithmetic replay.'),
    ('fresh_order11_order16_dimension60_20261018', '29dc73', '913e99',
     'Author and two full mathematical composition reads: C/B-first priorities, D105 whole/adaptive outputs, current216 obligations, finite420 plus A1/B10 absorption, complete terminal domain, positive absolute extension, factor fibers, one4650 correction and completed61 higher tail. No reviewer census, DP, convolution or2365 scalar replay.'),
    ('fresh_new5121_raw4_two_face_domain_20261018', '5193b3', '0a054d',
     'Author and peer full source/report and field/projection helper reads: current-domain exclusion, exact provenance, complete selected embeddings/orientations, immutable native field decoding, positive heights and native products. No reviewer embedding, field arithmetic or classification replay. No compatibility or cover update.'),
    ('fresh_joint30_33_hinge_native_pilot_20261018', '2b5dd4', '74735e',
     'Root and peer full producer/report/standalone-helper reads; peer also read inherited profile/native/corner/scoring/deletion interfaces. Corrected actual-side frame profiles, exact15/10/1 debit, bounded selection,13 provisional incidence diagnostics versus15 clean final pairs,14 failed closing duties and sole old restoration. Helper1c2af8 checks saved records; reviewers did not repeat normal, classifier or cone computation.'),
]
HELPERS = ['scripts/check_joint30_33_hinge_saved_20261018.py']
PREDECESSORS = {
    'results/depth8_finite425_and_joint_motion_frontier_20261018.md':
        '72b01763e658078930b5b310ec4ad2d241bfdd997bb84c39744312c79692bfe8',
    'results/depth8_finite425_and_joint_motion_frontier_20261018_audit.json':
        'a81a7991f2e36031896785e81932b2957acf8c747a60a3debca51f9df5a65ac8',
    'scripts/check_depth8_finite425_checkpoint_20261018.py':
        '2f26436462ce4dcbc867e5269b5333c4bf0a1833325079914cc863b9fe9d0345',
}
MUTABLE = ['README.md', 'RESEARCH_PLAN.md']


def digest(name):
    return sha256((ROOT / name).read_bytes()).hexdigest()


def link_exists(document, target):
    path = (ROOT / document).parent / unquote(target.split('#', 1)[0])
    assert path.exists() or path.resolve() == (ROOT / AUDIT).resolve(), (document, target)


def headlines(docs):
    d = docs[PACKAGES[0][0]]
    assert d['new_costs'] == dict(D16='111', D17='119')
    d = docs[PACKAGES[1][0]]
    assert d['finite_maximum'] == d['finite_curvature'] == '420'
    assert (d['complete_large_vectors'], d['interaction_keys']) == (1132, 854)
    assert d['new_residual_tables'] and d['optional_A1_B10_half_to_unit_absorption']
    d = docs[PACKAGES[2][0]]
    assert d['new_family_totals'] == dict(A='99/2', B='53', C='52', D='53')
    assert d['coupon_upgrade']['added_component_pair_cost'] == '1/2'
    d = docs[PACKAGES[3][0]]
    assert (d['dimension_bound'], d['excluded_dimension'], d['inherited_dimension_bound']) == (60, 61, 61)
    assert d['complete_terminal_ideal_profiles'] == 2365
    assert d['normalized_ideal_minimum'] == '3/250'
    assert d['actual_ideal_absolute_surplus'] == '64'
    assert d['finite_surplus'] == '0' and d['strict_global_constant'] == '4650'
    assert d['contradiction'] == '0 >= 0*v_f + 64*v_i + 4650 > 0'
    assert not any(d[k] for k in ('arithmeticity_required', 'compactness_required',
        'simple_cusps_required', 'facet_count_bound_required', 'new_finite_volume_example',
        'external_refereeing', 'formal_proof_assistant_verified'))
    d = docs[PACKAGES[4][0]]
    assert d['counts'] == dict(active_embeddings=56, actual_corner_occurrences=12,
        distinct_native_stars=28, faces=2, literal_orientation_occurrences=40, vertex_types=6)
    assert (d['baseline_represented_classes'], d['baseline_represented_forests']) == (322, 600)
    assert d['pair_compatibility_tests'] == 0 and not d['new_local_cover']
    assert d['new_dimension_bound'] is None
    d = docs[PACKAGES[5][0]]
    assert d['counts'] == dict(norm_equations=397, norm_branches=794)
    assert d['exact_deletion_debit']['lost_finite_count'] == 15
    assert d['exact_deletion_debit']['compact_faces_by_dimension'] == {'2': 10, '3': 1}
    assert len(d['first_bank']['caps']) == 12 and len(d['records']) == 9
    assert (len(d['accepted']), len(d['required_corner_failures']), len(d['unresolved'])) == (1, 14, 13)
    a = d['accepted'][0]
    assert a['old33_restored'] and a['old30_restored'] and a['open_ends'] == 10
    return dict(completed_unrestricted_dimension_bound=60, finite_maximum='420',
        finite_surplus='0', ideal_absolute_surplus='64', strict_global_constant='4650',
        compact_arithmetic_dimension_bound=21, arbitrary_compact_dimension_bound=29,
        arithmetic_cover_classes=322, arithmetic_cover_vertex_types=600,
        new_raw4_stars=28, new_raw4_stars_merged=False, new_finite_volume_example=False)


def build():
    for name, expected in PREDECESSORS.items():
        assert digest(name) == expected, ('predecessor', name)
    bound, packages, artifacts, docs = {}, [], [], {}
    bindings = 0
    for stem, generation, replay, review in PACKAGES:
        names = [f'scripts/{stem}.py', f'results/{stem}.json', f'results/{stem}.md']
        d = json.loads((ROOT / names[1]).read_text())
        prefix = 'BOUNDED' if stem.startswith('fresh_joint') else 'PASS'
        assert d['status'].startswith(prefix) and d['sources']
        for name, expected in d['sources'].items():
            assert digest(name) == expected, ('binding', stem, name)
            assert name not in bound or bound[name] == expected
            bound[name] = expected
        bindings += len(d['sources'])
        artifacts.extend(names)
        docs[stem] = d
        packages.append(dict(package=stem, status=d['status'], direct_bindings=len(d['sources']),
            generation=dict(identifier=generation, exit_code=0),
            fresh_replay=dict(identifier=replay, exit_code=0), review_scope=review))
    values = headlines(docs)
    artifacts += HELPERS + list(PREDECESSORS) + [CHECKPOINT, str(Path(__file__).relative_to(ROOT))]
    artifacts = list(dict.fromkeys(artifacts))
    asts = [name for name in artifacts if name.endswith('.py')]
    for name in asts:
        ast.parse((ROOT / name).read_text(), filename=name)
    reports = [f'results/{stem}.md' for stem, *_ in PACKAGES] + [CHECKPOINT] + MUTABLE
    links = []
    for name in reports:
        for raw in re.findall(r'\]\(([^\n]+?)\)', (ROOT / name).read_text()):
            target = raw.split(' "', 1)[0].strip('<>')
            if target.startswith(('http:', 'https:', 'mailto:', '#')):
                continue
            link_exists(name, target)
            links.append(dict(document=name, target=target))
    return dict(status=STATUS, date='2026-10-18', checkpoint=CHECKPOINT,
        saved_headlines=values, packages=packages,
        counts=dict(completed_packages=len(packages), direct_bindings=bindings,
            distinct_bound_files=len(bound), python_asts=len(asts), documents=len(reports), local_links=len(links)),
        artifact_sha256={name: digest(name) for name in artifacts}, bound_file_sha256=bound,
        mutable_document_snapshots={name: digest(name) for name in MUTABLE},
        python_ast_files=asts, local_links=links,
        scope='Integrity, saved conclusions, syntax and links. Prior execution identifiers are records, not executions by this auditor. No new mathematical, numerical, geometric, optimization or proof-assistant verification. Mutable document hashes are historical snapshots.')


def verify(data):
    assert data['status'] == STATUS
    for group in ('artifact_sha256', 'bound_file_sha256'):
        for name, expected in data[group].items():
            assert digest(name) == expected, (group, name)
    for name in data['python_ast_files']:
        ast.parse((ROOT / name).read_text(), filename=name)
    for row in data['local_links']:
        link_exists(row['document'], row['target'])
    docs = {stem: json.loads((ROOT / f'results/{stem}.json').read_text()) for stem, *_ in PACKAGES}
    assert headlines(docs) == data['saved_headlines']


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    path = ROOT / AUDIT
    if args.verify:
        result = json.loads(path.read_text())
        verify(result)
    else:
        assert not path.exists(), 'Frozen audit exists; use --verify.'
        result = build()
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['counts'])
