#!/usr/bin/env python3
"""Integrity, syntax and link audit for the dimension63 research checkpoint.

This is not a mathematical replay. Historical README/PLAN hashes are snapshots;
--verify checks the frozen artifacts and recorded links without requiring those
two current-status documents to remain unchanged in later research turns.
"""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/dimension63_order11_order17_and_qj_20261018.md'
AUDIT = CHECKPOINT.replace('.md', '_audit.json')
PACKAGES = [
    ('fresh_order11_bcd_half_pair_transfer_20261018', '35489d', '5544c5',
     'Two separate full source/report mathematical reviews of complete endpoint profiles, universal spectator alternatives, finite coupons and simultaneous short-path interfaces. No duplicate arithmetic or geometric census.'),
    ('fresh_order17_e6_half_coupon_transfer_20261018', 'a4fd87', 'f68e65',
     'Two full source/report mathematical reviews before freezing, plus a third after local report freeze recorded in checkpoint. Review of four whole-core alternatives, coupon increment and short8/triangle interfaces. Inherited140-pattern/96-map geometric classification not regenerated.'),
    ('fresh_e6_b10_half_coupon_forest_budget_20261018', 'dfea64', 'abf19c',
     'Root and two separate full source/report reviews of finite formulas, type domination, both new coupons and upper/attainment verification. No reviewer repeated Bellman arithmetic.'),
    ('fresh_ambient64_ideal_margin_closure_20261018', '621472', '902f1b',
     'Independent implementation reconstructs all2419 profiles and compares every scalar field to the earlier calculation. Two separate full source/report mathematical reviews of concentration domain, V monotonicity, envelopes and constants; second recorded after local report freeze. Reviewers did not rerun scalar arithmetic.'),
    ('fresh_order11_order17_dimension63_20261018', '753c2e', '0bdd94',
     'Two full source/report composition reviews of actual-corner priority split, all18 distance8 small-tree targets, cycle/deep-corridor interfaces, added coupons and global inequality. Separate saved-only tree inclusion check5502cf passed. No historical geometric census, Bellman or ideal scalar replay by reviewers.'),
    ('fresh_qj_connector_20261018', '14342', '26392',
     'Root and a separate full producer/report/helper review, second recorded after local report freeze. No-project-import saved checker ed7fb8 verifies855 recipe indices, disk topology, saved residual separation, rounded-coefficient signs and recorded cone/PD agreement. No independent normal reconstruction, numerical cone replay or spherical-classifier proof.'),
]
HELPERS = ['scripts/check_qj_connector_saved_20261018.py']
FROZEN_PREDECESSORS = {
    'results/finite464_cycle11_and_shared_edge_20261018.md':
        '3f71095d09029e38578ae0d937ac49bd21c0cb9e9f135b4a132decb580963208',
    'results/finite464_cycle11_and_shared_edge_20261018_audit.json':
        'a348f50964ae1bff115aa594b06e43e607392e8751d6450d7cd92617734e9fc1',
}
MUTABLE_DOCS = ['README.md', 'RESEARCH_PLAN.md']


def digest(name):
    return sha256((ROOT / name).read_bytes()).hexdigest()


def check_link(document, target):
    relative = unquote(target.split('#', 1)[0])
    path = (ROOT / document).parent / relative
    assert path.exists() or path.resolve() == (ROOT / AUDIT).resolve(), (document, target)


def build():
    for name, value in FROZEN_PREDECESSORS.items():
        assert digest(name) == value, ('changed predecessor', name)
    bound, packages, artifacts = {}, [], []
    bindings = 0
    for stem, generation, replay, review in PACKAGES:
        names = [f'scripts/{stem}.py', f'results/{stem}.json', f'results/{stem}.md']
        data = json.loads((ROOT / names[1]).read_text())
        if stem == 'fresh_qj_connector_20261018':
            assert data['status'] == 'NUMERICAL_QJ_CONNECTOR_NECESSARY_DOMAIN_SCOUT'
        else:
            assert data['status'].startswith('PASS'), stem
        for name, value in data['sources'].items():
            assert digest(name) == value, ('binding', stem, name)
            assert name not in bound or bound[name] == value
            bound[name] = value
        bindings += len(data['sources'])
        artifacts.extend(names)
        packages.append(dict(package=stem, status=data['status'],
            direct_bindings=len(data['sources']),
            generation=dict(identifier=generation, exit_code=0, status='PASS'),
            fresh_replay=dict(identifier=replay, exit_code=0, status='PASS'),
            review_scope=review))
    artifacts += HELPERS + list(FROZEN_PREDECESSORS) + [CHECKPOINT, str(Path(__file__).relative_to(ROOT))]
    artifacts = list(dict.fromkeys(artifacts))
    asts = [name for name in artifacts if name.endswith('.py')]
    for name in asts:
        ast.parse((ROOT / name).read_text(), filename=name)
    docs = [f'results/{stem}.md' for stem, *_ in PACKAGES] + [CHECKPOINT] + MUTABLE_DOCS
    links = []
    for name in docs:
        for raw in re.findall(r'\]\(([^\n]+?)\)', (ROOT / name).read_text()):
            target = raw.split(' "', 1)[0].strip('<>')
            if target.startswith(('http:', 'https:', 'mailto:', '#')):
                continue
            check_link(name, target)
            links.append(dict(document=name, target=target))
    return dict(status='PASS_DIMENSION63_CHECKPOINT_INTEGRITY_SYNTAX_AND_LINKS',
        date='2026-10-18', checkpoint=CHECKPOINT, packages=packages,
        counts=dict(completed_packages=len(packages), direct_bindings=bindings,
                    distinct_bound_files=len(bound), python_asts=len(asts),
                    documents=len(docs), local_links=len(links)),
        artifact_sha256={name: digest(name) for name in artifacts},
        bound_file_sha256=bound,
        mutable_document_snapshots={name: digest(name) for name in MUTABLE_DOCS},
        python_ast_files=asts, local_links=links,
        scope='Integrity, syntax and recorded local links only. Run identifiers are recorded execution evidence, not reexecuted by this audit. No additional mathematical, geometric, numerical or proof-assistant verification. Mutable current-status document hashes are historical snapshots.')


def verify(data):
    assert data['status'] == 'PASS_DIMENSION63_CHECKPOINT_INTEGRITY_SYNTAX_AND_LINKS'
    for group in ('artifact_sha256', 'bound_file_sha256'):
        for name, value in data[group].items():
            assert digest(name) == value, (group, name)
    for name in data['python_ast_files']:
        ast.parse((ROOT / name).read_text(), filename=name)
    for row in data['local_links']:
        check_link(row['document'], row['target'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = ROOT / AUDIT
    if args.verify:
        result = json.loads(path.read_text())
        verify(result)
    else:
        assert not path.exists(), 'Frozen audit exists; use --verify.'
        result = build()
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['counts'])
