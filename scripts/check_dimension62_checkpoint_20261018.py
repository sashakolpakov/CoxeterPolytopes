#!/usr/bin/env python3
"""Integrity, syntax and links for the dimension62 checkpoint, not proof replay.

Current README/RESEARCH_PLAN hashes are historical snapshots. Frozen artifacts
and recorded links remain checked when those current-status documents advance.
"""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/dimension62_sparse_coupons_and_wall66_20261018.md'
AUDIT = CHECKPOINT.replace('.md', '_audit.json')
PACKAGES = [
    ('fresh_order11_c_endpoint_zero_transfer_20261018', '08ec1b', '555738',
     'Root and a separate full source/report mathematical read of complete B10 attachment profiles, all-native bounds, R10 identity, nonnegative spectators and actual-corner priority. Second read arrived after the local report froze. No reviewer repeated the matrix arithmetic or finite/affine classification.'),
    ('fresh_order12_b_center_terminal_transfer_20261018', '5661ec', '494d67',
     'Separate full source/report mathematical read of complete principal extensions, R9 identity, native gap, coordinate Cauchy, whole A1 sealing, physical gauge and exact finite coupon roles. No reviewer matrix replay or finite optimization.'),
    ('fresh_order13_d_center_pair_half_transfer_20261018', '04dc72', '6a0f6d',
     'Root and two separate full source/report mathematical reads of complete D10 extensions, R8 identity, center extremes, whole A2/A3-middle alternatives and retained affine center units. Third read arrived after the local report froze. No reviewer formula, inertia or arbitrary-rank classification replay.'),
    ('fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018', '42537', '28025',
     'Root and a separate full source/report mathematical read, including inherited selected/literal routines and sealing. Reviewed proper completions, all-native extra contacts, radical whole-factor control, six coupons, restored-unit obligations and deep11 interfaces. No reviewer census or literal arithmetic replay.'),
    ('fresh_sparse_coupon_rank63_forest_budget_20261018', 'b7b486', '6c791d',
     'Root and two separate full source/report mathematical reads of nine exact coupon roles, simultaneous two-end type domination, distinct B9/E6/A17, complete multiplicity domain, conservative E6/D11 allowance, Bellman upper/attainment certificates and actual449 witness. Third read arrived after the local report froze. No reviewer Bellman replay.'),
    ('fresh_ambient63_ideal_margin_closure_20261018', 'd7483c', '14112e',
     'Independent implementation reconstructs all2204 profiles and compares every scalar field to the earlier calculation. Two separate full mathematical/source reads of concentration domain, cost preservation, V monotonicity and positive absolute bound; second arrived after local freeze. No reviewer scalar replay.'),
    ('fresh_sparse_coupon_order11_13_dimension62_20261018', '4b3b23', 'e48d4f',
     'Root and two separate full source/report mathematical composition reads of the exhaustive actual-corner split, protected outputs, inherited sparse coverage, triangle/deep11 interfaces, whole-component roles, nine coupon increments, exact finite/ideal comparison, single extremum constant and completed63 higher tail. No reviewer historical census, Bellman or2204-profile replay.'),
    ('fresh_wall66_single_edge_conic_20261018', '87e83b', 'd8fb39',
     'Root full producer/report/helper review and separate disk-topology review. Independent saved-record checker31406e verifies saved coefficients, native gates, zero deletion, disk links and all ledgers/accounting. No independent normal reconstruction or numerical cone replay; no exact new normal or finite-volume example.'),
]
HELPERS = ['scripts/check_wall66_single_edge_conic_saved_20261018.py']
PREDECESSORS = {
    'results/dimension63_order11_order17_and_qj_20261018.md':
        '21be7c1bdad7717af53c5f9382b643a73fdf8f10136d5e97f1c5a82e58e3ebf0',
    'results/dimension63_order11_order17_and_qj_20261018_audit.json':
        '113930b4f258be15d8a5dea1713bd3e2ec4da9064f0f83393303fe25960fde77',
}
MUTABLE = ['README.md', 'RESEARCH_PLAN.md']


def digest(name):
    return sha256((ROOT/name).read_bytes()).hexdigest()


def check_link(document, target):
    path = (ROOT/document).parent/unquote(target.split('#', 1)[0])
    assert path.exists() or path.resolve() == (ROOT/AUDIT).resolve(), (document, target)


def build():
    assert all('PENDING' not in str(row) for row in PACKAGES)
    for name, value in PREDECESSORS.items():
        assert digest(name) == value, ('changed predecessor', name)
    bound, packages, artifacts = {}, [], []
    bindings = 0
    for stem, generation, replay, review in PACKAGES:
        names = [f'scripts/{stem}.py', f'results/{stem}.json', f'results/{stem}.md']
        data = json.loads((ROOT/names[1]).read_text())
        if stem == 'fresh_wall66_single_edge_conic_20261018':
            assert data['status'] == 'NUMERICAL_WALL66_SINGLE_EDGE_CONIC'
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
    artifacts += HELPERS+list(PREDECESSORS)+[CHECKPOINT, str(Path(__file__).relative_to(ROOT))]
    artifacts = list(dict.fromkeys(artifacts))
    asts = [name for name in artifacts if name.endswith('.py')]
    for name in asts:
        ast.parse((ROOT/name).read_text(), filename=name)
    docs = [f'results/{stem}.md' for stem, *_ in PACKAGES]+[CHECKPOINT]+MUTABLE
    links = []
    for name in docs:
        for raw in re.findall(r'\]\(([^\n]+?)\)', (ROOT/name).read_text()):
            target = raw.split(' "', 1)[0].strip('<>')
            if target.startswith(('http:', 'https:', 'mailto:', '#')):
                continue
            check_link(name, target)
            links.append(dict(document=name, target=target))
    return dict(status='PASS_DIMENSION62_CHECKPOINT_INTEGRITY_SYNTAX_AND_LINKS',
        date='2026-10-18', checkpoint=CHECKPOINT, packages=packages,
        counts=dict(completed_packages=len(packages), direct_bindings=bindings,
                    distinct_bound_files=len(bound), python_asts=len(asts),
                    documents=len(docs), local_links=len(links)),
        artifact_sha256={name: digest(name) for name in artifacts},
        bound_file_sha256=bound,
        mutable_document_snapshots={name: digest(name) for name in MUTABLE},
        python_ast_files=asts, local_links=links,
        scope='Integrity, syntax and recorded local links only. Run identifiers record execution evidence, not reexecution by this audit. No additional mathematical, geometric, numerical or proof-assistant verification. Current-document hashes are historical snapshots.')


def verify(data):
    assert data['status'] == 'PASS_DIMENSION62_CHECKPOINT_INTEGRITY_SYNTAX_AND_LINKS'
    for group in ('artifact_sha256', 'bound_file_sha256'):
        for name, value in data[group].items():
            assert digest(name) == value, (group, name)
    for name in data['python_ast_files']:
        ast.parse((ROOT/name).read_text(), filename=name)
    for row in data['local_links']:
        check_link(row['document'], row['target'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = ROOT/AUDIT
    if args.verify:
        result = json.loads(path.read_text())
        verify(result)
    else:
        assert not path.exists(), 'Frozen audit exists; use --verify.'
        result = build()
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['counts'])
