#!/usr/bin/env python3
"""Integrity, syntax and links for the next ambient62 frontier checkpoint.

This audit does not replay mathematics or recorded runs. Current-document
hashes are historical snapshots; frozen artifact hashes remain enforced.
"""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/dimension62_three_transfers_raw39_and_wall59_20261018.md'
AUDIT = CHECKPOINT.replace('.md', '_audit.json')
PACKAGES = [
    ('fresh_ambient62_payment_frontier_20261018', '65d97d', '93eaa3',
     'Root and a separate full source/report mathematical read: simultaneous forest type domination, complete multiplicity/interaction domain, Bellman upper and attainment certificates, exact large-factor tail, low-only concentration, complete2641 domain and conditional scenarios. No reviewer DP or scalar replay.'),
    ('fresh_order13_d_center_pair_zero_transfer_20261018', '61a495', '892cd1',
     'Root full source/report mathematical read of inherited universal D13 branches, distinct full units, whole A2/A3-middle versus D10 roles, exact coupon increments and simultaneous priority. No reviewer formula, inertia or arbitrary-rank classification replay.'),
    ('fresh_order16_affine_b_singleton_projection_transfer_20261018', '428bb4', '809c4f',
     'Root and a separate full source/report and inherited contact-domain read: complete selected-map domain, forced simple/double contact, positive inverse end blocks, strict norm obstruction, preserved short8/deep11 and whole-face interfaces. No reviewer map or inverse replay.'),
    ('fresh_order17_e6_unit_coupon_zero_transfer_20261018', '0254a5', 'ee44e3',
     'Root and a separate full source/report mathematical read of complete inherited proper exclusions, radical whole-component sealing, distinct E6/D11 units, actual-corner priority, independent short8 route and pointwise absorption under existing allowance4. No reviewer census, matrix or Bellman replay.'),
    ('fresh_ambient62_three_transfer_frontier_20261018', '67777e', 'd53025',
     'Root and a separate full source/report mathematical composition read: independent actual-corner priority followed by B16 pointwise agreement, simultaneous D13/E6 budget absorption, complete exact-tail domain, independent direct-convolution comparisons and conditional-only order14 query. No reviewer census, DP or scalar replay.'),
    ('fresh_raw39_additional_target_frontier_20261018', 'e4e740', '55afaf',
     'Root and a separate full producer/report/helper mathematical read including inherited actual-face and native recovery interfaces, actual source-pair events, complete embedding/event/pair domains, unique unlabelled recovery, four-type equalities, conservative graphs and exact covering/charge certificates. Independent saved helper3f1159 checks all new arithmetic, graphs and certificates at its stated scope. Source reviewers did not rerun geometry, graphs, optimization or arithmetic.'),
    ('fresh_wall59_bond_replacement_20261018', 'e0a3b8', 'ef15b2',
     'Root full producer/report/helper read and a separate topology review. Independent saved helperc43eb7 checks deletion, complete sample domain, saved native products, disk links and ledgers. No independent physical-normal reconstruction or numerical cone replay; no exact new normal or improved chamber.'),
    ('fresh_wall59_low_bond_companions_20261018', '7cedd3', '85750a',
     'Root full producer/report/helper read and a separate topology review. Independent saved helperc43eb7 checks272 caps, all7200 pair domains, saved native products, disk links and cone-ledger bookkeeping. No finite pair junction occurs, so no dynamic third-center stage is reached. No independent normal or cone replay and no finite-volume example.'),
]
HELPERS = ['scripts/check_raw39_additional_targets_saved_20261018.py',
           'scripts/check_wall59_replacement_saved_20261018.py']
PREDECESSORS = {
    'results/dimension62_sparse_coupons_and_wall66_20261018.md':
        '03f290ac6c9713d12fac08e2649abdacf67580f4cf238331921686009a1f6e58',
    'results/dimension62_sparse_coupons_and_wall66_20261018_audit.json':
        'ba6274dfc75ff7fa69774587a92eec6db76f0fc271e176162126f3880b063175',
}
MUTABLE = ['README.md', 'RESEARCH_PLAN.md']
STATUS = 'PASS_DIMENSION62_FRONTIER_CHECKPOINT_INTEGRITY_SYNTAX_AND_LINKS'


def digest(name):
    return sha256((ROOT/name).read_bytes()).hexdigest()


def check_link(document, target):
    path = (ROOT/document).parent/unquote(target.split('#', 1)[0])
    assert path.exists() or path.resolve() == (ROOT/AUDIT).resolve(), (document, target)


def build():
    for name, expected in PREDECESSORS.items():
        assert digest(name) == expected, ('changed predecessor', name)
    bound, packages, artifacts = {}, [], []
    bindings = 0
    for stem, generation, replay, review in PACKAGES:
        names = [f'scripts/{stem}.py', f'results/{stem}.json', f'results/{stem}.md']
        data = json.loads((ROOT/names[1]).read_text())
        construction_status = {
            'fresh_wall59_bond_replacement_20261018': 'BOUNDED_NUMERICAL_WALL59_BOND_REPLACEMENT',
            'fresh_wall59_low_bond_companions_20261018': 'BOUNDED_LOW_BOND_PENTAGON_COMPANIONS',
        }
        if stem in construction_status:
            assert data['status'] == construction_status[stem], (stem, data['status'])
        else:
            assert data['status'].startswith('PASS'), (stem, data['status'])
        for name, expected in data['sources'].items():
            assert digest(name) == expected, ('binding', stem, name)
            assert name not in bound or bound[name] == expected
            bound[name] = expected
        bindings += len(data['sources'])
        artifacts.extend(names)
        packages.append(dict(package=stem, status=data['status'],
            direct_bindings=len(data['sources']),
            generation=dict(identifier=generation, exit_code=0, status='PASS'),
            fresh_replay=dict(identifier=replay, exit_code=0, status='PASS'),
            review_scope=review))
    artifacts += HELPERS + list(PREDECESSORS) + [CHECKPOINT, str(Path(__file__).relative_to(ROOT))]
    artifacts = list(dict.fromkeys(artifacts))
    asts = [name for name in artifacts if name.endswith('.py')]
    for name in asts:
        ast.parse((ROOT/name).read_text(), filename=name)
    docs = [f'results/{stem}.md' for stem, *_ in PACKAGES] + [CHECKPOINT] + MUTABLE
    links = []
    for name in docs:
        for raw in re.findall(r'\]\(([^\n]+?)\)', (ROOT/name).read_text()):
            target = raw.split(' "', 1)[0].strip('<>')
            if target.startswith(('http:', 'https:', 'mailto:', '#')):
                continue
            check_link(name, target)
            links.append(dict(document=name, target=target))
    return dict(status=STATUS, date='2026-10-18', checkpoint=CHECKPOINT,
        completed_unrestricted_bound=62, new_dimension_bound=None, new_finite_volume_example=False,
        packages=packages,
        counts=dict(completed_packages=len(packages), direct_bindings=bindings,
                    distinct_bound_files=len(bound), python_asts=len(asts),
                    documents=len(docs), local_links=len(links)),
        artifact_sha256={name: digest(name) for name in artifacts},
        bound_file_sha256=bound,
        mutable_document_snapshots={name: digest(name) for name in MUTABLE},
        python_ast_files=asts, local_links=links,
        scope='Integrity, syntax and recorded local links only. Run identifiers record execution evidence, not reexecution by this audit. No additional mathematical, geometric, numerical or proof-assistant verification. Current-document hashes are historical snapshots.')


def verify(data):
    assert data['status'] == STATUS
    for group in ('artifact_sha256', 'bound_file_sha256'):
        for name, expected in data[group].items():
            assert digest(name) == expected, (group, name)
    for name in data['python_ast_files']:
        ast.parse((ROOT/name).read_text(), filename=name)
    for row in data['local_links']:
        check_link(row['document'], row['target'])


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
