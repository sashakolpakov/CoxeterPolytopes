#!/usr/bin/env python3
"""Integrity, syntax and link audit for the finite464 research checkpoint.

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
CHECKPOINT = 'results/finite464_cycle11_and_shared_edge_20261018.md'
AUDIT = CHECKPOINT.replace('.md', '_audit.json')
PACKAGES = [
    ('fresh_ambient64_payment_frontier_20261018', '16c58d', '244610',
     'Full separate source/report review of current scalar domain, corrected2419-profile conditional domain and explicit finite coupon/Bellman formulas. No independent scalar/Bellman rerun.'),
    ('fresh_depth11_finite8_zero_transfer_20261018', '35f9a0', '29c166',
     'Two separate full source/report mathematical reviews of A27 corridor, native minimal-core bridge, cycle7/48 interface and actual-corner composition. Second review recorded in checkpoint after local report freeze; no independent144-value replay.'),
    ('fresh_finite_d16_d17_two_unit_exchange_20261018', 'c4c221', '7c4f9b',
     'Separate full source/report mathematical review of four radical cases, selected maps,216 literal interfaces and noncircular depth11 composition. No independent map or literal rerun.'),
    ('fresh_cycle11_half_pair_transfer_20261018', '95bef3', '97502e',
     'Two full source/report mathematical reviews, including inherited all-rank spectator inequalities, native profiles and triangle exclusion. No duplicate matrix or scalar rerun.'),
    ('fresh_ambient64_finite_closure_20261018', '96a342', 'f3ee6d',
     'Separate full binder/report composition review, including exact scenario match and conditional budget report. No inherited Bellman/geometric/ideal census replay. Finite464 complete; no63 claim.'),
    ('fresh_raw39_shared_edge_target_incidence_20261018', '64774', '52711',
     'Root and separate reviewer read complete source/report/helper. No-project-import helper91e5bc independently checks56 supports,10 frames,15 exclusions,29 cliques,23 vectors,46 inequalities,27 terms and2columns/repairs. No historical ownership or complete embedding recensus.'),
    ('fresh_split30_continuous_u_20261018', '6769', '4323',
     'Root full source/report/helper and inherited-interface review. Independent saved checkerdf01db checks topology/debit/saved native values and boundary ledgers; no physical normal/root/cone reconstruction or exactification.'),
]
HELPERS = ['scripts/check_raw39_shared_edge_target_saved_20261018.py',
           'scripts/check_split30_continuous_u_saved_20261018.py']
FROZEN_PREDECESSORS = {
    'results/dimension64_universal_corridors_and_changed_caps_20261018.md':
        'de8946b2ff0d4b2d410aa4f0a842c89ae7c558caf64b1de53785f46df33f4e8e',
    'results/dimension64_universal_corridors_and_changed_caps_20261018_audit.json':
        'c3f030b5179a8fde199558b4b33119a8bbf035c4ee5a1bed9a3aeff3a3e0cfe4',
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
        if stem == 'fresh_split30_continuous_u_20261018':
            assert data['status'] == 'BOUNDED_CONTINUOUS_U_NATIVE_FIRST_SCOUT'
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
    return dict(status='PASS_FINITE464_CHECKPOINT_INTEGRITY_SYNTAX_AND_LINKS',
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
    assert data['status'] == 'PASS_FINITE464_CHECKPOINT_INTEGRITY_SYNTAX_AND_LINKS'
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
