#!/usr/bin/env python3
"""Check frozen checkpoint integrity, syntax and links; no mathematical replay."""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/ambient62_ideal_closure_and_current_search_20261018.md'
AUDIT = CHECKPOINT.replace('.md', '_audit.json')
STATUS = 'PASS_AMBIENT62_CURRENT_SEARCH_CHECKPOINT_INTEGRITY_SYNTAX_AND_LINKS'
PACKAGES = [
    ('fresh_order14_c_inner_endpoint_projection_transfer_20261018', '7aff88', '1ebc2b',
     'Root and peer full source/report plus inherited contact-domain read. Complete selected support, forced physical norm-two projection, positive inverse end blocks, arbitrary supporting roots and full short8/current-priority interfaces. No reviewer embedding, inverse or census replay.'),
    ('fresh_order14_b_terminal_residual_transfer_20261018', '486d85', '3c04d0',
     'Root and peer full source/report plus inherited contact-domain read. Direct projection cases and doubled-terminal residual Cauchy cases, all extra contacts, full selected/nonselected and current-priority interfaces. No reviewer map, inverse or census replay.'),
    ('fresh_order14_ambient62_ideal_margin_closure_20261018', '3e7486', '0ca6d1',
     'Root and two peers full mathematical composition reads. Independent short8 priority, B/C14 selected patch, complete low9-only domain, exact rank13 savings, positive normalized and absolute margins and all-ideal higher tail. No reviewer convolution, scalar, DP or geometric census replay.'),
    ('fresh_finite_d17_five_unit_reduction_20261018', '7afe75', 'af8415',
     'Root and two peers full source/report and inherited selected/literal/coupon interfaces. Radical whole-component sealing, all current literal debits, distinct retained ideal halves, whole-small/D12 coupons and B14/C14 compatibility. No reviewer map, inertia, literal arithmetic, census or DP replay. The third review arrived after the local packet froze and is recorded here.'),
    ('fresh_finite_d17_cost120_rank62_budget_20261018', '5fe3ab', 'af30d9',
     'Root and peer full source/report plus inherited count engine. Current actual D13/E6 rule, simultaneous representative domination, complete large-count domain, residual upper/attainment certificates, exact438 witness and representative/type-tie distinction. No reviewer DP or arithmetic replay.'),
    ('fresh_finite_rank15_long_pair_support_scout_20261018', '63d147', 'd8cd33',
     'Root full source/report read within bounded diagnostic scope: complete536-pattern native-map domains for each rank15 type, all216 current literal records after D17=120, baseline-unit guard, invariant long orbits and separate hypothetical orbit trials. No new local theorem inferred and no reviewer map or arithmetic replay.'),
    ('fresh_rank16_existing_payment_exchange_20261018', '9b0cdf', '72c267',
     'Root and peer full source/report plus inherited contact-domain, sparse and literal interfaces. Six whole-D16 Lorentz seals, eight radical cases, whole-small/D13 halves and restored A17 halves, complete saved maps, all216 unchanged literal obligations and current-priority/nonselected/deep11 compatibility. No reviewer map, inertia, arithmetic or census replay.'),
    ('fresh_finite_rank15_existing_half_transfer_20261018', 'dcb9e6', '1182bb',
     'Root and peer full source/report plus scout and inherited contact-domain reads. Whole A2/A13 endpoint half, all remaining D15 contacts and retained D16 half, simultaneous216 literal debits including rank16, complete1691/83/304/536 partition and current ideal/nonselected routes. No reviewer arithmetic, map or census replay.'),
    ('fresh_finite_rank15_rank16_rank62_budget_20261018', 'c2ed7d', '746c65',
     'Root and two peers full source/report mathematical reads. Independent current base/cross formulas, conservative D16 envelope, separate A17 restoration, all35721 simultaneous type comparisons, complete1227 fixed vectors, all66258 unchanged residual items, exact fixed shifts and actual434 attainment. The921 complete Bellman tables are inherited without a DP rerun; reviewers did not replay arithmetic.'),
    ('fresh_rank15_rank16_dimension61_20261018', 'cd6064', 'caf2c7',
     'Author plus three peer full source/report composition reads. Independent actual ideal priority and short8/deep-zero routes, mutually retained rank15/rank16 alternatives, disjoint selected cases and simultaneous216 literal debits, counted A17 restoration, finite434/zero surplus, unchanged ideal absolute73, finite curvature434, single4805 extremum correction and completed bound62 higher tail. No reviewer census, DP or scalar convolution replay.'),
    ('fresh_cross_vertex_current_full_model_20261018', '872707', '749374',
     'Root author and peer full producer/helper/report reads, including inherited primal recovery. All-owner column costs, positive primitive scaling, monotone exact dual repair, full matrix replay and limited independent saved-helper scope. Preparation9636da, saved helper3389d5 and small display check4caa94 passed. Reviewers did not repeat the global optimization or matrix replay.'),
    ('fresh_wall66_seam_subdivision_20261018', 'c1bdab', '29cdba',
     'Root full producer/helper/report reads and separate peer disk-topology review. Saved helper8d5946 checks exact old199-support debit, complete representative-selection and low-recipe domains, saved native values, all2685 final failures and cone/PD ledgers. No independent physical-normal reconstruction, cone replay or interval signs.'),
]
HELPERS = ['scripts/check_cross_vertex_current_model_saved_20261018.py',
           'scripts/check_wall66_seam_subdivision_saved_20261018.py']
PREDECESSORS = {
    'results/dimension62_three_transfers_raw39_and_wall59_20261018.md':
        '173f60ced4bcdd62ff533279acc322fce7c3fa7050f05bb7d636ca8197b8bacb',
    'results/dimension62_three_transfers_raw39_and_wall59_20261018_audit.json':
        '10c78d4fffa2ab33f66eb1ce260d88e042d165756a67c1735068df52efd8cba8',
}
MUTABLE = ['README.md', 'RESEARCH_PLAN.md']


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
        if stem == 'fresh_wall66_seam_subdivision_20261018':
            assert data['status'] == 'BOUNDED_CHANGED_WALL66_SEAM_SUBDIVISION'
        else:
            assert data['status'].startswith('PASS'), (stem, data['status'])
        sources = data['prepared']['sources'] if 'prepared' in data else data['sources']
        for name, expected in sources.items():
            assert digest(name) == expected, ('binding', stem, name)
            assert name not in bound or bound[name] == expected
            bound[name] = expected
        bindings += len(sources)
        artifacts.extend(names)
        packages.append(dict(package=stem, status=data['status'],
            direct_bindings=len(sources),
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
        completed_unrestricted_bound=61, all_ideal_bound=61,
        new_finite_volume_example=False, packages=packages,
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
