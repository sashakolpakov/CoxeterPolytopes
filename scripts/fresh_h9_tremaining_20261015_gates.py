#!/usr/bin/env python3
"""New literal three-/two-high block binders for unresolved T representatives.

The frozen initial binder is reused as a read-only function. Each new metric
map checks its complete literal face premises; no target isomorphism is assumed.
"""
import argparse
import json
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
from fresh_h9_tnext_20261015_gates import run as initial_bind

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_tremaining_20261015_gates'
LEDGER = 'results/fresh_h9_tnext_20261015_ledger.json'
GEOMETRY = 'results/fresh_h9_tree_cube_wedges_20261014.json'
INITIAL_REVIEW = 'results/fresh_h9_tnext_20261015_scope_review.json'
HIGH = 'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks.json'
HIGH_REVIEW = 'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.json'
PREVIOUS_TRANSFER_REVIEW = 'results/fresh_h9_tnext_20261015_case31_transfers_review.json'
ENGINE = 'scripts/fresh_h9_tnext_20261015_search.py'


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def run(case):
    seen = set()
    def read(p):
        d = json.loads((ROOT / p).read_text())
        if p not in seen:
            seen.add(p)
            for q, h in d.get('sources', {}).items():
                assert sha(q) == h, q
                if q.endswith('.json'):
                    read(q)
        return d
    ledger, geo, theorem, review = map(read, (LEDGER, GEOMETRY, HIGH, HIGH_REVIEW))
    assert case in ledger['remaining_representatives']
    assert read(INITIAL_REVIEW)['status'].startswith('PASS_')
    assert read(PREVIOUS_TRANSFER_REVIEW)['status'].startswith('PASS_')
    assert theorem['status'] == 'PASS_T8_TWO_HIGH_BLOCKS_ALL_OTHER_LABELS_FREE_EXCLUDED'
    assert review['sources'][HIGH] == sha(HIGH) and review['status'].startswith('PASS_')
    initial = initial_bind(case)
    g = geo['representatives'][case]['geometry']
    vertices = {tuple(v['incident']) for v in g['vertices']}
    missing = {frozenset(m) for m in g['minimal_nonfaces']}
    blocks = [b for b in g['blocks'] if len(b) == 2]
    assert len(blocks) == 4
    outside = sorted(set(range(18)) - set(sum(blocks, [])))
    V = [set(v['incident']) for v in g['vertices']]
    # All within-block root flips really preserve this full nerve. They do
    # not change internal-label predicates, so one root ordering suffices.
    for flags in product((False, True), repeat=4):
        pm = list(range(18))
        for take, (a, b) in zip(flags, blocks):
            if take:
                pm[a], pm[b] = b, a
        assert {tuple(sorted(pm[i] for i in v)) for v in vertices} == vertices
        assert {frozenset(pm[i] for i in m) for m in missing} == missing
    supports = [('positive_support', theorem['positive_support'])]
    supports += [('meeting_pair', r['pair']) for r in theorem['all45_core_meeting_pairs']]
    for key in ('cross_right_guards', 'anchor_high_right_budgets',
                'third_other_blocks_right_budgets', 'third_anchor_upper_stars'):
        supports += [(key, r['support']) for r in theorem[key]]
    assert len(supports) == 80
    def actual(s):
        return next((i for i, v in enumerate(V) if set(s) <= v), None)
    maps = []
    for high_pair in combinations(range(4), 2):
        low = sorted(set(range(4)) - set(high_pair))
        for high in permutations(high_pair):
            for own in low:
                order = [*high, next(i for i in low if i != own), own]
                for bridge in outside:
                    if not all(frozenset(blocks[i] + [bridge]) in missing for i in high):
                        continue
                    for anchor in outside:
                        if anchor == bridge or frozenset(blocks[own] + [anchor]) not in missing:
                            continue
                        mp = dict(zip(theorem['positive_support'], sum((blocks[i] for i in order), [])))
                        mp.update({14: bridge, 11: anchor})
                        assert len(set(mp.values())) == len(mp) == 10
                        checks = []
                        for kind, source in supports:
                            target = [mp[x] for x in source]
                            vi = actual(target)
                            if vi is None:
                                break
                            checks.append(dict(kind=kind, source=source, target=target, actual_vertex=vi))
                        else:
                            triples = [sorted(mp[x] for x in s) for s in theorem['missing_triples']]
                            assert all(frozenset(s) in missing for s in triples)
                            bounds = sorted(sorted((mp[a], mp[b])) + [m]
                                            for a, b, m in theorem['finite_label_lower_bounds'])
                            assert not theorem['fixed_labels']
                            maps.append(dict(mapping=sorted(map(list, mp.items())),
                                             high_block_order=list(high), anchor_block=own,
                                             bridge=bridge, anchor=anchor, supports=checks,
                                             missing_triples=triples, bounds=bounds))
    families = initial['families'][:]
    keys = {(tuple(map(tuple, f['fixed'])), tuple(sorted(map(tuple, f['bounds'])))) for f in families}
    new = []
    for i, row in enumerate(maps):
        key = ((), tuple(map(tuple, row['bounds'])))
        if key not in keys:
            keys.add(key)
            f = dict(id=f'two_high_metric_{len(new)}',origin='reviewed_literal_two_high_negative_plane',
                     fixed=[],bounds=row['bounds'],witness_map=i)
            new.append(f);families.append(f)
    assert len(keys) == len(families)
    sources = dict(initial['sources'])
    paths = [LEDGER, INITIAL_REVIEW, HIGH, HIGH_REVIEW, PREVIOUS_TRANSFER_REVIEW,
             ENGINE, f'scripts/{TAG}.py']
    sources.update({p:sha(p) for p in paths})
    return dict(status='PASS_TNEXT_LITERAL_REVIEWED_CONDITIONS', package='tremaining', case=case,
                sources=sources, vertices=initial['vertices'],
                missing_triangles=initial['missing_triangles'], dotted_pairs=initial['dotted_pairs'],
                inherited_predicates=initial['inherited_predicates'],
                three_large_predicates=initial['new_predicates'], two_high_predicates=len(new),
                total_predicates=len(families), families=families,
                initial_literal_binder=initial,
                two_high_maps=maps, two_high_new_conditions=new,
                literal_two_high_support_occurrences=80*len(maps),
                simultaneous_twin_swaps=16,
                scope='All inherited local predicates, all four candidate three-block conditions and '
                      'all two-high block-order/bridge/anchor maps of the named theorem are considered '
                      'on this literal target. Each retained map checks every proof premise. Within-pair '
                      'orderings are reduced only by verified full-nerve twin automorphisms. '
                      'No whole-incidence exclusion is made by this binder and no other metric theorem '
                      'or arbitrary local embedding is claimed to have been exhaustively searched.')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',type=int,required=True)
    p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run(a.case);out=ROOT/f'results/{TAG}_case_{a.case:02d}.json'
    if a.verify:
        assert d==json.loads(out.read_text())
    else:
        out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],a.case,d['total_predicates'],d['three_large_predicates'],d['two_high_predicates'])
