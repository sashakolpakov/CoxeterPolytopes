#!/usr/bin/env python3
"""Fresh literal transfers of the independently proved negative-plane lemma.

No incidence automorphism is assumed. The generic scalar argument is a
frozen independent dependency; every new support is located directly.
"""
import argparse
import json
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_t8_two_high_transfers_20261014_review'
GEOMETRY = 'results/fresh_h9_tree_cube_wedges_20261014.json'
THEOREM = 'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.json'
PRIMARY = 'results/fresh_h9_t8_sparse_search_20261014_bridge_variants.json'

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    geometry = json.loads((ROOT/GEOMETRY).read_text())
    theorem = json.loads((ROOT/THEOREM).read_text())
    assert theorem['status'] == 'PASS_INDEPENDENT_TWO_HIGH_BLOCKS_DOMAIN_AND_INDEX_REVIEW'
    bound = set()
    def bind(p):
        if p in bound: return
        bound.add(p)
        for target,expected in json.loads((ROOT/p).read_text()).get('sources',{}).items():
            assert digest(target) == expected, target
            if target.endswith('.json'): bind(target)
    bind(GEOMETRY)
    bind(THEOREM)
    bind(PRIMARY)
    primary = json.loads((ROOT/PRIMARY).read_text())
    assert primary['status'] == 'PASS_FOUR_LITERAL_TWO_HIGH_BRIDGE_PROOF_TRANSFERS'
    variants = {(r['bridge'],r['anchor']):r for r in primary['variants']}
    assert set(variants) == {(b,a) for b in (14,15) for a in (11,17)}
    row = geometry['representatives'][8]
    assert row['id'] == 8
    g = row['geometry']
    vertices = [set(v['incident']) for v in g['vertices']]
    missing = {frozenset(v) for v in g['minimal_nonfaces']}
    S = [1,2,4,5,7,8,12,13]
    def face(s):
        ids = sorted(s)
        vi = next((i for i,v in enumerate(vertices) if set(ids) <= v),None)
        assert vi is not None, ids
        return dict(support=ids, actual_vertex=vi, incident=sorted(vertices[vi]))
    common = face(S)
    rows = []
    for bridge,anchor in product((14,15),(11,17)):
        high = ((1,2),(4,5) if bridge == 14 else (7,8))
        low = ((7,8) if bridge == 14 else (4,5),(12,13))
        core = S+[anchor,bridge]
        meeting = [face(pair) for pair in combinations(core,2)]
        assert len(meeting) == 45
        triples = [sorted(list(block)+[bridge]) for block in high]+[sorted(list(low[1])+[anchor])]
        assert all(frozenset(t) in missing for t in triples)
        cross = []
        blocks = high+low
        for i,j in combinations(range(4),2):
            for a,b in product(blocks[i],blocks[j]):
                if i < 2:
                    support = sorted(set(blocks[i])|{b})
                    witness = 'strict (>=7,>=3) path, relaxed to affine(6,3)'
                else:
                    support = sorted(set(blocks[i])|set(blocks[j]))
                    witness = 'affine(4,3,4) path after testing cross label>=3'
                cross.append(dict(pair=[a,b],witness=witness,**face(support)))
        assert len(cross) == 24
        budgets = [face([a,b]+list(low[0])+list(low[1])+[bridge]) for a,b in product(*high)]
        anchors = [face(list(high[0])+list(high[1])+list(low[0])+[r,anchor]) for r in low[1]]
        stars = [face([anchor,bridge,a,b]) for a,b in product(*high)]
        assert len(budgets) == len(stars) == 4 and len(anchors) == 2
        lower = sorted([list(b)+[7] for b in high]+[list(b)+[4] for b in low])
        variant = variants[bridge,anchor]
        assert sorted(variant['finite_label_lower_bounds']) == lower
        assert variant['fixed_labels'] == [] and sorted(variant['positive_support']) == S
        assert sorted(variant['missing_triples']) == sorted(triples)
        for guards in variant['incidence_checks'].values():
            for guard in guards:
                assert set(guard['support']) <= vertices[guard['actual_vertex']]
        rows.append(dict(bridge=bridge,anchor=anchor,meeting_pairs=meeting,
                         finite_label_lower_bounds=lower,
                         missing_triples=triples,cross_guards=cross,
                         third_budgets=budgets,anchor_budgets=anchors,affine_stars=stars))
    condition_sets = sorted({tuple(tuple(x) for x in r['finite_label_lower_bounds']) for r in rows})
    assert [[list(x) for x in row] for row in condition_sets] == primary['distinct_condition_sets']
    paths = [GEOMETRY,THEOREM,PRIMARY,f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_FOUR_LITERAL_TWO_HIGH_TRANSFERS',
        sources={p:digest(p) for p in paths},recursively_bound_manifests=sorted(bound),
        positive_support=common,distinct_condition_sets=[[list(x) for x in row] for row in condition_sets],
        transfers=rows,
        proof='The independently proved two-high negative-plane argument uses only the roles of two high blocks, two low blocks, bridge and anchor, and the displayed supports. Every support and missing triple is present literally for each new assignment. The same uniform Lanner cost>5/4 and disjoint projection argument give a negative residual plane; no other label is fixed and no global incidence automorphism is used.',
        scope='Only T8 with first plus second or third internal blocks>=7 and the other two>=4. Both bridges14/15 and both anchors11/17 are checked separately. No whole-incidence exclusion, new metric search or polynomial bank.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args = parser.parse_args()
    result = run()
    out = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(result['transfers']))
