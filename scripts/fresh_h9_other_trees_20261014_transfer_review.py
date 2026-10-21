#!/usr/bin/env python3
"""Independent literal-premise audit for local theorem transfers, no imports."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_other_trees_20261014_transfer_review'
FILES=['results/fresh_h9_other_trees_20261014_transfers.json',
       'results/fresh_h9_tree_cube_wedges_20261014.json',
       'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks.json',
       'results/fresh_h9_t8_asymmetric_20261014.json',
       'results/fresh_h9_t8_sparse_frontier_20261014.json',
       'results/fresh_h9_t8_sparse_search_20261014_all_blocks.json',
       'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.json',
       'results/fresh_h9_t8_all_blocks_20261014_review.json',
       'results/fresh_h9_t8_sparse_frontier_20261014_review.json']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    docs=[json.loads((ROOT/p).read_text())for p in FILES]
    for d in docs:
        for p,h in d['sources'].items():assert sha(p)==h,p
    transfer,geometry,high,domain,sparse,allblocks,*reviews=docs
    assert all('PASS'in d['status']for d in reviews)
    hp=[('positive_support',high['positive_support'])]
    hp +=[('meeting_pair',list(p))for p in combinations(high['core_roots'],2)] if 'core_roots'in high else [('meeting_pair',r['pair'])for r in high['all45_core_meeting_pairs']]
    assert sum(k=='meeting_pair'for k,s in hp)==45
    for key in ['cross_right_guards','anchor_high_right_budgets','third_other_blocks_right_budgets','third_anchor_upper_stars']:
        hp +=[(key,r['support'])for r in high[key]]
    sp=[('positive_support',domain['positive_support'])]
    sp +=[('meeting_pair',list(p))for p in combinations(sparse['core_roots'],2)]
    assert sum(k=='meeting_pair'for k,s in sp)==55
    for key in ['cross_right_guards','anchor_budgets','third_budgets']:
        sp +=[(key,r['support'])for r in domain[key]]
    for pattern in domain['own_support_patterns']:
        guard=pattern['anchor_pair_guard']
        if 'support'in guard:sp.append(('anchor_pair_guard',guard['support']))
        else:
            assert guard['domain']=='any finite label; delta in[0,1)'
            sp.append(('anchor_pair_guard',guard['pair']))
    sp +=[('third_upper_guard',r['support'])for r in sparse['replacement_third_contact_guards']]
    counts=Counter();checks=0;rows=[]
    for case in transfer['cases']:
        cid=case['case'];g=geometry['representatives'][cid]['geometry']
        missing={frozenset(m)for m in g['minimal_nonfaces']}
        vertices={frozenset(v)for v in combinations(range(18),9)if not any(m<=set(v)for m in missing)}
        assert vertices=={frozenset(v['incident'])for v in g['vertices']}
        assert len(vertices)==case['vertices']
        for rec in case['maps']:
            mp=dict(rec['mapping']);assert len(mp)==len(rec['mapping'])==len(set(mp.values()))
            kind=rec['kind'];is_high=kind=='two_high_negative_plane'
            assert is_high or kind=='sparse_all_blocks'
            premises=hp if is_high else sp;proof=high if is_high else allblocks
            assert len(premises)==len(rec['actual_support_checks'])
            for (name,original),saved in zip(premises,rec['actual_support_checks']):
                mapped=[mp[x]for x in original]
                assert saved['kind']==name and saved['source_support']==original and saved['target_support']==mapped
                target=set(mapped)
                assert len(target)==len(original) and not any(m<=target for m in missing)
                assert any(target<=v for v in vertices)
                assert target<=set(g['vertices'][saved['actual_vertex']]['incident'])
                checks+=1
            triangles=high['missing_triples'] if is_high else [[11,12,13],[12,13,17]]
            expected=[sorted(mp[x]for x in t)for t in triangles]
            assert expected==rec['required_missing_triples'] and all(frozenset(t)in missing for t in expected)
            def mapped(rows):return sorted([sorted([mp[a],mp[b]])+[m]for a,b,m in rows])
            assert mapped(proof['fixed_labels'])==rec['fixed_labels']
            assert mapped(proof['finite_label_lower_bounds'])==rec['finite_label_lower_bounds']
            counts[kind]+=1
        keys=set()
        for family in case['families']:
            witness=case['maps'][family['witness_map']]
            assert family['kind']==witness['kind']
            assert family['fixed']==witness['fixed_labels'] and family['bounds']==witness['finite_label_lower_bounds']
            key=(tuple(map(tuple,family['fixed'])),tuple(map(tuple,family['bounds'])))
            assert key not in keys;keys.add(key)
        rows.append(dict(case=cid,vertices=len(vertices),maps=len(case['maps']),distinct_conditions=len(case['families'])))
    assert dict(counts)=={'two_high_negative_plane':17,'sparse_all_blocks':12}
    paths=FILES+[f'scripts/{TAG}.py','scripts/fresh_h9_other_trees_20261014_transfers.py']
    return dict(status='PASS_INDEPENDENT_LOCAL_PREMISE_TRANSFERS_TO_FOUR_TREES',
        sources={p:sha(p)for p in paths},cases=rows,map_counts=dict(counts),literal_face_checks=checks,
        scope='Separate incidence reconstruction, injective maps, all proof faces/missing triples and every applied condition set. The analytic polynomial and negative-plane proofs are inherited from reviewed frozen sources. No isomorphism assumption and no new exhaustive map classification. Whole-incidence exclusions require the separate complete searches.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['map_counts'],d['literal_face_checks'],flush=True)
