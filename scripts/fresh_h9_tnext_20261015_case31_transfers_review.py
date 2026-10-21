#!/usr/bin/env python3
"""Independent literal application review; imports neither transfer producer nor search."""
import argparse,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_tnext_20261015_case31_transfers_review'
def run():
    src='results/fresh_h9_tnext_20261015_case31_transfers.json';d=json.loads((ROOT/src).read_text());sha=lambda p:sha256((ROOT/p).read_bytes()).hexdigest()
    for p,h in d['sources'].items():assert sha(p)==h,p
    geom=json.loads((ROOT/'results/fresh_h9_tree_cube_wedges_20261014.json').read_text())['representatives'][31]['geometry']
    V=[set(v['incident'])for v in geom['vertices']];M={frozenset(x)for x in geom['minimal_nonfaces']}
    tp='results/fresh_h9_t8_sparse_search_20261014_two_high_blocks.json';t=json.loads((ROOT/tp).read_text())
    rp='results/fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.json';r=json.loads((ROOT/rp).read_text());assert r['sources'][tp]==sha(tp)
    supports=[('positive_support',t['positive_support'])]+[('meeting_pair',x['pair'])for x in t['all45_core_meeting_pairs']]
    for name in ('cross_right_guards','anchor_high_right_budgets','third_other_blocks_right_budgets','third_anchor_upper_stars'):
        supports.extend((name,x['support'])for x in t[name])
    assert len(supports)==80 and t['fixed_labels']==[]
    for record in d['maps']:
        mp=dict(record['mapping']);assert len(mp)==len(set(mp.values()))==10
        assert len(record['supports'])==80
        for (kind,S),w in zip(supports,record['supports']):
            target=[mp[i]for i in S]
            assert w['kind']==kind and w['source']==S and w['target']==target
            assert set(target)<=V[w['actual_vertex']] and any(set(target)<=v for v in V)
        assert {frozenset(mp[i]for i in x)for x in t['missing_triples']}<=M
        expected=sorted(sorted([mp[a],mp[b]])+[m]for a,b,m in t['finite_label_lower_bounds'])
        assert expected==record['bounds']
    for c in d['new_conditions']:
        assert c['fixed']==[] and c['bounds']==d['maps'][c['witness_map']]['bounds']
    state=json.loads((ROOT/'results/fresh_h9_tnext_20261015_search_case_31_batch_001.json').read_text())
    labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
    hits=[i for i,c in enumerate(d['new_conditions'])if all(labels[a,b]>=m for a,b,m in c['bounds'])]
    assert hits==d['saved_survivor_hits'] and hits
    files=[src,tp,rp,'results/fresh_h9_tree_cube_wedges_20261014.json',f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_LITERAL_T31_NEGATIVE_PLANE_TRANSFER_REVIEW',sources={p:sha(p)for p in files},maps=len(d['maps']),actual_support_premises=80*len(d['maps']),missing_triple_occurrences=3*len(d['maps']),new_conditions=len(d['new_conditions']),saved_survivor_hits=hits,scope='Checks each listed map of the already independently reviewed all-label negative-plane theorem. No claim of complete T31 exclusion or exhaustive transfer enumeration.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
