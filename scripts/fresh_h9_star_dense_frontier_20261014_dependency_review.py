#!/usr/bin/env python3
"""Scope audit for raising ten noncore own contacts in five frozen proofs.

This binds prior independent arithmetic instead of redoing its polynomial
banks. The new checks concern literal core equality and monotone dependencies.
"""
import argparse,json
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_dense_frontier_20261014_dependency_review'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
GATES='results/fresh_h9_star_zero_20261014_gates.json'
NAMES=['fresh_h9_star_zero_20261014_own_three','fresh_h9_star_zero_20261014_opposite_single',
       'fresh_h9_star_lastblock_20261014','fresh_h9_star_shared_high_20261014','fresh_h9_star_opposite_high_20261014']
S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for f,h in d.get('sources',{}).items():
            assert sha(f)==h,f
            if f.endswith('.json'):read(f)
        return d
    geo,gates=read(GEOM),read(GATES)
    g=geo['representatives'][0]['geometry'];V=[set(v['incident'])for v in g['vertices']]
    core=set(S+[11,14,17]);assert len(core)==11 and any(set(S)<=v for v in V)
    assert all(any(set(pair)<=v for v in V)for pair in combinations(sorted(core),2))
    common={tuple(sorted((leaf,a))):3 for i in range(3)for leaf in(P[i],C[i])for a in S[2*i:2*i+2]}
    fixed={(3,14):3,(4,14):3};free={p:m for p,m in common.items()if p not in fixed}
    assert len(free)==10 and all(not set(p)<=core for p in free)
    assert all(set(p)<=core for p in fixed)
    by={f['id']:f for f in gates['families']}
    paths=[GEOM,GATES];proofs=[]
    def literal_scan(obj,counts):
        if isinstance(obj,dict):
            if 'actual_vertex'in obj:
                vi=obj['actual_vertex'];assert 0<=vi<len(V);counts['vertex_witness_occurrences']+=1
                if 'support'in obj:assert set(obj['support'])<=V[vi]
                if 'pair'in obj:assert set(obj['pair'])<=V[vi]
            for value in obj.values():literal_scan(value,counts)
        elif isinstance(obj,list):
            for x in obj:literal_scan(x,counts)
    for name in NAMES:
        p=f'results/{name}.json';review=f'results/{name}_review.json'
        d,r=read(p),read(review)
        assert d['status'].startswith('PASS_')and r['status'].startswith('PASS_')
        assert r['sources'][p]==sha(p)
        paths.extend([p,review,f'scripts/{name}.py'])
        ownrows=d.get('all_own_contacts',d.get('fixed_other_own_contacts',d.get('fixed_first_three_own',d.get('fixed_own_contacts'))))
        own={(a,b):m for a,b,m in ownrows}
        assert all(own[p]==3 for p in common)
        high=d.get('arbitrary_high_block_reduction',d.get('first_three_high_reductions'))
        bound_rows=[]
        for h in high:
            block=h.get('high_block')
            i=block if isinstance(block,int)else next(j for j in range(4)if S[2*j:2*j+2]==block)
            lab=own|{tuple(S[2*j:2*j+2]):6 if j==i else 4 for j in range(4)}
            if name.endswith('opposite_single'):lab|={(11,13):4,(12,17):4}
            if name.endswith(('shared_high','opposite_high')):lab[(12,13)]=7
            gate=by[h['reviewed_gate']];assert not gate['fixed']
            assert all(lab.get((a,b),0)>=m for a,b,m in gate['bounds'])
            bound_rows.append(dict(high_block=i,family=gate['id'],bounds=gate['bounds']))
        counts={'vertex_witness_occurrences':0};literal_scan(d,counts)
        rank_cases=0
        if name.endswith('own_three'):
            for record in d['records']:
                assert len(record['block_labels'])==4
                for row in record['cases']:
                    assert len(row['patterns'])==3 and all(0<=i<7 for i in row['patterns'])
                    assert len(row['contact_labels'])==2;rank_cases+=1
        elif name.endswith('lastblock_20261014'):
            for record in d['records']:
                assert len(record['block_labels'])==4
                assert len(record['anchor_profiles'])==2 and all(len(p)==2 for p in record['anchor_profiles'])
                assert len(record['pattern_order'])==3
                assert all(e is None or e in range(8)for row in record['pattern_order']for e in row)
                rank_cases+=len(record['determinants'])
        else:
            for row in d['cases']:
                if 'foreign_contacts'in row:
                    assert len(row['foreign_contacts'])==3 and all(e is None or e in range(6)for e in row['foreign_contacts'])
                else:assert row['foreign_third']is None or row['foreign_third']in range(2,8)
                assert len(row['third_labels'])==2;rank_cases+=1
        proofs.append(dict(name=name,rank_cases_inherited=rank_cases,high_bound_dependencies=bound_rows,**counts))
    assert [p['rank_cases_inherited']for p in proofs]==[21952,896,77952,7840,47040]
    # The literal first-state guards inherited by own_three are also
    # actual-vertex failures, not additional exact-core rank assumptions.
    first=read('results/fresh_h9_star_cube_metric_20261014.json')
    firstguards=first['forced_right_witnesses']+first['third_contact_upper_witnesses']
    assert len(firstguards)==30
    for row in firstguards:assert set(row['pair'])<=V[row['actual_vertex']]
    # Inspect the existing all-label composition's remaining high-last
    # reductions. They also have only lower-bound premises.
    composition=read('results/fresh_h9_star_lastblock_complete_20261014.json')
    for row in composition['strong_own_threshold_cases']+composition['dense_own_cases']:
        gate=by[row['family']];assert not gate['fixed']and row['bounds']==gate['bounds']
        assert all((a,b)not in free or m<=3 for a,b,m in gate['bounds'])
    paths+=['results/fresh_h9_star_cube_metric_20261014.json','results/fresh_h9_star_lastblock_complete_20261014.json',f'scripts/{TAG}.py']
    return dict(status='PASS_FIVE_FAMILY_NONCORE_MONOTONE_DEPENDENCY_REVIEW',
        sources={p:sha(p)for p in paths},recursive_manifests=len(seen),rank_core=sorted(core),
        exact_own_contacts=sorted([list(p)+[m]for p,m in fixed.items()]),
        noncore_own_lower_bounds=sorted([list(p)+[m]for p,m in free.items()]),
        internal_lower_bounds=[S[2*i:2*i+2]+[4]for i in range(4)],proof_dependencies=proofs,
        inherited_first_state_vertex_guards=30,composition_bound_only_guards=6,
        proof='The ten raised own entries are outside the eleven-root residual rank core. All rank formulas therefore remain literal equalities. Every use of these entries outside the core is either a bound-only reviewed gate or a nonpositive actual-vertex Z-matrix test, preserved by Perron monotonicity. Projection budgets concern core leaves or provide lower bounds only. The explicit internal>=4 assumptions replace the former missing-triple inference where required.',
        scope='Source/proof and literal-dependency review, not a new155680-case arithmetic replay. Inherits five hash-bound independent arithmetic reviews. Certifies widening only the ten listed own contacts toarbitrary finite labels>=3, with14own3/3 exact and four internals>=4. Last internal2/3, lowering an unused own contact toright, or changing14own contacts remain outside.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['proof_dependencies']))
