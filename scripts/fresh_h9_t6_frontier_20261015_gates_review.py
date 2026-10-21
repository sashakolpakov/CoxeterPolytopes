#!/usr/bin/env python3
"""Independent literal-premise/source audit only; does not execute a DFS."""
import argparse
import ast
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_t6_frontier_20261015_gates_review'
P='fresh_h9_t6_frontier_20261015'

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    paths=[f'scripts/{P}_gates.py',f'results/{P}_gates.json',f'scripts/{P}_independent.py',
           'results/fresh_h9_tree_cube_wedges_20261014.json',
           'results/fresh_h9_dotted_central_20261014_next_targets.json',
           'results/fresh_h9_star_three_large_blocks_20261015.json',
           'results/fresh_h9_star_three_large_blocks_20261015_review.json',
           'scripts/fresh_h9_star_large_block_resume_20261015_independent.py',
           'results/fresh_h9_star_closure_hostile_20261015.json',f'scripts/{TAG}.py']
    seen={}
    def load(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for q,h in d.get('sources',{}).items():
            assert digest(q)==h,q
            if q.endswith('.json'):load(q)
        return d
    for p in paths:
        if p.endswith('.json'):load(p)
    geometry=load(paths[3])['representatives'][6]['geometry']
    missing={frozenset(v) for v in geometry['minimal_nonfaces']}
    vertices={tuple(v['incident']) for v in geometry['vertices']}
    def face(v):return not any(m<=set(v) for m in missing)
    assert vertices=={v for v in combinations(range(18),9) if face(v)}
    assert len(vertices)==488
    assert sum(len(m)==3 for m in missing)==11 and sum(len(m)==2 for m in missing)==10
    def witness(row):
        s=set(row['support']);assert face(s)
        assert s<=set(geometry['vertices'][row['actual_vertex']]['incident'])
    target=load(paths[4])['cases'][6]
    binder=load(paths[1]);bounds_rebuilt=[];counts={};occurrences=0
    for name,item in target['gates'].items():
        for row in item['maps']:
            A,B,x,y=tuple(row['A']),tuple(row['B']),row['x'],row['y']
            assert len(set(A+B+(x,y)))==6
            assert {frozenset(A+(x,)),frozenset(B+(y,))}==set(map(frozenset,row['missing_own_triples']))
            assert all(frozenset(t) in missing for t in row['missing_own_triples'])
            premises=[A+B,B+(x,),A+(y,)]+[(x,y,a,b) for a in A for b in B]
            def threshold(u,v,k):return (*sorted((u,v)),k)
            if name=='dense_high':
                thresholds=[(*sorted(A),6),(*sorted(B),7)]+[threshold(x,a,3) for a in A]
            elif name=='strong_contact':
                choices=[b for b in B if list(threshold(y,b,4)) in row['bounds']]
                assert len(choices)==1
                b=choices[0];premises.append(A+(y,b))
                thresholds=[(*sorted(A),4),(*sorted(B),6)]+[threshold(x,a,3) for a in A]+[threshold(y,b,4)]
            else:
                assert name=='two_dense'
                premises.extend((x,y,b) for b in B)
                thresholds=[(*sorted(A),6),(*sorted(B),4)]+[threshold(x,a,3) for a in A]+[threshold(y,b,3) for b in B]
            assert sorted(tuple(sorted(p)) for p in premises)==sorted(tuple(w['support']) for w in row['supports'])
            assert sorted(thresholds)==sorted(map(tuple,row['bounds']))
            assert sorted(t[:2] for t in thresholds)==sorted(tuple(w['support']) for w in row['threshold_pair_witnesses'])
            for w in row['supports']+row['threshold_pair_witnesses']:witness(w);occurrences+=1
        counts[name]=len(item['conditions'])
        for cond in item['conditions']:
            assert cond['bounds']==item['maps'][cond['witness_map']]['bounds']
            bounds_rebuilt.append(tuple(map(tuple,cond['bounds'])))
    assert counts==dict(dense_high=23,strong_contact=88,two_dense=44)
    assert occurrences==2288==binder['inherited_support_occurrences']
    assert all(not f['fixed'] for f in binder['families'])
    assert bounds_rebuilt==[tuple(map(tuple,f['bounds'])) for f in binder['families'][:155]]
    blocks=[tuple(v) for v in geometry['blocks'] if len(v)==2]
    assert blocks==[(1,2),(6,7),(9,10),(12,13)]
    leaves=[[i for i in range(18) if frozenset(b+(i,)) in missing] for b in blocks]
    assert not set(sum(leaves,[]))&set(sum(blocks,()))
    all_choices={}
    for chosen in combinations(range(4),3):
        # Every cross-node and foreign-leaf pair must be forced right by
        # an actual three-node path with a >=6 internal edge.
        for i,j in combinations(chosen,2):
            for b in blocks[j]:assert face(blocks[i]+(b,))
            for a in blocks[i]:assert face(blocks[j]+(a,))
        valid=[]
        for ls in product(*(leaves[i] for i in chosen)):
            if len(set(ls))!=3 or not face(ls):continue
            if all(face(blocks[j]+(ls[k],)) for k,i in enumerate(chosen) for j in chosen if j!=i):valid.append(ls)
        all_choices[chosen]=valid
    assert all_choices=={(0,1,2):[(14,5,8),(14,5,17)],(0,1,3):[(14,5,11),(14,5,17)],(0,2,3):[],(1,2,3):[(15,16,11)]}
    for r in binder['three_large_transfers']:
        chosen=tuple(r['chosen_blocks']);ls=tuple(r['chosen_leaves'])
        assert ls==all_choices[chosen][0] and r['eligible_leaf_choices']==len(all_choices[chosen])
        assert r['family']['bounds']==[list(blocks[j])+[6] for j in chosen]
        assert sorted(map(frozenset,r['own_Lanner_triples']),key=lambda a:sorted(a))==sorted((frozenset(blocks[j]+(ls[k],)) for k,j in enumerate(chosen)),key=lambda a:sorted(a))
        assert len(r['cross_right_guards'])==len(r['foreign_right_guards'])==12
        for w in r['cross_right_guards']+r['foreign_right_guards']:
            witness(w);assert tuple(w['large_pair']) in blocks
            assert set(w['large_pair'])|set(w['trial_nonright_pair'])==set(w['support'])
            assert len(set(w['large_pair'])&set(w['trial_nonright_pair']))==1
        witness(r['leaf_triangle'])
    assert binder['uncovered_three_large_block_triples']==[[0,2,3]]
    assert len(binder['families'])==158
    assert len({tuple(map(tuple,f['bounds'])) for f in binder['families']})==158
    # The expensive recursion has exactly the previously audited implementation;
    # compare its nested functions structurally, permitting no changed pruning.
    old=ast.parse((ROOT/paths[7]).read_text());new=ast.parse((ROOT/paths[2]).read_text())
    names=['elliptic','can_prune_family','visit']
    # Find actual common nested function names rather than relying on prose.
    oldfunc={n.name:ast.dump(n,include_attributes=False) for n in ast.walk(old) if isinstance(n,ast.FunctionDef)}
    newfunc={n.name:ast.dump(n,include_attributes=False) for n in ast.walk(new) if isinstance(n,ast.FunctionDef)}
    checked=[n for n in oldfunc if n!='run']
    assert set(oldfunc)==set(newfunc)
    for n in checked:assert oldfunc[n]==newfunc[n],n
    # Independently verify the case-specific simultaneous symmetry guards.
    for choices in product((0,1),repeat=4):
        perm=list(range(18))
        for take,(a,b) in zip(choices,blocks):
            if take:perm[a],perm[b]=b,a
        assert {tuple(sorted(perm[i] for i in v)) for v in vertices}==vertices
        assert {frozenset(perm[i] for i in m) for m in missing}==missing
    return dict(status='PASS_T6_LITERAL_BINDER_AND_SOURCE_REVIEW',sources={p:digest(p) for p in paths},
                checked_recursive_manifests=len(seen),vertices=488,missing_triangles=11,dotted_pairs=10,
                inherited_counts=counts,inherited_witness_occurrences=occurrences,new_transfers=[dict(blocks=list(k),valid_leaves=[list(x) for x in v]) for k,v in all_choices.items()],
                inherited_predicates=155,new_predicates=3,total_predicates=158,
                unchanged_function_bodies=checked,simultaneous_symmetries=16,
                scope='Independent literal reconstruction and source review only. Analytic theorems and geometric realization are inherited. No DFS execution or new whole-case exclusion is claimed; the generic transfer does not cover block triple023.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['total_predicates'])
