#!/usr/bin/env python3
"""Literal existing-gate binder and three fresh generic transfers for T-case6."""
import argparse
import json
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_t6_frontier_20261015_gates'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
TARGETS='results/fresh_h9_dotted_central_20261014_next_targets.json'
THEOREM='results/fresh_h9_star_three_large_blocks_20261015.json'
REVIEW='results/fresh_h9_star_three_large_blocks_20261015_review.json'


def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for q,h in d.get('sources',{}).items():
            assert sha(q)==h,q
            if q.endswith('.json'):read(q)
        return d
    geo,targets,theorem,review=[read(p) for p in (GEOM,TARGETS,THEOREM,REVIEW)]
    assert targets['status']=='PASS_LITERAL_GATE_ELIGIBILITY_ONLY_NO_SEARCH'
    assert theorem['status']=='PASS_FOUR_THREE_LARGE_BLOCK_LORENTZ_OBSTRUCTIONS'
    assert review['status'].startswith('PASS_') and review['sources'][THEOREM]==sha(THEOREM)
    g=geo['representatives'][6]['geometry']
    V=[set(v['incident']) for v in g['vertices']]
    MN={frozenset(m) for m in g['minimal_nonfaces']}
    assert len(V)==488 and {len(m) for m in MN}=={2,3}
    assert {tuple(sorted(v)) for v in V}=={v for v in combinations(range(18),9) if not any(m<=set(v) for m in MN)}
    B=[b for b in g['blocks'] if len(b)==2]
    assert B==[[1,2],[6,7],[9,10],[12,13]]
    target=targets['cases'][6]
    assert target['case']==6 and target['total_unique_conditions']==155 and not target['closed']
    families=[];old_counts={};old_witnesses=0
    def actual(nodes):return next((i for i,v in enumerate(V) if set(nodes)<=v),None)
    for kind,rec in target['gates'].items():
        for r in rec['maps']:
            A,C,x,y=r['A'],r['B'],r['x'],r['y']
            assert len(set(A+C+[x,y]))==6
            assert {frozenset(m) for m in r['missing_own_triples']}=={frozenset(A+[x]),frozenset(C+[y])}
            assert all(frozenset(m) in MN for m in r['missing_own_triples'])
            expected=[A+C,C+[x],A+[y]]+[[x,y,a,b] for a in A for b in C]
            bounds=[]
            if kind=='dense_high':
                bounds=[sorted(A)+[6],sorted(C)+[7]]+[sorted((x,a))+[3] for a in A]
            elif kind=='strong_contact':
                hits=[b for b in C if sorted((y,b))+[4] in r['bounds']];assert len(hits)==1
                expected.append(A+[y,hits[0]])
                bounds=[sorted(A)+[4],sorted(C)+[6]]+[sorted((x,a))+[3] for a in A]+[sorted((y,hits[0]))+[4]]
            else:
                assert kind=='two_dense'
                expected += [[x,y,b] for b in C]
                bounds=[sorted(A)+[6],sorted(C)+[4]]+[sorted((x,a))+[3] for a in A]+[sorted((y,b))+[3] for b in C]
            assert sorted(map(sorted,expected))==sorted(s['support'] for s in r['supports'])
            assert sorted(bounds)==r['bounds']
            assert sorted(s['support'] for s in r['threshold_pair_witnesses'])==sorted(b[:2] for b in bounds)
            for s in r['supports']+r['threshold_pair_witnesses']:
                assert set(s['support'])<=V[s['actual_vertex']]
                old_witnesses+=1
        for i,r in enumerate(rec['conditions']):
            assert r['bounds']==rec['maps'][r['witness_map']]['bounds']
            families.append(dict(id=f'{kind}_{i}',origin='reviewed_literal_'+kind,fixed=[],bounds=r['bounds']))
        old_counts[kind]=len(rec['conditions'])
    assert old_counts==dict(dense_high=23,strong_contact=88,two_dense=44)
    leaves={i:[x for x in range(18) if frozenset(B[i]+[x]) in MN] for i in range(4)}
    assert all(x not in set(sum(B,[])) for choices in leaves.values() for x in choices)
    placements=[];uncovered=[]
    for chosen in combinations(range(4),3):
        cross=[]
        for i,j in combinations(chosen,2):
            for a in B[i]:
                for b in B[j]:
                    z=next(x for x in B[i] if x!=a)
                    vi=actual([z,a,b]);assert vi is not None
                    cross.append(dict(ordered_support=[z,a,b],support=sorted([z,a,b]),actual_vertex=vi,
                                      large_pair=sorted([z,a]),trial_nonright_pair=sorted([a,b])))
        successful=[]
        for L in product(*(leaves[i] for i in chosen)):
            if len(set(L))<3 or actual(L) is None:continue
            foreign=[]
            for k,i in enumerate(chosen):
                for j in chosen:
                    if i==j:continue
                    for a in B[j]:
                        z=next(x for x in B[j] if x!=a);vi=actual([z,a,L[k]])
                        if vi is None:break
                        foreign.append(dict(ordered_support=[z,a,L[k]],support=sorted([z,a,L[k]]),actual_vertex=vi,
                                            large_pair=sorted([z,a]),trial_nonright_pair=sorted([a,L[k]])))
            if len(foreign)==12:successful.append((L,foreign))
        if not successful:
            uncovered.append(list(chosen));continue
        L,foreign=successful[0]
        own=[sorted(B[i]+[L[k]]) for k,i in enumerate(chosen)]
        assert all(frozenset(m) in MN for m in own)
        assert all(actual(p) is not None for m in own for p in combinations(m,2))
        f=dict(id='three_large_blocks_'+'_'.join(map(str,chosen)),origin='literal_three_timelike_leaf_transfer',
               fixed=[],bounds=[B[i]+[6] for i in chosen])
        families.append(f)
        placements.append(dict(chosen_blocks=list(chosen),chosen_leaves=list(L),eligible_leaf_choices=len(successful),
                               cross_right_guards=cross,foreign_right_guards=foreign,own_Lanner_triples=own,
                               leaf_triangle=dict(support=sorted(L),actual_vertex=actual(L)),family=f))
    assert [r['chosen_blocks'] for r in placements]==[[0,1,2],[0,1,3],[1,2,3]]
    assert uncovered==[[0,2,3]]
    keys={tuple(tuple(sorted(map(tuple,f[k]))) for k in ('fixed','bounds')) for f in families}
    assert len(families)==len(keys)==158
    paths=[f'scripts/{TAG}.py',GEOM,TARGETS,THEOREM,REVIEW]
    return dict(status='PASS_T6_158_LITERAL_REVIEWED_CONDITIONS',case=6,sources={p:sha(p) for p in paths},
                vertices=488,missing_triangles=11,dotted_pairs=10,old_counts=old_counts,
                inherited_predicates=155,new_predicates=3,total_predicates=158,families=families,
                inherited_support_occurrences=old_witnesses,three_large_transfers=placements,
                uncovered_three_large_block_triples=uncovered,
                earlier_search=dict(status=target['earlier_status'],nodes=target['earlier_nodes']),
                scope='Necessary exact conditions for compact T-case6. The generic three-timelike-leaf proof is inherited; '
                      'all target premises are rechecked literally. Triple023 is not excluded by this transfer. '
                      'No search or whole-incidence exclusion in this binder.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['total_predicates'],d['uncovered_three_large_block_triples'])
