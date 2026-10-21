#!/usr/bin/env python3
"""Independent literal dependency removal and star involution review."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_single_right_20261014'
FREE={frozenset(x)for x in[(6,15),(7,15),(9,16),(10,16)]}
PERM=[1,0,17,12,13,16,9,10,15,6,7,14,3,4,11,8,5,2]
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
    primary=f'results/{TAG}.json';d=read(primary)
    geo=read('results/fresh_h9_tree_cube_wedges_20261014.json')['representatives'][0]['geometry']
    V=[tuple(v['incident'])for v in geo['vertices']];assert all(v==tuple(sorted(v))for v in V)
    M={frozenset(m)for m in geo['minimal_nonfaces']};S=[3,4,6,7,9,10,12,13]
    core=set(S+[11,14,17]);assert core==set(d['rank_core'])and all(not p<=core for p in FREE)
    assert all(any(set(pair)<=set(v)for v in V)for pair in combinations(sorted(core),2))
    first=read('results/fresh_h9_star_cube_metric_20261014.json')
    counts=[];dependencies=[]
    for audit in d['five_family_audits']:
        name=audit['name'];proof=read(f'results/{name}.json');review=read(f'results/{name}_review.json')
        assert review['status'].startswith('PASS_')and review['sources'][f'results/{name}.json']==sha(f'results/{name}.json')
        originals={}
        if name.endswith('own_three'):
            assert proof['zero_rank_cases']==[]and not proof['unresolved']
            for i,row in enumerate(proof['foreign_budgets']):
                if row['leaf']in(11,14,17):originals[('foreign_budgets',i)]=row
            needed=[r for r in first['forced_right_witnesses']if set(r['pair'])<=core]+first['third_contact_upper_witnesses']
            assert len(needed)==27
            for i,row in enumerate(needed):originals[('inherited_needed_guard',i)]=row
            assert audit['unused_own_pair_guards']==[[2,14],[5,15],[8,16]]
            assert sorted(audit['unused_leaf_budgets'])==[2,5,8,15,16]
        else:
            # Stack traversal is independent of the producer's recursive scanner.
            stack=[((),proof)]
            while stack:
                path,obj=stack.pop()
                if isinstance(obj,dict):
                    if'actual_vertex'in obj:originals[path]=obj
                    stack.extend((path+(k,),v)for k,v in obj.items())
                elif isinstance(obj,list):stack.extend((path+(i,),v)for i,v in enumerate(obj))
        assert len(originals)==audit['needed_witness_occurrences']
        checked=set();summ=Counter()
        for row in audit['literal_witnesses']:
            path=tuple(row['path']);assert path in originals and path not in checked;checked.add(path)
            old=originals[path];vi=old['actual_vertex'];assert row['actual_vertex']==vi and 0<=vi<len(V)
            if'support'in old:ids=tuple(sorted(old['support']));mode='literal_support'
            else:ids=V[vi][:old.get('prefix_order',old.get('leading_order',9))];mode='leading_principal_prefix'
            assert list(ids)==row['support']and row['mode']==mode and set(ids)<=set(V[vi])
            assert all(not p<=set(ids)for p in FREE)
            summ[ids]+=1
        assert checked==set(originals)
        assert audit['unique_witness_supports']==[dict(support=list(k),occurrences=v)for k,v in sorted(summ.items())]
        counts.append(len(checked));dependencies.append(name)
    assert counts==[30,96,708,479,3643]and sum(counts)==4956
    gates={r['id']:r for r in read('results/fresh_h9_star_zero_20261014_gates.json')['families']}
    for row in d['bound_only_arithmetic_reductions']+d['bound_only_composition_reductions']:
        g=gates[row['family']];assert not g['fixed']and g['bounds']==row['bounds']
        assert all(frozenset((a,b))not in FREE for a,b,m in g['bounds'])
    assert len(d['bound_only_arithmetic_reductions'])==17 and len(d['bound_only_composition_reductions'])==6
    assert d['fixed']==[[3,14,3],[4,14,3]]
    expected=[[2,3,3],[2,4,3],[3,4,4],[5,6,3],[5,7,3],[6,7,4],[8,9,3],[8,10,3],[9,10,4],[12,13,4]]
    assert d['bounds']==expected and d['inherited_residual_cases']==155680
    old=read('results/fresh_h9_star_dense_frontier_20261014.json')['families']
    assert d['families'][:-1]==old and len(d['families'])==161
    last=d['families'][-1];assert last['fixed']==d['fixed']and last['bounds']==expected
    # Independently verify the proposed literal reversal, without depending
    # on the forthcoming symmetry artifact (which may bind this review).
    assert sorted(PERM)==list(range(18))and all(PERM[PERM[i]]==i for i in range(18))
    assert{frozenset(PERM[i]for i in v)for v in V}=={frozenset(v)for v in V}
    assert{frozenset(PERM[i]for i in m)for m in M}==M
    def key(f):return tuple(tuple(sorted(map(tuple,f[k])))for k in('fixed','bounds'))
    original={key(f)for f in d['families']};images=set()
    for family in d['families']:
        image={k:[sorted((PERM[a],PERM[b]))+[m]for a,b,m in family[k]]for k in('fixed','bounds')}
        images.add(key(image))
    assert len(original)==161 and len(images-original)==12 and len(original|images)==173
    paths=[primary,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_TWO_FREE_CUT_LEAVES_DEPENDENCY_REVIEW',sources={p:sha(p)for p in paths},
        recursive_manifests=len(seen),five_proof_dependencies=dependencies,needed_witness_counts=counts,
        needed_witness_occurrences=4956,removed_pairs=[sorted(p)for p in sorted(FREE,key=lambda x:sorted(x))],
        bound_only_reductions=23,inherited_arithmetic_cases=155680,
        independent_reversal_check=dict(permutation=PERM,actual_vertices=len(V),minimal_nonfaces=len(M),original_predicates=161,new_images=12,total=173),
        scope='Independent source/proof and exact literal support audit. No repeated determinant arithmetic. The four cut-leaf own pairs can take any finite labels, includingright; six positive-leaf own contacts remain>=3,14own3/3 exact, fourinternals>=4. Independently confirms proposed literal reversal and all161predicate images, without claiming remaining-case exhaustion.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}_review.json'
    if args.verify:assert d==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['needed_witness_occurrences'])
