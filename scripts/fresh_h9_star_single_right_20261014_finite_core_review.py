#!/usr/bin/env python3
"""Independent rational-Laplacian review of the finite core domain."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_single_right_20261014_finite_core'
S=[3,4,6,7,9,10,12,13];LEAVES=[11,14,17]
Z=(F(0),F(0));ONE=(F(1),F(0))
COS={3:(F(1,2),F(0)),4:(F(0),F(1,2))}
def mul(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def dec(x):return tuple(map(F,x))
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
    geometry=read('results/fresh_h9_tree_cube_wedges_20261014.json')['representatives'][0]['geometry']
    V=[set(v['incident'])for v in geometry['vertices']]
    core=set(S+LEAVES);assert len(core)==11 and sorted(core)==d['rank_core']
    assert any(set(S)<=v for v in V)
    assert all(any(set(p)<=v for v in V)for p in combinations(sorted(core),2))
    own={tuple(sorted((leaf,r))):3 for leaf,roots in[(14,[3,4]),(11,[12,13]),(17,[12,13])]for r in roots}
    assert sorted([list(p)+[m]for p,m in own.items()])==d['core_own_equalities']
    blocks={tuple(S[2*i:2*i+2]):4 for i in range(4)}
    expected_cross={tuple(sorted((a,b)))for i,j in combinations(range(4),2)for a in S[2*i:2*i+2]for b in S[2*j:2*j+2]}
    checked_cross=set();checked_mutual=set();checked_upper=set();matrices=[]
    for r in d['positive_kernel_guards']:
        path=r['ordered_support'];n=len(path);assert len(set(path))==n and set(path)<=core
        assert set(path)<=V[r['actual_vertex']]
        pair=tuple(r['test_pair']);trial=r['test_label'];kind=r['kind']
        available=own|blocks|{pair:trial}
        bonds={tuple(sorted((path[i],path[i+1]))):m for i,m in enumerate(r['labels'])}
        if r['complete_triangle']:bonds[tuple(sorted((path[0],path[-1])))]=3
        assert all(available.get(p)==m for p,m in bonds.items())
        if kind=='cross_block_right':
            assert pair in expected_cross and trial==3 and r['labels']==[4,3,4]
            checked_cross.add(pair)
        elif kind=='mutual_anchor_right':
            assert pair==(11,17)and trial==3 and set(path)=={11,12,17}and r['complete_triangle']
            checked_mutual.add(pair)
        else:
            assert kind=='third_contact_upper'and pair in((11,14),(14,17))and trial==4
            assert path==([3,14,11,12,17]if pair==(11,14)else[3,14,17,12,11])
            assert r['labels']==[3,4,3,3]and r['actual_vertex']==294
            checked_upper.add(pair)
        B=[]
        for i,a in enumerate(path):
            row=[]
            for j,b in enumerate(path):
                p=tuple(sorted((a,b)))
                row.append(ONE if i==j else tuple(-x for x in COS[bonds[p]])if p in bonds else Z)
            B.append(row)
        assert B==[[dec(x)for x in row]for row in r['gram']]
        k=list(map(dec,r['positive_kernel']));assert len(k)==n
        assert all(a>=0 and b>=0 and(a or b)for a,b in k)
        L=[]
        for i in range(n):
            row=[]
            for j in range(n):
                x=mul(mul(k[i],B[i][j]),k[j]);assert not x[1]
                row.append(x[0])
            L.append(row)
        assert all(sum(row)==0 for row in L)
        assert all(L[i][i]>0 for i in range(n))and all(L[i][j]<=0 for i in range(n)for j in range(n)if i!=j)
        matrices.append(dict(kind=kind,test_pair=list(pair),ordered_support=path,
            actual_vertex=r['actual_vertex'],rational_congruent_laplacian=[[str(x)for x in row]for row in L]))
    assert checked_cross==expected_cross and checked_mutual=={(11,17)}and checked_upper=={(11,14),(14,17)}and len(matrices)==27
    got=set()
    for row in d['foreign_budgets']:
        leaf=row['leaf'];assert leaf in LEAVES
        roots=[3,4]if leaf==14 else[12,13]
        assert row['removed_own_root']==roots[0]and row['retained_own_root']==roots[1]
        support=(set(S)-{roots[0]})|{leaf}
        assert len(support)==8 and support==set(row['support'])and support<=core and support<=V[row['actual_vertex']]
        assert row['own_cost']=='1/4'and row['foreign_inverse_diagonal_lower_bound']=='2';got.add(leaf)
    assert got==set(LEAVES)and len(d['foreign_budgets'])==3
    # For every4/5block c>=sqrt2/2 andc<1, its inverse is
    # (1-c²)^-1[[1,c],[c,1]], entrywise positive with diagonal>=2.
    assert F(1,4)+2*F(1,4)*2==F(5,4)>1
    assert F(1,4)+F(1,2)*2==F(5,4)>1
    oldpath='results/fresh_h9_star_zero_20261014_own_three.json'
    old=read(oldpath);review=read('results/fresh_h9_star_zero_20261014_own_three_review.json')
    assert review['status'].startswith('PASS_')and review['sources'][oldpath]==sha(oldpath)
    assert old['zero_rank_cases']==[]and old['unresolved']==[]
    oldown={(a,b):m for a,b,m in old['all_own_contacts']}
    assert all(oldown[p]==m for p,m in own.items())
    expected_cases=set(product(range(7),range(7),range(7),(2,3),(2,3)))
    seen_blocks=set();count=0
    for record in old['records']:
        block=tuple(record['block_labels']);assert block not in seen_blocks;seen_blocks.add(block)
        keys=set()
        for row in record['cases']:
            key=tuple(row['patterns']+row['contact_labels']);assert key not in keys;keys.add(key)
            # Canonical degree<4 ANP values in the previously independently
            # verified primitive field. This is a coverage/nonzero check,
            # not a second arithmetic evaluation of these determinants.
            coeff=list(map(F,row['determinant']));assert len(coeff)<=4 and any(coeff)
            count+=1
        assert keys==expected_cases
    assert seen_blocks==set(product((4,5),repeat=4))and count==21952
    assert review['nonzero_rank11_determinants']==count
    assert d['inherited_rank_cases']==count and d['inherited_zero_cases']==0
    predicates=set()
    for f in d['families']:
        assert f['bounds']==[]
        fixed={(a,b):m for a,b,m in f['fixed']}
        assert set(fixed)==set(own)|set(blocks)and all(fixed[p]==m for p,m in own.items())
        predicates.add(tuple(fixed[p]for p in blocks))
    assert predicates==seen_blocks and len(d['families'])==16
    paths=[primary,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_FINITE_CORE_DOMAIN_REVIEW',sources={p:sha(p)for p in paths},
        arithmetic='Standalone Fraction rational-Laplacian congruences for27 new domain guards; no producer imports. Frozen21952 determinant arithmetic inherited through its prior independent hash-bound review.',
        recursive_manifests=len(seen),actual_vertices=len(V),rank_core=sorted(core),
        rational_kernel_guards=matrices,core_foreign_budgets=3,inherited_nonzero_rank_cases=count,
        internal_block_choices=16,projection_triples_per_block=343,third_label_pairs_per_block=4,
        scope='Literal star0, six core own equalities at3 and four internal labels independently4or5. Every noncore own contact is free. Checks new guard dependencies and complete inherited case coverage; no new bank arithmetic, high-internal theorem or whole-incidence traversal.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}_review.json'
    if args.verify:assert d==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['inherited_nonzero_rank_cases'])
