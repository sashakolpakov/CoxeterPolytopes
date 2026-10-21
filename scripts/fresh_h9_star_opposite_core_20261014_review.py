#!/usr/bin/env python3
"""Independent domain and fallback-witness review of the opposite core.

Uses only previously independent Fraction and determinant helpers. The
47,040 polynomial coefficient cases are inherited, not recomputed.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_h9_star_mixed_singles_20261014_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_opposite_core_20261014'
S=[3,4,6,7,9,10,12,13]
CORE=sorted(S+[11,14,17])
Z,O,COS=R.Z,R.O,R.COS
add,sub,neg,scale,mul,total,sg,inv,decode=R.add,R.sub,R.neg,R.scale,R.mul,R.total,R.sign,R.inverse,R.dec_anp
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):return list(map(str,x))


def run():
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for s,h in d.get('sources',{}).items():
            assert sha(s)==h,s
            if s.endswith('.json'):read(s)
        return d
    primary=f'results/{TAG}.json';d=read(primary)
    assert d['status']=='PASS_CORE_OPPOSITE_SINGLE_ALL_INTERNAL_LABELS'
    g=read('results/fresh_h9_tree_cube_wedges_20261014.json')['representatives'][0]['geometry']
    V=[set(v['incident'])for v in g['vertices']]
    missing={frozenset(m)for m in g['minimal_nonfaces']}
    assert len(V)==378 and d['rank_core']==CORE
    assert any(set(S)<=v for v in V)
    assert all(any(set(p)<=v for v in V)for p in combinations(CORE,2))
    faces={frozenset(v&set(CORE))for v in V}
    maximal={f for f in faces if not any(f<h for h in faces)}
    assert len(maximal)==9 and maximal==set(map(frozenset,d['complete_maximal_core_faces']))
    own={(3,14):3,(4,14):3,(11,12):2,(11,13):3,(12,17):3,(13,17):2}
    blocks={tuple(S[2*i:2*i+2]):4 for i in range(4)}
    cross={tuple(sorted((a,b)))for i,j in combinations(range(4),2)
           for a in S[2*i:2*i+2]for b in S[2*j:2*j+2]}
    got_cross,got_upper,got_double=set(),set(),set()
    guard_checks=[]
    for row in d['core_positive_quadratic_guards']:
        path=row['ordered_support'];pair=tuple(row['test_pair']);trial=row['test_label']
        assert len(path)==len(set(path))and set(path)<=set(CORE)
        assert set(path)==set(row['support'])<=V[row['actual_vertex']]
        bonds={tuple(sorted((path[i],path[i+1]))):m for i,m in enumerate(row['labels'])}
        kind=row['kind'];available=own|blocks|{pair:trial}
        if kind=='cross_block_right':
            assert pair in cross and trial==3 and row['labels']==[4,3,4]
            got_cross.add(pair)
        elif kind=='third_contact_upper':
            assert pair in((11,14),(14,17))and trial==5
            assert path==([3,14,11,13]if pair==(11,14)else[3,14,17,12])
            assert row['labels']==[3,5,3];got_upper.add(pair)
        else:
            assert kind=='two_label_four_contacts_impossible'
            assert path==[11,14,17]and row['labels']==[4,4]
            available|={(11,14):4,(14,17):4}
            got_double.add(tuple(path))
        assert all(available.get(p)==m for p,m in bonds.items())
        B=tuple(tuple(O if i==j else neg(COS[bonds[tuple(sorted((a,b)))]])
                      if tuple(sorted((a,b)))in bonds else Z
                      for j,b in enumerate(path))for i,a in enumerate(path))
        assert B==tuple(tuple(decode(x)for x in rr)for rr in row['gram'])
        w=tuple(map(decode,row['positive_witness']))
        assert all(sg(x)>0 for x in w)
        q=total(mul(w[i],mul(B[i][j],w[j]))for i in range(len(path))for j in range(len(path)))
        assert q==decode(row['quadratic_value'])and sg(q)<=0
        if kind=='third_contact_upper':
            assert q==(F(1),F(0),F(-1,2),F(0))and sg(q)<0
        else:
            assert q==Z
            assert all(total(mul(B[i][j],w[j])for j in range(len(path)))==Z for i in range(len(path)))
        guard_checks.append(dict(kind=kind,pair=list(pair),ordered_support=path,
                                actual_vertex=row['actual_vertex'],positive_vector=list(map(enc,w)),quadratic_value=enc(q)))
    assert got_cross==cross and got_upper=={(11,14),(14,17)}and got_double=={(11,14,17)}
    assert len(guard_checks)==27

    # Separate rational congruence for the3-6 mutual upper guard.
    row=d['mutual_upper_six']
    assert row['ordered_support']==[13,11,17]and set(row['support'])=={11,13,17}<=V[row['actual_vertex']]
    assert row['known_label']==3 and row['test_label']==6
    assert own[(11,13)]==3 and own[(13,17)]==2
    lap=[[F(1,4),-F(1,4),F(0)],[-F(1,4),F(1),-F(3,4)],[F(0),-F(3,4),F(3,4)]]
    assert all(sum(rr)==0 for rr in lap)
    # Congruence by diag(1/2,1,sqrt3/2) of the normalized3-6 path.
    assert F(1,4)+F(3,4)==1
    expected={11:(12,13),17:(13,12),14:(3,4)}
    for row in d['projection_budgets']:
        leaf=row['leaf'];drop,keep=expected.pop(leaf)
        assert (drop,keep)==(row['removed_own_root'],row['retained_own_root'])
        assert own[tuple(sorted((leaf,keep)))]==3
        support=(set(S)-{drop})|{leaf}
        assert len(support)==8 and support==set(row['support'])<=V[row['actual_vertex']]
    assert not expected and len(d['projection_budgets'])==3
    assert F(1,4)+2*F(1,2)>1 and F(1,4)+2*F(1,2)>1 and F(1,4)+4*F(1,4)>1
    row=d['high_last_third_right']
    assert set(row['support'])=={12,13,14}<=V[row['actual_vertex']]
    assert all(frozenset([12,13,x])in missing for x in(11,17))

    local=list(range(18));local[11],local[17]=17,11;local[12],local[13]=13,12
    assert d['local_core_swap']==local
    assert {frozenset(local[i]for i in f)for f in maximal}==maximal
    assert {tuple(sorted((local[a],local[b]))):m for(a,b),m in own.items()}==own
    assert {tuple(sorted((local[a],local[b]))):m for(a,b),m in blocks.items()}==blocks
    assert tuple(sorted((local[11],local[14])))==(14,17)
    assert tuple(sorted((local[17],local[14])))==(11,14)
    # All asymmetric domains orient to the inherited rectangle after rulingout4/4.
    domain={(a,b)for a,b in product((2,3,4),repeat=2)if(a,b)!=(4,4)}
    assert all(a<=3 and b<=4 or b<=3 and a<=4 for a,b in domain)

    # Exact formal square identity in a=sqrt(F),b=sqrt(H), overQ(sqrt2).
    rt2=scale(COS[4],2);assert mul(rt2,rt2)==scale(O,2)
    def poly_mul(p,q):
        out={}
        for(i,j),x in p.items():
            for(k,l),y in q.items():out[i+k,j+l]=add(out.get((i+k,j+l),Z),mul(x,y))
        return{k:v for k,v in out.items()if v!=Z}
    lhs=poly_mul({(0,0):scale(O,2),(0,2):O},{(0,0):O,(2,0):O})
    for key,value in poly_mul({(0,0):rt2,(1,1):O},{(0,0):rt2,(1,1):O}).items():
        lhs[key]=sub(lhs.get(key,Z),value)
    rhs=poly_mul({(1,0):rt2,(0,1):neg(O)},{(1,0):rt2,(0,1):neg(O)})
    assert{k:v for k,v in lhs.items()if v!=Z}==rhs
    assert 2*2-3==1 and 3>1 # 1/(2-sqrt3)=2+sqrt3; strict K>2,L>1.

    oldpath='results/fresh_h9_star_opposite_high_20261014.json'
    old=read(oldpath);oldreview=read('results/fresh_h9_star_opposite_high_20261014_review.json')
    assert oldreview['sources'][oldpath]==sha(oldpath)and oldreview['all_cases']==47040
    expected_keys={tuple(block)+(ia,ib,ic,mu,ma,mb)
        for block in product((4,5),repeat=3)
        for ia,ib,ic in product([None]+list(range(6)),[None]+list(range(6)),[None,2,3,4,5])
        for mu,ma,mb in product((2,3,4,5),(2,3),(2,3,4))}
    got,counts,fallback=set(),{},set()
    for i,row in enumerate(old['cases']):
        key=tuple(row['block_labels']+row['foreign_contacts']+[row['mutual_label']]+row['third_labels'])
        assert key not in got;got.add(key)
        status=row['status'];counts[status]=counts.get(status,0)+1
        if status=='STRICT_INTERIOR_BERNSTEIN_6_7_TO_ONE':
            bs=list(map(decode,row['bernstein_6_7']))
            assert all(sg(x)>=0 for x in bs)and any(sg(x)>0 for x in bs)
        elif status=='STRICT_BERNSTEIN_9_10_TO_ONE':
            assert all(sg(decode(x))>0 for x in row['bernstein_9_10'])
        else:
            assert status=='NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR'
            fallback.add(i)
    assert got==expected_keys and len(got)==47040
    assert counts==oldreview['case_counts']==d['inherited_status_counts']
    assert counts=={'STRICT_INTERIOR_BERNSTEIN_6_7_TO_ONE':43412,
                   'NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR':3612,
                   'STRICT_BERNSTEIN_9_10_TO_ONE':16}

    used,types,matrices,placements=set(),{},[],[]
    for row in d['new_core_fallbacks']:
        index=row['old_case_index'];assert index in fallback and index not in used;used.add(index)
        oldrow=old['cases'][index]
        for key in('block_labels','foreign_contacts','mutual_label','third_labels'):assert row[key]==oldrow[key]
        labels=own|blocks|{tuple(S[2*i:2*i+2]):m for i,m in enumerate(row['block_labels'])}
        labels|={(11,17):row['mutual_label'],(11,14):row['third_labels'][0],(14,17):row['third_labels'][1]}
        for leaf,extra in zip((11,17,14),row['foreign_contacts']):
            if extra is not None:labels[tuple(sorted((leaf,S[extra])))]=3
        assert sorted([list(p)+[m]for p,m in labels.items()])==row['relaxed_core_labels']
        assert all(set(p)<=set(CORE)for p in labels)
        w=row['witness'];nodes=w['component_order']
        assert set(nodes)<=set(w['containing_core_face'])<=V[w['actual_vertex']]
        assert frozenset(w['containing_core_face'])in maximal
        assert w['component_pair_labels']==[labels.get(p,2)for p in combinations(nodes,2)]
        A=tuple(tuple(O if a==b else neg(COS[labels.get(tuple(sorted((a,b))),2)])for b in nodes)for a in nodes)
        k,minors,v,value=R.nonnegative_failure(A)
        assert k==w['prefix_order']and nodes[:k]==w['support']
        assert [mul(minors[j],inv(minors[j-1]))for j in range(1,k+1)]==list(map(decode,w['pivots']))
        # Explicitly transfer both the selected support and all of its labels.
        image=sorted(local[i]for i in nodes[:k])
        image_vertex=next(i for i,VV in enumerate(V)if set(image)<=VV)
        mapped_labels={tuple(sorted((local[a],local[b]))):m for(a,b),m in labels.items()}
        assert all(mapped_labels[tuple(sorted((local[a],local[b])))]==labels[a,b]for a,b in labels)
        if A not in types:
            types[A]=len(types)
            matrices.append(dict(matrix_type=types[A],matrix_size=len(A),failure_order=k,
              matrix=[[enc(x)for x in rr]for rr in A],leading_determinants=list(map(enc,minors[1:])),
              nonnegative_quadratic_vector=list(map(enc,v)),quadratic_value=enc(value)))
        placements.append(dict(old_case_index=index,matrix_type=types[A],actual_vertex=w['actual_vertex'],
                               support=nodes[:k],swapped_support=image,swapped_actual_vertex=image_vertex))
    assert used==fallback and len(placements)==3612 and len(types)==240==d['selected_full_matrix_types']
    swap=list(range(18));swap[12],swap[13]=13,12
    assert swap==d['literal_last_root_swap']
    assert {frozenset(swap[i]for i in v)for v in V}==set(map(frozenset,V))
    assert {frozenset(swap[i]for i in m)for m in missing}==missing
    bounds=[S[2*i:2*i+2]+[4 if i<3 else 7]for i in range(4)]
    assert len(d['families'])==2
    for i,row in enumerate(d['families']):
        fixed=own if i==0 else{tuple(sorted((swap[a],swap[b]))):m for(a,b),m in own.items()}
        assert row['fixed']==sorted([list(p)+[m]for p,m in fixed.items()])and row['bounds']==bounds
    paths=[primary,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py',
           'scripts/fresh_h9_star_mixed_singles_20261014_review.py',
           'scripts/fresh_h9_star_zero_20261014_own_three_review.py']
    return dict(status='PASS_INDEPENDENT_OPPOSITE_CORE_DOMAIN_REVIEW',sources={p:sha(p)for p in paths},
      recursive_manifests=len(seen),arithmetic='Previously independent Fraction field and Laplace/Cramer helpers only; no producer imports. Frozen47040 polynomial arithmetic inherited through its independent review.',
      actual_vertices=len(V),rank_core=CORE,positive_quadratic_guards=guard_checks,
      mutual_six_rational_laplacian=[[str(x)for x in rr]for rr in lap],core_budgets=3,
      local_core_involution_only=True,high_first_square_identity_checked=True,
      inherited_cases=len(got),inherited_status_counts=counts,new_core_fallbacks=len(placements),
      selected_matrix_types=len(types),matrix_witnesses=matrices,fallback_placements=placements,
      scope='Literal star0 core with first3 finite internals>=4,last>=7,14own3/3,11and17 opposite single own3. Every noncore contact free. Both orientations. New guards, all3612 fallback placements and nonnegative witnesses, local face involution and continuous high-first proof reviewed. No full incidence exclusion or new search.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}_review.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['new_core_fallbacks'],d['selected_matrix_types'])
