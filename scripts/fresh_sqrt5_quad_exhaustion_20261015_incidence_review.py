#!/usr/bin/env python3
"""Independent union coverage replay: Dynkin recognition and reverse-order DFS.

No producer imports and no Gram arithmetic from the producer. The frozen,
independently reviewed Lanner catalogue supplies the minimal piece classes.
"""
import argparse,json
from collections import Counter,defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_sqrt5_quad_exhaustion_20261015';TAG=BASE+'_incidence_review'

def pair(a,b):return tuple(sorted((a,b)))



def finite(nodes,g):
    nodes=tuple(sorted(nodes))
    return finite_word(len(nodes),tuple(g.get(pair(i,j),2)for i,j in combinations(nodes,2)))

@lru_cache(None)
def finite_word(n,labels):
    nodes=range(n);g=dict(zip(combinations(nodes,2),labels))
    """Recognize the finite Dynkin forests, rather than eliminate the matrix."""
    adj={i:[]for i in nodes}
    for i,j in combinations(sorted(nodes),2):
        m=g.get((i,j),2)
        if m!=2:
            if m==0:return False
            adj[i].append(j);adj[j].append(i)
    unseen=set(nodes)
    while unseen:
        root=min(unseen);todo=[root];component={root}
        while todo:
            i=todo.pop()
            for j in adj[i]:
                if j not in component:component.add(j);todo.append(j)
        unseen-=component;n=len(component)
        edges=[(i,j,g[pair(i,j)])for i in component for j in adj[i]if i<j]
        if len(edges)!=n-1:return False
        if n<=2:continue
        if any(m>=6 for i,j,m in edges):return False
        five=[(i,j)for i,j,m in edges if m==5]
        if five:
            if len(five)!=1 or n not in(3,4) or any(m not in(3,5)for i,j,m in edges):return False
            if any(len(adj[i])>2 for i in component):return False
            i,j=five[0]
            if len(adj[i])==1 or len(adj[j])==1:continue
            return False
        four=[(i,j)for i,j,m in edges if m==4]
        if four:
            if len(four)!=1 or any(len(adj[i])>2 for i in component):return False
            i,j=four[0]
            if len(adj[i])==1 or len(adj[j])==1:continue
            if n==4:continue
            return False
        branches=[i for i in component if len(adj[i])>2]
        if not branches:continue
        if len(branches)!=1 or len(adj[branches[0]])!=3:return False
        center=branches[0];lengths=[]
        for first in adj[center]:
            prev,cur=center,first;length=1
            while len(adj[cur])==2:
                nxt=next(x for x in adj[cur]if x!=prev);prev,cur=cur,nxt;length+=1
            lengths.append(length)
        lengths.sort()
        if lengths[:2]==[1,1]:continue
        if lengths in[[1,2,2],[1,2,3],[1,2,4]]:continue
        return False
    return True


def corners(n,g):
    W=list(range(4,n))
    return all(finite(W+[i,j],g)for i,j in product((0,2),(1,3)))


def no_short(n,g):
    return all(g.get(pair(i,j),2)==2 and not any(g.get(pair(i,w),2)!=2 and g.get(pair(j,w),2)!=2 for w in range(4,n))for i,j in product((0,2),(1,3)))


def white_key(n,g):
    return min(tuple(g.get(pair(p[i],p[j]),2)for i,j in combinations(range(n),2))for tail in permutations(range(2,n))for p in[(0,1)+tail])


def independent_pieces(saved):
    cat=json.loads((ROOT/'results/fresh_sqrt5_lanner_catalogue_20261015.json').read_text())
    keys={(2,(0,))}
    for row in cat['order_three']['accepted']:
        keys.update((3,p)for p in permutations(row['labels']))
    for entry in cat['higher_orders']:
        n=entry['order']
        for row in entry['accepted_classes']:
            template=dict(zip(combinations(range(n),2),row['labels']))
            for p in permutations(range(n)):
                g={pair(i,j):template[pair(p[i],p[j])]for i,j in combinations(range(n),2)}
                keys.add((n,white_key(n,g)))
    literal=[(r['order'],tuple(r['labels']))for r in saved]
    assert len(keys)==len(literal)==153 and set(literal)==keys
    # Input ordering is solely a deterministic naming convention.
    return [(n,dict(zip(combinations(range(n),2),labels)))for n,labels in literal]

def identifications(left,right):
    for k in range(min(left,right)+1):
        for R in combinations(range(right),k):
            for L in permutations(range(4,4+left),k):
                d=dict(zip(R,L));nxt=4+left;mapping=[]
                for j in range(right):
                    if j in d:mapping.append(d[j])
                    else:mapping.append(nxt);nxt+=1
                yield tuple(mapping)


def union_enumeration(base):
    bank=independent_pieces(base['piece_catalogue']);result=set();count=Counter()
    for li,(ln,lg)in enumerate(bank):
        for rn,rg in bank[li:]:
            for rmw in identifications(ln-2,rn-2):
                count['identifications']+=1;lm=(0,2)+tuple(range(4,ln+2));rm=(1,3)+rmw;g={};consistent=True
                for old,mp in[(lg,lm),(rg,rm)]:
                    for (i,j),m in old.items():
                        e=pair(mp[i],mp[j])
                        if e in g and g[e]!=m:consistent=False;break
                        g[e]=m
                    if not consistent:break
                if not consistent:continue
                n=max(lm+rm)+1
                if not no_short(n,g)or not corners(n,g):continue
                for i,j in product((0,2),(1,3)):g[pair(i,j)]=2
                free=list(reversed([e for e in combinations(range(n),2)if e not in g]))
                def rec(k):
                    count['dfs_nodes']+=1
                    if k==len(free):
                        if not(set(lm)&set(rm))and not any(g.get(pair(i,j),2)!=2 for i,j in product(lm,rm)):
                            return
                        count['completions']+=1;result.add((n,tuple(g.get(e,2)for e in combinations(range(n),2))));return
                    e=free[k]
                    for m in(10,6,5,4,3,2):
                        g[e]=m
                        if no_short(n,g)and corners(n,g):rec(k+1)
                    del g[e]
                rec(0)
    expected={(r['nodes'],tuple(r['labels']))for r in base['labelled_union_supports']}
    assert count['identifications']==127875 and count['completions']==5502
    assert result==expected and len(result)==5440
    return dict(counts=dict(count),exact_literal_support_set_matches=True)



def canonical_marked(row):
    n=row['nodes'];g=dict(zip(combinations(range(n),2),row['labels']));keys=[]
    for p in permutations(range(4)):
        if {frozenset((p[0],p[2])),frozenset((p[1],p[3]))}!={frozenset((0,2)),frozenset((1,3))}:continue
        groups=defaultdict(list)
        for w in range(4,n):
            signature=tuple(g[pair(w,b)]for b in p)+tuple(sorted(g[pair(w,v)]for v in range(4,n)if v!=w))
            groups[signature].append(w)
        blocks=[groups[k]for k in sorted(groups)]
        for chosen in product(*(permutations(b)for b in blocks)):
            order=p+tuple(v for group in chosen for v in group)
            keys.append(tuple(g[pair(order[i],order[j])]for i,j in combinations(range(n),2)))
    return n,min(keys)


def class_check(base):
    rows=base['labelled_union_supports'];covered=[];keys=set();mapped_pairs=0
    for rec in base['marked_isomorphism_classes']:
        representative=rows[rec['representative']];key=canonical_marked(representative)
        assert key not in keys;keys.add(key)
        n=representative['nodes'];target=dict(zip(combinations(range(n),2),representative['labels']))
        for member in rec['members']:
            index=member['index'];source=rows[index];mp=member['to_representative'];covered.append(index)
            assert source['nodes']==n and sorted(mp)==list(range(n))
            assert all((i<4)==(mp[i]<4)for i in range(n))
            assert {frozenset((mp[0],mp[2])),frozenset((mp[1],mp[3]))}=={frozenset((0,2)),frozenset((1,3))}
            for (i,j),m in zip(combinations(range(n),2),source['labels']):
                assert m==target[pair(mp[i],mp[j])];mapped_pairs+=1
    assert sorted(covered)==list(range(len(rows))) and len(keys)==1180
    return dict(classes=len(keys),literal_maps=len(covered),mapped_pairs=mapped_pairs,
                separation='Independent canonicalization by opposite-pair-preserving black permutations and signature-grouped white permutations')


def run():
    base=json.loads((ROOT/'results'/f'{BASE}.json').read_text())
    enumeration=union_enumeration(base);classes=class_check(base)
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',
           'scripts/fresh_sqrt5_lanner_catalogue_20261015.py','results/fresh_sqrt5_lanner_catalogue_20261015.json',
           'scripts/fresh_sqrt5_lanner_catalogue_20261015_review.py','results/fresh_sqrt5_lanner_catalogue_20261015_review.json']
    return dict(status='PASS_INDEPENDENT_COMPLETE_QSQRT5_MINIMAL_PAIR_UNION_REPLAY',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                union_enumeration=enumeration,marked_classes=classes,
                arithmetic_scope='Independent finite Dynkin forest recognition instead of Qphi elimination; no producer/project imports',
                scope='All finite contacts on chosen minimal-piece union only. External white attachments, arithmetic full Gram constraints and rank completion remain open.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print(out['union_enumeration']);print(out['marked_classes'])
