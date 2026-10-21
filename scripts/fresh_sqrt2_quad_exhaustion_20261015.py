#!/usr/bin/env python3
"""Finite union-support enumeration for Q(sqrt2) quadrilateral stars.

This does not enumerate external white attachments and is not a dimension proof.
All shared white identifications and contacts on the union of the two chosen
minimal indefinite pieces are included. The finite-label alphabet and Lanner
piece list are separate published/catalogue inputs, bound in the final manifest.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import networkx as nx

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt2_quad_exhaustion_20261015'
TRIANGLES=((2,3,8),(2,4,8),(2,6,8),(2,8,8),(3,3,4),(3,8,8),(4,4,4),(4,6,6),(4,8,8))
SQ={3:Q(1,4),4:Q(1,2),6:Q(3,4)}

def pair(a,b):return tuple(sorted((a,b)))

def piece_key(n,g):
    # Local black nodes0,1 are fixed; only white names are quotiented.
    keys=[]
    for order in permutations(range(2,n)):
        perm=(0,1)+order
        keys.append(tuple(g[pair(perm[i],perm[j])] for i,j in combinations(range(n),2)))
    return min(keys)

def pieces():
    out=[dict(kind='L2',order=2,edges={(0,1):0})]
    for labels in TRIANGLES:
        for row in sorted(set(permutations(labels))):
            out.append(dict(kind='L3:'+','.join(map(str,labels)),order=3,
                            edges=dict(zip(combinations(range(3),2),row))))
    for n in (4,5):
        seen={}
        for order in permutations(range(1,n)):
            cyc=(0,)+order
            if cyc[1]>cyc[-1]:continue
            bonds=[pair(cyc[i],cyc[(i+1)%n]) for i in range(n)]
            for high in bonds:
                g={p:2 for p in combinations(range(n),2)}
                for p in bonds:g[p]=4 if p==high else 3
                key=piece_key(n,g)
                if key not in seen:
                    seen[key]=dict(kind='C'+str(n),order=n,
                                   edges=dict(zip(combinations(range(n),2),key)))
        out+=list(seen.values())
    return out

def maps(left_white,right_count):
    # Canonical names for new nodes; all identifications with old whites tried.
    def rec(values,unused,next_node):
        if len(values)==right_count:
            yield tuple(values);return
        for u in unused:
            yield from rec(values+[u],[v for v in unused if v!=u],next_node)
        yield from rec(values+[next_node],unused,next_node+1)
    yield from rec([],list(left_white),4+len(left_white))

def finite(nodes,g):
    """Exact positive definiteness of acute finite-alphabet Gram.

    A cycle has a nonpositive all-ones quadratic witness on its shortest cycle.
    A label>=6 cannot have a neighboring nonright edge in a positive diagram.
    Remaining forests are eliminated using rational squared entries.
    """
    adj={v:{} for v in nodes}
    for i,j in combinations(nodes,2):
        lab=g.get((i,j),2)
        if lab!=2:adj[i][j]=lab;adj[j][i]=lab
    for u in nodes:
        for v,lab in adj[u].items():
            if lab in (6,8) and (len(adj[u])!=1 or len(adj[v])!=1):return False
            if lab==0:return False
    values={v:Q(1) for v in nodes}
    while values:
        leaves=[v for v in values if len(adj[v])<=1]
        if not leaves:return False
        v=min(leaves)
        if values[v]<=0:return False
        if adj[v]:
            u,lab=next(iter(adj[v].items()))
            if lab==8:
                # This is an isolated finite I2(8) component.
                assert len(adj[u])==1
                del values[u],adj[u]
            else:
                values[u]-=SQ[lab]/values[v]
                del adj[u][v]
        del values[v],adj[v]
    return True

def no_short(W,g):
    for a,b in product((0,2),(1,3)):
        if g.get(pair(a,b),2)!=2:return False
        if any(g.get(pair(a,w),2)!=2 and g.get(pair(b,w),2)!=2 for w in W):return False
    return True

def corners_finite(W,g):
    return all(finite(sorted(W+[a,b]),g) for a,b in product((0,2),(1,3)))

def prepare(lp,rp,mapping):
    lm=[0,2]+list(range(4,lp['order']+2))
    rm=[1,3]+list(mapping)
    g={}
    for rec,mp in ((lp,lm),(rp,rm)):
        for (i,j),v in rec['edges'].items():
            p=pair(mp[i],mp[j])
            if p in g and g[p]!=v:return None
            g[p]=v
    n=max(lm+rm)+1;W=list(range(4,n))
    if not no_short(W,g) or not corners_finite(W,g):return None
    # A no-short candidate necessarily has right adjacent black pairs.
    for a,b in product((0,2),(1,3)):g[pair(a,b)]=2
    free=[p for p in combinations(range(n),2) if p not in g]
    return n,W,g,free,lm,rm

def key(n,g):return tuple(g.get(p,2) for p in combinations(range(n),2))

def connected_pieces(lm,rm,g):
    if set(lm)&set(rm):return True
    return any(g.get(pair(i,j),2)!=2 for i,j in product(lm,rm))

def encode_piece(p):return dict(kind=p['kind'],order=p['order'],labels=list(p['edges'].values()))

def isomorphism_classes(rows):
    graphs=[];classes=[]
    for index,row in enumerate(rows):
        G=nx.Graph()
        G.add_nodes_from((i,dict(black=i<4)) for i in range(row['nodes']))
        for p,label in zip(combinations(range(row['nodes']),2),row['labels']):
            G.add_edge(*p,label=label,opposite=p in ((0,2),(1,3)))
        for H,rec in zip(graphs,classes):
            match=nx.algorithms.isomorphism.GraphMatcher(G,H,node_match=lambda a,b:a==b,edge_match=lambda a,b:a==b)
            if match.is_isomorphic():
                mp=match.mapping
                assert set(mp)==set(range(row['nodes'])) and len(set(mp.values()))==row['nodes']
                assert all(G.nodes[v]==H.nodes[mp[v]] for v in G)
                assert all(G.edges[u,v]==H.edges[mp[u],mp[v]] for u,v in G.edges)
                rec['members'].append(dict(index=index,to_representative=[mp[i] for i in range(row['nodes'])]))
                break
        else:
            graphs.append(G);classes.append(dict(representative=index,members=[dict(index=index,to_representative=list(range(row['nodes'])))]))
    return classes

def enumerate_supports(limit=None):
    bank=pieces();counts=Counter();survivors={};checkpoints=[]
    for li,lp in enumerate(bank):
        for ri,rp in enumerate(bank):
            # Interchange opposite pairs, so unordered piece indices suffice.
            if ri<li:continue
            for mp in maps(list(range(4,lp['order']+2)),rp['order']-2):
                counts['base_identifications']+=1
                prep=prepare(lp,rp,mp)
                if prep is None:continue
                counts['base_corner_and_short_pass']+=1
                n,W,g,free,lm,rm=prep
                # Try right first. All later additions only worsen PD and
                # no-short conditions, so their failures are monotone prunes.
                def rec(k):
                    counts['dfs_nodes']+=1
                    if limit is not None and counts['dfs_nodes']>limit:
                        raise RuntimeError('diagnostic node limit reached; no completeness result')
                    if k==len(free):
                        if not connected_pieces(lm,rm,g):
                            counts['orthogonal_indefinite_pair']+=1;return
                        counts['labelled_completions']+=1
                        ky=(n,key(n,g))
                        if ky not in survivors:
                            survivors[ky]=dict(nodes=n,white=W,labels=list(ky[1]),
                                               left_piece=li,right_piece=ri,left_nodes=lm,right_nodes=rm,
                                               shared_white=sorted(set(lm)&set(rm)),
                                               added_contacts=[[a,b,g[(a,b)]] for a,b in free if g[(a,b)]!=2])
                        return
                    p=free[k]
                    g[p]=2;rec(k+1)
                    for label in (3,4,6,8):
                        g[p]=label
                        if no_short(W,g) and corners_finite(W,g):rec(k+1)
                    del g[p]
                rec(0)
        checkpoints.append(dict(left_piece=li,counts=dict(counts),distinct=len(survivors)))
        print('piece',li,'of',len(bank),'nodes',counts['dfs_nodes'],'survivors',len(survivors),flush=True)
    rows=list(survivors.values())
    return dict(piece_catalogue=list(map(encode_piece,bank)),counts=dict(counts),
                labelled_union_supports=rows,marked_isomorphism_classes=isomorphism_classes(rows),checkpoints=checkpoints)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--limit',type=int);args=ap.parse_args()
    data=enumerate_supports(args.limit)
    assert len(data['piece_catalogue'])==51 and len(data['labelled_union_supports'])==115
    assert len(data['marked_isomorphism_classes'])==52
    data['status']='PASS_FINITE_MINIMAL_PAIR_UNION_SUPPORT_ENUMERATION_EXTERNAL_WHITES_OPEN'
    data['scope']='No marked distance<=2 corner; all pairings on union of two minimal Lanner pieces. No external white attachments, rank completion or global dimension theorem claimed.'
    cat=json.loads((ROOT/'results/fresh_sqrt2_lanner_catalogue_20261015.json').read_text())
    assert cat['finite_labels']==[2,3,4,6,8]
    assert sorted(tuple(r['labels']) for r in cat['order_three']['accepted'])==sorted(TRIANGLES)
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt2_lanner_catalogue_20261015.py','results/fresh_sqrt2_lanner_catalogue_20261015.json']
    data['sources']={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    target=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert data==json.loads(target.read_text())
    else:target.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],'distinct',len(data['labelled_union_supports']))

if __name__=='__main__':main()
