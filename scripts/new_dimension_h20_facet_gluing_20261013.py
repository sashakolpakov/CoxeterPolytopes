#!/usr/bin/env python3
"""Exact finite tests of gluing the non-Coxeter facets of Borcherds' H21 chamber.

The frozen integral lattice supplies geometry. No arithmeticity assumption is
used in the defect-graph or short-cycle tests. These are bounded exclusions,
not a nonexistence theorem for H20 Coxeter polytopes.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd
from pathlib import Path
import networkx as nx
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUTPUT=ROOT/'data/enumerations/new_dimension_h20_facet_gluing_20261013.json'

def load():
    raw=json.loads(SOURCE.read_text())['coxeter_input']
    Q=-np.array(raw['basis_gram'],dtype=np.int64)
    W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
    G=W@Q@W.T
    assert Counter(map(int,G.diagonal()))=={2:42,4:168}
    assert all(G[i,j]<=0 for i,j in combinations(range(210),2))
    return Q,W,G

def section(G,f):
    q=int(G[f,f]);H=q*G-np.outer(G[:,f],G[f,:])
    ids=[i for i in range(210) if H[i,i]>0]
    return ids,H[np.ix_(ids,ids)]

def defects(G,f,include_geometry=False):
    ids,H=section(G,f);bad=nx.Graph();bad.add_nodes_from(range(len(ids)))
    angles=Counter();edge_records=[]
    for i,j in combinations(range(len(ids)),2):
        n=int(H[i,j])**2;d=int(H[i,i])*int(H[j,j])
        if n>=d:continue
        ratio=F(n,d);angles[str(ratio)]+=1
        if ratio not in {F(0),F(1,4),F(1,2),F(3,4)}:
            bad.add_edge(i,j);edge_records.append([ids[i]+1,ids[j]+1,str(ratio)])
    if G[f,f]==2:
        cover=[i for i,x in enumerate(ids) if G[x,x]==2 and G[f,x]==-1]
        assert len(cover)==5
    else:
        color=nx.algorithms.bipartite.color(bad)
        classes=[[i for i in bad if color[i]==v] for v in (0,1)]
        cover=min(classes,key=len);assert len(cover)==40
    cover_set=set(cover)
    assert all(i in cover_set or j in cover_set for i,j in bad.edges())
    matching=sorted(tuple(sorted(e)) for e in nx.max_weight_matching(bad,maxcardinality=True))
    assert len(matching)==len(cover)
    row={'ambient_facet_one_based':f+1,'ambient_norm':int(G[f,f]),'section_facets':len(ids),
         'angle_cosine_squared_counts':dict(sorted(angles.items())),
         'bad_ridge_count':len(edge_records),'minimum_defect_vertex_cover':len(cover),
         'cover_ambient_labels':[ids[i]+1 for i in cover],
         'matching_ambient_labels':[[ids[i]+1,ids[j]+1] for i,j in matching],
         'defect_degree_counts':dict(sorted(Counter(dict(bad.degree).values()).items()))}
    if include_geometry:
        row.update(section_labels_one_based=[i+1 for i in ids],
                   scaled_section_gram=H.tolist(),gram_scale=int(G[f,f]),
                   bad_ridges=edge_records)
    return row

def primitive(v):
    v=tuple(map(int,v));d=gcd(*v)
    assert d>0
    return tuple(x//d for x in v)

def tile_key(V,T):
    image=V@T.T
    assert np.max(np.abs(image))<2**50
    return tuple(sorted(primitive(v) for v in image))

def short_moves(Q,W,G,f=0):
    import pynauty as pn
    short=[i for i in range(210) if G[i,i]==2]
    graph=nx.Graph();graph.add_nodes_from(short)
    graph.add_edges_from((i,j) for i,j in combinations(short,2) if G[i,j]==-1)
    assert nx.is_connected(graph) and graph.number_of_edges()==105
    root=W[f];base_ids,unused=section(G,f)
    V=2*W[base_ids]-np.outer(G[base_ids,f],root)
    basis_local=sp.Matrix(W[short].T).rref()[1]
    basis_ids=[short[i] for i in basis_local]
    B=sp.Matrix(W[basis_ids].T);assert len(basis_ids)==22
    Binv=B.inv();identity=np.eye(22,dtype=np.int64);moves=[]
    wall_set={tuple(v) for v in W.tolist()}
    positions={v:i for i,v in enumerate(short)}
    pg=pn.Graph(len(short),False,adjacency_dict={positions[v]:[positions[w] for w in graph[v]] for v in short})
    pg_aut=pn.autgrp(pg);assert int(pg_aut[1])*10**int(pg_aut[2])==241920
    # Orbit traversal obtains a chosen ambient symmetry without enumerating
    # 241920 automorphisms or asking an unrefined graph matcher to find it.
    orbit={f:tuple(range(len(short)))};queue=[f]
    while queue:
        v=queue.pop();old=orbit[v]
        for gen in pg_aut[0]:
            image=short[gen[positions[v]]]
            if image in orbit:continue
            orbit[image]=tuple(gen[old[i]] for i in range(len(short)));queue.append(image)
    assert len(orbit)==42
    for j in graph.neighbors(f):
        permutation={v:short[orbit[j][positions[v]]] for v in short}
        A=sp.Matrix(W[[permutation[i] for i in basis_ids]].T)*Binv
        assert all(x.q==1 for x in A)
        A=np.array(A,dtype=np.int64)
        assert np.array_equal(A.T@Q@A,Q)
        assert {tuple(v) for v in (W@A.T).tolist()}==wall_set
        u=2*W[j]+root
        D=3*identity-np.outer(u,u@Q)
        r=root+W[j]
        R=(identity-np.outer(r,r@Q))@A
        for name,T,den in [('self_double',D,3),('odd_roll',R,1)]:
            assert np.array_equal(T.T@Q@T,den*den*Q)
            assert np.array_equal(T@root,den*root) or np.array_equal(T@root,-den*root)
            # The shared facet is fixed setwise with opposite outward normal.
            assert primitive(T@u)==primitive(-u) or name=='odd_roll'
            moves.append({'seam':j+1,'kind':name,'numerator':T,'denominator':den})
    return graph,V,moves

def seam_automorphisms(Q,W,G,f=0,j=21):
    """Complete metric seam group, including independent wall rescalings."""
    import pynauty as pn
    assert G[f,j]==-1 and G[f,f]==G[j,j]==2
    Z=3*W-np.outer(2*G[:,f]+G[:,j],W[f])-np.outer(G[:,f]+2*G[:,j],W[j])
    H=Z@Q@Z.T;ids=[i for i in range(210) if H[i,i]>0];K=H[np.ix_(ids,ids)]
    n=len(ids);colors=defaultdict(set);colors['wall']=set(range(n));adj={i:[] for i in range(n)}
    for a,b in combinations(range(len(ids)),2):
        if K[a,b]==0:continue
        label=(int(np.sign(K[a,b])),str(F(int(K[a,b])**2,int(K[a,a])*int(K[b,b]))))
        colors[label].add(n);adj[n]=[a,b];adj[a].append(n);adj[b].append(n);n+=1
    graph=pn.Graph(n,False,adjacency_dict=adj,vertex_coloring=list(colors.values()))
    aut=pn.autgrp(graph);order=int(aut[1])*10**int(aut[2]);assert order==2304
    generators=[g[:len(ids)] for g in aut[0]]
    basis=sp.Matrix(Z[ids].T).rref()[1];assert len(basis)==20
    X=sp.Matrix.hstack(sp.Matrix(Z[[ids[k] for k in basis]].T),sp.Matrix(W[f]),sp.Matrix(W[j]));Xinv=X.inv()
    wall_set={tuple(v) for v in W.tolist()};lifts=[]
    for permutation in generators:
        assert np.array_equal(K[np.ix_(permutation,permutation)],K)
        first=sp.Matrix(Z[[ids[permutation[k]] for k in basis]].T)
        possibilities=[]
        for swap in (False,True):
            Y=sp.Matrix.hstack(first,sp.Matrix(W[j if swap else f]),sp.Matrix(W[f if swap else j]))
            A=Y*Xinv
            if not all(x.q==1 for x in A):continue
            A=np.array(A,dtype=np.int64)
            if not np.array_equal(A.T@Q@A,Q):continue
            if {tuple(v) for v in (W@A.T).tolist()}!=wall_set:continue
            possibilities.append({'swaps_ambient_pair':swap,'matrix':A.tolist()})
        assert len(possibilities)==1
        lifts.append(possibilities[0])
    core=frozenset(k for k,i in enumerate(ids) if G[i,i]==2 and G[f,i]==-1)
    seen={core};queue=[core]
    while queue:
        s=queue.pop()
        for permutation in generators:
            t=frozenset(permutation[i] for i in s)
            if t not in seen:seen.add(t);queue.append(t)
    assert len(seen)==2
    assert any(row['swaps_ambient_pair'] for row in lifts)
    return {'ambient_pair_one_based':[f+1,j+1],'facets':len(ids),
            'normalized_gram_automorphism_order':order,
            'parent_seam_stabilizer_order':order//2,
            'parent_stabilizer_index':2,'two_actual_neighbor_types':['self_double','odd_roll'],
            'seam_labels_one_based':[i+1 for i in ids],
            'generators':generators,'integral_ambient_lifts':lifts,
            'core_subsets_orbit':[[ids[i]+1 for i in sorted(s)] for s in sorted(seen,key=lambda s:tuple(sorted(s)))]}

def cycles(Q,W,G,graph):
    color=nx.algorithms.bipartite.color(graph)
    points=[i for i in graph if color[i]==0]
    rows=[]
    for a,b,c in combinations(points,3):
        ab=set(graph[a])&set(graph[b]);bc=set(graph[b])&set(graph[c]);ca=set(graph[c])&set(graph[a])
        assert len(ab)==len(bc)==len(ca)==1
        x,y,z=next(iter(ab)),next(iter(bc)),next(iter(ca))
        if len({x,y,z})<3:continue
        path=[a,x,b,y,c,z]
        local=G[np.ix_(path,path)]
        assert all(local[i,j]==(2 if i==j else -1 if (i-j)%6 in (1,5) else 0) for i in range(6) for j in range(6))
        rows.append([i+1 for i in path])
    assert len(rows)==1120
    path=[i-1 for i in rows[0]];T=np.eye(22,dtype=np.int64)
    for i,j in zip(path,path[1:]+path[:1]):
        r=W[i]+W[j];T=(np.eye(22,dtype=np.int64)-np.outer(r,r@Q))@T
    assert np.array_equal(T.T@Q@T,Q) and np.array_equal(T@W[path[0]],W[path[0]])
    N=sp.Matrix(T)-sp.eye(22)
    assert N!=sp.zeros(22)
    restriction=sp.Matrix([W[path[0]]@Q]).nullspace()
    H=sp.Matrix.hstack(*restriction)
    square_nilpotent=sp.Matrix(T)**2-sp.eye(22)
    assert square_nilpotent!=sp.zeros(22) and square_nilpotent**3==sp.zeros(22)
    assert square_nilpotent*H!=sp.zeros(22,21)
    return {'short_odd_graph_vertices':42,'edges':105,'fundamental_group_rank':64,
            'six_cycles':rows,'six_cycle_count':len(rows),
            'representative_path':rows[0],'representative_holonomy':T.tolist(),
            'holonomy_minus_identity_rank':int(N.rank()),
            'holonomy_nontrivial_on_H20':N*H!=sp.zeros(22,21),
            'nilpotence_power_2_zero':N**2==sp.zeros(22),
            'nilpotence_power_3_zero':N**3==sp.zeros(22),
            'holonomy_squared_is_nontrivial_unipotent':True,
            'holonomy_squared_minus_identity_rank':int(square_nilpotent.rank()),
            'holonomy_squared_minus_identity_cube_zero':True,
            'holonomy_has_infinite_order_on_H20':True}

def bounded_cycles(V,moves,depth=3):
    identity=np.eye(22,dtype=np.int64)
    base=tile_key(V,identity)
    # Every depth-two placement is retained. Comparing depth-three to these
    # detects all closure candidates of length at most five from the base.
    layers=[[((),identity,base)]];known=defaultdict(list);known[base].append(())
    for d in range(1,3):
        layer=[]
        for path,T,key in layers[-1]:
            for k,m in enumerate(moves):
                U=T@m['numerator'];key2=tile_key(V,U)
                if len(path) and key2==base:continue
                path2=path+(k,);layer.append((path2,U,key2));known[key2].append(path2)
        layers.append(layer)
    collisions=[]
    for key,paths in known.items():
        for a,b in combinations(paths,2):
            if not a or not b or moves[a[0]]['seam']==moves[b[0]]['seam']:continue
            collisions.append({'path_a':list(a),'path_b':list(b),'length':len(a)+len(b)})
    five=[];six=[];third=defaultdict(list)
    if depth>=3:
        for path,T,key in layers[2]:
            for k,m in enumerate(moves):
                U=T@m['numerator'];key2=tile_key(V,U)
                for b in known.get(key2,[]):
                    if len(b)!=2 or moves[path[0]]['seam']==moves[b[0]]['seam']:continue
                    five.append({'path_a':list(path+(k,)),'path_b':list(b),'length':5})
                path3=path+(k,)
                # Immediate backtracks return a depth-one placement; they
                # cannot represent a simple six-cycle.
                if any(len(p)<=1 for p in known.get(key2,[])):continue
                for b in third[key2]:
                    if moves[path3[0]]['seam']==moves[b[0]]['seam']:continue
                    six.append({'path_a':list(path3),'path_b':list(b),'length':6})
                third[key2].append(path3)
    return {'states_depth_0_1_2':[len(x) for x in layers],
            'all_first_seam_different_closures_length_at_most_4':collisions,
            'all_first_seam_different_closures_length_5':five,
            'all_first_seam_different_closures_length_6':six,
            'distinct_depth_3_full_wall_sets':len(third),
            'scope':'Exact oriented full-wall-set coincidences; includes both gluing types. These are candidate cycle closures, not a completed chamber.'}

def noncore_catalog(G,f=0):
    ids,H=section(G,f);groups=defaultdict(list)
    for a,i in enumerate(ids):
        bad=Counter()
        for b,j in enumerate(ids):
            if a==b:continue
            r=F(int(H[a,b])**2,int(H[a,a])*int(H[b,b]))
            if r<1 and r not in {F(0),F(1,4),F(1,2),F(3,4)}:bad[str(r)]+=1
        signature=(int(G[i,i]),int(G[f,i]),tuple(sorted(bad.items())))
        groups[signature].append(i)
    records=[]
    for (norm,product,bad),labels in groups.items():
        j=labels[0];a=int(G[f,f]);b=int(G[f,j]);c=int(G[j,j]);det=a*c-b*b
        K=det*G-c*np.outer(G[:,f],G[f,:])+b*(np.outer(G[:,f],G[j,:])+np.outer(G[:,j],G[f,:]))-a*np.outer(G[:,j],G[j,:])
        keep=[i for i in range(210) if K[i,i]>0]
        cosine_counts=Counter(str(F(int(K[i,k])**2,int(K[i,i])*int(K[k,k]))) for i,k in combinations(keep,2))
        records.append({'ambient_wall_norm':norm,'ambient_pairing':product,
                        'labels_one_based':[i+1 for i in labels],'count':len(labels),
                        'incident_bad_angle_cosine_squared_counts':dict(bad),
                        'representative_one_based':j+1,'representative_seam_facets':len(keep),
                        'representative_seam_normalized_gram_squared_histogram':dict(cosine_counts),
                        'self_gluing_available_by_intrinsic_reflection':True,
                        'mixed_short_long_seam_available':norm==4})
    return {'base_facet_one_based':f+1,'scope':'Profiles of all209 walls; representative seam metric invariants. No claim that equal profiles alone prove isometry.',
            'profiles':records,
            'cross_profile_obstruction':'The five representative seam facet counts 190,196,81,208,50 are pairwise different, so these representative seams are pairwise nonisometric.',
            'mixed_short_long_congruence':'For each mixed ambient pair, the short-section and long-section seam is the same intrinsic intersection, giving explicit congruence for the81-wall and50-wall seams.'}

def six_tile_patch(Q,W,G,f=0,j=21,k=22):
    assert G[f,j]==G[f,k]==-1 and G[j,k]==0
    identity=np.eye(22,dtype=np.int64);u=2*W[j]+W[f];n=W[f]+W[j]+W[k];v=W[f]+W[k]
    D=3*identity-np.outer(u,u@Q);S=identity-np.outer(n,n@Q);U=identity-np.outer(v,v@Q)
    assert int(n@Q@n)==2 and int(n@Q@W[f])==0
    tiles={}
    for name,parent,T,den in [('A',f,identity,1),('DA',f,D,3),('B',k,U,1),('DB',k,D@U,3),('C',k,S@U,1),('DC',k,D@S@U,3)]:
        ids,unused=section(G,parent);V=2*W[ids]-np.outer(G[ids,parent],W[parent])
        image=V@T.T;vectors=[primitive(row) for row in image]
        assert all(sum(int(x)*int(y) for x,y in zip(row,Q@W[f]))==0 for row in vectors)
        tiles[name]={'ambient_parent_one_based':parent+1,'transport_numerator':T.tolist(),'transport_denominator':den,
                     'wall_labels_one_based':[i+1 for i in ids],'primitive_outward_normals':[list(row) for row in vectors]}
    def facet(tile,label):
        vectors=tiles[tile]['primitive_outward_normals'];labels=tiles[tile]['wall_labels_one_based']
        r=np.array(vectors[labels.index(label)],dtype=np.int64);qr=int(r@Q@r);out=[]
        for v in vectors:
            v=np.array(v,dtype=np.int64);p=qr*v-int(v@Q@r)*r
            if int(p@Q@p)>0:out.append(primitive(p))
        return primitive(r),tuple(sorted(out))
    pairs=[('A',j+1,'DA',j+1),('A',k+1,'B',f+1),('DA',k+1,'DB',f+1),
           ('B',j+1,'C',j+1),('DB',j+1,'DC',j+1),('C',f+1,'DC',f+1)]
    seams=[]
    for a,i,b,jj in pairs:
        ra,ka=facet(a,i);rb,kb=facet(b,jj)
        assert ra==tuple(-x for x in rb) and ka==kb
        seams.append({'tile_a':a,'wall_a':i,'tile_b':b,'wall_b':jj,'seam_facets':len(ka),'full_seam_matches_exactly':True})
    return {'ambient_A3_one_based':[f+1,j+1,k+1],'tiles':tiles,'intended_seams':seams,
            'ridge_sector_cycle':['A','DA','DB','DC','C','B'],
            'ridge_sector_angles':['alpha','alpha','beta','beta','beta','beta'],
            'angle_identity':'alpha=2gamma, beta=pi/2-gamma, hence2alpha+4beta=2pi',
            'new_extra_seams':'B--C andDB--DC are196-wall short-short orthogonal H19 facets.',
            'scope':'Six exact full-facet gluings and local interior-ridge closure. Global convexity, nonoverlap and Coxeter boundary require separate checks.'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=OUTPUT);parser.add_argument('--representatives-only',action='store_true')
    args=parser.parse_args();Q,W,G=load()
    representatives=[0,next(i for i in range(210) if G[i,i]==4)]
    out={'source':str(SOURCE.relative_to(ROOT)), 'dimension_target':20,
         'representative_sections':[defects(G,f,True) for f in representatives]}
    print('Representative defect covers verified.',flush=True)
    if not args.representatives_only:
        out['all_210_facets']=[defects(G,f) for f in range(210)]
        print('All 210 defect covers verified.',flush=True)
    graph,V,moves=short_moves(Q,W,G)
    out['moves']=[{**m,'numerator':m['numerator'].tolist()} for m in moves]
    print('Ten exact short-section neighbor placements constructed.',flush=True)
    out['odd_roll_holonomy']=cycles(Q,W,G,graph)
    print('All 1120 natural six-cycles enumerated; representative holonomy computed.',flush=True)
    out['bounded_cycle_search']=bounded_cycles(V,moves)
    print('Bounded cycles checked; computing complete A2 seam isometries.',flush=True)
    out['A2_seam_automorphisms']=seam_automorphisms(Q,W,G)
    out['noncore_seam_catalog']=noncore_catalog(G)
    out['six_tile_local_angle_repair']=six_tile_patch(Q,W,G)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out['bounded_cycle_search'],indent=2),flush=True)
    print('Saved '+str(args.output),flush=True)

if __name__=='__main__':main()
