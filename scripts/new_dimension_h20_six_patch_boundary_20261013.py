#!/usr/bin/env python3
"""Exact signed ridge census of the three/six-tile H20 angle-repair patches.

Global support/nonoverlap is independently certified in the h20_gluing audit.
This script recognizes actual boundary ridges by exact normal-plane equality
with a genuine ridge of an original acute tile. It does not assume the glued
polytope remains acute.
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/new_dimension_h20_facet_gluing_20261013.json'
OUTPUT=ROOT/'data/enumerations/new_dimension_h20_six_patch_boundary_20261013.json'

def oriented_primitive(row):
    row=tuple(map(int,row));d=gcd(*row)
    assert d>0
    if next(x for x in row if x)<0:d=-d
    return tuple(x//d for x in row)

def plane_key(a,b):
    # Two scaled rows of the exact reduced row echelon form. No divisions
    # occur until primitive integral normalization; row orientation is fixed.
    p=next(i for i,(x,y) in enumerate(zip(a,b)) if x or y)
    second=int(a[p])*b-int(b[p])*a
    q=next(i for i,x in enumerate(second) if x)
    first=int(b[q])*a-int(a[q])*b
    return oriented_primitive(first),oriented_primitive(second)

def boundary(archive,Q,names):
    patch=archive['six_tile_local_angle_repair'];tiles=patch['tiles']
    remove={name:set() for name in names}
    for row in patch['intended_seams']:
        if row['tile_a'] not in names or row['tile_b'] not in names:continue
        remove[row['tile_a']].add(row['wall_a']);remove[row['tile_b']].add(row['wall_b'])
    normal_sources={}
    for name in names:
        tile=tiles[name]
        for label,row in zip(tile['wall_labels_one_based'],tile['primitive_outward_normals']):
            if label not in remove[name]:normal_sources.setdefault(tuple(row),[]).append([name,label])
    normal_rows=sorted(normal_sources);V=np.array(normal_rows,dtype=np.int64)
    normal_gram=V@Q@V.T;planes={};tile_ridges=0
    for name in names:
        tile=tiles[name];A=np.array(tile['primitive_outward_normals'],dtype=np.int64);H=A@Q@A.T
        assert all(H[i,j]<=0 for i,j in combinations(range(len(A)),2))
        for i,j in combinations(range(len(A)),2):
            if int(H[i,j])**2>=int(H[i,i])*int(H[j,j]):continue
            tile_ridges+=1
            key=plane_key(A[i],A[j])
            planes.setdefault(key,[name,tile['wall_labels_one_based'][i],tile['wall_labels_one_based'][j]])
    angles=Counter();examples={};ridges=[];unrealized=Counter();boundary_planes={}
    for i,j in combinations(range(len(V)),2):
        ip=int(normal_gram[i,j]);den=int(normal_gram[i,i])*int(normal_gram[j,j]);num=ip*ip
        if num>=den:continue
        key=(('obtuse' if ip>0 else 'right' if ip==0 else 'acute'),str(F(num,den)))
        normal_plane=plane_key(V[i],V[j]);witness=planes.get(normal_plane)
        if witness is None:unrealized[str(key)]+=1;continue
        boundary_planes.setdefault(normal_plane,set()).update((i,j))
        angles[key]+=1;ridges.append([i,j,key[0],key[1],witness])
        examples.setdefault(key,{'boundary_indices_zero_based':[i,j],
                                'boundary_sources':[normal_sources[normal_rows[i]],normal_sources[normal_rows[j]]],
                                'original_tile_ridge_witness':witness,'normal_pairing':ip,
                                'normal_norms':[int(normal_gram[i,i]),int(normal_gram[j,j])]})
    assert all(len(indices)==2 for indices in boundary_planes.values())
    return {'tiles':list(names),'facet_count':len(V),'primitive_outward_normals':V.tolist(),
            'boundary_wall_sources':[normal_sources[row] for row in normal_rows],
            'genuine_tile_ridges_counted_with_multiplicity':tile_ridges,
            'distinct_genuine_tile_ridge_planes':len(planes),
            'each_boundary_ridge_has_exactly_two_facets':True,
            'actual_boundary_ridge_count':sum(angles.values()),
            'signed_angle_census':[{'kind':k[0],'cosine_squared':k[1],'count':n,'witness':examples[k]} for k,n in sorted(angles.items())],
            'raw_positive_definite_pairs_not_realized':dict(unrealized),
            'all_actual_ridges':ridges,
            'recognition':'A boundary-normal pair is counted only when its exact normal plane equals that of a genuine ridge in an original acute tile.'}

def collar(Q,W,G):
    f=0;neighbors=[j for j in range(210) if G[j,j]==2 and G[f,j]==-1]
    pairs=list(combinations(neighbors,2));normals=np.array([W[f]+W[i]+W[j] for i,j in pairs],dtype=np.int64)
    gram=normals@Q@normals.T
    for a,p in enumerate(pairs):
        for b,q in enumerate(pairs):
            expected=2 if a==b else 0 if set(p)&set(q) else -2
            assert int(gram[a,b])==expected
    a=0;b=next(i for i,p in enumerate(pairs) if not set(p)&set(pairs[a]))
    I=np.eye(22,dtype=np.int64)
    reflections=[I-np.outer(v,v@Q) for v in normals]
    T=reflections[a]@reflections[b];N=sp.Matrix(T)-sp.eye(22)
    assert N!=sp.zeros(22) and N**3==sp.zeros(22)
    assert np.array_equal(T@W[f],W[f]) and np.array_equal(T.T@Q@T,Q)
    return {'central_facet_one_based':1,'core_neighbors_one_based':[i+1 for i in neighbors],
            'pairs_one_based':[[i+1,j+1] for i,j in pairs],'normals':normals.tolist(),'gram':gram.tolist(),
            'reflection_matrices':[A.tolist() for A in reflections],
            'disjoint_pair_indices_zero_based':[a,b],'parabolic_product':T.tolist(),
            'product_minus_identity_rank':int(N.rank()),'product_minus_identity_square_rank':int((N**2).rank()),
            'product_minus_identity_cube_zero':True,'product_has_infinite_order':True,
            'scope':'Closure under every collar reflection contains this infinite cyclic subgroup and is infinite. This does not rule out selective finite patches or a finite chamber for the infinite reflection group.'}

def dehn_balance(G,f=0):
    short=[i for i in range(210) if G[i,i]==2];long=[i for i in range(210) if G[i,i]==4]
    rows=[];counts=Counter()
    H=2*G-np.outer(G[:,f],G[f,:])
    for i,j in combinations([i for i in range(210) if i!=f],2):
        r=F(int(H[i,j])**2,int(H[i,i])*int(H[j,j]))
        if r not in {F(1,9),F(1,3),F(2,3)}:continue
        if G[i,i]==G[j,j]==2:
            kind='A3_center_alpha' if r==F(1,9) else 'A3_endpoint_beta'
            assert r in {F(1,9),F(1,3)}
        else:
            assert sorted([int(G[f,f]),int(G[i,i]),int(G[j,j])])==[2,2,4]
            kind='B3_short_center_beta' if r==F(1,3) else 'B3_short_endpoint_gamma'
            assert r in {F(1,3),F(2,3)}
        counts[kind]+=1
        rows.append({'ambient_face_one_based':[f+1,i+1,j+1],'type':kind,
                     'gram':G[np.ix_([f,i,j],[f,i,j])].tolist(),'cosine_squared':str(r)})
    assert counts=={'A3_center_alpha':10,'A3_endpoint_beta':20,'B3_short_center_beta':240,'B3_short_endpoint_gamma':240}
    a3=[]
    for triple in combinations(short,3):
        if sorted(int(G[i,j]) for i,j in combinations(triple,2))==[-1,-1,0]:a3.append([i+1 for i in triple])
    b3=[]
    for i,j in combinations(short,2):
        if G[i,j]!=-1:continue
        for k in long:
            if sorted([int(G[i,k]),int(G[j,k])])==[-2,0]:b3.append([i+1,j+1,k+1])
    assert len(a3)==420 and len(b3)==10080
    return {'fixed_short_facet_one_based':f+1,'irrational_ridge_types':dict(counts),'fixed_facet_triples':rows,
            'all_A3_faces_one_based':a3,'all_relevant_B3_faces_one_based':b3,
            'A3_face_count':len(a3),'B3_face_count':len(b3),
            'A3_incidence_cancellation':'One alpha and two beta incidences have gamma coefficients2,-1,-1 on the identical H18 face.',
            'B3_incidence_cancellation':'The two short-facet incidences are beta andgamma, with coefficients-1,+1 on the identical H18 face.',
            'conclusion':'For the codimension-two Dehn component, summing over42 congruent short facets gives42D=0; the rational target is torsion-free, henceD=0. This also holds for the volume-valued projection. No assertion about other generalized Dehn components or compatible gluing.'}

def defect_cover(row):
    import networkx as nx
    from scipy.optimize import milp,LinearConstraint,Bounds
    edges=[(i,j) for i,j,kind,r,witness in row['all_actual_ridges'] if kind=='obtuse' or F(r) not in {F(0),F(1,4),F(1,2),F(3,4)}]
    graph=nx.Graph();graph.add_nodes_from(range(row['facet_count']));graph.add_edges_from(edges)
    matching=sorted(tuple(sorted(e)) for e in nx.max_weight_matching(graph,maxcardinality=True))
    assert len({i for edge in matching for i in edge})==2*len(matching)
    assert all(graph.has_edge(i,j) for i,j in matching)
    triangles=[sorted(c) for c in nx.connected_components(graph) if len(c)==3 and graph.subgraph(c).number_of_edges()==3]
    triangle_vertices={i for c in triangles for i in c}
    residual_matching=[e for e in matching if not set(e)&triangle_vertices]
    strengthened_lower=len(residual_matching)+2*len(triangles)
    A=np.zeros((len(edges),row['facet_count']))
    for k,(i,j) in enumerate(edges):A[k,i]=A[k,j]=1
    result=milp(np.ones(row['facet_count']),integrality=np.ones(row['facet_count']),bounds=Bounds(0,1),constraints=LinearConstraint(A,1,np.inf),options={'time_limit':15})
    assert result.x is not None
    cover=[i for i,v in enumerate(result.x) if v>0.5];C=set(cover)
    assert all(i in C or j in C for i,j in edges)
    return {'bad_ridges':len(edges),'disjoint_matching_lower_bound':len(matching),
            'matching_boundary_indices_zero_based':[list(e) for e in matching],
            'cover_upper_bound':len(cover),'cover_boundary_indices_zero_based':cover,
            'disjoint_triangle_components':triangles,
            'disjoint_edge_and_triangle_lower_bound':strengthened_lower,
            'exact_optimum_certified':strengthened_lower==len(cover),
            'scope':'MILP proposes only a cover; exact incidence verifies it. Disjoint edges and triangle components supply the rigorous lower bound.'}

def main():
    archive=json.loads(SOURCE.read_text())
    raw=json.loads((ROOT/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json').read_text())['coxeter_input']
    Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64);G=W@Q@W.T
    out={'source':str(SOURCE.relative_to(ROOT)),'two_tile':boundary(archive,Q,('A','DA')),
         'three_tile':boundary(archive,Q,('A','B','C'))}
    print('Three-tile census:',[(r['kind'],r['cosine_squared'],r['count']) for r in out['three_tile']['signed_angle_census']],flush=True)
    out['six_tile']=boundary(archive,Q,('A','DA','B','DB','C','DC'))
    print('Six-tile census:',[(r['kind'],r['cosine_squared'],r['count']) for r in out['six_tile']['signed_angle_census']],flush=True)
    out['ten_mirror_collar']=collar(Q,W,G)
    out['short_facet_dehn_balance']=dehn_balance(G)
    for name in ('two_tile','three_tile','six_tile'):
        out[name]['defect_cover_bounds']=defect_cover(out[name])
    OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print('Saved '+str(OUTPUT),flush=True)

if __name__=='__main__':main()
