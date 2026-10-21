#!/usr/bin/env python3
"""Bounded Coxeter-compatible wall subsets of the exact six-tile H20 patch."""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as F
import json
from math import gcd
from pathlib import Path
import random
import numpy as np
from scipy.optimize import milp,LinearConstraint,Bounds
from scipy.sparse import coo_matrix

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/new_dimension_h20_six_patch_boundary_20261013.json'
OUTPUT=ROOT/'data/enumerations/new_dimension_h20_compatible_subsets_20261013.json'

def primitive(row):
    row=tuple(map(int,row));d=gcd(*row);assert d>0
    return tuple(x//d for x in row)

def cyclic_nonintegrality_witness(ids,G):
    """Intrinsic obstruction, unlike nonintegrality in one chosen lattice."""
    for position,i in enumerate(ids):
        for j in ids[position+1:]:
            value=F(4*int(G[i,j])**2,int(G[i,i])*int(G[j,j]))
            if value.denominator!=1:
                return {'indices_zero_based':[i,j], 'cyclic_product':str(value),
                        'pairing':int(G[i,j]),'norms':[int(G[i,i]),int(G[j,j])]}
    return None

def main():
    archive=json.loads(SOURCE.read_text());V=np.array(archive['six_tile']['primitive_outward_normals'],dtype=np.int64)
    raw=json.loads((ROOT/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json').read_text())['coxeter_input']
    Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64);G=V@Q@V.T
    n=len(V);neighbors=[set() for _ in range(n)];edges=[]
    for i in range(n):
        for j in range(i):
            ip=int(G[i,j]);num=ip*ip;den=int(G[i,i])*int(G[j,j])
            if ip>0 or (num<den and F(num,den) not in {F(0),F(1,4),F(1,2),F(3,4)}):
                edges.append((i,j));neighbors[i].add(j);neighbors[j].add(i)
    # Saturated integral basis of f^perp: x_21=2*x_0.
    H=np.zeros((22,21),dtype=np.int64)
    for i in range(20):H[i+1,i]=1
    H[0,20]=1;H[21,20]=2
    nonintegral=[any(int(x)%int(G[i,i]) for x in 2*(V[i]@Q@H)) for i in range(n)]
    index={tuple(v):i for i,v in enumerate(V.tolist())}
    base=[]
    for w in W:
        if np.array_equal(w,W[0]):continue
        norm=int(w@Q@w);coupling=int(w@Q@W[0])
        if norm==2 and coupling==-1:continue
        v=primitive(2*w-coupling*W[0]);assert v in index;base.append(index[v])
    base=sorted(set(base));assert len(base)==204
    assert not any(j in neighbors[i] for k,i in enumerate(base) for j in base[:k])
    patch=json.loads((ROOT/'data/enumerations/new_dimension_h20_facet_gluing_20261013.json').read_text())['six_tile_local_angle_repair']
    D=np.array(patch['tiles']['DA']['transport_numerator'],dtype=np.int64)
    permutation=[index[primitive(D@v)] for v in V]
    assert all(permutation[permutation[i]]==i for i in range(n))
    out={'source':str(SOURCE.relative_to(ROOT)),
         'scope':'Fixed bounded selection round: 32 seeded greedy maximal independent sets, base204 extensions, and6 MILP objectives with15second limits. Pair compatibility is exact; finite volume is checked separately.',
         'wall_order':'six_tile.primitive_outward_normals in source',
         'wall_count':n,'complete_conflict_edges':[[i,j] for i,j in edges],
         'complete_conflict_edge_count':len(edges),'nonintegral_reflection_indices_zero_based':[i for i,b in enumerate(nonintegral) if b],
         'base204_indices_zero_based':base,'reflection_D_permutation':permutation,'candidates':[]}
    seen=set()
    def save(ids,label,metadata=None,complete=True):
        ids=sorted(set(ids));key=tuple(ids)
        assert not any(j in neighbors[i] for k,i in enumerate(ids) for j in ids[:k])
        # Complete each selected set deterministically: maximality is useful
        # here, although maximum cardinality is neither assumed nor needed.
        if complete:
            for i in range(n):
                if i not in ids and not any(j in neighbors[i] for j in ids):ids.append(i)
        ids=sorted(ids);key=tuple(ids)
        maximal=all(i in ids or any(j in neighbors[i] for j in ids) for i in range(n))
        assert maximal or not complete
        if key in seen:return
        seen.add(key)
        out['candidates'].append({'id':f'compatible_{len(out["candidates"]):03d}','selection':label,
                                  'indices_zero_based':ids,'walls':len(ids),
                                  'nonintegral_reflections':sum(nonintegral[i] for i in ids),
                                  'D_invariant':set(ids)=={permutation[i] for i in ids},
                                  'length_two_Vinberg_nonintegrality_witness':cyclic_nonintegrality_witness(ids,G),
                                  'exact_pair_compatibility':True,'maximal_under_single_wall_addition':maximal,
                                  'metadata':metadata or {}})
    def greedy(seed,forced=()):
        rng=random.Random(seed);chosen=list(forced);available=set(range(n))-set(chosen)
        for i in chosen:available-=neighbors[i]
        while available:
            scores={i:(1+0.5*rng.random()+((seed%3)*0.3 if nonintegral[i] else 0))/(1+len(neighbors[i]&available)) for i in available}
            i=max(available,key=lambda i:(scores[i],-i));chosen.append(i);available.remove(i);available-=neighbors[i]
        return chosen
    save(greedy(0,base),'base204_greedy_extension')
    for seed in range(32):save(greedy(seed),f'greedy_seed_{seed}')
    for seed in (1,7,19):save(greedy(seed,base),f'base204_extension_seed_{seed}')
    OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print('Initial candidates:',[(r['id'],r['walls'],r['nonintegral_reflections']) for r in out['candidates']],flush=True)
    # Three unrestricted objectives and three objectives preserving D.
    for trial in range(6):
        symmetric=trial>=3;rng=random.Random(1000+trial)
        weights=np.array([1+(trial%3)*0.5*nonintegral[i]+0.03*rng.random() for i in range(n)])
        entries=[];rows=[];cols=[];lower=[];upper=[]
        for k,(i,j) in enumerate(edges):rows.extend([k,k]);cols.extend([i,j]);entries.extend([1.,1.]);lower.append(-np.inf);upper.append(1.)
        if symmetric:
            for i,j in enumerate(permutation):
                if i>=j:continue
                k=len(lower);rows.extend([k,k]);cols.extend([i,j]);entries.extend([1.,-1.]);lower.append(0.);upper.append(0.)
        A=coo_matrix((entries,(rows,cols)),shape=(len(lower),n)).tocsr()
        r=milp(-weights,integrality=np.ones(n),bounds=Bounds(0,1),constraints=LinearConstraint(A,lower,upper),options={'time_limit':15})
        if r.x is not None:
            raw_ids=[i for i,x in enumerate(r.x) if x>.5]
            save(raw_ids,f'milp_{trial}',{'status':int(r.status),'D_constraint':symmetric,'time_limit_seconds':15})
            if symmetric:save(raw_ids,f'milp_{trial}_D_invariant_before_completion',{'status':int(r.status),'D_constraint':True,'time_limit_seconds':15},complete=False)
        OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print('MILP',trial,'status',r.status,'candidates',len(out['candidates']),flush=True)
    print('Saved '+str(OUTPUT),flush=True)

if __name__=='__main__':main()
