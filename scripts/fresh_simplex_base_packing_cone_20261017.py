#!/usr/bin/env python3
"""Exact full disjoint-base packing cones on two through five vertices.

No optimizer or project imports. The completeness proof is in the report:
the five-vertex cover dual reduces to a fractional stable-set polytope on
the Petersen graph, whose vertices are half-integral.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_simplex_base_packing_cone_20261017'


def rank(rows):
    basis = {}
    for row in rows:
        row = list(map(F, row))
        for p, b in sorted(basis.items()):
            if row[p]:
                z = row[p]
                row = [x-z*y for x,y in zip(row,b)]
        pivot = next((j for j,x in enumerate(row) if x), None)
        if pivot is not None:
            z = row[pivot]
            basis[pivot] = [x/z for x in row]
    return len(basis)


def is_vertex(v, adjacent):
    # At fractional coordinates, a tight-edge component must contain an odd
    # cycle; a bipartite component admits a two-sided alternating perturbation.
    half = {i for i,x in enumerate(v) if x == 1}
    left = set(half)
    while left:
        start = min(left)
        color = {start:0}
        todo = [start]
        odd = False
        for i in todo:
            for j in sorted(adjacent[i] & half):
                if j not in color:
                    color[j] = 1-color[i]
                    todo.append(j)
                elif color[j] == color[i]:
                    odd = True
        left -= color.keys()
        if not odd:
            return False
    return True


def petersen_vertices():
    pairs = [m for m in range(32) if m.bit_count() == 2]
    edges = [(i,j) for i,j in combinations(range(10),2)
             if not pairs[i] & pairs[j]]
    adjacent = {i:set() for i in range(10)}
    for i,j in edges:
        adjacent[i].add(j)
        adjacent[j].add(i)
    feasible = []
    vertices = []
    for v in product(range(3),repeat=10):
        if any(v[i]+v[j] > 2 for i,j in edges):
            continue
        feasible.append(v)
        if is_vertex(v,adjacent):
            active = []
            for i,x in enumerate(v):
                if x in (0,2):
                    row = [0]*10
                    row[i] = 1
                    active.append(row)
            for i,j in edges:
                if v[i]+v[j] == 2:
                    row = [0]*10
                    row[i] = row[j] = 1
                    active.append(row)
            assert rank(active) == 10
            vertices.append(v)
    assert len(feasible) == 1959 and len(vertices) == 294
    hist = Counter(v.count(1) for v in vertices)
    assert dict(hist) == {0:76,5:12,6:60,7:90,8:45,9:10,10:1}
    return pairs, vertices, dict(hist)


def cone(n, pairs5, vertices5):
    masks = [m for m in range(1<<n) if m.bit_count() >= 2]
    packings = [()] + [(m,) for m in masks]
    packings += list((a,b) for a,b in combinations(masks,2) if not a&b)
    assert n <= 5  # No three nontrivial disjoint bases fit.
    rows = []
    if n <= 3:
        rows = [[2]*len(masks)]
    elif n == 4:
        complementary = [(m,15^m) for m in masks if m.bit_count()==2 and m < (15^m)]
        assert len(complementary) == 3
        for choices in product((0,1),repeat=3):
            chosen = {ab[k] for ab,k in zip(complementary,choices)}
            rows.append([2*int(m.bit_count()>=3 or m in chosen) for m in masks])
    else:
        for v in vertices5:
            d = dict(zip(pairs5,v))
            rows.append([d[m] if m.bit_count()==2 else
                         2-d[31^m] if m.bit_count()==3 else 2 for m in masks])
    incidence = [[int(m in packing) for m in masks] for packing in packings]
    facet_ranks = []
    for coefficients in rows:
        assert all(0 <= x <= 2 for x in coefficients)
        loads = [sum(a*b for a,b in zip(coefficients,p)) for p in incidence]
        assert max(loads) == 2
        tight = [p for p,load in zip(incidence,loads) if load == 2]
        # These vectors lie on a nonzero affine hyperplane. Full linear
        # rank therefore means affine dimension len(masks)-1.
        r = rank(tight)
        assert r == len(masks)
        facet_ranks.append(r-1)
    assert len({tuple(r) for r in rows}) == len(rows)
    return dict(vertices=n, base_masks=masks, integer_packings=packings,
                doubled_capacity=2, doubled_inequality_coefficients=rows,
                facet_affine_ranks=facet_ranks,
                inequality_count=len(rows), packing_count=len(packings))


def run():
    pairs,vertices,hist = petersen_vertices()
    cones = [cone(n,pairs,vertices) for n in range(2,6)]
    assert [r['inequality_count'] for r in cones] == [1,1,8,294]
    assert [r['packing_count'] for r in cones] == [2,5,15,52]
    blocks = [0b00011,0b00101,0b01110,0b11001]
    assert all(a&b for a,b in combinations(blocks,2))
    floors = []
    for T in range(1,32):
        for r in range(1,T.bit_count()+1):
            load = sum((F((B&T).bit_count()//r,3) for B in blocks),F())
            capacity = T.bit_count()//r
            assert load <= capacity
            floors.append([T,r,str(load),capacity])
    five = cones[-1]
    witness_loads = [sum((F(coef[five['base_masks'].index(B)],6) for B in blocks),F())
                    for coef in five['doubled_inequality_coefficients']]
    assert max(witness_loads) == F(4,3)
    violating = [i for i,x in enumerate(witness_loads) if x>1]
    assert violating
    own = f'scripts/{TAG}.py'
    return dict(status='PASS_EXACT_COMPLETE_SIMPLEX_BASE_PACKING_CONES',
        sources={own:sha256((ROOT/own).read_bytes()).hexdigest()},
        cones=cones, petersen_feasible_half_grid_count=1959,
        petersen_vertex_half_coordinate_histogram=sorted(hist.items()),
        witness=dict(vertices=5,base_masks=blocks,mass_each='1/3',
            subset_floor_checks=floors,all_subset_floor_checks=len(floors),
            pairwise_intersecting=True,clique_load='4/3',clique_capacity=1,
            full_cone_violated_rows=violating,maximum_full_cone_load='4/3'),
        completeness='Nonnegative base marginals are coverable by mixtures of disjoint packings with total mass<=c iff every listed inequality has load<=c. Five-vertex LP dual reduces to the ten-pair Petersen fractional stable-set polytope; its vertices are half-integral by the tight-edge perturbation proof. All294 vertices are enumerated exactly. No optimizer is used.',
        scope='Complete abstract necessary disjoint-base packing domain on at most five simplex vertices. This does not assert geometric realizability, sufficient Coxeter compatibility, a violation at any repository witness, or a dimension bound.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify',action='store_true')
    args = p.parse_args()
    result = json.loads(json.dumps(run()))
    out = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],[r['inequality_count'] for r in result['cones']],
          'witness4/3>1;',result['witness']['all_subset_floor_checks'],'floor checks pass',flush=True)
