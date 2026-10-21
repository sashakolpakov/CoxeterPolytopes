#!/usr/bin/env python3
"""Full-norm bounded Vinberg search on <-1-sqrt2> + E8 + <1>.

The integral lattice differs at its chosen positive block from the frozen
F4^2 host. No global inequivalence or novelty follows from that observation.
"""
import argparse
import json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, groupby
from math import comb
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver as S, SQRT2, silver_determinant, silver_inertia
from coxeter_polytopes.search.silver_path_host import bounded_silver_integers
import fresh_silver_f4_h9_20261013 as H

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_silver_e8_h9_20261015'
Z, O = S(0), S(1)
EPS = 1+SQRT2
ETA = 2+SQRT2
EDGES = {(i, i+1) for i in range(6)} | {(2, 7)}
E8 = tuple(tuple(S(2) if i==j else -O if tuple(sorted((i, j))) in EDGES else Z
                 for j in range(8)) for i in range(8))
EI = H.inverse(E8)
SCHUR = tuple(H.inverse(tuple(row[:k] for row in E8[:k])) for k in range(1, 9))
GRAM = tuple(tuple(-EPS if i==j==0 else O if i==j==9 else E8[i-1][j-1]
                   if 1<=i<=8 and 1<=j<=8 else Z for j in range(10)) for i in range(10))
NORMS = (O, ETA, S(2))
assert silver_determinant(E8)==1 and silver_determinant(GRAM)==-EPS
assert all(v.is_integral and v.sign()>0 for row in EI for v in row)
assert silver_inertia(GRAM).positive==9 and silver_inertia(GRAM).negative==1
assert silver_inertia(GRAM, embedding=1).positive==10


@lru_cache(None)
def blocks(target, q):
    alphabet = tuple(x for x in bounded_silver_integers((2*target, 2*target))
                     if x.sign()>=0 and (2*x/q).is_integral)
    found = []
    visits = 0
    def rec(prefix):
        nonlocal visits
        visits += 1
        k = len(prefix)
        if k:
            lower = H.pair(prefix, SCHUR[k-1], prefix)
            if any((target-lower).sign(e)<0 for e in (0, 1)):
                return
            principal = tuple(row[:k] for row in EI[:k])
            if (target-H.pair(prefix, principal, prefix)).sign()<0:
                return
        if k==8:
            y = H.mv(EI, prefix)
            assert all(v.is_integral for v in y)
            found.append((lower, y))
            return
        for value in alphabet:
            rec((*prefix, value))
    rec(())
    return tuple(found), visits


def candidates(x, q):
    target = q+EPS*x*x
    choices, visits = blocks(target, q)
    scalars = tuple(z for z in bounded_silver_integers((target, target))
                    if z.sign()>=0 and (2*z/q).is_integral)
    by_square = {z*z: z for z in scalars}
    out = []
    for cost, y in choices:
        if target-cost not in by_square:
            continue
        r = (x, *y, by_square[target-cost])
        assert H.pair(r, GRAM, r)==q
        assert all((2*v/q).is_integral for v in H.mv(GRAM, r))
        out.append(r)
    return out, visits


def run(height):
    zero = []
    for q in NORMS:
        roots, visits = candidates(Z, q)
        zero.append(dict(norm=q.to_json(), visits=visits, roots=[H.enc(r) for r in roots]))
    assert [len(row['roots']) for row in zero] == [1, 0, 2]
    walls = [tuple(-O if i==j else Z for i in range(10)) for j in range(1, 10)]
    keys = {H.ray_key(r) for r in walls}
    controls = sorted((EPS*x*x/q, q, x) for q in NORMS
                      for x in bounded_silver_integers((height*q/EPS, -q/EPS))
                      if x.sign()>0 and (2*EPS*x/q).is_integral)
    levels = []
    print('zero census', [len(r['roots']) for r in zero], 'controls', len(controls), flush=True)
    for ht, packet in groupby(controls, key=lambda row: row[0]):
        accepted = []
        for _, q, x in packet:
            roots, visits = candidates(x, q)
            kept, duplicates = 0, 0
            for r in roots:
                if any(H.pair(r, GRAM, w).sign()>0 for w in walls):
                    continue
                key = H.ray_key(r)
                if key in keys:
                    duplicates += 1
                    continue
                keys.add(key)
                accepted.append(r)
                kept += 1
            levels.append(dict(height=ht.to_json(), norm=q.to_json(), time_coordinate=x.to_json(),
                               recursion_visits=visits, candidates=len(roots), roots=kept,
                               parallel_duplicates=duplicates))
        assert all(H.pair(r, GRAM, s).sign()<=0 for i, r in enumerate(accepted) for s in accepted[:i])
        walls.extend(accepted)
        print('height', str(ht), 'new', len(accepted), 'walls', len(walls), flush=True)
    spatial = H.mv(EI, (O,)*8)
    small = S(F(1, 2))
    while True:
        inside = (O, *(small*v for v in spatial), small)
        if H.pair(inside, GRAM, inside).sign()<0 and all(H.pair(inside, GRAM, r).sign()<0 for r in walls):
            break
        small /= 2
    angles = H.angle_gate(walls, GRAM)
    geometric = H.geometry(walls, GRAM, inside, 9, exhaust=comb(len(walls), 9)<=5000)
    paths = [f'scripts/{TAG}.py', 'scripts/fresh_silver_f4_h9_20261013.py',
             'src/coxeter_polytopes/exact/silver.py', 'src/coxeter_polytopes/search/silver_path_host.py']
    return dict(status='PASS_COMPLETE_BOUNDED_FULL_NORM_VINBERG_PREFIX',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                field='Q(sqrt2)', dimension=9, host='<-1-sqrt2>+E8+<1>', gram=[H.enc(r) for r in GRAM],
                determinant=(-EPS).to_json(), norm_classes=H.enc(NORMS), height=height.to_json(),
                zero_height_dominant_census=zero, time_controls=len(controls), levels=levels,
                all_walls=[H.enc(r) for r in walls], wall_count=len(walls), angle_counts=angles,
                strict_interior=H.enc(inside), strict_interior_norm=H.pair(inside, GRAM, inside).to_json(),
                geometry=geometric,
                norm_completeness='O=Z[sqrt2] is a PID and det(Gram) is a unit. Every primitive integral reflection norm divides2 up to units. Totally positive units are squares and the totally positive ramified generator is2+sqrt2, giving classes1,2+sqrt2,2. Nonprimitive parallel representatives are deduplicated.',
                enumeration='For dominant E8 covector w, each w_i>=0 physically and w_i^2<=2target in both embeddings. Leading-E8 inverse bounds minimize the remaining unconstrained coordinates and are valid in both embeddings. The entire inverse E8 matrix is entrywise positive, giving the additional physical prefix bound. Every resulting vector, norm, integral reflection and chronological halfspace condition is checked exactly.',
                scope='Complete finite height prefix in all primitive norm classes, with exact geometry only to its displayed extent. No claim that the full reflection chamber has finite volume, no global nonreflectivity claim, and no lattice-isometry or worldwide novelty claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--height', type=int, default=30)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = run(S(args.height))
    path = ROOT/'results'/f'{TAG}_h{args.height}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['wall_count'], out['geometry']['status'], flush=True)
