#!/usr/bin/env python3
"""Exact bounded Vinberg search on new cubic-seven <-a>+E8+<1> hosts.

Numerical LP is used only to suggest a spacelike witness; every reported
escape is checked in the cubic field, against every wall.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import time

from coxeter_polytopes.exact.cubic7 import Cubic7 as C, CUBIC7_ZERO as Z, CUBIC7_ONE as O
from coxeter_polytopes.search.cubic7_e8_cubes import SUPPORTED_E8_CUBE_SEEDS, _control_levels
from coxeter_polytopes.search.cubic7_root_blocks import (
    _block_data, _bounded_for_caps, _mat_vec, enumerate_dominant_block_weights,
)
from coxeter_polytopes.search.cubic7_bna1 import _kernel


def gram_for(a):
    block, _, _ = _block_data('E8')
    return tuple(tuple(-a if i == j == 0 else O if i == j == 9 else
                       C(block[i-1][j-1]) if 1 <= i <= 8 and 1 <= j <= 8 else Z
                       for j in range(10)) for i in range(10))


def pair(x, gram, y):
    return sum((x[i]*gram[i][j]*y[j] for i in range(10)
                for j in range(10) if gram[i][j]), Z)


def candidates(x, target, q):
    block = enumerate_dominant_block_weights('E8', target, q)
    scalar = {}
    for u in _bounded_for_caps((target, target, target)):
        if u.sign() >= 0 and (2*u/q).is_integral and all(
                (target-u*u).sign(p) >= 0 for p in range(3)):
            scalar.setdefault(u*u, []).append(u)
    vectors = tuple((x, *v.coordinates, u) for v in block.vectors
                    for u in scalar.get(target-v.square, ()))
    return vectors, sum(block.visited_by_depth)


def encode_vector(v):
    return [x.to_json() for x in v]


def exact_escape(walls, gram):
    import numpy as np
    from scipy.optimize import linprog
    theta = 2*np.cos(2*np.pi/7)
    real = lambda x: float(x.a)+float(x.b)*theta+float(x.c)*theta**2
    rows = [tuple(sum((gram[i][j]*w[j] for j in range(10)), Z)
                  for i in range(10)) for w in walls]
    # The artificial box only supplies additional possible support equations.
    # An accepted point is verified against the original chamber, not the box.
    for j in range(1, 10):
        for sign in (-1, 1):
            rows.append(tuple(C(-100) if i == 0 else C(sign) if i == j else Z
                              for i in range(10)))
    f = np.array([[real(x) for x in row] for row in rows])
    rng = np.random.default_rng(20261013)
    for sample in range(100):
        result = linprog(rng.normal(size=9), A_ub=f[:,1:], b_ub=-f[:,0],
                         bounds=[(None,None)]*9, method='highs')
        if not result.success:
            continue
        active = np.flatnonzero(abs(f[:,1:]@result.x+f[:,0]) < 1e-7)
        if len(active) < 9:
            continue
        v = _kernel(tuple(rows[int(i)] for i in active))
        if v is None or v[0] == Z:
            continue
        v = tuple(x/v[0] for x in v)
        norm = pair(v, gram, v)
        products = tuple(pair(v, gram, w) for w in walls)
        if norm.sign() > 0 and all(x.sign() <= 0 for x in products):
            return {'verdict':'EXACT_SPACELIKE_ESCAPE', 'vector':encode_vector(v),
                    'norm':norm.to_json(), 'wall_pairings':encode_vector(products),
                    'lp_sample':sample,
                    'scope':'The bounded prefix has infinite volume. No conclusion about the full reflection group.'}
    return {'verdict':'UNDECIDED', 'scope':'Failure to find an escape does not certify finite volume.'}


def run(seed_index, height, output):
    start = time.monotonic()
    seed = SUPPORTED_E8_CUBE_SEEDS[seed_index]
    a = seed.timelike_coefficient
    gram = gram_for(a)
    zero_census = []
    for q in seed.canonical_mirror_norms:
        vectors, nodes = candidates(Z, q, q)
        zero_census.extend((q, v) for v in vectors)
    # Dominant representatives of the E8 and A1 root orbits.
    assert len(zero_census) == 2
    assert sorted(q for q, _ in zero_census) == [C(1), C(2)]
    assert all((q == 1 and v[9] == 1 and all(x == 0 for x in v[1:9])) or
               (q == 2 and v[9] == 0) for q, v in zero_census)
    walls = [tuple(-O if i == j else Z for i in range(10)) for j in range(1,10)]
    levels = []
    payload = {'host':f'<-({a})>+E8+<1>', 'dimension':9, 'seed_index':seed_index,
               'timelike_coefficient':a.to_json(), 'height_bound':height,
               'norm_classes':[q.to_json() for q in seed.canonical_mirror_norms],
               'zero_height_dominant_census':[{'norm':q.to_json(),'vector':encode_vector(v)} for q,v in zero_census],
               'scope':'Bounded exact search. A nonclosed prefix is not a nonreflectivity proof.',
               'levels':levels}
    controls = _control_levels(seed, height)
    print(json.dumps({'phase':'start','seed':seed_index,'control_levels':len(controls)}),flush=True)
    for index, level in enumerate(controls):
        vectors, nodes = candidates(level.control, level.target, level.norm)
        accepted = [v for v in vectors if all(pair(v,gram,w).sign() <= 0 for w in walls)]
        assert all(pair(v,gram,w).sign() <= 0 for i,v in enumerate(accepted) for w in accepted[:i])
        for v in accepted:
            assert pair(v,gram,v) == level.norm
            assert all((2*sum((gram[i][j]*v[j] for j in range(10)),Z)/level.norm).is_integral
                       for i in range(10))
        walls.extend(accepted)
        row = {'index':index,'height':level.height.to_json(),'norm':level.norm.to_json(),
               'control':level.control.to_json(),'candidates':len(vectors),'new_walls':len(accepted),
               'walls':len(walls),'nodes':nodes}
        levels.append(row)
        payload.update({'all_walls':[encode_vector(v) for v in walls],
                        'wall_count':len(walls),'completed_control_levels':len(levels),
                        'total_control_levels':len(controls),'search_complete':False})
        output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
        print(json.dumps(row),flush=True)
    _, inverse, _ = _block_data('E8')
    spatial = _mat_vec(inverse, (O,)*8)
    epsilon = C(Fraction(1,2))
    while True:
        interior = (O,*(epsilon*x for x in spatial),epsilon)
        if pair(interior,gram,interior).sign() < 0 and all(pair(interior,gram,w).sign() < 0 for w in walls):
            break
        epsilon /= 2
    payload.update({'all_walls':[encode_vector(v) for v in walls], 'wall_count':len(walls),
                    'strict_interior':encode_vector(interior),'search_complete':True,
                    'geometry':exact_escape(walls,gram),'seconds':round(time.monotonic()-start,3)})
    output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'phase':'done','seed':seed_index,'walls':len(walls),
                      'geometry':payload['geometry']['verdict'],'seconds':payload['seconds']}),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=int,choices=range(3),required=True)
    p.add_argument('--height',type=int,default=10)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    run(args.seed,args.height,args.output)
