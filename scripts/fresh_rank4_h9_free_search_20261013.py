#!/usr/bin/env python3
"""Bounded free-ratio search in the surviving two-weight golden H9 slice.

The integer-array filters are exact with explicit overflow guards. Every
candidate square and wall pairing is then replayed in field arithmetic.
"""
from __future__ import annotations
import argparse
import json
from math import isqrt
import time
import numpy as np
from scripts.fresh_multibasepoint_h9_merge import (
    ROOT, G, read, matrix, inverse, mv, dot, enc, coxeter_pair_value,
)
from scripts.fresh_rank3_h9_search import alphabet
from coxeter_polytopes.geometric.golden_fusion import golden_square_root


def positive(a, b):
    """Exact positivity of integer arrays a+b*phi; reject possible overflow."""
    A, B = 2*a+b, b
    assert np.max(np.abs(A), initial=0) < 10**9
    assert np.max(np.abs(B), initial=0) < 10**9
    return (((A >= 0) & (B >= 0) & ((A != 0) | (B != 0))) |
            ((A >= 0) & (B < 0) & (A*A > 5*B*B)) |
            ((A < 0) & (B >= 0) & (5*B*B > A*A)))


def run(t_bound, u_bound, z_bound):
    if not (0 < t_bound <= 300 and 0 < u_bound <= 600 and 0 < z_bound <= 600):
        raise ValueError('This int64 search is certified only for positive bounds at most 300,600,600.')
    start = time.monotonic()
    source = 'results/fresh_rank3_h9_cones_20261013.json'
    d = read(source); gram = matrix(d['gram']); walls = matrix(d['cases'][0]['all_walls'])
    v = tuple(map(G.from_json, d['cases'][0]['hyperideal_witness']['vector']))
    inside = tuple(map(G.from_json, d['cases'][0]['strict_interior']))
    einv = inverse(tuple(tuple(row[1:9]) for row in gram[1:9]))
    normals = tuple(mv(gram, r) for r in walls)
    norms = tuple(dot(r, n) for r, n in zip(walls, normals))
    phi = G(0, 1)
    # Pairing of (0, t*omega4+u*omega7, z) with the target v.
    assert v[5] == 6 and v[8] == G(5, 3)/2 and v[9] == G(-1, 1)/2
    cuts = []; censuses = []
    for q in (1, 2):
        ts = tuple(t for t in alphabet(t_bound, G(8*q)/15) if t.sign() > 0)
        us = tuple(u for u in alphabet(u_bound, G(4*q)/5) if u.sign() > 0)
        zs = tuple(z for z in alphabet(z_bound, G(q)) if int(z.a)%2 == q%2 and int(z.b)%2 == 0)
        UA = np.repeat(np.array([int(u.a) for u in us], dtype=np.int64), len(zs))
        UB = np.repeat(np.array([int(u.b) for u in us], dtype=np.int64), len(zs))
        ZA = np.tile(np.array([int(z.a) for z in zs], dtype=np.int64), len(us))
        ZB = np.tile(np.array([int(z.b) for z in zs], dtype=np.int64), len(us))
        U2a, U2b = UA*UA+UB*UB, 2*UA*UB+UB*UB
        Z2a, Z2b = ZA*ZA+ZB*ZB, 2*ZA*ZB+ZB*ZB
        counts = {'norm': q, 't': len(ts), 'u': len(us), 'z_after_parity': len(zs),
                  'triples': 0, 'conjugate_and_cut_gate': 0,
                  'integral_square_radicands': 0, 'square_norm_radicands': 0,
                  'integral_norm_solutions': 0, 'compatible_roots': 0}
        for t in ts:
            ta, tb = int(t.a), int(t.b)
            ea = 12*(ta*ta+tb*tb)+18*(ta*UA+tb*UB)+8*U2a+Z2a-q
            eb = 12*(2*ta*tb+tb*tb)+18*(ta*UB+tb*UA+tb*UB)+8*U2b+Z2b
            la = 12*ta+5*UA+3*UB-ZA+ZB
            lb = 12*tb+3*UA+8*UB+ZA
            ga = la*la+lb*lb-8*eb
            gb = 2*la*lb+lb*lb-8*(ea+eb)
            # E' <= q; E>q; and L² > 2phi*(E-q), exactly.
            gate = (~positive(ea+eb, -eb)) & positive(ea, eb) & positive(ga, gb)
            counts['triples'] += len(UA); counts['conjugate_and_cut_gate'] += int(gate.sum())
            gate &= (ea%2 == 0) & (eb%2 == 0)
            ids = np.flatnonzero(gate)
            counts['integral_square_radicands'] += len(ids)
            aa, bb = (eb[ids]-ea[ids])//2, ea[ids]//2
            assert np.max(np.abs(aa), initial=0) < 10**8 and np.max(np.abs(bb), initial=0) < 10**8
            nn = aa*aa+aa*bb-bb*bb
            for idx, a, b, n in zip(ids, aa, bb, nn):
                n = int(n)
                if n < 0 or isqrt(n)**2 != n:
                    continue
                counts['square_norm_radicands'] += 1
                x = golden_square_root(G(int(a), int(b)))
                if x is None or not x.is_integral:
                    continue
                if x.sign() < 0: x = -x
                u, z = G(int(UA[idx]), int(UB[idx])), G(int(ZA[idx]), int(ZB[idx]))
                r = (x, *(t*einv[k][4]+u*einv[k][7] for k in range(8)), z)
                assert all(y.is_integral for y in r)
                assert dot(r, mv(gram, r)) == q
                assert dot(r, mv(gram, v)).sign() > 0
                counts['integral_norm_solutions'] += 1
                products = tuple(dot(r, nrm) for nrm in normals)
                if any(p.sign() > 0 for p in products): continue
                if not all(coxeter_pair_value(4*p*p/(q*nrm)) for p, nrm in zip(products, norms)): continue
                assert all((2*y/q).is_integral for y in mv(gram, r))
                assert dot(r, mv(gram, inside)).sign() < 0
                counts['compatible_roots'] += 1
                cuts.append({'root': enc(r), 'norm': q, 'parameters': enc((x,t,u,z)),
                             'escape_pairing': dot(r, mv(gram, v)).to_json(),
                             'basepoint_height': (2*phi*x*x/q).to_json()})
        censuses.append(counts); print(json.dumps(counts), flush=True)
    out = {'source': source, 'source_case': 0, 'fundamental_weight_indices': [4,7],
           'bounds': {'t': t_bound, 'u': u_bound, 'scalar': z_bound}, 'censuses': censuses,
           'cuts': cuts, 'seconds': round(time.monotonic()-start, 3),
           'scope': 'Finite box of integral parameters t,u>0 and scalar>=0 in span(time,omega4,omega7,scalar). Exact conjugate bounds, scalar parity, norm equation and all 382 Coxeter gates; no global slice or host conclusion.'}
    (ROOT/'results/fresh_rank4_h9_free_search_20261013.json').write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'cuts': len(cuts), 'seconds': out['seconds']}), flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--t-bound', type=int, default=300)
    parser.add_argument('--u-bound', type=int, default=600)
    parser.add_argument('--z-bound', type=int, default=600)
    args=parser.parse_args(); run(args.t_bound, args.u_bound, args.z_bound)
