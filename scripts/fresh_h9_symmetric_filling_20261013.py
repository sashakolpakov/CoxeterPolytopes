#!/usr/bin/env python3
"""Numerical scout of simultaneous symmetric angle changes of an H9 chamber.

Two center walls and eight interchangeable pairs of leaf walls. All normals
have square 2. Unlike the previous sparse surgery, whole orbits of old right
angles change at once. The S8-standard Gram blocks force w=2+a+b or b-2-a;
the remaining rank condition determines the center-center weight.
"""
import json
import time
from collections import Counter
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/enumerations/fresh_h9_symmetric_filling_20261013.json'


def weights():
    result = [(f'm{m}', float(2*np.cos(np.pi/m))) for m in (2,3,4,5,6,8,10,12)]
    result[0] = ('m2', 0.)
    result += [('sqrt5', np.sqrt(5)), ('1+sqrt2', 1+np.sqrt(2)),
               ('phi2', (3+np.sqrt(5))/2), ('3', 3.),
               ('2sqrt3', 2*np.sqrt(3)), ('2phi', 1+np.sqrt(5)), ('4', 4.)]
    return result


def active_orbits():
    rows = []
    # For a simple ray, nine independent active facets occur. Every nine-set
    # is equivalent under S8 to one of these patterns, even at nonsimple rays.
    for ca, cb in product((0,1), repeat=2):
        for neither in range(9):
            for only_u in range(9-neither):
                for only_v in range(9-neither-only_u):
                    both = 8-neither-only_u-only_v
                    if ca+cb+only_u+only_v+2*both != 9:
                        continue
                    chosen = ([0] if ca else []) + ([1] if cb else [])
                    for i in range(neither, neither+only_u):
                        chosen.append(2+2*i)
                    for i in range(neither+only_u, neither+only_u+only_v):
                        chosen.append(3+2*i)
                    for i in range(8-both, 8):
                        chosen.extend((2+2*i, 3+2*i))
                    rows.append({'pattern': [ca,cb,neither,only_u,only_v,both],
                                 'active': chosen})
    return rows


def matrix(a,b,p,q,w,c):
    G = 2*np.eye(18)
    G[0,1] = G[1,0] = -c
    for i in range(8):
        u,v = 2+2*i,3+2*i
        G[0,u] = G[u,0] = G[1,v] = G[v,1] = -p
        G[0,v] = G[v,0] = G[1,u] = G[u,1] = -q
        G[u,v] = G[v,u] = -w
        for j in range(i):
            s,t = 2+2*j,3+2*j
            G[u,s] = G[s,u] = G[v,t] = G[t,v] = -a
            G[u,t] = G[t,u] = G[v,s] = G[s,v] = -b
    return G


def gate(G, orbits):
    eig, U = np.linalg.eigh(G)
    if sum(eig < -1e-7) != 1 or sum(eig > 1e-7) != 9:
        return {'status': 'wrong_signature'}
    indices = [0, *np.flatnonzero(eig > 1e-7)]
    R = U[:,indices]*np.sqrt(abs(eig[indices]))
    if sum(R[:,0]) < 0:
        R[:,0] *= -1
    if min(R[:,0]) <= 1e-8:
        return {'status': 'no_strict_symmetric_center'}
    N = R.copy(); N[:,0] *= -1
    rays = []
    for row in orbits:
        M = N[row['active']]
        _, s, V = np.linalg.svd(M)
        if s[-1] < 1e-8:
            continue
        v = V[-1]
        if abs(v[0]) < 1e-9:
            continue
        v /= v[0]
        if max(N@v) > 1e-6:
            continue
        norm = float(-1 + v[1:]@v[1:])
        witness = {'active_pattern': row['pattern'], 'active': row['active'],
                   'norm_time_normalized': norm,
                   'maximum_pairing': float(max(N@v)), 'vector': v.tolist()}
        rays.append(witness)
        if norm > 1e-6:
            return {'status': 'numerical_infinite_volume', 'escape': witness}
    if not rays:
        return {'status': 'no_rays_resolved'}
    maximum = max(row['norm_time_normalized'] for row in rays)
    return {'status': 'numerical_compact_candidate' if maximum < -1e-6 else
            'numerical_finite_volume_noncompact_candidate',
            'ray_orbits': len(rays), 'maximum_norm': maximum, 'rays': rays}


def run():
    alphabet = weights(); orbits = active_orbits(); counts = Counter()
    survivors = []; escapes = []; start = time.monotonic()
    finite = [(m,2*np.cos(np.pi/m)) for m in range(2,61)]
    for branch, rank_block in product((1,-1), ('minus','plus')):
        for (an,a),(bn,b),(pn,p),(qn,q) in product(alphabet,repeat=4):
            # Exchanging the two leaf species relative to the centers swaps p,q.
            if p < q:
                continue
            counts['parameter_tuples'] += 1
            w = 2+a+b if branch == 1 else b-2-a
            if w < -1e-8 or abs(w-2) < 1e-8:
                counts['wrong_matched_weight'] += 1; continue
            if w < 2 and not any(abs(w-x)<1e-8 for _,x in finite):
                counts['noncoxeter_matched_weight'] += 1; continue
            lm = 4-6*a+8*b if branch == 1 else -8*a+8*b
            if lm <= 1e-8:
                counts['wrong_minus_block'] += 1; continue
            lp = -8*(a+b) if branch == 1 else 4-6*a-8*b
            if rank_block == 'plus' and lp >= -1e-8:
                counts['wrong_plus_block'] += 1; continue
            c = 8*(p-q)**2/lm-2 if rank_block == 'minus' else 2-8*(p+q)**2/lp
            if c < -1e-8 or abs(c-2) < 1e-8:
                counts['wrong_center_weight'] += 1; continue
            if c < 2 and not any(abs(c-x)<1e-8 for _,x in finite):
                counts['noncoxeter_center_weight'] += 1; continue
            counts['gram_proposals'] += 1
            G = matrix(a,b,p,q,w,c)
            result = gate(G,orbits); counts[result['status']] += 1
            record = {'branch': branch, 'rank_block':rank_block,
                      'parameters': dict(zip('abpq',(an,bn,pn,qn))),
                      'matched_weight':float(w), 'center_weight':float(c), 'geometry':result}
            if 'candidate' in result['status']:
                record['gram'] = G.tolist(); survivors.append(record)
                print(json.dumps(record),flush=True)
            elif result['status'] == 'numerical_infinite_volume':
                escapes.append(record)
    data = {'schema':'H9-symmetric-S8-filling-scout-v1','alphabet':alphabet,
            'active_set_orbits':orbits,'counts':dict(counts),'survivors':survivors,
            'escapes':escapes,'seconds':round(time.monotonic()-start,3),
            'scope':'Floating point scout of explicitly listed symmetric parameters; '
                    'no exact absence, no all-real-parameter, and no worldwide novelty claim.'}
    OUT.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'counts':dict(counts),'seconds':data['seconds'],
                      'active_orbits':len(orbits)}),flush=True)


if __name__ == '__main__':
    run()
