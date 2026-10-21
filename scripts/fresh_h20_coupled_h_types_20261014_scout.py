#!/usr/bin/env python3
"""Find coupled H3/H4 local profiles outside the one-contact obstructions.

Numerical profile discovery and actual-neighbor checks only.
"""
from pathlib import Path
from itertools import combinations, permutations, product
from collections import Counter
import hashlib, json, math
import numpy as np
from scipy.optimize import least_squares
import fresh_h20_h4_edge_sections_20261013_scout as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_coupled_h_types_20261014'


def run():
    Q, W, G, h = S.load()
    C = W@Q
    corners = json.loads((ROOT/S.WALK).read_text())['corner_orbits']
    catalog = {}
    joins = {}
    seen = set()
    for corner in corners:
        for ei, edge in enumerate(corner['edge_orbits']):
            E = edge['edge_walls']
            if tuple(E) in seen:
                continue
            seen.add(tuple(E))
            cc = S.comps(E, G)
            base = {'corner': corner['id'], 'edge_index': ei, 'edge_walls': E}
            for comp in cc:
                if len(comp) in [3, 4]:
                    key, order = min((tuple(G[np.ix_(p, p)].ravel()), p) for p in permutations(comp))
                    catalog.setdefault(key, dict(base, component=list(order)))
                if len(comp) == 2:
                    for iso in [c[0] for c in cc if len(c) == 1]:
                        for a, b in [comp, comp[::-1]]:
                            key = (int(G[iso, iso]), int(G[a, a]), int(G[b, b]), int(G[a, b]))
                            joins.setdefault(key, dict(base, component=[iso, a, b]))
    profiles = []
    c5 = math.cos(math.pi/5)
    for key, seed in joins.items():
        ids = seed['component']
        B = G[np.ix_(ids, ids)].astype(float)
        cold = B[1, 2]**2/(B[1, 1]*B[2, 2])
        y = 1-cold/c5**2
        x = (1-y)/(1+3*y)
        a = np.array([math.sqrt(B[0, 0]*x), math.sqrt(B[1, 1]*y), 0])
        H = B-np.outer(a, a)
        assert np.linalg.eigvalsh(H)[0] > 1e-6
        profiles.append(dict(seed, kind='A1_plus_rank2_to_H3', pairings=a.tolist(),
                             squared_normalized_contacts=[x, y], old_Gram=B.tolist()))
    for key, seed in catalog.items():
        ids = seed['component']
        B = G[np.ix_(ids, ids)].astype(float)
        n = len(ids)
        target = np.eye(n)
        for i in range(n-1):
            target[i, i+1] = target[i+1, i] = -c5 if i == 0 else -.5
        solutions = []
        for center in range(n):
            adj = [i for i in range(n) if i != center and B[center, i] < 0]
            for leaf1, leaf2 in combinations(adj, 2):
                for fixed, free in [(leaf1, leaf2), (leaf2, leaf1)]:
                    low = abs(B[center, fixed])/math.sqrt(B[center, center])
                    high = math.sqrt(B[fixed, fixed])
                    if low >= high:
                        continue
                    for order in permutations(range(n)):
                        T = target[np.ix_(order, order)]
                        if T[center, fixed] != 0:
                            continue
                        pairs = list(combinations(range(n), 2))
                        def make(v):
                            a = np.zeros(n)
                            a[fixed] = low+(high-low)*v[0]
                            a[center] = B[center, fixed]/a[fixed]
                            a[free] = math.sqrt(B[free, free])*v[1]
                            return a
                        def residual(v):
                            a = make(v)
                            H = B-np.outer(a, a)
                            d = np.sqrt(np.diag(H))
                            return np.array([H[i, j]/(d[i]*d[j])-T[i, j] for i, j in pairs])
                        for start in product([.2, .5, .8], repeat=2):
                            opt = least_squares(residual, start, bounds=(1e-7, 1-1e-7),
                                                max_nfev=150, ftol=1e-11, xtol=1e-11, gtol=1e-11)
                            if max(abs(residual(opt.x))) > 1e-7:
                                continue
                            a = make(opt.x)
                            if np.linalg.eigvalsh(B-np.outer(a, a))[0] < 1e-7:
                                continue
                            if any(max(abs(a-old)) < 1e-6 for old in solutions):
                                continue
                            solutions.append(a)
                            profiles.append(dict(seed, kind=f'three_contact_H{n}',
                                                 pairings=a.tolist(), target_order=list(order),
                                                 old_Gram=B.tolist(), canceled_pair=[ids[center], ids[fixed]]))
        print('old block', n, 'solutions', len(solutions), flush=True)
    allowed = np.array([math.cos(math.pi/m)**2 for m in range(2, 101)])
    cases = []
    for si, seed in enumerate(profiles):
        E = seed['edge_walls']
        I = np.linalg.inv(G[np.ix_(E, E)])
        a = np.zeros(20)
        for j, val in zip(seed['component'], seed['pairings']):
            a[E.index(j)] = val
        v = I@a
        spatial = a@v
        gamma = 1-spatial
        assert gamma > 1e-8
        corner = corners[seed['corner']]
        edge = corner['edge_orbits'][seed['edge_index']]
        left = np.array(corner['primitive_ray'], float)
        right = np.array(edge['endpoint_primitive_ray'], float)
        left /= -(left@Q@h)
        right /= -(right@Q@h)
        w = v@W[E]
        for t, sign in product([1/3, .5, 2/3], [-1, 1]):
            P = (1-t)*left+t*right
            pp = C@P
            qp = P@Q@P
            z0 = right-(right@Q@P/qp)*P
            z = sign*math.sqrt(gamma/(z0@Q@z0))*z0
            e = w+z
            ep = C@e
            row = {'profile': si, 'position': t, 'branch': sign, 'normal': e.tolist()}
            for ri, release in enumerate(E):
                D = I[ri]@W[E]-v[ri]*z/gamma
                dp = C@D
                ids = [j for j in range(210) if j not in E and dp[j] < -1e-8]
                if not ids:
                    continue
                step, new = min((pp[j]/dp[j], j) for j in ids)
                xp = pp-step*dp
                qx = qp+step*step*(D@Q@D)
                active = np.flatnonzero(abs(xp) < 1e-7).tolist()
                if qx >= -1e-8 or len(active) != 20 or max(xp) > 1e-7:
                    continue
                H = G[np.ix_(active, active)]-np.outer(ep[active], ep[active])
                if np.linalg.eigvalsh(H)[0] < 1e-7:
                    continue
                defects = []
                for i, j in combinations(range(20), 2):
                    p = H[i, j]
                    cs = p*p/(H[i, i]*H[j, j])
                    if p > 1e-7:
                        defects.append((float(p), 'obtuse', active[i], active[j], float(cs)))
                    elif p < -1e-7 and min(abs(allowed-cs)) > 1e-5:
                        defects.append((float(min(abs(allowed-cs))), 'forbidden_numerical', active[i], active[j], float(cs)))
                if defects:
                    defect = max(defects)
                    row.update(status='NUMERICAL_ACTUAL_RIDGE_DEFECT',
                               witness={'release': release, 'new_wall': new, 'active': active,
                                        'pair': list(defect[2:4]), 'kind': defect[1], 'cosine_squared': defect[4]})
                    break
            row.setdefault('status', 'UNRESOLVED')
            cases.append(row)
        print('profile', si, seed['kind'], 'unresolved', sum(r['status']=='UNRESOLVED' for r in cases), flush=True)
    files = [S.RAW, S.WALK, 'scripts/'+Path(__file__).name, 'scripts/fresh_h20_h4_edge_sections_20261013_scout.py']
    out = {'status': 'NUMERICAL_COUPLED_NONCRYSTALLOGRAPHIC_PROFILE_SCOUT',
           'sources': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
           'old_connected_block_types': len(catalog), 'old_join_types': len(joins),
           'profiles': profiles, 'cases': cases, 'counts': dict(Counter(r['status'] for r in cases)),
           'scope': 'Bounded numerical local-profile discovery and three positions per selected source edge. Not a complete classification or exact geometric rejection.'}
    (ROOT/'results'/f'{TAG}_scout.json').write_text(json.dumps(out, indent=2)+'\n')
    print(len(profiles), out['counts'], flush=True)


if __name__ == '__main__':
    run()
