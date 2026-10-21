#!/usr/bin/env python3
"""New bounded A1+A1 -> A2 section scout; numerical proposals only."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from pathlib import Path
import numpy as np
import fresh_h20_h4_edge_sections_20261013_scout as S

ROOT = S.ROOT
TAG = 'fresh_h20_compact_section_coverage_20261014'
PAIRS = [('1/2', '6/7'), ('2/3', '2/3'), ('1', '2/5')]
PARAMETERS = ['1/5', '1/2', '4/5']


def run():
    Q, W, G, h = S.load()
    C = W @ Q
    walk = json.loads((ROOT / S.WALK).read_text())['corner_orbits']
    seen = set()
    catalogue = []
    for c in walk:
        for ei, edge in enumerate(c['edge_orbits']):
            E = edge['edge_walls']
            if tuple(E) in seen:
                continue
            seen.add(tuple(E))
            cc = S.comps(E, G)
            single = [a[0] for a in cc if len(a) == 1 and G[a[0], a[0]] == 2]
            if len(single) < 2 or max(map(len, cc)) == 1:
                continue
            catalogue.append(dict(corner=c['id'], edge_index=ei, edge_walls=E,
                                  contact_walls=single[:2],
                                  component_sizes=sorted(map(len, cc)),
                                  endpoint_types=sorted([c['type'], walk[edge['endpoint_orbit']]['type']])))
    selected = []
    for types in [('finite', 'finite'), ('finite', 'ideal'), ('ideal', 'ideal')]:
        signatures = set()
        for row in catalogue:
            sig = tuple(row['component_sizes'])
            if tuple(row['endpoint_types']) != types or sig in signatures:
                continue
            selected.append(row)
            signatures.add(sig)
            if len(signatures) == 3:
                break
    assert len(selected) == 9
    cases = []
    allowed = np.array([0., .25, .5, .75])
    for si, seed in enumerate(selected):
        E = seed['edge_walls']
        edge = walk[seed['corner']]['edge_orbits'][seed['edge_index']]
        I = np.linalg.inv(G[np.ix_(E, E)].astype(float))
        left = np.array(walk[seed['corner']]['primitive_ray'], float)
        right = np.array(edge['endpoint_primitive_ray'], float)
        left /= -(left @ Q @ h)
        right /= -(right @ Q @ h)
        for ar, br in PAIRS:
            a = np.zeros(20)
            for j, v in zip(seed['contact_walls'], [ar, br]):
                a[E.index(j)] = np.sqrt(float(F(v)))
            coeff = I @ a
            gamma = 1 - a @ coeff
            w = coeff @ W[E]
            H = G[np.ix_(E, E)] - np.outer(a, a)
            assert gamma > 0 and np.linalg.eigvalsh(H)[0] > 1e-8
            k, l = map(E.index, seed['contact_walls'])
            assert abs(H[k, l]**2/(H[k, k]*H[l, l])-.25) < 1e-10
            for ts in PARAMETERS:
                t = float(F(ts))
                p = (1-t)*left+t*right
                qp = p @ Q @ p
                z0 = right - (right @ Q @ p/qp)*p
                qz = z0 @ Q @ z0
                pp = C @ p
                assert qp < 0 and qz > 0
                assert max(abs(pp[E])) < 1e-7
                assert max(pp[i] for i in range(210) if i not in E) < -1e-8
                for sign in [-1, 1]:
                    z = sign*np.sqrt(gamma/qz)*z0
                    e = w+z
                    ep = C @ e
                    dual = I @ W[E] - np.outer(coeff/gamma, z)
                    assert abs(e @ Q @ e-1) < 1e-7
                    bad = []
                    for k, release in enumerate(E):
                        dp = C @ dual[k]
                        choices = [(pp[j]/dp[j], j) for j in range(210)
                                   if j not in E and dp[j] < -1e-9]
                        if not choices:
                            continue
                        step, new = min(choices)
                        x = p-step*dual[k]
                        xp = C @ x
                        active = np.flatnonzero(abs(xp) < 1e-7).tolist()
                        if x @ Q @ x >= -1e-8 or len(active) != 20 or max(xp) > 1e-7:
                            continue
                        A = G[np.ix_(active, active)]-np.outer(ep[active], ep[active])
                        if np.linalg.eigvalsh(A)[0] < 1e-7:
                            continue
                        defects = []
                        for u, v in combinations(range(20), 2):
                            pair = A[u, v]
                            cs = pair*pair/(A[u, u]*A[v, v])
                            if pair > 1e-7:
                                defects.append((pair, 'obtuse', active[u], active[v], cs))
                            elif pair < -1e-7 and min(abs(allowed-cs)) > 1e-5:
                                if 0 < cs < .25:
                                    defects.append((min(cs, .25-cs), 'between_0_and_1/4', active[u], active[v], cs))
                                elif .25 < cs < .5:
                                    defects.append((min(cs-.25, .5-cs), 'between_1/4_and_1/2', active[u], active[v], cs))
                        if defects:
                            score, kind, i, j, cs = max(defects)
                            bad.append(dict(release=release, new_wall=new, active=active,
                                            bad_pair=[i, j], kind=kind, score=float(score),
                                            cosine_squared=float(cs)))
                    rec = dict(seed=si, contact_squares=[ar, br], parameter=ts,
                               radical_sign=sign, numerical_normal=e.tolist())
                    if bad:
                        rec.update(status='NUMERICAL_ACTUAL_VERTEX_DEFECT', witness=max(bad, key=lambda b:b['score']))
                    else:
                        rec.update(status='UNRESOLVED')
                    cases.append(rec)
        print('seed', si, 'cases', len(cases), 'unresolved', sum(r['status']=='UNRESOLVED' for r in cases), flush=True)
    paths = [S.RAW, S.WALK, 'scripts/fresh_h20_h4_edge_sections_20261013_scout.py', f'scripts/{TAG}_scout.py']
    out = dict(status='BOUNDED_NEW_A1_A1_TO_A2_SECTION_SCOUT',
               sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
               unique_literal_edges=len(seen), eligible_literal_edges=len(catalogue),
               selected_edges=selected, contact_square_pairs=PAIRS, parameters=PARAMETERS,
               cases=cases, counts=dict(Counter(r['status'] for r in cases)),
               scope='Nine recorded old edge supports; three same-sign contact ratios, three positions, both residual signs. '
                     'All-right edge supports omitted. Numerical proposals only; no complete continuous family or all-profile exclusion.')
    (ROOT/f'results/{TAG}_scout.json').write_text(json.dumps(out, indent=2)+'\n')
    print(out['counts'])


if __name__ == '__main__':
    run()
