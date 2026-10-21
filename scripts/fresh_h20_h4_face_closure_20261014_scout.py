#!/usr/bin/env python3
"""Numerical finite-candidate reduction of continuous A4-to-H4 sections.

The necessary H4-face closure equation is mathematical, but all predicates
in this scout are floating point. No exact exclusion is asserted here.
"""
from pathlib import Path
import collections
import hashlib
import json
import math
import numpy as np
from fresh_h20_h4_edge_sections_20261013_scout import load, comps, A4

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_h4_face_closure_20261014'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK = 'results/fresh_h21_source_corner_walk_20261013.json'


def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()


def run():
    Q, W, G, h = load()
    C = W @ Q
    corners = json.loads((ROOT / WALK).read_text())['corner_orbits']
    qe = (1 + math.sqrt(5)) / 4
    gamma = qe - .8
    seen = set()
    families = []
    candidates = []
    counts = collections.Counter()
    allowed = np.array([math.cos(math.pi/m)**2 for m in range(2, 101)])
    for corner in corners:
        for ei, edge in enumerate(corner['edge_orbits']):
            E = edge['edge_walls']
            if tuple(E) in seen:
                continue
            seen.add(tuple(E))
            components = [cc for cc in comps(E, G) if A4(cc, G)]
            if not components:
                continue
            I = np.linalg.inv(G[np.ix_(E, E)].astype(float))
            dual = I @ G[E]
            v = np.array(corner['primitive_ray'], float)
            b = np.array(edge['endpoint_primitive_ray'], float)
            v /= -(v @ Q @ h)
            b /= -(b @ Q @ h)
            pv, pb = C @ v, C @ b
            plane = np.array([[v @ Q @ v, v @ Q @ b],
                              [b @ Q @ v, b @ Q @ b]])
            assert np.linalg.det(plane) < 0
            assert max(abs(pv[E])) < 1e-8 and max(abs(pb[E])) < 1e-8
            outside = [l for l in range(210) if l not in E]
            for cc in components:
                for anchor in cc:
                    if sum(G[anchor, j] != 0 for j in cc) != 2:
                        continue
                    ai = E.index(anchor)
                    eligible = [l for l in outside
                                if all(G[l, i] == 0 for i in cc if i != anchor)]
                    family = {'corner': corner['id'], 'edge_index': ei,
                              'edge_walls': E, 'component': cc, 'anchor': anchor,
                              'eligible_closure_walls': eligible}
                    fid = len(families)
                    families.append(family)
                    wc = dual[ai]
                    unique = []
                    for l in eligible:
                        p = np.array([pv[l], pb[l]])
                        assert np.dot(p, p) > 1e-12
                        k = gamma * G[l, anchor]
                        z0 = k * p / np.dot(p, p)
                        y = np.array([p[1], -p[0]])
                        aa = y @ plane @ y
                        bb = 2 * z0 @ plane @ y
                        dd = z0 @ plane @ z0 - gamma
                        if abs(aa) < 1e-13:
                            if abs(bb) < 1e-13:
                                counts['degenerate_no_solution'] += 1
                                continue
                            roots = [-dd / bb]
                        else:
                            disc = bb*bb - 4*aa*dd
                            if disc < -1e-13:
                                counts['negative_discriminant'] += 1
                                continue
                            if abs(disc) <= 1e-13:
                                counts['near_discriminant_boundary'] += 1
                            roots = [(-bb + s*math.sqrt(max(0, disc))) / (2*aa)
                                     for s in [-1, 1]]
                        for branch, root in enumerate(roots):
                            zc = z0 + root*y
                            zv, zb = plane @ zc
                            if zv*zb >= -1e-12:
                                counts['not_transverse_interior'] += 1
                                continue
                            if any(np.max(abs(zc - old)) < 1e-7 for old in unique):
                                counts['duplicate_within_family'] += 1
                                continue
                            unique.append(zc)
                            t = -zv / (zb-zv)
                            pp = (1-t)*pv + t*pb
                            qp = np.array([1-t, t]) @ plane @ np.array([1-t, t])
                            zp = zc[0]*pv + zc[1]*pb
                            ep = wc + zp
                            assert qp < 0 and 0 < t < 1
                            assert max(pp[outside]) < -1e-8
                            assert abs(ep[l] - qe*G[l, anchor]) < 1e-7
                            row = {'family': fid, 'closure_wall': l,
                                   'quadratic_branch': branch,
                                   'edge_plane_coefficients': zc.tolist(),
                                   'position': float(t)}
                            checked = 0
                            for release in [i for i in E if i not in cc] + cc:
                                ri = E.index(release)
                                dp = dual[ri] - I[ri, ai]*zp/gamma
                                ids = [j for j in outside if dp[j] < -1e-9]
                                if not ids:
                                    continue
                                step, j = min((pp[j]/dp[j], j) for j in ids)
                                xp = pp - step*dp
                                qd = I[ri, ri] + I[ri, ai]**2/gamma
                                qx = qp + step*step*qd
                                act = np.flatnonzero(abs(xp) < 1e-7).tolist()
                                if qx >= -1e-8 or len(act) != 20 or max(xp) > 1e-7:
                                    continue
                                H = G[np.ix_(act, act)] - np.outer(ep[act], ep[act])/qe
                                if np.linalg.eigvalsh(H)[0] < 1e-7:
                                    continue
                                checked += 1
                                defects = []
                                for ia in range(20):
                                    for ib in range(ia):
                                        pval = H[ia, ib]
                                        cs = pval*pval/(H[ia, ia]*H[ib, ib])
                                        if pval > 1e-7:
                                            defects.append((float(pval), 'obtuse', act[ia], act[ib], float(cs)))
                                        elif pval < -1e-7 and np.min(abs(allowed-cs)) > 1e-5:
                                            kind = ('between_0_and_1/4' if cs < .25 else
                                                    'between_1/4_and_1/2' if cs < .5 else 'other_numeric_bad')
                                            defects.append((float(np.min(abs(allowed-cs))), kind, act[ia], act[ib], float(cs)))
                                if defects:
                                    defect = max(defects)
                                    row['status'] = 'NUMERICAL_ACTUAL_RIDGE_DEFECT'
                                    row['witness'] = {'release': release, 'new_wall': j,
                                                      'active': act, 'vertex_norm': float(qx),
                                                      'pair': list(defect[2:4]), 'kind': defect[1],
                                                      'cosine_squared': defect[4]}
                                    break
                            row.setdefault('status', 'UNRESOLVED')
                            row['checked_finite_neighbors'] = checked
                            counts[row['status']] += 1
                            candidates.append(row)
                    if fid % 100 == 0:
                        print('families', len(families), 'candidates', len(candidates),
                              'unresolved', counts['UNRESOLVED'], flush=True)
    out = {'status': 'NUMERICAL_H4_FACE_CLOSURE_CENSUS',
           'sources': {p: sha(p) for p in [RAW, WALK, 'scripts/'+Path(__file__).name,
                       'scripts/fresh_h20_h4_edge_sections_20261013_scout.py']},
           'unique_literal_edge_supports': len(seen),
           'families': families, 'candidates': candidates, 'counts': dict(counts),
           'scope': 'All saved literal source edges with a norm2 A4 component, both endpoints of each A4. Closure equation is necessary for a Coxeter section. Floating-point enumeration and rejection only; unresolved or near-boundary cases require exact analysis.'}
    (ROOT/'results'/f'{TAG}_scout.json').write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], len(families), len(candidates), dict(counts), flush=True)


if __name__ == '__main__':
    run()
