#!/usr/bin/env python3
"""Numerically close the three-cap patch after releasing the old34 constraint.

Delete old34 and old52, retain the other61 frozen walls, and move52 on the
norm-two circle with its contacts to S6,8,30 fixed. This continuous degree
of freedom is used to impose the formerly missing A--B Coxeter equation.
"""
import argparse
import json
import math
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import brentq

import fresh_nonarithmetic_s6_split_bridge_scout_20261017 as F

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_s6_release34_patch_scout_20261017'
S, M = F.S, F.METRIC


def theta(m):
    return 0. if m == 2 else 2*math.cos(math.pi/m)


def run():
    d = F.read(F.TAG)
    old = np.array(d['old_normals'])
    h = np.array([1.] + [0.]*9)
    retained = [i for i in range(63) if i not in [34, 52]]
    frame = old[S+[8, 30]]
    contacts = np.array([0.]*7+[-theta(4)])
    inverse = np.linalg.inv(frame @ M @ frame.T)
    p = contacts @ inverse @ frame
    u = old[52]-p
    rho = float(u @ M @ u)
    assert rho > 0 and max(abs(frame @ M @ u)) < 1e-6
    v = null_space(np.vstack([frame, u]) @ M)[:, 0]
    v *= math.sqrt(rho/(v @ M @ v))
    if v @ M @ old[34] < 0:
        v *= -1
    b = float(old[8] @ M @ old[30])
    target = -theta(27)

    def evaluate(t):
        replacement = p+((1-t*t)*u+2*t*v)/(1+t*t)
        x = float(old[44] @ M @ replacement)
        f = np.vstack([old[8], old[44], replacement])
        g3 = np.array([[2., 0, 0], [0, 2., x], [0, x, 2.]])
        j3 = np.linalg.inv(g3)
        g30 = np.array([b, 0, -theta(4)])
        cc = np.array([0., -theta(15), -theta(4)])
        r = 2-cc @ j3 @ cc
        w2 = 2-g30 @ j3 @ g30
        if min(r, w2) <= 0:
            return None
        c = cc @ j3 @ f-math.sqrt(r/w2)*(old[30]-g30 @ j3 @ f)
        c30 = float(cc @ j3 @ g30-math.sqrt(r*w2))
        g4 = np.array([[2., 0, -theta(4), c30],
                       [0, 2., x, -theta(15)],
                       [-theta(4), x, 2., -theta(4)],
                       [c30, -theta(15), -theta(4), 2.]])
        j4 = np.linalg.inv(g4)
        if min(j4[1, 1], j4[2, 2]) <= 0:
            return None
        basis = np.vstack([old[30], old[44], replacement, c])
        a = -math.sqrt(2/j4[2, 2])*(j4[2] @ basis)
        bb = -math.sqrt(2/j4[1, 1])*(j4[1] @ basis)
        ab = 2*j4[1, 2]/math.sqrt(j4[1, 1]*j4[2, 2])
        return dict(parameter=float(t), AB_product=float(ab), residual=float(ab-target),
                    distance44_52=x, new52=replacement, A=a, B=bb, C=c)

    initial = evaluate(0.)
    assert initial is not None and abs(initial['residual']) < 3e-6
    grid = sorted(set([0.]+[s*10.**e for s in [-1, 1] for e in range(-10, -1)]))
    samples = [evaluate(t) for t in grid]
    brackets = []
    for l, r in zip(samples, samples[1:]):
        if l is not None and r is not None and l['residual']*r['residual'] < 0:
            brackets.append([l['parameter'], r['parameter']])
    hits = []
    for lo, hi in brackets:
        t = brentq(lambda t:evaluate(t)['residual'], lo, hi, xtol=1e-15)
        rec = evaluate(t)
        N = np.vstack([old[retained], rec['new52'], rec['A'], rec['B'], rec['C']])
        labels = retained+[52, 63, 64, 65]
        prescribed = {
            52:{**{i:2 for i in S+[8]}, 30:4},
            63:{**{i:2 for i in S+[30, 44, 65]}, 64:27},
            64:{**{i:2 for i in S+[30, 52, 65]}, 63:27},
            65:{**{i:2 for i in S+[8, 63, 64]}, 44:15, 52:4},
        }
        products = N @ M @ N.T
        gates, failures = Counter(), []
        for i, j in combinations(range(len(labels)), 2):
            li, lj = labels[i], labels[j]
            if li in retained and lj in retained:
                continue
            m = prescribed.get(li, {}).get(lj, prescribed.get(lj, {}).get(li))
            g = float(products[i, j])
            if m is not None:
                error = abs(g+theta(m))
                if error > 2e-6:
                    failures.append(dict(pair=[li, lj], product=g, expected_label=m, error=error))
                gates['finite_prescribed'] += 1
            elif g < -2-1e-7:
                gates['strict_dotted'] += 1
            else:
                failures.append(dict(pair=[li, lj], product=g, reason='not_prescribed_or_dotted'))
        # The eight new finite corners bound the closed left patch; every other
        # wall is checked, so this test does not rely on isolated diagrams.
        corner_sides = [[8, 44, 65], [8, 52, 65], [30, 44, 63],
                        [30, 52, 64], [44, 63, 65], [52, 64, 65],
                        [30, 63, 64], [63, 64, 65]]
        corners = []
        for sides in corner_sides:
            ids = S+sides
            f = N[[labels.index(i) for i in ids]]
            g = f @ M @ f.T
            eig = float(np.linalg.eigvalsh(g)[0])
            vertex = h-(h @ M @ f.T) @ np.linalg.inv(g) @ f
            values = vertex @ M @ N.T
            corners.append(dict(sides=sides, minimum_gram_eigenvalue=eig,
                                norm=float(vertex @ M @ vertex), maximum_wall_product=float(max(values)),
                                vertex=vertex.tolist()))
        rec = {k:val.tolist() if isinstance(val, np.ndarray) else val for k, val in rec.items()}
        rec.update(bracket=[lo, hi], gates=dict(gates), failures=failures,
                   removed34_new52_product=float(N[labels.index(52)] @ M @ old[34]),
                   common_old_center_maximum=float(max(h @ M @ N.T)),
                   maximum_norm_error=float(max(abs(np.diag(products)-2))), corners=corners,
                   wall_labels=labels, normals=N.tolist())
        hits.append(rec)
    paths = [f'scripts/{TAG}.py', f'scripts/{F.TAG}.py', f'results/{F.TAG}.json']
    return dict(status='COMPLETED_NUMERICAL_RELEASE34_THREE_CAP_PATCH_SCOUT',
                sources={q:sha256((ROOT/q).read_bytes()).hexdigest() for q in paths},
                removed_old_walls=[34, 52], retained_old_wall_count=len(retained),
                new_wall_labels=[52, 63, 64, 65], replacement_fixed_contacts=S+[8, 30],
                target_AB_label=27, initial_residual=initial['residual'],
                circle=dict(p=p.tolist(), u=u.tolist(), v=v.tolist(), squared_radius=rho),
                samples=[None if r is None else {k:r[k] for k in ['parameter', 'residual', 'distance44_52']} for r in samples],
                hits=hits, metric=[-1]+[1]*9,
                scope='Numerical root discovery after explicit deletion of old34 and replacement of old52. Prescribed angle equality, whole local-face closure and finite volume require separate exact certificates. No old-corner persistence is assumed across deleted walls.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = json.loads(json.dumps(run()))
    out = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], 'hits', len(result['hits']), flush=True)
    for r in result['hits']:
        print('parameter', r['parameter'], 'angle_residual', r['residual'], 'gates', r['gates'],
              'failures', r['failures'], 'corner_min', min(c['minimum_gram_eigenvalue'] for c in r['corners']), flush=True)
