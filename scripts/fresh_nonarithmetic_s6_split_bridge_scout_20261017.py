#!/usr/bin/env python3
"""Numerical discovery of two separately terminating S6 endpoint walls.

The new topology asks for endpoints on30/44 and30/52, then a mutual finite
corner on30. No wall is required to terminate both old edges. This first
scout restricts all unprescribed old contacts to be strictly dotted.
"""
import argparse
import json
import math
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import null_space

import fresh_nonarithmetic_escape62_20261017_cap as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_s6_split_bridge_scout_20261017'
ENDPOINT = 'fresh_nonarithmetic_h3_endpoint_domains_20261017'
S = [0, 1, 2, 3, 5, 6]
EDGES = [(30, 44), (30, 52)]
TAIL_LABELS = [6, 8, 10, 16, 32, 64, 256]
MUTUAL_LABELS = list(range(2, 17)) + [20, 24, 32, 64]
METRIC = np.diag([-1.] + [1.] * 9)


def read(tag):
    result = json.loads((ROOT / f'results/{tag}.json').read_text())
    for path, digest in result['sources'].items():
        assert sha256((ROOT / path).read_bytes()).hexdigest() == digest, path
    return result


def high(x):
    q = (x.l+x.h)/2 if hasattr(x, 'l') else (Fraction(x[0])+Fraction(x[1]))/2
    return mp.mpf(q.numerator)/q.denominator


def roots(a, b, c):
    scale = max(abs(a), abs(b), abs(c))
    if scale == 0:
        return []
    if abs(a) < 1e-13*scale:
        return [] if abs(b) < 1e-13*scale else [-c/b]
    disc = b*b-4*a*c
    if disc < 0:
        return []
    q = -0.5*(b+math.copysign(math.sqrt(disc), b))
    if q == 0:
        return [0.]
    return [q/a, c/q]


def intervals(polynomials):
    points = sorted(x for p in polynomials for x in roots(*p)
                    if math.isfinite(x) and x > 0)
    distinct = []
    for x in points:
        if not distinct or x-distinct[-1] > 1e-9*max(1., x):
            distinct.append(x)
    result = []
    for lo, hi in zip([0.] + distinct, distinct + [math.inf]):
        t = max(lo+1., 1.2*lo) if math.isinf(hi) else (lo+hi)/2
        vals = polynomials @ np.array([t*t, t, 1.])
        scales = abs(polynomials) @ np.array([t*t, t, 1.])
        if np.all(vals <= 1e-10*np.maximum(scales, 1.)):
            if result and abs(result[-1][1]-lo) < 1e-9*max(1., lo):
                result[-1][1] = hi
            else:
                result.append([lo, hi])
    return result


def samples(lo, hi):
    if math.isinf(hi):
        values = [max(lo, 0.01)*f + 1e-7 for f in
                  [1.00001, 1.001, 1.01, 1.05, 1.1, 1.25, 1.5, 2, 3, 5,
                   8, 13, 25, 50, 100, 500, 2000, 10000]]
    else:
        values = [lo+(hi-lo)*f for f in
                  [0.00001, .001, .01, .05, .1, .2, .35, .5, .65, .8, .9, .95, .99, .999, .99999]]
    return values


def evaluate(f, t):
    return f['p']+t*f['a']+f['b']/t


def run():
    endpoint = read(ENDPOINT)
    cap = read(C.TAG)
    z = C.geometry()
    mp.mp.dps = 80
    q = mp.matrix([[high(v) for v in row] for row in z['exact_gram']])
    eig, vec = mp.eigsy(q)
    assert eig[0] < 0 and all(v > 0 for v in eig[1:])
    transform = vec*mp.diag([mp.sqrt(abs(v)) for v in eig])
    h = mp.matrix([[high(v) for v in cap['common_center']]])*transform
    h = h/mp.sqrt(h[0]**2-sum(h[i]**2 for i in range(1, 10)))
    if h[0] < 0:
        for i in range(10):
            transform[i, 0] *= -1
        h[0] *= -1
    boost = mp.eye(10)
    boost[0, 0] = h[0]
    for i in range(1, 10):
        boost[0, i] = boost[i, 0] = -h[i]
        for j in range(1, 10):
            boost[i, j] += h[i]*h[j]/(h[0]+1)
    transform = transform*boost
    normals = np.array((mp.matrix([[high(v) for v in row] for row in z['normals']])*transform).tolist(), dtype=float)
    h = np.array([1.] + [0.]*9)
    assert np.max(abs(np.diag(normals @ METRIC @ normals.T)-2)) < 2e-5
    vertex_ids = S+[7, 8, 11]
    vn = normals[vertex_ids]
    vertex = h-(h @ METRIC @ vn.T) @ np.linalg.inv(vn @ METRIC @ vn.T) @ vn
    assert vertex @ METRIC @ vertex < 0
    groups, census = [], []
    for edge in EDGES:
        domain = next(r for r in endpoint['endpoint_domains'] if r['sides'] == list(edge))
        profiles = [r['labels'] for r in domain['profiles']['finite_label_profiles']]
        for tail in domain['profiles']['tails']:
            for m in TAIL_LABELS:
                ls = [2]*4
                ls[tail['slot']] = m
                profiles.append(ls)
        frame_ids = S+list(edge)
        frame = normals[frame_ids]
        gram = frame @ METRIC @ frame.T
        inverse = np.linalg.inv(gram)
        k0 = h-(h @ METRIC @ frame.T) @ inverse @ frame
        k0 /= math.sqrt(-k0 @ METRIC @ k0)
        kernel = null_space(frame @ METRIC)
        k1 = kernel[:, 0]+(kernel[:, 0] @ METRIC @ k0)*k0
        if k1 @ METRIC @ k1 < 1e-10:
            k1 = kernel[:, 1]+(kernel[:, 1] @ METRIC @ k0)*k0
        k1 /= math.sqrt(k1 @ METRIC @ k1)
        assert k0[0] > 0
        group = []
        for profile_index, labels in enumerate(profiles):
            contacts = np.array([0.]*4+[-2*math.cos(math.pi/m) if m > 2 else 0. for m in labels])
            p = contacts @ inverse @ frame
            residual = 2-p @ METRIC @ p
            assert residual > 0
            for sign in [-1, 1]:
                a = math.sqrt(residual)/2*(k0+sign*k1)
                b = math.sqrt(residual)/2*(-k0+sign*k1)
                other = [j for j in range(63) if j not in frame_ids]
                pa, aa, ba = [v @ METRIC @ normals[other].T for v in [p, a, b]]
                polys = np.array([[av, pv+2, bv] for av, pv, bv in zip(aa, pa, ba)]
                                 + [[av, 0, -bv] for av, bv in zip(aa, ba)]
                                 + [[a @ METRIC @ vertex, p @ METRIC @ vertex, b @ METRIC @ vertex],
                                    [-a[0], 0, b[0]]])
                domains = intervals(polys)
                if domains:
                    group.append(dict(profile_index=profile_index, labels=labels, branch=sign,
                                      frame=frame_ids, p=p, a=a, b=b, intervals=domains,
                                      residual=residual))
        groups.append(group)
        census.append(dict(edge=list(edge),profiles=len(profiles),branches=len(profiles)*2,
                           surviving_branches=len(group),intervals=sum(len(f['intervals']) for f in group)))
        print('SPLIT_BRIDGE_DOMAINS', census[-1], flush=True)
    hits, count = [], Counter()
    fixed = normals[S+[30]]
    fi = np.linalg.inv(fixed @ METRIC @ fixed.T)
    for ia, a in enumerate(groups[0]):
        for lo, hi in a['intervals']:
            for t in samples(lo, hi):
                va = evaluate(a, t)
                for ib, b in enumerate(groups[1]):
                    aa, bb, cc = [float(va @ METRIC @ b[key]) for key in ['a', 'b', 'p']]
                    for m in MUTUAL_LABELS:
                        target = -2*math.cos(math.pi/m) if m > 2 else 0.
                        for u in roots(aa, cc-target, bb):
                            count['mutual_equation_roots'] += 1
                            if not any(u > x+1e-8*max(1., x) and u < y-1e-8*max(1., y)
                                       if math.isfinite(y) else u > x+1e-8*max(1., x)
                                       for x, y in b['intervals']):
                                continue
                            vb = evaluate(b, u)
                            whole = np.vstack([fixed, va, vb])
                            eigen = float(np.linalg.eigvalsh(whole @ METRIC @ whole.T)[0])
                            if eigen < 1e-7:
                                count['mutual_corner_not_positive_definite'] += 1
                                continue
                            newvertex = h-(h @ METRIC @ whole.T) @ np.linalg.inv(whole @ METRIC @ whole.T) @ whole
                            products = newvertex @ METRIC @ normals.T
                            if newvertex @ METRIC @ newvertex >= -1e-9 or max(products) > 1e-7:
                                count['mutual_corner_not_in_chamber'] += 1
                                continue
                            for normal, family in [(va, a), (vb, b)]:
                                other = [j for j in range(63) if j not in family['frame']]
                                assert max(normal @ METRIC @ normals[other].T) < -2+2e-6
                            count['finite_bridge_hits'] += 1
                            if len(hits) < 100:
                                hits.append(dict(family_indices=[ia, ib], parameters=[t, u], mutual_label=m,
                                                 mutual_product=float(va @ METRIC @ vb),
                                                 minimum_corner_eigenvalue=eigen,
                                                 normals=[va.tolist(), vb.tolist()],
                                                 new_vertex=newvertex.tolist(),
                                                 new_vertex_norm=float(newvertex @ METRIC @ newvertex)))
        print('SPLIT_BRIDGE_FAMILY', ia, len(groups[0]), dict(count), flush=True)
    serial_groups = []
    for group in groups:
        serial_groups.append([{**f, **{key:f[key].tolist() for key in ['p','a','b']},
                               'intervals':[[lo, None if math.isinf(hi) else hi] for lo,hi in f['intervals']]}
                              for f in group])
    paths = [f'scripts/{TAG}.py', f'scripts/{C.TAG}.py', f'results/{C.TAG}.json', f'results/{ENDPOINT}.json']
    return dict(status='COMPLETED_NUMERICAL_SPLIT_ENDPOINT_BRIDGE_SCOUT',
                sources={path:sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
                retained_walls=63,edges=[list(e) for e in EDGES],tail_labels=TAIL_LABELS,
                mutual_labels=MUTUAL_LABELS,domain_census=census,families=serial_groups,
                counts=dict(count),hits=hits,metric=[-1]+[1]*9,
                old_normals=normals.tolist(),
                root_basis_to_lorentz=np.array(transform.tolist(),dtype=float).tolist(),
                scope='Floating-point discovery, with all unprescribed old contacts strictly dotted. Bounded tail labels and samples in the first parameter; the second parameter solves a mutual Coxeter equation. No exact acceptance, continuous nonexistence, or finite-volume claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    out = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'], result['counts'], 'saved hits', len(result['hits']), flush=True)
