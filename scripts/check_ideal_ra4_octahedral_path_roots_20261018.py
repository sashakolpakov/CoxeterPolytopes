#!/usr/bin/env python3
"""Exact finite D4 facts used by the restricted octahedral-chain proof.

This does not verify the universal cusp, surface, or cutting arguments.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import fcntl
import json
from pathlib import Path
import signal
import run_lean_batches_signal_safe_20261018 as safe

ROOT = Path(__file__).resolve().parents[1]


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def lorentz(x, y):
    return dot(x[1:], y[1:])-x[0]*y[0]


def reflected(x, n):
    coefficient = 2*lorentz(x, n)
    return tuple(a-coefficient*b for a, b in zip(x, n))


def run():
    roots = set()
    for i, j in combinations(range(4), 2):
        for a, b in product((-1, 1), repeat=2):
            r = [0]*4
            r[i], r[j] = a, b
            roots.add(tuple(r))
    roots = sorted(roots)
    vertices = []
    for i in range(4):
        for a in (-1, 1):
            v = [F(0)]*4
            v[i] = F(a)
            vertices.append(tuple(v))
    vertices += [tuple(F(a, 2) for a in v)
                 for v in product((-1, 1), repeat=4)]
    assert len(roots) == len(vertices) == 24
    assert all(dot(v, v) == 1 for v in vertices)
    faces = {r: {v for v in vertices if dot(r, v) == 1} for r in roots}
    assert all(len(vs) == 6 for vs in faces.values())
    assert all(dot(r, v) <= 1 for r in roots for v in vertices)
    counts, reflected_oct_vertices = Counter(), 0
    for r, s in combinations(roots, 2):
        p = dot(r, s)
        assert p in (-2, -1, 0, 1)
        assert len(faces[r] & faces[s]) == {-2: 0, -1: 0, 0: 1, 1: 3}[p]
        counts[p] += 1
        if p == -2:
            assert tuple(-a for a in r) == s
        if p != -1:
            continue
        c = tuple(a+b for a, b in zip(r, s))
        assert c in faces and dot(c, r) == dot(c, s) == 1
        nr, ns, nc = (1,)+r, (1,)+s, (1,)+c
        assert lorentz(nr, nc) == lorentz(ns, nc) == 0
        for cut, other in ((nr, ns), (ns, nr)):
            images = [reflected((F(1),)+v, cut) for v in faces[c]]
            assert all(lorentz(x, x) == lorentz(x, nc) == 0 for x in images)
            assert all(lorentz(x, cut) >= 0 and lorentz(x, other) <= 0
                       for x in images)
            center = tuple(sum(x[i] for x in images)/6 for i in range(5))
            assert lorentz(center, center) < 0
            assert lorentz(center, cut) > 0 and lorentz(center, other) < 0
            reflected_oct_vertices += len(images)
    assert dict(counts) == {-2: 12, -1: 96, 0: 72, 1: 96}
    assert reflected_oct_vertices == 1152
    return {'status': 'PASS_OCTAHEDRAL_PATH_D4_ROOT_FACTS',
            'roots': 24, 'vertices': 24, 'root_pairs': 276,
            'dot_product_counts': dict(sorted(counts.items())),
            'common_neighbor_cases': 96,
            'reflected_octahedral_vertices': reflected_oct_vertices,
            'universal_geometry_verified': False}


if __name__ == '__main__':
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        safe.resources()
        signal.alarm(59)
        result = run()
        signal.alarm(0)
        safe.resources(start=False)
    print(json.dumps(result, sort_keys=True))
