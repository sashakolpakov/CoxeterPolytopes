#!/usr/bin/env python3
"""Shared intrinsic A6+A6+A7 cusp geometry; floating proposals, exact base data.

Coordinates use a rational unit-root basis. Numerical vectors are mapped to
Euclidean orthonormal coordinates by a Cholesky factor. No arithmetic lattice
restriction is imposed on prospective caps.
"""
from fractions import Fraction as F
from itertools import product
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_balanced_20261014'
BENCHMARK = 'results/fresh_h20_intrinsic_cusp_benchmark_20261014.json'


def sha(p):
    return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def exact_data():
    # A_n finite unit Gram inverse: 2 min(i,j)(n+1-max(i,j))/(n+1).
    ranks = [6, 6, 7]
    gram = [[F(0) for _ in range(19)] for _ in range(19)]
    normals, factor_vertices, factors = [], [], []
    basis_start = wall_start = 0
    for n in ranks:
        for i in range(n):
            gram[basis_start+i][basis_start+i] = F(1)
            if i+1 < n:
                gram[basis_start+i][basis_start+i+1] = F(-1, 2)
                gram[basis_start+i+1][basis_start+i] = F(-1, 2)
        inverse = [[F(2*min(i+1, j+1)*(n-max(i, j)), n+1) for j in range(n)] for i in range(n)]
        block = [[F(-1)]*n]+[[F(i == j) for j in range(n)] for i in range(n)]
        normals += [[F(0)]*basis_start+v+[F(0)]*(19-basis_start-n) for v in block]
        vv = []
        for omit in range(n+1):
            rhs = [F(1)-(n+1 if omit == i+1 else 0) for i in range(n)]
            z = [sum(inverse[i][j]*rhs[j] for j in range(n)) for i in range(n)]
            vv.append([F(0)]*basis_start+z+[F(0)]*(19-basis_start-n))
        factor_vertices.append(vv)
        factors.append(list(range(wall_start, wall_start+n+1)))
        basis_start += n
        wall_start += n+1
    return dict(ranks=ranks, spatial_Gram=gram, seed_vectors=normals,
                factor_vertices=factor_vertices, factors=factors,
                marks=[[1]*(n+1) for n in ranks])


def data(scales=(1., 1., 1.)):
    e = exact_data()
    B = np.array(e['spatial_Gram'], dtype=float)
    L = np.linalg.cholesky(B)
    coeff = np.array(e['seed_vectors'], dtype=float)
    V = coeff@L
    offsets = np.concatenate([np.full(n+1, s) for n, s in zip(e['ranks'], scales)])
    fv = [np.array(v, dtype=float)@L*scales[k] for k, v in enumerate(e['factor_vertices'])]
    vertices = []
    for local in product(*[range(n+1) for n in e['ranks']]):
        z = sum(fv[k][a] for k, a in enumerate(local))
        omitted = [e['factors'][k][a] for k, a in enumerate(local)]
        retained = [i for i in range(22) if i not in omitted]
        gaps = offsets-V@z
        assert np.min(gaps) > -1e-9 and np.max(abs(gaps[retained])) < 1e-9
        assert np.max(abs(gaps[omitted]-np.array([7, 7, 8])*scales)) < 1e-9
        vertices.append(dict(local_omitted=list(local), omitted=omitted, omitted_walls=omitted,
                             retained=retained, point=z, heights=gaps[omitted]))
    assert len(vertices) == 392
    return dict(factors=e['factors'], walls=list(range(22)), ranks=e['ranks'],
                marks=e['marks'], unit_Gram=V@V.T, vectors=V, offsets=offsets,
                vertices=vertices, scales=list(scales), spatial_Gram=B,
                factor_vertices=fv, cholesky=L)


if __name__ == '__main__':
    d = data()
    verts = np.array([v['point'] for v in d['vertices']])
    assert np.max(abs(np.sum(verts**2, axis=1)-196)) < 1e-9
    print('BALANCED_GEOMETRY', len(verts), 'vertices; common square196; no volume claim')
