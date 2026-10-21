#!/usr/bin/env python3
"""Numerical proposals with exact strict escape certificates for selected cores."""
from fractions import Fraction as Q
import numpy as np
import mpmath as mp
from scipy.spatial import HalfspaceIntersection
from scripts.fresh_h9_replaced_core_20261013_base import (
    K, O, scale, value, strict_escape)


def volume_gate(I, center, selected, max_escapes=12):
    mp.mp.dps = 100
    C = np.array([[float(value(r['field'], t)) for t in r['c']] for r in selected])
    In = np.array([[float(value(K, t)) for t in row] for row in I])
    hn = np.array([float(value(K, t)) for t in center])
    hn /= hn.sum()
    assert np.min(C @ hn) > 0
    half = np.c_[C[:, 9, None]-C[:, :9], -C[:, 9]]
    half /= np.linalg.norm(half[:, :9], axis=1)[:, None]
    vertices = HalfspaceIntersection(half, hn[:9], qhull_options='Qx').intersections
    ys = np.c_[vertices, 1-vertices.sum(axis=1)]
    norms = np.einsum('bi,ij,bj->b', ys, In, ys)
    escapes, supports = [], set()
    for index in sorted(np.flatnonzero(norms > 1e-8), key=lambda i: -norms[i]):
        support = tuple(int(i) for i in np.flatnonzero(ys[index] > 1e-7))
        if support in supports:
            continue
        cert = None
        for epsilon in (1e-3, 1e-5, 1e-7):
            yn = ys[index]+epsilon*hn
            for denominator in (10**9, 10**12, 10**15):
                y = tuple(scale(O, Q(round(float(t)*denominator), denominator)) for t in yn)
                try:
                    cert = strict_escape(y, I, center, selected)
                except AssertionError:
                    continue
                break
            if cert:
                break
        if cert:
            escapes.append({'certificate': cert, 'proposal_support': list(support),
                            'proposal_vertex_norm': float(norms[index])})
            supports.add(support)
        if len(escapes) >= max_escapes:
            break
    return {'status': 'EXACT_STRICT_ESCAPES' if escapes else 'UNRESOLVED_VOLUME',
            'numeric_vertices': len(ys), 'numeric_spacelike_vertices': int((norms > 1e-8).sum()),
            'maximum_numeric_vertex_norm': float(norms.max()), 'escapes': escapes,
            'scope': 'Only the exact feasible future spacelike certificates prove '
                     'infinite volume. No absence-of-escape volume conclusion.'}
