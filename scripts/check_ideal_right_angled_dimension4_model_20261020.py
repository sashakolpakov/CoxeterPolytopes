#!/usr/bin/env python3
"""Exact finite model checks for the cusp-height dimension argument.

This checks the standard 24-cell and the scalar identity used in its
spherical-code reduction. It does not prove Musin's theorem, Kolpakov's
classification, or the universal geometric reduction in the manuscript.
Run through run_bounded_research_check_20261018.py --parent-lock.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
import json


def dot(x, y):
    return -x[0] * y[0] + sum(a * b for a, b in zip(x[1:], y[1:]))


def spatial_dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def main():
    normals = []
    for i, j in combinations(range(4), 2):
        for si, sj in product((-1, 1), repeat=2):
            r = [Q(0)] * 4
            r[i], r[j] = Q(si), Q(sj)
            normals.append((Q(1), *r))
    assert len(set(normals)) == 24
    assert all(dot(n, n) == 1 for n in normals)
    gram = Counter(dot(a, b) for a, b in combinations(normals, 2))
    assert gram == Counter({Q(0): 96, Q(-1): 72, Q(-2): 96, Q(-3): 12})
    assert all(spatial_dot(a[1:], b[1:]) <= 1
               for a, b in combinations(normals, 2))

    opposite = []
    for i, n in enumerate(normals):
        partner = (Q(1), *(-x for x in n[1:]))
        j = normals.index(partner)
        assert tuple(a + b for a, b in zip(n, partner)) == (2, 0, 0, 0, 0)
        if i < j:
            opposite.append((i, j))
    assert len(opposite) == 12

    vertices = []
    for i in range(4):
        for sign in (-1, 1):
            x = [Q(0)] * 4
            x[i] = Q(sign)
            vertices.append((Q(1), *x))
    vertices.extend((Q(1), *(Q(s, 2) for s in signs))
                    for signs in product((-1, 1), repeat=4))
    assert len(set(vertices)) == 24
    links = []
    for v in vertices:
        assert dot(v, v) == 0
        assert all(dot(v, n) <= 0 for n in normals)
        incident = [i for i, n in enumerate(normals) if dot(v, n) == 0]
        assert len(incident) == 6
        link_gram = Counter(dot(normals[i], normals[j])
                            for i, j in combinations(incident, 2))
        assert link_gram == Counter({Q(0): 12, Q(-1): 3})
        links.append(incident)
    assert all(sum(i in link for link in links) == 6 for i in range(24))

    # Every displayed wall is an actual facet: this positive multiple of
    # its unit facet-center is feasible and strict at every other wall.
    for i, n in enumerate(normals):
        facet_point = (Q(2), *n[1:])
        assert dot(facet_point, facet_point) < 0
        assert dot(facet_point, n) == 0
        assert all(dot(facet_point, m) < 0 for j, m in enumerate(normals) if j != i)

    # Polynomial coefficient checks, variables A=a^2 and B=b^2:
    # (1+A)(1+B)-4AB = (1-A)(1+B)+2A(1-B).
    left = {(0, 0): 1, (1, 0): 1, (0, 1): 1, (1, 1): -3}
    right = Counter({(0, 0): 1, (0, 1): 1})
    right.subtract({(1, 0): 1, (1, 1): 1})
    right.update({(1, 0): 2})
    right.subtract({(1, 1): 2})
    assert dict(right) == left
    print(json.dumps({
        'status': 'PASS_EXACT_24_CELL_MODEL',
        'arithmetic': 'integer and Fraction, no floating point',
        'unit_normals': 24, 'unordered_gram_pairs': 276,
        'opposite_normal_pairs': 12, 'displayed_null_vertices': 24,
        'cubical_cusp_links': 24, 'actual_facet_witnesses': 24,
        'scalar_polynomial_identity': True,
        'external_and_universal_geometry_verified': False,
        'vertex_enumeration_completeness_claimed': False,
    }, sort_keys=True))


if __name__ == '__main__':
    main()
