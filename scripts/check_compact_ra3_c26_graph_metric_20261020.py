#!/usr/bin/env python3
"""Fixed C26 graph and Q(sqrt17) identities, not a fullerene enumeration.

The count-one classification and actual Andreev symmetry bridge are written
premises, not checked here. Parent wrapper owns lock/resources. Root runs jobs.
"""
import argparse
import hashlib
import json
import os
import signal
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'results/compact_ra3_c26_graph_metric_20261020.json'
CAP = 65536
PINS = {
    'results/compact_ra3_twenty_six_vertices_no_regular_pentagon_20261020.md':
        '1bbe6adc630ff2f44424f6ef3c573b3c71dcef4f25496def9b13fd35b2e4a03f',
    'results/compact_ra3_c26_graph_metric_preparation_20261020.md':
        '09560b0d19d188bff51b5acd26c852256cbd34f6483da92a303ea4ea6aa3b35a',
    'manuscripts/compact_ra4_volume_constraints.tex':
        'afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
}


class K:
    """Exact a+b*sqrt17; no floating arithmetic or number-field helper."""
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        if isinstance(a, K):
            assert b == 0
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = F(a), F(b)

    def __add__(self, other):
        q = K(other)
        return K(self.a + q.a, self.b + q.b)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-K(other))

    def __rsub__(self, other):
        return K(other) + (-self)

    def __mul__(self, other):
        q = K(other)
        return K(self.a*q.a + 17*self.b*q.b,
                 self.a*q.b + self.b*q.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        q = K(other)
        norm = q.a*q.a - 17*q.b*q.b
        assert norm != 0
        return self * K(q.a/norm, -q.b/norm)

    def __rtruediv__(self, other):
        return K(other) / self

    def __pow__(self, exponent):
        assert isinstance(exponent, int) and exponent >= 0
        ans, factor = K(1), self
        while exponent:
            if exponent & 1:
                ans = ans * factor
            factor = factor * factor
            exponent //= 2
        return ans

    def __eq__(self, other):
        if not isinstance(other, (K, int, F)):
            return False
        q = K(other)
        return self.a == q.a and self.b == q.b

    def sign(self):
        a, b = self.a, self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0 or (a > 0) == (b > 0):
            return 1 if b > 0 else -1
        d = a*a - 17*b*b
        assert d != 0  # sqrt17 is irrational; neither a nor b vanishes here.
        return ((d > 0) - (d < 0)) * (1 if a > 0 else -1)

    def encoded(self):
        return [str(self.a), str(self.b)]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def connected(vertices, adjacency):
    assert vertices
    seen, todo = set(), [min(vertices)]
    while todo:
        v = todo.pop()
        if v not in seen:
            seen.add(v)
            todo.extend(adjacency[v] - seen)
    return seen == set(vertices)


def check():
    sources = {path: digest(ROOT/path) for path in PINS}
    assert sources == PINS, 'pinned proof/control changed'
    own = Path(__file__).resolve().relative_to(ROOT).as_posix()
    sources[own] = digest(ROOT/own)
    families = ['A', 'Ap', 'B', 'Bp', 'H']
    names = [f'{family}{i}' for family in families for i in range(3)]
    index = {v: i for i, v in enumerate(names)}

    def name(family, i):
        return f'{family}{i % 3}'

    triangles = set()

    def triangle(*vertices):
        assert len(set(vertices)) == 3
        item = tuple(sorted(vertices))
        assert item not in triangles
        triangles.add(item)

    for a, b in [('A', 'B'), ('Ap', 'Bp')]:
        triangle(*(name(a, i) for i in range(3)))
        for i in range(3):
            triangle(name(a, i), name(a, i+1), name(b, i))
            triangle(name(a, i), name(b, i), name('H', i))
            triangle(name(a, i), name('H', i), name(b, i-1))
    for i in range(3):
        triangle(name('B', i), name('Bp', i), name('H', i))
        triangle(name('B', i), name('Bp', i), name('H', i+1))
    assert len(names) == 15 and len(triangles) == 26
    assert set().union(*(set(t) for t in triangles)) == set(names)
    edge_triangles, adjacency = {}, {v: set() for v in names}
    for t in sorted(triangles):
        for e in combinations(t, 2):
            edge_triangles.setdefault(e, []).append(t)
            adjacency[e[0]].add(e[1])
            adjacency[e[1]].add(e[0])
    assert len(edge_triangles) == 39
    assert all(len(t) == 2 for t in edge_triangles.values())
    assert connected(set(names), adjacency)
    assert 15 - len(edge_triangles) + len(triangles) == 2
    links = {}
    for v in names:
        link_edges = [tuple(x for x in t if x != v)
                      for t in sorted(triangles) if v in t]
        link = {x: set() for x in adjacency[v]}
        for x, y in link_edges:
            link[x].add(y)
            link[y].add(x)
        degree = 6 if v.startswith('H') else 5
        assert len(link) == len(link_edges) == degree
        assert all(len(neighbors) == 2 for neighbors in link.values())
        assert connected(set(link), link)
        links[v] = [list(e) for e in link_edges]

    def permute(kind):
        out = []
        for f in families:
            for i in range(3):
                if kind == 'rho':
                    g, j = f, i+1
                elif kind == 'sigma':
                    g = {'A': 'Ap', 'Ap': 'A', 'B': 'Bp', 'Bp': 'B', 'H': 'H'}[f]
                    j = i
                else:
                    assert kind == 'tau'
                    g, j = f, -i-1 if f in ('B', 'Bp') else -i
                out.append(index[name(g, j)])
        assert sorted(out) == list(range(15))
        return tuple(out)

    generators = {kind: permute(kind) for kind in ('rho', 'sigma', 'tau')}
    identity = tuple(range(15))
    group, todo = {identity}, [identity]
    while todo:
        p = todo.pop()
        for q in generators.values():
            pq = tuple(p[q[i]] for i in range(15))
            if pq not in group:
                group.add(pq)
                todo.append(pq)
                assert len(group) <= 12
    assert len(group) == 12
    for p in group:
        assert {tuple(sorted(names[p[index[v]]] for v in t))
                for t in triangles} == triangles
    orbits = {v: sorted({names[p[index[v]]] for p in group})
              for v in ('A0', 'B0', 'H0')}
    assert set(orbits['A0']) == {name(f, i) for f in ('A', 'Ap') for i in range(3)}
    assert set(orbits['B0']) == {name(f, i) for f in ('B', 'Bp') for i in range(3)}
    assert set(orbits['H0']) == {name('H', i) for i in range(3)}

    def thirds(x, y):
        edge = tuple(sorted((x, y)))
        return {next(v for v in t if v not in edge) for t in edge_triangles[edge]}

    for i in range(3):
        for a, b in [('A', 'B'), ('Ap', 'Bp')]:
            assert thirds(name(a, i), name(a, i+1)) == {name(a, i-1), name(b, i)}
        assert thirds(name('B', i), name('Bp', i)) == {name('H', i), name('H', i+1)}

    root = K(0, 1)
    assert root*root == 17 and (root-4).sign() > 0 and (5-root).sign() > 0
    x = (root-1)/6
    assert 9*x*x+3*x-4 == 0 and (x-F(1, 3)).sign() > 0
    r2, u2 = 1+x, K(F(2, 3))
    a2 = 2*(1+x)/(3*x)
    upper_a2 = a2+F(1, 3)
    b2 = (1+x)/(2*(3*x-1))
    v2, upper_b2 = b2+F(1, 2), K(F(1, 2))
    values = {'h_squared': x, 'R_squared': r2, 'U_squared': u2,
              'a_squared': a2, 'A_squared': upper_a2, 'b_squared': b2,
              'V_squared': v2, 'B_squared': upper_b2}
    assert all(v.sign() > 0 for v in values.values())
    assert -x+r2 == 1
    assert -a2+u2+upper_a2 == 1
    assert -b2+v2+upper_b2 == 1
    assert u2*r2 == a2*x
    assert v2*r2 == 4*b2*x
    assert upper_a2 == a2+u2/2
    assert upper_b2 == v2-b2
    assert upper_a2*upper_b2 == a2*b2/(1+x)**2
    edge_a2 = 9*a2*b2*x*x/(1+x)**2
    edge_b = (1+3*x)/2
    assert edge_a2 == 3*x/(3*x-1) == (7+root)/4
    assert edge_b == (1+root)/4
    assert (edge_a2-F(11, 4)).sign() > 0
    assert (K(F(3, 2))-edge_b).sign() > 0
    assert F(4) < F(5) < F(25, 4)  # 3/2<phi; phi^2<11/4.
    assert all(digest(ROOT/path) == pin for path, pin in sources.items()), 'source changed'
    return {
        'status': 'EXACT_FIXED_C26_GRAPH_AND_METRIC_IDENTITIES_PASS',
        'scope': 'Fixed15-vertex dual graph and exact metric identities; NOT a fullerene enumeration or Andreev geometry check.',
        'external_premises': ['Primary BGM2012v2 Table2 count-one fullerene at26',
                              'Actual compact RA3 metric uniqueness and symmetry realization'],
        'owns_lock': False, 'sources': sources,
        'field_encoding': '[a,b] means a+b*sqrt17, exact rational strings',
        'vertices': names, 'edges': [list(e) for e in sorted(edge_triangles)],
        'triangles': [list(t) for t in sorted(triangles)], 'vertex_links': links,
        'euler_characteristic': 2,
        'D3h_generators': {k: list(v) for k, v in generators.items()},
        'D3h_order': len(group), 'face_orbits': orbits,
        'positive_squared_parameters': {k: v.encoded() for k, v in values.items()},
        'A_orbit_witness_edge_cosh_squared': edge_a2.encoded(),
        'B_orbit_witness_edge_cosh': edge_b.encoded(),
        'strict_comparisons': ['Aedge_cosh_squared>11/4>phi_squared',
                               'Bedge_cosh<3/2<phi'],
        'enumerations_or_hyperbolic_realization_searches': 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.write and CERT.exists():
        parser.error('refusing to overwrite certificate')
    assert int(os.environ.get('COXETER_CHECK_OWNER', '-1')) == os.getppid(), 'bounded parent wrapper required'
    signal.alarm(45)
    result = check()
    payload = json.dumps(result, sort_keys=True, indent=2)+'\n'
    assert len(payload.encode()) <= CAP, 'certificate exceeds64KiB'
    if args.write:
        with CERT.open('x') as stream:
            stream.write(payload)
    else:
        assert CERT.stat().st_size <= CAP
        assert result == json.loads(CERT.read_text()), 'saved certificate differs'
    print(json.dumps({'status': 'PASS', 'mode': 'write' if args.write else 'verify',
                      'vertices': 15, 'edges': 39, 'triangles': 26, 'group_order': 12,
                      'certificate_bytes': len(payload.encode()),
                      'certificate_sha256': hashlib.sha256(payload.encode()).hexdigest()}))


if __name__ == '__main__':
    main()
