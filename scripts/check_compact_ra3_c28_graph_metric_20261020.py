#!/usr/bin/env python3
"""Fixed C28 graph/Q(sqrt5) identities; written geometry remains an input.

No fullerene enumeration, realization search or old-worker import. The
bounded parent owns lock/resources; root alone runs --write and --verify.
"""
import argparse
import hashlib
import json
import os
import signal
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'results/compact_ra3_c28_graph_metric_20261020.json'
CAP = 131072
PINS = {
    'results/compact_ra3_twenty_eight_vertices_no_regular_pentagon_20261020.md':
        '076864a5d836456e56ac377dffbcaa581a4738ec977dbf58a7905580d9838afc',
    'results/compact_ra3_c28_graph_metric_preparation_20261020.md':
        'fa34cba3412e73685cf89901bcc6162df37ce40aadb17be584890fdec4e8ca56',
    'results/review_compact_ra3_c28_proof_root_20261020.md':
        '200863417476eaa93dce3f83040b06dfc684f995b03f24eb51c118d6ebadaa33',
    'results/review_compact_ra3_c28_finite57_20261020.md':
        'c1082bf6cc5040e4db426eca2c4af3c6c7a7d5b4fa258119253392400b8a4e33',
    'manuscripts/compact_ra4_volume_constraints.tex':
        'afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f',
    'results/compact_ra3_twenty_six_vertices_no_regular_pentagon_20261020.md':
        '1bbe6adc630ff2f44424f6ef3c573b3c71dcef4f25496def9b13fd35b2e4a03f',
    'scripts/check_compact_ra3_c26_graph_metric_20261020.py':
        '47108a5ff3f732bf4b54d7893cc925f8e2a9423106c0ac812f690d355a9976e4',
    'results/compact_ra3_c26_graph_metric_20261020.json':
        '3f5567187bdaa7f2f42c6d4064a8cb46d09f30c99cdf12f667be050b3261cd1b',
    'results/review_compact_ra3_c26_proof_and_checker_root_20261020.md':
        'c618874b8293d89560aa38eca03a9a751459a32fd07bb82c7e9207511d7fac61',
    'results/compact_ra4_farthest_dodecahedron_wall_radius_20261020.md':
        'd082e4a918d4504a59ac200f32bbb141f1f2cb52d701006eaea7ba5b0e5bcc22',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
}


class K:
    """Exact a+b*sqrt5 with rational coefficients and exact physical sign."""
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        if isinstance(a, K):
            assert b == 0
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = F(a), F(b)

    def __add__(self, other):
        q = K(other)
        return K(self.a+q.a, self.b+q.b)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-K(other))

    def __rsub__(self, other):
        return K(other) + (-self)

    def __mul__(self, other):
        q = K(other)
        return K(self.a*q.a+5*self.b*q.b, self.a*q.b+self.b*q.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        q = K(other)
        norm = q.a*q.a-5*q.b*q.b
        assert norm != 0
        return self*K(q.a/norm, -q.b/norm)

    def __rtruediv__(self, other):
        return K(other)/self

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        answer, factor = K(1), self
        while n:
            if n & 1:
                answer *= factor
            factor *= factor
            n //= 2
        return answer

    def __eq__(self, other):
        if not isinstance(other, (K, int, F)):
            return False
        q = K(other)
        return self.a == q.a and self.b == q.b

    def sign(self):
        a, b = self.a, self.b
        if b == 0:
            return (a > 0)-(a < 0)
        if a == 0 or (a > 0) == (b > 0):
            return 1 if b > 0 else -1
        delta = a*a-5*b*b
        assert delta != 0
        return ((delta > 0)-(delta < 0))*(1 if a > 0 else -1)

    def encoded(self):
        return [str(self.a), str(self.b)]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def snapshot():
    sources = {path: digest(ROOT/path) for path in PINS}
    assert sources == PINS, 'pinned proof/control changed'
    own = Path(__file__).resolve().relative_to(ROOT).as_posix()
    sources[own] = digest(ROOT/own)
    return sources


def fresh(sources):
    assert all(digest(ROOT/path) == pin for path, pin in sources.items()), 'source changed'


def connected(vertices, adjacency):
    assert vertices
    seen, todo = set(), [min(vertices)]
    while todo:
        v = todo.pop()
        if v not in seen:
            seen.add(v)
            todo.extend(adjacency[v]-seen)
    return seen == set(vertices)


def complex_check(names, raw):
    triangles = {tuple(sorted(t)) for t in raw}
    assert len(names) == len(set(names)) == 16 and len(raw) == len(triangles) == 28
    assert all(len(t) == len(set(t)) == 3 for t in triangles)
    assert set().union(*(set(t) for t in triangles)) == set(names)
    edges, adjacency = {}, {v: set() for v in names}
    for t in sorted(triangles):
        for edge in combinations(t, 2):
            edges.setdefault(edge, []).append(t)
            x, y = edge
            adjacency[x].add(y)
            adjacency[y].add(x)
    assert len(edges) == 42 and all(len(ts) == 2 for ts in edges.values())
    assert connected(set(names), adjacency) and 16-42+28 == 2
    links = {}
    for v in names:
        link_edges = [tuple(x for x in t if x != v) for t in sorted(triangles) if v in t]
        link = {x: set() for x in adjacency[v]}
        for x, y in link_edges:
            link[x].add(y)
            link[y].add(x)
        degree = 6 if v.startswith('H') else 5
        assert len(link_edges) == len(link) == degree
        assert all(len(ns) == 2 for ns in link.values()) and connected(set(link), link)
        links[v] = [list(e) for e in link_edges]
    assert sum(len(ns) == 5 for ns in adjacency.values()) == 12
    assert sum(len(ns) == 6 for ns in adjacency.values()) == 4
    result = {'vertices': names, 'triangles': [list(t) for t in sorted(triangles)],
              'edges': [list(e) for e in sorted(edges)], 'link_edges': links,
              'euler_characteristic': 2}
    return triangles, edges, adjacency, result


def automorphism(p, names, triangles):
    assert set(p) == set(names) == set(p.values())
    assert {tuple(sorted(p[v] for v in t)) for t in triangles} == triangles


def graph_checks():
    p = lambda i, j: f'P{i}{j}'
    names = [f'H{i}' for i in range(4)]+[p(i, j) for i in range(4) for j in range(4) if i != j]
    raw = [tuple(p(i, j) for j in range(4) if j != i) for i in range(4)]
    for k in range(4):
        i, j, l = [v for v in range(4) if v != k]
        boundary = [p(i, j), p(j, i), p(j, l), p(l, j), p(l, i), p(i, l)]
        raw.extend((f'H{k}', boundary[t], boundary[(t+1) % 6]) for t in range(6))
    triangles, edges, adj, tetra = complex_check(names, raw)
    group = []
    for perm in permutations(range(4)):
        mapping = {f'H{i}': f'H{perm[i]}' for i in range(4)}
        mapping.update({p(i, j): p(perm[i], perm[j]) for i in range(4) for j in range(4) if i != j})
        automorphism(mapping, names, triangles)
        group.append(mapping)
    assert len(group) == 24 and len({g['P01'] for g in group}) == 12
    for i in range(4):
        for j in range(i+1, 4):
            edge = tuple(sorted((p(i, j), p(j, i))))
            thirds = {next(v for v in t if v not in edge) for t in edges[edge]}
            assert thirds == {f'H{k}' for k in range(4) if k not in (i, j)}
    assert not any(x.startswith('H') and y.startswith('H') for x, y in edges)
    tetra['symmetry_order'] = 24
    tetra['all_pentagons_have_edge_cosh_4_over_3'] = True

    rows = ['H0 H1 A', 'H0 H1 B', 'H0 B C', 'H0 C D', 'H0 D E', 'H0 E A',
            'A E F', 'A F G', 'A G H1', 'H1 G J', 'H1 J K', 'H1 K B',
            'B K M', 'B M C', 'C M N', 'C N D', 'D N H2', 'D H2 E',
            'E H2 F', 'H2 F P', 'F G P', 'G J P', 'J H3 P', 'J K H3',
            'K M H3', 'M N H3', 'N H2 H3', 'H2 H3 P']
    names = ['H0', 'H1', 'H2', 'H3']+list('ABCDEFGJKMNP')
    triangles, edges, adj, second = complex_check(names, [tuple(r.split()) for r in rows])
    cycles = ['H1 B C D E A', 'H0 B K J G A', 'D E F P H3 N', 'J K M N H2 P',
              'H0 E F G H1', 'H0 C M K H1', 'H0 D N M B', 'H0 E H2 N C',
              'H0 A F H2 D', 'A G P H2 E', 'A H1 J P F', 'H1 K H3 P G',
              'H1 B M H3 J', 'B C N H3 K', 'C D H2 H3 M', 'F G J H3 H2']
    for v, row in zip(names, cycles):
        cyc = row.split()
        expected = {tuple(sorted((cyc[i], cyc[(i+1) % len(cyc)]))) for i in range(len(cyc))}
        assert expected == {tuple(e) for e in second['link_edges'][v]}
    hex_edges = [list(e) for e in edges if all(v.startswith('H') for v in e)]
    assert sorted(hex_edges) == [['H0', 'H1'], ['H2', 'H3']]

    def involution(pairs):
        result = {}
        for x, y in pairs:
            result[x], result[y] = y, x
        automorphism(result, names, triangles)
        assert all(result[result[v]] == v for v in names)
        return result

    alpha = involution([('H0', 'H1'), ('H2', 'H3'), ('A', 'B'), ('C', 'G'),
                        ('D', 'J'), ('E', 'K'), ('F', 'M'), ('N', 'P')])
    beta = involution([('H0', 'H2'), ('H1', 'H3'), ('A', 'N'), ('B', 'P'),
                       ('C', 'F'), ('D', 'E'), ('G', 'M'), ('J', 'K')])
    assert all(alpha[beta[v]] == beta[alpha[v]] for v in names)
    identity = {v: v for v in names}
    group = [identity, alpha, beta, {v: alpha[beta[v]] for v in names}]
    assert len({tuple(g[v] for v in names) for g in group}) == 4
    orbits = [{g[v] for g in group} for v in ('A', 'D', 'C')]
    assert orbits == [set('ABNP'), set('DEJK'), set('CFGM')]
    pentagons = set('ABCDEFGJKMNP')
    propagation = {}
    for representative, initial in zip(('A', 'D', 'C'), orbits):
        regular = set(initial)
        while True:
            enlarged = regular | {v for v in pentagons if len(adj[v] & regular) >= 2}
            enlarged |= set().union(*(o for o in orbits if o & enlarged))
            if enlarged == regular:
                break
            regular = enlarged
        assert regular == (set('DEJK') if representative == 'D' else pentagons)
        propagation[representative] = sorted(regular)
    assert len(adj['H0'] & pentagons) == 5
    assert alpha['C'] == 'G' and len({'A', 'C', alpha['C']}) == 3
    assert all(v in adj['D'] and v in adj['E'] for v in ('H0', 'H2'))
    assert 'H1' in adj['A'] and 'H1' not in adj['C']
    second.update({'alpha': alpha, 'beta': beta, 'group_order': 4,
                   'pentagon_orbits': [sorted(o) for o in orbits],
                   'two_regular_neighbor_closure': propagation,
                   'hexagon_induced_edges': sorted(hex_edges)})
    return {'tetrahedral': tetra, 'second': second,
            'nonisomorphism_witness': 'degree-six induced graphs have zero versus two edges'}


def determinant(matrix):
    if not matrix:
        return K(1)
    return sum(((-1)**j)*matrix[0][j]*determinant(
        [row[:j]+row[j+1:] for row in matrix[1:]]) for j in range(len(matrix)))


def metric_checks():
    phi = K(F(1, 2), F(1, 2))
    u = phi-1
    assert phi*phi == phi+1 and phi*u == 1
    assert (phi-F(4, 3)).sign() > 0 and (2-phi).sign() > 0
    # Tetrahedral equations use squares/products, avoiding auxiliary radicals.
    B2, A2, AB, a2, h2, R2 = F(3, 4), F(3), F(3, 2), F(7, 4), F(3, 4), F(7, 4)
    assert A2 == 4*B2 and AB == 2*B2 and (-B2/3) < 0
    assert a2 == A2-(2*AB+B2)/3 == 2*AB-(A2+B2)/3
    assert -a2+A2+B2-2*AB/3 == 1 and -h2+R2 == 1
    assert (A2+2*AB+B2)*R2 == 9*a2*h2
    assert h2+R2/3 == F(4, 3)
    vgram = [[K(1 if i == j else F(-1, 3)) for j in range(4)] for i in range(4)]
    assert determinant(vgram) == 0 and determinant([row[:3] for row in vgram[:3]]).sign() > 0
    assert all(sum(row) == 0 for row in vgram)

    basis = [[K(1), K(0), K(0), K(0)], [K(0), K(1), K(0), K(0)],
             [K(0), K(0), K(1), -phi], [K(0), K(0), -phi, K(1)]]
    assert determinant(basis) == -phi

    def dot(x, y):
        return sum(x[i]*basis[i][j]*y[j] for i in range(4) for j in range(4))

    D, E, H0, H2 = [[K(int(i == j)) for i in range(4)] for j in range(4)]
    A, C = [-phi, K(0), phi, K(1)], [K(0), -phi, phi, K(1)]
    L1, L3 = [-u, K(0), K(1), u], [K(0), -u, u, K(1)]
    assert [dot(A, v) for v in (E, H0, D, H2)] == [0, 0, -phi, -phi]
    assert [dot(C, v) for v in (D, H0, E, H2)] == [0, 0, -phi, -phi]
    assert dot(A, A) == dot(C, C) == 1
    assert dot(A, L1) == dot(A, E) == 0
    # These contractions give the two equations for arbitrary real k,e.
    assert dot(L1, L1) == dot(L3, L3) == -(u**3)
    assert dot(L1, E) == dot(L3, D) == 0
    assert dot(L1, L3) == -1 and dot(L1, D)+dot(E, L3) == -2*u
    coefficients = [phi**2, -(8-2*phi), 33-20*phi]
    assert coefficients == [phi**2, -2*phi-4*u**2, 1-4*u**5]
    roots = [phi**4, (6+7*phi)/29]
    leading = coefficients[2]
    assert leading.sign() > 0
    assert coefficients == [leading*roots[0]*roots[1], -leading*sum(roots), leading]
    assert (1-roots[1]).sign() > 0 and (roots[0]-1).sign() > 0
    k = phi**2
    e = (phi-k*k)/(2*u*k)
    assert e == -phi and e*e == 1+u**3*k*k
    H1 = [k*L1[i]+e*E[i] for i in range(4)]
    H3 = [k*L3[i]+e*D[i] for i in range(4)]
    H = [H0, H1, H2, H3]
    gram = [[dot(x, y) for y in H] for x in H]
    assert gram == [[1, 0, -phi, -phi**2], [0, 1, -phi**2, -phi],
                    [-phi, -phi**2, 1, 0], [-phi**2, -phi, 0, 1]]
    assert determinant([row[:3] for row in gram[:3]]) == 1-phi**2-phi**4
    assert determinant([row[:3] for row in gram[:3]]).sign() < 0
    assert determinant(gram) == 0
    relation = [1, -1, -1, 1]
    assert all(sum(relation[j]*H[j][i] for j in range(4)) == 0 for i in range(4))
    assert (2-2*phi).sign() < 0
    product_row = [K(0), K(0), -phi, -phi]
    assert [dot(A, h) for h in H] == [dot(C, h) for h in H] == product_row
    # Alpha swaps both zero entries and both -phi entries, so G has this row.
    assert [product_row[j] for j in (1, 0, 3, 2)] == product_row
    difference, midpoint = [A[i]-C[i] for i in range(4)], [(A[i]+C[i])/2 for i in range(4)]
    assert dot(difference, difference) == 2*phi**2
    assert all(dot(difference, h) == 0 for h in H) and dot(midpoint, difference) == 0
    assert determinant([H0, H1, H2, difference]) != 0
    assert dot(midpoint, midpoint)+dot(difference, difference)/4 == 1
    assert (1-dot(midpoint, midpoint))/dot(difference, difference) == F(1, 4)
    # Any unit vector with the common H row is midpoint +/- difference/2.
    assert F(65+17+16, 12) == F(49, 6) and F(120+0+16, 12) > F(49, 6)
    encode = lambda xs: [K(x).encoded() for x in xs]
    return {'field_encoding': '[a,b] means a+b*sqrt5, rational strings',
            'phi': phi.encoded(), 'tetrahedral_squared_parameters':
            {'B2': str(B2), 'A2': str(A2), 'AB': str(AB), 'a2': str(a2),
             'h2': str(h2), 'R2': str(R2), 'edge_cosh': '4/3'},
            'D_E_H0_H2_basis_gram': [encode(row) for row in basis],
            'A_C_coordinates': [encode(A), encode(C)],
            'H1_H3_k_linear_vectors': [encode(L1), encode(L3)],
            'quartic_coefficients_in_k_squared': encode(coefficients),
            'roots_in_k_squared': encode(roots), 'selected_k_e': encode([k, e]),
            'H_coordinates': [encode(h) for h in H], 'H_gram': [encode(row) for row in gram],
            'H_null_relation': relation, 'unit_normal_fiber_squared_parameter': '1/4',
            'duplicate_labels': ['A', 'C', 'G'], 'conditional_RA4_volume_coefficient': '49/6'}


def check(sources):
    graphs, algebra = graph_checks(), metric_checks()
    fresh(sources)
    return {'status': 'EXACT_FIXED_C28_GRAPHS_AND_METRIC_IDENTITIES_PASS',
            'scope': 'Two fixed dual spheres and exact algebra, not graph classification or geometric realization.',
            'external_and_written_premises': [
                'BGM arXiv1207.7010v2 Table2 count-two fullerene at28 (not fetched by worker)',
                'Actual compact RA3 metric uniqueness and isometric symmetry realization',
                'Polygon capacity/full L7 equality pattern and regular base length comparison',
                'Actual edge-normal length identity, pentagon regularity and three-edge closure',
                'Positive one-dimensional unit-normal fiber and distinct actual facet normals',
                'Inherited RA4 capacity, dodecahedral-facet existence, W gap and volume formula'],
            'primary_classification_url': 'https://arxiv.org/pdf/1207.7010v2',
            'owns_lock': False, 'sources': sources, 'graphs': graphs, 'algebra': algebra,
            'fullerene_enumerations_or_realization_searches': 0}


def read_certificate():
    assert CERT.stat().st_size <= CAP
    with CERT.open('rb') as stream:
        raw = stream.read(CAP+1)
    assert len(raw) <= CAP
    return raw, json.loads(raw)


def main():
    assert __debug__, 'assertions must be enabled'
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    mode.add_argument('--verify-code', action='store_true', help='provenance only; no proof replay')
    args = parser.parse_args()
    if args.write and CERT.exists():
        parser.error('refusing to overwrite certificate')
    assert int(os.environ.get('COXETER_CHECK_OWNER', '-1')) == os.getppid(), 'bounded parent wrapper required'
    signal.alarm(45)
    sources = snapshot()
    saved_raw, saved = (None, None) if args.write else read_certificate()
    if args.verify_code:
        assert saved['sources'] == sources
        fresh(sources)
        assert read_certificate()[0] == saved_raw, 'saved certificate changed'
        print(json.dumps({'status': 'PROVENANCE_ONLY_PASS', 'proof_replay': False,
                          'source_paths': len(sources)}))
        return
    result = check(sources)
    payload = (json.dumps(result, sort_keys=True, indent=2)+'\n').encode()
    assert len(payload) <= CAP, 'certificate exceeds128KiB'
    fresh(sources)
    if args.write:
        with CERT.open('xb') as stream:
            stream.write(payload)
    else:
        assert saved == result, 'saved certificate differs'
        assert read_certificate()[0] == saved_raw, 'saved certificate changed'
    print(json.dumps({'status': 'PASS', 'mode': 'write' if args.write else 'verify',
                      'fixed_graphs': 2, 'vertices_edges_triangles_each': [16, 42, 28],
                      'certificate_bytes': len(payload),
                      'certificate_sha256': hashlib.sha256(payload).hexdigest()}))


if __name__ == '__main__':
    main()
