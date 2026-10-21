#!/usr/bin/env python3
"""Separate sparse-radical review of the two finite-edge caps.

Imports only an earlier independent arithmetic implementation, not the cap
producer or its arithmetic. The finite-support census uses exhaustive subsets
and FLINT principal minors, independently of the producer's clique/LDL code.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from flint import fmpq, fmpq_mat
import fresh_h20_e7_replacement_escape_20261014_review as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_finite_edgewalk_20261014'
PRIMARY = f'results/{TAG}_exact.json'
FRONTIER = f'results/{TAG}_frontier.json'
BRIDGE = 'results/fresh_h20_e7_causal_residuals_20261014_bridge.json'
FACE = 'results/fresh_h20_e7_intrinsic_caps_20261014.json'
ARITH = 'scripts/fresh_h20_e7_replacement_escape_20261014_review.py'
RAD = [1, 2, 79, 158, 11, 22, 869, 1738, 3, 6, 237, 474, 33, 66, 2607, 5214]

def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()

def read(p):
    d = json.loads((ROOT / p).read_text())
    for s, h in d.get('sources', {}).items():
        assert sha(s) == h, s
    return d

def run():
    d, old, face = map(read, [PRIMARY, BRIDGE, FACE])
    assert d['normals'][:88] == old['normals']
    assert d['form'] == old['form'] and d['field'] == old['field']
    def base(a):
        assert len(a) in [8, 16]
        return {r: F(x) for r, x in zip(RAD, a) if F(x)}
    D = base(d['field']['discriminant'])
    assert R.sign(D) > 0
    zero = ({}, {})
    def num(q): return (R.c(q), {})
    def add(a, b): return (R.add(a[0], b[0]), R.add(a[1], b[1]))
    def scale(a, q): return (R.scale(a[0], q), R.scale(a[1], q))
    def mul(a, b):
        return (R.add(R.mul(a[0], b[0]), R.mul(D, R.mul(a[1], b[1]))),
                R.add(R.mul(a[0], b[1]), R.mul(a[1], b[0])))
    def sign(a):
        sa, sb = R.sign(a[0]), R.sign(a[1])
        if not sb: return sa
        if not sa or sa == sb: return sb
        return sa * R.sign(R.sub(R.mul(a[0], a[0]), R.mul(D, R.mul(a[1], a[1]))))
    def total(xs):
        z = zero
        for x in xs: z = add(z, x)
        return z
    def parse(a):
        assert len(a) == 2
        return (base(a[0]), base(a[1]))
    def vec(a): return tuple(map(parse, a))
    def bvec(a): return tuple((base(x), {}) for x in a)
    Q = [[F(x) for x in row] for row in d['form']]
    assert len(Q) == 21 and all(len(row) == 21 for row in Q)
    def cov(v): return tuple(total(scale(v[j], q) for j, q in enumerate(row) if q) for row in Q)
    def rawdot(v, w): return total(mul(a, b) for a, b in zip(v, w))
    def dot(v, w): return rawdot(v, cov(w))
    def lin(vs, cs): return tuple(total(mul(v[k], c) for v, c in zip(vs, cs)) for k in range(21))
    def rat(a):
        assert not a[1]
        return R.rational(a[0])
    N = list(map(vec, d['normals']))
    CN = list(map(cov, N))
    norms = list(map(F, d['norms']))
    assert len(N) == len(norms) == 90
    H = vec(d['strict_center']); CH = cov(H)
    assert H == vec(old['strict_center']) and sign(rawdot(H, CH)) < 0
    assert all(dot(n, n) == num(q) for n, q in zip(N, norms))
    assert all(sign(rawdot(H, n)) < 0 for n in CN)
    pair_counts = Counter()
    pair_records = {}
    for i in [88, 89]:
        for j in range(i):
            g = rawdot(N[i], CN[j]); p = norms[i] * norms[j]
            assert sign(g) <= 0
            sq = mul(g, g)
            if g == zero: kind = 'right'
            elif sq == num(p / 4): kind = 'pi3'
            elif sq == num(p / 2): kind = 'pi4'
            elif sq == num(3 * p / 4): kind = 'pi6'
            elif sq == num(p): kind = 'parallel'
            else:
                assert sign(add(sq, num(-p))) > 0
                kind = 'dotted'
            pair_counts[kind] += 1
            pair_records[i, j] = (g, kind)
    for rec in d['new_pairs']:
        assert pair_records[tuple(rec['walls'])] == (parse(rec['pairing']), rec['relation'])
    assert len(d['new_pairs']) == len(pair_records) == 177
    assert dict(pair_counts) == d['new_pair_counts']
    assert rawdot(N[88], CN[89]) == zero
    # Independent reconstruction of the shared-support projection and E7 dual.
    S = [12, 13, 14, 15, 18, 31, 65]
    I = R.inverse([[rat(rawdot(N[i], CN[j])) for j in S] for i in S])
    omega = lin([N[i] for i in S], [num(q) for q in I[S.index(13)]])
    factors = [{2: F(1)} if i in face['face_facet_walls'][11:13] else R.c(1) for i in range(90)]
    for index, (vi, contacts, boost) in enumerate([(4, [0, 14], F(7, 4)), (10, [10, 14], F(3, 2))]):
        rec = d['recipes'][index]
        support = [i for i in face['face_finite_vertices'][vi]['ambient_incident'] if i != 13]
        rows = [lin([N[i]], [(factors[i], {})]) for i in support]
        G = [[rat(dot(a, b)) for b in rows] for a in rows]
        GI = R.inverse(G)
        rhs = [(R.scale(factors[i], -int(i in contacts)), {}) for i in support]
        coeff = [total(scale(rhs[j], GI[k][j]) for j in range(19)) for k in range(19)]
        p = lin(rows, coeff); V = bvec(face['face_finite_vertices'][vi]['ambient_point'])
        cap = lin([p, omega, V], [num(1), num(F(-2, 3)), num(boost)])
        assert cap == N[88 + index] == bvec(rec['normal'])
        assert p == bvec(rec['projection'])
        assert dot(p, p) == num([F(47, 30), F(22, 15)][index])
    # Facet projection is timelike. Every other inequality is strict by
    # H*n_j<0, nonpositive n_i*n_j, and positive projection coefficient.
    for i, n in enumerate(N):
        c = scale(rawdot(H, CN[i]), -1 / norms[i])
        p = lin([H, n], [num(1), c])
        assert p == vec(d['facet_points'][i])
        assert sign(c) > 0 and rawdot(p, CN[i]) == zero
        assert sign(dot(p, p)) < 0 and sign(rawdot(p, CH)) < 0
    old_vertices = d['preserved_old_vertices']
    new_vertices = d['new_cap_finite_vertices']
    all_vertices = old_vertices + new_vertices
    assert len(old_vertices) == 44 and len(new_vertices) == 40
    supports = []
    for rec in all_vertices:
        p = bvec(rec['point']); active = sorted(rec['incident'])
        assert len(active) == len(set(active)) == 20
        assert sign(dot(p, p)) < 0 and sign(rawdot(p, CH)) < 0
        signs = [sign(rawdot(p, n)) for n in CN]
        assert all(s == 0 if i in active else s < 0 for i, s in enumerate(signs))
        supports.append(tuple(active))
    assert len(set(supports)) == 84
    for k, rec in enumerate(old_vertices):
        assert rec['point'] == face['face_finite_vertices'][k]['ambient_point']
    print('review: all90 walls and84 finite vertices pass', flush=True)
    # Re-enumerate every possible finite vertex incident to a new wall.
    # A finite Coxeter vertex is simple; its20 facets have pairwise finite
    # intersections. Full subsets rather than the producer's clique recursion.
    neighbours = {}
    for i in [88, 89]:
        neighbours[i] = [j for j in range(90) if j != i and
            sign(add(mul(rawdot(N[i], CN[j]), rawdot(N[i], CN[j])), num(-norms[i] * norms[j]))) < 0]
        assert neighbours[i] == d['new_cap_finite_neighbor_sets'][str(i)]
    universe = sorted(set(neighbours[88] + neighbours[89] + [88, 89]))
    gram = {}
    finite = {}
    for i in universe:
        for j in universe:
            g = rawdot(N[i], CN[j])
            gram[i, j] = rat(mul(g, (R.mul(factors[i], factors[j]), {})))
            finite[i, j] = i != j and sign(add(mul(g, g), num(-norms[i] * norms[j]))) < 0
    candidates = set()
    for i in [88, 89]:
        for sub in combinations(neighbours[i], 19):
            if all(finite[a, b] for a, b in combinations(sub, 2)):
                candidates.add(tuple(sorted(sub + (i,))))
    positives, negatives = set(), set()
    for ids in sorted(candidates):
        # Sylvester's criterion, using an unrelated rational matrix library.
        pd = all(fmpq_mat([[fmpq(str(gram[i, j])) for j in ids[:k]] for i in ids[:k]]).det() > 0
                 for k in range(1, 21))
        (positives if pd else negatives).add(ids)
    assert positives == {tuple(sorted(v['incident'])) for v in new_vertices}
    assert negatives == set(map(tuple, d['new_cap_nonPD_cliques']))
    assert len(candidates) == d['new_cap_candidate_cliques'] == 3080
    assert len(positives) == 40 and len(negatives) == 3040
    assert sum(88 in v and 89 in v for v in positives) == 9
    endpoints = d['selected_edge_endpoints']
    for rec in endpoints:
        vi, released, wall = rec['old_vertex'], rec['released'], rec['new_wall']
        V = bvec(old_vertices[vi]['point']); p = bvec(rec['point']); t = F(rec['parameter'])
        w = lin([N[i] for i in S], [num(q) for q in I[S.index(released)]])
        assert t > 0 and p == lin([V, w], [num(1), num(-t)])
        assert dot(p, p) == (base(rec['norm']), {}) and sign(dot(p, p)) < 0
        assert sign(rawdot(p, CH)) < 0 and rawdot(p, CN[released]) == num(-t)
        active = tuple(sorted(rec['incident']))
        assert active in positives and wall in active
        assert set(rec['edge_support']) == set(old_vertices[vi]['incident']) - {released} == set(active) - {wall}
        assert all(sign(rawdot(p, n)) == 0 if i in active else sign(rawdot(p, n)) < 0 for i, n in enumerate(CN))
    assert len(endpoints) == 31 and len({tuple(sorted(v['incident'])) for v in endpoints}) == 31
    assert Counter(v['new_wall'] for v in endpoints) == {88: 19, 89: 12}
    edges = []
    for i, j in combinations(range(84), 2):
        common = sorted(set(supports[i]) & set(supports[j]))
        if j >= 44 and len(common) == 19: edges.append(dict(vertices=[i, j], support=common))
    assert edges == d['finite_edges_with_new_endpoint'] and len(edges) == 145
    for rec in d['preserved_cusps']:
        p = bvec(rec['point'])
        assert dot(p, p) == zero and sign(rawdot(p, CH)) < 0
        assert all(sign(rawdot(p, n)) <= 0 for n in CN)
        assert all(sign(rawdot(p, CN[i])) < 0 for i in [87, 88, 89])
    E = vec(d['strict_escape']['point'])
    assert E == vec(old['strict_escape']['point'])
    q = dot(E, E)
    assert q == parse(d['strict_escape']['norm']) and sign(q) > 0 and sign(rawdot(E, CH)) < 0
    assert [rawdot(E, n) for n in CN] == list(map(parse, d['strict_escape']['all90_pairings']))
    assert all(sign(rawdot(E, n)) < 0 for n in CN)
    frontier = read(FRONTIER)
    directions = []
    for i in range(44, 84):
        for j in supports[i]:
            support = set(supports[i]) - {j}
            ends = [k for k, active in enumerate(supports) if k != i and support <= set(active)]
            assert len(ends) <= 1
            if not support & {88, 89}: assert ends and ends[0] < 44
            directions.append(dict(vertex=i, released=j, support=sorted(support),
                                   other_finite_vertex=ends[0] if ends else None))
    assert directions == frontier['new_corner_directions']
    assert len(directions) == 800 and sum(x['other_finite_vertex'] is not None for x in directions) == 259
    cusp = frontier['new_parallel_cusp']
    nu = bvec(cusp['point']); z = bvec(cusp['lineality_vector'])
    assert nu == lin([N[9], N[88]], [num(1), num(1)])
    assert dot(nu, nu) == zero and sign(rawdot(nu, CH)) < 0
    pairings = [rawdot(nu, n) for n in CN]
    assert all(sign(g) <= 0 for g in pairings)
    active = [i for i, g in enumerate(pairings) if g == zero]
    assert active == cusp['incident'] and len(active) == 20
    assert dot(z, nu) == zero and dot(z, z) == num(F(98, 27))
    assert all(rawdot(z, CN[i]) == zero for i in active)
    cg = [[rat(mul(rawdot(N[i], CN[j]), (R.mul(factors[i], factors[j]), {}))) for j in active] for i in active]
    assert cg == [[F(x) for x in row] for row in cusp['scaled_incident_gram']]
    assert fmpq_mat([[fmpq(str(x)) for x in row] for row in cg]).rank() == 18 == cusp['quotient_normal_rank']
    covered = []
    for component in cusp['components']:
        ids = component['walls']; covered.extend(ids)
        M = [[cg[active.index(i)][active.index(j)] for j in ids] for i in ids]
        assert M == [[F(x) for x in row] for row in component['gram']]
        assert all(not cg[active.index(i)][active.index(j)] for i in ids for j in active if j not in ids)
        assert fmpq_mat([[fmpq(str(x)) for x in row] for row in M]).rank() == component['rank']
        for kernel in component['kernel']:
            v = list(map(F, kernel))
            assert all(x > 0 for x in v) and all(sum(a * b for a, b in zip(row, v)) == 0 for row in M)
    assert sorted(covered) == active
    paths = [PRIMARY, FRONTIER, BRIDGE, FACE, ARITH, f'scripts/{TAG}_review.py']
    return dict(status='PASS_SEPARATE_SPARSE_RADICAL_90_WALL_REVIEW', sources={p: sha(p) for p in paths},
        new_pair_counts=dict(pair_counts), actual_facets=90, distinct_finite_vertices=84,
        complete_new_cap_finite_supports=40, candidate_supports=3080, rejected_supports=3040,
        adjacent_endpoints=31, finite_edges_with_new_endpoint=145, preserved_bounded_cusps=2,
        new_corner_germs=800, germs_with_second_finite_vertex=259, germs_without_second_finite_vertex=541,
        new_parallel_cusp=dict(incident_facets=20, quotient_normal_rank=18, lineality_norm='98/27'),
        surviving_escape_norm=str(rat(q)),
        scope='Independent sparse squarefree arithmetic with exact outward sign bounds; rebuild both cap recipes,177 pairs,90 facet points,84 distinct finite points,31 forward adjacent endpoints,145 finite edges and the strict escape. Exhaustive19-subsets and FLINT principal minors independently reproduce the complete new-cap finite-vertex census and800 finite-endpoint statuses. Independently check the new null cusp,all90 inequalities,rank18 and positive-norm tangent lineality. The frozen88 Coxeter proof and two old cusp-link completeness proofs are inherited. No full ambient old-wall vertex census or finite-volume assertion.')

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); a = ap.parse_args()
    d = run(); p = ROOT / f'results/{TAG}_review.json'
    if a.verify: assert d == json.loads(p.read_text())
    else: p.write_text(json.dumps(d, indent=2) + '\n')
    print(d['status'], flush=True)
