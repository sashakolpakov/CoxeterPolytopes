#!/usr/bin/env python3
"""Exact local maximal-simplex packing test of one frozen rank-13 dual.

No optimizer or new diagram census.  The graph helper is frozen and source
bound; all incidence masses and the 22 exported successor columns are rebuilt
from literal stars with Fraction arithmetic.
"""
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path

import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_simplex_port_incidence_20261017'
LP = 'fresh_sqrt5_rank13_full_flexible_simplex_lp_20261017'
STARS = 'fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017'
FLAGS = 'fresh_sqrt5_rank13_triangle_simplex_flags_20261017'


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def read(tag):
    d = json.loads((ROOT / 'results' / (tag + '.json')).read_text())
    for p, h in d.get('sources', {}).items():
        assert sha(p) == h, (tag, p)
    return d


def freeze(x):
    return tuple(map(freeze, x)) if isinstance(x, list) else x


def code(pieces, marks=()):
    return tuple(sorted(C.code(p, marks) for p in pieces))


def forest_graph(forest, models):
    g = {}
    for name in forest:
        h = models[name]
        offset = len(g)
        for i, row in h.items():
            g[offset + i] = {offset + j: m for j, m in row.items()}
    assert len(g) == 13
    return g


def star_corners(star, finite):
    g = C.graph(range(14), {(i, j): m for i, j, m in star['edges']})
    S = {r['omitted'] for r in star['spherical_deletions']}
    assert len(S) - 1 == star['simplex_dimension']
    out = []
    for row in star['spherical_deletions']:
        o = row['omitted']
        h = {i: {j: m for j, m in rr.items() if j != o}
             for i, rr in g.items() if i != o}
        pieces = C.pieces(h)
        forest = tuple(sorted(finite[C.code(p)] for p in pieces))
        assert forest == tuple(row['forest'])
        ports = tuple(sorted(S - {o}))
        out.append((o, forest, ports, pieces, code(pieces, ports)))
    return out


def run():
    d = read(LP)
    stars = read(STARS)['stars']
    dual = next(c for c in d['certificates']
                if c['kind'] == 'EXACT_RESTRICTED_MODEL_DUAL')
    alpha = {i: F(q) for i, q in dual['face_multipliers']}
    mu = {i: F(q) for i, q in dual['vertex_multipliers']}
    assert len(alpha) == 346 and len(mu) == 203
    models = dict(C.models())
    finite = {C.code(g): name for name, g in models.items() if len(g) <= 13}
    assert len(finite) == sum(len(g) <= 13 for g in models.values())
    forest_vertex = {tuple(v['forest']): v['index'] for v in d['vertices']}
    pair_index = {(tuple(v['forest']), freeze(v['marked_forest_code'])): v['index']
                  for v in d['variables']}
    mass = defaultdict(lambda: defaultdict(F))
    origins = defaultdict(lambda: defaultdict(list))
    loads = defaultdict(F)
    dimensions = Counter()

    # Independently rebuild every positive-dual row's entire pair incidence.
    for fi, a in sorted(alpha.items()):
        face = d['faces'][fi]
        assert a > 0 and len(face['provenance']) == 1
        si = face['provenance'][0]['star']
        star = stars[si]
        assert star['index'] == si and not star['fixed_discriminant']
        corners = star_corners(star, finite)
        assert face['required'] == 2 * comb(len(corners), 3)
        dimensions[len(corners) - 1] += 1
        coefficients = Counter()
        for omitted, forest, ports, pieces, bc in corners:
            v = forest_vertex[forest]
            assert v in mu
            mass[v][bc] += a
            origins[v][bc].append(dict(face=fi, star=si, omitted=omitted,
                                       ports=list(ports), multiplier=str(a)))
            for pair in combinations(ports, 2):
                coefficients[pair_index[forest, code(pieces, pair)]] += 1
        assert sorted(coefficients.items()) == [tuple(z) for z in face['coefficients']]
        for j, n in coefficients.items():
            loads[j] += a * n
    assert set(mass) == set(mu)
    for var in d['variables']:
        assert loads[var['index']] <= mu.get(var['vertex'], F()) * var['multiplicity']
    numerator = sum((a * d['faces'][i]['required'] for i, a in alpha.items()), F())
    normalizer = sum(mu.values(), F())
    assert numerator == F(dual['objective_numerator'])
    assert normalizer == F(dual['normalizer'])
    assert numerator / normalizer == F(dual['bound']) > 12

    summaries = []
    violations = []
    all_orbits = {}
    tests = 0
    literal_pair_checks = 0
    for serial, v in enumerate(sorted(mu)):
        forest = d['vertices'][v]['forest']
        g = forest_graph(forest, models)
        pieces = C.pieces(g)
        orbits = defaultdict(list)
        for size in (2, 3, 4):
            for B in combinations(range(13), size):
                bc = code(pieces, B)
                if bc in mass[v]:
                    orbits[bc].append(B)
        assert set(orbits) == set(mass[v])
        literal = {B: mass[v][bc] / len(bs)
                   for bc, bs in orbits.items() for B in bs}
        pair_load = defaultdict(F)
        for B, value in literal.items():
            for p in combinations(B, 2):
                pair_load[p] += value
        for p in combinations(range(13), 2):
            j = pair_index[tuple(forest), code(pieces, p)]
            assert pair_load[p] == loads[j] / d['variables'][j]['multiplicity']
            assert pair_load[p] <= mu[v]
            literal_pair_checks += 1
        # Needed for the fixed-alpha minimum-mass consequence: the old beta
        # is the least mass allowed by its pair capacities, not a loose bound.
        assert max(pair_load.values()) == mu[v]

        best = None
        for usize, k in ((4, 3), (5, 4), (6, 4)):
            assert 2 * k - usize >= 2
            for U in combinations(range(13), usize):
                value = sum((literal.get(B, F()) for size in range(k, 5)
                             for B in combinations(U, size)), F())
                tests += 1
                if best is None or value > best[0]:
                    best = value, usize, k, U
        value, usize, k, U = best
        summary = dict(vertex=v, forest=forest, vertex_mass=str(mu[v]),
                       block_orbits=len(orbits), literal_blocks=len(literal),
                       subset=list(U), minimum_block_size=k,
                       packing_mass=str(value), ratio=str(value / mu[v]),
                       pair_capacity_maximum=str(max(pair_load.values())),
                       violates=value > mu[v])
        summaries.append(summary)
        if value > mu[v]:
            Ucode = code(pieces, U)
            Uorbit = [V for V in combinations(range(13), usize)
                      if code(pieces, V) == Ucode]
            terms = []
            for bc, bs in sorted(orbits.items()):
                hits = [B for B in bs if len(B) >= k and set(B) <= set(U)]
                if not hits:
                    continue
                coefficient = F(len(hits), len(bs))
                # Double-count incidences to verify the orbit-average charge.
                for B in bs:
                    assert F(sum(set(B) <= set(V) for V in Uorbit), len(Uorbit)) == coefficient
                terms.append(dict(literal_blocks=[list(B) for B in bs],
                                  selected_blocks=[list(B) for B in hits],
                                  block_orbit_mass=str(mass[v][bc]),
                                  averaged_coefficient=str(coefficient),
                                  contribution=str(coefficient * mass[v][bc]),
                                  origins=origins[v][bc]))
            assert sum((F(t['contribution']) for t in terms), F()) == value
            violations.append(dict(**summary, subset_orbit=[list(V) for V in Uorbit],
                                   subset_orbit_size=len(Uorbit), terms=terms))
            all_orbits[v] = (pieces, U, k, Uorbit)
        if serial % 40 == 0:
            print('literal packing forests', serial, '/', len(mu), flush=True)

    assert len(violations) == 22
    assert tests == len(mu) * (comb(13, 4) + comb(13, 5) + comb(13, 6))
    by_v = {r['vertex']: r for r in violations}
    assert F(by_v[390]['packing_mass']) == F(107171663, 10**10)
    assert F(by_v[390]['ratio']) == F(107171663, 58638750)
    assert F(by_v[951]['packing_mass']) == F(300792247, 10**10)
    assert F(by_v[951]['ratio']) == F(300792247, 283103766)

    # Export valid extra primal columns over EVERY old face, not just dual rows.
    # One column = automorphism-average of a clique indicator.  Its entire
    # vertex cost is1; its coefficient on a simplex is the sum over corners.
    column_coeff = {v: defaultdict(F) for v in all_orbits}
    block_coefficient_cache = {}
    for fi, face in enumerate(d['faces']):
        assert len(face['provenance']) == 1
        si = face['provenance'][0]['star']
        star = stars[si]
        assert star['index'] == si
        corners = star_corners(star, finite)
        assert face['required'] == 2 * comb(len(corners), 3)
        for _, forest, ports, _, bc in corners:
            v = forest_vertex[forest]
            if v not in all_orbits:
                continue
            key = v, bc
            if key not in block_coefficient_cache:
                pieces, U, k, Uorbit = all_orbits[v]
                size = len(ports)
                if size < k or size > len(U):
                    coefficient = F()
                else:
                    orbit = [B for B in combinations(range(13), size)
                             if code(pieces, B) == bc]
                    assert orbit
                    coefficient = F(sum(set(B) <= set(U) for B in orbit), len(orbit))
                    assert F(sum(set(orbit[0]) <= set(V) for V in Uorbit), len(Uorbit)) == coefficient
                block_coefficient_cache[key] = coefficient
            coefficient = block_coefficient_cache[key]
            if coefficient:
                column_coeff[v][fi] += coefficient

    columns = []
    for v in sorted(all_orbits):
        record = by_v[v]
        coefficients = column_coeff[v]
        load = sum((alpha.get(fi, F()) * c for fi, c in coefficients.items()), F())
        assert load == F(record['packing_mass']) > mu[v]
        columns.append(dict(index=len(columns), vertex=v, forest=record['forest'],
                            vertex_cost='1', subset=record['subset'],
                            minimum_block_size=record['minimum_block_size'],
                            subset_orbit=record['subset_orbit'],
                            face_coefficients=[[i, str(q)] for i, q in sorted(coefficients.items())],
                            fixed_dual_load=str(load), fixed_dual_capacity=str(mu[v])))
    extra = sum((F(r['packing_mass']) - mu[r['vertex']] for r in violations), F())
    assert extra == F(7706660539, 300000000000)
    repaired_upper = numerator / (normalizer + extra)
    assert repaired_upper == F(2104176880435, 179495552776) < 12
    paths = [f'scripts/{TAG}.py', f'scripts/{C.TAG}.py', f'scripts/{LP}.py',
             f'results/{LP}.json', f'results/{STARS}.json',
             f'scripts/{FLAGS}.py', f'results/{FLAGS}.json', f'results/{FLAGS}.md']
    return dict(
        status='PASS_EXACT_MAXIMAL_SIMPLEX_PORT_PACKING_VIOLATION',
        sources={p: sha(p) for p in paths},
        theorem=dict(
            ambient='Simple compact acute-angled hyperbolic polytope; in particular a compact Coxeter polytope.',
            ports='A d-simplex face through a vertex selects its d incident edge ports.',
            intersection='Distinct maximal simplex faces have port-set intersection of size at most1, by unique maximal-simplex extension of each triangular two-face.',
            clique='For U of size u and k with2k-u>=2, at most one maximal-simplex port block B with |B|>=k and B subset U occurs.',
            orbit_average='Average the clique indicator over the full labelled-diagram automorphism orbit of U; its total at every vertex remains at most1.',
            completeness='The theorem concerns every maximal simplex, without an observed dimension or label cutoff. The diagnostic tests only the frozen dual and three small (u,k) clique types.'),
        arithmetic='Exact Fraction arithmetic; canonical marked-tree graphs from a bound frozen helper. No optimizer and no new spherical or arithmetic-star census.',
        counts=dict(dual_faces=len(alpha), positive_vertex_forests=len(mu),
                    dual_simplex_dimensions={str(k): v for k, v in sorted(dimensions.items())},
                    literal_pair_checks=literal_pair_checks,
                    tested_cliques=tests, violated_forests=len(violations),
                    successor_columns=len(columns),
                    successor_faces=len(d['faces']),
                    successor_nonzero_coefficients=sum(len(c['face_coefficients']) for c in columns)),
        old_dual=dict(objective_numerator=str(numerator), normalizer=str(normalizer),
                      bound=str(numerator / normalizer), all_old_pair_capacities_pass=True),
        forest_tests=summaries, violations=violations,
        fixed_face_mass_consequence=dict(
            minimum_additional_vertex_mass=str(extra),
            minimum_total_vertex_mass=str(normalizer + extra),
            demand_per_vertex_at_most=str(repaired_upper),
            scope='Necessary only for realization with these exact fixed face masses. Increasing vertex masses is not asserted sufficient for a realizable distribution.'),
        successor_columns=columns,
        successor_scope='Every old face and budget remains. These22 valid clique charges strengthen the partial model; no augmented optimization, feasible rule, complete face coverage or dimension bound is claimed.',
        diagnosis='The frozen dual certifies an obstruction to its pair-charge model, but its face/vertex masses violate actual maximal-simplex incidence. The contradiction is to that hypothetical global incidence distribution, not to the arithmetic validity of any individual local star.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    result = run()
    path = ROOT / 'results' / (TAG + '.json')
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['counts'], result['fixed_face_mass_consequence'], flush=True)
