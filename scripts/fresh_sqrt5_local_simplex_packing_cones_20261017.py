#!/usr/bin/env python3
"""Fourteen small complete local packing cones; exact replay never solves.

Ten forests maximize the prior all-subset mass increase.  Four further forests
are selected by new violated intersecting-orbit clique inequalities.  This is
not a repeat of the global LP and does not inspect other forests' full cones.
"""
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_local_simplex_packing_cones_20261017'
SUBSETS = 'fresh_sqrt5_all_subset_port_packing_20261017'
CLIQUES = 'fresh_sqrt5_intersecting_block_cliques_20261017'
POINTED = 'fresh_sqrt5_pointed_simplex_base_capacity_20261017'
N = 13
FULL = (1 << N) - 1


def dump(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'))


def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()


def read(tag):
    obj = json.loads((ROOT / f'results/{tag}.json').read_text())
    for p, h in obj['sources'].items():
        assert sha(p) == h, p
    return obj


class Packings:
    def __init__(self, columns):
        self.columns = columns
        self.ids = [c['index'] for c in columns]
        self.local = {j: i for i, j in enumerate(self.ids)}
        self.bybit = defaultdict(list)
        self.masks = {}
        for i, c in enumerate(columns):
            masks = c['literal_masks']
            assert masks and len(set(masks)) == len(masks)
            for mask in masks:
                assert 0 < mask <= FULL and mask.bit_count() == c['block_size']
                assert mask not in self.masks
                self.masks[mask] = i
                for p in range(N):
                    if mask >> p & 1:
                        self.bybit[p].append((mask, i))

        @lru_cache(None)
        def count_vectors(mask):
            if mask == 0:
                return frozenset((0,))
            p = (mask & -mask).bit_length() - 1
            out = set(count_vectors(mask ^ (1 << p)))
            for block, i in self.bybit[p]:
                if block & mask == block:
                    # Four bits suffice: there are only thirteen ports.
                    out.update(x + (1 << (4 * i)) for x in count_vectors(mask ^ block))
            return frozenset(out)
        self.dp = count_vectors
        self.vectors = sorted(count_vectors(FULL))
        self.counts = [[(x >> (4 * i)) & 15 for i in range(len(columns))]
                       for x in self.vectors]
        assert all(sum(row[i] * columns[i]['block_size'] for i in range(len(row))) <= N
                   for row in self.counts)

    def witness(self, vector):
        def go(mask, x):
            if mask == 0:
                assert x == 0
                return []
            p = (mask & -mask).bit_length() - 1
            remaining = mask ^ (1 << p)
            if x in self.dp(remaining):
                return go(remaining, x)
            for block, i in self.bybit[p]:
                unit = 1 << (4 * i)
                if block & mask == block and ((x >> (4 * i)) & 15) and x - unit in self.dp(mask ^ block):
                    return [[self.ids[i], block]] + go(mask ^ block, x - unit)
            raise AssertionError((mask, x))
        assert vector in self.dp(FULL)
        return go(FULL, vector)

    def verify_witness(self, blocks):
        used = 0
        row = [0] * len(self.columns)
        for j, mask in blocks:
            assert mask in self.masks and self.ids[self.masks[mask]] == j
            assert used & mask == 0
            used |= mask
            row[self.local[j]] += 1
        encoded = sum(m << (4 * i) for i, m in enumerate(row))
        assert encoded in self.dp(FULL)
        return row


def rational_linear_solution(rows, rhs, guess):
    """Exact Gauss-Jordan, retaining supplied rational free coordinates."""
    n = len(guess)
    a = [[F(v) for v in row] + [F(b)] for row, b in zip(rows, rhs)]
    pivot = []
    r = 0
    for col in range(n):
        pos = next((i for i in range(r, len(a)) if a[i][col]), None)
        if pos is None:
            continue
        a[r], a[pos] = a[pos], a[r]
        scale = a[r][col]
        a[r] = [v / scale for v in a[r]]
        for i in range(len(a)):
            if i != r and a[i][col]:
                scale = a[i][col]
                a[i] = [v - scale * w for v, w in zip(a[i], a[r])]
        pivot.append(col)
        r += 1
    if any(not any(row[:n]) and row[n] for row in a):
        return None
    answer = list(guess)
    free = set(range(n)) - set(pivot)
    for row, col in zip(a, pivot):
        answer[col] = row[n] - sum((row[j] * answer[j] for j in free), F(0))
    assert all(sum((F(v) * z for v, z in zip(row, answer)), F(0)) == F(b)
               for row, b in zip(rows, rhs))
    return answer


def solve_local(pack, weights, scale):
    import numpy as np
    from scipy.optimize import linprog
    mat = np.asarray(pack.counts, dtype=float)
    w = np.asarray([float(a / scale) for a in weights])
    result = linprog(-w, A_ub=mat, b_ub=np.ones(len(mat)), bounds=(0, None),
                     method='highs', options={'time_limit': 20})
    assert result.success, (result.status, result.message)
    charges = [max(F(0), F(float(a)).limit_denominator(10**6)) for a in result.x]
    cost = max(sum((m * a for m, a in zip(row, charges)), F(0)) for row in pack.counts)
    assert cost > 0
    charges = [a / cost for a in charges]
    target = sum((w * a for w, a in zip(weights, charges)), F(0))

    active = [i for i, a in enumerate(result.ineqlin.marginals) if -a > 1e-12]
    guesses = [F(float(-result.ineqlin.marginals[i])).limit_denominator(10**12) * scale for i in active]
    covering = [[pack.counts[i][j] for i in active] for j in range(len(weights))]
    primal = None
    for tolerance in (0, 1e-11, 1e-9, 1e-7):
        tight = [j for j in range(len(weights)) if charges[j] > 0 or (
            tolerance and abs(sum(pack.counts[i][j] * -result.ineqlin.marginals[i] for i in active) - w[j]) < tolerance)]
        rows = [[1] * len(active)] + [covering[j] for j in tight]
        rhs = [target] + [weights[j] for j in tight]
        proposed = rational_linear_solution(rows, rhs, guesses)
        if proposed is None or any(a < 0 for a in proposed):
            continue
        if all(sum((m * a for m, a in zip(row, proposed)), F(0)) >= weights[j]
               for j, row in enumerate(covering)):
            primal = proposed
            break
    assert primal is not None, 'Exact local primal recovery failed; no second solve is permitted.'
    mixture = [dict(mass=str(a), blocks=pack.witness(pack.vectors[i]))
               for i, a in zip(active, primal) if a]
    return dict(charges=[[j, str(a)] for j, a in zip(pack.ids, charges) if a],
                packing_mixture=mixture,
                solver=dict(method='highs', time_limit_seconds=20, success=bool(result.success),
                    status=int(result.status), message=str(result.message),
                    normalized_objective=float(-result.fun), iterations=int(result.nit)),
                exact_local_optimum=str(target))


def verify_local(pack, weights, solution):
    q = {j: F(a) for j, a in solution['charges']}
    assert len(q) == len(solution['charges']) and set(q) <= set(pack.ids)
    assert all(a > 0 for a in q.values())
    charges = [q.get(j, F(0)) for j in pack.ids]
    maxcost = max(sum((m * a for m, a in zip(row, charges)), F(0)) for row in pack.counts)
    assert maxcost == 1
    value = sum((w * a for w, a in zip(weights, charges)), F(0))
    cover = [F(0)] * len(weights)
    mass = F(0)
    for item in solution['packing_mixture']:
        a = F(item['mass'])
        assert a > 0
        counts = pack.verify_witness(item['blocks'])
        for j, m in enumerate(counts):
            cover[j] += a * m
        mass += a
    assert all(a >= w for a, w in zip(cover, weights))
    assert mass == value == F(solution['exact_local_optimum'])
    assert all((a - w) * charge == 0 for a, w, charge in zip(cover, weights, charges))
    return value, [[j, str(a - w)] for j, a, w in zip(pack.ids, cover, weights) if a != w]


def run(saved=None):
    data = read(SUBSETS)
    clique = read(CLIQUES)
    pointed = read(POINTED)
    assert data['face_multipliers'] == clique['face_multipliers']
    assert next(t for t in pointed['tests'] if t['name'] == 'full_optimized')['repaired_dual']['face_multipliers'] == data['face_multipliers']
    alpha = {i: F(a) for i, a in data['face_multipliers']}
    subset_mu = {v: F(a) for v, a in data['vertex_multipliers']}
    mu = {v: F(a) for v, a in clique['vertex_multipliers']}
    assert all(mu[v] >= a for v, a in subset_mu.items())
    weights = defaultdict(F)
    for f in data['positive_simplex_features']:
        for j, m in f['block_orbits']:
            weights[j] += alpha[f['face']] * m
    increases = {}
    forests = {}
    for r in data['subset_maxima']:
        v = r['vertex']
        forests[v] = r['forest']
        inc = F(r['required_vertex_mass']) - F(r['old_vertex_mass'])
        increases[v] = max(increases.get(v, inc), inc)
    top = sorted(increases, key=lambda v: (-increases[v], v))[:10]
    extra = sorted({r['vertex'] for r in clique['critical_columns']})
    selected = top + [v for v in extra if v not in top]
    assert len(top) == 10 and len(extra) == 4 and len(selected) == 14
    literals = defaultdict(list)
    for c in data['literal_block_orbits']:
        literals[c['vertex']].append(c)
    for cs in literals.values():
        cs.sort(key=lambda c: c['index'])
    saved_rows = {} if saved is None else {r['vertex']: r for r in saved['local_cones']}
    records, columns = [], []
    repaired = dict(mu)
    for v in selected:
        pack = Packings(literals[v])
        w = [weights[j] for j in pack.ids]
        assert all(a > 0 for a in w)
        if saved is None:
            solution = solve_local(pack, w, subset_mu[v])
        else:
            solution = saved_rows[v]['certificate']
        value, slack = verify_local(pack, w, solution)
        repaired[v] = max(repaired[v], value)
        q = {j: F(a) for j, a in solution['charges']}
        facecoeff = []
        for f in data['positive_simplex_features']:
            a = sum((m * q.get(j, F(0)) for j, m in f['block_orbits']), F(0))
            if a:
                facecoeff.append([f['face'], str(a)])
        assert sum((alpha[i] * F(a) for i, a in facecoeff), F(0)) == value
        columns.append(dict(index=len(columns), vertex=v, vertex_cost=1,
                            positive_simplex_coefficients=facecoeff, exact_load=str(value)))
        records.append(dict(vertex=v, forest=forests[v],
            selection='TOP10_ALL_SUBSET_MASS_INCREASE' if v in top else 'NEW_VIOLATED_INTERSECTING_ORBIT_CLIQUE',
            counts=dict(block_orbits=len(pack.ids), literal_blocks=len(pack.masks),
                complete_packing_count_vectors=len(pack.vectors), cached_available_masks=pack.dp.cache_info().currsize,
                positive_mixture_packings=len(solution['packing_mixture'])),
            block_orbit_marginals=[[j, str(a)] for j, a in zip(pack.ids, w)],
            literal_block_orbits=literals[v],
            count_vector_digest=hashlib.sha256(dump(pack.vectors).encode()).hexdigest(),
            certificate=solution, primal_cover_slacks=slack,
            all_subset_vertex_mass=str(subset_mu[v]), clique_repaired_vertex_mass=str(mu[v]),
            local_optimum=str(value), additional_mass_over_subsets=str(max(F(0), value - subset_mu[v])),
            additional_mass_over_cliques=str(max(F(0), value - mu[v]))))
        print('EXACT_LOCAL_CONE', v, len(pack.vectors), str(value),
              'over_subsets', str(max(F(0), value - subset_mu[v])),
              'over_cliques', str(max(F(0), value - mu[v])), flush=True)
    numerator = F(data['objective_numerator'])
    normalizer = sum(repaired.values(), F(0))
    paths = [f'scripts/{TAG}.py', f'scripts/{SUBSETS}.py', f'results/{SUBSETS}.json',
             f'scripts/{CLIQUES}.py', f'results/{CLIQUES}.json',
             f'scripts/{POINTED}.py', f'results/{POINTED}.json',
             'results/fresh_simplex_disjoint_ports_20261017.md']
    return dict(status='PASS_EXACT_FOURTEEN_LOCAL_PACKING_CONE_CERTIFICATES',
        sources={p: sha(p) for p in paths}, dimension=13, target=12,
        generation_local_solver_invocations=14, saved_replay_solver_invocations=0,
        global_solver_invocations=0, generation_time_limit_per_local_solver_seconds=20,
        counts=dict(selected_forests=len(selected), top_increase_forests=len(top), new_clique_forests=len(extra),
            complete_packing_count_vectors=sum(r['counts']['complete_packing_count_vectors'] for r in records),
            newly_increased_vertex_masses=sum(F(r['additional_mass_over_cliques']) > 0 for r in records)),
        top_increase_forests=top, new_clique_forests=extra, local_cones=records,
        exact_charge_columns=columns, objective_numerator=str(numerator), normalizer=str(normalizer),
        bound=str(numerator / normalizer), excludes_target12=numerator / normalizer > 12,
        face_multipliers=data['face_multipliers'], vertex_multipliers=[[v, str(a)] for v, a in sorted(repaired.items())],
        theorem='The selected maximal-simplex blocks are pairwise disjoint at each finite vertex. For arbitrary nonnegative orbit-invariant block charges, the vertex cost is their maximum sum over disjoint literal packings. The DP enumerates all achievable orbit-count vectors by either skipping the least available port or selecting its unique block in the packing. The dual charge has exact maximum cost1. A nonnegative mixture of explicit disjoint packings covers every required orbit marginal with the same total mass, proving the exact local optimum. Averaging each packing over the finite-forest automorphism group gives uniform marginals on each complete orbit, so orbit-count mixtures suffice.',
        scope='Only14 justified local cones of the existing positive simplex-block domain are solved; the remaining726 forests are not solved. Unknown blocks can be assigned zero charge and omitted from packing, so the stated charges are valid restricted indicators. All quad coefficients are zero. Every old face mass, actual-provenance split, inherited signed transfer and pointed capacity remains unchanged; vertex masses only increase. No global optimization, complete-face coverage or dimension bound is claimed.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--verify', action='store_true')
    mode.add_argument('--refresh', action='store_true', help='reconstruct from saved exact certificates, without solving')
    a = p.parse_args()
    out = ROOT / f'results/{TAG}.json'
    old = json.loads(out.read_text()) if a.verify or a.refresh else None
    result = json.loads(json.dumps(run(old)))
    if a.verify:
        assert result == old
    else:
        out.write_text(dump(result) + '\n')
    print(result['status'], result['counts'], result['bound'], flush=True)
