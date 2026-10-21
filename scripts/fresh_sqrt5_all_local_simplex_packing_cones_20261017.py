#!/usr/bin/env python3
"""All positive-simplex vertex packing cones, with exact saved certificates.

No global LP. Fourteen frozen local certificates are reused. Fresh replay
performs no numerical optimization. Explicit resource stops remain open.
"""
import argparse
import hashlib
import json
import time
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

import fresh_sqrt5_local_simplex_packing_cones_20261017 as H

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_all_local_simplex_packing_cones_20261017'
PROGRESS = Path('/private/tmp') / (TAG + '_progress.jsonl')
ENTRY_LIMIT = 1000000
ENUM_SECONDS = 30
SOLVE_SECONDS = 20


class ResourceStop(Exception):
    pass


class Packings(H.Packings):
    def __init__(self, columns):
        self.columns = columns
        self.ids = [c['index'] for c in columns]
        self.local = {j: i for i, j in enumerate(self.ids)}
        self.bybit = defaultdict(list)
        self.masks = {}
        self.entries = 0
        start = time.monotonic()
        for i, c in enumerate(columns):
            values = c['literal_masks']
            assert values and len(set(values)) == len(values)
            for mask in values:
                assert 0 < mask <= H.FULL and mask.bit_count() == c['block_size']
                assert mask not in self.masks
                self.masks[mask] = i
                for p in range(H.N):
                    if mask >> p & 1:
                        self.bybit[p].append((mask, i))

        def guard(size=0):
            if self.entries + size > ENTRY_LIMIT:
                raise ResourceStop('COUNT_VECTOR_ENTRY_LIMIT')
            if time.monotonic() - start > ENUM_SECONDS:
                raise ResourceStop('ENUMERATION_TIME_LIMIT')

        @lru_cache(None)
        def count_vectors(mask):
            guard()
            if mask == 0:
                self.entries += 1
                return frozenset((0,))
            p = (mask & -mask).bit_length() - 1
            out = set(count_vectors(mask ^ (1 << p)))
            for block, i in self.bybit[p]:
                if block & mask == block:
                    out.update(x + (1 << (4 * i)) for x in count_vectors(mask ^ block))
                    guard(len(out))
            self.entries += len(out)
            guard()
            return frozenset(out)
        self.dp = count_vectors
        self.vectors = sorted(count_vectors(H.FULL))
        self.counts = [[(x >> (4 * i)) & 15 for i in range(len(columns))]
                       for x in self.vectors]
        assert all(sum(row[i] * columns[i]['block_size'] for i in range(len(row))) <= H.N
                   for row in self.counts)
        guard()


def solve_and_recover(pack, weights, scale):
    """One local solve, then exact recovery; fallback is an honest bracket."""
    import numpy as np
    from scipy.optimize import linprog
    mat = np.asarray(pack.counts, dtype=float)
    w = np.asarray([float(a / scale) for a in weights])
    result = linprog(-w, A_ub=mat, b_ub=np.ones(len(mat)), bounds=(0, None),
                     method='highs', options={'time_limit': SOLVE_SECONDS})
    solver = dict(method='highs', time_limit_seconds=SOLVE_SECONDS,
        success=bool(result.success), status=int(result.status), message=str(result.message),
        normalized_objective=None if result.fun is None else float(-result.fun),
        iterations=int(getattr(result, 'nit', 0)))
    if result.x is None:
        charges = [F(0)] * len(weights)
        charges[max(range(len(weights)), key=lambda j: weights[j])] = F(1)
    else:
        charges = [max(F(0), F(float(a)).limit_denominator(10**6)) for a in result.x]
    cost = max(sum((m * a for m, a in zip(row, charges)), F(0)) for row in pack.counts)
    if not cost:
        charges[0] = F(1)
        cost = max(row[0] for row in pack.counts)
    charges = [a / cost for a in charges]
    lower = sum((w * a for w, a in zip(weights, charges)), F(0))
    primal = None
    active = []
    marginals = getattr(getattr(result, 'ineqlin', None), 'marginals', None)
    if result.success and marginals is not None:
        active = [i for i, a in enumerate(marginals) if -a > 1e-12]
        guess = [F(float(-marginals[i])).limit_denominator(10**12) * scale for i in active]
        covering = [[pack.counts[i][j] for i in active] for j in range(len(weights))]
        for tolerance in (0, 1e-12, 1e-11, 1e-9, 1e-7):
            tight = [j for j in range(len(weights)) if charges[j] > 0 or (
                tolerance and abs(sum(pack.counts[i][j] * -marginals[i] for i in active) - w[j]) < tolerance)]
            proposed = H.rational_linear_solution(
                [[1] * len(active)] + [covering[j] for j in tight],
                [lower] + [weights[j] for j in tight], guess)
            if proposed is None or any(a < 0 for a in proposed):
                continue
            if all(sum((m * a for m, a in zip(row, proposed)), F(0)) >= weights[j]
                   for j, row in enumerate(covering)):
                primal = proposed
                break
    mixture = []
    if primal is not None:
        mixture = [dict(mass=str(a), blocks=pack.witness(pack.vectors[i]))
                   for i, a in zip(active, primal) if a]
        status = 'EXACT_LOCAL_OPTIMUM'
    else:
        # No second optimization. Rationally retain a proposed mixture, then
        # fill any exact marginal deficits with explicit singleton packings.
        cover = [F(0)] * len(weights)
        if active:
            for i in active:
                a = max(F(0), F(round(float(-marginals[i]) * 10**15), 10**15) * scale)
                if not a:
                    continue
                mixture.append(dict(mass=str(a), blocks=pack.witness(pack.vectors[i])))
                for j, m in enumerate(pack.counts[i]):
                    cover[j] += a * m
        for j, weight in enumerate(weights):
            if cover[j] < weight:
                a = weight - cover[j]
                mixture.append(dict(mass=str(a), blocks=[[pack.ids[j], pack.columns[j]['literal_masks'][0]]]))
        status = 'EXACT_COVER_BRACKET' if result.success else 'OPEN_LOCAL_SOLVER_WITH_EXACT_COVER'
    upper = sum((F(item['mass']) for item in mixture), F(0))
    return dict(status=status, charges=[[j, str(a)] for j, a in zip(pack.ids, charges) if a],
        packing_mixture=mixture, solver=solver, exact_lower_bound=str(lower),
        exact_cover_mass=str(upper))


def verify_certificate(pack, weights, certificate):
    q = {j: F(a) for j, a in certificate['charges']}
    assert len(q) == len(certificate['charges']) and set(q) <= set(pack.ids)
    assert all(a > 0 for a in q.values())
    charges = [q.get(j, F(0)) for j in pack.ids]
    cost = max(sum((m * a for m, a in zip(row, charges)), F(0)) for row in pack.counts)
    assert cost == 1
    lower = sum((w * a for w, a in zip(weights, charges)), F(0))
    cover = [F(0)] * len(weights)
    upper = F(0)
    for item in certificate['packing_mixture']:
        a = F(item['mass'])
        assert a > 0
        counts = pack.verify_witness(item['blocks'])
        for j, m in enumerate(counts):
            cover[j] += a * m
        upper += a
    assert all(a >= w for a, w in zip(cover, weights))
    assert lower == F(certificate['exact_lower_bound'])
    assert upper == F(certificate['exact_cover_mass']) and upper >= lower
    if certificate['status'] == 'EXACT_LOCAL_OPTIMUM':
        assert upper == lower
    return lower, upper, [[j, str(a - w)] for j, a, w in zip(pack.ids, cover, weights) if a != w]


def run(saved=None):
    data, clique, pilot = H.read(H.SUBSETS), H.read(H.CLIQUES), H.read(H.TAG)
    alpha = {i: F(a) for i, a in data['face_multipliers']}
    assert data['face_multipliers'] == clique['face_multipliers'] == pilot['face_multipliers']
    mu = {v: F(a) for v, a in pilot['vertex_multipliers']}
    assert pilot['vertex_multipliers'] == clique['vertex_multipliers']
    weights = defaultdict(F)
    for f in data['positive_simplex_features']:
        for j, m in f['block_orbits']:
            weights[j] += alpha[f['face']] * m
    literals = defaultdict(list)
    for c in data['literal_block_orbits']:
        literals[c['vertex']].append(c)
    for cs in literals.values():
        cs.sort(key=lambda c: c['index'])
    forests = {r['vertex']: r['forest'] for r in data['subset_maxima']}
    reused = {r['vertex']: r for r in pilot['local_cones']}
    assert len(literals) == 740 and len(reused) == 14
    order = sorted(literals, key=lambda v: (-mu[v], -len(literals[v]), v))
    inputs = {f'results/{tag}.json': H.sha(f'results/{tag}.json')
              for tag in (H.SUBSETS, H.CLIQUES, H.TAG)}
    progress_key = hashlib.sha256(H.dump(inputs).encode()).hexdigest()
    existing = {}
    if saved is None and PROGRESS.exists():
        for line in PROGRESS.read_text().splitlines():
            p = json.loads(line)
            assert p['inputs_digest'] == progress_key
            assert p['vertex'] not in existing
            existing[p['vertex']] = p['record']
    saved_rows = {} if saved is None else {r['vertex']: r for r in saved['local_cones']}
    repaired = dict(mu)
    records, cuts, open_cases = [], [], []
    solver_calls = 0
    start = time.monotonic()
    for serial, v in enumerate(order):
        w = [weights[c['index']] for c in literals[v]]
        assert all(a > 0 for a in w)
        source = 'REUSED_FROZEN_FOURTEEN_CERTIFICATE' if v in reused else 'NEW_COMPLETE_LOCAL_CONE'
        try:
            pack = Packings(literals[v])
        except ResourceStop as exc:
            rec = dict(vertex=v, forest=forests[v], status='OPEN_ENUMERATION',
                       reason=str(exc), source=source, block_orbits=len(literals[v]))
            open_cases.append(v)
            records.append(rec)
            if saved is not None:
                assert rec == saved_rows[v]
            elif v not in existing:
                with PROGRESS.open('a') as stream:
                    stream.write(H.dump(dict(inputs_digest=progress_key, vertex=v, record=rec)) + '\n')
            continue
        if saved is not None:
            certificate = saved_rows[v]['certificate']
        elif v in existing:
            assert existing[v]['status'] != 'OPEN_ENUMERATION', 'Previously open forest requires a fresh named successor.'
            certificate = existing[v]['certificate']
        elif v in reused:
            old = reused[v]['certificate']
            certificate = dict(status='EXACT_LOCAL_OPTIMUM', charges=old['charges'],
                packing_mixture=old['packing_mixture'], solver=old['solver'],
                exact_lower_bound=old['exact_local_optimum'], exact_cover_mass=old['exact_local_optimum'])
        else:
            certificate = solve_and_recover(pack, w, mu[v])
        if source == 'NEW_COMPLETE_LOCAL_CONE':
            solver_calls += 1
        lower, upper, slack = verify_certificate(pack, w, certificate)
        repaired[v] = max(mu[v], upper)
        if certificate['status'].startswith('OPEN'):
            open_cases.append(v)
        q = {j: F(a) for j, a in certificate['charges']}
        if lower > mu[v]:
            terms = []
            for f in data['positive_simplex_features']:
                value = sum((m * q.get(j, F(0)) for j, m in f['block_orbits']), F(0))
                if value:
                    terms.append([f['face'], str(value)])
            assert sum((alpha[i] * F(a) for i, a in terms), F(0)) == lower
            cuts.append(dict(index=len(cuts), vertex=v, vertex_cost=1,
                positive_simplex_coefficients=terms, exact_load=str(lower), old_vertex_mass=str(mu[v]),
                excess=str(lower - mu[v]), orbit_charges=certificate['charges']))
        rec = dict(vertex=v, forest=forests[v], status=certificate['status'], source=source,
            counts=dict(block_orbits=len(pack.ids), literal_blocks=len(pack.masks),
                packing_count_vectors=len(pack.vectors), cached_vector_entries=pack.entries,
                cached_available_masks=pack.dp.cache_info().currsize,
                positive_mixture_packings=len(certificate['packing_mixture'])),
            block_orbit_marginals=[[j, str(a)] for j, a in zip(pack.ids, w)],
            literal_block_orbits=literals[v],
            count_vector_digest=hashlib.sha256(H.dump(pack.vectors).encode()).hexdigest(),
            certificate=certificate, primal_cover_slacks=slack, old_vertex_mass=str(mu[v]),
            repaired_vertex_mass=str(repaired[v]), lower_excess=str(max(F(0), lower - mu[v])),
            cover_mass_increase=str(repaired[v] - mu[v]))
        records.append(rec)
        if saved is None and v not in existing:
            with PROGRESS.open('a') as stream:
                stream.write(H.dump(dict(inputs_digest=progress_key, vertex=v, record=rec)) + '\n')
        if serial % 20 == 0 or lower > mu[v]:
            print('ALL_LOCAL_CONES', serial + 1, '/', len(order), 'vertex', v,
                  'cuts', len(cuts), 'open', len(open_cases),
                  'seconds', round(time.monotonic() - start, 2), flush=True)
        pack.dp.cache_clear()
    complete = not open_cases
    numerator = F(data['objective_numerator'])
    normalizer = sum(repaired.values(), F(0))
    paths = [f'scripts/{TAG}.py', f'scripts/{H.TAG}.py', f'results/{H.TAG}.json',
        f'scripts/{H.SUBSETS}.py', f'results/{H.SUBSETS}.json',
        f'scripts/{H.CLIQUES}.py', f'results/{H.CLIQUES}.json',
        'results/fresh_simplex_disjoint_ports_20261017.md',
        f'results/{H.POINTED}.json']
    return dict(status='PASS_EXACT_ALL_LOCAL_PACKING_CONE_COVER' if complete else 'PASS_PARTIAL_LOCAL_PACKING_CONES_OPEN_CASES_RECORDED',
        sources={p: H.sha(p) for p in paths}, dimension=13, target=12,
        complete_positive_simplex_forest_coverage=complete,
        resource_limits=dict(cached_vector_entries_per_forest=ENTRY_LIMIT,
            enumeration_seconds_per_forest=ENUM_SECONDS, local_solver_seconds=SOLVE_SECONDS),
        generation_new_local_solver_invocations=solver_calls, reused_frozen_local_certificates=14,
        saved_replay_solver_invocations=0, global_solver_invocations=0,
        counts=dict(positive_simplex_forests=len(order), processed_forests=len(records), open_forests=len(open_cases),
            exact_local_optima=sum(r['status'] == 'EXACT_LOCAL_OPTIMUM' for r in records),
            exact_cover_brackets=sum(r['status'] == 'EXACT_COVER_BRACKET' for r in records),
            packing_count_vectors=sum(r.get('counts', {}).get('packing_count_vectors', 0) for r in records),
            positive_mixture_packings=sum(r.get('counts', {}).get('positive_mixture_packings', 0) for r in records),
            new_exact_cut_forests=len(cuts), increased_vertex_masses=sum(repaired[v] > a for v, a in mu.items())),
        priority_order=order, open_forests=open_cases, local_cones=records, new_exact_charge_columns=cuts,
        objective_numerator=str(numerator), normalizer=str(normalizer), bound=str(numerator / normalizer),
        excludes_target12=complete and numerator / normalizer > 12,
        face_multipliers=data['face_multipliers'], vertex_multipliers=[[v, str(a)] for v, a in sorted(repaired.items())],
        theorem='Each explicit rational mixture of disjoint literal simplex-port blocks covers all required orbit marginals. Therefore every nonnegative orbit charge of maximum packing cost1 has load at most the mixture mass. Raising each vertex mass to that exact cover mass proves every such local charge inequality simultaneously. Matching charge and mixture values certify the local optimum when recorded. The frozen pilot proves the orbit symmetrization and complete least-port recursion. All current positive simplex rows and all740 contributing forests are included when complete is true.',
        scope='No global LP and no new geometric face census. Arbitrary nonnegative local simplex-block charges are checked on all current positive simplex features, while quad coefficients are zero. Additional geometric incompatibilities, charges coupling simplex and quad configurations, or signed inter-face relations are outside this local packing cone. Face masses, actual quad provenances, all inherited signed transfers and pointed-base constraints are unchanged. Resource stops remain explicitly open. No ambient realization or dimension bound is asserted.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--verify', action='store_true')
    mode.add_argument('--refresh', action='store_true')
    a = p.parse_args()
    out = ROOT / f'results/{TAG}.json'
    saved = json.loads(out.read_text()) if a.verify or a.refresh else None
    result = json.loads(json.dumps(run(saved)))
    if a.verify:
        assert result == saved
    else:
        out.write_text(H.dump(result) + '\n')
    print(result['status'], result['counts'], result['bound'], flush=True)
