#!/usr/bin/env python3
"""Exact no-optimizer repair of the cross-coface cap capacities.

Rebuilds the inherited full matrix; rechecks all local packing mixtures and
updates every affected joint three/four/five-port orbit at the changed alpha.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from math import comb
import json
from pathlib import Path

import fresh_sqrt5_restricted_flag_local_cones_20261017 as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_cross_coface_cap_repair_20261017'
CAP = 'fresh_sqrt5_cross_coface_cap_capacity_20261017'
BASE = 'fresh_sqrt5_joint_five_port_diagnostic_20261017'
THREE = 'fresh_sqrt5_joint_triangle_square_ports_20261017'
FOUR = 'fresh_sqrt5_joint_four_port_diagnostic_20261017'
POINTED = 'fresh_sqrt5_pointed_simplex_base_capacity_20261017'


def read(tag):
    x = json.loads((ROOT/f'results/{tag}.json').read_text())
    for obj in (x, x.get('prepared', {})):
        for p, h in obj.get('sources', {}).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    return x


def run():
    cap, base, prior, three, four, pointed = map(read, [CAP, BASE, R.TAG, THREE, FOUR, POINTED])
    t = next(x for x in prior['tests'] if x['name'] == 'legacy_export_seed')
    assert t['complete_positive_simplex_forest_coverage'] and not t['open_forests']
    assert base['face_multipliers'] == three['face_multipliers'] == four['face_multipliers'] == t['face_multipliers']
    olda = {i: F(a) for i, a in base['face_multipliers']}
    oldmu = {i: F(a) for i, a in base['vertex_multipliers']}
    alpha = dict(olda)
    for c in cap['columns']:
        load = sum((m*olda.get(i, F()) for i, m in c['positive_face_coefficients']), F())
        assert load == F(c['load'])
        s = c['simplex_face']; alpha[s] = max(alpha.get(s, F()), load/c['capacity'])
    alpha = {i: a for i, a in alpha.items() if a}
    delta = {i: a-olda.get(i, F()) for i, a in alpha.items() if a > olda.get(i, F())}
    assert len(delta) == 53 and set(alpha) == set(olda)
    print('REPAIR_ALPHA', len(alpha), len(delta), flush=True)
    d, prepared = R.N.build()
    inherited = read(R.N.TAG)
    assert json.loads(json.dumps(prepared)) == inherited['prepared']
    print('REPAIR_MODEL_REBUILT', prepared['counts'], flush=True)
    profiles = {r['simplex_face']: r for r in pointed['simplex_profiles']}
    columns = {c['index']: c for c in t['literal_block_orbits']}
    features = {f['face']: f for f in t['positive_simplex_features']}
    assert set(delta) <= set(features) and len(features) == 1864
    oldw, neww = defaultdict(F), defaultdict(F)
    dw = defaultdict(F); corner_increments = defaultdict(F)
    for fi, f in features.items():
        assert d['faces'][fi]['required'] == 2*comb(f['simplex_dimension']+1, 3)
        count = Counter()
        for j, m in f['block_orbits']:
            assert m > 0
            oldw[j] += m*olda[fi]; neww[j] += m*alpha[fi]
            count[columns[j]['vertex']] += m
            if fi in delta:
                dw[j] += m*delta[fi]; corner_increments[columns[j]['vertex']] += m*delta[fi]
        assert sum(count.values()) == f['simplex_dimension']+1
        if fi in delta:
            assert count == Counter(v for root, v in profiles[fi]['corners'])
    assert set(oldw) == set(neww) == set(columns)
    assert all(neww[j] == oldw[j]+dw[j] for j in columns)
    uniform = {v: oldmu[v]+corner_increments[v] for v in oldmu}
    assert set(corner_increments) <= set(oldmu)
    mu = dict(oldmu)
    # Every inherited positive and signed column is recounted at the new alpha.
    loads = defaultdict(F)
    for fi, a in alpha.items():
        for j, m in d['faces'][fi]['coefficients']: loads[j] += a*m
    pos = signed = 0
    for c in d['variables']:
        if c['vertex'] is None:
            assert loads[c['index']] <= 0; signed += 1
        else:
            v = c['vertex']; need = loads[c['index']]/c['cost']
            mu[v] = max(mu.get(v, F()), need); pos += 1
    assert signed == 94481
    model_mu = dict(mu)
    # The old exact mixtures are revalidated, then actual orbit deficits are
    # filled with singleton packings. No old cover is used at changed marginals.
    covers = []; inherited_terms = added_terms = 0
    for rec in sorted(t['local_cones'], key=lambda r: r['vertex']):
        v = rec['vertex']; cs = {c['index']: c for c in rec['literal_block_orbits']}
        assert all(c == columns[j] for j, c in cs.items())
        covered = {j: F() for j in cs}; mass = F()
        for term in rec['certificate']['packing_mixture']:
            a = F(term['mass']); assert a > 0
            used = 0
            for j, mask in term['blocks']:
                assert j in cs and mask in cs[j]['literal_masks'] and not used & mask
                used |= mask; covered[j] += a
            mass += a; inherited_terms += 1
        assert mass == F(rec['certificate']['exact_cover_mass'])
        assert all(covered[j] >= oldw[j] for j in cs)
        oldmass = mass; additions = []
        for j in sorted(cs):
            a = max(F(), neww[j]-covered[j])
            if a:
                additions.append(dict(mass=str(a), blocks=[[j, cs[j]['literal_masks'][0]]]))
                mass += a; covered[j] += a; added_terms += 1
        assert all(covered[j] >= neww[j] for j in cs)
        assert mass <= oldmass+corner_increments[v] <= uniform[v]
        mu[v] = max(mu[v], mass)
        covers.append(dict(vertex=v, old_cover_mass=str(oldmass), new_cover_mass=str(mass),
            new_required_marginals=[[j, str(neww[j])] for j in sorted(cs)],
            covered_marginals=[[j, str(covered[j])] for j in sorted(cs)],
            singleton_packing_increments=additions))
    print('REPAIR_LOCAL_COVERS', len(covers), added_terms, flush=True)
    cover_mu = dict(mu)
    changed_by_vertex = defaultdict(list)
    for j, a in dw.items():
        if a: changed_by_vertex[columns[j]['vertex']].append((j, a))
    joint_checks = []; joint_stages = []
    # Frozen complete orbit domains plus the unchanged Q footprints are used.
    # Only changed simplex histograms need reconstruction.
    for size, src, orbit_key in [(3, three, 'finite_forest_triple_orbits'),
                                 (4, four, 'finite_forest_four_port_orbits'),
                                 (5, base, 'finite_forest_five_port_orbits')]:
        domains = {r['vertex']: r['orbits'] for r in src[orbit_key]}
        assert set(domains) == set(oldmu)
        savedloads = {(v, oi): vals for v, oi, *vals in src[
            'all_orbit_histogram_loads' if size == 5 else 'all_orbit_loads_and_slacks']}
        changed = checks = 0; identities = []
        for v, rows in sorted(domains.items()):
            assert sum(r['multiplicity'] for r in rows) == comb(13, size)
            sums = defaultdict(F)
            for row in rows:
                oi = row['orbit']; old = savedloads[v, oi]
                if size == 5:
                    hist = list(map(F, old[0])); assert len(hist) == 8
                else:
                    hist = list(map(F, old[:2])); assert len(hist) == 2
                increment = [F() for _ in range(4 if size == 5 else 1)]
                if v in changed_by_vertex:
                    U = sum(1 << i for i in row['representative'])
                    for j, a in changed_by_vertex[v]:
                        masks = columns[j]['literal_masks']; hs = Counter((m & U).bit_count() for m in masks)
                        if size == 3: co = [F(2*sum(n for h, n in hs.items() if h >= 2)+hs[3], len(masks))]
                        elif size == 4: co = [F(sum(max(0, h-1)*n for h, n in hs.items()), len(masks))]
                        else: co = [F(hs[h], len(masks)) for h in range(2, 6)]
                        for k, value in enumerate(co): increment[k] += a*value
                    changed += 1
                for k, a in enumerate(increment): sums[k] += a*row['multiplicity']
                if size == 5:
                    oldhist = list(hist)
                    for k in range(4): hist[k] += increment[k]
                    for facet in src['facets']:
                        co, cost = facet['coefficients'], facet['capacity']
                        oldload = sum((x*y for x, y in zip(co, oldhist)), F())
                        load = sum((x*y for x, y in zip(co, hist)), F())
                        assert oldload <= cost*oldmu[v]
                        assert load-oldload <= cost*corner_increments[v]
                        mu[v] = max(mu[v], load/cost); checks += 1
                else:
                    oldload = sum(hist); hist[0] += increment[0]; load = sum(hist)
                    assert oldload <= 3*oldmu[v]
                    assert load-oldload <= 3*corner_increments[v]
                    mu[v] = max(mu[v], load/3); checks += 1
                if any(increment):
                    joint_checks.append(dict(size=size, vertex=v, orbit=oi,
                        simplex_increment=list(map(str, increment)), new_loads=list(map(str, hist))))
            if v in changed_by_vertex:
                expected = [F() for _ in range(4 if size == 5 else 1)]
                for j, a in changed_by_vertex[v]:
                    k = columns[j]['block_size']
                    if size == 3: values = [2*comb(k, 2)*(13-k)+3*comb(k, 3)]
                    elif size == 4: values = [sum((h-1)*comb(k, h)*comb(13-k, 4-h) for h in range(2, min(k, 4)+1))]
                    else: values = [comb(k, h)*comb(13-k, 5-h) if h <= k else 0 for h in range(2, 6)]
                    for z, value in enumerate(values): expected[z] += a*value
                assert [sums[k] for k in range(len(expected))] == expected
                identities.append([v, list(map(str, expected))])
        joint_stages.append(dict(subset_size=size, total_inequalities=checks,
            recomputed_changed_orbits=changed, binomial_sum_identities=identities))
        print('REPAIR_JOINT', size, checks, changed, flush=True)
    mu = {v: a for v, a in mu.items() if a}
    assert set(mu) == set(oldmu) and all(oldmu[v] <= mu[v] <= uniform[v] for v in mu)
    for c in d['variables']:
        if c['vertex'] is None: assert loads[c['index']] <= 0
        else: assert loads[c['index']] <= c['cost']*mu.get(c['vertex'], F())
    cap_checks = []
    for c in cap['columns']:
        load = sum((m*alpha.get(i, F()) for i, m in c['positive_face_coefficients']), F())
        assert load <= c['capacity']*alpha.get(c['simplex_face'], F())
        cap_checks.append([c['index'], str(load), str(c['capacity']*alpha.get(c['simplex_face'], F()))])
    num = sum((d['faces'][i]['required']*a for i, a in alpha.items()), F())
    oldnum = sum((d['faces'][i]['required']*a for i, a in olda.items()), F())
    assert oldnum == F(base['objective_numerator'])
    den = sum(mu.values(), F()); uden = sum(uniform.values(), F())
    paths = [f'scripts/{TAG}.py', f'scripts/{R.N.TAG}.py', f'results/{R.N.TAG}.json',
             f'scripts/{R.TAG}.py', f'results/{R.TAG}.json']
    for tag in (CAP, BASE, THREE, FOUR, POINTED):
        paths += [f'scripts/{tag}.py', f'results/{tag}.json']
    return dict(status='PASS_EXACT_CROSS_COFACE_CAP_DUAL_REPAIR',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        inherited_model_sha256=prepared['model_sha256'], optimizer_invocations=0,
        dimension=13, target=12, initial_bound=base['bound'],
        uniform_singleton_bound=str(num/uden), bound=str(num/den),
        objective_numerator=str(num), normalizer=str(den),
        face_multipliers=[[i, str(a)] for i, a in sorted(alpha.items())],
        vertex_multipliers=[[v, str(a)] for v, a in sorted(mu.items())],
        simplex_mass_raises=[[i, str(olda[i]), str(alpha[i])] for i in sorted(delta)],
        uniform_vertex_increments=[[v, str(a)] for v, a in sorted(corner_increments.items()) if a],
        full_model_minimum_vertex_masses=[[v, str(a)] for v, a in sorted(model_mu.items()) if a],
        post_cover_vertex_masses=[[v, str(a)] for v, a in sorted(cover_mu.items()) if a],
        counts=dict(positive_faces=len(alpha), simplex_faces=len(features), quad_faces=355,
            raised_simplex_masses=len(delta), new_positive_simplex_faces=0,
            changed_forests=len(changed_by_vertex), changed_block_orbits=sum(bool(a) for a in dw.values()),
            local_forests=len(covers), block_orbits=len(columns), inherited_packing_terms=inherited_terms,
            added_singleton_packing_terms=added_terms, inherited_faces=len(d['faces']),
            inherited_positive_columns=pos, inherited_signed_columns=signed, new_cap_columns=len(cap_checks)),
        local_packing_covers=covers, changed_joint_orbit_loads=joint_checks,
        joint_scan_stages=joint_stages, cap_checks=cap_checks,
        monotone_other_transfers='All quad masses and common actual provenances are fixed. Only simplex masses increase. Consequently every other inherited Q-positive/S-negative transfer, including the complete coface-base cone capacities, remains valid.',
        scope='Exact repair of one fixed witness, without optimization. Rebuilds the full inherited matrix; revalidates all old literal local packing mixtures and adds actual singleton deficits; recomputes every changed joint3/4/5 orbit and all corresponding inequalities. The repaired ratio is not an optimum, primal feasibility, or dimension theorem. No old packing cover is asserted valid at changed marginals without the saved explicit additions.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); a = p.parse_args()
    result = json.loads(json.dumps(run())); out = ROOT/f'results/{TAG}.json'
    if a.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n')
    print(result['status'], result['counts'], result['bound'], float(F(result['bound'])), flush=True)
