#!/usr/bin/env python3
"""Complete joint 3/4/5-port separation at the final tightened-cap masses."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb
import json
from pathlib import Path

import fresh_sqrt5_joint_triangle_square_ports_20261017 as J
import fresh_sqrt5_cross_coface_cap_full_model_20261017 as N

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_profile_clique_joint_ports_20261017'
BASE = 'fresh_sqrt5_profile_clique_local_cover_update_20261017'
FIVE_DOMAIN = 'fresh_sqrt5_joint_five_port_domain_20261017'
FOUR_DOMAIN = 'fresh_sqrt5_joint_four_port_domain_20261017'
C = J.C


def dump(x): return json.dumps(x, sort_keys=True, separators=(',', ':'))
def vec(): return [F() for _ in range(4)]


def run():
    base, model = N.read(BASE), N.read(N.TAG)
    assert base['complete_positive_simplex_forest_coverage']
    alpha = {i: F(a) for i, a in base['face_multipliers']}
    mu = {v: F(a) for v, a in base['vertex_multipliers']}; updates = dict(mu)
    inverse = {fi: p for p, fi in model['prepared']['provenance_to_face']}
    qfaces = sorted(set(alpha)&set(inverse)); assert len(qfaces) == 336
    features = base['positive_simplex_features']; cols = {c['index']: c for c in base['literal_block_orbits']}
    assert set(alpha) == set(qfaces)|{f['face'] for f in features}
    weights = defaultdict(F); features_by_v = defaultdict(list)
    for f in features:
        touched = set()
        for j, m in f['block_orbits']: weights[j] += m*alpha[f['face']]; touched.add(cols[j]['vertex'])
        for v in touched: features_by_v[v].append(f)
    assert set(weights) == set(cols)
    forests = N.read(J.FORESTS)['vertices']; vlookup = {tuple(v['forest']): v['index'] for v in forests}
    models = dict(C.models()); finite = {C.code(g): name for name, g in models.items()}
    actual, quad, raw, spectators = map(N.read, [J.ACTUAL, J.B.Q.TAG, J.B.Q.CORES, J.B.Q.SPECTATORS])
    qcorners = []; qrecords = []
    for fi in qfaces:
        p = inverse[fi]; rec = actual['provenance_masses'][p]; ci = rec['canonical_prism_core']; core = quad['canonical_cores'][ci]
        g0, _ = J.B.Q.graph(raw['records'][core['source_core']]); g = {i: {k: m[0] for k, m in row.items()} for i, row in g0.items()}
        sp = spectators['contexts'][rec['source_context']]['spectator_forest']
        for name in sp:
            off = len(g)
            for i, row in models[name].items(): g[off+i] = {off+j: m for j, m in row.items()}
        assert len(g) == 15
        cs = J.X.quad_features(dict(quadrilateral_black_partitions=core['quad_partitions']))
        assert len(cs) == 2*rec['prism_dimension']; literal = []
        for c in cs:
            omitted = set(c['vertex_omitted'])
            h = {i: {j: m for j, m in row.items() if j not in omitted} for i, row in g.items() if i not in omitted}
            assert len(h) == 13 and all(m > 0 for row in h.values() for m in row.values())
            pieces = C.pieces(h); forest = tuple(sorted(finite[C.code(p)] for p in pieces)); v = vlookup[forest]
            edges = set(map(tuple, c['occupied_pairs'])); e = len(edges)
            assert 1 <= e == rec['prism_dimension']-1 <= 4
            hubs = set.intersection(*(set(edge) for edge in edges)); assert hubs
            hub = min(hubs); leaves = set().union(*(set(edge) for edge in edges))-{hub}
            assert len(leaves) == e
            qcorners.append(dict(face=fi, vertex=v, graph=h, pieces=pieces, hub=hub, leaves=leaves, edges=e, omitted=sorted(omitted)))
            literal.append(dict(vertex=v, omitted=sorted(omitted), occupied_pairs=c['occupied_pairs']))
        qrecords.append(dict(face=fi, provenance=p, canonical_core=ci, spectator_forest=sp,
            prism_dimension=rec['prism_dimension'], corners=literal))
    vertices = {c['vertex'] for c in cols.values()}|{c['vertex'] for c in qcorners}
    assert vertices <= set(mu)
    graphs = {}
    for v in sorted(vertices):
        g = {}
        for name in forests[v]['forest']:
            off = len(g)
            for i, row in models[name].items(): g[off+i] = {off+j: m for j, m in row.items()}
        assert set(g) == set(range(13)); graphs[v] = g
    five = N.read(FIVE_DOMAIN)['necessary_domain']['facets']; assert len(five) == 15
    facet_families = {
        3: [dict(index=0, coefficients=[2, 3, 0, 0, 1, 1, 1, 1], capacity=3)],
        4: [dict(index=0, coefficients=[1, 2, 3, 0, 0, 1, 1, 1], capacity=3)],
        5: five}
    assert N.read(FOUR_DOMAIN)['status'] == 'PASS_EXACT_SELECTED_DOTTED_PRISM_FOUR_PORT_FACET'
    stages = []; cuts = []
    for size in (3, 4, 5):
        domains = {}; lookup = {}; domain_hash = sha256(); load_hash = sha256()
        for v, g in sorted(graphs.items()):
            pieces = C.pieces(g); groups = defaultdict(list)
            for U in combinations(range(13), size): groups[tuple(sorted(C.code(p, U) for p in pieces))].append(U)
            rows = []
            for oi, (code, sets) in enumerate(sorted(groups.items())):
                rows.append(dict(orbit=oi, representative=sets[0], multiplicity=len(sets)))
                lookup[v, code] = oi
                domain_hash.update((dump([v, oi, code, sets[0], len(sets)])+'\n').encode())
            assert sum(r['multiplicity'] for r in rows) == comb(13, size)
            domains[v] = rows
        scoeff = {}; sload = defaultdict(vec); sident = sha256()
        for serial, (j, c) in enumerate(sorted(cols.items())):
            v = c['vertex']; masks = c['literal_masks']; k = c['block_size']; totals = vec()
            assert masks and len(masks) == len(set(masks)) and all(m.bit_count() == k for m in masks)
            for row in domains[v]:
                oi = row['orbit']; U = sum(1 << i for i in row['representative']); hs = Counter((m & U).bit_count() for m in masks)
                co = tuple(F(hs[h], len(masks)) for h in range(2, 6))
                if any(co):
                    scoeff[j, oi] = co
                    for z, a in enumerate(co): sload[v, oi][z] += weights[j]*a; totals[z] += a*row['multiplicity']
            expected = [comb(k, h)*comb(13-k, size-h) if h <= k and 0 <= size-h <= 13-k else 0 for h in range(2, 6)]
            assert totals == expected; sident.update((dump([j, list(map(str, totals))])+'\n').encode())
            if serial % 1000 == 0: print('JOINT_NEW_S', size, serial+1, '/', len(cols), flush=True)
        qload = defaultdict(vec); qcoeff = defaultdict(lambda: defaultdict(vec)); qident = sha256()
        for serial, corner in enumerate(qcorners):
            fi, v, h, pieces, hub, leaves, e = (corner[k] for k in ('face', 'vertex', 'graph', 'pieces', 'hub', 'leaves', 'edges'))
            hits = defaultdict(Counter)
            for rest in combinations(sorted(set(h)-{hub}), size-1):
                count = len(leaves & set(rest))
                if count:
                    U = (hub,)+rest; code = tuple(sorted(C.code(p, U) for p in pieces)); oi = lookup[v, code]
                    hits[oi][count] += 1
            totals = vec()
            for oi, counts in hits.items():
                mult = domains[v][oi]['multiplicity']; assert sum(counts.values()) <= mult
                for z in range(1, 5):
                    a = F(counts[z], mult)
                    if a:
                        qload[v, oi][z-1] += alpha[fi]*a; qcoeff[v, oi][fi][z-1] += a; totals[z-1] += a*mult
            expected = [comb(e, z)*comb(12-e, size-1-z) if z <= e and 0 <= size-1-z <= 12-e else 0 for z in range(1, 5)]
            assert totals == expected
            qident.update((dump([fi, corner['omitted'], list(map(str, totals))])+'\n').encode())
            if serial % 500 == 0: print('JOINT_NEW_Q', size, serial+1, '/', len(qcorners), flush=True)
        startcuts = len(cuts); check_count = 0; stage_max = defaultdict(F)
        for v, rows in sorted(domains.items()):
            for row in rows:
                oi = row['orbit']; hist = sload[v, oi]+qload[v, oi]
                load_hash.update((dump([v, oi, list(map(str, hist))])+'\n').encode())
                for facet in facet_families[size]:
                    co = facet['coefficients']; cost = facet['capacity']
                    assert len(co) == 8 and min(co) >= 0 and cost > 0
                    load = sum((a*b for a, b in zip(co, hist)), F()); need = load/cost
                    updates[v] = max(updates[v], need); stage_max[v] = max(stage_max[v], need); check_count += 1
                    if need <= mu[v]: continue
                    terms = Counter()
                    for fi, qq in qcoeff[v, oi].items():
                        a = sum((x*y for x, y in zip(co[4:], qq)), F())
                        if a: terms[fi] += a
                    for f in features_by_v[v]:
                        a = F()
                        for j, m in f['block_orbits']:
                            if cols[j]['vertex'] != v: continue
                            cc = scoeff.get((j, oi))
                            if cc: a += m*sum((x*y for x, y in zip(co[:4], cc)), F())
                        if a: terms[f['face']] += a
                    assert sum((alpha[fi]*a for fi, a in terms.items()), F()) == load
                    cuts.append(dict(index=len(cuts), subset_size=size, facet=facet['index'], vertex=v,
                        orbit=oi, representative=list(row['representative']), multiplicity=row['multiplicity'],
                        vertex_cost=cost, load=str(load), old_vertex_mass=str(mu[v]), required_vertex_mass=str(need),
                        positive_face_coefficients=[[fi, str(a)] for fi, a in sorted(terms.items())]))
        assert all(updates[v] >= a for v, a in stage_max.items())
        stages.append(dict(subset_size=size, forests=len(domains), literal_subsets=comb(13, size)*len(domains),
            subset_orbits=sum(map(len, domains.values())), inequalities=check_count, violated_columns=len(cuts)-startcuts,
            complete_orbit_domain_digest=domain_hash.hexdigest(), exact_histogram_digest=load_hash.hexdigest(),
            simplex_binomial_identity_digest=sident.hexdigest(), quad_binomial_identity_digest=qident.hexdigest(),
            vertex_required_masses=[[v, str(a)] for v, a in sorted(stage_max.items())]))
        print('JOINT_NEW_STAGE_PASS', size, check_count, len(cuts)-startcuts, flush=True)
    assert all(updates[v] >= a for v, a in mu.items())
    num = F(base['objective_numerator']); den = sum(updates.values(), F())
    paths = [f'scripts/{TAG}.py', f'scripts/{BASE}.py', f'results/{BASE}.json', f'results/{N.TAG}.json',
        f'scripts/{J.TAG}.py', f'scripts/{J.X.TAG}.py', f'scripts/{J.B.Q.TAG}.py',
        f'results/{J.ACTUAL}.json', f'results/{J.FORESTS}.json', f'results/{J.B.Q.TAG}.json',
        f'results/{J.B.Q.CORES}.json', f'results/{J.B.Q.SPECTATORS}.json',
        f'results/{FIVE_DOMAIN}.json', f'results/{FOUR_DOMAIN}.json']
    return dict(status='PASS_COMPLETE_JOINT_PORT_SCAN_AFTER_TIGHTENED_CAPS',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}, optimizer_invocations=0,
        initial_bound=base['bound'], bound=str(num/den), objective_numerator=str(num), normalizer=str(den),
        counts=dict(positive_faces=len(alpha), simplex_faces=len(features), actual_quad_faces=len(qfaces),
            block_orbits=len(cols), actual_quad_vertex_occurrences=len(qcorners), tested_forests=len(vertices),
            total_inequalities=sum(r['inequalities'] for r in stages), violated_columns=len(cuts),
            increased_vertex_masses=sum(updates[v] > mu[v] for v in mu)),
        face_multipliers=base['face_multipliers'], initial_vertex_multipliers=base['vertex_multipliers'],
        vertex_multipliers=[[v, str(a)] for v, a in sorted(updates.items())], actual_quad_features=qrecords,
        complete_family_stages=stages, critical_columns=cuts,
        scope='Same actual Q provenances and unchanged final tightened-cap face masses. Reconstructs every S block and Q occupied-star contribution on every3/4/5 subset orbit; all strengthened3-port,4-port and15five-port inequalities checked. Exact digests bind full histogram streams; all violated finite columns are exported. Raising only mu preserves every inherited full-model, cap and complete local-packing certificate. These17 families do not exhaust all joint geometric constraints; no global optimizer or dimension theorem.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); a = p.parse_args()
    result = json.loads(dump(run())); out = ROOT/f'results/{TAG}.json'
    if a.verify: assert result == json.loads(out.read_text())
    else: out.write_text(dump(result)+'\n')
    print(result['status'], result['counts'], result['bound'], float(F(result['bound'])), flush=True)
