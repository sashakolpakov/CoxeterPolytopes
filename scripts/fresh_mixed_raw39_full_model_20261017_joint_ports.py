#!/usr/bin/env python3
"""Joint 3/4/5-port separation at the new enlarged-model masses.

Selected actual Q classes only; coefficients on unexpanded merged Q rows
are explicitly zero. No optimizer and no independent provenance selection.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb
import json
from pathlib import Path

import fresh_sqrt5_joint_triangle_square_ports_20261017 as J
import fresh_mixed_raw39_full_model_20261017 as N

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_mixed_raw39_full_model_20261017_joint_ports'
BASE = 'fresh_mixed_raw39_full_model_20261017_local_cones'
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
    features = base['positive_simplex_features']; cols = {c['index']: c for c in base['literal_block_orbits']}
    simplex_faces={f['face']for f in features}
    assert len(alpha)==5059
    full_dual = next(c for c in model['certificates'] if c['kind'] == 'EXACT_PARTIAL_MODEL_DUAL')
    assert base['face_multipliers'] == full_dual['face_multipliers']
    weights = defaultdict(F); features_by_v = defaultdict(list)
    for f in features:
        touched = set()
        for j, m in f['block_orbits']: weights[j] += m*alpha[f['face']]; touched.add(cols[j]['vertex'])
        for v in touched: features_by_v[v].append(f)
    assert set(weights) == set(cols)
    forest_data=N.read(J.FORESTS);forests=forest_data['vertices'];oldfaces=forest_data['faces']; vlookup = {tuple(v['forest']): v['index'] for v in forests}
    models = dict(C.models()); finite = {C.code(g): name for name, g in models.items()}
    actual, quad, raw, spectators = map(N.read, [J.ACTUAL, J.B.Q.TAG, J.B.Q.CORES, J.B.Q.SPECTATORS])
    qchoices={}
    for fi in sorted(alpha):
        if fi in inverse:
            p=inverse[fi];rec=dict(actual['provenance_masses'][p]);rec['provenance_kind']='existing_explicit_map'
        elif fi not in simplex_faces:
            assert fi<len(oldfaces)and oldfaces[fi]['kind']=='canonical_dotted_quad_aggregate'
            contexts=oldfaces[fi]['source_contexts']
            if len(contexts)!=1:continue
            ctx,ci=contexts[0];core=quad['canonical_cores'][ci]
            if core['source_core']not in(23,25,39):continue
            assert oldfaces[fi]['required']==core['required']
            assert quad['source_core_maps'][spectators['contexts'][ctx]['core']]['canonical_core']==ci
            p=None;rec=dict(canonical_prism_core=ci,source_context=ctx,
                prism_dimension={23:5,25:4,39:2}[core['source_core']],
                provenance_kind='unique_saved_raw23_raw25_raw39_context',source_contexts=contexts)
        else:continue
        qchoices[fi]=(p,rec)
    qfaces=sorted(qchoices);explicit=sum(p is not None for p,r in qchoices.values())
    assert explicit==461
    outside=set(alpha)-set(qfaces)-simplex_faces
    assert all(fi<len(oldfaces)and oldfaces[fi]['kind']=='canonical_dotted_quad_aggregate'for fi in outside)
    assert len(qfaces)+len(outside)==2173 and len(features)==2886
    print('JOINT_NEW_DOMAIN',len(features),'S',explicit,'explicit Q',len(qfaces)-explicit,'unique raw23/raw25/raw39 Q',len(outside),'unexpanded Q',flush=True)
    qcorners = []; qrecords = []
    for fi in qfaces:
        p,rec=qchoices[fi];ci=rec['canonical_prism_core'];core=quad['canonical_cores'][ci]
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
        qrecords.append(dict(face=fi, provenance=p, provenance_kind=rec['provenance_kind'],
            source_context=rec['source_context'],source_contexts=rec.get('source_contexts'),canonical_core=ci, spectator_forest=sp,
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
        f'scripts/{N.TAG}.py',
        f'scripts/{J.TAG}.py', f'scripts/{J.X.TAG}.py', f'scripts/{J.B.Q.TAG}.py',
        f'results/{J.ACTUAL}.json', f'results/{J.FORESTS}.json', f'results/{J.B.Q.TAG}.json',
        f'results/{J.B.Q.CORES}.json', f'results/{J.B.Q.SPECTATORS}.json',
        f'results/{FIVE_DOMAIN}.json', f'results/{FOUR_DOMAIN}.json']
    return dict(status='PASS_NEW_MIXED_RAW39_DUAL_SELECTED_Q_JOINT_PORT_SCAN',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}, optimizer_invocations=0,
        initial_bound=base['bound'], bound=str(num/den), objective_numerator=str(num), normalizer=str(den),
        counts=dict(positive_faces=len(alpha), simplex_faces=len(features), actual_quad_faces=len(qfaces),explicit_provenance_Q_rows=explicit,unique_saved_raw23_raw25_raw39_Q_rows=len(qfaces)-explicit,
            unselected_Q_rows_with_zero_new_coefficients=len(outside),
            block_orbits=len(cols), actual_quad_vertex_occurrences=len(qcorners), tested_forests=len(vertices),
            total_inequalities=sum(r['inequalities'] for r in stages), violated_columns=len(cuts),
            increased_vertex_masses=sum(updates[v] > mu[v] for v in mu)),
        face_multipliers=base['face_multipliers'], initial_vertex_multipliers=base['vertex_multipliers'],
        zero_coefficient_unselected_Q_rows=sorted(outside),
        vertex_multipliers=[[v, str(a)] for v, a in sorted(updates.items())], actual_quad_features=qrecords,
        complete_family_stages=stages, critical_columns=cuts,
        scope='Full newly optimized face masses remain fixed. Every positive S block, every positive explicitly split actual Q class, and every positive raw23/raw25/raw39 row with exactly one saved context are reconstructed on every3/4/5 subset orbit. Every remaining merged Q row has zero new coefficients. This is valid restricted nonnegative charging, not a claim to have reconstructed their geometry. All17 selected joint families are checked; exact histogram digests and violated columns are saved. Only mu increases, preserving every inherited finite-model column, including52 vertex-cost signed incidences, and the new complete local simplex covers. No optimizer, all-Q geometric coverage, optimality or dimension theorem.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); a = p.parse_args()
    result = json.loads(dump(run())); out = ROOT/f'results/{TAG}.json'
    if a.verify: assert result == json.loads(out.read_text())
    else: out.write_text(dump(result)+'\n')
    print(result['status'], result['counts'], result['bound'], float(F(result['bound'])), flush=True)
