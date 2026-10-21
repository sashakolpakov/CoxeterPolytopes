#!/usr/bin/env python3
"""Exact (simplex vertex, immediate coface) capacities; no optimizer."""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

import fresh_sqrt5_pointed_simplex_base_capacity_20261017 as P
import fresh_sqrt5_prism_simplex_base_incidence_20261017 as B

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_flag_refined_base_capacity_20261017'
M, I, C = P.M, P.I, B.B.C
THEOREM = f'results/{TAG}.md'


def dump(x):
    return json.dumps(x, separators=(',', ':'), sort_keys=True)


def read(tag):
    return M.R.P.read(tag)


def restrict(g, ids):
    ids = set(ids)
    return {i: {j: m for j, m in g[i].items() if j in ids} for i in sorted(ids)}


def flag_orbits(core):
    """Ordered flags, with every accepted map verified on the full graph."""
    g = P.graph(core)
    colors, _ = I.refinement(g)
    black = core['spherical_deletions']
    white = sorted(set(g)-set(black))
    flags = list(product(black, white))
    orbits, tests = [], []
    left = set(flags)
    while left:
        p, q = min(left)
        orbit = []
        for a, b in sorted(left):
            mp = I.isomorphism(g, g,
                {i: (colors[i], i == p, i == q) for i in g},
                {i: (colors[i], i == a, i == b) for i in g})
            if mp is not None:
                assert mp[p] == a and mp[q] == b
                assert {mp[r] for r in black} == set(black)
                assert all(g[i].get(j, 2) == g[mp[i]].get(mp[j], 2)
                           for i, j in combinations(g, 2))
                orbit.append([a, b])
            tests.append(dict(source=[p, q], target=[a, b], equivalent=mp is not None,
                automorphism=None if mp is None else sorted(map(list, mp.items()))))
        assert orbit
        orbits.append(orbit)
        left -= set(map(tuple, orbit))
    assert sum(map(len, orbits)) == len(flags)
    return dict(core=core['index'], black=black, white=white, orbits=orbits, tests=tests)


def run():
    saved = json.loads((ROOT/f'results/{M.TAG}.json').read_text())
    for path, digest in saved['prepared']['sources'].items():
        assert M.R.P.sha(path) == digest, path
    support = read(P.SUPPORT)
    pointed = read(P.TAG)
    vertex = read(M.V.TAG)
    incidence = read(P.INCIDENCE)
    actual = read(M.A.TAG)
    quad = read(B.Q.TAG)
    raw = read(B.Q.CORES)
    profiles = {p['base_type']: p for p in vertex['simplex_vertex_profiles']}
    cores = {ci: flag_orbits(incidence['cores'][ci])
             for ci in sorted({p['core'] for p in profiles.values()})}
    pc = {r['core']: r for r in pointed['core_pointed_orbits']}
    point_index = {ci: {p: oi for oi, orbit in enumerate(r['orbits']) for p in orbit}
                   for ci, r in pc.items()}
    flag_index = {ci: {tuple(f): oi for oi, orbit in enumerate(r['orbits']) for f in orbit}
                  for ci, r in cores.items()}
    models = dict(C.models())
    finite = {C.code(g): name for name, g in models.items()}
    rooted = {}
    for name in sorted({n for p in profiles.values() for n in p['spectator_forest']}):
        g = models[name]
        bycode = defaultdict(list)
        for q in g:
            bycode[C.code(g, [q])].append(q)
        rooted[name] = [dict(roots=roots, code=json.loads(dump(code)))
                        for code, roots in sorted(bycode.items())]
        assert sum(len(r['roots']) for r in rooted[name]) == len(g)
    rooted_index = {(name, dump(r['code'])): oi
                    for name, rows in rooted.items() for oi, r in enumerate(rows)}
    templates = {t['canonical_core']: t for t in incidence['prism_templates']}
    coeff = defaultdict(Counter)
    occurrence_rows, finite_maps = [], []
    recovered = defaultdict(Counter)
    counts = Counter()
    location_cache = {}
    for row in vertex['actual_provenance_base_vertices']:
        j = row['provenance']
        ci = row['canonical_prism_core']
        assert actual['provenance_masses'][j]['canonical_prism_core'] == ci
        tmpl = templates[ci]
        by_deleted = {b['deleted_normal']: b for b in tmpl['bases']}
        g, _ = B.Q.graph(raw['records'][quad['canonical_cores'][ci]['source_core']])
        prism_black = set().union(*(set(left+right) for left, right in quad['canonical_cores'][ci]['quad_partitions']))
        for base in row['bases']:
            t, b, q = base['base_type'], base['deleted_normal'], base['opposite_base_facet']
            p = profiles[t]
            temp = by_deleted[b]
            mp = dict(temp['core_isomorphism'])
            assert q == temp['opposite_base_facet'] and q != b
            assert q not in temp['full_spherical_deletions']
            # S is supported by all surviving nonblack roots.  Joining with
            # Q removes precisely q, hence this literal q identifies U.
            supportS = set(g)-{b}-set(temp['full_spherical_deletions'])
            assert supportS & prism_black == {q}
            assert {mp[r] for r in temp['full_spherical_deletions']} == set(cores[p['core']]['black'])
            key = ci, b
            if key not in location_cache:
                if q in mp:
                    assert mp[q] in cores[p['core']]['white']
                    loc = ('core', mp[q])
                else:
                    h = restrict(g, set(g)-{b})
                    piece = next(c for c in C.pieces(h) if q in c)
                    assert set(piece).isdisjoint(mp)
                    z = {i: {k: m[0] for k, m in rr.items()} for i, rr in piece.items()}
                    assert all(m[0] > 0 for rr in piece.values() for m in rr.values())
                    name = finite[C.code(z)]
                    oi = rooted_index[name, dump(C.code(z, [q]))]
                    target = rooted[name][oi]['roots'][0]
                    gg = models[name]
                    cz, _ = I.refinement(z)
                    cg, _ = I.refinement(gg)
                    iso = I.isomorphism(z, gg,
                        {i: (cz[i], i == q) for i in z},
                        {i: (cg[i], i == target) for i in gg})
                    assert iso is not None and iso[q] == target
                    finite_maps.append(dict(prism_core=ci, deleted=b, coface_root=q,
                        component=name, rooted_orbit=oi, roots=sorted(z),
                        edges=[[i,k,m] for i,rr in z.items() for k,m in rr.items() if i<k],
                        isomorphism=sorted(map(list, iso.items()))))
                    loc = ('finite', name, oi)
                location_cache[key] = loc
            loc = location_cache[key]
            if loc[0] == 'finite':
                assert Counter(p['spectator_forest'])[loc[1]] >= 1
            seen = set()
            for pair, r, forest in base['vertices']:
                assert b in pair and len(pair) == 2
                oldp = next(k for k in pair if k != b)
                assert mp[oldp] == r and r not in seen
                seen.add(r)
                poi = point_index[p['core']][r]
                if loc[0] == 'core':
                    label = ('core', flag_index[p['core']][r, loc[1]])
                else:
                    label = ('finite', poi, loc[1], loc[2])
                coeff[t, label][j] += 1
                recovered[t, poi][j] += 1
                occurrence_rows.append(dict(provenance=j, base_type=t, deleted=b,
                    coface_root=q, omitted_pair=pair, canonical_vertex=r,
                    coface_location=list(loc), flag_label=list(label)))
                counts[loc[0]+'_base_vertex_occurrences'] += 1
            assert len(seen) == tmpl['prism_dimension']
            counts['base_occurrences'] += 1
    oldcolumns = {(c['base_type'], c['pointed_orbit']): c for c in pointed['pointed_transfer_columns']}
    assert set(recovered) <= set(oldcolumns)
    for key, col in oldcolumns.items():
        assert recovered[key] == dict(col['provenance_coefficients'])

    columns, lookup = [], {}
    for t, p in sorted(profiles.items()):
        ci = p['core']
        caps_by_point = Counter()
        labels = []
        for oi, orbit in enumerate(cores[ci]['orbits']):
            points = {point_index[ci][a] for a, b in orbit}
            assert len(points) == 1
            poi = next(iter(points))
            labels.append((('core', oi), poi, len(orbit)))
        for name, multiplicity in sorted(Counter(p['spectator_forest']).items()):
            for poi, orbit in enumerate(pc[ci]['orbits']):
                for oi, rooted_orbit in enumerate(rooted[name]):
                    capacity = len(orbit)*multiplicity*len(rooted_orbit['roots'])
                    labels.append((('finite', poi, name, oi), poi, capacity))
        for label, poi, capacity in labels:
            assert capacity > 0
            caps_by_point[poi] += capacity
            col = dict(index=len(columns), base_type=t, simplex_face=p['simplex_face'],
                core=ci, simplex_dimension=p['simplex_dimension'], flag_label=list(label),
                pointed_orbit=poi, capacity=capacity,
                provenance_coefficients=sorted(map(list, coeff[t,label].items())))
            lookup[t, label] = col['index']
            columns.append(col)
        for poi, orbit in enumerate(pc[ci]['orbits']):
            assert caps_by_point[poi] == (13-p['simplex_dimension'])*len(orbit)
            assert caps_by_point[poi] == oldcolumns[t,poi]['capacity']
        assert sum(caps_by_point.values()) == (13-p['simplex_dimension'])*(p['simplex_dimension']+1)
    assert set(coeff) <= set(lookup)
    for rec in occurrence_rows:
        rec['column'] = lookup[rec['base_type'], tuple(rec['flag_label'])]
    pmap = dict(saved['prepared']['provenance_to_face'])
    assert pmap == dict(pointed['provenance_to_face'])
    for col in columns:
        col['positive_face_coefficients'] = [[pmap[j],m] for j,m in col['provenance_coefficients']]
    initial = support['exported_only_repaired_dual']
    alpha = {i:F(a) for i,a in initial['face_multipliers']}
    assert len(alpha) == len(initial['face_multipliers']) and all(a>0 for a in alpha.values())
    updated, checks, violations = dict(alpha), [], []
    for col in columns:
        load = sum((m*alpha.get(fi,F(0)) for fi,m in col['positive_face_coefficients']),F(0))
        need = load/col['capacity']
        a = alpha.get(col['simplex_face'],F(0))
        checks.append([col['index'],str(load),str(load-col['capacity']*a)])
        if need>a:
            violations.append(col['index'])
            updated[col['simplex_face']] = max(updated.get(col['simplex_face'],F(0)),need)
    print('FLAG_COUNTS',len(columns),dict(counts),'VIOLATIONS',len(violations),
          'RAISED',sum(a>alpha.get(i,F(0)) for i,a in updated.items()),flush=True)
    d, prepared = M.build()
    assert json.loads(json.dumps(prepared)) == saved['prepared']
    M.verify_dual(d,initial)
    repaired = M.dual_certificate(d,updated)
    M.verify_dual(d,repaired)
    aa = {i:F(a) for i,a in repaired['face_multipliers']}
    assert all(aa.get(i,F(0))>=a for i,a in alpha.items())
    for col,(_,load,_) in zip(columns,checks):
        assert F(load) <= col['capacity']*aa.get(col['simplex_face'],F(0))
    assert all(aa.get(fi,F(0))==alpha.get(fi,F(0)) for fi in pmap.values())
    ranked = sorted(violations,key=lambda i: F(checks[i][2]),reverse=True)
    full_initial = next(c for c in saved['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    M.verify_dual(d,full_initial)
    full_alpha = {i:F(a) for i,a in full_initial['face_multipliers']}
    full_updated = dict(full_alpha)
    full_checks, full_violations = [], []
    for col in columns:
        load = sum((m*full_alpha.get(fi,F(0)) for fi,m in col['positive_face_coefficients']),F(0))
        a = full_alpha.get(col['simplex_face'],F(0))
        full_checks.append([col['index'],str(load),str(load-col['capacity']*a)])
        if load>col['capacity']*a:
            full_violations.append(col['index'])
            fi = col['simplex_face']
            full_updated[fi] = max(full_updated.get(fi,F(0)),load/col['capacity'])
    full_repaired = M.dual_certificate(d,full_updated)
    M.verify_dual(d,full_repaired)
    full_aa = {i:F(a) for i,a in full_repaired['face_multipliers']}
    assert all(full_aa.get(i,F(0))>=a for i,a in full_alpha.items())
    assert all(full_aa.get(fi,F(0))==full_alpha.get(fi,F(0)) for fi in pmap.values())
    for col,(_,load,_) in zip(columns,full_checks):
        assert F(load)<=col['capacity']*full_aa.get(col['simplex_face'],F(0))
    print('FLAG_FULL',len(full_violations),'RAISED',
          sum(a>full_alpha.get(i,F(0)) for i,a in full_aa.items()),full_repaired['bound'],flush=True)
    paths = [f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',f'results/{P.TAG}.json',
        f'scripts/{M.TAG}.py',f'results/{M.TAG}.json',f'results/{P.SUPPORT}.json',
        f'scripts/{B.TAG}.py',f'results/{B.TAG}.json',f'results/{M.V.TAG}.json',
        f'results/{M.A.TAG}.json',f'scripts/{B.Q.TAG}.py',f'results/{B.Q.TAG}.json',
        f'results/{B.Q.CORES}.json',f'scripts/{I.TAG}.py',
        'results/fresh_prism_simplex_shared_capacity_20261017.md',
        'results/fresh_vertex_weighted_prism_base_capacity_20261017.md']
    return dict(status='PASS_EXACT_FLAG_REFINED_BASE_CAPACITY_DIAGNOSTIC',
        sources={p:M.R.P.sha(p) for p in paths},dimension=13,target=12,optimizer_invocations=0,
        inherited_model_sha256=prepared['model_sha256'],
        counts=dict(counts,simplex_core_types=len(cores),simplex_profiles=len(profiles),
            actual_provenances=len(vertex['actual_provenance_base_vertices']),
            flag_columns=len(columns),pointed_columns=len(oldcolumns),
            complete_core_flag_isomorphism_tests=sum(len(r['tests']) for r in cores.values()),
            external_root_isomorphism_maps=len(finite_maps),violating_flags=len(violations),
            raised_simplex_rows=sum(a>alpha.get(i,F(0)) for i,a in aa.items()),
            initial_positive_faces=len(alpha),repaired_positive_faces=len(aa)),
        core_flag_orbits=list(cores.values()),finite_root_orbits=rooted,
        external_coface_maps=finite_maps,literal_flag_occurrences=occurrence_rows,
        flag_transfer_columns=columns,initial_flag_loads_and_slacks=checks,
        violating_columns=violations,violations_by_excess=ranked,
        initial_bound=initial['bound'],repaired_dual=repaired,
        full_optimized_test=dict(initial_bound=full_initial['bound'],
            flag_loads_and_slacks=full_checks,violating_columns=full_violations,
            raised_simplex_rows=sum(a>full_alpha.get(i,F(0)) for i,a in full_aa.items()),
            repaired_dual=full_repaired,
            simplex_mass_changes=[[i,str(full_alpha.get(i,F(0))),str(a)]
                                  for i,a in sorted(full_aa.items()) if a>full_alpha.get(i,F(0))]),
        simplex_mass_changes=[[i,str(alpha.get(i,F(0))),str(a)] for i,a in sorted(aa.items())
                              if a>alpha.get(i,F(0))],
        full_star_orbit_rule='The unique connected Lorentz component is preserved. Ordered (black vertex,white coface root) core flags use simultaneous graph automorphisms. For a coface root in a finite orthogonal component, its rooted finite orbit combines independently with the pointed core orbit and all identical component copies; capacities include that multiplicity. Summing flag capacities and positive coefficients exactly recovers every pointed capacity.',
        scope='Exact separation and monotone repair of the frozen exported-only dual. All selected prism provenances use one common mass. Unselected positive coefficients remain zero, yielding valid weaker columns. No optimizer, new face classification, ambient realization or dimension claim. A repaired ratio is a witnessed partial-model lower bound, not an optimum.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()))
    path=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result==json.loads(path.read_text())
    else:
        path.write_text(dump(result)+'\n')
    print(result['status'],result['counts'],result['initial_bound'],result['repaired_dual']['bound'],flush=True)
