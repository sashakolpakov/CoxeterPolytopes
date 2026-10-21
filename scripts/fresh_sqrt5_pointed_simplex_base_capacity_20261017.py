#!/usr/bin/env python3
"""Pointed complete-simplex-star capacity: exact separation, no optimizer.

Use the already expanded actual quad provenance rows.  New positive transfer
coefficients are zero on unselected merged quad rows, exactly as in the frozen
forest-labelled model.  Test both its full and exported-only exact duals.
"""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

import fresh_sqrt5_labelled_packing_full_model_20261017 as M
import fresh_sqrt5_rank14_core_spectator_domain_20261017 as I

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_pointed_simplex_base_capacity_20261017'
SUPPORT = 'fresh_sqrt5_labelled_packing_dual_support_20261017'
INCIDENCE = 'fresh_sqrt5_prism_simplex_base_incidence_20261017'
THEOREM = 'results/fresh_vertex_weighted_prism_base_capacity_20261017.md'


def graph(core):
    g = {i: {} for i in range(core['nodes'])}
    for i, j, label in core['edges']:
        assert i < j and j not in g[i]
        g[i][j] = g[j][i] = tuple(label)
    # There is exactly one connected Lorentz component in the complete star.
    reached, todo = set(), [0]
    while todo:
        p = todo.pop()
        if p in reached:
            continue
        reached.add(p)
        todo.extend(g[p])
    assert reached == set(g)
    return g


def pointed_orbits(core):
    g = graph(core)
    colors, _ = I.refinement(g)
    roots = core['spherical_deletions']
    assert roots == sorted(set(roots)) and len(roots) == core['simplex_dimension'] + 1
    equivalent = {(p, p): True for p in roots}
    tests = []
    for p, q in combinations(roots, 2):
        mp = I.isomorphism(g, g,
            {i: (colors[i], i == p) for i in g},
            {i: (colors[i], i == q) for i in g})
        equivalent[p, q] = equivalent[q, p] = mp is not None
        if mp is not None:
            assert mp[p] == q and set(mp) == set(g) and set(mp.values()) == set(g)
            assert {mp[r] for r in roots} == set(roots)
            assert all(g[i].get(j, 2) == g[mp[i]].get(mp[j], 2)
                       for i, j in combinations(g, 2))
        tests.append(dict(source_root=p, target_root=q, equivalent=mp is not None,
                          automorphism=None if mp is None else sorted(map(list, mp.items()))))
    orbits = []
    left = set(roots)
    while left:
        p = min(left)
        orbit = sorted(q for q in roots if equivalent[p, q])
        assert set(orbit) <= left
        assert all(equivalent[a, b] for a in orbit for b in orbit)
        orbits.append(orbit)
        left -= set(orbit)
    return dict(core=core['index'], nodes=core['nodes'], edges=core['edges'],
                spherical_deletions=roots, simplex_dimension=core['simplex_dimension'],
                orbits=orbits, all_pointed_pair_tests=tests)


def run():
    saved = json.loads((ROOT / f'results/{M.TAG}.json').read_text())
    for path, digest in saved['prepared']['sources'].items():
        assert M.R.P.sha(path) == digest, path
    support = M.R.P.read(SUPPORT)
    vertex = M.R.P.read(M.V.TAG)
    incidence = M.R.P.read(INCIDENCE)
    actual = M.R.P.read(M.A.TAG)
    d, prepared = M.build()
    assert json.loads(json.dumps(prepared)) == saved['prepared']
    pmap = dict(prepared['provenance_to_face'])
    assert set(pmap) == set(range(1858)) and len(set(pmap.values())) == 1858
    profiles = {p['base_type']: p for p in vertex['simplex_vertex_profiles']}
    cores = {ci: pointed_orbits(incidence['cores'][ci]) for ci in sorted({p['core'] for p in profiles.values()})}
    templates = {p['canonical_core']: p for p in incidence['prism_templates']}
    root_to_orbit = {ci: {p: oi for oi, orbit in enumerate(c['orbits']) for p in orbit}
                     for ci, c in cores.items()}
    forest_columns = {(c['base_type'], c['vertex']): c for c in vertex['labelled_transfer_columns']}
    coeff = defaultdict(Counter)
    reconstructed_forest = defaultdict(Counter)
    base_checks = vertex_checks = 0
    for row in vertex['actual_provenance_base_vertices']:
        j = row['provenance']
        assert actual['provenance_masses'][j]['source_context'] == row['source_context']
        assert actual['provenance_masses'][j]['canonical_prism_core'] == row['canonical_prism_core']
        template = templates[row['canonical_prism_core']]
        by_deleted = {b['deleted_normal']: b for b in template['bases']}
        pc = Counter()
        for base in row['bases']:
            t = base['base_type']
            p = profiles[t]
            b = base['deleted_normal']
            tmpl = by_deleted[b]
            assert tmpl['base_core'] == p['core']
            assert base['opposite_base_facet'] == tmpl['opposite_base_facet']
            mp = dict(tmpl['core_isomorphism'])
            corners = dict(p['corners'])
            assert {mp[r] for r in tmpl['full_spherical_deletions']} == set(corners)
            assert set(corners) == set(cores[p['core']]['spherical_deletions'])
            seen = set()
            for pair, r, v in base['vertices']:
                assert pair == sorted(set(pair)) and len(pair) == 2 and b in pair
                assert tuple(pair) not in seen
                seen.add(tuple(pair))
                other = next(q for q in pair if q != b)
                assert mp[other] == r and corners[r] == v
                oi = root_to_orbit[p['core']][r]
                coeff[t, oi][j] += 1
                reconstructed_forest[t, v][j] += 1
                vertex_checks += 1
            assert len(seen) == template['prism_dimension']
            pc[t] += len(seen)
            base_checks += 1
        assert dict(pc) == dict(actual['provenance_masses'][j]['base_coefficients'])
    for key, c in forest_columns.items():
        assert reconstructed_forest[key] == dict(c['provenance_coefficients'])
    assert set(reconstructed_forest) <= set(forest_columns)

    columns, profile_rows, split_profiles = [], [], []
    for t, p in sorted(profiles.items()):
        core = cores[p['core']]
        corners = dict(p['corners'])
        per_forest = defaultdict(list)
        orows = []
        for oi, roots in enumerate(core['orbits']):
            forests = {corners[r] for r in roots}
            assert len(forests) == 1
            v = next(iter(forests))
            per_forest[v].append(oi)
            cap = (13 - p['simplex_dimension']) * len(roots)
            c = dict(index=len(columns), base_type=t, simplex_face=p['simplex_face'],
                     core=p['core'], simplex_dimension=p['simplex_dimension'],
                     pointed_orbit=oi, spherical_deletion_roots=roots, vertex_forest=v,
                     simplex_vertex_count=len(roots), capacity=cap,
                     provenance_coefficients=sorted(map(list, coeff[t, oi].items())),
                     positive_face_coefficients=[[pmap[j], m] for j, m in sorted(coeff[t, oi].items())])
            columns.append(c)
            orows.append(dict(orbit=oi, roots=roots, vertex=v, column=c['index']))
        profile_rows.append(dict(base_type=t, simplex_face=p['simplex_face'], core=p['core'],
            simplex_dimension=p['simplex_dimension'], spectator_forest=p['spectator_forest'],
            corners=p['corners'], pointed_orbits=orows))
        for v, ids in sorted(per_forest.items()):
            old = forest_columns[t, v]
            pieces = [columns[r['column']] for r in orows if r['orbit'] in ids]
            sums = Counter()
            for c in pieces:
                sums.update(dict(c['provenance_coefficients']))
            assert sums == dict(old['provenance_coefficients'])
            assert sum(c['capacity'] for c in pieces) == old['capacity']
            if len(ids) > 1:
                split_profiles.append(dict(base_type=t, core=p['core'], vertex_forest=v,
                    old_column=old['index'], pointed_columns=[c['index'] for c in pieces]))

    # Certify implication by an existing column when capacities agree and
    # every actual positive coefficient is bounded.  This is stronger than
    # checking one particular dual, and includes all zero-mass provenances.
    implications = []
    for c in columns:
        old = forest_columns[c['base_type'], c['vertex_forest']]
        own = dict(c['provenance_coefficients'])
        if c['capacity'] == old['capacity']:
            assert own == dict(old['provenance_coefficients'])
            witness, kind = old, 'IDENTICAL_OLD_FOREST_COLUMN'
        else:
            candidates = [x for (t, v), x in forest_columns.items()
                          if t == c['base_type'] and x['capacity'] == c['capacity']
                          and all(m <= dict(x['provenance_coefficients']).get(j, 0) for j, m in own.items())]
            assert candidates, c
            witness, kind = min(candidates, key=lambda x: x['index']), 'DOMINATED_BY_OTHER_FOREST_COLUMN'
        assert witness['simplex_face'] == c['simplex_face']
        upper = dict(witness['provenance_coefficients'])
        slack = [[j, upper.get(j, 0) - own.get(j, 0)] for j in sorted(set(upper) | set(own))
                 if upper.get(j, 0) != own.get(j, 0)]
        assert all(m > 0 for j, m in slack)
        implications.append(dict(pointed_column=c['index'], old_forest_column=witness['index'],
                                 kind=kind, positive_coefficient_slack=slack))

    certs = [
        ('full_optimized', next(c for c in saved['certificates'] if c['kind'] == 'EXACT_PARTIAL_MODEL_DUAL')),
        ('exported_only', support['exported_only_repaired_dual']),
    ]
    tests = []
    for name, initial in certs:
        alpha = {i: F(a) for i, a in initial['face_multipliers']}
        assert len(alpha) == len(initial['face_multipliers']) and all(a > 0 for a in alpha.values())
        M.verify_dual(d, initial)
        mu = {v: F(a) for v, a in initial['vertex_multipliers']}
        updated = dict(alpha)
        checks, violations = [], []
        for c in columns:
            load = sum((m * alpha.get(fi, F(0)) for fi, m in c['positive_face_coefficients']), F(0))
            a = alpha.get(c['simplex_face'], F(0))
            value = load - c['capacity'] * a
            checks.append([c['index'], str(load), str(value)])
            if value > 0:
                violations.append(c['index'])
            updated[c['simplex_face']] = max(updated.get(c['simplex_face'], F(0)), load / c['capacity'])
        repaired = M.dual_certificate(d, updated)
        M.verify_dual(d, repaired)
        newalpha = {i: F(a) for i, a in repaired['face_multipliers']}
        newmu = {v: F(a) for v, a in repaired['vertex_multipliers']}
        assert all(newalpha.get(i, F(0)) >= a for i, a in alpha.items())
        assert all(newmu.get(v, F(0)) >= a for v, a in mu.items())
        for c, (_, load, _) in zip(columns, checks):
            assert F(load) <= c['capacity'] * newalpha.get(c['simplex_face'], F(0))
        # The exact implication above predicts no changes; verify that claim
        # against every inherited positive and signed column, not a sample.
        assert not violations and repaired == initial
        tests.append(dict(name=name, initial_bound=initial['bound'],
            counts=dict(positive_faces=len(alpha), positive_forests=len(mu),
                        positive_selected_actual_provenances=sum(alpha.get(fi, F(0)) > 0 for fi in pmap.values()),
                        pointed_violations=len(violations), raised_simplex_rows=0, increased_vertex_masses=0),
            pointed_column_loads_and_signed_slacks=checks, violating_columns=violations,
            repaired_dual=repaired))
        print('POINTED_TEST', name, tests[-1]['counts'], repaired['bound'], flush=True)

    split_core_ids = sorted({p['core'] for p in split_profiles})
    paths = [f'scripts/{TAG}.py', f'scripts/{M.TAG}.py', f'results/{M.TAG}.json',
             f'scripts/{SUPPORT}.py', f'results/{SUPPORT}.json',
             f'scripts/{M.V.TAG}.py', f'results/{M.V.TAG}.json',
             f'scripts/{INCIDENCE}.py', f'results/{INCIDENCE}.json',
             f'scripts/{I.TAG}.py', f'results/{M.A.TAG}.json', THEOREM]
    return dict(status='PASS_EXACT_POINTED_BASE_REFINEMENT_DOMINATED',
        sources={p: M.R.P.sha(p) for p in paths}, dimension=13, target=12,
        inherited_model_sha256=prepared['model_sha256'], optimizer_invocations=0,
        counts=dict(simplex_core_types=len(cores), simplex_profiles=len(profiles),
            all_pointed_root_pairs=sum(len(c['all_pointed_pair_tests']) for c in cores.values()),
            old_forest_labelled_columns=len(forest_columns), pointed_columns=len(columns),
            genuinely_split_core_types=len(split_core_ids), split_simplex_profiles=len(split_profiles),
            identical_old_columns=sum(r['kind'] == 'IDENTICAL_OLD_FOREST_COLUMN' for r in implications),
            dominated_refined_columns=sum(r['kind'] == 'DOMINATED_BY_OTHER_FOREST_COLUMN' for r in implications),
            actual_provenances=len(pmap), actual_base_occurrences=base_checks,
            literal_base_vertex_checks=vertex_checks),
        core_pointed_orbits=list(cores.values()), simplex_profiles=profile_rows,
        split_core_indices=split_core_ids, split_simplex_profiles=split_profiles,
        provenance_to_face=prepared['provenance_to_face'], pointed_transfer_columns=columns,
        old_column_implications=implications, tests=tests,
        complete_star_orbit_argument='The complete normalized S-star is an orthogonal union of one connected Lorentz core and finite positive-definite components. Any weighted graph automorphism preserves components and their Gram signature, so preserves that unique core. Conversely every core automorphism extends by the identity on the spectators. Every vertex of S is a spherical deletion of the core. Thus its complete-star vertex orbit is exactly the computed core-root orbit; spectator permutations cannot merge two distinct orbits. Exact (m,a,b) edge labels are preserved, including dotted values.',
        implication='For each new pointed label, either its signed transfer column equals an inherited forest-labelled column, or its negative simplex coefficient agrees and every positive actual-provenance coefficient is smaller. Hence the existing nonpositive old signed load implies the new load is nonpositive for every nonnegative distribution on the selected actual domain. The two exact dual recounts independently check this implication.',
        scope='Exact structural diagnostic on the frozen rank13 partial model. All1858 selected quad provenances are retained simultaneously, including zero-mass ones. New pointed positive coefficients are zero on unselected merged quad classes; no actual split or omitted base features are inferred for those classes. Every containing simplex uses all spherical one-deletions, while prism vertices are the saved union of genuine K2,2 quadrilateral vertices. No optimizer, new geometric domain, global realization or dimension conclusion.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    result = json.loads(json.dumps(run()))
    out = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(M.R.B.dump(result) + '\n')
    print(result['status'], result['counts'], flush=True)
