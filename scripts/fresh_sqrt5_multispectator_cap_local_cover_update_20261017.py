#!/usr/bin/env python3
"""No-optimizer update after the new two-spectator geometric cap family."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import fresh_sqrt5_cross_coface_cap_full_model_20261017 as N
import fresh_sqrt5_restricted_flag_local_cones_20261017 as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_multispectator_cap_local_cover_update_20261017'
OLD = 'fresh_sqrt5_profile_clique_local_cover_update_20261017'
JOINT = 'fresh_sqrt5_profile_clique_joint_ports_20261017'


def dump(x): return json.dumps(x, sort_keys=True, separators=(',', ':'))
def key(c): return (c['vertex'], c['block_size'], R.A.freeze(c['marked_forest_code']))


def run(cap_tag):
    old, joint, caps = N.read(OLD), N.read(JOINT), N.read(cap_tag)
    assert old['complete_positive_simplex_forest_coverage']
    assert joint['face_multipliers'] == old['face_multipliers']
    assert len(caps['columns']) == 555
    inherited_caps = N.read(old['cap_package'])
    assert len(inherited_caps['columns']) == 851
    assert caps['status'] == 'PASS_EXACT_TWO_A1_MULTISPECTATOR_CAP_CAPACITIES'
    d, prepared = N.build(); assert json.loads(json.dumps(prepared)) == N.read(N.TAG)['prepared']
    olda = {i: F(a) for i, a in old['face_multipliers']}; alpha = dict(olda)
    oldmu = {v: F(a) for v, a in joint['vertex_multipliers']}
    assert all(oldmu.get(v, F()) >= F(a) for v, a in old['vertex_multipliers'])
    pmap = dict(prepared['provenance_to_face'])
    for c in caps['columns']:
        assert [[pmap[p], m] for p, m in c['provenance_coefficients']] == c['positive_face_coefficients']
        assert c['capacity'] > 0
        load = sum((m*olda.get(fi, F()) for fi, m in c['positive_face_coefficients']), F())
        fi = c['simplex_face']; alpha[fi] = max(alpha.get(fi, F()), load/c['capacity'])
    alpha = {i: a for i, a in alpha.items() if a}
    initial = N.N.M.dual_certificate(d, dict(alpha)); N.N.verify(d, initial)
    assert {i: F(a) for i, a in initial['face_multipliers']} == alpha
    qrows = {f['index'] for f in d['faces'] if f['kind'] == 'canonical_dotted_quad_aggregate'}
    assert all(alpha.get(i, F()) == olda.get(i, F()) for i in qrows)
    assert sum(i in qrows for i in alpha) == 336
    mu = dict(oldmu)
    for v, a in initial['vertex_multipliers']: mu[v] = max(mu.get(v, F()), F(a))
    print('MULTISPECTATOR_ALPHA', len(alpha), sum(a > olda.get(i, F()) for i, a in alpha.items()), flush=True)
    columns, features, dimensions, corners = R.A.corner_features(d, alpha)
    columns, features = json.loads(dump(columns)), json.loads(dump(features))
    lookup = {key(c): c['index'] for c in columns}; assert len(lookup) == len(columns)
    oldcolumns = {c['index']: c for c in old['block_orbit_columns']}
    remap = {j: lookup[key(c)] for j, c in oldcolumns.items()}
    assert len(set(remap.values())) == len(remap)
    oldliteral = {c['index']: c for c in old['literal_block_orbits']}
    literals = {remap[j]: dict(c, index=remap[j]) for j, c in oldliteral.items()}
    missing = [c for c in columns if c['index'] not in literals]
    if missing:
        for c in R.literal_orbits(d, missing): literals[c['index']] = c
    assert set(literals) == set(range(len(columns)))
    newfeatures = {f['face']: f for f in features}
    for f in old['positive_simplex_features']:
        assert sorted([remap[j], m] for j, m in f['block_orbits']) == newfeatures[f['face']]['block_orbits']
    weights = defaultdict(F); oldweights = defaultdict(F)
    for f in features:
        assert sum(m for j, m in f['block_orbits']) == f['simplex_dimension']+1
        for j, m in f['block_orbits']: weights[j] += m*alpha[f['face']]
    for f in old['positive_simplex_features']:
        for j, m in f['block_orbits']: oldweights[remap[j]] += m*olda[f['face']]
    assert all(weights[j] >= a for j, a in oldweights.items())
    byvertex = defaultdict(list)
    for c in columns: byvertex[c['vertex']].append(c['index'])
    oldcovers = {r['vertex']: r for r in old['local_covers']}
    records = []; inherited_terms = additions_count = 0
    for v, js in sorted(byvertex.items()):
        covered = {j: F() for j in js}; mass = F(); remapped_terms = []
        if v in oldcovers:
            rec = oldcovers[v]
            assert {remap[c['index']] for c in rec['literal_block_orbits']} <= set(js)
            terms = rec['inherited_remapped_packing_mixture'] + rec['singleton_packing_increments']
            for term in terms:
                a = F(term['mass']); assert a > 0
                used = 0; blocks = []
                for oldj, mask in term['blocks']:
                    j = remap[oldj]
                    assert j in covered and mask in literals[j]['literal_masks'] and not used & mask
                    used |= mask; covered[j] += a; blocks.append([j, mask])
                mass += a; inherited_terms += 1
                remapped_terms.append(dict(mass=str(a), blocks=blocks))
            assert mass == F(rec['new_cover_mass'])
            assert all(F(a) == oldweights[remap[j]] for j, a in rec['required_marginals'])
            assert all(covered[j] >= oldweights[j] for j in covered)
        oldmass = mass; added = []
        for j in sorted(js):
            a = max(F(), weights[j]-covered[j])
            if a:
                block = [j, literals[j]['literal_masks'][0]]
                added.append(dict(mass=str(a), blocks=[block])); covered[j] += a; mass += a; additions_count += 1
        assert all(covered[j] >= weights[j] for j in js)
        mu[v] = max(mu.get(v, F()), mass)
        records.append(dict(vertex=v, forest=d['budgets'][v]['forest'], old_cover_mass=str(oldmass),
            new_cover_mass=str(mass), literal_block_orbits=[literals[j] for j in sorted(js)],
            required_marginals=[[j, str(weights[j])] for j in sorted(js)],
            covered_marginals=[[j, str(covered[j])] for j in sorted(js)],
            inherited_remapped_packing_mixture=remapped_terms, singleton_packing_increments=added))
    # Direct check at the final masses, rather than only a monotonicity assertion.
    assert inherited_terms == 3507+61
    loads = defaultdict(F)
    for fi, a in alpha.items():
        for j, m in d['faces'][fi]['coefficients']: loads[j] += a*m
    for c in d['variables']:
        if c['vertex'] is None: assert loads[c['index']] <= 0
        else: assert loads[c['index']] <= c['cost']*mu.get(c['vertex'], F())
    cap_checks = []
    for source, data in [(old['cap_package'], inherited_caps), (cap_tag, caps)]:
        for c in data['columns']:
            load = sum((m*alpha.get(i, F()) for i, m in c['positive_face_coefficients']), F())
            rhs = c['capacity']*alpha.get(c['simplex_face'], F())
            assert load <= rhs
            cap_checks.append([source, c['index'], str(load), str(rhs)])
    num = sum((d['faces'][fi]['required']*a for fi, a in alpha.items()), F()); den = sum(mu.values(), F())
    assert num == F(initial['objective_numerator'])
    paths = [f'scripts/{TAG}.py', f'scripts/{N.TAG}.py', f'results/{N.TAG}.json',
        f'scripts/{R.TAG}.py', f'scripts/{R.A.TAG}.py', f'scripts/{OLD}.py', f'results/{OLD}.json',
        f'scripts/{JOINT}.py', f'results/{JOINT}.json',
        f'scripts/{cap_tag}.py', f'results/{cap_tag}.json', f'results/{old["cap_package"]}.json']
    return dict(status='PASS_EXACT_MULTISPECTATOR_CAP_INCREMENTAL_COMPLETE_LOCAL_COVER',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        cap_package=cap_tag, inherited_model_sha256=prepared['model_sha256'], optimizer_invocations=0,
        initial_bound=joint['bound'], bound=str(num/den), objective_numerator=str(num), normalizer=str(den),
        face_multipliers=[[i, str(a)] for i, a in sorted(alpha.items())],
        vertex_multipliers=[[v, str(a)] for v, a in sorted(mu.items()) if a],
        simplex_mass_raises=[[i, str(olda.get(i, F())), str(a)] for i, a in sorted(alpha.items()) if a > olda.get(i, F())],
        counts=dict(positive_faces=len(alpha), simplex_faces=len(features), actual_quad_faces=336,
            forests=len(records), block_orbits=len(columns), new_block_orbits=len(missing),
            old_packing_terms=inherited_terms, added_singleton_terms=additions_count,
            raised_simplex_masses=sum(a > olda.get(i, F()) for i, a in alpha.items()),
            increased_vertex_masses=sum(a > oldmu.get(v, F()) for v, a in mu.items()),
            inherited_signed_columns=len(d['zeros']), cap_checks=len(cap_checks)),
        complete_positive_simplex_forest_coverage=True, positive_simplex_features=features,
        block_orbit_columns=columns, literal_block_orbits=[literals[j] for j in sorted(literals)],
        old_to_new_block_orbit_map=sorted(map(list, remap.items())), local_covers=records, cap_checks=cap_checks,
        inherited_joint_status='The previous final joint-port vertex masses are retained as a lower bound, but changed simplex masses require a fresh full joint3/4/5 scan; that is a separate successor.',
        scope='No local or global optimization. Canonical block keys transport all3507 original literal mixtures and61 prior singleton increments, revalidated at their old required marginals; exact new deficits receive singleton packings. This proves complete nonnegative local-simplex charge coverage at the new masses. Every inherited model column and all851old+555new cap columns are checked. Complete joint3/4/5 families at changed alpha remain a separate obligation; no optimality or dimension claim.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--caps'); p.add_argument('--verify', action='store_true'); a = p.parse_args()
    out = ROOT/f'results/{TAG}.json'; saved = json.loads(out.read_text()) if a.verify else None
    cap_tag = saved['cap_package'] if saved else a.caps
    assert cap_tag and (not a.caps or a.caps == cap_tag)
    result = json.loads(dump(run(cap_tag)))
    if a.verify: assert result == saved
    else: out.write_text(dump(result)+'\n')
    print(result['status'], result['counts'], result['bound'], float(F(result['bound'])), flush=True)
