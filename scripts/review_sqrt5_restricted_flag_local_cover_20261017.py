#!/usr/bin/env python3
"""Independent exact cover check; no producer imports, optimizer or packing DP."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'review_sqrt5_restricted_flag_local_cover_20261017'
COVER = 'fresh_sqrt5_restricted_flag_local_cones_20261017'
INPUT = 'fresh_sqrt5_flag_local_packing_dual_support_20261017'


def read(tag):
    obj = json.loads((ROOT / f'results/{tag}.json').read_text())
    for name, expected in obj['sources'].items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    return obj


def run():
    top, old = read(COVER), read(INPUT)
    assert len(top['tests']) == 1
    d = s = top['tests'][0]
    assert d['name'] == 'legacy_export_seed'
    c = next(t['repaired_dual'] for t in old['tests'] if t['name'] == d['name'])
    assert d['complete_positive_simplex_forest_coverage'] and not d['open_forests']
    assert d['face_multipliers'] == c['face_multipliers']
    assert d['initial_vertex_multipliers'] == c['vertex_multipliers']
    assert F(d['objective_numerator']) == F(c['objective_numerator'])
    assert F(d['input_bound']) == F(c['bound'])
    alpha = {i: F(a) for i, a in s['face_multipliers']}
    assert len(alpha) == len(s['face_multipliers']) and all(a > 0 for a in alpha.values())
    support = {r['face']: r for r in old['positive_face_support']}
    assert set(alpha) <= set(support)
    expected_s = {i for i in alpha if support[i]['kind'] in
                  ('maximal_simplex_aggregate', 'maximal_edge_incidence')}
    feature_ids = [f['face'] for f in s['positive_simplex_features']]
    assert len(feature_ids) == len(set(feature_ids)) == 1864
    assert set(feature_ids) == expected_s
    required = {i: F(r['objective_contribution']) / F(r['mass']) for i, r in support.items()}
    recomputed_demand = sum((alpha[i] * required[i] for i in alpha), F())
    assert recomputed_demand == F(s['objective_numerator'])
    weights = defaultdict(F)
    for feature in s['positive_simplex_features']:
        assert sum(m for _, m in feature['block_orbits']) == feature['simplex_dimension'] + 1
        for j, multiplicity in feature['block_orbits']:
            assert multiplicity > 0
            weights[j] += alpha[feature['face']] * multiplicity
    columns = {r['index']: r for r in s['literal_block_orbits']}
    assert len(columns) == len(s['literal_block_orbits']) and set(columns) == set(weights)
    byvertex = defaultdict(list)
    for j, r in columns.items():
        assert weights[j] > 0
        assert len(r['literal_masks']) == len(set(r['literal_masks']))
        assert all(0 < mask < 1 << 13 and mask.bit_count() == r['block_size']
                   for mask in r['literal_masks'])
        byvertex[r['vertex']].append(j)
    records = {r['vertex']: r for r in d['local_cones']}
    assert len(records) == len(d['local_cones']) == len(byvertex) == 593
    assert set(records) == set(byvertex)
    mu = {v: F(a) for v, a in c['vertex_multipliers']}
    repaired = dict(mu)
    checked, packing_terms, literal_blocks = [], 0, 0
    for v, js in sorted(byvertex.items()):
        rec = records[v]
        expected = [columns[j] for j in sorted(js)]
        assert rec['literal_block_orbits'] == expected
        assert rec['block_orbit_marginals'] == [[j, str(weights[j])] for j in sorted(js)]
        masks = set()
        for j in js:
            assert not masks.intersection(columns[j]['literal_masks'])
            masks.update(columns[j]['literal_masks'])
        covered = {j: F() for j in js}
        mass = F()
        mixture = rec['certificate']['packing_mixture']
        for item in mixture:
            a = F(item['mass'])
            assert a > 0
            used = 0
            for j, mask in item['blocks']:
                assert j in covered and mask in columns[j]['literal_masks']
                assert not used & mask
                used |= mask
                covered[j] += a
                literal_blocks += 1
            mass += a
            packing_terms += 1
        assert all(covered[j] >= weights[j] for j in js)
        assert mass == F(rec['certificate']['exact_cover_mass'])
        assert F(rec['old_vertex_mass']) == mu[v]
        repaired[v] = max(mu[v], mass)
        assert repaired[v] == F(rec['repaired_vertex_mass'])
        checked.append(dict(vertex=v, orbits=len(js), packing_terms=len(mixture),
                            exact_cover_mass=str(mass), vertex_mass=str(repaired[v])))
    assert [[v, str(a)] for v, a in sorted(repaired.items())] == d['vertex_multipliers']
    numerator = F(s['objective_numerator'])
    denominator = sum(repaired.values(), F())
    bound = numerator / denominator
    assert numerator == F(d['objective_numerator'])
    assert denominator == F(d['normalizer']) and bound == F(d['bound']) > 12
    paths = [f'scripts/{TAG}.py', f'scripts/{COVER}.py', f'results/{COVER}.json',
             f'results/{INPUT}.json']
    return dict(status='PASS_INDEPENDENT_RESTRICTED_593_FOREST_LOCAL_CHARGE_COVER',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        counts=dict(forests=len(checked), block_orbit_marginals=len(weights),
                    explicit_packing_terms=packing_terms, literal_blocks_checked=literal_blocks),
        forests=checked, bound=str(bound), normalizer=str(denominator),
        scope='The exact input face/vertex masses are checked against the frozen legacy restriction. All explicit literal packing mixtures and their rational orbit marginals are checked with an independent standard-library implementation. No producer module, numerical optimizer or count-vector enumeration is used. This proves the upper cover certificate sufficient for every nonnegative charge depending only on the finite forest and block orbit. It does not independently establish complete orbit reconstruction, the exact local lower optima, inherited demand coefficients or signed transfers, realizability, or an ambient dimension theorem.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    result = run()
    out = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(result['status'], result['counts'], result['bound'], flush=True)
