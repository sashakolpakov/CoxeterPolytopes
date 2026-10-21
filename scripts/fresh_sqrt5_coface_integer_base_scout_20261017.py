#!/usr/bin/env python3
"""Exact, nonoptimizing coface-summed integer base-packing diagnostic."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_coface_integer_base_scout_20261017'
FLAGS = 'fresh_sqrt5_flag_refined_base_capacity_20261017'
ACTUAL = 'fresh_sqrt5_actual_provenance_dual_20261017'
MODEL = 'fresh_sqrt5_flag_local_packing_full_model_20261017'
BASE = 'fresh_sqrt5_joint_five_port_diagnostic_20261017'


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def read(stem):
    return json.loads((ROOT / 'results' / (stem + '.json')).read_text())


def run():
    flags, actual, model, base = map(read, [FLAGS, ACTUAL, MODEL, BASE])
    alpha = {i: F(a) for i, a in base['face_multipliers']}
    assert len(alpha) == len(base['face_multipliers']) == 2219
    assert all(a > 0 for a in alpha.values())
    pf = dict(model['prepared']['provenance_to_face'])
    assert len(pf) == len(set(pf.values())) == 1858
    assert len(set(pf.values()) & set(alpha)) == 355
    profiles = {}
    for c in flags['flag_transfer_columns']:
        key = c['base_type']
        value = (c['simplex_face'], c['simplex_dimension'])
        assert profiles.setdefault(key, value) == value
    assert len(profiles) == len({si for si, d in profiles.values()}) == 4969
    assert not ({si for si, d in profiles.values()} & set(pf.values()))
    grouped = defaultdict(list)
    for o in flags['literal_flag_occurrences']:
        grouped[o['provenance'], o['deleted']].append(o)
    bases = []
    counts_by_q = Counter()
    loads = defaultdict(F)
    positive_coefficients = defaultdict(Counter)
    for (j, deleted), rows in sorted(grouped.items()):
        first = rows[0]
        t, q = first['base_type'], first['coface_root']
        assert all(o['base_type'] == t and o['coface_root'] == q for o in rows)
        vertices = sorted(o['canonical_vertex'] for o in rows)
        omitted = [tuple(o['omitted_pair']) for o in rows]
        k = len(vertices)
        assert k == len(set(vertices)) == len(set(omitted))
        rec = actual['provenance_masses'][j]
        assert rec['index'] == j and k == rec['prism_dimension']
        assert 2 <= k <= 5
        si, d = profiles[t]
        assert 1 <= d <= 4 and k <= d + 1
        fi, mass = pf[j], alpha.get(pf[j], F())
        bases.append(dict(provenance=j, face=fi, deleted=deleted, base_type=t,
                          coface_root=q, canonical_vertices=vertices,
                          prism_dimension=k, simplex_face=si, simplex_dimension=d))
        counts_by_q[j] += 1
        for divisor in range(2, d + 2):
            coefficient = k // divisor
            loads[t, divisor] += coefficient * mass
            if coefficient:
                positive_coefficients[t, divisor][j] += coefficient
    assert len(bases) == 6016 and set(counts_by_q) == set(pf)
    # A square has four invariant edge bases; higher prisms have two bases.
    for j, n in counts_by_q.items():
        assert n == (4 if actual['provenance_masses'][j]['prism_dimension'] == 2 else 2)
    checks, violations = [], []
    histogram = Counter()
    positive_histogram = Counter()
    largest_ratio = (F(), None)
    for (t, divisor), load in sorted(loads.items()):
        si, d = profiles[t]
        capacity = (13 - d) * ((d + 1) // divisor)
        assert capacity > 0
        terms = positive_coefficients[t, divisor]
        assert load == sum((n * alpha.get(pf[j], F()) for j, n in terms.items()), F())
        old = alpha.get(si, F())
        slack = capacity * old - load
        row = dict(base_type=t, simplex_face=si, simplex_dimension=d,
                   divisor=divisor, capacity=capacity,
                   provenance_coefficients=[[j, n] for j, n in sorted(terms.items())],
                   load=str(load), simplex_mass=str(old), slack=str(slack))
        checks.append(row)
        histogram[d, divisor] += 1
        if load:
            positive_histogram[d, divisor] += 1
            assert old > 0, ('positive load at absent simplex', t, divisor)
            ratio = load / (capacity * old)
            if ratio > largest_ratio[0]:
                largest_ratio = (ratio, [t, divisor, si])
        if slack < 0:
            violations.append(row)
    assert len(checks) == 7504 and not violations
    assert sum(load > 0 for load in loads.values()) == 1604
    # No alpha or mu is changed: all inherited exact constraints remain as saved.
    sources = {}
    for stem in [FLAGS, ACTUAL, MODEL, BASE]:
        for directory, suffix in [('scripts', '.py'), ('results', '.json')]:
            p = f'{directory}/{stem}{suffix}'
            sources[p] = sha(p)
    for p in [f'scripts/{TAG}.py',
              'results/fresh_prism_simplex_shared_capacity_20261017.md',
              'results/fresh_vertex_weighted_prism_base_capacity_20261017.md']:
        sources[p] = sha(p)
    return dict(status='PASS_EXACT_COFACE_SUMMED_INTEGER_BASE_SCOUT',
                sources=sources, dimension=13, optimizer_invocations=0,
                input_bound=base['bound'], output_bound=base['bound'],
                alpha_changed=False, mu_changed=False,
                counts=dict(positive_faces=len(alpha), positive_actual_quads=355,
                            actual_quad_provenances=len(pf), base_types=len(profiles),
                            actual_bases=len(bases), positive_mass_bases=sum(alpha.get(b['face'], F()) > 0 for b in bases),
                            columns=len(checks), positive_load_columns=1604, violations=0),
                column_counts_by_dimension_and_divisor=[[d, t, n, positive_histogram[d, t]]
                                                       for (d, t), n in sorted(histogram.items())],
                maximum_capacity_utilization=dict(ratio=str(largest_ratio[0]),
                                                   base_type_divisor_simplex=largest_ratio[1]),
                actual_bases=bases, column_checks=checks, violations=violations,
                scope='All coface-summed floor(k/t) columns on the 1858 saved actual prism provenances and 4969 represented maximal-simplex types, at the unchanged 2219-face restricted dual. This does not test labelled coface/subset orbits or the full subset-packing cone. No dimension conclusion.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT / 'results' / (TAG + '.json')
    if args.verify:
        assert json.loads(path.read_text()) == result
    else:
        path.write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(result['status'], result['counts'], result['maximum_capacity_utilization'], flush=True)


if __name__ == '__main__':
    main()
