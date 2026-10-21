#!/usr/bin/env python3
"""Independent exact check of the first cross-coface cap witness."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'review_sqrt5_cross_coface_full_base_witness_20261017'
TAGS = dict(raw='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017',
            incidence='fresh_sqrt5_prism_simplex_base_incidence_20261017',
            pointed='fresh_sqrt5_pointed_simplex_base_capacity_20261017',
            actual='fresh_sqrt5_actual_provenance_dual_20261017',
            labelled='fresh_sqrt5_labelled_coface_base_packing_20261017',
            weights='fresh_sqrt5_joint_five_port_diagnostic_20261017')


def K(a=0, b=0): return (F(a), F(b))
def add(x, y): return (x[0]+y[0], x[1]+y[1])
def neg(x): return (-x[0], -x[1])
def sub(x, y): return add(x, neg(y))
def mul(x, y): return (x[0]*y[0]+5*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def div(x, y):
    den = y[0]*y[0]-5*y[1]*y[1]
    assert den
    return mul(x, (y[0]/den, -y[1]/den))
def sign(x):
    a, b = x
    if not a: return (b > 0)-(b < 0)
    if not b or (a > 0) == (b > 0): return (a > 0)-(a < 0)
    return ((a > 0)-(a < 0))*((a*a > 5*b*b)-(a*a < 5*b*b))
def encode(x): return list(map(str, x))
def total(xs):
    ans = K()
    for x in xs: ans = add(ans, x)
    return ans


def inverse_det(G):
    n = len(G)
    a = [list(row)+[K(int(i == j)) for j in range(n)] for i, row in enumerate(G)]
    determinant = K(1)
    for j in range(n):
        i = next(i for i in range(j, n) if a[i][j] != K())
        if i != j:
            a[i], a[j] = a[j], a[i]; determinant = neg(determinant)
        pivot = a[j][j]; determinant = mul(determinant, pivot)
        a[j] = [div(x, pivot) for x in a[j]]
        for i in range(n):
            if i != j:
                scale = a[i][j]
                a[i] = [sub(x, mul(scale, y)) for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a], determinant


def positive_definite(G):
    a = [list(row) for row in G]
    for j in range(len(a)):
        if sign(a[j][j]) <= 0: return False
        for i in range(j+1, len(a)):
            for k in range(j+1, len(a)):
                a[i][k] = sub(a[i][k], div(mul(a[i][j], a[j][k]), a[j][j]))
    return True


def run():
    data = {k: json.loads((ROOT/f'results/{tag}.json').read_text()) for k, tag in TAGS.items()}
    rec = data['raw']['records'][23]
    assert rec['index'] == 23 and rec['nodes'] == 7
    labels = dict(zip(combinations(range(7), 2), rec['row']['labels']))
    nonright = {e: m for e, m in labels.items() if m != 2}
    assert nonright == {(0, 2): 0, (1, 6): 5, (2, 3): 3,
                        (3, 5): 3, (4, 5): 3, (4, 6): 3}
    tau = tuple(map(F, rec['parameters']['0,2']))
    assert tau == K(F(7, 2), F(1, 2))
    template = next(t for t in data['incidence']['prism_templates'] if t['canonical_core'] == 19)
    assert template['source_core'] == 23 and template['prism_dimension'] == 5
    base = next(b for b in template['bases'] if b['deleted_normal'] == 2)
    mp = dict(base['core_isomorphism'])
    assert mp == {1: 4, 3: 3, 4: 1, 5: 2, 6: 0}
    assert base['opposite_base_facet'] == 0 and base['finite_extra'] == ['A1']
    assert base['base_core'] == 1 and base['full_spherical_deletions'] == sorted(mp)
    core = data['incidence']['cores'][1]
    edges = {(i, j): lab[0] for i, j, lab in core['edges']}
    assert core['nodes'] == 5 and edges == {(0, 1): 3, (0, 4): 5, (1, 2): 3, (2, 3): 3}
    def label(i, j): return edges.get(tuple(sorted((i, j))), 2)
    assert all(labels[tuple(sorted((i, j)))] == label(mp[i], mp[j]) for i, j in combinations(mp, 2))
    assert all(labels[tuple(sorted((0, i)))] == 2 for i in mp)
    automorphisms = [list(p) for p in permutations(range(5))
                     if all(label(i, j) == label(p[i], p[j]) for i, j in combinations(range(5), 2))]
    assert automorphisms == [list(range(5))]
    G = [[K(2) if i == j else {2: K(), 3: K(-1), 5: K(F(-1, 2), F(-1, 2))}[label(i, j)]
          for j in range(5)] for i in range(5)]
    inv, det = inverse_det(G)
    assert det == K(4, -2) and sign(det) < 0
    assert all(positive_definite([[G[i][j] for j in range(5) if j != r]
                                 for i in range(5) if i != r]) for r in range(5))
    # A positive definite codimension-one restriction and negative determinant give inertia (4,1).
    p = [K(-int(i == 3)) for i in range(5)]
    assert all(labels[tuple(sorted((2, old)))] == (3 if new == 3 else 2) for old, new in mp.items())
    v = [total(mul(x, y) for x, y in zip(row, p)) for row in inv]
    assert [total(mul(x, y) for x, y in zip(row, v)) for row in G] == p
    v2 = total(mul(x, y) for x, y in zip(p, v))
    assert v2 == K(F(1, 4), F(-1, 4))
    assert add(v2, div(tau, K(2))) == K(2)
    assert sign(v2) < 0 and sign(add(v2, K(1))) > 0
    profile = next(p for p in data['pointed']['simplex_profiles'] if p['base_type'] == 54999)
    assert profile['core'] == 1 and profile['simplex_face'] == 80143
    assert profile['spectator_forest'] == ['A1']*9 and profile['simplex_dimension'] == 4
    provenance = data['actual']['provenance_masses'][440]
    assert provenance['canonical_prism_core'] == 19 and provenance['merged_face'] == 26259
    bases = [b for b in data['labelled']['actual_bases'] if b['provenance'] == 440]
    assert len(bases) == 2 and sorted(b['base_type'] for b in bases) == [54998, 54999]
    B = next(b for b in bases if b['base_type'] == 54999)
    assert B['deleted'] == 2 and B['coface_root'] == 0 and B['canonical_vertices'] == list(range(5))
    assert B['coface_location'] == ['finite', 'A1', 0] and B['face'] == 26259
    alpha = {i: F(a) for i, a in data['weights']['face_multipliers']}
    qmass, smass = alpha[26259], alpha[80143]
    assert qmass == F(7551, 10**10) and smass == F(839, 10**10) and qmass == 9*smass
    paths = [f'scripts/{TAG}.py']+[f'results/{tag}.json' for tag in TAGS.values()]
    return dict(status='PASS_INDEPENDENT_EXACT_CROSS_COFACE_WITNESS',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                arithmetic='Independent stdlib Fraction Q(sqrt(5)); no producer arithmetic imports',
                source_core=23, canonical_prism_core=19, provenance=440, quad_face=26259,
                base_type=54999, simplex_face=80143, coface_A1_copies=9,
                core_gram=[[encode(x) for x in row] for row in G],
                determinant=encode(det), physical_inertia=[4, 1], all_five_deletions_PD=True,
                all_120_permutations_checked=True, core_automorphisms=automorphisms,
                core_contact_vector=[encode(x) for x in p], core_projection=[encode(x) for x in v],
                squared_dotted_product=encode(tau), distinct_copy_product=encode(v2),
                distinct_copy_product_strictly_between_minus_one_and_zero=True,
                signed_column=dict(positive_quad_face=26259, coefficient=1,
                                   negative_simplex_face=80143, capacity=1),
                current_quad_mass=str(qmass), current_simplex_mass=str(smass),
                violation_factor=str(qmass/smass), excess=str(qmass-smass),
                scope='Checks the exact first cross-coface incidence obstruction and saved mass violation. Prior local-star realization and canonical aggregation are inherited. No mass repair, broader family census, ambient realization, or dimension bound is asserted.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--verify', action='store_true'); args = parser.parse_args()
    result = run(); output = ROOT/f'results/{TAG}.json'
    if args.verify: assert json.loads(output.read_text()) == result
    else: output.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(result['status'], result['violation_factor'], result['excess'], flush=True)
