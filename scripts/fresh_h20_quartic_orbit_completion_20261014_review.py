#!/usr/bin/env python3
"""Independent integer-bank and archive audit; geometry reviewed at source level.

Does not import either scout, rerun its floating geometry, or certify its
numerical exclusions. NumPy is used only to read the saved NPZ arrays.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import product
import hashlib
import json
from math import gcd
from pathlib import Path
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_quartic_orbit_completion_20261014'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
LOCAL = 'results/fresh_h20_rank3_omitted_sections_20261014.json'
WALK = 'results/fresh_h21_source_corner_walk_20261013.json'
BANK = f'data/enumerations/{TAG}_bank.npz'
SCOUT = f'results/{TAG}.json'
EDGES = f'results/{TAG}_edges.json'
OUT = ROOT / 'results' / f'{TAG}_review.json'


def read(p):
    return json.loads((ROOT / p).read_text())


def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def run():
    scout, edges, local, walk = map(read, [SCOUT, EDGES, LOCAL, WALK])
    for saved in [scout, edges]:
        for path, digest in saved['sources'].items():
            assert sha(path) == digest, path
    raw = read(RAW)['coxeter_input']
    Q = [[-int(x) for x in row] for row in raw['basis_gram']]
    W = [tuple(map(int, row)) for row in raw['wall_coordinates_as_columns_transposed']]
    WQ = [[dot(r, col) for col in zip(*Q)] for r in W]
    norms = [dot(r, qr) for r, qr in zip(W, WQ)]
    G = [[dot(qr, s) for s in W] for qr in WQ]
    h = list(map(sum, zip(*W)))
    hgcd = reduce(gcd, h)
    h = [x // hgcd for x in h]
    hQ = [dot(h, col) for col in zip(*Q)]
    assert dot(h, hQ) < 0 and all(dot(w, hQ) < 0 for w in W)
    bank = np.load(ROOT / BANK)
    roots = [tuple(map(int, row)) for row in bank['roots']]
    recipes = [tuple(map(int, row)) for row in bank['source_reflection']]
    expected = set(W)
    for j, r in enumerate(W):
        for i, w in enumerate(W):
            quotient, remainder = divmod(2*G[j][i], norms[j])
            assert remainder == 0
            v = tuple(x-quotient*y for x, y in zip(w, r))
            if dot(v, hQ) > 0:
                v = tuple(-x for x in v)
            expected.add(v)
    assert roots == sorted(expected) and len(roots) == 29043
    counts = Counter()
    for r, (i, j) in zip(roots, recipes):
        if j < 0:
            assert j == -1
            v = W[i]
        else:
            num = 2*G[j][i]
            assert num % norms[j] == 0
            v = tuple(x-(num//norms[j])*y for x, y in zip(W[i], W[j]))
            if dot(v, hQ) > 0:
                v = tuple(-x for x in v)
        assert r == v
        qr = [dot(r, col) for col in zip(*Q)]
        n = dot(r, qr)
        assert n == norms[i]
        counts[n] += 1
    assert dict(counts) == {2: 2163, 4: 26880}
    assert {str(k): v for k, v in counts.items()} == scout['bank_root_norm_counts']

    # If phi(m) divides 8, every prime divisor p of m has p-1 dividing 8.
    # Thus p is 2, 3, or 5; the respective exponents are at most 4, 2, 1.
    labels = sorted({2**a * 3**b * 5**c
                     for a, b, c in product(range(5), range(3), range(2))
                     if 2**a * 3**b * 5**c >= 2
                     and 8 % int(sp.totient(2**a * 3**b * 5**c)) == 0})
    assert labels == scout['angle_label_superset']
    quartic = {i for i, r in enumerate(local['directions'])
               if len(r['minimal_polynomial']) == 5}
    assert len(quartic) == 44
    x = sp.Symbol('x')
    for i in quartic:
        cs = [sp.Rational(Fraction(a)) for a in local['directions'][i]['minimal_polynomial']]
        assert sp.Poly.from_list(cs, x).is_irreducible
    pairs = {(i, u) for i in quartic for u in scout['tested_parameters']}
    assert len(pairs) == 220
    assert {(s['local_record'], s['parameter_u']) for s in scout['sections']} == pairs
    assert len(scout['sections']) == len(edges['sections']) == 220
    sums = Counter()
    for s, e in zip(scout['sections'], edges['sections']):
        assert (s['local_record'], s['parameter_u']) == (e['local_record'], e['parameter_u'])
        spec = local['systems'][local['directions'][s['local_record']]['system']]
        cusp = walk['corner_orbits'][spec['cusp']]
        seed = [i for i in cusp['active_walls'] if i != spec['omitted_wall']]
        assert s['seed_walls'] == e['seed_walls'] == seed
        factors = [[w for w in f['walls'] if w in seed] for f in cusp['affine_components']]
        assert sorted(sum(factors, [])) == sorted(seed)
        assert sum(len(f)-1 for f in factors) == 19
        omissions = set(product(*factors))
        archived = [tuple(a) for a in e['uncovered_edge_omissions']]
        assert len(archived) == len(omissions) and set(archived) == omissions
        assert not e['covered_edges']
        ids = s['compatible_cap_indices'] + s['near_angle_cap_indices']
        assert len(set(ids)) == len(ids) and all(0 <= i < len(roots) for i in ids)
        assert len(ids) == e['bank_caps_including_near_angle_remainders']
        sums['seed_edges'] += len(omissions)
        sums['compatible_cap_proposals'] += len(s['compatible_cap_indices'])
        sums['near_angle_cap_proposals'] += len(s['near_angle_cap_indices'])
        sums['seed_acute_cap_proposals'] += s['positive_acute_caps']
    assert dict(sums) == dict(seed_edges=42420, compatible_cap_proposals=8136,
                             near_angle_cap_proposals=96, seed_acute_cap_proposals=89242)
    assert edges['counts'] == dict(tested_sections=220, seed_edges=42420, covered_seed_edges=0)
    paths = [RAW, LOCAL, WALK, BANK, SCOUT, EDGES,
             f'scripts/{TAG}.py', f'scripts/{TAG}_edges.py', f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_INTEGER_BANK_AND_SOURCE_SCOPE_REVIEW',
                sources={p: sha(p) for p in paths}, bank_size=len(roots),
                bank_norm_counts={str(k): v for k, v in counts.items()},
                quartic_records=44, sections=220, finite_label_superset=labels,
                archive_counts=dict(sums),
                exact_replay=['complete one-reflection integer bank and every recipe',
                              'all bank norms', 'quartic irreducibility and degree label pool',
                              'source manifests and complete edge-omission tuples'],
                source_proof_review=['projection formula', 'cusp-edge support reduction',
                                     'finite-or-ideal Schur necessary inequality'],
                limits='No independent floating geometry rerun or exact endpoint exclusion. '
                       'The 8136 accepted and 96 near-angle entries remain numerical seed-compatible proposals; '
                       'pairwise compatibility of the added caps and finite volume were not certified.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    if args.verify:
        assert json.loads(OUT.read_text()) == out
    else:
        OUT.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'])
    print(out['archive_counts'])
