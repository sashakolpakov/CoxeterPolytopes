#!/usr/bin/env python3
"""Exact finite input for the sparse-support and compact A1+A1/G2 gates."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_compact_section_coverage_20261014'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUT = ROOT / f'results/{TAG}_gates.json'


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def forest(n, es):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for a, b in es:
        a, b = find(a), find(b)
        if a == b:
            return False
        parent[a] = b
    return True


def run():
    raw = json.loads((ROOT/RAW).read_text())['coxeter_input']
    Q = [[-int(x) for x in row] for row in raw['basis_gram']]
    W = [list(map(int, row)) for row in raw['wall_coordinates_as_columns_transposed']]
    C = [[dot(row, w) for row in Q] for w in W]
    G = [[dot(w, c) for c in C] for w in W]
    norms = [G[i][i] for i in range(210)]
    assert set(norms) == {2, 4}
    for i in range(210):
        assert all(2*x % norms[i] == 0 for x in C[i])
    short = [i for i, q in enumerate(norms) if q == 2]
    pairs = [(a,b) for a,b in combinations(short,2) if G[a][b] == 0]
    table = Counter()
    modified = Counter()
    for a,b in pairs:
        for l in range(210):
            if l in [a,b]:
                continue
            ga,gb,d = G[a][l],G[b][l],norms[l]
            # Schur complement against the orthogonal norm-two pair.
            if 2*d-ga*ga-gb*gb <= 0:
                continue
            table[(d,ga,gb)] += 1
            if ga and gb:
                modified[(d,ga,gb)] += 1
    assert set(modified) == {(2,-1,-1)}
    assert {(norms[i], norms[j], G[i][j]) for i in short for j in range(210)
            if i != j and G[i][j] < 0} == {(2,2,-1),(2,4,-2)}

    rt = s.sqrt(3)
    contact_square = 4*rt-6
    N = s.simplify(2-1/contact_square)
    cos2 = s.simplify(1/(2*N))
    assert s.simplify(contact_square**2/(2-contact_square)**2) == s.Rational(3,4)
    assert bool(contact_square > 0) and bool(contact_square < 1)
    assert bool(N > 0) and bool(N < 1)
    assert bool(cos2 > s.Rational(1,2)) and bool(cos2 < s.Rational(5,8))
    assert bool(s.Rational(5,8) < (3+s.sqrt(5))/8)
    assert bool(1/contact_square > N)  # modified/modified are dotted even at old pairing zero
    assert bool(1/N > 1)              # modified/unchanged long nonright pairs are dotted

    counts = Counter()
    admissible_three = []
    for n in [3,4]:
        all_edges = list(combinations(range(n),2))
        for mask in range(1 << len(all_edges)):
            es = {p for k,p in enumerate(all_edges) if mask >> k & 1}
            if not forest(n, es):
                continue
            for tail in product([-1,1], repeat=n-1):
                signs = (1,)+tail  # quotient global sign reversal
                counts[f'forest_signed_{n}'] += 1
                if max(signs.count(1), signs.count(-1)) >= 3:
                    counts[f'same_sign_triangle_{n}'] += 1
                    continue
                if any(signs[i] != signs[j] and (i,j) not in es for i,j in all_edges):
                    counts[f'obtuse_cross_pair_{n}'] += 1
                    continue
                assert n == 3
                assert len(es) == 2
                majority = [i for i in range(3) if signs.count(signs[i]) == 2]
                center = next(i for i in range(3) if i not in majority)
                assert es == {tuple(sorted((center,j))) for j in majority}
                admissible_three.append(dict(old_edges=sorted(es), signs=signs,
                                             mandatory_new_leaf_edge=majority,
                                             at_least_one_old_edge_must_cancel=True))
    assert len(admissible_three) == 3
    # Three negative finite Coxeter contacts cannot be PD: evaluating on
    # the all-one vector is <=0 after normalizing the three diagonal entries.
    minimum_triangle = s.Matrix([[1,-s.Rational(1,2),-s.Rational(1,2)],
                                 [-s.Rational(1,2),1,-s.Rational(1,2)],
                                 [-s.Rational(1,2),-s.Rational(1,2),1]])
    assert minimum_triangle*s.ones(3,1) == s.zeros(3,1)
    sources = [RAW, f'scripts/{TAG}_gates.py']
    return dict(status='PASS_SPARSE_SUPPORT_AND_COMPACT_ORTHOGONAL_SHORT_G2_INPUT',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                orthogonal_short_pairs=len(pairs),
                spherical_boundary_table=[dict(norm=k[0],pair_a=k[1],pair_b=k[2],count=v)
                                           for k,v in sorted(table.items())],
                unique_modified_boundary=[2,-1,-1],
                equal_contact_square=str(contact_square), modified_norm=str(N),
                modified_unchanged_short_cosine_squared=str(cos2),
                strict_forbidden_interval=['1/2','5/8'],
                sparse_signed_forest_counts=dict(counts),
                three_contact_patterns=admissible_three,
                scope='Sparse support <=3 applies at a finite simple section vertex. '
                      'The compact-only G2 theorem assumes an actual G2 ridge from two old orthogonal norm-two walls. '
                      'It is not a no-section theorem, does not cover noncompact G2 faces, and does not classify all sparse profiles.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = json.loads(json.dumps(run()))
    if args.verify:
        assert out == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'])
    print(out['sparse_signed_forest_counts'])
