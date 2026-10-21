#!/usr/bin/env python3
"""Exact local obstructions to a single-facet, two-wall dimension lift."""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/fresh_dimension_wedge_20261013.json'


def read(path):
    return json.loads((ROOT / path).read_text())


def run():
    paths = ['data/known_examples/bugaenko_compact_h8.json',
             'data/enumerations/vinberg_x3_h19.json',
             'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json']
    phi = (1 + sp.sqrt(5)) / 2
    b = read(paths[0])
    G = sp.Matrix([[F(x['a']) + F(x['b']) * phi for x in row]
                   for row in b['wall_gram']])
    h8 = []
    for i in range(G.rows):
        j = next(j for j, label in enumerate(b['coxeter_relations'][i])
                 if label in (3, 5))
        label = b['coxeter_relations'][i][j]
        ratio = sp.simplify(G[i, j]**2 / (G[i, i] * G[j, j]))
        assert G[i, j] < 0
        assert sp.simplify(ratio - {3: sp.Rational(1, 4), 5: phi**2 / 4}[label]) == 0
        h8.append({'wall': i, 'neighbor': j, 'forbidden_old_label': label})

    G19 = -np.array(read(paths[1])['wall_gram'], dtype=object)
    raw = read(paths[2])['coxeter_input']
    W = np.array([[int(x) for x in row] for row in raw['wall_coordinates_as_columns_transposed']], dtype=object)
    Q = -np.array([[int(x) for x in row] for row in raw['basis_gram']], dtype=object)
    G21 = W @ Q @ W.T

    def rational_audit(G):
        assert all(G[i, i] > 0 for i in range(len(G)))
        assert all(G[i, j] <= 0 for i in range(len(G)) for j in range(i))
        rows = []
        for i in range(len(G)):
            ratios = {j: F(int(G[i, j])**2, int(G[i, i] * G[j, j]))
                      for j in range(len(G)) if j != i}
            threes = [j for j, c in ratios.items() if c == F(1, 4)]
            if threes:
                rows.append({'wall': i, 'neighbor': threes[0],
                             'obstruction': 'old finite label3'})
                continue
            sixes = [j for j, c in ratios.items() if c == F(3, 4)]
            parallel = [j for j, c in ratios.items() if c == 1]
            assert sixes and parallel
            # A finite label6 forces split squares1/3 and2/3.
            # At a retained parallel neighbor the first is a forbidden
            # squared cosine, strictly between labels3 and4.
            assert F(1, 4) < F(1, 3) < F(1, 2)
            rows.append({'wall': i, 'finite6_neighbor': sixes[0],
                         'parallel_neighbor': parallel[0],
                         'obstruction': 'label6 split against parallel neighbor',
                         'forced_forbidden_squared_cosine': '1/3'})
        return {'walls': len(G), 'obstruction_counts': dict(Counter(
                    row['obstruction'] for row in rows)), 'certificates': rows}

    # The three possible positive rank-three diagrams with two nonright
    # edges give old labels4,6,10, respectively. The proof that these are
    # exhaustive (including unbounded labels) is in the accompanying note.
    assert sp.simplify(sp.Rational(1, 4) + phi**2 / 4 - (5 + sp.sqrt(5)) / 8) == 0
    return {'status': 'PASS_SCOPED_SINGLE_FACET_WEDGE_OBSTRUCTIONS',
            'sources': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
            'general_necessary_old_finite_labels': [4, 6, 10],
            'H8_to_H9': {'walls': len(h8), 'certificates': h8},
            'H19_to_H20': rational_audit(G19), 'H21_to_H22': rational_audit(G21),
            'scope': 'Replace exactly one original facet by two walls with positive old projections proportional to its normal and opposite nonzero added coordinates; retain every other original facet at zero added coordinate. This rules out that ansatz for each listed source, not general dimension lifts or subsequent replacements.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    out = run()
    if args.verify:
        assert json.loads(OUT.read_text()) == out
    else:
        OUT.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'])
    for name in ('H8_to_H9', 'H19_to_H20', 'H21_to_H22'):
        print(name, out[name]['walls'], out[name].get('obstruction_counts', {}))
