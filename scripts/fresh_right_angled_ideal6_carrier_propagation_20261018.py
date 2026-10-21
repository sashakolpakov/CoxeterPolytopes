#!/usr/bin/env python3
"""Exact necessary ideal6 and mixed-six-face bounds; no existence/exclusion census."""
from fractions import Fraction
from hashlib import sha256
from math import comb
from pathlib import Path
import argparse
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/fresh_right_angled_ideal6_carrier_propagation_20261018.json'
PINS = {
    'results/right_angled_five_face_carrier_20261018.md':
        'ce90d9a3fa3c217e6811163a48f84d8ccadaafbe9c99eb4d62ef86de36dc1e0d',
    'results/ideal_ra4_second_volume_equality_20261018.md':
        '4d227e9f7f4577f0a61b30a3ebf2790c763feeb4f1acddf7799a83f96fc1d7e2',
}


def build():
    for p, h in PINS.items():
        assert sha256((ROOT / p).read_bytes()).hexdigest() == h
    # The report proves each incidence premise; these are their exact consequences.
    carrier_facets, carrier_vertices = 25, 30
    facets5 = 1 + carrier_facets + carrier_vertices
    vertices5 = 3 * facets5 + 2
    facets6 = 1 + facets5 + vertices5
    vertices6 = Fraction(vertices5 * facets6, 10)
    assert (facets5, vertices5, facets6, vertices6) == (56, 170, 227, 3859)
    rows = []
    for n in (11, 12):
        a, b = n // 2, (n + 1) // 2
        assert 0 <= 5 < 6 <= b
        upper = Fraction(comb(n - 5, n - 6) * (comb(a, 5) + comb(b, 5)),
                         comb(a, 6) + comb(b, 6))
        fraction = (upper - 27) / (facets6 - 27)
        assert upper == 42 and fraction == Fraction(3, 40)
        rows.append({'n': n, 'six_face_facet_average_strict_upper': str(upper),
                     'ideal_six_face_fraction_strict_upper': str(fraction)})
    return {
        'status': 'NECESSARY_IDEAL6_SIZE_AND_MIXED_SIX_FACE_DENSITY',
        'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'sources': PINS,
        'internal_incidence_premises': {
            'ideal5_euler': 'f3 = 3*v + 3*f4 - 6',
            'ideal5_facet_incidence': '2*f3 >= 24*f4',
            'actual_ideal_facet_injection': 'facets(P) >= 1 + facets(Q) + vertices(Q)',
            'ideal6_vertex_facet_incidence': '10*v6 = sum(vertices(Q5))',
        },
        'ideal5_minima': {'facets': facets5, 'vertices': vertices5},
        'ideal6_minima': {'facets': facets6, 'vertices': int(vertices6)},
        'external_inputs': {
            'ideal_RA4_facet_minimum': 24,
            'general_finite_volume_RA6_facet_minimum': 27,
            'general_RA6_derivation': 'Dufour nu5 >= 26 and facets >= 1 + nu5; no compact RA6',
            'NK_range': 'simple at edges; k=6,l=5,n=11,12',
            'primary_statements': 'https://arxiv.org/html/2303.09533v1#S6',
        },
        'mixed_rows': rows,
        'existence_asserted': False,
        'new_dimension_bound': False,
        'optimality_asserted': False,
        'scope': 'exact algebra of the universal written incidence arguments and explicit external bounds',
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    data = build()
    if a.verify:
        assert json.loads(OUT.read_text()) == data
    else:
        assert not OUT.exists()
        OUT.write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'verify': a.verify,
                      'ideal5': data['ideal5_minima'], 'ideal6': data['ideal6_minima'],
                      'mixed_six_fraction': '3/40'}))
