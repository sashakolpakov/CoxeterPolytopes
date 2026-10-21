#!/usr/bin/env python3
"""Exact algebra for coupled H3 sections of the saved Borcherds chamber.

No numerical profile discovery is used. Geometric and arithmetic-group
steps, and the precise scope, are in the accompanying proof note.
"""
from pathlib import Path
import argparse
import hashlib
import json
from collections import Counter
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_coupled_h3_obstruction_20261014'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUT = ROOT / f'results/{TAG}.json'
SQ = s.sqrt(5)
C5 = (3 + SQ) / 8
X, U, R = s.symbols('X U R')
CASES = [(2, 2, 2, -1, -1), (2, 4, 2, -2, -1),
         (2, 2, 4, -1, -2), (4, 2, 4, -2, -2),
         (4, 4, 2, -2, -2)]


def rat(z):
    return s.cancel(z, extension=SQ)


def monic(z, v):
    return s.Poly(z, v, extension=SQ).monic().as_expr()


def norm_certificate(z, v, expected_degree):
    p = s.Poly(s.expand(z * z.xreplace({SQ: -SQ})), v, domain=s.QQ)
    factors = s.factor_list(p)[1]
    assert len(factors) == 1
    f, multiplicity = factors[0]
    assert multiplicity == 1 and f.degree() == expected_degree > 1
    assert f.is_irreducible
    return {'degree': f.degree(), 'irreducible_over_Q': True,
            'coefficients_descending': list(map(str, f.monic().all_coeffs()))}


def run():
    raw = json.loads((ROOT / RAW).read_text())['coxeter_input']
    Q = -s.Matrix(raw['basis_gram']).applyfunc(s.Rational)
    W = s.Matrix(raw['wall_coordinates_as_columns_transposed']).applyfunc(s.Rational)
    C = W * Q
    G = C * W.T
    norms = [G[i, i] for i in range(210)]
    assert Counter(norms) == {2: 42, 4: 168}
    pools = {}
    for d, e in [(2, 2), (2, 4), (4, 4)]:
        pool = sorted({G[i, j] for i in range(210) for j in range(i)
                       if sorted([norms[i], norms[j]]) == [d, e]})
        pools[f'{d},{e}'] = list(map(int, pool))
    assert pools == {'2,2': [-1, 0], '2,4': [-2, 0], '4,4': [-8, -6, -4, -2]}
    assert all((2 * C[i, j] / norms[i]).q == 1
               for i in range(210) for j in range(22))

    joined = []
    for cold in [s.Rational(1, 4), s.Rational(1, 2)]:
        y = rat(1 - cold / C5)
        x = rat((1 - y) / (1 + 3 * y))
        ratio = rat(x / y)
        norm = s.simplify(ratio * ratio.xreplace({SQ: -SQ}))
        assert norm in [s.Rational(1, 11), s.Rational(1, 5)]
        assert s.sqrt(norm).is_Rational is False
        assert s.simplify((1 - x) / (1 + 3 * x) - y) == 0
        assert s.expand(s.radsimp(y*y)).coeff(SQ) != 0
        joined.append({'old_cosine_squared': str(cold), 'X': str(x), 'Y': str(y),
                       'ratio_X_over_Y': str(ratio), 'field_norm_of_ratio': str(norm)})

    two = []
    three = []
    for dc, df, dl, gf, gl in CASES:
        old = [dc, df, dl, gf, gl]
        A = rat(dc - gl*gl / (dl*C5))
        P = monic(A*(3*A+dc)*R*R - 8*gf*A*R + 4*gf*gf - (dc-A)*df, R)
        B = monic(s.resultant(P, U-A*R*R, R), U)
        two.append({'old_dc_df_dl_gcf_gcl': old, 'center_contact_squared': str(A),
                    'ratio_polynomial_over_Qsqrt5': str(P),
                    'contact_ratio': norm_certificate(P, R, 4),
                    'leaf_contact_squared': norm_certificate(B, U, 4)})
        for five_on_leaf_pair in [True, False]:
            c1, c2 = ((C5, s.Rational(1, 4)) if five_on_leaf_pair
                      else (s.Rational(1, 4), C5))
            r2 = rat(c1*dl*(df-X)/(X*(c1*df+(1-c1)*X)))
            c0 = rat(gl*gl + gf*gf*r2 - c2*(dc-gf*gf/X)*(dl-X*r2))
            P = monic(rat(c0*c0 - 4*gl*gl*gf*gf*r2).as_numer_denom()[0], X)
            quantity = r2 if five_on_leaf_pair else rat(X*X*r2)
            num, den = quantity.as_numer_denom()
            B = monic(s.resultant(P, U*den-num, X), U)
            three.append({'old_dc_df_dl_gcf_gcl': old,
                          'five_on_leaf_pair': five_on_leaf_pair,
                          'R_squared': str(r2), 'quartic_over_Qsqrt5': str(P),
                          'fixed_leaf_contact_squared': norm_certificate(P, X, 8),
                          'quantity': 'R_squared' if five_on_leaf_pair else 'X_squared_R_squared',
                          'five_edge_contact_ratio_gate': norm_certificate(B, U, 8)})

    A = rat(2 - 1/(2*C5))
    T = 1/A
    B = rat(2*(2-T)/(2+3*T))
    S = rat(((2-A)*(2-B)/4 - 1 - A*B)/2)
    defect = rat(S*S-A*B)
    assert defect != 0
    # Old B3/B4 endpoint contact: normalize its pairing to1. A label4
    # becomes5 with this normal norm. Any nonzero integral pairing with
    # the changed endpoint makes an orthogonal boundary root nonspacelike.
    b_endpoint = []
    for anchor_norm in [2, 4]:
        q = rat((2+SQ)/anchor_norm)
        for wall_norm in [2, 4]:
            minimum = 2 if 4 in [anchor_norm, wall_norm] else 1
            assert s.sign(wall_norm-q*minimum**2) == -1
            b_endpoint.append({'anchor_norm': anchor_norm, 'wall_norm': wall_norm,
                               'normal_norm': str(q), 'smallest_nonzero_pairing': minimum,
                               'projected_norm_upper_bound': str(rat(wall_norm-q*minimum**2))})

    files = [RAW, f'scripts/{TAG}.py',
             'results/fresh_h20_h4_face_closure_20261014.json',
             'results/fresh_h20_h4_face_closure_20261014_review.json']
    return {'status': 'PASS_EXACT_COUPLED_H3_ALGEBRA',
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
            'source_pair_pools': pools, 'all_source_reflections_integral': True,
            'joined_cases': joined, 'two_contact_cases': two,
            'three_contact_cases': three, 'D4_extension_defect': str(defect),
            'B_endpoint_bounds': b_endpoint,
            'scope': 'Uniform local coupled H3 algebra. Compact-face and lattice arguments are in the proof note; no claim of a complete classification of all H20 sections.'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    z = run()
    if args.verify:
        assert z == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(z, indent=2) + '\n')
    print(z['status'])
    print('joined:', len(z['joined_cases']), 'two-contact:', len(z['two_contact_cases']),
          'three-contact:', len(z['three_contact_cases']))
