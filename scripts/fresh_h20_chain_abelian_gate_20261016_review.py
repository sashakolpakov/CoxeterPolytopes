#!/usr/bin/env python3
"""Independent symbolic/modular review of the nonabelian rotation gate."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import sympy as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_chain_abelian_gate_20261016_review'
PRIMARY = 'results/fresh_h20_chain_abelian_gate_20261016.json'
CHAIN = 'results/fresh_h20_s61_multicap_20261014_exact.json'
OLD = 'results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
BRIDGE = 'results/fresh_h20_e7_bridge_face_20261014_exact.json'
UNIFORM = 'results/fresh_h20_uniform_chain_cap_gate_20261016.json'
UNIFORM_REVIEW = 'results/fresh_h20_uniform_chain_cap_gate_20261016_review.json'


def sha(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def read(path):
    data = json.loads((ROOT/path).read_text())
    for name, expected in data.get('sources', {}).items():
        assert sha(name) == expected, name
    return data


def rational_mod(value, p):
    value = F(value)
    assert value.denominator % p
    return value.numerator*pow(value.denominator, -1, p) % p


def residue(value, p, base_roots, extra_roots=()):
    if extra_roots:
        assert len(value) == 2
        return (residue(value[0], p, base_roots, extra_roots[:-1])+
                extra_roots[-1]*residue(value[1], p, base_roots, extra_roots[:-1])) % p
    assert len(value) == 16
    total = 0
    for mask, coefficient in enumerate(value):
        term = rational_mod(coefficient, p)
        for bit, root in enumerate(base_roots):
            if mask & (1 << bit):
                term = term*root % p
        total = (total+term) % p
    return total


def check_base(p, roots):
    assert S.isprime(p) and len(roots) == 4
    for value, root in zip((2, 79, 11, 3), roots):
        assert root*root % p == value % p and 2*root % p


def run():
    primary, chain, old, bridge, uniform = map(read, (PRIMARY, CHAIN, OLD, BRIDGE, UNIFORM))
    uniform_review = read(UNIFORM_REVIEW)
    assert uniform_review['status'] == 'PASS_INDEPENDENT_ALL_LABEL_UNIFORM_CHAIN_CAP_GATE'
    assert uniform_review['sources'][UNIFORM] == sha(UNIFORM)
    assert primary['status'] == 'PASS_EXACT_NONABELIAN_CHAIN_ROTATION_AND_NONDOTTED_CAP_OBSTRUCTION'
    radD, radH, w = old['field']['discriminant'], bridge['field']['new_positive_sqrt_radical'], chain['frame']['W_norm']
    assert old['field']['base_square_root_basis'] == [1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]
    roots = [193, 198, 295, 230]
    check_base(313, roots)
    d = residue(radD, 313, roots)
    assert d == 244 and S.legendre_symbol(d, 313) == -1
    roots = [340, 325, 332, 196]
    check_base(359, roots)
    d = residue(radD, 359, roots)
    assert d == 246 and 167*167 % 359 == d and 2*167 % 359
    h = residue(radH, 359, roots, (167,))
    assert h == 342 and S.legendre_symbol(h, 359) == -1
    # The four independent rational prime squareclasses give degree16;
    # these two nonsquare residue witnesses double it twice to degree64.
    c, ww = S.symbols('c w')
    A, B = S.chebyshevu(12, c)-S.chebyshevu(11, c), S.chebyshevu(12, c)+S.chebyshevu(11, c)
    # Write b=-z*A, z=sqrt((1+c)/(1-c)); one inverse rotation
    # sends (A,B) to (c*A-(1-c)*B,(1+c)*A+c*B).
    ar, br = S.Integer(1), S.Integer(1)
    for _ in range(12):
        ar, br = S.expand(c*ar-(1-c)*br), S.expand((1+c)*ar+c*br)
    assert S.expand(ar-A) == 0 and S.expand(br-B) == 0
    R, T = (1+c)*A*A, (1-c)*B*B
    assert S.expand(R+T) == 2
    # Multiplying the actual endpoint equation by sqrt(1-c) gives
    # sqrt(2c(w-5/6)) = sqrt((1+c)/2) A + sqrt(w(1-c)) B.
    # Squaring twice yields this polynomial; no converse is needed.
    assert chain['parameter']['equation'] == 'a*sqrt(w-5/6)+b_endpoint/sqrt(2)-c_endpoint*sqrt(w)=0'
    P = S.Poly(S.expand(36*((2*c-T)*ww-S.Rational(5, 3)*c-R/2)**2-72*ww*R*T), c, ww, domain=S.ZZ)
    assert P.degree(c) == 50 and P.degree(ww) == 2
    coeffs = []
    for power in range(3):
        polynomial = S.Poly(S.expand(P.as_expr()).coeff(ww, power), c)
        coeffs.append([int(polynomial.nth(i)) for i in range(polynomial.degree()+1)])
    assert coeffs == primary['endpoint_polynomial_coefficients_by_w_degree']
    certificates = []
    for data in primary['modular_reductions']:
        p, roots, rd, rh = data['p'], data['roots'], data['rd'], data['rh']
        check_base(p, roots)
        d = residue(radD, p, roots)
        h = residue(radH, p, roots, (rd,))
        assert d and h and rd*rd % p == d and rh*rh % p == h
        assert 2*rd % p and 2*rh % p
        wr = residue(w, p, roots, (rd, rh))
        assert wr == data['w']
        f = S.Poly(P.as_expr().subs(ww, wr), c, modulus=p)
        assert f.degree() == 50 and int(f.LC()) % p
        f = f.monic()
        assert list(reversed([int(x) % p for x in f.all_coeffs()])) == data['reduced_monic_polynomial']
        assert S.gcd(f, f.diff()).degree() == 0
        _, factors = S.factor_list(f)
        assert all(multiplicity == 1 and factor.is_irreducible for factor, multiplicity in factors)
        degrees = sorted(factor.degree() for factor, multiplicity in factors)
        assert degrees == ([50] if p == 2017 else [13, 37])
        certificates.append(dict(prime=p, factor_degrees=degrees,
                                 simple_tower_roots=True, leading_unit=True, squarefree=True))
        print('independent symbolic factorization', p, degrees, flush=True)
    assert sorted(row['prime'] for row in certificates) == [359, 2017]
    lo, hi = map(F, chain['parameter']['isolating_interval'])
    assert lo > 0 and lo*lo > 280
    assert 2+4/(lo*lo) < F(141, 70)
    assert list(map(F, uniform['time_interval'])) == [1/hi, 2/lo]
    assert F(uniform['finite_support_cost_census']['strict_above_two_minimum']) > F(141, 70)
    assert all(F(row['whole_parameter_positive_pairing'][0]) > 0 for row in uniform['exact_blockers'])
    # A nondotted block of five has four adjacent differences. If its two
    # internal differences vanish, the recurrence forces the outer two to
    # vanish too; otherwise it expresses c as a ratio in Q^ab.
    assert all(5*dots+4 < 26 for dots in range(5)) and 5*5+4 >= 26
    paths = [f'scripts/{TAG}.py', PRIMARY, CHAIN, OLD, BRIDGE, UNIFORM, UNIFORM_REVIEW]
    return dict(status='PASS_INDEPENDENT_NONABELIAN_ROTATION_AND_FIVE_DOTTED_GATE',
                sources={p: sha(p) for p in paths}, inherited_tower_degree=64,
                independently_reconstructed_endpoint_polynomial_degree=50,
                modular_factorizations=certificates, minimum_dotted_chain_contacts=5,
                proof_dependencies='Frozen exact chain parameter/geometry and uniform-axis gate, including the separate independent all-label support-cost and continuous-blocker review.',
                scope='Independent symbolic Chebyshev polynomial construction, raw radical-tower residues, nonsquare witnesses and SymPy finite-field factorizations; separate mathematical review of Hensel/irreducibility/non-Galois and recurrence arguments. No producer imports. Proves a necessary condition on new walls retaining this133-wall chamber, not finite volume or a dimension bound.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    output = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(output.read_text())
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
