#!/usr/bin/env python3
"""Delete the C11 endpoint half using every native heavy-end profile.

The universal geometry is proved in the accompanying report. No project
arithmetic code or bounded label catalogue is imported.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_order11_c_endpoint_zero_transfer_20261018'
PREVIOUS = 'fresh_order11_bcd_half_pair_transfer_20261018'
BASE = 'fresh_order11_order17_dimension63_20261018'
STEMS = (PREVIOUS, BASE)


def run():
    docs = {tag: json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag, packet in docs.items():
        assert packet['status'].startswith('PASS'), tag
        for path, digest in packet['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, (tag, path)
    old = next(r for r in docs[PREVIOUS]['factor_rules'] if r['kind'] == 'C')
    assert old['changed_pairs'] == [[0, 10]] and old['total_cost'] == '109/2'
    assert docs[BASE]['dimension_bound'] == 63
    g = 2*s.eye(11)
    for i in range(10):
        g[i, i+1] = g[i+1, i] = -s.sqrt(2) if i in (0, 9) else -1
    nu = s.Matrix([1/s.sqrt(2)]+[1]*9+[1/s.sqrt(2)])
    assert g*nu == s.zeros(11, 1) and g.rank() == 10
    W = list(range(1, 10))
    J = g.extract(W, W).inv()
    va, vb = [-nu[i]*g.extract([i], W).T for i in (0, 10)]
    S = s.diag(nu[0], nu[10])*g.extract([0, 10], [0, 10])*s.diag(nu[0], nu[10])
    assert (S-s.Matrix.hstack(va, vb).T*J*s.Matrix.hstack(va, vb)).applyfunc(s.simplify) == s.Matrix([[1, -1], [-1, 1]])/10
    for mark in (0, 10):
        ids = [i for i in range(11) if i != mark]
        deleted = g.extract(ids, ids)
        assert all(deleted[:i, :i].det().is_positive for i in range(1, 11))
        inverse = deleted.inv()
        allowed = [(ids[i], s.simplify(inverse[i, i])) for i in range(10) if inverse[i, i] <= 2]
        assert allowed == ([(1, 1), (2, 2)] if mark == 0 else [(8, 2), (9, 1)])
    t, u = s.symbols('t u', nonnegative=True)
    profile_data = [
        ('isolated', {}, 0, 2, t/s.sqrt(2), 1, 4),
        ('end3', {1: 1}, 1, 1, t/s.sqrt(2)+1, 0, 3),
        ('end4', {1: s.sqrt(2)}, s.sqrt(2), 0, t/s.sqrt(2)+s.sqrt(2), 0, 2),
        ('next3', {2: 1}, 2, 0, t/s.sqrt(2)+1, 0, 4),
    ]
    profiles, native = [], {}
    for name, extra, moment, c0, height, minimum, maximum in profile_data:
        vector = s.Matrix([extra.get(i, 0) for i in W])
        assert s.simplify(nu[0]*t+sum(nu[i]*q for i, q in extra.items())-height) == 0
        assert s.simplify((vb.T*J*vector)[0]-moment/10) == 0
        assert s.simplify(2-(vector.T*J*vector)[0]-c0-moment**2/s.Integer(10)) == 0
        cost = s.factor(2*moment/height+c0/height**2)
        assert s.simplify(cost.subs(t, minimum)) == maximum
        # t>=1 for an isolated positive-height native opposite contact;
        # every other profile allows all t>=0. Positive polynomial
        # numerator and denominator prove the bound on the entire ray.
        difference = s.factor((maximum-cost).subs(t, u+minimum))
        numerator, denominator = s.fraction(difference)
        np = s.Poly(s.expand(numerator), u)
        dp = s.Poly(s.expand(denominator), u)
        assert all(c.is_nonnegative for c in np.all_coeffs())
        assert all(c.is_nonnegative for c in dp.all_coeffs()) and dp.eval(0).is_positive
        profiles.append(dict(profile=name, height=str(height), moment=str(moment), c0=str(c0),
                             contact_minimum=minimum, maximum_cost=maximum,
                             exact_cost=str(cost), upper_gap=str(difference),
                             upper_gap_numerator_coefficients=list(map(str, np.all_coeffs())),
                             upper_gap_denominator_coefficients=list(map(str, dp.all_coeffs()))))
        native[name] = vector
    pair_bounds = []
    for left, right in product(profiles, repeat=2):
        x, y = native[left['profile']], s.Matrix(list(reversed(native[right['profile']])))
        moment_product = s.sympify(left['moment'])*s.sympify(right['moment'])
        assert s.simplify((x.T*J*y)[0]-moment_product/10) == 0
        bound = 1-F(left['maximum_cost']+right['maximum_cost'], 10)
        assert bound >= F(1, 5)
        pair_bounds.append(dict(left=left['profile'], right=right['profile'], lower_bound=str(bound)))
    h, k, l, m, c0, d0, p, q, r = s.symbols('h k l m c0 d0 p q r', nonzero=True)
    a = s.Rational(1, 10)
    T = s.Matrix([[a, -a, -h+l*a, -m*a], [-a, a, -l*a, -k+m*a],
                  [-h+l*a, -l*a, c0+l*l*a-p, -l*m*a-q],
                  [-m*a, -k+m*a, -l*m*a-q, d0+m*m*a-r]])
    identity = 1-(2*l/h+c0/h**2+2*m/k+d0/k**2)/10+(p/h**2-2*q/(h*k)+r/k**2)/10
    assert s.factor(T.det()/(h*h*k*k)-identity) == 0
    weights = {tuple(row['pair']): F(row['weight']) for row in old['weighted_pairs']}
    assert weights[(0, 10)] == F(1, 2)
    weights[(0, 10)] = F()
    assert set(weights) == set(combinations(range(11), 2)) and sum(weights.values(), F()) == 54
    assert all(weights[i, j] == 1 for i, j in weights if j-i <= 7)
    assert min(F(row['lower_bound']) for row in pair_bounds) == F(1, 5)
    paths = [Path(__file__)]+[ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS
                             for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md',
              ROOT/'results/fresh_finite7_ideal5_triangle_payments_20261017.md']
    return dict(status='PASS_ORDER11_C_ENDPOINT_ZERO_TRANSFER',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
                target_gram=[[str(x) for x in row] for row in g.tolist()], positive_kernel=list(map(str, nu)),
                complete_profiles=profiles, all_ordered_profile_bounds=pair_bounds,
                determinant_identity=str(identity), strict_minimum='1/5',
                changed_pairs=[[0, 10]], old_weight='1/2', new_weight='0', old_cost='109/2', new_cost='54',
                weighted_pairs=[dict(pair=list(pair), weight=str(value)) for pair, value in sorted(weights.items())],
                new_order11_costs=dict(A='99/2', B='54', C='54', D='53'), new_order11_envelope='54/11',
                finite_coupon_increment='0', other_weights_unchanged=True,
                actual_alternative='An adjacent black contact or shared spectator coordinate gives a distinct actual corner a retained length1 or2 unit.',
                simultaneous_priority='First handle any changed C endpoint corner by this independent unit alternative. If none exists, all actual weights equal the completed baseline and its whole-face coverage is inherited.',
                triangle_minimum=2, quadrilateral_minimum=1, new_dimension_bound=None,
                scope='Universal local zero transfer using all native profiles and arbitrary supporting spectators. All16 profile-pair bounds are strict; affine alternative halves are never spent. No new finite budget or ambient dimension theorem.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = ROOT/f'results/{TAG}.json'
    out = json.loads(json.dumps(run()))
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], 'C cost', out['new_cost'], 'uniform gap', out['strict_minimum'], flush=True)
