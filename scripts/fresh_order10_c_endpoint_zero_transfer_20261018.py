#!/usr/bin/env python3
"""Delete the order10 C endpoint unit using every native heavy-end profile.

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
TAG = 'fresh_order10_c_endpoint_zero_transfer_20261018'
PREVIOUS = 'fresh_order11_c_endpoint_zero_transfer_20261018'
BASE = 'fresh_a15_d15_dimension58_20261018'
FRONT = 'fresh_ambient58_payment_frontier_20261018'
ORIGIN = 'fresh_deep_finite8_dimension64_20261018'
STEMS = (PREVIOUS, BASE, FRONT, ORIGIN)
M, R = 10, 9



def run():
    docs = {tag: json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag, packet in docs.items():
        assert packet['status'].startswith('PASS'), tag
        for path, digest in packet['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, (tag, path)
    assert docs[PREVIOUS]['strict_minimum'] == '1/5'
    assert docs[BASE]['dimension_bound'] == 58
    assert docs[FRONT]['ideal']['current_envelopes']['9'] == '9/2'
    assert docs[ORIGIN]['ideal_weights_unchanged'] is True
    # The original all-pairs-unit rule at orders <=11 is inherited.
    # Every subsequent changed ideal role has order >=11; the completed
    # rank9 envelope is still 45/10. The report records that interface.
    g = 2*s.eye(M)
    for i in range(M-1):
        g[i, i+1] = g[i+1, i] = -s.sqrt(2) if i in (0, M-2) else -1
    nu = s.Matrix([1/s.sqrt(2)]+[1]*(M-2)+[1/s.sqrt(2)])
    assert g*nu == s.zeros(M, 1) and g.rank() == R
    W = list(range(1, M-1))
    J = g.extract(W, W).inv()
    va, vb = [-nu[i]*g.extract([i], W).T for i in (0, M-1)]
    S = s.diag(nu[0], nu[M-1])*g.extract([0, M-1], [0, M-1])*s.diag(nu[0], nu[M-1])
    assert (S-s.Matrix.hstack(va, vb).T*J*s.Matrix.hstack(va, vb)).applyfunc(s.simplify) == s.Matrix([[1, -1], [-1, 1]])/R
    inverse_records = []
    for mark in (0, M-1):
        ids = [i for i in range(M) if i != mark]
        deleted = g.extract(ids, ids)
        determinants = [s.factor(deleted[:i, :i].det()) for i in range(1, M)]
        assert all(value.is_positive for value in determinants)
        inverse = deleted.inv()
        assert all(value.is_nonnegative for value in inverse)
        diagonal = [s.simplify(inverse[i, i]) for i in range(R)]
        allowed = [(ids[i], diagonal[i]) for i in range(R) if diagonal[i] <= 2]
        assert allowed == ([(1, 1), (2, 2)] if mark == 0 else [(M-3, 2), (M-2, 1)])
        positions = [ids.index(i) for i, _ in allowed]
        two_contact_projection = s.simplify(sum(inverse[i, j] for i in positions for j in positions))
        assert two_contact_projection >= 3 > 2
        inverse_records.append(dict(deleted_mark=mark, old_nodes=ids,
            leading_determinants=list(map(str, determinants)), inverse_diagonal=list(map(str, diagonal)),
            admissible_contacts=[dict(node=i, inverse_diagonal=str(v)) for i, v in allowed],
            two_contact_projection_minimum=str(two_contact_projection)))
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
        assert s.simplify((vb.T*J*vector)[0]-moment/R) == 0
        assert s.simplify(2-(vector.T*J*vector)[0]-c0-moment**2/s.Integer(R)) == 0
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
        assert s.simplify((x.T*J*y)[0]-moment_product/R) == 0
        bound = 1-F(left['maximum_cost']+right['maximum_cost'], R)
        assert bound >= F(1, 9)
        pair_bounds.append(dict(left=left['profile'], right=right['profile'], lower_bound=str(bound)))
    h, k, l, m, c0, d0, p, q, r = s.symbols('h k l m c0 d0 p q r', nonzero=True)
    a = s.Rational(1, R)
    T = s.Matrix([[a, -a, -h+l*a, -m*a], [-a, a, -l*a, -k+m*a],
                  [-h+l*a, -l*a, c0+l*l*a-p, -l*m*a-q],
                  [-m*a, -k+m*a, -l*m*a-q, d0+m*m*a-r]])
    identity = 1-(2*l/h+c0/h**2+2*m/k+d0/k**2)/R+(p/h**2-2*q/(h*k)+r/k**2)/R
    assert s.factor(T.det()/(h*h*k*k)-identity) == 0
    old_weights = {pair: F(1) for pair in combinations(range(M), 2)}
    weights = dict(old_weights)
    weights[(0, M-1)] = F()
    assert len(weights) == 45 and sum(old_weights.values(), F()) == 45
    assert sum(weights.values(), F()) == 44
    assert all(weights[i, j] == 1 for i, j in weights if j-i <= 8)
    assert all(weights[i, j] == weights[M-1-j, M-1-i] for i, j in weights)
    assert min(F(row['lower_bound']) for row in pair_bounds) == F(1, 9)
    paths = [Path(__file__)]+[ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS
                             for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md',
              ROOT/'results/fresh_finite7_ideal5_triangle_payments_20261017.md']
    return dict(status='PASS_ORDER10_C_ENDPOINT_ZERO_TRANSFER',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
                target_gram=[[str(x) for x in row] for row in g.tolist()], positive_kernel=list(map(str, nu)),
                deleted_principal_inverse_checks=inverse_records,
                complete_profiles=profiles, all_ordered_profile_bounds=pair_bounds,
                determinant_identity=str(identity), strict_minimum='1/9',
                changed_pairs=[[0, 9]], old_weight='1', new_weight='0', old_cost='45', new_cost='44',
                old_weighted_pairs=[dict(pair=list(pair), weight=str(value)) for pair, value in sorted(old_weights.items())],
                weighted_pairs=[dict(pair=list(pair), weight=str(value)) for pair, value in sorted(weights.items())],
                new_order10_costs=dict(A='45', B='45', C='44', D='45'),
                order10_envelope_unchanged='9/2', finite_coupon_increment='0',
                new_cross_component_coupons=[], other_weights_unchanged=True,
                all_short8_weights_unchanged=True, reflection_invariant=True,
                actual_alternative='An adjacent black contact or shared spectator coordinate gives a distinct actual corner a retained length1 or2 unit.',
                simultaneous_priority='Handle any actual changed C10 endpoint corner first with its independent short unit. If absent, all actual payments equal the completed58 rule.',
                triangle_minimum=2, quadrilateral_minimum=1, current_dimension_bound=58, new_dimension_bound=None,
                low9_envelope_unchanged=True, new_finite_budget=False,
                scope='Universal local zero transfer using complete native profiles and arbitrary supporting spectators. All16 profile-pair bounds are strict. Only the C family improves; the all-family rank9 envelope, finite costs and completed dimension bound remain unchanged.')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = ROOT/f'results/{TAG}.json'
    out = json.loads(json.dumps(run()))
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        assert not path.exists()
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], 'C cost', out['new_cost'], 'uniform gap', out['strict_minimum'], flush=True)
