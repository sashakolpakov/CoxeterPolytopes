#!/usr/bin/env python3
"""Exact all-rank scalar composition; written geometric proofs remain inputs."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_cusp_minimum_factor_rank17_dimension29_20261017'
PACKAGES = (
    'fresh_cusp_factor_payment_framework_20261017',
    'fresh_cusp_spectator_role_payments_20261017',
    'fresh_simple_cusp_vinberg_weight_frontier_20261017',
    'fresh_large_affine_factor_quad_payment_20261017',
    'fresh_cusp_link_shape_curvature_20261017',
    'fresh_cusp_link_rank_monotonicity_20261017',
    'fresh_private_affine_tree_exclusion_20261017',
    'fresh_affine_tree_carrier_exclusion_20261017',
)
WRITTEN = (
    'fresh_simple_finite_quad_vinberg_transfer_20261017',
    'fresh_private_affine_cycle_exclusion_20261017',
    'fresh_affine_cycle_carrier_exclusion_20261017',
)


def digest(p):
    h = sha256()
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def beta(n):
    return F(n*(n-6),8) if n % 2 == 0 else F((n-1)*(n-5),8)


def run():
    paths = [Path(__file__)]
    input_statuses = {}
    for tag in PACKAGES:
        p = ROOT / f'results/{tag}.json'
        data = json.loads(p.read_text())
        assert data['status'].startswith('PASS'), (tag,data['status'])
        input_statuses[tag] = data['status']
        for name, expected in data.get('sources',{}).items():
            assert isinstance(expected,str) and digest(ROOT/name) == expected, (tag,name)
        paths += [p,ROOT/f'results/{tag}.md',ROOT/f'scripts/{tag}.py']
    paths += [ROOT/f'results/{tag}.md' for tag in WRITTEN]
    assert all(p.is_file() for p in paths)

    t,s,m = sp.symbols('t s m', nonnegative=True)
    even_n,odd_n = 30+2*t,31+2*t
    even = sp.expand(even_n*(even_n-6)-24*even_n)
    odd = sp.expand((odd_n-1)*(odd_n-5)-24*odd_n)
    assert even == 4*t*t+60*t
    assert odd == 4*t*t+64*t+36
    assert all(c>=0 for c in sp.Poly(even,t).all_coeffs())
    assert all(c>=0 for c in sp.Poly(odd,t).all_coeffs())
    assert sp.Poly(odd,t).TC() > 0

    # Odd r=17+2t and even r=18+2t have the same pointwise cap.
    odd_r,even_r = 17+2*t,18+2*t
    cap_odd=(odd_r+3)/(4*(odd_r+2))
    cap_even=(even_r+2)/(4*(even_r+1))
    assert sp.simplify(cap_odd-cap_even)==0
    cap_slack=sp.factor(sp.Rational(5,19)-cap_odd)
    assert sp.simplify(cap_slack-t/(38*(19+2*t)))==0
    increment=F(3,4)-F(5,2)*F(5,19)
    assert increment==F(7,76)>0

    # Compare 39s/4 with (5/2)*sqrt(323s/12); both are positive.
    squared_gap=sp.expand(48*((sp.Rational(39,4)*s)**2-sp.Rational(25,4)*323*s/12))
    assert sp.expand(squared_gap-s*(4563*s-8075))==0
    assert sp.expand((4563*s-8075).subs(s,2+t))==1051+4563*t

    affine_costs={
        'cycle':(3*m,0),
        'path':(3*m-6,6),
        'single_fork':(3*m-5,5),
        'double_fork':(3*m-4,4),
    }
    for near,end in affine_costs.values():
        assert sp.expand(near+end)==3*m
    checks=[dict(n=n,beta=str(beta(n)),finite_slack=str(beta(n)-3*n),
                 extremum_slack=5*(n-1)) for n in range(30,201)]
    assert all(F(c['finite_slack'])>=0 and c['extremum_slack']>0 for c in checks)
    assert beta(30)==90

    return dict(
        status='PASS_ALL_RANK_MINIMUM_CUSP_FACTOR_RANK17_DIMENSION29_COMPOSITION',
        sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
        input_statuses=input_statuses,
        hypothesis='Every irreducible affine cusp factor has rank at least17.',
        dimension_bound=29,arithmetic_hypothesis=False,factor_count_limit=None,
        affine_order_exclusion_through=17,
        finite_cost_coefficient=3,ideal_cost_per_factor_per_link_vertex=3,
        even_finite_slack_times8=str(even),odd_finite_slack_times8=str(odd),
        pointwise_cap_slack=str(cap_slack),curvature_increment=str(increment),
        multiple_factor_squared_gap_times48=str(squared_gap),
        multiple_factor_shifted_positive_coefficient='1051+4563*t',
        affine_cost_identities={k:[str(a),str(b),str(sp.expand(a+b))] for k,(a,b) in affine_costs.items()},
        scalar_spot_checks=checks,
        all_rank_argument='Nonnegative parity polynomials; uniform centered-uniform increment; strict multiple-factor squared comparison; positive global5(n-1)slack.',
        geometric_inputs='Actual factor fibers, triangle transfer, large ideal-corner payments, cycle/tree complement exclusions, and the audited published finite-corner proof are written proof inputs, not formalized here.',
        enumeration_cutoff_used=False,unrestricted163_tail_used=False,new_example=None,
        scope='Exact scalar composition and direct dependency hashes. Does not independently rederive published diagram tables or formally verify the written geometric lemmas.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['dimension_bound'],len(result['sources']),flush=True)
