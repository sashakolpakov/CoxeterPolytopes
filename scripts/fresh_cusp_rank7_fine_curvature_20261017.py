#!/usr/bin/env python3
"""All-rank fine-budget curvature: monotonicity and 21 minimal profiles."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_cusp_rank7_fine_curvature_20261017'
INPUTS = (
    'results/fresh_cusp_link_shape_curvature_20261017.md',
    'results/fresh_cusp_link_shape_curvature_20261017.json',
    'scripts/fresh_cusp_link_shape_curvature_20261017.py',
    'results/fresh_cusp_link_rank_monotonicity_20261017.md',
    'results/fresh_cusp_link_rank_monotonicity_20261017.json',
    'scripts/fresh_cusp_link_rank_monotonicity_20261017.py',
    'results/fresh_cusp_minimum_factor_rank17_dimension29_20261017_review.json',
    'scripts/fresh_cusp_minimum_factor_rank17_dimension29_20261017_review.py',
    'results/fresh_cusp_spectator_role_payments_20261017.md',
)


def digest(path):
    h = sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def partitions(total, count, minimum):
    if count == 1:
        if total >= minimum:
            yield (total,)
    else:
        for first in range(minimum, total // count + 1):
            for tail in partitions(total-first, count-1, first):
                yield (first,) + tail


def convolution(ranks):
    h = [1]
    for r in ranks:
        nxt = [0] * (len(h)+r)
        for i, a in enumerate(h):
            for j in range(r+1):
                nxt[i+j] += a
        h = nxt
    return h


def inclusion_exclusion(ranks):
    """Independent bounded-composition coefficient formula."""
    s = len(ranks)
    result = []
    for j in range(sum(ranks)+1):
        coefficient = 0
        for size in range(s+1):
            for subset in combinations(range(s), size):
                remaining = j-sum(ranks[i]+1 for i in subset)
                if remaining >= 0:
                    coefficient += (-1)**size * comb(remaining+s-1, s-1)
        result.append(coefficient)
    return result


def profile(ranks):
    h = convolution(ranks)
    assert h == inclusion_exclusion(ranks)
    d = sum(ranks)
    v = prod(r+1 for r in ranks)
    assert h == h[::-1] and sum(h) == v
    g = sum(a*min(j,d-j) for j,a in enumerate(h))
    kappa_over_v = F(5*g-d*v, 2*v)
    deviation = sum((F(a)*abs(F(j)-F(d,2)) for j,a in enumerate(h)), F())/v
    assert kappa_over_v == F(3*d,4)-F(5,2)*deviation
    return dict(ranks=list(ranks), dimension=d+1, link_vertices=v,
                h_coefficients=h, G=g, mean_absolute_deviation=str(deviation),
                normalized_curvature=str(kappa_over_v),
                margin_over_3s=str(kappa_over_v-3*len(ranks)))


def run():
    for name in INPUTS:
        if name.endswith('.json'):
            data = json.loads((ROOT/name).read_text())
            assert data['status'].startswith('PASS')
            for source, expected in data.get('sources', {}).items():
                assert isinstance(expected,str) and digest(ROOT/source)==expected

    # For odd r=7+2t or even r=8+2t the cap is (10+2t)/(4(9+2t)).
    # Clearing denominators of 5/18 minus that cap leaves 4t.
    assert [5*4*9-18*10, 5*4*2-18*2] == [0,4]
    increment = F(3,4)-F(5,2)*F(5,18)
    assert increment == F(1,18)>0
    # At s>=7, 16[(9s/4)^2-(25/4)(21s/4)] = s(81s-525).
    assert 16*F(9,4)**2==81
    assert 16*F(25,4)*F(21,4)==525
    assert 81*7-525==42>0
    # Independent exhaustive generation of the minimal rank profiles.
    rows=[]
    for s in range(2,7):
        total=max(29,7*s)
        group=list(partitions(total,s,7))
        # Shift by7 and use a second unrestricted partition recursion.
        def increments(n,k,upper):
            if k==0:
                if n==0:yield ()
            else:
                for last in range(min(n,upper)+1):
                    for rest in increments(n-last,k-1,last):
                        yield rest+(last,)
        alternate={tuple(7+x for x in p) for p in increments(total-7*s,s,total-7*s)}
        assert set(group)==alternate
        for ranks in group:
            row=profile(ranks)
            assert F(row['margin_over_3s'])>0
            rows.append(row)
    assert len(rows)==21
    minimum=min(F(r['margin_over_3s']) for r in rows)
    assert minimum==F(41,256)
    assert [r['ranks'] for r in rows if F(r['margin_over_3s'])==minimum]==[[7,7,7,8]]
    # Single-factor and finite-vertex coefficients, coefficientwise for t>=0.
    assert [30*24-24*30,2*30+2*24-48,4]==[0,60,4]
    assert [30*26-24*31,2*30+2*26-48,4]==[36,64,4]
    obstruction=profile((6,)*5)
    assert F(obstruction['margin_over_3s'])==F(-49795,33614)<0
    paths=[Path(__file__)]+[ROOT/p for p in INPUTS]
    return dict(status='PASS_ALL_RANK_MINIMUM_RANK7_FINE_CURVATURE_COMPARISON',
                sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
                minimum_factor_rank=7,ambient_dimension_at_least=30,
                pointwise_increment_cap='5/18',pointwise_cap_slack='t/(18*(9+2*t))',
                curvature_increment=str(increment),
                many_factor_threshold=7,many_factor_squared_margin_times16='s*(81*s-525)',
                many_factor_shifted_positive_factor=[42,81],
                minimal_profiles=rows,minimum_base_margin=str(minimum),
                coefficients_8_beta_minus_3n_even=[0,60,4],
                coefficients_8_beta_minus_3n_odd=[36,64,4],
                rank6_counterexample=obstruction,
                finite_rank_cutoff_used=False,unrestricted163_tail_used=False,
                scope='All-rank curvature inequality only. A dimension29 theorem still requires a complete geometric payment rule under the proposed cusp hypothesis. The rank6 counterexample is a legitimate local cusp-link shape, not an ambient finite-volume counterexample.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(result['minimal_profiles']),'minimal profiles',flush=True)
