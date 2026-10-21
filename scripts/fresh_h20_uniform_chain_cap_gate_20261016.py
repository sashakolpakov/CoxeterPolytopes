#!/usr/bin/env python3
"""All-label obstruction to a cap having one finite angle to all26 chain walls.

The finite support-cost census uses standalone Fraction multiquadratic
arithmetic. Physical frame signs reuse the frozen inherited field helper.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path

from fresh_h20_coupled_roles_20261014_field import Field, parse
from fresh_h20_s61_multicap_20261014_exact import I, FILES, BITS

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_uniform_chain_cap_gate_20261016'
PRIMARY = 'results/fresh_h20_s61_multicap_20261014_exact.json'
REVIEW = 'results/fresh_h20_s61_multicap_20261014_exact_review.json'
RAD = (1, 2, 3, 6, 5, 10, 15, 30)
ZERO = (F(0),)*8


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def scalar(x):
    return (F(x),)+(F(0),)*7


def add(x, y):
    return tuple(a+b for a, b in zip(x, y))


def scale(x, t):
    return tuple(a*t for a in x)


def mul(x, y):
    out = [F(0)]*8
    for i, a in enumerate(x):
        if not a:
            continue
        for j, b in enumerate(y):
            if b:
                out[i ^ j] += a*b*RAD[i & j]
    return tuple(out)


DEN = 1 << 160
ROOT_BOUNDS = []
for d in RAD:
    r = isqrt(d*DEN*DEN)
    ROOT_BOUNDS.append((F(r, DEN), F(r + (r*r != d*DEN*DEN), DEN)))


def sign(x):
    if x == ZERO:
        return 0
    lo = sum(a*(l if a >= 0 else h) for a, (l, h) in zip(x, ROOT_BOUNDS))
    hi = sum(a*(h if a >= 0 else l) for a, (l, h) in zip(x, ROOT_BOUNDS))
    assert lo > 0 or hi < 0, x
    return 1 if lo > 0 else -1


def cost_census(J, S):
    one = scalar(1)
    sqrt2 = tuple(F(i == 1) for i in range(8))
    sqrt3 = tuple(F(i == 2) for i in range(8))
    phi = tuple(F(1, 2) if i in (0, 4) else F(0) for i in range(8))
    values = {3: one, 4: sqrt2, 5: phi, 6: sqrt3}
    pair_values = {(i, j): mul(values[i], values[j]) for i in values for j in values}
    assert mul(sqrt2, sqrt2) == scalar(2)
    assert mul(sqrt3, sqrt3) == scalar(3)
    assert mul(phi, phi) == add(phi, one)
    counts = Counter()
    witnesses = []
    digest = sha256()
    total = 0
    for k in (1, 2, 3):
        for ids in combinations(range(18), k):
            for labels in product((3, 4, 5, 6), repeat=k):
                q = ZERO
                for a, i in enumerate(ids):
                    for b, j in enumerate(ids):
                        q = add(q, scale(pair_values[labels[a], labels[b]], J[i][j]))
                side = sign(add(q, scalar(-2)))
                counts[{-1: 'below_two', 0: 'equal_two', 1: 'above_gap'}[side]] += 1
                if side > 0:
                    compare = sign(add(q, scalar(-F(137, 66))))
                    assert compare >= 0
                    if not compare:
                        witnesses.append(dict(walls=[S[i] for i in ids], labels=list(labels)))
                digest.update(json.dumps([list(ids), list(labels), list(map(str, q))], separators=(',', ':')).encode())
                total += 1
    assert total == 54744 and len(witnesses) == 8
    assert counts == {'below_two': 97, 'equal_two': 8, 'above_gap': 54639}
    return dict(profiles=total, labels=[3, 4, 5, 6], maximum_nonzero_contacts=3,
                counts=dict(counts), strict_above_two_minimum='137/66',
                minimizing_profiles=witnesses, ordered_profile_digest=digest.hexdigest())


def interval(x):
    return [str(x.l), str(x.h)]


def run():
    docs = [json.loads((ROOT / p).read_text()) for p in FILES]
    old, bridge, edge, polar, cur = docs
    base = json.loads((ROOT / PRIMARY).read_text())
    review = json.loads((ROOT / REVIEW).read_text())
    for doc in docs+[base, review]:
        for p, h in doc['sources'].items():
            assert sha(p) == h, p
    assert review['sources'][PRIMARY] == sha(PRIMARY)
    S = base['support']
    G = [[F(x) for x in row] for row in base['support_gram']]
    J = [[F(x) for x in row] for row in base['support_inverse']]
    assert len(S) == 18
    assert all(sum(G[i][k]*J[k][j] for k in range(18)) == int(i == j)
               for i in range(18) for j in range(18))
    assert all(x >= 0 for row in J for x in row)
    diagonals = [J[i][i] for i in range(18)]
    assert min(diagonals) == F(1, 2)
    assert [S[i] for i, q in enumerate(diagonals) if q == F(1, 2)] == [17, 31]
    assert min(q for q in diagonals if q > F(1, 2)) == F(2, 3)
    assert sum(sorted(diagonals)[:4]) == F(7, 3)
    for wall in (17, 31):
        i = S.index(wall)
        assert all(J[i][j] == (F(1, 2) if i == j else 0) for j in range(18))

    # z=2cos(2pi/7)>1 is the root of z^3+z^2-2z-1.  The polynomial
    # is increasing for z>=1 and is negative at31/25.  Hence
    # (2cos(pi/7))^2=2+z>81/25.
    z = F(31, 25)
    cubic_value = z**3+z*z-2*z-1
    assert cubic_value < 0
    ceiling = F(141, 70)
    assert F(7, 3) > ceiling
    assert F(2, 3)*F(81, 25) > ceiling
    assert F(1, 2)*F(81, 25)+F(1, 2) > ceiling
    assert F(8, 3) > ceiling and F(5, 2) > ceiling
    assert F(137, 66) > ceiling
    census = cost_census(J, S)
    print('PASS complete finite-support cost gap', census['profiles'], census['counts'], flush=True)

    A = Field(old, bridge, edge)
    M = [tuple(parse(x, 1) for x in v) for v in old['normals']]
    assert all(F(old['norms'][i]) == 2 for i in S)
    assert [[A.rational(A.dot(M[i], M[j], 1)) for j in S] for i in S] == G
    X = tuple(parse(x, 2) for x in base['frame']['X'])
    qX = parse(base['frame']['A'], 2)
    assert A.dot(X, X, 2) == A.scale(qX, -1, 2)
    assert all(A.dot(X, A.vlift(M[i], 1, 2), 2) == A.z[2] for i in S)
    alo, ahi = map(F, base['parameter']['isolating_interval'])
    assert alo**2 > 280 and ahi > alo > 0
    # The interval includes the closed endpoints, so covers every m>=3.
    time = I(1/ahi, 2/alo)
    transverse = (4+2*time*time).sqrt()
    assert 2+4/(alo*alo) < ceiling
    blocker_rows = []
    for support_wall, blocker in ((17, 83), (31, 11)):
        qold = F(old['norms'][blocker])
        normden = I(qold/2).sqrt()
        xp = I(*A.bounds(A.dot(X, A.vlift(M[blocker], 1, 2), 2), 2, BITS))
        xp = xp/(I(*A.bounds(qX, 2, BITS)).sqrt()*normden)
        wp = I(*A.bounds(A.dot(M[support_wall], M[blocker], 1), 1, BITS))/(2*normden)
        physical = time*xp-transverse*wp
        assert physical.l > 0
        blocker_rows.append(dict(single_nonfinite_support_wall=support_wall,
                                 old_blocker=blocker, unit_time_pairing=interval(xp),
                                 dual_weight_pairing=interval(wp),
                                 whole_parameter_positive_pairing=interval(physical)))
    # Three consecutive spatial normals on the saved rotation circle are
    # noncollinear: determinant=2 sin(delta)(1-cos(delta)) times a positive
    # radius factor.  Therefore equal pairings force spatial residual zero.
    aa = I(alo, ahi)
    cosine = aa*aa/(aa*aa+2)
    sine = 2*(aa*aa+1).sqrt()/(aa*aa+2)
    assert cosine.l > 0 and cosine.h < 1 and sine.l > 0
    noncollinear = 2*sine*(1-cosine)
    assert noncollinear.l > 0
    paths = [f'scripts/{TAG}.py', PRIMARY, REVIEW] + FILES + [
        'scripts/fresh_h20_coupled_roles_20261014_field.py',
        'scripts/fresh_h20_s61_multicap_20261014_exact.py']
    return dict(status='PASS_ALL_LABEL_UNIFORM_FINITE_CHAIN_ANGLE_CAP_OBSTRUCTION',
                sources={p: sha(p) for p in paths},
                support=S, support_inverse=base['support_inverse'],
                finite_support_cost_census=census,
                all_label_reduction=dict(cost_ceiling=str(ceiling),
                    label7_squared_contact_lower_bound='81/25',
                    cubic_at_31_over25=str(cubic_value),
                    four_contact_lower_bound='7/3',
                    remaining_single_nonfinite_contacts=[17, 31]),
                chain_parameter_interval=[str(alo), str(ahi)],
                time_interval=interval(time), transverse_contact_interval=interval(transverse),
                spatial_noncollinearity_bound=interval(noncollinear),
                exact_blockers=blocker_rows,
                scope='Every real spacelike additional normal, positively normalized to square2, forming Coxeter pairs with all133 retained walls and meeting all26 saved new chain walls at one common finite Coxeter angle. All18 retained support contacts may be right, finite Coxeter, parallel or dotted. Coherence requires a common strict interior point. Every finite label is covered by analytic reduction. This does not exclude unequal chain angles, caps dotted to all26, changed retained walls or arbitrary H20 completions.',
                no_new_wall_constructed=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    output = ROOT / 'results' / f'{TAG}.json'
    if args.verify:
        assert result == json.loads(output.read_text())
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
