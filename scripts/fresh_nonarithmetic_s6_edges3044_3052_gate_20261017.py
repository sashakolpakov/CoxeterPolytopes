#!/usr/bin/env python3
"""All-label obstruction to one cap terminating both S6 edges30/44,30/52.

Every rejection uses outward rational intervals for frozen algebraic recipes.
The two infinite angle tails are covered by continuous derivative bounds.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_h3_region_truncations_20261017 as R
import fresh_nonarithmetic_h3_edge730_domain_20261017 as D

E, P, U, I = R.E, R.P, R.U, R.I
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_s6_edges3044_3052_gate_20261017'
PATCH = dict(sides=[30, 44, 52], edges=[[30, 44], [30, 52]])


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def subdivide(lo, hi, predicate):
    pending, accepted = [(lo, hi)], []
    while pending:
        a, b = pending.pop()
        record = predicate(I(a, b))
        if record is not None:
            accepted.append(dict(interval=[str(a), str(b)], **record))
        else:
            assert b-a > F(1, 4096), ('unresolved continuous gate', a, b)
            mid = (a+b)/2
            pending.extend([(mid, b), (a, mid)])
    accepted.sort(key=lambda r: F(r['interval'][0]))
    assert F(accepted[0]['interval'][0]) == lo
    assert F(accepted[-1]['interval'][1]) == hi
    assert all(a['interval'][1] == b['interval'][0] for a, b in zip(accepted, accepted[1:]))
    return accepted


def run():
    endpoint = E.B.read(f'results/{E.TAG}.json')
    domains = {tuple(r['sides']): r for r in endpoint['endpoint_domains']}
    profiles = R.profile_intersection(PATCH, domains)
    low = [r for r in profiles if 6 not in r.values()]
    tails = [r for r in profiles if 6 in r.values()]
    assert len(low) == 38 and len(tails) == 2
    assert all(r[44] == 6 and r[30] == r[5] == r[6] == 2 for r in tails)
    assert sorted(r[52] for r in tails) == [2, 3]
    _, qi, normals, h, _, _, _, _ = P.physical_data()
    pair = lambda x, y: U.idot(x, U.imv(qi, y))
    ids = P.S + PATCH['sides']
    frame = [normals[i] for i in ids]
    gram = [[pair(x, y) for y in frame] for x in frame]
    inv, pivots = P.C53.inverse(gram)
    kernel = P.sub(h, P.combine(U.imv(inv, [pair(h, n) for n in frame]), frame))
    square = pair(kernel, kernel)
    assert square.l > 0
    k = [U.rnd(-x/U.sqrt(square)) for x in kernel]
    assert pair(k, h).h < 0
    # All candidates are p(g) +/- sqrt(2-p(g)^2)*k: rank9 frame, k^2=1.
    # These exact identities follow from projection, not from interval zero tests.
    decisions = []
    for index, labels in enumerate(low):
        contact = [-D.theta(labels.get(i, 2)) for i in ids]
        p = P.combine(U.imv(inv, contact), frame)
        residual = U.rnd(2-pair(p, p))
        assert residual.l > 0
        for sign in (-1, 1):
            r = [U.rnd(x+sign*U.sqrt(residual)*y) for x, y in zip(p, k)]
            gates = [(j, pair(r, normals[j])) for j in range(55) if j not in ids]
            positive = [(j, value) for j, value in gates if value.l > 0]
            if positive:
                wall, value = positive[0]
                result = dict(status='OBTUSE_OLD_CONTACT', wall=wall, product=value.enc())
            else:
                result = None
                for wall, value in gates:
                    if value.h < 0 and value.l > -2:
                        for m in range(2, 65):
                            if (value+D.theta(m)).h < 0 and (value+D.theta(m+1)).l > 0:
                                result = dict(status='CONSECUTIVE_LABEL_GAP', wall=wall,
                                              product=value.enc(), consecutive_labels=[m, m+1])
                                break
                    if result is not None:
                        break
                assert result is not None, (index, labels, sign)
            decisions.append(dict(profile=index, labels={str(i): m for i, m in labels.items()},
                                  branch=sign, residual=residual.enc(), **result))
    print('LOW_PROFILES', len(low), dict(Counter(r['status'] for r in decisions)), flush=True)

    tail_records = []
    for labels in tails:
        contact = [-D.theta(labels.get(i, 2)) if i != 44 else I(0) for i in ids]
        p0 = P.combine(U.imv(inv, contact), frame)
        p1 = P.combine(U.imv(inv, [I(-int(i == 44)) for i in ids]), frame)
        a, b, d = -pair(p1, p1), U.rnd(-2*pair(p0, p1)), U.rnd(2-pair(p0, p0))
        coefficients = {}
        for j in (4, 8):
            coefficients[j] = (pair(p1, normals[j]), pair(p0, normals[j]), pair(k, normals[j]))

        def rho(c):
            return U.rnd(a*c*c+b*c+d)

        def value(c, wall=8, sign=1):
            alpha, beta, gamma = coefficients[wall]
            return U.rnd(alpha*c+beta+sign*gamma*U.sqrt(rho(c)))

        def bound(c):
            residual = rho(c)
            if residual.l <= 0:
                return None
            alpha, beta, gamma = coefficients[8]
            derivative = U.rnd(alpha+gamma*(2*a*c+b)/(2*U.sqrt(residual)))
            other = value(c, wall=4, sign=-1)
            if derivative.h < 0 and other.l > 0:
                return dict(residual=residual.enc(), positive_branch_wall8_derivative=derivative.enc(),
                            negative_branch_wall4_product=other.enc())
            return None

        continuous = subdivide(D.theta(6).l, F(2), bound)
        endpoints = {str(m): value(D.theta(m)) for m in (6, 10, 11)}
        endpoints['infinity'] = value(I(2))
        intervals = []
        if labels[52] == 2:
            assert (endpoints['6']+D.theta(3)).h < 0
            assert (endpoints['infinity']+D.theta(4)).l > 0
            intervals.append(dict(labels_to44='all integers m>=6', gap=[3, 4]))
        else:
            assert (endpoints['6']+D.theta(4)).h < 0
            assert (endpoints['10']+D.theta(5)).l > 0
            assert (endpoints['11']+D.theta(5)).h < 0
            assert (endpoints['infinity']+D.theta(6)).l > 0
            intervals.extend([dict(labels_to44='6<=m<=10', gap=[4, 5]),
                              dict(labels_to44='all integers m>=11', gap=[5, 6])])
        tail_records.append(dict(labels={str(i): m for i, m in labels.items()},
                                 parameter='c=2cos(pi/m), m>=6', p0=P.enc(p0), p1=P.enc(p1),
                                 residual_polynomial=P.enc([a, b, d]),
                                 coefficients={str(j): P.enc(z) for j, z in coefficients.items()},
                                 continuous_domain='sqrt3<=c<=2', continuous_cover=continuous,
                                 negative_branch='strictly obtuse to old4 throughout',
                                 positive_branch_wall8_values={j: z.enc() for j, z in endpoints.items()},
                                 positive_branch_all_label_gaps=intervals))
        print('TAIL', labels[52], 'continuous intervals', len(continuous), flush=True)
    paths = [f'scripts/{TAG}.py', f'scripts/{R.TAG}.py', f'scripts/{D.TAG}.py',
             f'scripts/{E.TAG}.py', f'results/{E.TAG}.json', f'scripts/{P.TAG}.py',
             f'results/{P.TAG}.json', P.U.CAP, P.C53.CL, P.CAP52, P.CAPS55, P.CONE55,
             'results/fresh_nonarithmetic_nonright_s6_20261016_cone_review.json']
    return dict(status='PASS_ALL_LABEL_TWO_EDGE_COMMON_CAP_OBSTRUCTION',
                sources={p: sha(p) for p in paths}, retained_walls=55, face_support=P.S,
                target_edges=PATCH['edges'], full_frame=ids,
                frame_gram=[P.enc(row) for row in gram], inverse_pivots=pivots,
                projected_center_square=square.enc(), positive_kernel=P.enc(k),
                finite_profiles=len(low), finite_branch_decisions=decisions,
                tails=tail_records, decision_counts=dict(Counter(r['status'] for r in decisions)),
                exact_recipe='p=F*G^-1*g; k=-proj_ker(F)(h)/sqrt(proj_ker(F)(h)^2); r=p+/-sqrt(2-p^2)*k.',
                conclusion='No norm-two Coxeter wall acute to the retained55 can simultaneously give finite endpoints to S6+[30,44] and S6+[30,52]. All fields and all finite labels are covered. The conclusion also holds for every extension retaining these55walls.',
                scope='A necessary local obstruction. It excludes the proposed two-wall bridge around old44/52 whose second cap must terminate both selected edges; different multiwall incidence patterns remain open. No finite-volume example or global dimension claim.',
                proof='Each finite endpoint forces one of the frozen complete spherical profiles. Their intersection has38 finite profiles and2 isolated high-label tails. The independent nine-normal frame makes the norm-two solutions exactly two roots. Each low branch has a strictly obtuse old contact or lies strictly between consecutive Coxeter values. Each tail has a continuously obtuse negative branch; the positive branch wall8 product is strictly decreasing. Exact values at6,10,11 and the c=2 limit place every integer label in the displayed consecutive-label gaps. The continuous cover includes every c in[sqrt3,2], so no high-label cutoff is assumed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
        print('PASS_SAVED_ALL_LABEL_TWO_EDGE_COMMON_CAP_REPLAY')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], flush=True)
