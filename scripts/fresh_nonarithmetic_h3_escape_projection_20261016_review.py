#!/usr/bin/env python3
"""Independent finite-vertex projection obstruction using physical coordinates."""
import argparse
import json
from fractions import Fraction as F
import fresh_nonarithmetic_escape52_20261016_cap_review as C

B, ROOT, I = C.B, C.ROOT, C.I
TAG = 'fresh_nonarithmetic_h3_escape_projection_20261016_review'
PRIMARY = 'results/fresh_nonarithmetic_h3_escape_projection_20261016.json'
CONE = 'results/fresh_nonarithmetic_escape52_20261016_cone.json'
CONE_REVIEW = 'results/fresh_nonarithmetic_escape52_20261016_cone_review.json'
S = [0, 1, 2, 3, 5, 6]
SUPPORTS = [S + [7, 8, 11], S + [7, 10, 33], S + [8, 10, 44], S + [8, 12, 34]]


def pd(matrix):
    a = [row[:] for row in matrix]
    pivots = []
    for k in range(len(a)):
        p = a[k][k]
        assert p.lo > 0
        pivots.append(p.saved())
        for i in range(k + 1, len(a)):
            for j in range(i, len(a)):
                a[j][i] = a[i][j] = a[i][j] - a[i][k] * a[k][j] / p
    return pivots


def overlap(value, raw):
    lo, hi = map(F, raw)
    assert value.lo <= hi and lo <= value.hi


def run():
    target, cone, reviewed = map(B.read, (PRIMARY, CONE, CONE_REVIEW))
    assert target['status'] == 'PASS_EXACT_RENEWED53_ESCAPE_POSITIVE_FINITE_VERTEX_PROJECTION'
    assert reviewed['status'] == 'PASS_INDEPENDENT_COMPLETE_53_WALL_CONE_AND_RESIDUAL_H3_FACE'
    assert target['right_support'] == S
    g = C.geometry()
    N, cov, H, hc, Q = g['NI'], g['covI'], g['HI'], g['hI'], g['QI']
    pair = lambda a, b: B.idot(a, B.imv(Q, b))
    E = [B.iv(B.dec(x)) for x in cone['strict_escape']['point']]
    assert target['escape'] == cone['strict_escape']['point']
    assert pair(E, E).lo > 0 and B.idot(E, hc).hi < 0
    assert all(B.idot(E, v).hi < 0 for v in cov)
    SG = [[B.idot(N[i], cov[j]) for j in S] for i in S]
    support_pd = pd(SG)
    SJ, support_inverse_pivots = C.inverse_interval(SG)
    coefficients = B.imv(SJ, [B.idot(E, cov[j]) for j in S])
    projection = [x - y for x, y in zip(E, C.combine(coefficients, [N[j] for j in S]))]
    L = -B.idot(projection, hc)
    assert L.lo > 0
    projection = [v / L for v in projection]
    assert pair(projection, projection).hi < 0
    V, corners = [], []
    for support in SUPPORTS:
        assert support in reviewed['complete_census']['finite_vertices']
        G = [[B.idot(N[i], cov[j]) for j in support] for i in support]
        gp = pd(G)
        J, pivots = C.inverse_interval(G)
        coefficients = B.imv(J, [B.idot(H, cov[j]) for j in support])
        point = [x - y for x, y in zip(H, C.combine(coefficients, [N[j] for j in support]))]
        norm = pair(point, point)
        scale = -B.idot(point, hc)
        assert norm.hi < 0 and scale.lo > 0
        point = [v / scale for v in point]
        products = [B.idot(point, v) for v in cov]
        for j, v in enumerate(products):
            if j in support:
                assert v.lo <= 0 <= v.hi
            else:
                assert v.hi < 0, (support, j)
        V.append(point)
        corners.append(dict(support=support,positive_Gram_pivots=gp,
            inverse_pivots=pivots,unscaled_norm=norm.saved(),normalization=scale.saved(),
            point=[v.saved() for v in point],all53_pairs=[v.saved() for v in products]))
    # Use four literal coordinate functionals selected by complete
    # pivoting, rather than the producer's four old-root covectors.
    minor = B.R.rank_minor(V)
    assert minor['rank_lower_bound'] == 4
    indices = minor['original_columns']
    M = [[v[i] for v in V] for i in indices]
    inverse, inverse_pivots = C.inverse_interval(M)
    weights = B.imv(inverse, [projection[i] for i in indices])
    assert all(v.lo > 0 for v in weights)
    for value, raw in zip(weights, target['positive_weights']):
        overlap(value, raw)
    # All four corner vectors and projection are exactly in S-perp by
    # their Gram-projection definitions. The nonsingular coordinate
    # minor proves they span this4space; solving those4equations thus
    # proves the full identity. Remaining overlaps are sanity checks.
    reconstructed = C.combine(weights, V)
    residual = [x - y for x, y in zip(projection, reconstructed)]
    assert all(v.lo <= 0 <= v.hi for v in residual)
    total = sum(weights, I())
    assert total.lo <= 1 <= total.hi
    paths = [f'scripts/{TAG}.py', PRIMARY, CONE, CONE_REVIEW,
             f'scripts/{C.TAG}.py', f'results/{C.TAG}.json']
    return dict(status='PASS_INDEPENDENT_ALL_REAL_H3_RIGHT_SUPPORT_ESCAPE_OBSTRUCTION',
        sources={p:B.digest(p) for p in paths},right_support=S,
        support_positive_pivots=support_pd,support_inverse_pivots=support_inverse_pivots,
        positive_projection_scale=L.saved(),projected_escape=[v.saved() for v in projection],
        projected_norm=pair(projection,projection).saved(),vertices=corners,
        physical_coordinate_minor=minor,coordinate_matrix=[[v.saved() for v in row] for row in M],
        coordinate_inverse_pivots=inverse_pivots,positive_weights=[v.saved() for v in weights],
        weight_sum=total.saved(),identity_residual_intervals=[v.saved() for v in residual],
        strict_nonincident_vertex_checks=176,
        persistence='For anycommonfuturestrictinteriorX, eachfinitecorner v has X=a*v+sum b_i*n_i, witha>0 andb_i<=0 byinversepositivityofitsPDZGram. Acuteness gives X.r>=a*v.r, so v.r<0. No oldcenterretention isneeded.',
        proof='Independent53physicalgeometry; fourliteralPD9cornerprojections withallstrictsigns; independentcoordinate4minor establishesbasis ofSperp. Positive solution givesprojectedE asstrictconvexcombination. Finitevertexpersistence thenforces everycoherentacuteSrightwall tomissEstrictly.',
        scope='Allrealcoefficientsandlabels, forspacelikeaddednormalsrighttoS andadmittingacommonstrictfutureinteriorwiththe53oldwalls. Anyfinitefamilyofsuchwallsretainsstrictescape. ArbitrarynonrightScompletionsremainopen; noclaimforincoherenthalfspaces.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
