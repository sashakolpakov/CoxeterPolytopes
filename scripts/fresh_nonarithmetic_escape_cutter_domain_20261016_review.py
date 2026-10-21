#!/usr/bin/env python3
"""Independent all-real cutter gate from the reviewed complete 52-wall cone.

The extra inequality is r.E >= 0. Every crossing is certified by an
independent rank-eight minor. Strict future timelike projections of every
generator imply the same property for every nonzero conic combination.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
import fresh_nonarithmetic_affine_filtered_20261016_cluster_review as C

B, ROOT, I = C.B, C.ROOT, C.I
TAG = 'fresh_nonarithmetic_escape_cutter_domain_20261016_review'
PRIMARY = 'results/fresh_nonarithmetic_escape_cutter_domain_20261016.json'
CONE = 'results/fresh_nonarithmetic_affine_filtered_20261016_cone_review.json'
LOCAL = f'results/{C.TAG}.json'
ESCAPE = 'results/fresh_nonarithmetic_affine_filtered_20261016_cone.json'


def run():
    target, cone, local, source = map(B.read, (PRIMARY, CONE, LOCAL, ESCAPE))
    assert cone['status'] == 'PASS_INDEPENDENT_COMPLETE_FILTERED52_CONE_AND_ESCAPE'
    assert cone['sources'][ESCAPE] == B.digest(ESCAPE)
    assert cone['sources'][LOCAL] == B.digest(LOCAL)
    assert local['sources'][C.PRIMARY] == B.digest(C.PRIMARY)
    g = C.geometry()
    cov, Q, h = g['covI'], g['QI'], g['hI']
    E = [B.iv(B.dec(x)) for x in source['strict_escape']['point']]
    assert source['strict_escape']['point'] == target['cut_escape']
    e_cov = B.imv(Q, E)
    cut_cov = [-x for x in e_cov]
    rays = [dict(active=row['active'], mask=sum(1 << i for i in row['active']),
                 v=[I(*map(F, x)) for x in row['normalized_coordinate_bounds']])
            for row in cone['final_rays']]
    signs = [B.idot(row['v'], cut_cov).strict_sign() for row in rays]
    neg = [i for i, s in enumerate(signs) if s < 0]
    pos = [i for i, s in enumerate(signs) if s > 0]
    assert len(neg) == 146 and len(pos) == 2311
    out = {rays[i]['mask']: dict(rays[i]) for i in neg}
    crossings = []
    for i in neg:
        for j in pos:
            a, b = rays[i], rays[j]
            common = a['mask'] & b['mask']
            support = B.ids(common)
            if len(support) < 8:
                continue
            minor = B.R.rank_minor([cov[k] for k in support])
            # Eight common independent equations in dimension ten give a
            # two-dimensional face. Its two extreme generators are adjacent.
            assert minor['rank_lower_bound'] == 8, support
            ga, gb = B.idot(a['v'], cut_cov), B.idot(b['v'], cut_cov)
            assert ga.hi < 0 < gb.lo
            ca, cb = gb / (gb-ga), -ga / (gb-ga)
            assert ca.lo > 0 and cb.lo > 0
            v = [ca*x + cb*y for x, y in zip(a['v'], b['v'])]
            mask = common | (1 << 52)
            assert mask not in out
            out[mask] = dict(active=B.ids(mask), mask=mask, v=v)
            crossings.append(dict(negative_support=a['active'], positive_support=b['active'],
                                  new_support=B.ids(mask), rank8_minor=minor,
                                  convex_coefficients=[ca.saved(), cb.saved()]))
    assert len(crossings) == 196 and len(out) == 342
    expected = target['clipping_stages'][-1]
    keys = ('negative_support', 'positive_support', 'new_support')
    normalize = lambda row: tuple(tuple(row[k]) for k in keys)
    assert sorted(map(normalize, crossings)) == sorted(map(normalize, expected['crossings']))
    all_cov = cov + [cut_cov]
    records = []
    for row in sorted(out.values(), key=lambda x: x['active']):
        v, mask = row['v'], row['mask']
        q, hp = B.idot(v, B.imv(Q, v)), B.idot(v, h)
        assert hp.hi < 0
        kind = 'finite_vertices' if q.strict_sign() < 0 else 'hyperideal_rays'
        for k, c in enumerate(all_cov):
            z = B.idot(v, c)
            if mask >> k & 1:
                assert z.lo <= 0 <= z.hi
            else:
                assert z.hi < 0, (row['active'], k)
        projections = []
        for cap in (49, 50, 51):
            pair = I() if mask >> cap & 1 else B.idot(v, cov[cap])
            assert pair.hi <= 0
            ch = B.idot(g['NI'][cap], h)
            projected_norm = q - pair*pair/2
            future_chart = hp - pair*ch/2
            assert projected_norm.hi < 0 and future_chart.hi < 0
            projections.append(dict(cap=cap, projected_norm=projected_norm.saved(),
                                    future_chart=future_chart.saved()))
        records.append(dict(active=row['active'], kind=kind, norm=q.saved(),
                            projections=projections))
    census = {kind: sorted(r['active'] for r in records if r['kind'] == kind)
              for kind in ('finite_vertices', 'cusps', 'hyperideal_rays', 'opposite_sheet_rays')}
    assert census == target['complete_census']
    assert target['certified_dotted_caps'] == [49, 50, 51]
    counts = dict(Counter(r['kind'] for r in records))
    assert counts == dict(finite_vertices=1, hyperideal_rays=341)
    paths = [f'scripts/{TAG}.py', PRIMARY, CONE, LOCAL, ESCAPE,
             f'scripts/{C.TAG}.py', f'scripts/{B.TAG}.py', C.PRIMARY]
    return dict(status='PASS_INDEPENDENT_ALL_REAL_ESCAPE_CUTTER_DOTTED_TO_THREE_CAPS',
                sources={p: B.digest(p) for p in paths}, counts=counts,
                crossings=crossings, complete_census=census, generators=records,
                checked_inequalities=53*len(records), checked_projections=3*len(records),
                certified_dotted_caps=[49, 50, 51],
                proof='Independent complete52cone clipped by -E, with196 rank8-certified crossings. All342 generators satisfy every inequality and project strictly future timelike to each cap hyperplane. Positive conic combinations preserve that property, with no zero projection. For any spacelike candidate r, r^2-(r.cap)^2/2<0 and r.cap<=0 imply r.cap<-sqrt(2*r^2).',
                scope='All real acute spacelike normals with r.E>=0 are strictly dotted to caps49,50,51. No field, angle, normalization or retained-center constraint. This forbids a single cutter that also meets any of these caps; simultaneous separate ridge closures and escape cutters remain open.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], flush=True)
