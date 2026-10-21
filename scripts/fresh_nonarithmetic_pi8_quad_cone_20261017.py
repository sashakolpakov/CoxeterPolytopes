#!/usr/bin/env python3
"""Complete ambient cone after the two-wall r8 quadrilateral construction."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_a2_escape_partition_20261017 as P

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_pi8_quad_cone_20261017'
CAP = 'fresh_nonarithmetic_pi8_quad_completion_20261017'
TETRA = 'fresh_nonarithmetic_pi8_tetrahedron_completion_20261017'
PRIOR = 'fresh_nonarithmetic_pi8_tetrahedron_cone_20261017'
U, I = P.U, P.I


def read(tag):
    obj = json.loads((ROOT/f'results/{tag}.json').read_text())
    for path, digest in obj['sources'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    return obj


def sign(x):
    if x.l > 0:
        return 1
    if x.h < 0:
        return -1
    raise AssertionError(('unresolved sign', x.enc()))


def ids(mask):
    return [i for i in range(mask.bit_length()) if (mask >> i) & 1]


def vector(raw):
    return [U.rnd(I(*x)) for x in raw]


def rank8(rows):
    a = [row[:] for row in rows]
    row_ids, col_ids = list(range(len(a))), list(range(10))
    pivots, selected = [], []
    for k in range(8):
        candidates = [(abs((a[i][j].l+a[i][j].h)/2), i, j)
                      for i in range(k, len(a)) for j in range(k, 10)
                      if a[i][j].l > 0 or a[i][j].h < 0]
        assert candidates, ('rank below8', k)
        _, i, j = max(candidates)
        a[k], a[i] = a[i], a[k]
        row_ids[k], row_ids[i] = row_ids[i], row_ids[k]
        for row in a:
            row[k], row[j] = row[j], row[k]
        col_ids[k], col_ids[j] = col_ids[j], col_ids[k]
        pivot = a[k][k]
        pivots.append(pivot.enc())
        selected.append([row_ids[k], col_ids[k]])
        for i in range(k+1, len(a)):
            factor = U.rnd(a[i][k]/pivot)
            for j in range(k+1, 10):
                a[i][j] = U.rnd(a[i][j]-factor*a[k][j])
    return dict(pivots=pivots, pivot_positions=selected)


def tighten(v, support, cov):
    """Same exact projective ray, reconstructed from nine active equations."""
    anchor = max(range(10), key=lambda j: abs((v[j].l+v[j].h)/2))
    orientation = sign(v[anchor])
    columns = [j for j in range(10) if j != anchor]
    rows = [[cov[i][j] for j in columns] for i in support]
    work = [r[:] for r in rows]
    row_ids = list(range(len(support)))
    for k in range(9):
        candidates = [(abs((work[i][j].l+work[i][j].h)/2), i, j)
                      for i in range(k, len(work)) for j in range(k, 9)
                      if work[i][j].l > 0 or work[i][j].h < 0]
        assert candidates, ('ray reconstruction rank', support, k)
        _, i, j = max(candidates)
        work[k], work[i] = work[i], work[k]
        row_ids[k], row_ids[i] = row_ids[i], row_ids[k]
        for row in work:
            row[k], row[j] = row[j], row[k]
        pivot = work[k][k]
        for i in range(k+1, len(work)):
            factor = U.rnd(work[i][k]/pivot)
            for j in range(k+1, 9):
                work[i][j] = U.rnd(work[i][j]-factor*work[k][j])
    selected = [support[i] for i in row_ids[:9]]
    matrix = [[cov[i][j] for j in columns] for i in selected]
    inverse, pivots = P.C53.inverse(matrix)
    solution = U.imv(inverse, [-orientation*cov[i][anchor] for i in selected])
    out = [I(orientation) if j == anchor else solution[columns.index(j)] for j in range(10)]
    # The clipping point has the same active equations and nonzero anchor.
    # The invertible selected minor proves uniqueness, hence exact equality
    # after positive normalization; interval overlap is only a sanity check.
    for x, y in zip(out, v):
        previous = U.rnd(orientation*y/v[anchor])
        assert x.l <= previous.h and previous.l <= x.h
    return out, dict(anchor_coordinate=anchor, anchor_value=orientation,
                     selected_active_rows=selected, inverse_pivots=pivots)


def run():
    old, cap, tetra = read(PRIOR), read(CAP), read(TETRA)
    assert old['status'] == 'PASS_COMPLETE_58_WALL_CONE_AND_STRICT_ESCAPE'
    assert cap['status'] == 'PASS_EXACT_60_WALL_CHAMBER_WITH_COMPACT_R8_QUADRILATERAL_FACE'
    _, q, normals, _, _, _, _, _ = P.physical_data()
    normals.extend(vector(row['physical']) for row in tetra['new_normals'])
    normals.extend(vector(row['physical']) for row in cap['new_normals'])
    assert len(normals) == 60
    cov = [U.imv(q, v) for v in normals]
    h = vector(cap['physical_common_center'])
    hcov = U.imv(q, h)
    assert U.idot(h, hcov).h < 0
    assert all(U.idot(h, c).h < 0 for c in cov)
    rays = [dict(active=r['active'], mask=sum(1 << j for j in r['active']),
                 kind=r['kind'], v=vector(r['normalized_coordinate_bounds'])) for r in old['final_rays']]
    assert len(rays) == 2798
    stages = []
    for wall in (58, 59):
        products = [U.idot(r['v'], cov[wall]) for r in rays]
        signs = [sign(x) for x in products]
        negative = [i for i, s in enumerate(signs) if s < 0]
        positive = [i for i, s in enumerate(signs) if s > 0]
        assert all(rays[i]['kind'] == 'hyperideal_rays' for i in positive)
        out = {rays[i]['mask']: dict(rays[i]) for i in negative}
        crossings = []
        for ni in negative:
            a = rays[ni]
            for pi in positive:
                b = rays[pi]
                common = a['mask'] & b['mask']
                if common.bit_count() < 8:
                    continue
                support = ids(common)
                minor = rank8([cov[k] for k in support])
                ga, gb = products[ni], products[pi]
                ca, cb = U.rnd(gb/(gb-ga)), U.rnd(-ga/(gb-ga))
                assert ca.l > 0 and cb.l > 0
                v = [U.rnd(ca*x+cb*y) for x, y in zip(a['v'], b['v'])]
                mask = common | (1 << wall)
                v, reconstruction = tighten(v, ids(mask), cov)
                norm = U.idot(v, U.imv(q, v))
                if norm.l <= 0 <= norm.h:
                    print('UNRESOLVED_NEW_NORM', wall, ids(common | (1 << wall)),
                          'parent_supports', a['active'], b['active'],
                          'float_norm_interval', float(norm.l), float(norm.h),
                          'max_coordinate', max(abs(float(x.l)) for x in v), flush=True)
                kind = 'finite_vertices' if sign(norm) < 0 else 'hyperideal_rays'
                if kind == 'finite_vertices':
                    assert U.idot(v, hcov).h < 0
                assert mask not in out, ('duplicate ray', ids(mask))
                out[mask] = dict(active=ids(mask), mask=mask, kind=kind, v=v)
                crossings.append(dict(negative_support=a['active'], positive_support=b['active'],
                                      new_support=ids(mask), rank8_minor=minor, kind=kind,
                                      positive_coefficients=P.enc([ca, cb]), norm=norm.enc(),
                                      direct_ray_reconstruction=reconstruction))
        rays = sorted(out.values(), key=lambda r: r['active'])
        stage = dict(wall=wall, negative=len(negative), positive=len(positive),
                     crossings=crossings, counts=dict(Counter(r['kind'] for r in rays)))
        stages.append(stage)
        print('COMPLETE60_CLIP', wall, stage['negative'], stage['positive'], len(crossings), stage['counts'], flush=True)
    final, counts = [], Counter()
    for number, ray in enumerate(rays, 1):
        v = ray['v']
        norm, time = U.idot(v, U.imv(q, v)), U.idot(v, hcov)
        kind = 'finite_vertices' if sign(norm) < 0 else 'hyperideal_rays'
        assert kind == ray['kind']
        if kind == 'finite_vertices':
            assert time.h < 0
        for i, c in enumerate(cov):
            value = U.idot(v, c)
            if (ray['mask'] >> i) & 1:
                assert value.l <= 0 <= value.h
            else:
                assert value.h < 0, (ray['active'], i, value.enc())
        final.append(dict(active=ray['active'], kind=kind, norm=norm.enc(),
                          center_product=time.enc(), normalized_coordinate_bounds=P.enc(v)))
        counts[kind] += 1
        if number % 500 == 0:
            print('COMPLETE60_FINAL_GATES', number, flush=True)
    old_finite = {tuple(r['active']) for r in old['final_rays'] if r['kind'] == 'finite_vertices'}
    new_finite = {tuple(r['active']) for r in final if r['kind'] == 'finite_vertices'}
    assert old_finite <= new_finite
    assert all(tuple(sorted(v['incident'])) in new_finite for v in tetra['tetrahedron_vertices'])
    assert all(tuple(sorted(v['incident'])) in new_finite for v in cap['complete_six_pair_intersections'] if v['kind'] == 'finite_vertex')
    faces = []
    for name, support in [('original_H4_A2', P.S), ('new_H4_I2_8', tetra['new_face_support']), ('r8_quadrilateral', cap['face_support'])]:
        rows = [r for r in final if set(support) <= set(r['active'])]
        faces.append(dict(name=name, support=support, counts=dict(Counter(r['kind'] for r in rows)),
                          ray_supports=[r['active'] for r in rows]))
    assert faces[1]['counts'] == faces[2]['counts'] == {'finite_vertices': 4}
    exterior = next(r for r in rays if r['kind'] == 'hyperideal_rays' and U.idot(r['v'], hcov).h < 0)
    for power in range(128):
        epsilon = I(F(1, 2**power))
        escape = [U.rnd(x+epsilon*y) for x, y in zip(exterior['v'], h)]
        norm, future = U.idot(escape, U.imv(q, escape)), U.idot(escape, hcov)
        if norm.l > 0 and future.h < 0:
            products = [U.idot(escape, c) for c in cov]
            if all(x.h < 0 for x in products):
                break
    else:
        raise AssertionError('strict escape perturbation failed')
    paths = [f'scripts/{TAG}.py', f'results/{CAP}.json', f'scripts/{CAP}.py',
             f'results/{PRIOR}.json', f'scripts/{PRIOR}.py', f'results/{TETRA}.json', f'scripts/{TETRA}.py', f'scripts/{P.TAG}.py',
             f'results/{P.TAG}.json', P.U.CAP, P.C53.CL, P.CAP52, P.CAPS55]
    return dict(status='PASS_COMPLETE_60_WALL_CONE_AND_STRICT_ESCAPE',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                initial_counts=old['counts'], counts=dict(counts), clipping_stages=stages,
                final_rays=final, new_finite_supports=[list(x) for x in sorted(new_finite-old_finite)],
                final_inequality_checks=60*len(final), faces=faces,
                strict_escape=dict(ray_support=exterior['active'], center_coefficient_power=power,
                                   point=P.enc(escape), norm=norm.enc(), future=future.enc(),
                                   all60_products=P.enc(products)),
                proof='Start from the exact complete58-wall cone, whose source chain starts from the independently checked55-wall cone. For each new halfspace every old ray has a strict sign. Every opposite-sign pair with at least8 common active normals is tested; an exact rank8 minor proves an edge, since two distinct rays force rank at most8. Positive combinations give every new extreme ray. Each new point is tightened by solving9 independent active equations with one nonzero coordinate fixed to its sign; uniqueness proves it is the same exact ray after positive normalization. All final strict/nonincident gates and norms are checked. Active equalities follow from clipping identities. The old pointed cone and common center fix the halfspace orientation.',
                scope='Complete ambient cone after the two prescribed additions, using inherited exact root definitions and outward rational intervals. This is a new primary clipping certificate, not an independent reconstruction of the new algebraic roots. A strict spacelike escape proves infinite volume; no new finite-volume polytope is claimed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], data['counts'], flush=True)
