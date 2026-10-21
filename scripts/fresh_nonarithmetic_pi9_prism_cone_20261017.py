#!/usr/bin/env python3
"""Complete ambient cone after the two-wall H4+I2(9) prism completion."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_pi8_quad_cone_20261017 as Q
import fresh_nonarithmetic_pi9_prism_completion_20261017 as B

P, U, I = Q.P, Q.U, Q.I
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_pi9_prism_cone_20261017'
CAP = 'fresh_nonarithmetic_pi9_prism_completion_20261017'
PRIOR = Q.TAG


def run():
    old, cap = Q.read(PRIOR), Q.read(CAP)
    tetra, quad = Q.read(Q.TETRA), Q.read(Q.CAP)
    assert old['status'] == 'PASS_COMPLETE_60_WALL_CONE_AND_STRICT_ESCAPE'
    assert cap['status'] == 'PASS_EXACT_62_WALL_CHAMBER_WITH_COMPACT_PI9_TRIANGULAR_PRISM'
    geometry = B.geometry()
    q, normals = geometry['exact_gram'], geometry['normals']
    assert len(normals) == 62
    cov = [U.imv(q, n) for n in normals]
    h = Q.vector(cap['common_center'])
    hcov = U.imv(q, h)
    assert U.idot(h, hcov).h < 0 and all(U.idot(h, c).h < 0 for c in cov)
    rays = []
    reconstruction_digest = sha256()
    for i, r in enumerate(old['final_rays']):
        converted = geometry['coordinates'](Q.vector(r['normalized_coordinate_bounds']))
        v, reconstruction = Q.tighten(converted, r['active'], cov)
        reconstruction_digest.update(json.dumps(reconstruction, sort_keys=True).encode())
        rays.append(dict(active=r['active'], mask=sum(1 << j for j in r['active']), kind=r['kind'], v=v))
        if i % 500 == 0:
            print('COMPLETE62_RECONSTRUCT_OLD', i, flush=True)
    assert len(rays) == 2868
    stages = []
    for wall in (60, 61):
        products = [U.idot(r['v'], cov[wall]) for r in rays]
        signs = []
        for r, x in zip(rays, products):
            if x.l <= 0 <= x.h:
                print('UNRESOLVED_CLIP_SIGN', wall, r['active'], x.enc(), flush=True)
            signs.append(Q.sign(x))
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
                support = Q.ids(common)
                minor = Q.rank8([cov[k] for k in support])
                ga, gb = products[ni], products[pi]
                ca, cb = U.rnd(gb/(gb-ga)), U.rnd(-ga/(gb-ga))
                assert ca.l > 0 and cb.l > 0
                v = [U.rnd(ca*x+cb*y) for x, y in zip(a['v'], b['v'])]
                mask = common | (1 << wall)
                v, reconstruction = Q.tighten(v, Q.ids(mask), cov)
                norm = U.idot(v, U.imv(q, v))
                if norm.l <= 0 <= norm.h:
                    print('UNRESOLVED_CLIP_NORM', wall, Q.ids(mask), norm.enc(), flush=True)
                kind = 'finite_vertices' if Q.sign(norm) < 0 else 'hyperideal_rays'
                if kind == 'finite_vertices':
                    assert U.idot(v, hcov).h < 0
                assert mask not in out, ('duplicate ray', Q.ids(mask))
                out[mask] = dict(active=Q.ids(mask), mask=mask, kind=kind, v=v)
                crossings.append(dict(negative_support=a['active'], positive_support=b['active'],
                                      new_support=Q.ids(mask), rank8_minor=minor, kind=kind,
                                      positive_coefficients=P.enc([ca, cb]), norm=norm.enc(),
                                      direct_ray_reconstruction=reconstruction))
        rays = sorted(out.values(), key=lambda r: r['active'])
        stage = dict(wall=wall, negative=len(negative), positive=len(positive),
                     crossings=crossings, counts=dict(Counter(r['kind'] for r in rays)))
        stages.append(stage)
        print('COMPLETE62_CLIP', wall, len(negative), len(positive), len(crossings), stage['counts'], flush=True)
    final, counts = [], Counter()
    for number, ray in enumerate(rays, 1):
        v = ray['v']
        norm, time = U.idot(v, U.imv(q, v)), U.idot(v, hcov)
        kind = 'finite_vertices' if Q.sign(norm) < 0 else 'hyperideal_rays'
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
                          center_product=time.enc(), normalized_computational_coordinate_bounds=P.enc(v)))
        counts[kind] += 1
        if number % 500 == 0:
            print('COMPLETE62_FINAL_GATES', number, flush=True)
    old_finite = {tuple(r['active']) for r in old['final_rays'] if r['kind'] == 'finite_vertices'}
    new_finite = {tuple(r['active']) for r in final if r['kind'] == 'finite_vertices'}
    assert old_finite <= new_finite
    supplied_vertices = [r for r in cap['complete_triple_intersections'] if r['kind'] == 'finite_vertex']
    assert len(supplied_vertices) == 6
    assert all(tuple(sorted(r['incident'])) in new_finite for r in supplied_vertices)
    faces = []
    targets = [('original_H4_A2', P.S), ('compact_H4_I2_8_tetrahedron', tetra['new_face_support']),
               ('r8_quadrilateral', quad['face_support']), ('new_H4_I2_9_prism', cap['face_support'])]
    for name, support in targets:
        rows = [r for r in final if set(support) <= set(r['active'])]
        faces.append(dict(name=name, support=support, counts=dict(Counter(r['kind'] for r in rows)),
                          ray_supports=[r['active'] for r in rows]))
    assert faces[1]['counts'] == faces[2]['counts'] == {'finite_vertices': 4}
    assert faces[3]['counts'] == {'finite_vertices': 6}
    escape = None
    exterior = next((r for r in rays if r['kind'] == 'hyperideal_rays' and U.idot(r['v'], hcov).h < 0), None)
    if exterior is not None:
        for power in range(128):
            point = [U.rnd(x+I(F(1, 2**power))*y) for x, y in zip(exterior['v'], h)]
            norm, future = U.idot(point, U.imv(q, point)), U.idot(point, hcov)
            if norm.l > 0 and future.h < 0:
                products = [U.idot(point, c) for c in cov]
                if all(x.h < 0 for x in products):
                    escape = dict(ray_support=exterior['active'], center_coefficient_power=power,
                                  computational_point=P.enc(point), physical_point=P.enc(P.combine(point, geometry['basis'])),
                                  norm=norm.enc(), future=future.enc(), all62_products=P.enc(products))
                    break
        assert escape is not None
    else:
        assert counts['hyperideal_rays'] == 0
    paths = [f'scripts/{TAG}.py', f'results/{CAP}.json', f'scripts/{CAP}.py',
             f'results/{PRIOR}.json', f'scripts/{PRIOR}.py', f'results/{Q.TETRA}.json',
             f'results/{Q.CAP}.json', f'scripts/{P.TAG}.py', f'results/{P.TAG}.json',
             P.U.CAP, P.C53.CL, P.CAP52, P.CAPS55]
    return dict(status='PASS_COMPLETE_62_WALL_CONE_AND_STRICT_ESCAPE' if escape else 'PASS_COMPLETE_62_WALL_COMPACT_CONE',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                initial_counts=old['counts'], counts=dict(counts), clipping_stages=stages,
                coordinate_system='Exact computational root basis, with old physical rays converted and re-solved from active equations.',
                computational_root_basis=geometry['basis_ids'], computational_gram=[P.enc(r) for r in q],
                computational_common_center=P.enc(h), initial_ray_reconstruction_sha256=reconstruction_digest.hexdigest(),
                final_rays=final, new_finite_supports=[list(x) for x in sorted(new_finite-old_finite)],
                final_inequality_checks=62*len(final), faces=faces, strict_escape=escape,
                proof='Start from the frozen complete60-wall cone. Convert its physical rays to the exact root basis, then solve nine independent active equations to obtain the same exact oriented rays with narrower interval bounds. Strict signs partition every old ray at each cut. Opposite-sign pairs with8independent common active normals are exactly crossed edges; positive combinations give all new rays. Nine active equations reconstruct each ray exactly after positive scaling, using the frozen rank8/tighten helpers. All62final inequalities and causal signs are checked. All old finite corners persist; the retained tetrahedron/quadrilateral and new prism are identified by their complete ambient active-ray sets.',
                scope='Complete ambient cone for these two prescribed exact additions. Root definitions and prior cone completeness are inherited; this is a new primary clipping certificate, not an independent arithmetic reconstruction. If a strict escape is exported, it proves infinite volume.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['counts'], flush=True)
