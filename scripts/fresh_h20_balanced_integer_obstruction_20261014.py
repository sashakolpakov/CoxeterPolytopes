#!/usr/bin/env python3
"""Exact exhaustive radius-two integer-gap scout at two specified probes.

This verifies a finite cap family, not unrestricted Coxeter completion.
All distance decisions use integers, scaled by 28 times the probe denominator.
For an endpoint cap of radius sqrt(28), the contact cosine is rational/sqrt(7).
Its square is rational; the finite-angle possibilities are 0,1/4,1/2,3/4.
Only zero can be a rational square divided by7. Thus right/dotted is the
complete gate, not just a sufficient selection criterion.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_balanced_integer_obstruction_20261014'


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for i in range(total+1):
            for rest in compositions(total-i, length-1):
                yield (i,)+rest


def factor(n):
    inverse28 = [[56*min(i, j)*(n+1-max(i, j))//(n+1)
                  for j in range(1, n+1)] for i in range(1, n+1)]
    assert all(56*min(i, j)*(n+1-max(i, j)) % (n+1) == 0
               for i in range(1, n+1) for j in range(1, n+1))
    rows = []
    for gaps in compositions(n+1, n+1):
        scalar = [1-g for g in gaps[1:]]
        coeff28 = [sum(a*b for a, b in zip(row, scalar)) for row in inverse28]
        rows.append(dict(gaps=gaps, scalar=scalar, coeff28=coeff28,
                         norm28=sum(a*b for a, b in zip(scalar, coeff28)),
                         dot_vertex28=sum(coeff28)))
    vertices = []
    for omitted in range(n+1):
        scalar = [1-(n+1)*(j == omitted) for j in range(1, n+1)]
        coeff28 = [sum(a*b for a, b in zip(row, scalar)) for row in inverse28]
        assert all(x % 28 == 0 for x in coeff28)
        vertices.append([x//28 for x in coeff28])
    return rows, vertices, sum(map(sum, inverse28))


def run():
    factors = [factor(n) for n in (6, 6, 7)]
    vertices = [sum(parts, []) for parts in product(*[f[1] for f in factors])]
    assert len(vertices) == 392
    probes = []
    for p, q in ((2, 5), (3, 5)):
        selected = []
        for rows, _, vnorm28 in factors:
            candidates = []
            for r in rows:
                distance = q*q*r['norm28']-2*p*q*r['dot_vertex28']+p*p*vnorm28
                if distance <= 112*q*q:
                    candidates.append((r, distance))
            selected.append(candidates)
        records, tangent_records = [], []
        for parts in product(*selected):
            power_num = sum(x[1] for x in parts)-112*q*q
            if power_num > 0:
                continue
            rows = [x[0] for x in parts]
            norm28 = sum(r['norm28'] for r in rows)
            scalar = sum([r['scalar'] for r in rows], [])
            central_ok = norm28 in (224, 336) or norm28 >= 448
            bad = None
            for vi, vertex in enumerate(vertices):
                distance28 = norm28+196*28-56*sum(a*b for a, b in zip(scalar, vertex))
                t = distance28-32*28
                if t == 0 or (t > 0 and t*t >= 448*28*28):
                    continue
                bad = dict(vertex=vi, distance28=distance28)
                break
            record = dict(gaps=[list(r['gaps']) for r in rows],
                          probe_power_numerator=power_num,
                          norm28=norm28, central_ok=central_ok,
                          first_incompatible_endpoint=bad)
            (records if power_num < 0 else tangent_records).append(record)
        accepted = sum(r['central_ok'] and r['first_incompatible_endpoint'] is None
                       for r in records)
        assert len(records) == 555
        assert accepted == (555 if p == 2 else 0)
        tangent_accepted = sum(r['central_ok'] and r['first_incompatible_endpoint'] is None
                               for r in tangent_records)
        assert p != 3 or tangent_accepted == 0
        probes.append(dict(multiple=f'{p}/{q}', strict_cutters=len(records),
                           core_compatible=accepted, power_denominator=28*q*q,
                           tangent_caps=len(tangent_records),
                           core_compatible_tangents=tangent_accepted,
                           tangent_records=tangent_records, records=records))
    paths = [f'scripts/{TAG}.py', 'results/fresh_h20_balanced_integer_caps_20261014.json']
    old = json.loads((ROOT/paths[1]).read_text())
    assert {tuple(map(tuple, r['gaps'])) for r in old['recipes']} == {
        tuple(map(tuple, r['gaps'])) for r in probes[0]['records']}
    return dict(status='PASS_EXACT_INTEGER_CAP_CENSUS_AND_WHOLE_FAMILY_ESCAPE',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                factor_composition_counts=[len(f[0]) for f in factors],
                probes=probes,
                whole_family_escape=dict(spatial_multiple='3/5', vertex=0,
                                         D='98783/2800', U=1, norm='1/1400',
                                         integer_cap_power_lower_bound='1/700',
                                         central_power='1664/25',
                                         endpoint_power_lower_bound='84/25',
                                         minimum_seed_gap='2/5'),
                scope='All radius2 caps with nonnegative integer seed gaps summing7,7,8. The .4 census exactly matches the numerical scout. No such cap cutting or tangent to .6z_000 is compatible with the central cap and all392 Rsqrt28 endpoint caps. Every cap power lies in(1/700)Z, so all core-compatible caps have power at least1/700. The displayed norm1/1400 point escapes the entire family, even without mutual compatibility. Cyclic seed symmetry transports this statement to every old vertex. This excludes finite-volume completion using only this family while retaining the baseline. It does not exclude other radii, noninteger gaps, or changed endpoint layers.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == out
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], [(p['multiple'], p['strict_cutters'], p['core_compatible'])
                          for p in out['probes']], flush=True)
