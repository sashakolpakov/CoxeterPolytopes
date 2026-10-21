#!/usr/bin/env python3
"""Necessary compatibility for the explicit {24-cell, double, Q16} facet spectrum."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/fresh_ideal_ra5_known_carrier_cusp_compatibility_20261018.json'
PINS = {
 'results/fresh_ideal_ra4_double24_volume_20261018.json': '386d9650b8f669e3edcb3cdee56ed53d1c9ea4a90af2b4ba5c8c1b83f9d7bfa3',
 'results/fresh_ideal_ra4_positive_kappa_doubling_20261018.json': 'cf5e71eb3720d16070852c19d3274f0a62b052ed38035d3f9d51adce666deb1c',
 'results/ideal_ra4_positive_kappa_doubling_20261018.md': '3ce0bfd400920c7fc6c4fd9b5b5e2cfab3b0e33fa8f4018a34ce34a84bfed724',
 'results/right_angled_five_face_carrier_20261018.md': 'ce90d9a3fa3c217e6811163a48f84d8ccadaafbe9c99eb4d62ef86de36dc1e0d',
}
NODES = tuple((a, s) for a in range(4) for s in range(2))


def dot(a, b):
    return -a[0]*b[0] + sum(x*y for x, y in zip(a[1:], b[1:]))


def star_records(stage):
    ns = [tuple(map(F, x)) for x in stage['walls']]
    vs = [tuple(map(F, x)) for x in stage['vertices']]
    fv = [[j for j, v in enumerate(vs) if dot(n, v) == 0] for n in ns]
    records = []
    for vi, v in enumerate(vs):
        ids = [i for i, n in enumerate(ns) if dot(n, v) == 0]
        pairs = [(i, j) for i, j in combinations(ids, 2) if dot(ns[i], ns[j]) == -1]
        assert len(ids) == 6 and len(pairs) == 3
        assert sorted(x for p in pairs for x in p) == ids
        axes = []
        for i, j in pairs:
            width = ns[i][0] + ns[j][0]
            assert width > 0
            assert tuple(a+b for a, b in zip(ns[i], ns[j])) == tuple(width*x for x in v)
            axes.append({'width': str(width), 'ridge_vertices': [len(fv[i]), len(fv[j])],
                         'facet_indices': [i, j]})
        records.append({'vertex': vi, 'axes': axes})
    hist = Counter(tuple(sorted(F(a['width']) for a in r['axes'])) for r in records)
    normalized = Counter()
    for w, count in hist.items():
        normalized[tuple(str(x/w[0]) for x in w)] += count
    return {'vertices': len(vs), 'facets': len(ns), 'ridge_size_histogram': dict(sorted(Counter(map(len, fv)).items())),
            'cusp_shape_histogram': [{'shape': list(k), 'count': v} for k, v in sorted(normalized.items())],
            'records': records}


def domains(records, new_width, opposite_type):
    ww = [new_width, 4, 4, 8]
    out = {}
    for node in NODES:
        rem = [a for a in range(4) if a != node[0]]
        labels = set()
        typ = opposite_type if node == (0, 1) else 16
        for r in records[typ]['records']:
            for perm in permutations(range(3)):
                axes = [r['axes'][j] for j in perm]
                if len({F(q['width'])/ww[a] for a, q in zip(rem, axes)}) != 1:
                    continue
                for flips in product((False, True), repeat=3):
                    labs = tuple(x for q, flip in zip(axes, flips)
                                 for x in (q['ridge_vertices'][::-1] if flip else q['ridge_vertices']))
                    labels.add(labs)
        out[node] = sorted(labels)
    out[(0, 0)] = [(24, 24, 24, 24, 14, 14)]
    return out


def solve(dom):
    neighbors = {n: [m for m in NODES if m[0] != n[0]] for n in NODES}
    index = {n: {m: j for j, m in enumerate(neighbors[n])} for n in NODES}
    visits = 0
    def go(chosen):
        nonlocal visits
        visits += 1
        if len(chosen) == 8:
            return chosen
        candidates = {}
        for n in NODES:
            if n in chosen:
                continue
            candidates[n] = [p for p in dom[n] if all(n[0] == m[0] or p[index[n][m]] == q[index[m][n]]
                                                       for m, q in chosen.items())]
            if not candidates[n]:
                return None
        n = min(candidates, key=lambda x: (len(candidates[x]), x))
        for p in candidates[n]:
            answer = go({**chosen, n: p})
            if answer is not None:
                return answer
        return None
    answer = go({})
    return visits, answer


def build():
    for p, h in PINS.items():
        assert sha256((ROOT/p).read_bytes()).hexdigest() == h
    base = json.loads((ROOT/list(PINS)[0]).read_text())
    doubled = json.loads((ROOT/list(PINS)[1]).read_text())
    records = {1: star_records(base['original']), 2: star_records(doubled['stages'][0]),
               16: star_records(doubled['stages'][-1])}
    assert records[16]['ridge_size_histogram'] == {6: 128, 9: 48, 14: 18, 24: 4}
    assert records[16]['cusp_shape_histogram'] == [
        {'shape': ['1','1','1'], 'count':195}, {'shape':['1','1','2'], 'count':49},
        {'shape':['1','2','2'], 'count':14}]
    special = [r for r in records[16]['records']
               if sorted(x for a in r['axes'] for x in a['ridge_vertices']) == [14,14,24,24,24,24]]
    assert len(special) == 1 and special[0]['vertex'] == 85
    assert sorted((F(a['width']), sorted(a['ridge_vertices'])) for a in special[0]['axes']) == [
        (F(4), [24,24]), (F(4), [24,24]), (F(8), [14,14])]
    cases = []
    for width, typ in product((4,8), (2,16)):
        dom = domains(records, width, typ)
        visits, answer = solve(dom)
        assert (answer is not None) == (width == 4 and typ == 16)
        cases.append({'new_width': width, 'opposite_type': typ,
                      'domains': [{'node':list(n), 'labels':dom[n]} for n in NODES],
                      'search_nodes': visits, 'compatible': answer is not None,
                      'witness': None if answer is None else [{'node':list(n), 'labels':answer[n]} for n in NODES]})
    # An abstract arithmetic witness; no geometric chamber or full boundary is inferred.
    x,y,z = 361,66,30
    a,b,c,d = 2051,241,105,0
    assert 6*z-2*y == 48 and z%6 == 0
    assert (8*a+2*b+2*d, 6*b+4*c, 4*c+6*d) == (24*x+36*y+195*z, 6*y+49*z, 14*z)
    assert x+y+z == 457 and a+b+c+d == 2397
    assert 24*x+38*y+198*z == 2*(3*2397+3*457-6)
    return {'status':'RESTRICTED_FACET_SPECTRUM_COUPLED_CUSP_NECESSARY_CONSTRAINTS',
            'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(), 'sources':PINS,
            'facet_spectrum':[1,2,16], 'facet_types':records, 'special_vertex':85,
            'local_cases':cases, 'necessary_counts':{'Q16_multiple':6,'Q16_minimum':30,'double24_minimum':66},
            'abstract_incidence_witness':{'facet_counts':[x,y,z],'ambient_cusp_types':[a,b,c,d],
                                         'facets':457,'vertices':2397,'is_geometric_realization':False},
            'global_claims':{'ideal_RA5_example':False,'spectrum_excluded':False,'unrestricted_RA5_excluded':False},
            'scope':'exact cusp stars, necessary width/ridge-size compatibility, and universal counts only under the stated three-type facet restriction'}


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    data=build()
    if args.verify:
        assert json.loads(OUT.read_text()) == json.loads(json.dumps(data))
    else:
        assert not OUT.exists();OUT.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'status':'PASS','verify':args.verify,'stars':24+42+258,
                      'local_cases':[(q['new_width'],q['opposite_type'],q['compatible'],q['search_nodes']) for q in data['local_cases']],
                      'Q16_minimum':30,'double24_minimum':66}))
