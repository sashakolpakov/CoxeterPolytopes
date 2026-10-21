#!/usr/bin/env python3
"""Exact right-angle vertex caps with unrestricted dotted polar pairing."""
import json
from fractions import Fraction
import numpy as np
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h4_blocks_independent_verify import inverse
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K


def nullspace_K(matrix):
    a = [list(row) for row in matrix]
    pivots = []
    for j in range(len(a[0])):
        i = next((i for i in range(len(pivots), len(a)) if K.sign(a[i][j])), None)
        if i is None:
            continue
        a[len(pivots)], a[i] = a[i], a[len(pivots)]
        i = len(pivots)
        p = a[i][j]
        a[i] = [div(t, p) for t in a[i]]
        for k in range(len(a)):
            if k != i:
                factor = a[k][j]
                a[k] = [sub(x, mul(factor, y)) for x, y in zip(a[k], a[i], strict=True)]
        pivots.append(j)
        if len(pivots) == len(a):
            break
    free = [j for j in range(len(a[0])) if j not in pivots]
    basis = []
    for j in free:
        v = [Z]*len(a[0])
        v[j] = O
        for i, p in enumerate(pivots):
            v[p] = K.neg(a[i][j])
        basis.append(tuple(v))
    return tuple(free), tuple(basis)


def decode(record, name):
    d = parse(record['radicand'])
    f = Quadratic(K, d)
    raw = tuple(map(parse, record['raw_seed_pairings']))
    coeff = tuple(map(parse, record['raw_seed_coefficients']))
    return {'name': name, 'field': f, 'p': tuple((Z, t) for t in raw),
            'c': tuple((Z, t) for t in coeff), 'q': scale(O, 2)}


def finite_vertex(block, raw, ps, I, h):
    G = [[dot(K, ps[a]['c'], ps[b]['p']) for b in block] for a in block]
    assert inertia_K(G) == (8, 0, 0)
    gh = tuple(dot(K, ps[n]['c'], h) for n in block)
    alpha = matvec(K, inverse(G), gh)
    y = list(h)
    for n, t in zip(block, alpha, strict=True):
        y = [sub(a, mul(t, b)) for a, b in zip(y, ps[n]['p'], strict=True)]
    c = matvec(K, I, raw)
    t = div(dot(K, c, h), dot(K, c, raw))
    y = tuple(sub(a, mul(t, b)) for a, b in zip(y, raw, strict=True))
    norm, cp = pairing(K, y, I, y), pairing(K, y, I, h)
    assert K.sign(norm) < 0 and K.sign(cp) < 0
    signs = {n: r['field'].sign(dot(r['field'], r['c'],
              tuple(lift_K(r['field'], t) for t in y))) for n, r in ps.items()}
    assert all(s == (0 if n in block else 1) for n, s in signs.items())
    assert K.sign(dot(K, c, y)) == 0
    return {'y': list(map(enc, y)), 'norm': enc(norm), 'center_pairing': enc(cp),
            'old_wall_signs': signs, 'new_cap_sign': 0,
            'block_gram': [[enc(t) for t in row] for row in G]}


def run():
    G, I, h, ps, source, case = configuration('48')
    data = read(PREFIX+'_vertex_proposals.json')
    pool = angle_pool(ANGLE_SOURCE)
    caps, rejected = [], []
    for index, proposal in enumerate(data['proposals']):
        if proposal['kind'] != 'eight_right_others_dotted':
            continue
        block = proposal['block_names']
        if any(ps[n]['field'] is not K for n in block):
            continue
        free, basis = nullspace_K([ps[n]['c'] for n in block])
        assert len(free) == len(basis) == 2
        pn = proposal['seed_pairings_numeric']
        if abs(pn[free[1]]) < 1e-9:
            continue
        ratio = Fraction(float(pn[free[0]]/pn[free[1]])).limit_denominator(1000)
        orientation = 1 if pn[free[1]] > 0 else -1
        raw = tuple(scale(add(scale(x, ratio), y), orientation)
                    for x, y in zip(*basis, strict=True))
        qraw = pairing(K, raw, I, raw)
        if K.sign(qraw) <= 0:
            rejected.append({'proposal': index, 'reason': 'nonpositive_raw_norm'})
            continue
        d = div(scale(O, 2), qraw)
        row = {'kind': 'eight_right_others_dotted', 'proposal_index': index,
               'block_names': block, 'free_coordinates': list(free), 'rational_ratio': str(ratio),
               'orientation': orientation, 'radicand': enc(d),
               'raw_seed_pairings': list(map(enc, raw)),
               'raw_seed_coefficients': list(map(enc, matvec(K, I, raw)))}
        cap = decode(row, f'vertex{len(caps)}')
        checks = {n: relation(cap, r, pool) for n, r in ps.items()}
        if any(r != ({'kind': 'finite', 'label': 2} if n in block else {'kind': 'dotted'})
               for n, r in checks.items()):
            rejected.append({'proposal': index, 'reason': 'exact_angle_gate',
                             'bad': {n: r for n, r in checks.items()
                                     if r != ({'kind': 'finite', 'label': 2} if n in block else {'kind': 'dotted'})}})
            continue
        f = cap['field']
        assert f.sign(f.sub(dot(f, cap['c'], cap['p']), f.rational(2))) == 0
        assert f.sign(dot(f, cap['c'], tuple(lift_K(f, t) for t in h))) > 0
        row['pair_checks'] = checks
        row['finite_vertex'] = finite_vertex(block, raw, ps, I, h)
        caps.append(row)
        print({'cap': len(caps)-1, 'block': block, 'ratio': str(ratio)}, flush=True)
    result = {'status': 'EXACT_FINITE_VERTEX_CREATING_CAPS', 'source': source,
              'source_case': case, 'proposal_source': PREFIX+'_vertex_proposals.json',
              'caps': caps, 'rejected': rejected,
              'scope': 'Accepted caps have norm2, eight exact right angles, all remaining '
                       'old pairs dotted, and an exact actual finite vertex. This alone '
                       'does not prove finite volume or joint compatibility.'}
    (ROOT/(PREFIX+'_vertex_caps.json')).write_text(json.dumps(result, indent=2)+'\n')
    print({'accepted': len(caps), 'rejected': len(rejected)}, flush=True)


if __name__ == '__main__':
    run()
