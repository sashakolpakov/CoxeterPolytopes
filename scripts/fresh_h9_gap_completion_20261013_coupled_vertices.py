#!/usr/bin/env python3
"""Finite-angle coupling between two rank-eight orthogonal cap families."""
import json
from itertools import combinations
import numpy as np
import mpmath as mp
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h9_gap_completion_20261013_vertex_caps import decode, nullspace_K


def exact(block, first, label, branch, G, I, h, old):
    free, basis = nullspace_K([old[n]['c'] for n in block])
    assert len(basis) == 2
    rawc = tuple(map(parse, first['raw_seed_coefficients']))
    d = parse(first['radicand'])
    E = Quadratic(K, d)
    a, b = [dot(K, rawc, v) for v in basis]
    if not K.sign(b):
        basis = basis[::-1]
        a, b = b, a
    assert K.sign(b)
    k = {2: Z, 3: O, 5: P}[label]
    target = (Z, K.neg(div(k, d)))
    line0 = tuple(E.mul(E.lift(div(v, b)), target) for v in basis[1])
    line1 = tuple(E.lift(sub(v, mul(div(a, b), w))) for v, w in zip(*basis, strict=True))
    IE = tuple(tuple(E.lift(t) for t in row) for row in I)
    A = pairing(E, line1, IE, line1)
    B = E.mul(E.rational(2), pairing(E, line0, IE, line1))
    C = E.sub(pairing(E, line0, IE, line0), E.rational(2))
    D = E.sub(E.mul(B, B), E.mul(E.rational(4), E.mul(A, C)))
    assert E.sign(D) > 0 and E.sign(A)
    F = Quadratic(E, D)
    x = (E.div(E.neg(B), E.mul(E.rational(2), A)),
         E.div(E.rational(branch), E.mul(E.rational(2), A)))
    p = tuple(F.add(F.lift(v), F.mul(x, F.lift(w))) for v, w in zip(line0, line1, strict=True))
    c = matvec(F, tuple(tuple(lift_K(F, t) for t in row) for row in I), p)
    assert F.sign(F.sub(dot(F, c, p), F.rational(2))) == 0
    assert F.sign(dot(F, c, tuple(lift_K(F, t) for t in h))) > 0
    row = {'field': F, 'p': p, 'c': c, 'q': scale(O, 2)}
    return row, {'kind': 'coupled_eight_right', 'block_names': block,
        'anchor_label': label, 'branch': branch, 'base_radicand': enc(d),
        'line_constant': list(map(E.encode, line0)), 'line_direction': list(map(E.encode, line1)),
        'norm_polynomial': list(map(E.encode, (A, B, C))), 'discriminant': E.encode(D),
        'parameter': F.encode(x), 'seed_pairings': list(map(F.encode, p)),
        'seed_coefficients': list(map(F.encode, c))}


def run():
    G, I, h, old, source, case = configuration('48')
    bank = read(PREFIX+'_vertex_caps.json')
    mp.mp.dps = 80
    C = np.array([[float(value(r['field'], t)) for t in r['c']] for r in old.values()])
    Qn = np.array([float(value(K, r['q'])) for r in old.values()])
    In = np.array([[float(value(K, t)) for t in row] for row in I])
    hn = np.array([float(value(K, t)) for t in h])
    names = list(old)
    proposals = []
    for anchor_index, first in enumerate(bank['caps']):
        r = decode(first, 'anchor')
        cn = np.array([float(value(r['field'], t)) for t in r['c']])
        for second_index, second in enumerate(bank['caps']):
            block = second['block_names']
            if len(set(block)&set(first['block_names'])) != 7:
                continue
            free, basis = nullspace_K([old[n]['c'] for n in block])
            N = np.array([[float(value(K, v[i])) for v in basis] for i in range(10)])
            aa = cn@N
            if abs(aa[1]) < 1e-9:
                N = N[:, ::-1]
                aa = aa[::-1]
            direction = N@np.array([1, -aa[0]/aa[1]])
            A = direction@In@direction
            outside = [i for i, n in enumerate(names) if n not in block]
            for label in (2, 3, 5):
                origin = N[:, 1]*(-2*np.cos(np.pi/label)/aa[1])
                B = 2*origin@In@direction
                D = B*B-4*A*(origin@In@origin-2)
                if D <= 1e-9 or abs(A) < 1e-9:
                    continue
                for branch in (-1, 1):
                    p = origin+(-B+branch*np.sqrt(D))/(2*A)*direction
                    pair = C@p
                    if p@In@hn <= 1e-8 or not np.all(pair[outside] < -np.sqrt(2*Qn[outside])-1e-7):
                        continue
                    proposals.append((anchor_index, second_index, label, branch))
    print({'numeric_coupled_proposals': len(proposals),
           'anchors': sorted(set(p[0] for p in proposals))}, flush=True)
    # Pick the anchor admitting the largest star of distinct block attachments.
    anchors = sorted(set(p[0] for p in proposals),
                     key=lambda i: (-len(set(p[1] for p in proposals if p[0] == i)), i))
    chosen_anchors = anchors[:2]
    pool = angle_pool(ANGLE_SOURCE)
    rows = []
    for anchor_index, second_index, label, branch in proposals:
        if anchor_index not in chosen_anchors:
            continue
        first, second = bank['caps'][anchor_index], bank['caps'][second_index]
        try:
            cap, record = exact(second['block_names'], first, label, branch, G, I, h, old)
        except AssertionError:
            continue
        checks = {n: relation(cap, r, pool) for n, r in old.items()}
        assert all(r == ({'kind': 'finite', 'label': 2} if n in second['block_names'] else {'kind': 'dotted'})
                   for n, r in checks.items())
        anchor = decode(first, 'anchor')
        assert relation(cap, anchor, pool) == {'kind': 'finite', 'label': label}
        record.update({'anchor_cap_index': anchor_index, 'second_block_index': second_index,
                       'pair_checks': checks})
        rows.append(record)
        print({'accepted': len(rows), 'anchor': anchor_index, 'block': second_index, 'label': label}, flush=True)
    out = {'source': source, 'source_case': case, 'anchor_source': PREFIX+'_vertex_caps.json',
           'numeric_proposals': [list(p) for p in proposals], 'chosen_anchors': chosen_anchors,
           'caps': rows, 'scope': 'Exact norm2 caps orthogonal to eight old walls and at a '
           'specified finite angle to a retained new vertex cap. Only the displayed finite '
           'rational anchor samples were screened; no continuous-family exclusion.'}
    (ROOT/(PREFIX+'_coupled_vertices.json')).write_text(json.dumps(out, indent=2)+'\n')
    print({'accepted': len(rows)}, flush=True)


if __name__ == '__main__':
    run()
