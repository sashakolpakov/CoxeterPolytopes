#!/usr/bin/env python3
"""Independent rational Gram-matrix review of the418-wall bridge chamber."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from flint import fmpq as Q, fmpq_mat as M

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_balanced_bridge_review_20261014'
SOURCE = 'results/fresh_h20_balanced_20261014_bridge_exact.json'


def run():
    src = json.loads((ROOT/SOURCE).read_text())
    for p, h in src['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest() == h
    B = M(21, 21)
    seed, factor_vertices = [], []
    start = 0
    for n in (6, 6, 7):
        for i in range(n):
            B[start+i, start+i] = 1
            if i+1 < n:
                B[start+i, start+i+1] = B[start+i+1, start+i] = Q(-1, 2)
        C = M([[B[start+i, start+j] for j in range(n)] for i in range(n)])
        inv = C.inv()
        rows = [[Q(-1)]*n]+[[Q(i == j) for j in range(n)] for i in range(n)]
        seed += [[Q(0)]*start+r+[Q(0)]*(19-start-n)+[Q(1), Q(0)] for r in rows]
        fv = []
        for omit in range(n+1):
            rhs = M([[1-(n+1 if omit == i+1 else 0)] for i in range(n)])
            sol = inv*rhs
            fv.append([Q(0)]*start+[sol[i, 0] for i in range(n)]+[Q(0)]*(19-start-n))
        factor_vertices.append(fv)
        start += n
    B[19, 20] = B[20, 19] = -1
    vertices = [[sum(v[i] for v in chosen) for i in range(19)] for chosen in product(*factor_vertices)]
    def cap(c, r2):
        v = M([c+[Q(0), Q(0)]])
        norm = (v*B*v.transpose())[0, 0]
        return c+[(norm-r2)/2, Q(1)]
    roots = seed+[[Q(0)]*19+[Q(-2), Q(1)]]
    roots += [cap(v, Q(28)) for v in vertices]
    roots += [cap(list(map(Q, r['center'])), Q(r['radius'])**2) for r in src['new_caps']]
    W = M(roots)
    G = W*B*W.transpose()
    assert len(roots) == 418
    counts = Counter()
    for i in range(418):
        assert G[i, i] > 0
        for j in range(i):
            g, den2 = G[i, j], G[i, i]*G[j, j]
            if not g:
                typ = 'right'
            elif g < 0 and 4*g*g == den2:
                typ = 'pi3'
            else:
                assert g < 0 and g*g >= den2, (i, j, g, den2)
                typ = 'parallel' if g*g == den2 else 'dotted'
            counts[typ] += 1
    H = M([[Q(0)]*19+[Q(200), Q(1)]])
    HW = H*B*W.transpose()
    assert (H*B*H.transpose())[0, 0] == -400
    assert all(HW[0, i] < 0 for i in range(418))
    # Projection ofH to every wall yields a timelike point strict in all others.
    for i in range(418):
        assert -400-HW[0, i]**2/G[i, i] < 0
        assert all(HW[0, j]-HW[0, i]*G[i, j]/G[i, i] < 0 for j in range(418) if i != j)
    P = M([v+[Q(112), Q(1)] for v in vertices])
    PP, PW = P*B*P.transpose(), P*B*W.transpose()
    for i in range(392):
        assert PP[i, i] == -28
        zeros = [j for j in range(418) if PW[i, j] == 0]
        assert len(zeros) == 20 and 23+i in zeros
        assert all(PW[i, j] <= 0 for j in range(418))
        assert all(j < 22 or j == 23+i for j in zeros)
    e = src['strict_escape']
    X = M([list(map(Q, e['spatial']))+[Q(e['D']), Q(e['U'])]])
    assert (X*B*X.transpose())[0, 0] == 1
    assert (X*B*H.transpose())[0, 0] < 0
    XW = X*B*W.transpose()
    assert all(XW[0, j] < 0 for j in range(418))
    paths = [SOURCE, str(Path(__file__).resolve().relative_to(ROOT))]
    return dict(status='PASS_INDEPENDENT_418_WALL_BRIDGE_MATRIX_REVIEW',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                wall_count=418, pair_count=418*417//2, pair_types=dict(counts),
                actual_facet_witnesses=418, finite_vertices=392, strict_escape_norm=1,
                scope='Separate FLINT matrix reconstruction of all418 roots, all pairs, every actual-facet projection,392 finite endpoints and one strict spacelike escape. Confirms infinite volume. Does not independently check the radius-five general theorem or numerical scouts.')


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
    print(out['status'], out['pair_types'], flush=True)
