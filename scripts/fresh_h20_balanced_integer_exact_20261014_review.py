#!/usr/bin/env python3
"""Independent integer/Fraction replay of the frozen 810-wall chamber.

No producer, FLINT, NumPy, or project imports.  The exact ambient form and
type-A inverse formula are reconstructed here; denominator clearing makes
the full pair and endpoint checks ordinary integer arithmetic.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from math import lcm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_balanced_integer_exact_20261014_review'
BASE = 'fresh_h20_balanced_integer_exact_20261014'
LINKS = list(range(5)) + list(range(6, 11)) + list(range(12, 18))


def dual2(v):
    """2B v for diag(C_A6/2,C_A6/2,C_A7/2) plus DU=-1."""
    out = [2*x for x in v[:19]] + [-2*v[20], -2*v[19]]
    for i in LINKS:
        out[i] -= v[i+1]
        out[i+1] -= v[i]
    return out


def dot(a, b):
    return sum(x*y for x, y in zip(a, dual2(b))) / F(2)


def inverse_A(n, rhs):
    return [sum(F(2*min(i, j)*(n+1-max(i, j)), n+1)*rhs[j-1]
                for j in range(1, n+1)) for i in range(1, n+1)]


def rank(rows):
    a = [[F(x) for x in row] for row in rows]
    k = 0
    for j in range(len(a[0])):
        p = next((i for i in range(k, len(a)) if a[i][j]), None)
        if p is None:
            continue
        a[k], a[p] = a[p], a[k]
        q = a[k][j]
        a[k] = [x/q for x in a[k]]
        for i in range(k+1, len(a)):
            q = a[i][j]
            a[i] = [x-q*y for x, y in zip(a[i], a[k])]
        k += 1
        if k == len(a):
            break
    return k


def run():
    inp = ROOT/f'results/{BASE}.json'
    data = json.loads(inp.read_text())
    for p, h in data['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    for p in list(data['sources'])[:3]:
        doc = json.loads((ROOT/p).read_text())
        for q, h in doc['sources'].items():
            assert sha256((ROOT/q).read_bytes()).hexdigest() == h, q
    walls = [[F(x) for x in row] for row in data['walls']]
    assert len(walls) == 810 and all(len(r) == 21 for r in walls)

    seed = []
    start = 0
    for n in (6, 6, 7):
        for omitted in range(n+1):
            v = [F(0)]*21
            if omitted == 0:
                v[start:start+n] = [F(-1)]*n
            else:
                v[start+omitted-1] = F(1)
            v[19] = F(1)
            seed.append(v)
        start += n
    assert walls[:22] == seed
    assert walls[22] == [F(0)]*19+[F(-2), F(1)]
    vertices = []
    for a in range(7):
        for b in range(7):
            for c in range(8):
                spatial = []
                for n, o in zip((6, 6, 7), (a, b, c)):
                    spatial += inverse_A(n, [F(1)-(n+1 if i == o else 0)
                                             for i in range(1, n+1)])
                v = spatial+[F(0), F(0)]
                assert dot(v, v) == 196
                vertices.append(spatial+[F(112), F(1)])
                assert walls[22+len(vertices)] == spatial+[F(84), F(1)]
    motif = json.loads((ROOT/'results/fresh_h20_balanced_20261014_bridge_exact.json').read_text())
    for i, item in enumerate(motif['new_caps'], 415):
        c = list(map(F, item['center']))+[F(0), F(0)]
        assert walls[i] == c[:19]+[(dot(c, c)-F(item['radius'])**2)/2, F(1)]
    assert len(data['new_caps']) == 392
    for i, rec in enumerate(data['new_caps'], 418):
        gaps = list(map(F, rec['seed_gaps']))
        center = inverse_A(6, [1-x for x in gaps[1:7]])
        center += inverse_A(6, [1-x for x in gaps[8:14]])
        center += inverse_A(7, [1-x for x in gaps[15:22]])
        v = center+[F(0), F(0)]
        assert center == list(map(F, rec['center']))
        assert [1-dot(v, r) for r in seed] == gaps
        assert walls[i] == center+[(dot(v, v)-4)/2, F(1)]
        probe = [F(2, 5)*x for x in vertices[rec['probe']][:19]]+[F(0), F(0)]
        power = dot(probe, probe)+dot(v, v)-2*dot(probe, v)-4
        assert power == F(rec['probe_power']) < 0

    H = [F(0)]*19+[F(200), F(1)]
    X = list(map(F, data['strict_escape']['point']))
    assert X[:19] == [F(3, 5)*v for v in vertices[0][:19]]
    assert dot(H, H) == -400 and dot(X, X) == 1 and dot(H, X) < 0
    scale = lcm(*(x.denominator for r in walls+vertices+[H, X] for x in r))
    W = [[int(scale*x) for x in r] for r in walls]
    D = [dual2(r) for r in W]
    den = 2*scale**2
    def pair(a, d):
        return sum(x*y for x, y in zip(a, d))
    norms = [pair(r, d) for r, d in zip(W, D)]
    assert all(q > 0 for q in norms)
    hc = [pair([int(scale*x) for x in H], d) for d in D]
    xc = [pair([int(scale*x) for x in X], d) for d in D]
    assert max(hc) < 0 and max(xc) < 0
    counts = Counter()
    for i, r in enumerate(W):
        for j in range(i):
            p = pair(r, D[j])
            qq = norms[i]*norms[j]
            if p == 0:
                typ = 'right'
            elif p < 0 and 4*p*p == qq:
                typ = 'pi3'
            else:
                assert p < 0 and p*p >= qq, (i, j)
                typ = 'parallel' if p*p == qq else 'dotted'
            counts[typ] += 1
    assert dict(counts) == data['pair_types']
    # Strict facet projections follow immediately from these signs:
    # H_j - (H_i/q_i)G_ij < 0 and H^2-H_i^2/q_i < 0.
    for i, p in enumerate(vertices):
        assert dot(p, p) == -28 and dot(p, H) < 0
        pp = [pair([int(scale*x) for x in p], d) for d in D]
        zeros = [j for j, x in enumerate(pp) if not x]
        assert len(zeros) == 20 and max(pp) == 0
        assert 23+i in zeros and all(j < 22 or j == 23+i for j in zeros)
        # The remaining 19 form A6+A6+A7; the cap is orthogonal to all.
        assert all(pair(W[23+i], D[j]) == 0 for j in zeros if j != 23+i)

    delta = [walls[3][j]+walls[4][j]+walls[418][j]/2 for j in range(21)]
    assert dot(delta, delta) == 0 and dot(delta, H) < 0
    vals = [dot(delta, r) for r in walls]
    incident = [i for i, x in enumerate(vals) if x == 0]
    assert incident == [3, 4, 9, 12, 15, 17, 20, 418]
    assert max(x for x in vals if x != 0) == -F(1, 2)
    normalized = [walls[3], walls[4], [x/2 for x in walls[418]]]
    assert [[dot(a, b) for b in normalized] for a in normalized] == [
        [F(1), -F(1, 2), -F(1, 2)],
        [-F(1, 2), F(1), -F(1, 2)],
        [-F(1, 2), -F(1, 2), F(1)]]
    a1 = [9, 12, 15, 17, 20]
    for i in a1:
        assert dot(walls[i], walls[i]) == 1
        assert all(dot(walls[i], walls[j]) == 0 for j in incident if i != j)
    G = [[dot(walls[i], walls[j]) for j in incident] for i in incident]
    assert rank(G) == 7 and rank([walls[i] for i in incident]) == 8
    recession = [-x for x in walls[9]]
    assert dot(delta, recession) == 0 and dot(recession, recession) == 1
    assert all(dot(recession, walls[i]) <= 0 for i in incident)
    # All five spherical A1 components can be completed already at shift 100.
    # This only checks their individual compatibility with all old walls;
    # the existence lemma supplies the twelve residual orthogonal pairs.
    a1_mates = []
    for i in a1:
        r = [-x+100*y for x, y in zip(walls[i], delta)]
        assert dot(r, r) == 1 and dot(r, H) < 0
        for j, w in enumerate(walls):
            p = dot(r, w)
            if j in incident:
                assert p == (-1 if j == i else 0)
            else:
                assert p < 0 and p*p > dot(w, w)
        a1_mates.append(r)
    assert all(dot(a1_mates[i], a1_mates[j]) == 0
               for i in range(5) for j in range(i))
    sources = {f'scripts/{BASE}.py', f'results/{BASE}.json',
               f'scripts/{TAG}.py', *data['sources']}
    return dict(status='PASS_INDEPENDENT_810_GEOMETRY_AND_LOCAL_CUSP_DIAGNOSIS',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(sources)},
        arithmetic='stdlib integers and Fraction; no producer/project imports',
        walls=810, pairs=sum(counts.values()), pair_types=dict(counts),
        exact_root_reconstructions=810, retained_finite_endpoints=392,
        strict_escape_norm='1', strict_escape_max_pair=str(F(max(xc), den)),
        strict_center_max_pair=str(F(max(hc), den)),
        representative_cusp=dict(delta=list(map(str, delta)), incident=incident,
            components='affine A2 + 5 spherical A1', quotient_rank=7,
            missing_quotient_dimensions=12, strict_nonincident_max='-1/2',
            positive_kernel=[1, 1, 1], normalized_A2=[3, 4, '418/2'],
            positive_norm_recession='-wall9',
            five_A1_mates_checked_shift=100,
            local_repair='5 A1 mates and 12 orthogonal opposite pairs, with sufficiently large positive null shifts'),
        scope='Independent complete literal geometry replay; no complete cusp/vertex census. The representative cusp admits local repair, which makes no global finite-volume claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    output = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert output == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(output, indent=2)+'\n')
    print(output['status'], output['pairs'], output['representative_cusp']['components'])
