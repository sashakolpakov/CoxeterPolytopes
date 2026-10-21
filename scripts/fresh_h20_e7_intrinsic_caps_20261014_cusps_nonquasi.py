#!/usr/bin/env python3
"""Four-wall conditional non-quasi-arithmetic witness, independent Q(sqrt2)."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_intrinsic_caps_20261014_cusps_nonquasi'
SOURCE = 'results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
ZERO = (F(0), F(0))


def add(x, y): return (x[0]+y[0], x[1]+y[1])
def mul(x, y): return (x[0]*y[0]+2*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def scale(x, a): return (a*x[0], a*x[1])
def sub(x, y): return add(x, scale(y, -1))


def sign(x):
    a, b = x
    if not a: return (b > 0)-(b < 0)
    if not b or (a > 0) == (b > 0): return (a > 0)-(a < 0)
    n = a*a-2*b*b
    assert n
    return ((a > 0)-(a < 0))*((n > 0)-(n < 0))


def vec(rows):
    out = []
    for row in rows:
        assert len(row) == 8 and not any(map(F, row[2:]))
        out.append(tuple(map(F, row[:2])))
    return out


def digest(path): return sha256((ROOT/path).read_bytes()).hexdigest()


def encode(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, dict): return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [encode(v) for v in x]
    return x


def run():
    d = json.loads((ROOT/SOURCE).read_text())
    for path, value in d['sources'].items(): assert digest(path) == value, path
    assert len(d['normals']) == len(d['facet_points']) == 87
    ids = [19, 20, 22, 23]
    Q = [[F(v) for v in row] for row in d['form']]
    assert len(Q) == 21 and all(len(row) == 21 for row in Q)
    assert Q == list(map(list, zip(*Q)))
    sparse = [(i, j, q) for i, row in enumerate(Q) for j, q in enumerate(row) if q]

    def dot(x, y):
        out = ZERO
        for i, j, q in sparse: out = add(out, scale(mul(x[i], y[j]), q))
        return out

    normals = [vec(d['normals'][i]) for i in ids]
    H = vec(d['strict_center'])
    E = vec(d['strict_escape']['point'])
    assert sign(dot(H, H)) < 0 and sign(dot(E, H)) < 0
    assert all(dot(r, r) == (F(2), F(0)) for r in normals)
    assert all(F(d['norms'][i]) == 2 for i in ids)
    assert all(sign(dot(H, r)) < 0 and sign(dot(E, r)) < 0 for r in normals)
    G = [[dot(x, y) for y in normals] for x in normals]
    assert G[0][1] == (F(-2), F(-2))
    cycle = mul(G[0][1], G[1][0])
    assert cycle == (F(12), F(8))
    # All four roots have norm 2, so 2*(unit-normal Gram) is G.
    # Thus cycle is exactly the cyclic product in 2*(unit Gram).
    recovered_sqrt2 = scale(sub(cycle, (F(12), F(0))), F(1, 8))
    assert recovered_sqrt2 == (F(0), F(1))
    assert mul(recovered_sqrt2, recovered_sqrt2) == (F(2), F(0))
    block = [row[2:] for row in G[2:]]
    assert block == [[(F(2), F(0)), (F(-6), F(0))],
                     [(F(-6), F(0)), (F(2), F(0))]]
    determinant = sub(mul(block[0][0], block[1][1]), mul(block[0][1], block[1][0]))
    assert determinant == (F(-32), F(0))
    normalized = [[scale(v, F(1, 2)) for v in row] for row in block]
    assert normalized == [[(F(1), F(0)), (F(-3), F(0))],
                          [(F(-3), F(0)), (F(1), F(0))]]
    assert sub(mul(normalized[0][0], normalized[1][1]), mul(normalized[0][1], normalized[1][0])) == (F(-8), F(0))
    # Verify the saved facet-point recipes on the four-wall subsystem;
    # actuality against all 87 walls is inherited from SOURCE and its audit.
    for i, r in zip(ids, normals):
        p = vec(d['facet_points'][i])
        recipe = [sub(h, scale(v, F(1, 2))) for h, v in zip(H, [mul(dot(H, r), a) for a in r])]
        assert p == recipe and dot(p, r) == ZERO
        assert sign(dot(p, p)) < 0 and sign(dot(p, H)) < 0
        assert all(sign(dot(p, s)) < 0 for s in normals if s != r)
    en = dot(E, E)
    assert en == (F(1336028658157091827, 250000000000000000), F(0)) and sign(en) > 0
    paths = [SOURCE, f'scripts/{TAG}.py']
    return encode(dict(
        status='PASS_CONDITIONAL_NON_QUASI_ARITHMETIC_COMPLETION_WITNESS',
        sources={p: digest(p) for p in paths},
        arithmetic='Standalone Fraction arithmetic in basis [1,sqrt2]; no producer imports.',
        wall_indices=ids, four_wall_gram=G, norms=[F(2)]*4,
        irrational_cycle=dict(walls=[19, 20], pairing=G[0][1], cartan_two_cycle=cycle,
                              sqrt2_recovery='(cycle-12)/8'),
        invariant_indefinite_block=dict(walls=[22, 23], raw_gram=block,
                                       raw_determinant=determinant,
                                       unit_gram=normalized, unit_determinant=F(-8)),
        current_strict_escape_norm=en,
        criterion_url='https://arxiv.org/html/2607.14715v1#S3.SS1',
        scope='Every finite-volume Coxeter completion retaining these four walls as actual facets is non-quasi-arithmetic. The current 87-wall chamber has a certified strict spacelike escape and is not a lattice. No finite-volume existence claim.'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(path.read_text())
    else: path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
