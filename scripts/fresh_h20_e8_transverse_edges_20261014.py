#!/usr/bin/env python3
"""Exact all-parameter escape gate for 280 specified transverse-edge cap families.

No numerical optimization or finite label cutoff occurs in this certificate.
Each cap is orthogonal to the 19 retained supports of one transverse edge.
This is not an obstruction to arbitrary E8-coupled caps or to completion.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_contact_closure_20261014_replacement_review as K

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e8_transverse_edges_20261014'
BASE = 'results/fresh_h20_e8_cusp_completion_20261014.json'
VERTICES = 'results/fresh_h20_contact_closure_20261014_h12_replacement.json'
Z = K.ZERO


def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def vec(row): return [tuple(map(F, x)) for x in row]
def sub(a, b): return K.add(a, K.neg(b))


def radical_sign(a, b, h):
    """Sign of a+b sqrt(h), h rational nonnegative, using exact squaring."""
    assert h >= 0
    sa, sb = K.sign(a), K.sign(b)
    if not h or not sb: return sa
    if not sa or sa == sb: return sb
    d = sub(K.mul(a, a), K.scale(K.mul(b, b), h))
    return sa*K.sign(d)


def enc(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, (tuple, list)): return [enc(y) for y in x]
    if isinstance(x, dict): return {k: enc(v) for k, v in x.items()}
    return x


def run():
    d = json.loads((ROOT/BASE).read_text())
    b = json.loads((ROOT/VERTICES).read_text())
    for doc in (d, b):
        for p, h in doc['sources'].items(): assert digest(p) == h, p
    Q = [[F(x) for x in row] for row in d['form']]
    sparse = [[(j, x) for j, x in enumerate(row) if x] for row in Q]
    cache = {}

    def dot(v, w):
        key = tuple(w)
        if key not in cache:
            cache[key] = [tuple(sum((a*w[j][k] for j, a in row), F(0)) for k in range(4)) for row in sparse]
        out = Z
        for x, y in zip(v, cache[key]):
            if x != Z and y != Z: out = K.add(out, K.mul(x, y))
        return out

    N = list(map(vec, d['normals'])); assert len(N) == 79
    E = vec(d['strict_escape']['point'])
    support = d['E8_support']
    inverse = [[F(x) for x in row] for row in b['E8_inverse']]
    weights = [[tuple(sum(inverse[i][j]*N[support[j]][a][k] for j in range(8)) for k in range(4)) for a in range(21)] for i in range(8)]
    assert all(dot(weights[i], N[support[j]]) == K.scalar(int(i == j)) for i in range(8) for j in range(8))
    A = [K.neg(dot(w, E)) for w in weights]
    wall_a = [[K.neg(dot(w, n)) for n in N] for w in weights]
    assert all(K.sign(a) > 0 for a in A)
    rows = []
    for vi, row in enumerate(b['H12_face_finite_vertices']):
        V = vec(row['point']); q = K.neg(dot(V, V))
        assert q == K.neg(tuple(map(F, row['norm']))) and q[1:] == (0, 0, 0) and q[0] > 0
        q = q[0]; beta = K.neg(dot(V, E)); assert K.sign(beta) > 0
        wall_b = [K.neg(dot(V, n)) for n in N]
        assert all(K.sign(x) >= 0 for x in wall_b)
        assert [i for i, x in enumerate(wall_b) if x == Z] == row['incident']
        for j in range(8):
            degree = inverse[j][j]
            assert dot(weights[j], weights[j]) == K.scalar(degree)
            assert dot(weights[j], V) == Z
            h1 = (degree-2)/q
            rec = dict(vertex=vi, E8_position=j, E8_wall=support[j], weight_norm=degree, vertex_negative_norm=q, escape_A=A[j], escape_beta=beta)
            # r/c = -omega_j + u V, u=sqrt((degree-2/c^2)/q).
            # Coxeter contact with e_j implies c=1 or c>=sqrt(2).
            # More weakly c>=1, and u is increasing in c.
            at1 = radical_sign(A[j], K.neg(beta), h1)
            if at1 <= 0:
                rec.update(reason='minimum_coefficient_escape_bound', h=h1, sign_at_c1=at1)
            else:
                blockers = []
                for i, (a, bb) in enumerate(zip(wall_a[j], wall_b)):
                    if K.sign(a) <= 0: continue
                    if bb == Z:
                        blockers.append((i, None)); break
                    gap = sub(K.mul(a, beta), K.mul(A[j], bb))
                    if K.sign(gap) >= 0:
                        blockers.append((i, gap)); break
                if blockers:
                    i, gap = blockers[0]
                    rec.update(reason='acute_wall_escape_bound', blocking_wall=i, wall_a=wall_a[j][i], wall_b=wall_b[i], cross_product_gap=gap)
                else:
                    # Everything c>=sqrt(2) is excluded at the left endpoint.
                    h4 = (degree-1)/q
                    at4 = radical_sign(A[j], K.neg(beta), h4)
                    assert at4 < 0, (vi, j)
                    rec.update(reason='Coxeter_coefficient_gap', h_at_sqrt2=h4, sign_at_sqrt2=at4)
                    # The only remaining possible c is exactly1.
                    blocked = False
                    for i, (a, bb) in enumerate(zip(wall_a[j], wall_b)):
                        sg = radical_sign(a, K.neg(bb), h1)
                        if sg > 0:
                            rec.update(c1_reason='nonacute', blocking_wall=i, wall_a=a, wall_b=bb, c1_h=h1)
                            blocked = True; break
                        if sg == 0: continue
                        # g=a-b sqrt(h1). Exclude 0 < -g/sqrt(2 norm_i) < 1/2.
                        square_a = K.add(K.mul(a, a), K.scale(K.mul(bb, bb), h1))
                        square_b = K.scale(K.mul(a, bb), -2)
                        threshold = F(d['norms'][i])/2
                        gapsg = radical_sign(sub(square_a, K.scalar(threshold)), square_b, h1)
                        if gapsg < 0:
                            rec.update(c1_reason='forbidden_angle_between_pi3_and_pi2', blocking_wall=i, wall_a=a, wall_b=bb, c1_h=h1, squared_pairing_threshold=threshold)
                            blocked = True; break
                    assert blocked, (vi, j)
            rows.append(rec)
    assert len(rows) == 280
    counts = dict(Counter(row['reason'] for row in rows))
    paths = [BASE, VERTICES, f'scripts/{TAG}.py', 'scripts/fresh_h20_contact_closure_20261014_replacement_review.py']
    return enc(dict(status='PASS_ALL_280_SPECIFIED_TRANSVERSE_EDGE_FAMILIES_RETAIN_ESCAPE', sources={p: digest(p) for p in paths}, family='r=-c omega_j+sqrt((c*c*d_j-2)/(-V*V)) V; norm2; positive square root; c>=1', counts=counts, records=rows, scope='Every cap in these 280 single-contact families that is Coxeter with all79 walls has r.escape<=0, for all real parameters and all Coxeter labels. Their arbitrary union cannot remove the saved strict positive-norm escape (boundary equality allowed). General nonorthogonal endpoint profiles, multi-node E8 contacts, and other H12 residual directions are outside this statement. No new finite-volume example or unrestricted noncompletion result.'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    result = run(); p = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(p.read_text())
    else: p.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['counts'])
