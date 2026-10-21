#!/usr/bin/env python3
"""Standalone rational replay of saved physical high-S interval certificates.

No project imports, numerical root solving, native-bank or cone reconstruction.
The same Schur identities are implemented separately, with directed rational
rounding. This is not an alternative interval algorithm or scaffold proof.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_oblique_high_s_polar_interval_20261018'
PRIMARY = 'fresh_nonarithmetic_release49_62_replacement_20261017'
HIGH = 'fresh_oblique_high_n_right_r_join_20261018'
WHITE, OLD = [0, 1, 2, 3, 5, 6], [8, 10, 52, 63, 64]
SCALE = 10**70


def iv(x):
    return tuple(map(F, x)) if isinstance(x, (list, tuple)) else (F(x), F(x))


def roundout(v):
    a, b = v
    assert a <= b
    return F(a.numerator*SCALE//a.denominator, SCALE), F(-((-b.numerator*SCALE)//b.denominator), SCALE)


def add(a, b): return roundout((a[0]+b[0], a[1]+b[1]))
def neg(a): return -a[1], -a[0]
def sub(a, b): return add(a, neg(b))
def mul(a, b):
    z = [x*y for x in a for y in b]
    return roundout((min(z), max(z)))
def div(a, b):
    assert b[0]*b[1] > 0
    return mul(a, (1/b[1], 1/b[0]))
def total(values):
    out = iv(0)
    for v in values: out = add(out, v)
    return out
def sqrt(a):
    assert a[0] >= 0
    lo = isqrt(a[0].numerator*SCALE*SCALE//a[0].denominator)
    hi = isqrt(a[1].numerator*SCALE*SCALE//a[1].denominator)
    if F(hi, SCALE)**2 < a[1]: hi += 1
    return F(lo, SCALE), F(hi, SCALE)
def dot(x, g, y): return total(mul(mul(a, g[i][j]), b) for i, a in enumerate(x) for j, b in enumerate(y))
def enc(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, dict): return {str(k): enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [enc(v) for v in x]
    return x
def determinant(g):
    n, terms = len(g), []
    for p in permutations(range(n)):
        value = iv((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n)))
        for i in range(n): value = mul(value, g[i][p[i]])
        terms.append(value)
    return total(terms)
def inverse(g):
    n, d = len(g), determinant(g)
    assert d[0]*d[1] > 0
    return [[div(mul(iv((-1)**(i+j)), determinant([[g[a][b] for b in range(n) if b != i]
        for a in range(n) if a != j])), d) for j in range(n)] for i in range(n)], d
def gauss(g):
    n = len(g)
    a = [r[:]+[iv(int(i == j)) for j in range(n)] for i, r in enumerate(g)]
    pivots = []
    for i in range(n):
        p = a[i][i]
        assert p[0]*p[1] > 0
        pivots.append(p)
        a[i] = [div(x, p) for x in a[i]]
        for j in range(n):
            if i != j:
                q = a[j][i]
                a[j] = [sub(x, mul(q, y)) for x, y in zip(a[j], a[i])]
    return [r[n:] for r in a], pivots
def theta(m):
    if m in (2, 3): return iv(m-2)
    if m == 4: return sqrt(iv(2))
    if m == 5: return div(add(iv(1), sqrt(iv(5))), iv(2))
    if m == 6: return sqrt(iv(3))
    if m == 8: return sqrt(add(iv(2), sqrt(iv(2))))
    if m == 10: return sqrt(div(add(iv(5), sqrt(iv(5))), iv(2)))
    if m == 12: return div(add(sqrt(iv(2)), sqrt(iv(6))), iv(2))
    assert m == 7
    def atan(x):
        z = sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(100))
        return z, z+x**201/F(201)
    pi = sub(mul(iv(16), atan(F(1, 5))), mul(iv(4), atan(F(1, 239))))
    x = div(pi, iv(m))
    xx, term, value = mul(x, x), iv(1), iv(1)
    for k in range(1, 60):
        term = div(mul(neg(xx), term), iv((2*k-1)*(2*k)))
        value = add(value, term)
    e = xx[1]*max(abs(t) for t in term)/F(119*120)
    return mul(iv(2), (value[0]-e, value[1]+e))


def check():
    d = json.loads((ROOT/f'results/{TAG}.json').read_text())
    for p, h in d['sources'].items(): assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    pri = json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    old = json.loads((ROOT/f'results/{HIGH}.json').read_text())
    assert d['status'] == 'PASS_FIXED_HIGH_N_HIGH_S_HIGH_RM_PHYSICAL_INTERVAL_OBSTRUCTION'
    assert d['interval_scale'] == SCALE and len(d['parents']) == 2 and len(d['records']) == 20
    ids = WHITE+OLD
    metric, center = [list(map(iv, r)) for r in pri['computational_gram']], list(map(iv, pri['common_center']))
    roots = {}
    for row in pri['facet_witnesses']:
        if row['wall'] in ids:
            hp = iv(pri['center_products'][str(row['wall'])])
            assert hp[1] < 0
            roots[row['wall']] = [div(mul(iv(2), sub(x, y)), hp) for x, y in zip(center, map(iv, row['point']))]
    native = {a: {b: dot(roots[a], metric, roots[b]) for b in ids} for a in ids}
    wi, pivots = gauss([[native[a][b] for b in WHITE] for a in WHITE])
    assert all(p[0] > 0 for p in pivots)
    wc0 = {a: [native[a][w] for w in WHITE] for a in OLD}
    g0 = {a: {b: sub(native[a][b], dot(wc0[a], wi, wc0[b])) for b in OLD} for a in OLD}
    assert d['primary_projection'] == enc(dict(walls=ids, physical_Gram=native, white_inverse=wi,
        white_pivots=pivots, old_white_covectors=wc0, projected_Gram=g0))
    profiles = [[iv(0)]*6]
    for pos in (4, 5):
        for m in (3, 4, 5):
            v = [iv(0)]*6
            v[pos] = neg(theta(m))
            profiles.append(v)
    assert enc(profiles) == d['profiles']

    def extend(g, wc, row):
        wall, frame, pi, labels = (row[k] for k in ('wall', 'frame', 'profile', 'labels'))
        gram = [[g[a][b] for b in frame] for a in frame]
        inv, det = inverse(gram)
        assert det[1] < 0 and gram[0][0][0] > 0
        rhs = [sub(neg(theta(m)), dot(profiles[pi], wi, wc[i])) for i, m in zip(frame, labels)]
        square = sub(iv(2), dot(profiles[pi], wi, profiles[pi]))
        rem = sub(square, dot(rhs, inv, rhs))
        w, sign = row['witness'], row['witness_sign']
        assert w in g and w not in frame and sign in (-1, 1)
        cv = [g[i][w] for i in frame]
        q = sub(g[w][w], dot(cv, inv, cv))
        assert q[0] > 0 and rem[0] > 0
        linear, radical = dot(rhs, inv, cv), sqrt(mul(rem, q))
        branches = [sub(linear, radical), add(linear, radical)]
        factor = mul(iv(sign), sqrt(div(rem, q)))
        out = {}
        for j in g:
            cj = [g[i][j] for i in frame]
            out[j] = rhs[frame.index(j)] if j in frame else add(dot(rhs, inv, cj),
                mul(factor, sub(g[w][j], dot(cv, inv, cj))))
        expected = dict(frame_Gram=gram, inverse=inv, determinant=det, rhs=rhs,
            projected_square=square, remainder=rem, witness_residual=q,
            branch_witness_products=branches, new_products=out)
        assert all(row[k] == enc(v) for k, v in expected.items())
        distances = [abs(row['numerical_witness_product']-float(sum(v)/2)) for v in branches]
        assert distances == row['branch_comparison'] and distances[(sign+1)//2] < 2e-6 and distances[(1-sign)//2] > 1e-4
        assert row['all_product_comparison_error'] < 2e-6
        for j, v in out.items(): g[j][wall] = v
        g[wall], wc[wall] = {**out, wall: square}, profiles[pi]

    bounds, gaps, old8 = [], [], []
    index = 0
    for parent, source in zip(d['parents'], old['records']):
        assert parent['index'] == source['index'] and parent['origins'] == source['origins']
        assert parent['instruction_key'] == source['instruction_key']
        chain = json.loads(source['instruction_key'])
        assert len(chain) == len(parent['steps']) == 5 and [z[0] for z in chain] == [67, 68, 69, 70, 71]
        g, wc = {a: r.copy() for a, r in g0.items()}, wc0.copy()
        for step, instruction in zip(parent['steps'], chain):
            assert [step[k] for k in ('wall', 'frame', 'profile', 'labels', 'saved_numerical_branch')] == instruction
            extend(g, wc, step)
        assert len(source['caps']) == 10
        for cap in source['caps']:
            row = d['records'][index]
            assert (row['index'], row['parent'], row['N_index'], row['N_branch']) == (index, source['index'], cap['index'], cap['branch'])
            index += 1
            assert [z for z in cap['S_bank']['symbolic_high_families'] if z[1][2] == 2] == [[0, [6, 2, 2], [0]]]
            g1, wc1 = {a: r.copy() for a, r in g.items()}, wc.copy()
            assert [row['N_extension'][k] for k in ('wall', 'frame', 'profile', 'labels')] == [73, [10, 67, 69], 0, cap['labels']]
            extend(g1, wc1, row['N_extension'])
            c = row['continuous']
            frame = [52, 68, 71]
            inv, det = inverse([[g1[a][b] for b in frame] for a in frame])
            rhs, cv = [iv(-1), iv(0), iv(0)], [g1[i][8] for i in frame]
            alpha = neg(dot(rhs, inv, rhs))
            q8 = sub(g1[8][8], dot(cv, inv, cv))
            assert det[1] < 0 and alpha[0] > 0 and q8[0] > 0
            aa = {j: dot(rhs, inv, [g1[i][j] for i in frame]) for j in g1}
            bb = {j: div(sub(g1[8][j], dot(cv, inv, [g1[i][j] for i in frame])), sqrt(q8)) for j in g1}
            z = (div(theta(6), sqrt(add(iv(2), mul(iv(3), alpha))))[0], div(iv(2), sqrt(add(iv(2), mul(iv(4), alpha))))[1])
            positive8 = add(mul(z, aa[8]), bb[8])
            xn = [69, 73]
            gram = [[g1[a][b] for b in xn] for a in xn]
            xi, xd = inverse(gram)
            assert gram[0][0][0] > 0 and xd[0] > 0 and positive8[0] > 0
            mm, ar, br = [g1[a][71] for a in xn], [aa[a] for a in xn], [neg(bb[a]) for a in xn]
            m = sub(g1[71][71], dot(mm, xi, mm))
            q = [sub(iv(1), dot(br, xi, br)), neg(mul(iv(2), dot(ar, xi, br))), sub(neg(alpha), dot(ar, xi, ar))]
            cross = [sub(neg(bb[71]), dot(br, xi, mm)), sub(aa[71], dot(ar, xi, mm))]
            factor = sub(iv(3), mul(iv(2), m))
            h = [sub(neg(mul(factor, q[0])), mul(iv(2), mul(cross[0], cross[0]))),
                 sub(neg(mul(factor, q[1])), mul(iv(4), mul(cross[0], cross[1]))),
                 sub(neg(mul(factor, q[2])), mul(iv(2), mul(cross[1], cross[1])))]
            evaluate = lambda p: add(p[0], mul(z, add(p[1], mul(z, p[2]))))
            qr, hr = evaluate(q), evaluate(h)
            assert qr[1] < 0 and hr[0] > 0
            bound = F(3)-hr[0]/(-qr[0])
            assert 0 < bound < 3
            expected = dict(N_labels=cap['labels'], S_frame=frame, S_frame_inverse=inv,
                S_frame_determinant=det, alpha=alpha, residual8_square=q8, A_products=aa, B_products=bb,
                theta_domain=[theta(6)[0], F(2)], transformed_z_domain=z, positive_S_sheet_old8=positive8,
                XN_Gram=gram, XN_inverse=xi, XN_determinant=xd, M_residual_square=m,
                negative_S_residual_square_coefficients=q, M_cross_coefficients=cross,
                strict_high_RM_obstruction_coefficients=h, S_residual_square_range=qr,
                high_RM_obstruction_range=hr, forced_RM_squared_upper=bound)
            assert c == enc(expected)
            bounds.append(bound); gaps.append(hr[0]); old8.append(positive8[0])
    assert index == 20
    summary = dict(parent_classes=2, original_parent_origins=4, parent_extensions=10, N_geometries=20,
        continuous_S_families=20, positive_old8_S_sheets=20, polar_square_exclusions=20,
        minimum_old8_normalized_positive=min(old8), minimum_polynomial_margin=min(gaps),
        maximum_forced_RM_square=max(bounds), unresolved=0, accepted=0, new_cones=0)
    assert d['summary'] == enc(summary)
    assert d['exact_deletion_debit'] == old['exact_deletion_debit'] and d['deleted_walls'] == old['deleted_walls']
    assert d['topology'] == old['topology']
    assert d['necessary_conditions'] == dict(N_profile=0, S_profile=0, R_profile=0,
        NX_label_at_least=6, S52_label_at_least=6, RM_label_at_least=6, SL=2, SM=2, RX=2, RN=2, RS=2)
    print('PASS_STANDALONE_PHYSICAL_HIGH_S_INTERVAL_SAVED', 'extensions30', 'continuous20',
          'maximum_RM_square', float(max(bounds)), 'new_cones0', flush=True)


if __name__ == '__main__':
    check()
