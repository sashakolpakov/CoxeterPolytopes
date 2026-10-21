#!/usr/bin/env python3
"""Standalone rational replay of the complete low-R cube / continuous high-S proof.

No project imports, numerical root solving, native-bank or cone reconstruction.
The same Schur identities are implemented separately, with directed rational
rounding. This is not an alternative interval algorithm or scaffold proof.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations, product
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_oblique_high_s_complete_low_r_interval_20261018'
PRIOR = 'fresh_oblique_high_s_polar_interval_20261018'
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
    for p, h in d['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    pri = json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    old = json.loads((ROOT/f'results/{HIGH}.json').read_text())
    prior = json.loads((ROOT/f'results/{PRIOR}.json').read_text())
    assert d['status'] == 'PASS_COMPLETE_FIXED_LOW_R_CUBE_CONTINUOUS_HIGH_S_INTERVAL_EXCLUSION'
    assert d['interval_scale'] == SCALE and len(d['records']) == 20
    assert len(prior['parents']) == len(old['records']) == 2 and len(prior['records']) == 20
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
    assert d['inherited_physical_projection'] == enc(dict(walls=ids, physical_Gram=native, white_inverse=wi,
        white_pivots=pivots, old_white_covectors=wc0, projected_Gram=g0))
    profiles = [[iv(0)]*6]
    for pos in (4, 5):
        for m in (3, 4, 5):
            v = [iv(0)]*6
            v[pos] = neg(theta(m))
            profiles.append(v)
    assert enc(profiles) == d['inherited_profiles']

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

    remainders, positives, bounds, derivatives, keys = [], [], [], [], []
    for parent, source in zip(prior['parents'], old['records']):
        assert parent['index'] == source['index'] and parent['origins'] == source['origins']
        assert parent['instruction_key'] == source['instruction_key']
        chain = json.loads(source['instruction_key'])
        assert len(chain) == len(parent['steps']) == 5
        assert [z[0] for z in chain] == [67, 68, 69, 70, 71]
        g, wc = {a: r.copy() for a, r in g0.items()}, wc0.copy()
        for step, instruction in zip(parent['steps'], chain):
            assert [step[k] for k in ('wall', 'frame', 'profile', 'labels', 'saved_numerical_branch')] == instruction
            extend(g, wc, step)
        assert len(source['caps']) == 10
        for cap in source['caps']:
            match = [r for r in prior['records'] if (r['parent'], r['N_index']) == (source['index'], cap['index'])]
            assert len(match) == 1
            inherited = match[0]
            rec = d['records'][len(keys)]
            key = [source['index'], cap['index']]
            assert (rec['index'], rec['parent'], rec['N_index'], rec['N_labels']) == (len(keys), *key, cap['labels'])
            keys.append(key)
            assert [r for r in cap['S_bank']['symbolic_high_families'] if r[1][2] == 2] == [[0, [6, 2, 2], [0]]]
            gn, wn = {a: r.copy() for a, r in g.items()}, wc.copy()
            nstep = inherited['N_extension']
            assert [nstep[k] for k in ('wall', 'frame', 'profile', 'labels')] == [73, [10, 67, 69], 0, cap['labels']]
            extend(gn, wn, nstep)
            del gn[70], wn[70]
            for row in gn.values(): del row[70]
            # Independently certify the positive S sheet before considering R.
            sf = [52, 68, 71]
            si, sd = inverse([[gn[a][b] for b in sf] for a in sf])
            sr, s8 = [iv(-1), iv(0), iv(0)], [gn[i][8] for i in sf]
            alpha = neg(dot(sr, si, sr))
            q8 = sub(gn[8][8], dot(s8, si, s8))
            assert sd[1] < 0 and alpha[0] > 0 and q8[0] > 0
            z = (div(theta(6), sqrt(add(iv(2), mul(iv(3), alpha))))[0],
                 div(iv(2), sqrt(add(iv(2), mul(iv(4), alpha))))[1])
            positiveS8 = add(mul(z, dot(sr, si, s8)), div(q8, sqrt(q8)))
            assert positiveS8[0] > 0
            assert rec['positive_S_sheet_old8'] == inherited['continuous']['positive_S_sheet_old8']
            assert iv(rec['positive_S_sheet_old8'])[0] > 0
            # No numerical R selection: all28 recipes and both exact signs.
            frame, witness = [69, 71, 73], 10
            gram = [[gn[a][b] for b in frame] for a in frame]
            inv, det = inverse(gram)
            cv = [gn[i][witness] for i in frame]
            sigma = sub(gn[witness][witness], dot(cv, inv, cv))
            assert det[1] < 0 and gram[0][0][0] > 0 and sigma[0] > 0
            expected = dict(R_frame=frame, frame_Gram=gram, inverse=inv, determinant=det,
                witness=witness, witness_covector=cv, witness_residual=sigma)
            assert all(rec[k] == enc(v) for k, v in expected.items())
            assert len(rec['recipes']) == 28
            assert [(r['profile'], r['RM_label']) for r in rec['recipes']] == list(product(range(7), range(2, 6)))
            for index, (profile, label) in enumerate(product(range(7), range(2, 6))):
                row = rec['recipes'][index]
                labels = [2, label, 2]
                rhs = [sub(neg(theta(m)), dot(profiles[profile], wi, wn[i])) for i, m in zip(frame, labels)]
                square = sub(iv(2), dot(profiles[profile], wi, profiles[profile]))
                rho = sub(square, dot(rhs, inv, rhs))
                assert rho[0] > 0
                factor = sqrt(div(rho, sigma))
                sheets = {}
                for sign in (-1, 1):
                    contacts = {}
                    for j in gn:
                        cj = [gn[i][j] for i in frame]
                        contacts[j] = rhs[frame.index(j)] if j in frame else add(dot(rhs, inv, cj),
                            mul(mul(iv(sign), factor), sub(gn[witness][j], dot(cv, inv, cj))))
                    sheets[sign] = contacts
                positive8 = add(sheets[1][8], dot(profiles[profile], wi, wn[8]))
                assert positive8[0] > 0
                cc = [sheets[-1][i] for i in sf]
                a = dot(sr, si, cc)
                b = div(sub(sheets[-1][8], dot(s8, si, cc)), sqrt(q8))
                derivative = sub(a, mul(mul(alpha, z), b))
                assert derivative[1] < 0
                def evaluate(t):
                    return sub(mul(t, a), mul(sqrt(add(iv(2), mul(alpha, mul(t, t)))), b))
                at6, limit = evaluate(theta(6)), evaluate(iv(2))
                product_range = limit[0], at6[1]
                assert product_range[1] < -2
                continuous = dict(S_frame=sf, S_frame_inverse=si, S_frame_determinant=sd,
                    alpha=alpha, residual8_square=q8, A_dot_R=a, B_dot_R=b,
                    theta_domain=[theta(6)[0], F(2)], transformed_z_domain=z,
                    derivative_range=derivative, at_label6=at6, at_limit=limit,
                    full_product_range=product_range, monotonicity='strictly_decreasing')
                expected = dict(index=index, profile=profile, RM_label=label, R_labels=labels,
                    rhs=rhs, projected_square=square, norm_remainder=rho,
                    positive_R_sheet_old8=positive8,
                    negative_R_sheet_products={i: sheets[-1][i] for i in (8, 52, 68, 71)},
                    continuous=continuous)
                assert row == enc(expected)
                remainders.append(rho[0]); positives.append(positive8[0])
                bounds.append(product_range[1]); derivatives.append(derivative[1])
    assert len(keys) == len({tuple(k) for k in keys}) == 20 and len(bounds) == 560
    summary = dict(parent_classes=2, parent_origins=4, N_geometries=20, complete_R_recipes=560,
        exact_R_sheets=1120, positive_norm_remainders=560, excluded_positive_S_sheets=20,
        positive_old8_R_sheets=560, continuously_dotted_negative_R_sheets=560,
        minimum_norm_remainder=min(remainders), minimum_positive_R_old8=min(positives),
        maximum_product_upper=max(bounds), maximum_derivative_upper=max(derivatives),
        unresolved=0, accepted=0, new_cones=0)
    assert d['summary'] == enc(summary)
    assert d['exact_deletion_debit'] == old['exact_deletion_debit'] and d['deleted_walls'] == old['deleted_walls']
    assert d['exact_deletion_debit']['lost_finite_count'] == 19
    assert d['exact_deletion_debit']['compact_faces_by_dimension'] == {'2': 12, '3': 1}
    assert d['topology'] == old['topology']
    assert d['necessary_conditions'] == dict(N_profile=0, NX_label_at_least=6,
        S_profile=0, S52_label_at_least=6, SL=2, SM=2, RX=2, RN=2,
        R_profile_domain=list(range(7)), RM_labels=[2, 3, 4, 5])
    print('PASS_STANDALONE_COMPLETE_LOW_R_HIGH_S_INTERVAL', 'parent_N_extensions30',
          'recipes560', 'R_sheets1120', 'maximum_RS', float(max(bounds)), 'new_cones0', flush=True)


if __name__ == '__main__':
    check()
