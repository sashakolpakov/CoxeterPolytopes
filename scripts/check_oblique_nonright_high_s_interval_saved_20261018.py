#!/usr/bin/env python3
"""Standalone rational replay of eight nonright-SM symbolic high-S families.

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
TAG = 'fresh_oblique_nonright_high_s_interval_20261018'
CUBE = 'fresh_oblique_high_s_complete_low_r_interval_20261018'
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
    for p, h in d['sources'].items(): assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    pri = json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    old = json.loads((ROOT/f'results/{HIGH}.json').read_text())
    prior = json.loads((ROOT/f'results/{PRIOR}.json').read_text())
    cube = json.loads((ROOT/f'results/{CUBE}.json').read_text())
    assert d['status'] == 'PASS_EIGHT_NONRIGHT_HIGH_S_FAMILIES_PHYSICAL_INTERVAL_EXCLUSION'
    assert d['interval_scale'] == SCALE and len(d['records']) == 20
    assert cube['summary']['exact_R_sheets'] == 1120 and cube['summary']['positive_old8_R_sheets'] == 560
    assert len(old['records']) == len(prior['parents']) == 2
    assert len(cube['records']) == len(prior['records']) == 20
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
    assert prior['primary_projection'] == enc(dict(walls=ids, physical_Gram=native, white_inverse=wi,
        white_pivots=pivots, old_white_covectors=wc0, projected_Gram=g0))
    profiles = [[iv(0)]*6]
    for pos in (4, 5):
        for m in (3, 4, 5):
            v = [iv(0)]*6
            v[pos] = neg(theta(m))
            profiles.append(v)
    assert enc(profiles) == prior['profiles']

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

    families = [[0, [k, 2, 6], [2]] for k in range(2, 6)] + [[0, [6, 2, k], [0]] for k in range(3, 6)] + [[0, [6, 2, 6], [0, 2]]]
    assert d['symbolic_families'] == families
    keys, rho_lowers, positive_lowers, upper_products = [], [], [], []
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
            key = source['index'], cap['index']
            pm = [r for r in prior['records'] if (r['parent'], r['N_index']) == key]
            cm = [r for r in cube['records'] if (r['parent'], r['N_index']) == key]
            assert len(pm) == len(cm) == 1
            row = d['records'][len(keys)]
            assert (row['index'], row['parent'], row['N_index'], row['N_labels']) == (len(keys), *key, cap['labels'])
            keys.append(key)
            assert [f for f in cap['S_bank']['symbolic_high_families'] if f[1] != [6, 2, 2]] == families
            gn, wn = {a: r.copy() for a, r in g.items()}, wc.copy()
            nstep = pm[0]['N_extension']
            assert [nstep[k] for k in ('wall', 'frame', 'profile', 'labels')] == [73, [10, 67, 69], 0, cap['labels']]
            extend(gn, wn, nstep)
            frame = [52, 68, 71]
            gram = [[gn[a][b] for b in frame] for a in frame]
            inv, det = inverse(gram)
            av, cv, old8 = [iv(-1), iv(0), iv(0)], [iv(0), iv(0), iv(-1)], [gn[i][8] for i in frame]
            sigma = sub(gn[8][8], dot(old8, inv, old8))
            assert det[1] < 0 and gram[0][0][0] > 0 and sigma[0] > 0
            aa, ac, cc = dot(av, inv, av), dot(av, inv, cv), dot(cv, inv, cv)
            a8, c8 = dot(av, inv, old8), dot(cv, inv, old8)
            basis = dict(frame=frame, Gram=gram, inverse=inv, determinant=det,
                old8_covector=old8, old8_square=gn[8][8], residual8_square=sigma,
                AA=aa, AC=ac, CC=cc, A8=a8, C8=c8)
            assert row['S_basis'] == enc(basis)
            assert len(cm[0]['recipes']) == len(row['R_coefficients']) == 28
            assert [(r['profile'], r['RM_label']) for r in cm[0]['recipes']] == list(product(range(7), range(2, 6)))
            coeffs = []
            for index, r in enumerate(cm[0]['recipes']):
                assert r['index'] == index and iv(r['positive_R_sheet_old8'])[0] > 0
                contacts = [iv(r['negative_R_sheet_products'][str(i)]) for i in frame]
                ar, cr = dot(av, inv, contacts), dot(cv, inv, contacts)
                br = div(sub(iv(r['negative_R_sheet_products']['8']), dot(old8, inv, contacts)), sqrt(sigma))
                assert enc(ar) == r['continuous']['A_dot_R'] and enc(br) == r['continuous']['B_dot_R']
                coeffs.append(dict(index=index, profile=r['profile'], RM_label=r['RM_label'], AR=ar, CR=cr, BR=br))
            assert row['R_coefficients'] == enc(coeffs)
            assert len(row['families']) == 8
            for family, saved in zip(families, row['families']):
                profile, labels, high = family
                theta_domain = (theta(6)[0], F(2)) if 0 in high else theta(labels[0])
                tau_domain = (theta(6)[0], F(2)) if 2 in high else theta(labels[2])
                th, ta = theta_domain, tau_domain
                rho = sub(iv(2), add(add(mul(aa, mul(th, th)), mul(iv(2), mul(ac, mul(th, ta)))), mul(cc, mul(ta, ta))))
                assert rho[0] > 0
                radical = sqrt(rho)
                pos8 = add(add(mul(th, a8), mul(ta, c8)), mul(radical, sqrt(sigma)))
                assert pos8[0] > 0
                joins = []
                for r in coeffs:
                    value = sub(add(mul(th, r['AR']), mul(ta, r['CR'])), mul(radical, r['BR']))
                    assert value[1] < -2
                    joins.append(dict(R_index=r['index'], full_SR_range=value))
                    upper_products.append(value[1])
                expected = dict(family=family, theta_domain=th, tau_domain=ta,
                    norm_remainder_range=rho, positive_S_sheet_old8=pos8, negative_S_R_joins=joins)
                assert saved == enc(expected)
                rho_lowers.append(rho[0]); positive_lowers.append(pos8[0])
    assert len(keys) == len(set(keys)) == 20 and len(rho_lowers) == 160 and len(upper_products) == 4480
    summary = dict(parent_classes=2, parent_origins=4, N_geometries=20, symbolic_family_types=8,
        continuous_S_families=160, exact_R_recipes=560, inherited_R_sheets=1120,
        independent_two_high_families=20, positive_S_old8_sheets=160,
        continuously_dotted_rectangles=4480, combined_R_S_sheet_pairs=17920,
        minimum_S_norm_remainder=min(rho_lowers), minimum_positive_S_old8=min(positive_lowers),
        maximum_product_upper=max(upper_products), unresolved=0, accepted=0, new_cones=0)
    assert d['summary'] == enc(summary)
    assert 8-2*(3+1) == 0
    assert d['high_RM_nonright_SM_obstruction'] == dict(RM_square_at_least=3,
        SM_square_at_least=1, maximum_three_root_determinant=0, strict_positive_definiteness_impossible=True)
    assert d['exact_deletion_debit'] == old['exact_deletion_debit'] and d['deleted_walls'] == old['deleted_walls']
    assert d['exact_deletion_debit']['lost_finite_count'] == 19
    assert d['exact_deletion_debit']['compact_faces_by_dimension'] == {'2': 12, '3': 1}
    assert d['topology'] == old['topology']
    print('PASS_STANDALONE_NONRIGHT_HIGH_S_INTERVAL', 'parent_N_extensions30', 'families160',
          'continuous_rectangles4480', 'maximum_SR', float(max(upper_products)), 'new_cones0', flush=True)


if __name__ == '__main__':
    check()
