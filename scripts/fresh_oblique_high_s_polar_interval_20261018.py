#!/usr/bin/env python3
"""Physical interval obstruction for high S52 and high RM on twenty fixed N roots.

The original certificate supplies physical Gram enclosures. Five parent steps
and each selected N are reconstructed by exact interval Schur formulas. The
continuous S parameter is then eliminated by a quadratic inequality; no new
cone or native-wall bank is run.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import numpy as np
import fresh_compact10_cross_bridge_high_interval_20261018 as I
import fresh_oblique_high_n_right_r_join_20261018 as H

ROOT = I.ROOT
TAG = 'fresh_oblique_high_s_polar_interval_20261018'
PREVIOUS = 'fresh_oblique_high_rm_right_polar_20261018'
OLD = [8, 10, 52, 63, 64]
I.PREC = 10**70
V, add, sub, mul, div, neg, total, sqrt = I.iv, I.add, I.sub, I.mul, I.div, I.neg, I.total, I.sqrt


def determinant(g):
    n = len(g)
    terms = []
    for p in permutations(range(n)):
        value = V((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n)))
        for i in range(n):
            value = mul(value, g[i][p[i]])
        terms.append(value)
    return total(terms)


def inverse(g):
    # The small second principal pivot in the M frame makes interval Gaussian
    # elimination unnecessarily wide. Adjugates avoid that intermediate pivot.
    n = len(g)
    d = determinant(g)
    assert d[0]*d[1] > 0
    return [[div(mul(V((-1)**(i+j)), determinant(
        [[g[a][b] for b in range(n) if b != i] for a in range(n) if a != j])), d)
        for j in range(n)] for i in range(n)], d


def theta(m):
    if m == 2:
        return V(0)
    if m == 3:
        return V(1)
    if m == 4:
        return sqrt(V(2))
    if m == 5:
        return div(add(V(1), sqrt(V(5))), V(2))
    if m == 6:
        return sqrt(V(3))
    if m == 8:
        return sqrt(add(V(2), sqrt(V(2))))
    if m == 10:
        return sqrt(div(add(V(5), sqrt(V(5))), V(2)))
    if m == 12:
        return div(add(sqrt(V(2)), sqrt(V(6))), V(2))
    assert m == 7
    def atan(x):
        z = sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(100))
        return z, z+x**201/F(201)
    pi = sub(mul(V(16), atan(F(1, 5))), mul(V(4), atan(F(1, 239))))
    x = div(pi, V(m))
    xx, term, value = mul(x, x), V(1), V(1)
    for k in range(1, 60):
        term = div(mul(neg(xx), term), V((2*k-1)*(2*k)))
        value = add(value, term)
    error = xx[1]*max(abs(t) for t in term)/F(119*120)
    return mul(V(2), (value[0]-error, value[1]+error))


def initial(primary):
    ids = I.WHITE+OLD
    metric = [list(map(V, r)) for r in primary['computational_gram']]
    center = list(map(V, primary['common_center']))
    roots = {}
    for row in primary['facet_witnesses']:
        if row['wall'] not in ids:
            continue
        hp = V(primary['center_products'][str(row['wall'])])
        assert hp[1] < 0
        roots[row['wall']] = [div(mul(V(2), sub(x, y)), hp)
                              for x, y in zip(center, map(V, row['point']))]
    assert set(roots) == set(ids)
    native = {a: {b: I.dot(roots[a], metric, roots[b]) for b in ids} for a in ids}
    wi, pivots = I.inverse([[native[a][b] for b in I.WHITE] for a in I.WHITE])
    assert all(p[0] > 0 for p in pivots)
    wc = {a: [native[a][w] for w in I.WHITE] for a in OLD}
    g = {a: {b: sub(native[a][b], I.dot(wc[a], wi, wc[b])) for b in OLD} for a in OLD}
    profiles = [[V(0)]*6]
    for pos in (4, 5):
        for m in (3, 4, 5):
            v = [V(0)]*6
            v[pos] = neg(theta(m))
            profiles.append(v)
    return g, wc, wi, profiles, dict(walls=ids, physical_Gram=native,
        white_inverse=wi, white_pivots=pivots, old_white_covectors=wc, projected_Gram=g)


def extension(g, wc, wi, profiles, numerical, wall, frame, profile, labels, reference):
    gram = [[g[a][b] for b in frame] for a in frame]
    inv, det = inverse(gram)
    assert det[1] < 0 and gram[0][0][0] > 0
    rhs = [sub(neg(theta(m)), I.dot(profiles[profile], wi, wc[i]))
           for i, m in zip(frame, labels)]
    square = sub(V(2), I.dot(profiles[profile], wi, profiles[profile]))
    remainder = sub(square, I.dot(rhs, inv, rhs))
    assert remainder[0] > 0
    candidates = []
    for w in g:
        if w in frame:
            continue
        cv = [g[i][w] for i in frame]
        q = sub(g[w][w], I.dot(cv, inv, cv))
        if q[0] > 0:
            candidates.append((q[0], w, q, cv))
    _, witness, q, cv = max(candidates)
    linear = I.dot(rhs, inv, cv)
    radical = sqrt(mul(remainder, q))
    branch_values = [sub(linear, radical), add(linear, radical)]
    observed = float(reference@H.T.M@numerical[witness])
    distances = [abs(observed-float(sum(v)/2)) for v in branch_values]
    index = int(np.argmin(distances))
    assert distances[index] < 2e-6 and distances[1-index] > 1e-4
    sign = (-1, 1)[index]
    factor = mul(V(sign), sqrt(div(remainder, q)))
    products = {}
    for j in g:
        if j in frame:
            # These are the defining physical equations, not inferred zeros.
            products[j] = rhs[frame.index(j)]
        else:
            cj = [g[i][j] for i in frame]
            products[j] = add(I.dot(rhs, inv, cj), mul(factor,
                sub(g[witness][j], I.dot(cv, inv, cj))))
    error = max(abs(float(sum(v)/2)-float(reference@H.T.M@numerical[j]))
                for j, v in products.items())
    assert error < 2e-6
    proof = dict(wall=wall, frame=frame, profile=profile, labels=labels,
        frame_Gram=gram, inverse=inv, determinant=det, rhs=rhs,
        projected_square=square, remainder=remainder, witness=witness,
        witness_residual=q, branch_witness_products=branch_values, witness_sign=sign,
        numerical_witness_product=observed, branch_comparison=distances,
        all_product_comparison_error=error, new_products=products)
    for j, v in products.items():
        g[j][wall] = v
    g[wall] = {**products, wall: square}
    wc[wall], numerical[wall] = profiles[profile], reference
    return proof


def continuous(g, nlabels):
    frame = [52, 68, 71]
    inv, det = inverse([[g[a][b] for b in frame] for a in frame])
    assert det[1] < 0
    rhs = [V(-1), V(0), V(0)]
    alpha = neg(I.dot(rhs, inv, rhs))
    assert alpha[0] > 0
    witness = 8
    cv = [g[i][witness] for i in frame]
    q8 = sub(g[witness][witness], I.dot(cv, inv, cv))
    assert q8[0] > 0
    # A is the projected solution for rhs(-1,0,0). B is the normalized
    # residual of old8. Exactly A²=-alpha, B²=1, A.B=0.
    aa = {j: I.dot(rhs, inv, [g[i][j] for i in frame]) for j in g}
    bb = {j: div(sub(g[witness][j], I.dot(cv, inv, [g[i][j] for i in frame])), sqrt(q8)) for j in g}
    z = (div(theta(6), sqrt(add(V(2), mul(V(3), alpha))))[0],
         div(V(2), sqrt(add(V(2), mul(V(4), alpha))))[1])
    positive8 = add(mul(z, aa[8]), bb[8])
    assert positive8[0] > 0
    # Divide S by its positive sqrt(2+alpha*theta²). Only V=zA-B survives.
    xn = [69, 73]
    xn_gram = [[g[a][b] for b in xn] for a in xn]
    xi, xd = inverse(xn_gram)
    assert xn_gram[0][0][0] > 0 and xd[0] > 0
    mm, ar, br = [g[a][71] for a in xn], [aa[a] for a in xn], [neg(bb[a]) for a in xn]
    m = sub(g[71][71], I.dot(mm, xi, mm))
    q = [sub(V(1), I.dot(br, xi, br)), neg(mul(V(2), I.dot(ar, xi, br))),
         sub(neg(alpha), I.dot(ar, xi, ar))]
    d = [sub(neg(bb[71]), I.dot(br, xi, mm)), sub(aa[71], I.dot(ar, xi, mm))]
    factor = sub(V(3), mul(V(2), m))
    h = [sub(neg(mul(factor, q[0])), mul(V(2), mul(d[0], d[0]))),
         sub(neg(mul(factor, q[1])), mul(V(4), mul(d[0], d[1]))),
         sub(neg(mul(factor, q[2])), mul(V(2), mul(d[1], d[1])))]
    evaluate = lambda c: add(c[0], mul(z, add(c[1], mul(z, c[2]))))
    q_range, h_range = evaluate(q), evaluate(h)
    assert q_range[1] < 0 and h_range[0] > 0
    squared_upper = F(3)-h_range[0]/(-q_range[0])
    assert 0 < squared_upper < 3
    return dict(N_labels=nlabels, S_frame=frame, S_frame_inverse=inv,
        S_frame_determinant=det, alpha=alpha, residual8_square=q8,
        A_products=aa, B_products=bb, theta_domain=[theta(6)[0], F(2)],
        transformed_z_domain=z, positive_S_sheet_old8=positive8,
        XN_Gram=xn_gram, XN_inverse=xi, XN_determinant=xd, M_residual_square=m,
        negative_S_residual_square_coefficients=q, M_cross_coefficients=d,
        strict_high_RM_obstruction_coefficients=h, S_residual_square_range=q_range,
        high_RM_obstruction_range=h_range, forced_RM_squared_upper=squared_upper)


def run():
    old = json.loads((ROOT/f'results/{H.TAG}.json').read_text())
    prev = json.loads((ROOT/f'results/{PREVIOUS}.json').read_text())
    primary = json.loads((ROOT/f'results/{I.PRIMARY}.json').read_text())
    for doc in (old, prev, primary):
        for path, digest in doc['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    assert old['summary']['native_N'] == 20 and prev['summary']['right_S_polars'] == 192
    assert primary['status'] == 'PASS_EXACT_JOINT49_62_REPLACEMENT_63_WALLS_199_FINITE'
    g0, wc0, wi, profiles, physical = initial(primary)
    numerical_old = {int(k): np.array(v) for k, v in old['numerical_snapshot']['old_normals'].items()}
    parents, rows = [], []
    for source in old['records']:
        g, wc = {a: r.copy() for a, r in g0.items()}, wc0.copy()
        numerical = {a: numerical_old[a] for a in OLD}
        chain = json.loads(source['instruction_key'])
        assert len(chain) == 5 and [r[0] for r in chain] == [67, 68, 69, 70, 71]
        parent = dict(index=source['index'], origins=source['origins'], instruction_key=source['instruction_key'], steps=[])
        parents.append(parent)
        for wall, frame, profile, labels, branch in chain:
            rhs = [-H.T.B.theta(m)-float(sum(I.dot(profiles[profile], wi, wc[i]))/2)
                   for i, m in zip(frame, labels)]
            square = float(sum(sub(V(2), I.dot(profiles[profile], wi, profiles[profile])))/2)
            numerical_roots, info = H.T.N.norm_roots([numerical[i] for i in frame], rhs, square)
            assert info['status'] == 'two_real_norm_branches'
            proof = extension(g, wc, wi, profiles, numerical, wall, frame, profile, labels, dict(numerical_roots)[branch])
            proof['saved_numerical_branch'] = branch
            parent['steps'].append(proof)
        for i in (67, 68, 69, 71):
            assert max(abs(numerical[i]-np.array(source['prefix_normals'][str(i)]))) < 2e-6
        for cap in source['caps']:
            high = cap['S_bank']['symbolic_high_families']
            remaining = [r for r in high if r[1][2] == 2]
            assert remaining == [[0, [6, 2, 2], [0]]]
            g1, wc1, nums = {a: r.copy() for a, r in g.items()}, wc.copy(), numerical.copy()
            proof = extension(g1, wc1, wi, profiles, nums, 73, [10, 67, 69], 0, cap['labels'], np.array(cap['normal']))
            row = dict(index=len(rows), parent=source['index'], N_index=cap['index'], N_branch=cap['branch'],
                       N_extension=proof, continuous=continuous(g1, cap['labels']))
            rows.append(row)
            print('CONTINUOUS', row['parent'], cap['labels'], 'PASS', flush=True)
    assert len(parents) == 2 and len(rows) == 20
    debit = old['exact_deletion_debit']
    assert debit['lost_finite_count'] == 19 and debit['compact_faces_by_dimension'] == {'2': 12, '3': 1}
    paths = {Path(__file__), Path(I.__file__), Path(H.__file__)}
    paths.update(ROOT/p for doc in (old, prev, primary) for p in doc['sources'])
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in (H.TAG, PREVIOUS, I.TAG, I.PRIMARY)
                 for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md')))
    return I.enc(dict(status='PASS_FIXED_HIGH_N_HIGH_S_HIGH_RM_PHYSICAL_INTERVAL_OBSTRUCTION',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        interval_scale=I.PREC, primary_projection=physical, profiles=profiles, parents=parents, records=rows,
        exact_deletion_debit=debit, deleted_walls=old['deleted_walls'], topology=old['topology'],
        necessary_conditions=dict(N_profile=0, S_profile=0, R_profile=0, NX_label_at_least=6,
            S52_label_at_least=6, RM_label_at_least=6, SL=2, SM=2, RX=2, RN=2, RS=2),
        summary=dict(parent_classes=2, original_parent_origins=4, parent_extensions=10, N_geometries=20,
            continuous_S_families=20, positive_old8_S_sheets=20, polar_square_exclusions=20,
            minimum_old8_normalized_positive=min(r['continuous']['positive_S_sheet_old8'][0] for r in rows),
            minimum_polynomial_margin=min(r['continuous']['high_RM_obstruction_range'][0] for r in rows),
            maximum_forced_RM_square=max(r['continuous']['forced_RM_squared_upper'] for r in rows),
            unresolved=0, accepted=0, new_cones=0),
        scope='Physical rational interval exclusion, conditional on frozen original exact control and seven-profile finite classifier, for two canonical parent chains and twenty declared high-N recipes. Covers both S sheets and both R sheets with all finite S52>=6 and RM>=6 simultaneously. No rounded-coefficient surrogate, all-N claim, low-RM/high-S exclusion, arbitrary-incidence exclusion, new normal accepted as a chamber, or ambient-volume conclusion.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = ROOT/f'results/{TAG}.json'
    assert args.verify or not path.exists()
    result = run()
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['summary'], flush=True)
