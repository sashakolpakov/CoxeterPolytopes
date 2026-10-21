#!/usr/bin/env python3
"""Review the complete first-cap branches without importing their producer.

The old octic field and rigorous sign evaluator are shared. Linear algebra,
radical comparisons, and branch completeness are checked here directly.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import fresh_h9_joint_corner_20261015_bridge_review as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_two_caps_20261016_review'
PRIMARY = 'results/fresh_nonarithmetic_two_caps_20261016.json'
BASE = 'results/fresh_h9_joint_pair_20261015_two_extension.json'
BRIDGE = 'results/fresh_h9_joint_corner_20261015_bridge.json'
TARGET = 'results/fresh_nonarithmetic_escape_target_20261016.json'
K, Z, O = R.K, R.Z, R.O


def digest(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def read(p):
    doc = json.loads((ROOT / p).read_text())
    for name, value in doc['sources'].items():
        assert digest(name) == value, name
    return doc


def dec(v):
    return K.new([K.dom.convert(F(a)) for a in v])


def dot(a, b):
    assert len(a) == len(b)
    return sum((x*y for x, y in zip(a, b)), Z)


def mv(a, v):
    return [dot(row, v) for row in a]


def rank(a):
    a = [row[:] for row in a]
    rr = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rr, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rr], a[pivot] = a[pivot], a[rr]
        scale = a[rr][col]
        a[rr] = [x/scale for x in a[rr]]
        for i in range(rr+1, len(a)):
            scale = a[i][col]
            a[i] = [x-scale*y for x, y in zip(a[i], a[rr])]
        rr += 1
        if rr == len(a):
            break
    return rr


def sign_radical_sum(a, s, b, t):
    """Independent sign of a sqrt(s)+b sqrt(t), s,t nonnegative."""
    assert R.sg(s) >= 0 and R.sg(t) >= 0
    x, y = (R.sg(a) if s else 0), (R.sg(b) if t else 0)
    if x == 0:
        return y
    if y == 0 or x == y:
        return x
    squared_difference = R.sg(a*a*s - b*b*t)
    if squared_difference == 0:
        return 0
    return x if squared_difference > 0 else y


def positive_pivots(a):
    a = [row[:] for row in a]
    result = []
    for i in range(len(a)):
        q = a[i][i]
        assert R.sg(q) > 0
        result.append(q)
        for j in range(i+1, len(a)):
            for k in range(i+1, len(a)):
                a[j][k] -= a[j][i]*a[i][k]/q
    return result


def run():
    primary, base, bridge, target = map(read, (PRIMARY, BASE, BRIDGE, TARGET))
    assert primary['status'] == 'PASS_COMPLETE_BOUNDED_TWO_CAP_FIRST_STAGE_OBSTRUCTION'
    Q = R.m4(base['form'])
    normals = R.m4(base['normals']) + [R.v8(bridge['cap']['normal'])]
    H, E = R.v4(base['center']), R.v4(base['strict_escape']['point'])
    pair = lambda a, b: dot(a, mv(Q, b))
    cov = [mv(Q, n) for n in normals]
    norms = [pair(n, n) for n in normals]
    T = [0, 1, 2, 5, 6, 10, 30]
    assert primary['target_support'] == T
    assert len(normals) == 33 and all(R.sg(q) > 0 for q in norms)
    assert all(R.sg(pair(H, n)) < 0 for n in normals)
    assert R.sg(pair(H, H)) < 0 and R.sg(pair(E, E)) > 0
    endpoint = [[pair(normals[i], normals[j]) for j in T+[31]] for i in T+[31]]
    assert R.enc(positive_pivots(endpoint)) == primary['endpoint_positive_pivots']
    frames = {}
    for frame in primary['linear_frames']:
        anchor = frame['anchor']
        assert anchor in (3, 7) and anchor not in frames
        ids = T+[31, anchor]
        rows = [cov[i] for i in ids]
        p, k = list(map(dec, frame['particular'])), list(map(dec, frame['kernel']))
        assert len(p) == len(k) == 10 and rank(rows) == 9
        assert mv(rows, p) == [Z]*8+[-O] and mv(rows, k) == [Z]*9
        pp, kk = pair(p, p), pair(k, k)
        assert R.sg(pp) < 0 < R.sg(kk) and pair(p, k) == Z
        # The nine equations leave exactly a line; these identities identify
        # its full affine solution for each prescribed anchor magnitude.
        assert rank(rows+[mv(Q, k)]) == 10
        a, b = [pair(p, n) for n in normals], [pair(k, n) for n in normals]
        assert R.enc(a) == frame['old_pairing_a'] and R.enc(b) == frame['old_pairing_b']
        assert R.enc(pp) == frame['particular_norm'] and R.enc(kk) == frame['kernel_norm']
        frames[anchor] = p, k, pp, kk, a, b
    assert set(frames) == {3, 7}
    phi = (O+R.BASIS[1])/2
    prescribed = {2: Z, 3: O, 4: 2*O, 5: phi*phi}
    finite_cycles = {**prescribed, 6: 3*O}
    seen, counts, rows = set(), Counter(), []
    for recipe in primary['recipes']:
        anchor, label = recipe['anchor'], recipe['label']
        key = anchor, label
        assert key not in seen and anchor in frames and label in prescribed
        seen.add(key)
        p, k, pp, kk, a, b = frames[anchor]
        s = prescribed[label]
        d = (2*O-s*pp)/kk
        assert R.sg(d) > 0 and s*pp+d*kk == 2*O
        assert dec(recipe['s']) == s and dec(recipe['d']) == d
        assert sorted(br['sign'] for br in recipe['branches']) == [-1, 1]
        for branch in recipe['branches']:
            sign = branch['sign']
            signs = [sign_radical_sum(x, s, sign*y, d) for x, y in zip(a, b)]
            assert signs == branch['old_wall_pairing_signs']
            assert all(signs[i] == 0 for i in T+[31])
            assert signs[anchor] == (0 if label == 2 else -1)
            for name, vector in (('center', H), ('escape', E)):
                assert sign_radical_sum(pair(p, vector), s, sign*pair(k, vector), d) == branch[name+'_sign']
            positive = [i for i, value in enumerate(signs) if value > 0]
            assert positive == branch['positive_old_pairings']
            blocker = branch['blocker_wall']
            if positive:
                assert blocker in positive and branch['status'] == 'REJECT_OBTUSE'
                counts['obtuse'] += 1
            else:
                assert branch['status'] == 'REJECT_ALL_LABEL_ANGLE_GAP' and signs[blocker] == -1
                aa = 2*(s*a[blocker]**2+d*b[blocker]**2)/norms[blocker]
                bb = 4*sign*a[blocker]*b[blocker]/norms[blocker]
                cycle = branch['cartan_cycle']
                assert [dec(cycle[n]) for n in ('a', 'b', 'radicand')] == [aa, bb, s*d]
                lo_label, hi_label = branch['consecutive_labels']
                assert hi_label == lo_label+1
                if (lo_label, hi_label) == (7, 8):
                    lo, hi = K.convert(F(33, 10)), K.convert(F(10, 3))
                else:
                    lo, hi = finite_cycles[lo_label], finite_cycles[hi_label]
                assert dec(branch['lower_bound']) == lo and dec(branch['upper_bound']) == hi
                assert sign_radical_sum(aa-lo, O, bb, s*d) == 1
                assert sign_radical_sum(aa-hi, O, bb, s*d) == -1
                assert sign_radical_sum(aa, O, bb, s*d) == 1
                assert sign_radical_sum(aa-4*O, O, bb, s*d) == -1
                counts['angle_gap'] += 1
            rows.append(dict(anchor=anchor, label=label, branch=sign,
                             blocker=blocker, status=branch['status']))
            counts['wall_signs'] += 33
    assert seen == {(j, m) for j in (3, 7) for m in (2, 3, 4, 5)}
    assert counts == dict(obtuse=9, angle_gap=7, wall_signs=528)
    # Independently evaluate the exact rational bracket for c7 and c8.
    poly = lambda x: x**3-5*x*x+6*x-1
    assert poly(F(3)) == -1 and poly(F(33, 10)) == F(287, 1000)
    assert 3*3**2-10*3+6 == 3 and 6*3-10 == 8
    assert F(2) > (F(10, 3)-2)**2
    experiment = target['proposed_bounded_experiment']
    assert experiment['algebraic_recipe_count'] == 128
    paths = [f'scripts/{TAG}.py', PRIMARY, BASE, BRIDGE, TARGET,
             'results/fresh_nonarithmetic_two_caps_20261016.md',
             'scripts/fresh_h9_joint_corner_20261015_bridge_review.py',
             'scripts/fresh_h9_h8_all_high_20261014_review.py',
             'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_BOUNDED_TWO_CAP_BRANCH_REVIEW',
                sources={p: digest(p) for p in paths}, counts=dict(counts),
                first_stage_recipes=8, complete_isolated_branches=16,
                excluded_ordered_two_cap_recipes=128, branch_rejections=rows,
                scope='Independent linear rank, formal radical branches and all528 old-wall signs and seven universal finite-angle gaps. Shares the frozen octic arithmetic and rigorous field-sign evaluator. No producer import; no claim beyond the literal bounded family; no finite-volume construction.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    output = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(output.read_text())
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['counts'])
