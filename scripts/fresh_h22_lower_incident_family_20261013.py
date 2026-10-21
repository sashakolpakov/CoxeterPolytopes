#!/usr/bin/env python3
"""Exact all-real-shift incident-wall census at the lower H22 cusp delta4.

The complete transported affine rays are hash-anchored to the earlier exact
orbit certificate. Both height signs and both normal orientations are covered.
No floating arithmetic, optimizer, or producer arithmetic is used. Every native
integer contraction/product is bounded before evaluation. The conclusion is
about acute-compatible retained-facet additions, not arbitrary inequalities
which could make a retained wall redundant.
"""
import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction as F
from functools import reduce
from math import gcd, isqrt
from pathlib import Path

import numpy as np
import sympy as sp

P = Path(__file__).resolve().parents[1]
PRE = 'fresh_h22_lower_incident_family_20261013'
OUT = 'results/' + PRE + '.json'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
BASE = 'fresh_h22_new_cusp_completion_20261013'
COMPLETE = 'results/' + BASE + '_verify.json'
LOWER = 'results/' + BASE + '_lower.json'
LOCAL = 'results/fresh_h22_lower_new_cusp_20261013.json'
FRAMES = 'data/enumerations/' + BASE + '_frames.npz'
HIGHEST = 'data/enumerations/' + BASE + '_highest_roots.npz'
LIMIT = 2**63


def read(path):
    return json.loads((P / path).read_text())


def sha(path):
    return hashlib.sha256((P / path).read_bytes()).hexdigest()


def mx(a):
    assert a.dtype == np.int64
    assert not np.any(a == -LIMIT)
    return int(np.abs(a).max(initial=0))


def matmul(a, b):
    assert a.shape[-1] == b.shape[0]
    assert a.shape[-1] * mx(a) * mx(b) < LIMIT
    return a @ b


def mul(a, b):
    assert mx(a) * mx(b) < LIMIT
    return a * b


def subtract(a, b):
    assert mx(a) + mx(b) < LIMIT
    return a - b


def dot(a, b):
    return sum((x * y for x, y in zip(a, b, strict=True)), F(0))


def qmv(Q, x):
    return [dot(row, x) for row in Q]


def sqrt_rational(q):
    assert q >= 0
    a, b = isqrt(q.numerator), isqrt(q.denominator)
    return F(a, b) if a*a == q.numerator and b*b == q.denominator else None


def rational_pair_alphabet(q, oldnorm=F(2)):
    # For a rational cosine square, 4cos^2-2 is a rational algebraic integer.
    # The four finite possibilities and the parallel endpoint suffice for ALL m.
    out = {F(0)}
    for ratio in (F(1, 4), F(1, 2), F(3, 4), F(1)):
        root = sqrt_rational(q * oldnorm * ratio)
        if root is not None:
            out.add(-root)
    return sorted(out)


def run():
    complete, lower, local = read(COMPLETE), read(LOWER), read(LOCAL)
    assert all(sha(path) == digest for path, digest in complete['sources'].items())
    assert lower['source'] == COMPLETE and lower['source_sha256'] == sha(COMPLETE)
    assert all(sha(path) == digest for path, digest in local['sources'].items())
    source = read(RAW)['coxeter_input']
    Q = np.array([[-int(x) for x in row] for row in source['basis_gram']], dtype=np.int64)
    W = np.array(source['wall_coordinates_as_columns_transposed'], dtype=np.int64)
    h0 = [sum(int(row[j]) for row in W) for j in range(22)]
    div = reduce(gcd, h0)
    assert div > 0
    h = np.array([x // div for x in h0], dtype=np.int64)
    assert np.array_equal(matmul(Q, h), matmul(Q.T, h))
    hcov = matmul(Q, h)
    assert int(matmul(h, hcov)) == -14
    assert sp.Matrix(W.tolist()).rank() == 22
    assert max(matmul(W, hcov)) < 0
    oldterms = mul(W, matmul(W, Q))
    assert 22 * mx(oldterms) < LIMIT
    oldnorm = np.sum(oldterms, axis=1)
    assert set(map(int, oldnorm)) == {2, 4}
    frames, highest = np.load(P / FRAMES), np.load(P / HIGHEST)
    D, V, B = (frames[k] for k in ('six_cusp_old', 'first_tangent_numerator', 'second_tangent_numerator'))
    A, ids = highest['old_bases'], highest['cusp_indices']
    assert D.shape == V.shape == B.shape == (120960, 22)
    assert A.shape == (846720, 22)
    assert np.array_equal(ids, np.repeat(np.arange(len(D)), 7))
    for a in (D, V, B, A):
        mx(a)
    DD = mul(D, matmul(D, Q))
    assert 22 * mx(DD) < LIMIT and np.all(np.sum(DD, axis=1) == -1584)
    # c4 has old coordinate C/6 and new coordinate 2sqrt11.
    assert 6 * mx(W[32]) + mx(D[0]) < LIMIT
    C = D[0] - 6 * W[32]
    cov26, covC = matmul(Q, W[26]), matmul(Q, C)
    L, M0, Dh = matmul(D, cov26), matmul(D, covC), matmul(D, hcov)
    assert mx(M0) + 1584 < LIMIT
    assert all(not np.any((L == 0) & (M0 + 1584*sg == 0)) for sg in (1, -1))
    assert 6 * mx(W[26]) + mx(C) < LIMIT
    delta_six = C + 6 * W[26]
    ip = matmul(W, matmul(Q, delta_six))
    incident = np.flatnonzero(ip == 0).tolist()
    assert max(ip) == 0
    assert incident == [1, 3, 9, 10, 14, 21, 22, 23, 26, 28, 34, 39, 45, 46, 48, 50, 63]
    assert all(int(oldnorm[i]) == 2 for i in incident)
    archived = local['representatives'][0]
    assert archived['representative_type'] == 4 and archived['incident_old_wall_indices'] == incident
    assert archived['quotient_gram_rank'] == 17 and archived['lineality_dimension'] == 4
    Qf = [[F(int(x)) for x in row] + [F(0)] for row in Q] + [[F(0)]*22 + [F(11)]]
    c4 = [F(int(x), 6) for x in C] + [F(2)]
    delta = [F(int(x), 6) for x in delta_six] + [F(2)]
    old26 = [F(int(x)) for x in W[26]] + [F(0)]
    hf = [F(int(x)) for x in h] + [F(0)]
    assert list(map(str, delta)) == archived['null_ray']
    assert dot(delta, qmv(Qf, delta)) == 0 and dot(delta, qmv(Qf, hf)) == -72
    assert dot(c4, qmv(Qf, c4)) == 2 and dot(c4, qmv(Qf, old26)) == -2
    covs = [matmul(Q, W[i]) for i in incident]
    d_pairs = [matmul(D, cov) for cov in covs]
    stats, new_literals, duplicate_totals = [], [], Counter()
    angle_alphabets = {}
    families = [('highest_roots', A, ids, 1, None),
                ('first+', V, np.arange(len(D)), 28, 15624),
                ('first-', -V, np.arange(len(D)), 28, 15624),
                ('second+', B, np.arange(len(D)), 5208, 23870),
                ('second-', -B, np.arange(len(D)), 5208, 23870)]
    for name, base, ix, den, expected in families:
        bQ = matmul(base, Q)
        terms = mul(base, bQ)
        assert 22 * mx(terms) < LIMIT
        norms = np.sum(terms, axis=1)
        assert set(map(int, np.unique(norms))) == ({2, 4} if expected is None else {expected})
        terms = mul(bQ, D[ix])
        assert 22 * mx(terms) < LIMIT and np.all(np.sum(terms, axis=1) == 0)
        alphabets = {int(n): rational_pair_alphabet(F(int(n), den*den)) for n in np.unique(norms)}
        angle_alphabets[name] = {str(F(n, den*den)): list(map(str, pa)) for n, pa in alphabets.items()}
        # A positive-proportional old-wall copy is redundant, not an incompatible
        # new wall. Since a copy has zero added coordinate, it must have s=0.
        oldcopies = {}
        for n in alphabets:
            q = F(n, den*den)
            for i, w in enumerate(W):
                scale = sqrt_rational(q / int(oldnorm[i]))
                if scale is None:
                    continue
                coeff = den * scale
                row = [coeff * int(x) for x in w]
                if all(x.denominator == 1 for x in row):
                    key = tuple(int(x) for x in row)
                    assert key not in oldcopies or oldcopies[key] == i
                    oldcopies[key] = i
        U, VV, Bh = matmul(base, cov26), matmul(base, covC), matmul(base, hcov)
        for sg in (1, -1):
            ll, mm = L[ix], M0[ix] + 1584*sg
            j = np.flatnonzero(mul(U, mm) == mul(VV, ll))
            row = {'family': name, 'height_sign': sg, 'base_entries': len(base),
                   'incidence_solutions': len(j), 'degenerate_cusp_equation_pairs': 0}
            T = np.where(ll[j] != 0, ll[j], mm[j])
            a = np.where(ll[j] != 0, U[j], VV[j])
            assert np.all(T != 0)
            nh = subtract(mul(T, Bh[j]), mul(a, Dh[ix[j]]))
            nz = nh != 0
            row['zero_center_excluded'] = int(np.sum(~nz))
            # This selects the unique orientation with (r,h)<0, rather than
            # silently discarding the opposite normal orientation.
            o = -np.sign(nh[nz]) * np.sign(T[nz])
            j, T, a = j[nz], T[nz], a[nz]
            row['nonzero_center_solutions'] = len(j)
            duplicate_mask = np.zeros(len(j), dtype=bool)
            counts = Counter()
            z = np.flatnonzero(a == 0)
            if len(z) and oldcopies:
                signed = base[j[z]] * o[z, None]
                unique, inverse, multiplicity = np.unique(signed, axis=0, return_inverse=True, return_counts=True)
                flags = np.zeros(len(unique), dtype=bool)
                for k, normal in enumerate(unique):
                    old_id = oldcopies.get(tuple(map(int, normal)))
                    if old_id is not None:
                        assert old_id in incident
                        flags[k] = True
                        counts[old_id] += int(multiplicity[k])
                duplicate_mask[z] = flags[inverse]
            row['zero_shift_solutions'] = int(len(z))
            row['redundant_old_wall_copies'] = {str(k): v for k, v in sorted(counts.items())}
            row['redundant_literal_count'] = int(np.sum(duplicate_mask))
            duplicate_totals.update(counts)
            keep = ~duplicate_mask
            j, T, a, o = j[keep], T[keep], a[keep], o[keep]
            rejects = {}
            for wall, cov, dp in zip(incident, covs, d_pairs, strict=True):
                p = mul(o, subtract(mul(T, matmul(base[j], cov)), mul(a, dp[ix[j]])))
                good = np.zeros(len(j), dtype=bool)
                for n, alphabet in alphabets.items():
                    mask = norms[j] == n
                    for value in alphabet:
                        # p/(den*T) == value; use integer cross multiplication.
                        assert den * abs(value.numerator) * mx(T) < LIMIT
                        assert value.denominator * mx(p) < LIMIT
                        good |= mask & (p*value.denominator == T*(den*value.numerator))
                rejects[str(wall)] = int(np.sum(~good))
                j, T, a, o = j[good], T[good], a[good], o[good]
                if not len(j):
                    break
            row['first_old17_rejections'] = rejects
            row['remaining_new_literal_entries'] = j.tolist()
            stats.append(row)
            for entry, tt, aa, orient in zip(j, T, a, o, strict=True):
                entry, tt, aa, orient = map(int, (entry, tt, aa, orient))
                g = int(ix[entry])
                shift = F(-6*aa, den*tt)
                r = [orient * (F(int(x), den) + shift*F(int(y), 6))
                     for x, y in zip(base[entry], D[g], strict=True)] + [orient*2*sg*shift]
                rc = qmv(Qf, r)
                q = dot(r, rc)
                assert q == F(int(norms[entry]), den*den)
                assert dot(r, qmv(Qf, hf)) < 0
                assert dot(c4, rc) == dot(old26, rc) == dot(delta, rc) == 0
                assert all(dot([F(int(x)) for x in W[i]]+[F(0)], rc) in rational_pair_alphabet(q) for i in incident)
                bad = [(i, dot([F(int(x)) for x in w]+[F(0)], rc)) for i, w in enumerate(W)]
                bad = [(i, p) for i, p in bad if p > 0]
                assert bad
                new_literals.append({'family': name, 'entry': entry, 'cusp': g,
                                     'height_sign': sg, 'orientation': orient, 'shift': str(shift),
                                     'normal': list(map(str, r)), 'norm': str(q),
                                     'center_pairing': str(dot(hf, rc)),
                                     'positive_old_wall_witness': {'index': bad[0][0], 'pairing': str(bad[0][1])}})
    distinct = {}
    for row in new_literals:
        key = tuple(row['normal'])
        if key not in distinct:
            distinct[key] = {'normal': row['normal'], 'norm': row['norm'],
                             'center_pairing': row['center_pairing'],
                             'positive_old_wall_witness': row['positive_old_wall_witness'], 'sources': []}
        distinct[key]['sources'].append({k: row[k] for k in ('family', 'entry', 'cusp', 'height_sign', 'orientation', 'shift')})
    assert len(new_literals) == 5 and len(distinct) == 3
    distinct_rows = list(distinct.values())
    nonsquare_ratios = []
    for i, row in enumerate(distinct_rows):
        for previous in distinct_rows[:i]:
            ratio = F(row['norm']) / F(previous['norm'])
            # These vectors have rational coefficients in the fixed coordinate
            # basis, so proportionality would have a rational scalar.
            assert sqrt_rational(ratio) is None
            nonsquare_ratios.append(str(ratio))
    return {'status': 'PASS_ALL_REAL_SIGNED_AFFINE_RAY_FAMILY_CANNOT_REPAIR_FIXED_H22_CUSP',
            'sources': {p: sha(p) for p in (RAW, COMPLETE, LOWER, LOCAL, FRAMES, HIGHEST)},
            'retained_core': {'original_walls': 210, 'additional_wall': 'lower type4 c4',
                              'null_ray': list(map(str, delta)), 'incident_original_indices': incident,
                              'old_wall_span_rank': 22, 'positive_sum_center_divisor': div},
            'trivial_pre_reduction_copies': 'Positive-proportional copies of c4 or old26 are redundant and excluded before imposing orthogonality to those two normals; no nonpositive self-pairing is required for a repeated mirror.',
            'families': 'Seven highest-root bases and four signed residual tangent bases over all120960 cusp frames, both added-coordinate signs, all real shifts, both normal orientations.',
            'equations': {'C': 'D0-6W32', 'U': 'Bnum.Q.W26', 'V': 'Bnum.Q.C',
                          'L': 'Dg.Q.W26', 'M': 'Dg.Q.C+1584*height_sign',
                          'incidence': '6U+s*bden*L=0; 6V+s*bden*M=0',
                          'crossfilter': 'U*M=V*L', 'degenerate_L_M_zero': 0},
            'rational_allowed_incident_pairings': angle_alphabets,
            'enumeration': stats,
            'redundant_old_wall_copy_counts': {str(k): v for k, v in sorted(duplicate_totals.items())},
            'redundant_old_wall_literal_count': sum(duplicate_totals.values()),
            'new_literal_survivors_of_incident_gate': len(new_literals),
            'distinct_new_survivors_of_incident_gate': list(distinct.values()),
            'distinct_survivor_norm_ratios_are_rational_nonsquares': nonsquare_ratios,
            'acute_compatible_genuinely_new_incident_walls': 0,
            'scope': 'No finite addition from these complete signed all-real-shift affine-ray families can give finite volume while retaining c4 and all original210 facets and requiring nonpositive pairings with every distinct retained mirror. Every compatible incident candidate is a redundant positive-proportional copy of a retained incident wall. Nonincident additions leave the deficient cusp link unchanged sufficiently far out. Arbitrary inequalities making retained facets redundant, replacement/deletion of retained walls, and different normal families are not excluded.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.verify:
        assert result == read(OUT)
    else:
        (P / OUT).write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'])
    print('redundant literal copies:', result['redundant_old_wall_literal_count'],
          'new mirrors rejected:', len(result['distinct_new_survivors_of_incident_gate']))
