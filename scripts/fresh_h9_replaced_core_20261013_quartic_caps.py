#!/usr/bin/env python3
"""Three-distance caps using the actual newly retained quartic anchor."""
import argparse
from itertools import combinations, product
import json
import numpy as np
import mpmath as mp
from scripts.fresh_h9_replaced_core_20261013_base import *

CORE = 'data/enumerations/fresh_h9_replaced_core_20261013_cores.json'
TARGET = 'data/enumerations/fresh_h9_replaced_core_20261013_round2_selection.json'
REQUIRED = 'data/enumerations/fresh_h9_replaced_core_20261013_required2.json'
OUT = ROOT/'data/enumerations/fresh_h9_replaced_core_20261013_quartic_caps.json'


def exact(spec, G, I, center, ps, pool, anchor_ref):
    u, v, w, fixed_indices, fixed_labels, magnitude_label, polar_label, branch = spec
    E, a = collapsed_anchor(anchor_ref)
    b = tuple(E.lift(t) for t in ps['old348']['c'])
    M = E.lift({2: Z, 3: O, 5: P}[magnitude_label])
    N = E.lift(Z if polar_label == 2 else O)
    fixed = [E.zero]*10
    for j, label in zip(fixed_indices, fixed_labels):
        fixed[j] = E.lift(O if label == 3 else P)
    determinant = E.sub(E.mul(a[v], b[w]), E.mul(a[w], b[v]))
    assert E.sign(determinant)
    r0, r1 = E.sub(M, dot(E, a, fixed)), E.sub(N, dot(E, b, fixed))
    line0, line1 = list(fixed), [E.zero]*10
    line0[v] = E.div(E.sub(E.mul(r0, b[w]), E.mul(a[w], r1)), determinant)
    line0[w] = E.div(E.sub(E.mul(a[v], r1), E.mul(r0, b[v])), determinant)
    line1[u] = E.one
    line1[v] = E.div(E.sub(E.mul(a[w], b[u]), E.mul(a[u], b[w])), determinant)
    line1[w] = E.div(E.sub(E.mul(a[u], b[v]), E.mul(a[v], b[u])), determinant)
    IE = tuple(tuple(E.lift(t) for t in row) for row in I)
    A = pairing(E, line1, IE, line1)
    B = E.mul(E.rational(2), pairing(E, line0, IE, line1))
    C = E.sub(pairing(E, line0, IE, line0), E.rational(2))
    D = E.sub(E.mul(B, B), E.mul(E.rational(4), E.mul(A, C)))
    assert E.sign(D) > 0 and E.sign(A)
    F = Quadratic(E, D)
    x = (E.div(E.neg(B), E.mul(E.rational(2), A)),
         E.div(E.rational(branch), E.mul(E.rational(2), A)))
    pa = tuple(F.neg(F.add(F.lift(t), F.mul(x, F.lift(s))))
               for t, s in zip(line0, line1, strict=True))
    cf = matvec(F, tuple(tuple(F.lift(t) for t in row) for row in IE), pa)
    assert F.sign(F.sub(dot(F, cf, pa), F.rational(2))) == 0
    assert all(F.sign(t) <= 0 for t in pa)
    assert F.sign(dot(F, cf, tuple(lift_K(F, t) for t in center))) > 0
    assert F.sign(F.add(dot(F, pa, tuple(F.lift(t) for t in a)), F.lift(M))) == 0
    assert F.sign(F.add(dot(F, pa, tuple(F.lift(t) for t in b)), F.lift(N))) == 0
    row = {'field': F, 'c': cf, 'p': pa, 'q': scale(O, 2)}
    seedchecks = {str(i): relation(row, ps[f'seed{i}'], pool) for i in range(10)}
    assert all(r['kind'] in ('finite', 'dotted') for r in seedchecks.values())
    return row, {'kind': 'full_quartic_anchor', 'anchor_reference': anchor_ref,
        'anchor_label': magnitude_label, 'free_seed_indices': [u, v, w],
        'fixed_seed_indices': list(fixed_indices), 'fixed_seed_labels': list(fixed_labels),
        'polar_label': polar_label, 'branch': branch, 'base_radicand': enc(E.d),
        'line_constant': list(map(E.encode, line0)), 'line_direction': list(map(E.encode, line1)),
        'norm_polynomial': list(map(E.encode, (A, B, C))), 'discriminant': E.encode(D),
        'parameter': F.encode(x), 'seed_pairings': list(map(F.encode, pa)),
        'seed_coefficients': list(map(F.encode, cf)), 'seed_relations': seedchecks}


def run(limit=30, maximum_fixed=2, focus=0):
    G, I, h, ps = profiles()
    core, required = read(CORE), read(REQUIRED)
    add_required(ps, core, required)
    anchor_ref = required[0]
    E, a = collapsed_anchor(anchor_ref)
    mp.mp.dps = 90
    an = np.array([float(value(E, t)) for t in a])
    bn = np.array([float(value(K, t)) for t in ps['old348']['c']])
    In = np.array([[float(value(K, t)) for t in row] for row in I])
    hn = np.array([float(value(K, t)) for t in h])
    ys = [tuple(map(parse, t['certificate']['y'])) for t in read(TARGET)['cases'][0]['geometry']['escapes']]
    Yn = np.array([[float(value(K, t)) for t in y] for y in ys])
    names = list(dict.fromkeys(n for case in core['cases'] for n in case['walls']))
    cn = np.array([[float(value(ps[n]['field'], t)) for t in ps[n]['c']] for n in names])
    qn = np.array([float(value(K, ps[n]['q'])) for n in names])
    ci = [[names.index(n) for n in case['walls']] for case in core['cases']]
    allowed = np.array([np.cos(np.pi/m)**2 for m in (2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 24, 30)])
    proposals = []
    counts = {'linear_families': 0, 'seed_center_roots': 0, 'focus_cutters': 0, 'coherent': 0}
    for free in combinations(range(10), 3):
        u, v, w = free
        # Choose a nonzero 2x2 minor for the two anchor equations.
        best = max(((abs(an[j]*bn[k]-an[k]*bn[j]), i, j, k)
                    for i in free for j, k in combinations([z for z in free if z != i], 2)))
        det, u, v, w = best
        if det < 1e-9:
            continue
        det = an[v]*bn[w]-an[w]*bn[v]
        remaining = [i for i in range(10) if i not in free]
        specs = [((), ())]
        for count in range(1, maximum_fixed+1):
            specs += [(s, labels) for s in combinations(remaining, count)
                      for labels in product((3, 5), repeat=count)]
        fixed = np.zeros((len(specs), 10))
        for j, (positions, labels) in enumerate(specs):
            for i, label in zip(positions, labels):
                fixed[j, i] = 1 if label == 3 else float(value(K, P))
        line1 = np.zeros(10)
        line1[u] = 1
        line1[v] = (an[w]*bn[u]-an[u]*bn[w])/det
        line1[w] = (an[u]*bn[v]-an[v]*bn[u])/det
        A = line1 @ In @ line1
        if abs(A) < 1e-9:
            continue
        for am, pm in product((2, 3, 5), (2, 10)):
            M, N = {2: 0, 3: 1, 5: float(value(K, P))}[am], 0 if pm == 2 else 1
            r0, r1 = M-fixed@an, N-fixed@bn
            line0 = fixed.copy()
            line0[:, v] = (r0*bn[w]-an[w]*r1)/det
            line0[:, w] = (an[v]*r1-r0*bn[v])/det
            B = 2*line0 @ In @ line1
            C = np.einsum('bi,ij,bj->b', line0, In, line0)-2
            D = B*B-4*A*C
            counts['linear_families'] += len(specs)
            for branch in (-1, 1):
                t = (-B+branch*np.sqrt(np.maximum(D, 0)))/(2*A)
                ws = line0+t[:, None]*line1
                ratio = ws*ws/4
                near = np.min(abs(ratio[:, :, None]-allowed), axis=2) < 1e-7
                ok = (D > 1e-8) & np.all((ws >= -1e-7) & ((ws > 2+1e-7) | near), axis=1)
                ok &= ws @ In @ hn < -1e-7
                counts['seed_center_roots'] += int(ok.sum())
                cuts = ws @ In @ Yn.T
                ok &= cuts[:, focus] > 1e-7
                counts['focus_cutters'] += int(ok.sum())
                for j in np.flatnonzero(ok):
                    pairs = -cn @ ws[j]
                    rr = pairs*pairs/(2*qn)
                    good = (pairs <= 1e-7) & ((rr > 1+1e-7) |
                        (np.min(abs(rr[:, None]-allowed), axis=1) < 1e-7))
                    retained = [sum(good[ids]) for ids in ci]
                    if max(retained) < len(core['cases'][0]['walls']):
                        continue
                    counts['coherent'] += 1
                    positions, labels = specs[j]
                    score = (int((cuts[j] > 1e-7).sum()), cuts[j, focus], max(cuts[j]))
                    proposals.append((score, (u, v, w, positions, labels, am, pm, branch)))
    proposals.sort(reverse=True)
    print({'proposal_counts': counts, 'scored': len(proposals)}, flush=True)
    pool = angle_pool(ANGLE_SOURCE)
    caps, rejected, seen = [], [], set()
    for score, spec in proposals:
        if len(caps)+len(rejected) >= limit:
            break
        try:
            row, record = exact(spec, G, I, h, ps, pool, anchor_ref)
        except AssertionError:
            rejected.append({'spec': spec, 'reason': 'exact_norm_seed_gate'})
            continue
        nk = tuple(round(float(value(row['field'], p)), 10) for p in row['p'])
        if nk in seen:
            continue
        seen.add(nk)
        checks = {n: relation(row, ps[n], pool) for n in names}
        retained = [[n for n in case['walls'] if checks[n]['kind'] in ('finite', 'dotted')]
                    for case in core['cases']]
        if max(map(len, retained)) < len(core['cases'][0]['walls']):
            rejected.append({'spec': spec, 'reason': 'exact_core_gate'})
            continue
        F = row['field']
        signs = [F.sign(dot(F, row['c'], tuple(lift_K(F, t) for t in y))) for y in ys]
        assert signs[focus] < 0
        record.update({'pair_checks': checks, 'retained_by_core': retained,
                       'target_signs': signs, 'cuts': [i for i, s in enumerate(signs) if s < 0]})
        caps.append(record)
        print({'accepted': len(caps), 'retained': list(map(len, retained)),
               'anchor_label': record['anchor_label'], 'polar_label': record['polar_label'],
               'cuts': record['cuts']}, flush=True)
    result = {'status': 'EXACT_CHANGED_CORE_QUARTIC_CAPS', 'core_source': CORE,
              'target_source': TARGET, 'target_case': 0, 'focus_target': focus,
              'required_caps': required, 'angle_source': ANGLE_SOURCE,
              'proposal_counts': counts, 'caps': caps, 'rejected': rejected,
              'scope': 'A bounded full-quartic-anchor construction targeted at an actual '
                       'leading escape of the changed41-wall core. No exhaustive family claim.'}
    OUT.write_text(json.dumps(result, indent=2)+'\n')
    print({'accepted': len(caps), 'rejected': len(rejected)}, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=30)
    p.add_argument('--max-fixed', type=int, default=2)
    p.add_argument('--focus', type=int, default=0)
    a = p.parse_args()
    run(a.limit, a.max_fixed, a.focus)
