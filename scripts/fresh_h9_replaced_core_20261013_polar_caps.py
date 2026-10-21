#!/usr/bin/env python3
"""Targeted finite-polar-angle caps on changed, affine-free H9 cores."""
import argparse
from itertools import combinations, product
import json
import numpy as np
import mpmath as mp
from scripts.fresh_h9_replaced_core_20261013_base import *
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import finite_label

CORE = 'data/enumerations/fresh_h9_replaced_core_20261013_cores.json'
TARGET = 'data/enumerations/fresh_h9_replaced_core_20261013_base_geometry.json'
OUT = ROOT/'data/enumerations/fresh_h9_replaced_core_20261013_polar_caps.json'
POLAR_LABELS = (2, 3, 4, 5, 6, 10)


def specs_for(free, maximum=2):
    remaining = [i for i in range(10) if i not in free]
    yield (), ()
    for count in range(1, maximum+1):
        for positions in combinations(remaining, count):
            for labels in product((3, 5), repeat=count):
                yield positions, labels


def polar_square(q, label):
    factor = {2: Z, 3: scale(O, Q(1, 2)), 4: O,
              5: scale(add(O, P), Q(1, 2)), 6: scale(O, Q(3, 2)),
              10: scale(add(scale(O, 2), P), Q(1, 2))}[label]
    return mul(q, factor)


def exact_cap(spec, G, I, center, profiles, pool):
    u, v, positions, labels, polar_label, branch = spec
    polar = profiles['old348']
    assert polar['field'] is K
    R = polar_square(polar['q'], polar_label)
    E = Quadratic(K, R)
    kappa = (Z, O) if polar_label != 2 else E.zero
    w0, w1 = [E.zero]*10, [E.zero]*10
    for i, label in zip(positions, labels, strict=True):
        w0[i] = E.lift(O if label == 3 else P)
    pc = polar['c']
    w0[v] = E.div(E.sub(kappa, dot(E, tuple(E.lift(x) for x in pc), w0)),
                   E.lift(pc[v]))
    w1[u], w1[v] = E.one, E.lift(K.neg(div(pc[u], pc[v])))
    IE = tuple(tuple(E.lift(t) for t in row) for row in I)
    A = pairing(E, w1, IE, w1)
    B = E.mul(E.rational(2), pairing(E, w0, IE, w1))
    C = E.sub(pairing(E, w0, IE, w0), E.rational(2))
    D = E.sub(E.mul(B, B), E.mul(E.rational(4), E.mul(A, C)))
    assert E.sign(D) > 0 and E.sign(A)
    F = Quadratic(E, D)
    x = (E.div(E.neg(B), E.mul(E.rational(2), A)),
         E.div(E.rational(branch), E.mul(E.rational(2), A)))
    pa = tuple(F.neg(F.add(F.lift(a), F.mul(x, F.lift(b)))) for a, b in zip(w0, w1, strict=True))
    IF = tuple(tuple(F.lift(t) for t in row) for row in IE)
    c = matvec(F, IF, pa)
    assert F.sign(F.sub(dot(F, pa, c), F.rational(2))) == 0
    assert all(F.sign(p) <= 0 for p in pa)
    assert F.sign(dot(F, c, tuple(lift_K(F, t) for t in center))) > 0
    row = {'field': F, 'c': c, 'p': pa, 'q': scale(O, 2)}
    # Complete finite seed gate; parallel seed pairs are rejected immediately.
    seed_relations = {}
    for i in range(10):
        r = relation(row, profiles[f'seed{i}'], pool)
        assert r['kind'] in ('finite', 'dotted')
        seed_relations[str(i)] = r
    anchor = relation(row, polar, pool)
    assert anchor == {'kind': 'finite', 'label': polar_label}
    return row, {
        'free_seed_indices': [u, v], 'fixed_seed_indices': list(positions),
        'fixed_seed_labels': list(labels), 'polar_label': polar_label,
        'branch': branch, 'base_radicand': enc(R),
        'line_constant': list(map(E.encode, w0)), 'line_direction': list(map(E.encode, w1)),
        'norm_polynomial': list(map(E.encode, (A, B, C))),
        'discriminant': E.encode(D), 'parameter': F.encode(x),
        'seed_pairings': list(map(F.encode, pa)), 'seed_coefficients': list(map(F.encode, c)),
        'seed_relations': seed_relations}


def run(limit=80, max_fixed=2, minimum_retained=None, required=(), target_case=0):
    G, I, h, ps = profiles()
    core = read(CORE)
    add_required(ps, core, required)
    if minimum_retained is None:
        minimum_retained = 39+len(required)
    pool = angle_pool(ANGLE_SOURCE)
    target = read(TARGET)
    ys = [tuple(map(parse, r['certificate']['y'])) for r in target['cases'][target_case]['geometry']['escapes']]
    mp.mp.dps = 80
    In = np.array([[float(value(K, x)) for x in row] for row in I])
    hn = np.array([float(value(K, x)) for x in h])
    Yn = np.array([[float(value(K, x)) for x in y] for y in ys])
    names = list(dict.fromkeys(n for case in core['cases'] for n in case['walls']))
    cn = np.array([[float(value(ps[n]['field'], x)) for x in ps[n]['c']] for n in names])
    qn = np.array([float(value(K, ps[n]['q'])) for n in names])
    coreids = [[names.index(n) for n in case['walls']] for case in core['cases']]
    pc = np.array([float(value(K, t)) for t in ps['old348']['c']])
    allowed = np.array([np.cos(np.pi/m)**2 for m in (2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 24, 30)])
    proposals, counts = [], {'linear_families': 0, 'real_seed_center': 0, 'target_cutters': 0,
                              'coherent_proposals': 0}
    for uv in combinations(range(10), 2):
        u, v = uv
        if abs(pc[v]) < 1e-9:
            u, v = v, u
        if abs(pc[v]) < 1e-9:
            continue
        specs = list(specs_for((u, v), max_fixed))
        fixed = np.zeros((len(specs), 10))
        for j, (positions, labels) in enumerate(specs):
            for i, label in zip(positions, labels):
                fixed[j, i] = 1 if label == 3 else float(value(K, P))
        w1 = np.zeros(10)
        w1[u], w1[v] = 1, -pc[u]/pc[v]
        A = w1 @ In @ w1
        if abs(A) < 1e-10:
            continue
        for label in POLAR_LABELS:
            kappa = np.sqrt(float(value(K, polar_square(ps['old348']['q'], label))))
            w0 = fixed.copy()
            w0[:, v] = (kappa-fixed @ pc)/pc[v]
            B = 2*w0 @ In @ w1
            C = np.einsum('bi,ij,bj->b', w0, In, w0)-2
            D = B*B-4*A*C
            counts['linear_families'] += len(specs)
            for branch in (-1, 1):
                t = (-B+branch*np.sqrt(np.maximum(D, 0)))/(2*A)
                ws = w0+t[:, None]*w1
                ratios = ws*ws/4
                near = np.min(abs(ratios[:, :, None]-allowed), axis=2) < 1e-7
                ok = (D > 1e-9) & np.all((ws >= -1e-8) & ((ws > 2+1e-7) | near), axis=1)
                ok &= ws @ In @ hn < -1e-7
                counts['real_seed_center'] += int(ok.sum())
                cuts = ws @ In @ Yn.T
                ok &= np.any(cuts > 1e-7, axis=1)
                counts['target_cutters'] += int(ok.sum())
                for j in np.flatnonzero(ok):
                    vals = -cn @ ws[j]
                    ratio = vals*vals/(2*qn)
                    gate = (vals <= 1e-7) & ((ratio > 1+1e-7) |
                        (np.min(abs(ratio[:, None]-allowed), axis=1) < 1e-7))
                    retained = [sum(gate[ids]) for ids in coreids]
                    if max(retained) < minimum_retained:
                        continue
                    counts['coherent_proposals'] += 1
                    positions, labels = specs[j]
                    mask = tuple(i for i in range(len(ys)) if cuts[j, i] > 1e-7)
                    score = (max(retained), len(mask), min(cuts[j]), max(cuts[j]), -sum(ws[j]))
                    proposals.append((score, (u, v, positions, labels, label, branch), mask))
    proposals.sort(reverse=True)
    # Avoid repeated numeric mirrors before the exact source-anchored replay.
    seen, accepted, rejected = set(), [], []
    for score, spec, mask in proposals:
        if len(accepted)+len(rejected) >= limit:
            break
        try:
            row, record = exact_cap(spec, G, I, h, ps, pool)
        except AssertionError:
            rejected.append({'spec': spec, 'reason': 'exact_root_seed_gate'})
            continue
        numeric_key = tuple(round(float(value(row['field'], p)), 10) for p in row['p'])
        if numeric_key in seen:
            continue
        seen.add(numeric_key)
        checks = {name: relation(row, ps[name], pool) for name in names}
        retained = [[n for n in case['walls'] if checks[n]['kind'] in ('finite', 'dotted')]
                    for case in core['cases']]
        if max(map(len, retained)) < minimum_retained:
            rejected.append({'spec': spec, 'reason': 'exact_core_gate'})
            continue
        F = row['field']
        signs = [F.sign(dot(F, row['c'], tuple(lift_K(F, t) for t in y))) for y in ys]
        assert any(s < 0 for s in signs)
        record.update({'pair_checks': checks, 'retained_by_core': retained,
                       'target_signs': signs, 'cuts': [i for i, s in enumerate(signs) if s < 0],
                       'numeric_score': list(map(float, score))})
        accepted.append(record)
        print({'accepted': len(accepted), 'retained': list(map(len, retained)),
               'polar_label': record['polar_label'], 'cuts': record['cuts']}, flush=True)
    result = {'status': 'EXACT_TARGETED_POLAR_CAPS', 'core_source': CORE, 'target_source': TARGET,
              'target_case': target_case, 'required_caps': list(required),
              'angle_source': ANGLE_SOURCE, 'maximum_fixed_seed_support': max_fixed,
              'minimum_retained': minimum_retained, 'proposal_counts': counts,
              'numeric_proposals': len(proposals), 'caps': accepted, 'rejected': rejected,
              'scope': 'Bounded numerical proposals, exact accepted roots and pair relations. '
                       'Finite-label acceptance uses the saved26-label pool; no claim that '
                       'this pool exhausts arbitrary pair fields or the cap family.'}
    OUT.write_text(json.dumps(result, indent=2)+'\n')
    print({'status': result['status'], 'counts': counts, 'accepted': len(accepted)}, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=80)
    p.add_argument('--max-fixed', type=int, default=2)
    p.add_argument('--min-retained', type=int)
    p.add_argument('--required-file')
    p.add_argument('--target-source')
    p.add_argument('--target-case', type=int, default=0)
    p.add_argument('--output')
    a = p.parse_args()
    if a.target_source:
        TARGET = a.target_source
    if a.output:
        OUT = ROOT/a.output
    required = read(a.required_file) if a.required_file else []
    run(a.limit, a.max_fixed, a.min_retained, required, a.target_case)
