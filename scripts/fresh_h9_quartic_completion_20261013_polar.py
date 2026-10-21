#!/usr/bin/env python3
"""Two-distance polar-angle cap. --verify uses only Fraction arithmetic.

Numerical halfspace vertices propose a strict rational spacelike certificate;
all cap, pairing, center and escape conclusions are checked exactly.
"""
import argparse
from fractions import Fraction as Q
import json
from scripts.fresh_h9_quartic_completion_20261013_verify import (
    ROOT, read, matrix, load_old, node)
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import (
    K, Quadratic, Z, O, P, parse, enc, add, sub, mul, div, scale,
    dot, matvec, pairing)

SOURCE = 'data/enumerations/fresh_h9_quartic_completion_20261013_second_geometry.json'
OUTPUT = ROOT/'data/enumerations/fresh_h9_quartic_completion_20261013_polar.json'


def lift_K(f, x):
    if f is K:
        return x
    return f.lift(lift_K(f.base, x))


def attach(field, radicand_K, radical_E, coefficients):
    """Embed a saved K(sqrt(d0))(sqrt(d1)) cap over an arbitrary K-tower."""
    t = Quadratic(field, lift_K(field, radicand_K))
    embed_E = lambda x: (lift_K(field, x[0]), lift_K(field, x[1]))
    u = Quadratic(t, embed_E(radical_E))
    return u, tuple((embed_E(x[0]), embed_E(x[1])) for x in coefficients)


def build():
    geometry = read(SOURCE)
    case = geometry['cases'][0]
    assert case['facet_count'] == 44
    caps = read(geometry['source'])
    base = read(caps['base_source'])
    G, I = matrix(base['seed_gram']), matrix(base['seed_gram_inverse'])
    center = tuple(map(parse, base['strict_center_y']))
    assert K.sign(pairing(K, center, I, center)) < 0
    old = load_old(G)
    A, B, C = (add(scale(O, -20), scale(P, Q(25, 2))),
               add(scale(O, -25), scale(P, 15)),
               add(scale(O, Q(-42, 5)), scale(P, Q(24, 5))))
    D = sub(mul(B, B), scale(mul(A, C), 4))
    assert D == add(scale(O, -62), scale(P, 39))
    assert K.sign(A) > 0 and K.sign(B) < 0 and K.sign(C) < 0 and K.sign(D) > 0
    E = Quadratic(K, D)
    s = (div(K.neg(B), scale(A, 2)), div(O, scale(A, 2)))
    assert E.sign(s) > 0
    assert E.sign(E.sub(s, E.rational(3))) > 0
    assert E.sign(E.sub(s, E.rational(4))) < 0
    assert E.sign(E.add(E.add(E.mul(E.lift(A), E.mul(s, s)),
                             E.mul(E.lift(B), s)), E.lift(C))) == 0
    delta2 = add(scale(O, 2), P)
    F = Quadratic(E, E.lift(delta2))
    w = [F.zero]*10
    w[4] = E.zero, s
    w[8] = E.zero, E.add(s, E.one)
    pa = tuple(F.neg(x) for x in w)
    IF = tuple(tuple(lift_K(F, x) for x in row) for row in I)
    c = matvec(F, IF, pa)
    norm = dot(F, c, pa)
    assert F.sign(F.sub(norm, F.rational(2))) == 0
    cp = dot(F, c, tuple(lift_K(F, x) for x in center))
    assert F.sign(cp) > 0
    assert all(F.sign(x) <= 0 for x in pa)
    seedchecks = []
    for i, p in enumerate(pa):
        if i in (4, 8):
            assert F.sign(F.sub(F.mul(p, p), F.rational(4))) > 0
            seedchecks.append({'index': i, 'kind': 'dotted'})
        else:
            assert F.sign(p) == 0
            seedchecks.append({'index': i, 'kind': 'finite', 'label': 2})
    oldchecks = []
    for index in case['retained_old_caps']:
        n = old[index]
        T = Quadratic(F, lift_K(F, n['d']))
        value = (dot(F, pa, tuple(lift_K(F, x) for x in n['a'])),
                 dot(F, pa, tuple(lift_K(F, x) for x in n['b'])))
        assert T.sign(value) < 0
        square = T.mul(value, value)
        if index == 348:
            assert T.sign(T.sub(square, lift_K(T, div(n['q'], scale(O, 2))))) == 0
            kind = {'kind': 'finite', 'label': 3}
        else:
            assert T.sign(T.sub(square, lift_K(T, scale(n['q'], 2)))) > 0
            kind = {'kind': 'dotted'}
        oldchecks.append({'index': index, 'pairing': T.encode(value), **kind})
    Eold = Quadratic(K, parse(caps['anchor_radicand']))
    newchecks, prior = [], []
    for index in case['selected_new_caps']:
        ref = geometry['cap_references'][index]
        src = read(ref['source'])
        r = src['caps'][ref['index']]
        f = Quadratic(Eold, Eold.parse(r['discriminant_over_E']))
        cs = tuple(map(f.parse, r['seed_coefficients']))
        ps = tuple(map(f.parse, r['seed_pairings']))
        prior.append((ref, f, cs, ps))
        T, csT = attach(F, Eold.d, f.d, cs)
        value = dot(T, tuple(T.lift(T.base.lift(x)) for x in pa), csT)
        assert T.sign(value) < 0
        assert T.sign(T.sub(T.mul(value, value), T.rational(4))) > 0
        newchecks.append({'reference': ref, 'kind': 'dotted', 'pairing': T.encode(value)})
    target_data = read(caps['target_source'])
    targets = [tuple(map(parse, r['certificate']['y'])) for r in target_data['targets']]
    target_signs = [F.sign(dot(F, c, tuple(lift_K(F, t) for t in y))) for y in targets]
    assert target_signs[3] < 0
    record = {
        'status': 'EXACT_POLAR_PI3_CAP',
        'geometry_source': SOURCE, 'base_source': caps['base_source'],
        'A': enc(A), 'B': enc(B), 'C': enc(C), 'D': enc(D),
        'delta_squared': enc(delta2), 'positive_s': E.encode(s),
        'field_degree_upper_bound': 8,
        'seed_pairings': list(map(F.encode, pa)),
        'seed_coefficients': list(map(F.encode, c)), 'norm': F.encode(norm),
        'center_inequality': F.encode(cp),
        'seed_checks': seedchecks, 'old_checks': oldchecks, 'new_checks': newchecks,
        'target_source': caps['target_source'], 'target_signs': target_signs,
        'target_cuts': [i for i, s in enumerate(target_signs) if s < 0],
        'facet_count': 45,
        'scope': 'One exact cap preserving all 44 existing facets. Its only new finite '
                 'nonright angle is pi/3 with polar348; remaining nonzero pairings are dotted.'}
    return record, (G, I, center, old, case, F, pa, c, prior)


def mirror_check(context):
    G, I, center, old, case, F, pa, c, prior = context
    previous = old + [node(r, G) for r in read(
        'data/enumerations/fresh_h9_joint_completion_20261013.json')['caps']]
    sq = tuple(F.mul(p, p) for p in pa)
    support = [i for i, p in enumerate(pa) if F.sign(p)]
    counts = {'preceding_entries': len(previous)+165, 'matching_zero_pattern': 0}
    for n in previous:
        q = Quadratic(K, n['d'])
        ps = tuple(zip(n['pa'], n['pb']))
        if [i for i, p in enumerate(ps) if q.sign(p)] != support:
            continue
        counts['matching_zero_pattern'] += 1
        T = Quadratic(F, lift_K(F, n['d']))
        normalized = [T.mul(T.mul((lift_K(F, p[0]), lift_K(F, p[1])),
                                  (lift_K(F, p[0]), lift_K(F, p[1]))),
                            lift_K(T, div(scale(O, 2), n['q']))) for p in ps]
        assert any(T.sign(T.sub(T.lift(a), b)) for a, b in zip(sq, normalized, strict=True))
    for suffix in ('targeted', 'adaptive', 'second'):
        src = read(f'data/enumerations/fresh_h9_quartic_completion_20261013_{suffix}.json')
        e = Quadratic(K, parse(src['anchor_radicand']))
        for r in src['caps']:
            f = Quadratic(e, e.parse(r['discriminant_over_E']))
            ps = tuple(map(f.parse, r['seed_pairings']))
            if [i for i, p in enumerate(ps) if f.sign(p)] != support:
                continue
            counts['matching_zero_pattern'] += 1
            T, psT = attach(F, e.d, f.d, ps)
            assert any(T.sign(T.sub(T.lift(T.base.lift(a)), T.mul(b, b)))
                       for a, b in zip(sq, psT, strict=True))
    return {'status': 'DISTINCT_FROM_ALL_PRECEDING_ENTRIES', **counts}


def check_escape(cert, context):
    G, I, center, old, case, F, pa, c, prior = context
    y = tuple(map(parse, cert['y']))
    assert len(y) == 10 and all(K.sign(t) > 0 for t in y)
    norm, cp = pairing(K, y, I, y), pairing(K, y, I, center)
    assert K.sign(norm) > 0 and K.sign(cp) < 0
    assert norm == parse(cert['norm']) and cp == parse(cert['center_pairing'])
    signs = []
    for i in case['retained_old_caps']:
        n = old[i]
        f = Quadratic(K, n['d'])
        signs.append(f.sign((dot(K, n['a'], y), dot(K, n['b'], y))))
    for ref, f, cs, ps in prior:
        signs.append(f.sign(dot(f, cs, tuple(lift_K(f, t) for t in y))))
    signs.append(F.sign(dot(F, c, tuple(lift_K(F, t) for t in y))))
    assert all(s > 0 for s in signs)
    assert signs == cert['nonseed_inequality_signs']
    caps = read(read(SOURCE)['source'])
    targets = read(caps['target_source'])['targets']
    assert cert['identical_saved_target_indices'] == [
        i for i, target in enumerate(targets)
        if tuple(map(parse, target['certificate']['y'])) == y]
    projective = []
    for i, target in enumerate(targets):
        t = tuple(map(parse, target['certificate']['y']))
        if all(mul(y[0], b) == mul(t[0], a) for a, b in zip(y, t, strict=True)):
            projective.append(i)
    assert cert['projectively_equal_saved_target_indices'] == projective


def propose_escape(context):
    # The replay path does not import numerical or producer modules.
    import mpmath as mp
    import numpy as np
    from scipy.spatial import HalfspaceIntersection
    mp.mp.dps = 100
    def numeric(f, x):
        if f is K:
            return (mp.mpf(x[0].numerator)/x[0].denominator+
                    mp.mpf(x[1].numerator)/x[1].denominator*(1+mp.sqrt(5))/2)
        return numeric(f.base, x[0])+mp.sqrt(numeric(f.base, f.d))*numeric(f.base, x[1])
    G, I, center, old, case, F, pa, c, prior = context
    cn = []
    for i in case['retained_old_caps']:
        n = old[i]
        f = Quadratic(K, n['d'])
        cn.append([float(numeric(f, x)) for x in zip(n['a'], n['b'])])
    cn += [[float(numeric(f, x)) for x in cs] for ref, f, cs, ps in prior]
    cn.append([float(numeric(F, x)) for x in c])
    A = np.vstack((np.eye(10), np.array(cn)))
    hn = np.array([float(numeric(K, x)) for x in center])
    hn /= hn.sum()
    In = np.array([[float(numeric(K, x)) for x in row] for row in I])
    half = np.c_[A[:, 9, None]-A[:, :9], -A[:, 9]]
    half /= np.linalg.norm(half[:, :9], axis=1)[:, None]
    vertices = HalfspaceIntersection(half, hn[:9], qhull_options='Qx').intersections
    ys = np.c_[vertices, 1-vertices.sum(axis=1)]
    norms = np.einsum('bi,ij,bj->b', ys, In, ys)
    stats = {'numeric_vertices': len(ys), 'numeric_spacelike_vertices': int((norms > 1e-8).sum())}
    for index in sorted(np.flatnonzero(norms > 1e-8), key=lambda i: -norms[i]):
        for epsilon in (1e-3, 1e-5, 1e-7):
            yn = ys[index]+epsilon*hn
            for denominator in (10**9, 10**12, 10**15):
                y = tuple(scale(O, Q(round(float(t)*denominator), denominator)) for t in yn)
                cert = {'y': list(map(enc, y)), 'norm': enc(pairing(K, y, I, y)),
                        'center_pairing': enc(pairing(K, y, I, center)),
                        'nonseed_inequality_signs': [1]*(len(case['retained_old_caps'])+len(prior)+1)}
                caps = read(read(SOURCE)['source'])
                targets = read(caps['target_source'])['targets']
                cert['identical_saved_target_indices'] = [
                    i for i, target in enumerate(targets)
                    if tuple(map(parse, target['certificate']['y'])) == y]
                cert['projectively_equal_saved_target_indices'] = [
                    i for i, target in enumerate(targets)
                    if all(mul(y[0], b) == mul(parse(target['certificate']['y'][0]), a)
                           for a, b in zip(y, map(parse, target['certificate']['y']), strict=True))]
                try:
                    check_escape(cert, context)
                except AssertionError:
                    continue
                return {**stats, 'status': 'EXACT_STRICT_SPACELIKE_ESCAPE', 'certificate': cert}
    return {**stats, 'status': 'UNRESOLVED_VOLUME'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    record, context = build()
    record['mirror_check'] = mirror_check(context)
    if args.verify:
        saved = json.loads(OUTPUT.read_text())
        volume = saved.pop('volume')
        assert saved == record
        assert volume['status'] == 'EXACT_STRICT_SPACELIKE_ESCAPE'
        check_escape(volume['certificate'], context)
        record['volume'] = volume
    else:
        # Save the exact positive root before numerical volume work.
        OUTPUT.write_text(json.dumps(record, indent=2)+'\n')
        print({'status': record['status'], 'pairings': 44,
               'target_cuts': record['target_cuts'], 'mirror_check': record['mirror_check']}, flush=True)
        record['volume'] = propose_escape(context)
        OUTPUT.write_text(json.dumps(record, indent=2)+'\n')
    print({'status': 'PASS_EXACT_POLAR_COUPLING', 'facets': record['facet_count'],
           'volume_status': record['volume']['status']}, flush=True)


if __name__ == '__main__':
    main()
