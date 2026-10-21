#!/usr/bin/env python3
"""Independent exact reconstruction of the proposed two-distance polar cap.

No search/producer imports. New-old tower fields are embedded in explicit
composita, so the new quadratic radicand need not equal the old anchor's.
"""
import copy
import json
from fractions import Fraction as Q
from pathlib import Path

from scripts.fresh_h9_quartic_completion_20261013_verify import read, matrix, load_old, node
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import (
    K, Quadratic, Z, O, parse, dot, matvec, pairing)

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = 'data/enumerations/fresh_h9_quartic_completion_20261013_second_geometry.json'
SOURCE = 'data/enumerations/fresh_h9_quartic_completion_20261013_polar.json'
OUT = ROOT/'results/fresh_h9_quartic_completion_20261013_polar_independent.json'


def relation(field, value, other_norm):
    """Only claim the explicit zero/pi3/dotted possibilities of this cap."""
    s = field.sign(value)
    assert s <= 0
    if s == 0:
        return 'right'
    square = field.mul(value, value)
    if field.sign(field.sub(square, field.mul(field.rational(Q(1, 2)), other_norm))) == 0:
        return 'pi/3'
    assert field.sign(field.sub(square, field.mul(field.rational(2), other_norm))) > 0
    return 'dotted'


def reconstruct():
    geometry = read(GEOMETRY)
    case = geometry['cases'][0]
    capdata = read(geometry['source'])
    base = read(capdata['base_source'])
    G, I = matrix(base['seed_gram']), matrix(base['seed_gram_inverse'])
    center = tuple(map(parse, base['strict_center_y']))
    assert len(G) == len(I) == len(center) == 10
    for i, row in enumerate(G):
        assert matvec(K, I, row) == tuple(O if j == i else Z for j in range(10))
    assert all(K.sign(x) > 0 for x in center)
    assert K.sign(pairing(K, center, I, center)) < 0

    A = (Q(-20), Q(25, 2))
    B = (Q(-25), Q(15))
    C = (Q(-42, 5), Q(24, 5))
    D = (Q(-62), Q(39))
    assert K.sub(K.mul(B, B), K.mul(K.rational(4), K.mul(A, C))) == D
    assert tuple(K.sign(x) for x in (A, B, C, D)) == (1, -1, -1, 1)
    E = Quadratic(K, D)
    F = Quadratic(E, E.lift((Q(2), Q(1))))
    lift = lambda x: F.lift(E.lift(x))
    s = E.div((K.neg(B), O), E.lift(K.mul(K.rational(2), A)))
    assert E.sign(E.sub(s, E.rational(3))) > 0
    assert E.sign(E.sub(s, E.rational(4))) < 0
    assert E.sign(E.add(E.add(E.mul(E.lift(A), E.mul(s, s)), E.mul(E.lift(B), s)), E.lift(C))) == 0
    # A>0,C<0 prove opposite signs of the two roots, hence branch uniqueness.
    pa = [F.zero]*10
    pa[4] = F.neg((E.zero, s))
    pa[8] = F.neg((E.zero, E.add(s, E.one)))
    pa = tuple(pa)
    cf = matvec(F, tuple(tuple(lift(x) for x in row) for row in I), pa)
    assert matvec(F, tuple(tuple(lift(x) for x in row) for row in G), cf) == pa
    assert F.sign(F.sub(dot(F, cf, pa), F.rational(2))) == 0
    assert F.sign(dot(F, cf, tuple(lift(x) for x in center))) > 0
    rows = [{'kind': 'seed', 'index': j, 'relation': relation(F, x, lift(G[j][j]))}
            for j, x in enumerate(pa)]

    old = load_old(G)
    old_fields = []
    for j in case['retained_old_caps']:
        n = old[j]
        X = Quadratic(F, lift(n['d']))
        cc = tuple((lift(a), lift(b)) for a, b in zip(n['a'], n['b'], strict=True))
        pp = tuple(X.lift(x) for x in pa)
        value = dot(X, pp, cc)
        rows.append({'kind': 'old', 'index': j, 'relation': relation(X, value, X.lift(lift(n['q'])))})
        old_fields.append((X, cc))

    new_fields = []
    for j in case['selected_new_caps']:
        ref = geometry['cap_references'][j]
        source = read(ref['source'])
        oldE = Quadratic(K, parse(source['anchor_radicand']))
        c = source['caps'][ref['index']]
        oldF = Quadratic(oldE, oldE.parse(c['discriminant_over_E']))
        # F <- X=F[sqrt(old anchor)] <- Y=X[sqrt(old cap discriminant)].
        X = Quadratic(F, lift(oldE.d))
        embedE = lambda x: (lift(x[0]), lift(x[1]))
        Y = Quadratic(X, embedE(oldF.d))
        cc = tuple((embedE(x[0]), embedE(x[1])) for x in map(oldF.parse, c['seed_coefficients']))
        value = dot(Y, tuple(Y.lift(X.lift(x)) for x in pa), cc)
        rows.append({'kind': 'new', 'index': j, 'relation': relation(Y, value, Y.rational(2))})
        new_fields.append((Y, cc))
    assert len(rows) == case['facet_count'] == 44
    assert [r for r in rows if r['relation'] == 'pi/3'] == [
        {'kind': 'old', 'index': 348, 'relation': 'pi/3'}]

    targets = read(capdata['target_source'])
    target = tuple(map(parse, targets['targets'][3]['certificate']['y']))
    assert all(K.sign(x) > 0 for x in target)
    assert K.sign(pairing(K, target, I, target)) > 0
    assert K.sign(pairing(K, target, I, center)) < 0
    target_value = dot(F, cf, tuple(lift(x) for x in target))
    assert F.sign(target_value) < 0
    return {'status': 'PASS_EXACT_NEW_CAP', 'source': GEOMETRY,
            'positive_parameter_interval': [3, 4], 'norm': 2,
            'all_44_pairings': rows, 'cuts_persistent_target3': True,
            'positive_common_center': True}, (F, cf, lift, G, I, center, old_fields, new_fields)


def embed_K(field, value):
    if field is K:
        return value
    return field.lift(embed_K(field.base, value))


def check_archive(saved, result, context):
    F, cf, lift, G, I, center, old_fields, new_fields = context
    E = F.base
    assert saved['geometry_source'] == GEOMETRY
    geometry = read(GEOMETRY)
    source = read(geometry['source'])
    case = geometry['cases'][0]
    assert saved['base_source'] == source['base_source']
    assert saved['target_source'] == source['target_source']
    assert saved['facet_count'] == 45
    assert [parse(saved[k]) for k in ('A', 'B', 'C')] == [
        (Q(-20), Q(25, 2)), (Q(-25), Q(15)), (Q(-42, 5), Q(24, 5))]
    assert tuple(map(F.parse, saved['seed_coefficients'])) == cf
    pa = matvec(F, tuple(tuple(lift(x) for x in row) for row in G), cf)
    assert tuple(map(F.parse, saved['seed_pairings'])) == pa
    assert E.parse(saved['positive_s']) == F.neg(pa[4])[1]
    assert F.parse(saved['norm']) == dot(F, cf, pa)
    assert F.parse(saved['center_inequality']) == dot(F, cf, tuple(lift(x) for x in center))
    assert parse(saved['D']) == E.d
    assert parse(saved['delta_squared']) == F.d[0] and F.d[1] == K.zero
    labels = lambda row: 'right' if row.get('label') == 2 else 'pi/3' if row.get('label') == 3 else row['kind']
    saved_relations = [labels(row) for key in ('seed_checks', 'old_checks', 'new_checks') for row in saved[key]]
    assert saved_relations == [r['relation'] for r in result['all_44_pairings']]
    assert [r['index'] for r in saved['seed_checks']] == list(range(10))
    assert [r['index'] for r in saved['old_checks']] == case['retained_old_caps']
    assert [r['reference'] for r in saved['new_checks']] == [
        geometry['cap_references'][j] for j in case['selected_new_caps']]
    targets = read(saved['target_source'])['targets']
    signs = [F.sign(dot(F, cf, tuple(lift(parse(x)) for x in row['certificate']['y']))) for row in targets]
    assert signs == saved['target_signs']
    assert saved['target_cuts'] == [i for i, sign in enumerate(signs) if sign < 0] == [3]

    assert saved['volume']['status'] == 'EXACT_STRICT_SPACELIKE_ESCAPE'
    cert = saved['volume']['certificate']
    y = tuple(map(parse, cert['y']))
    assert len(y) == 10 and all(K.sign(x) > 0 for x in y)
    norm, cp = pairing(K, y, I, y), pairing(K, y, I, center)
    assert K.sign(norm) > 0 and K.sign(cp) < 0
    assert norm == parse(cert['norm']) and cp == parse(cert['center_pairing'])
    ss = [f.sign(dot(f, cc, tuple(embed_K(f, x) for x in y))) for f, cc in old_fields + new_fields]
    ss.append(F.sign(dot(F, cf, tuple(lift(x) for x in y))))
    assert ss == cert['nonseed_inequality_signs'] == [1]*35

    # Independently count the literal preceding entries and verify zero-pattern
    # incompatibility with this cap; this implies distinct mirrors at any scale.
    old = load_old(G)
    old += [node(r, G) for r in read('data/enumerations/fresh_h9_joint_completion_20261013.json')['caps']]
    count = len(old)
    for n in old:
        f = Quadratic(K, n['d'])
        assert [j for j, x in enumerate(zip(n['pa'], n['pb'])) if f.sign(x)] != [4, 8]
    for suffix in ('targeted', 'adaptive', 'second'):
        src = read(f'data/enumerations/fresh_h9_quartic_completion_20261013_{suffix}.json')
        e = Quadratic(K, parse(src['anchor_radicand']))
        count += len(src['caps'])
        for row in src['caps']:
            f = Quadratic(e, e.parse(row['discriminant_over_E']))
            assert [j for j, x in enumerate(map(f.parse, row['seed_pairings'])) if f.sign(x)] != [4, 8]
    assert count == saved['mirror_check']['preceding_entries'] == 3030
    assert saved['mirror_check']['matching_zero_pattern'] == 0
    return {'strict_future_spacelike_escape': True, 'wall_count': 45,
            'distinct_from_literal_preceding_entries': count,
            'volume_conclusion': 'INFINITE_VOLUME'}


def main():
    result, context = reconstruct()
    saved = read(SOURCE)
    result.update(check_archive(saved, result, context))
    for kind in ('cap', 'escape'):
        bad = copy.deepcopy(saved)
        if kind == 'cap':
            bad['seed_coefficients'][0][0][0]['a'] = '1'
        else:
            bad['volume']['certificate']['y'] = [{'a': '0', 'b': '0'}]*10
        try:
            check_archive(bad, result, context)
        except AssertionError:
            result[f'forged_{kind}_rejected'] = True
        else:
            raise AssertionError('Forged certificate accepted')
    result['status'] = 'PASS_EXACT_NEW_CAP_AND_INFINITE_VOLUME'
    OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'all_44_pairings'}))


if __name__ == '__main__':
    main()
