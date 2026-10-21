#!/usr/bin/env python3
"""Producer-free exact replay of changed-core chambers and volume witnesses."""
import argparse
import json
from itertools import combinations
from scripts.fresh_h9_replaced_core_20261013_base import *
from scripts.fresh_h9_replaced_core_20261013_cores import normalized_field, inertia
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K


def cap_check(record, G, I, center, ps, pool, name):
    row = cap_profile(record, name)
    F, E = row['field'], row['field'].base
    R = parse(record['base_radicand'])
    polar = ps['old348']
    factors = {2: Z, 3: scale(O, Q(1, 2)), 4: O,
               5: scale(add(O, P), Q(1, 2)), 6: scale(O, Q(3, 2)),
               10: scale(add(scale(O, 2), P), Q(1, 2))}
    full_anchor = record.get('kind') == 'full_quartic_anchor'
    if full_anchor:
        anchor_field, anchor_c = collapsed_anchor(record['anchor_reference'])
        assert R == anchor_field.d
    else:
        assert R == mul(polar['q'], factors[record['polar_label']])
    line0 = tuple(map(E.parse, record['line_constant']))
    line1 = tuple(map(E.parse, record['line_direction']))
    x = F.parse(record['parameter'])
    free = record['free_seed_indices']
    assert len(free) == (3 if full_anchor else 2)
    assert len(set(free)) == len(free) and all(0 <= i < 10 for i in free)
    u = free[0]
    assert E.sign(E.sub(line1[u], E.one)) == 0 and E.sign(line0[u]) == 0
    fixed = dict(zip(record['fixed_seed_indices'], record['fixed_seed_labels'], strict=True))
    assert set(fixed).isdisjoint(free)
    for i in range(10):
        if i in free:
            continue
        expected = E.lift(O if fixed.get(i) == 3 else P if fixed.get(i) == 5 else Z)
        assert E.sign(E.sub(line0[i], expected)) == 0
        assert E.sign(line1[i]) == 0
    if full_anchor:
        assert record['polar_label'] in (2, 10)
        target = E.lift(O if record['polar_label'] == 10 else Z)
        anchor_target = E.lift({2: Z, 3: O, 5: P}[record['anchor_label']])
        assert E.sign(E.sub(dot(E, anchor_c, line0), anchor_target)) == 0
        assert E.sign(dot(E, anchor_c, line1)) == 0
    else:
        target = (Z, O) if record['polar_label'] != 2 else E.zero
    pc = tuple(E.lift(t) for t in polar['c'])
    assert E.sign(E.sub(dot(E, pc, line0), target)) == 0
    assert E.sign(dot(E, pc, line1)) == 0
    IE = tuple(tuple(E.lift(t) for t in r) for r in I)
    A = pairing(E, line1, IE, line1)
    B = E.mul(E.rational(2), pairing(E, line0, IE, line1))
    C = E.sub(pairing(E, line0, IE, line0), E.rational(2))
    saved_poly = tuple(map(E.parse, record['norm_polynomial']))
    assert all(E.sign(E.sub(a, b)) == 0 for a, b in zip((A, B, C), saved_poly, strict=True))
    discriminant = E.sub(E.mul(B, B), E.mul(E.rational(4), E.mul(A, C)))
    assert E.sign(E.sub(discriminant, F.d)) == 0 and E.sign(F.d) > 0
    for i in range(10):
        expected = F.neg(F.add(F.lift(line0[i]), F.mul(x, F.lift(line1[i]))))
        assert F.sign(F.sub(row['p'][i], expected)) == 0
    GG = tuple(tuple(lift_K(F, t) for t in r) for r in G)
    assert all(F.sign(F.sub(a, b)) == 0
               for a, b in zip(row['p'], matvec(F, GG, row['c']), strict=True))
    assert F.sign(F.sub(dot(F, row['c'], row['p']), F.rational(2))) == 0
    assert all(F.sign(t) <= 0 for t in row['p'])
    assert F.sign(dot(F, row['c'], tuple(lift_K(F, t) for t in center))) > 0
    actual_seed = {str(i): relation(row, ps[f'seed{i}'], pool) for i in range(10)}
    assert actual_seed == record['seed_relations']
    assert all(r['kind'] in ('finite', 'dotted') for r in actual_seed.values())
    assert relation(row, polar, pool) == {'kind': 'finite', 'label': record['polar_label']}
    if full_anchor:
        ref = record['anchor_reference']
        anchor = cap_profile(read(ref['source'])['caps'][ref['index']], 'full_anchor')
        assert relation(row, anchor, pool) == {'kind': 'finite', 'label': record['anchor_label']}
    return row


def census_check(names, pairs, saved):
    """Check an archived frontier for exact inertia and exhaustive closure."""
    f, weights = normalized_field({r['label'] for r in pairs.values() if r['kind'] == 'finite'})
    def rel(a, b):
        return pairs['|'.join(sorted((a, b)))]
    rows = {tuple(r['walls']): tuple(r['inertia']) for r in saved['exact_inertias']}
    assert len(rows) == saved['connected_tested_subsets'] == len(saved['exact_inertias'])
    positive, affine = set(), []
    for ids, sig in rows.items():
        assert tuple(sorted(set(ids))) == ids and set(ids) <= set(names)
        seen = {ids[0]}
        while True:
            added = {i for i in ids if i not in seen and any(
                rel(i, j)['kind'] == 'finite' and rel(i, j)['label'] != 2 for j in seen)}
            if not added:
                break
            seen |= added
        assert seen == set(ids)
        matrix = [[f.rational(2) if a == b else weights[rel(a, b)['label']]
                   for b in ids] for a in ids]
        assert inertia(f, matrix) == sig
        if sig == (len(ids), 0, 0):
            positive.add(ids)
        elif not sig[1]:
            assert sig == (len(ids)-1, 0, 1)
            affine.append(list(ids))
    for name in names:
        assert (name,) in rows
    for ids in positive:
        for name in names:
            if name in ids:
                continue
            rs = [rel(name, other) for other in ids]
            if all(r['kind'] == 'finite' for r in rs) and any(r['label'] != 2 for r in rs):
                assert tuple(sorted((*ids, name))) in rows
    for ids in rows:
        if len(ids) > 1:
            assert any(tuple(n for n in ids if n != removed) in positive for removed in ids)
    assert sorted(affine) == saved['affine_circuits']
    return len(rows)


def run(source):
    data = read(source)
    capdata = read(data['source'])
    G, I, h, ps = profiles()
    assert inertia_K(G) == (9, 1, 0)
    assert all(matvec(K, I, row) == tuple(O if i == j else Z for j in range(10))
               for i, row in enumerate(G))
    assert K.sign(pairing(K, h, I, h)) < 0 and all(K.sign(t) > 0 for t in h)
    pool = angle_pool(ANGLE_SOURCE)
    verified = set()
    for i, ref in enumerate(capdata.get('required_caps', [])):
        r = read(ref['source'])['caps'][ref['index']]
        ps[f'required{i}'] = cap_check(r, G, I, h, ps, pool, f'required{i}')
        verified.add((ref['source'], ref['index']))
    counts = {'configurations': 0, 'all_pair_checks': 0, 'affine_inertia_checks': 0,
              'strict_escape_checks': 0}
    for case in data['cases']:
        for i in case['selected_caps']:
            ps[f'cap{i}'] = cap_check(capdata['caps'][i], G, I, h, ps, pool, f'cap{i}')
            verified.add((data['source'], i))
        names = case['walls']
        assert len(names) == len(set(names)) == case['facet_count']
        for name in names:
            r = ps[name]
            f = r['field']
            assert f.sign(f.sub(dot(f, r['c'], r['p']), lift_K(f, r['q']))) == 0
            assert f.sign(dot(f, r['c'], tuple(lift_K(f, t) for t in h))) > 0
        pairs = {'|'.join(sorted((a, b))): relation(ps[a], ps[b], pool)
                 for a, b in combinations(names, 2)}
        assert pairs == case['pair_checks']
        assert all(r['kind'] in ('finite', 'dotted') for r in pairs.values())
        counts['all_pair_checks'] += len(pairs)
        counts['affine_inertia_checks'] += census_check(names, pairs, case['affine_census'])
        if 'geometry' in case:
            assert not case['affine_census']['affine_circuits']
            assert case['status'] == 'EXACT_STRICT_ESCAPES'
            assert case['geometry']['escapes']
            for escape in case['geometry']['escapes']:
                y = tuple(map(parse, escape['certificate']['y']))
                assert strict_escape(y, I, h, [ps[n] for n in names]) == escape['certificate']
                counts['strict_escape_checks'] += 1
        else:
            assert case['status'] == 'REJECTED_POSITIVE_AFFINE_COMPACTNESS'
            assert case['affine_census']['affine_circuits']
        counts['configurations'] += 1
    return {'status': 'PASS_EXACT_CHANGED_CORE_GEOMETRY', 'source': source,
            'verified_cap_sources': [{'source': s, 'index': i} for s, i in sorted(verified)],
            **counts, 'scope': 'Independent of numerical proposal/selection/volume modules. '
            'Checks selected cap equations, all pairs, exhaustive affine-census closure, '
            'common strict timelike center, and all strict future spacelike escapes.'}


def mutation_checks(source):
    from copy import deepcopy
    data = read(source)
    bank = read(data['source'])
    case = next(c for c in data['cases'] if 'geometry' in c)
    G, I, h, ps = profiles()
    pool = angle_pool(ANGLE_SOURCE)
    record = deepcopy(bank['caps'][case['selected_caps'][0]])
    row = cap_profile(record, 'forged')
    field = row['field']
    record['seed_coefficients'][0] = field.encode(field.add(row['c'][0], field.one))
    try:
        cap_check(record, G, I, h, ps, pool, 'forged')
    except AssertionError:
        pass
    else:
        raise AssertionError('Forged cap was accepted')
    try:
        strict_escape(h, I, h, [])
    except AssertionError:
        pass
    else:
        raise AssertionError('Timelike escape replacement was accepted')
    return {'forged_cap': 'REJECTED', 'timelike_escape_replacement': 'REJECTED'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', default='data/enumerations/fresh_h9_replaced_core_20261013_selection.json')
    p.add_argument('--verify', action='store_true')
    p.add_argument('--self-test', action='store_true')
    a = p.parse_args()
    result = run(a.source)
    output = ROOT/a.source.replace('.json', '_verify.json')
    if a.verify:
        assert json.loads(output.read_text()) == result
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print({k: v for k, v in result.items() if k != 'verified_cap_sources'}, flush=True)
    if a.self_test:
        print(mutation_checks(a.source), flush=True)
