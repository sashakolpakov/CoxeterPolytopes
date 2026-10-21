#!/usr/bin/env python3
"""Exact replay without numerical proposal, selection, or geometry imports."""
import argparse
from copy import deepcopy
from itertools import combinations
import json
from scripts.fresh_h9_gap_completion_20261013_base import (
    current, configuration, profiles, SOURCE, PREFIX, canonical_anchor, cap_profile, read, ROOT, K,
    O, Z, P, parse, add, scale, mul, enc, matvec, dot, pairing, lift_K,
    relation, angle_pool, ANGLE_SOURCE, strict_escape)
from scripts.fresh_h9_replaced_core_20261013_verify import census_check


def verify_cap(record, name, G, I, h, ps, pool):
    row = cap_profile(record, name)
    F, E = row['field'], row['field'].base
    anchor_field, anchor = canonical_anchor(record['anchor_reference'])
    assert E.d == anchor_field.d
    line0, line1 = [tuple(map(E.parse, record[k])) for k in ('line_constant', 'line_direction')]
    free = record['free_seed_indices']
    assert len(free) == len(set(free)) == 3 and all(0 <= i < 10 for i in free)
    assert E.sign(line0[free[0]]) == 0 and E.sign(E.sub(line1[free[0]], E.one)) == 0
    fixed = dict(zip(record['fixed_seed_indices'], record['fixed_seed_labels'], strict=True))
    assert set(fixed).isdisjoint(free) and set(fixed) <= set(range(10))
    assert all(m in (3, 5) for m in fixed.values())
    for i in range(10):
        if i in free:
            continue
        expected = E.lift(O if fixed.get(i) == 3 else P if fixed.get(i) == 5 else Z)
        assert E.sign(E.sub(line0[i], expected)) == E.sign(line1[i]) == 0
    M = E.lift({2: Z, 3: O, 5: P}[record['anchor_label']])
    assert record['polar_label'] in (2, 10)
    N = E.lift(O if record['polar_label'] == 10 else Z)
    polar = tuple(E.lift(t) for t in ps['old348']['c'])
    for c, target in ((anchor, M), (polar, N)):
        assert E.sign(E.sub(dot(E, c, line0), target)) == 0
        assert E.sign(dot(E, c, line1)) == 0
    IE = tuple(tuple(E.lift(t) for t in r) for r in I)
    A = pairing(E, line1, IE, line1)
    B = E.mul(E.rational(2), pairing(E, line0, IE, line1))
    C = E.sub(pairing(E, line0, IE, line0), E.rational(2))
    saved = tuple(map(E.parse, record['norm_polynomial']))
    assert all(E.sign(E.sub(x, y)) == 0 for x, y in zip((A, B, C), saved, strict=True))
    D = E.sub(E.mul(B, B), E.mul(E.rational(4), E.mul(A, C)))
    assert E.sign(E.sub(F.d, D)) == 0 and E.sign(D) > 0 and E.sign(A)
    x = F.parse(record['parameter'])
    branch = record['branch']
    assert branch in (-1, 1)
    expected = (E.div(E.neg(B), E.mul(E.rational(2), A)),
                E.div(E.rational(branch), E.mul(E.rational(2), A)))
    assert F.sign(F.sub(x, expected)) == 0
    for i in range(10):
        p = F.neg(F.add(F.lift(line0[i]), F.mul(x, F.lift(line1[i]))))
        assert F.sign(F.sub(p, row['p'][i])) == 0
    GF = tuple(tuple(lift_K(F, t) for t in r) for r in G)
    assert all(F.sign(F.sub(x, y)) == 0 for x, y in
               zip(matvec(F, GF, row['c']), row['p'], strict=True))
    assert F.sign(F.sub(dot(F, row['c'], row['p']), F.rational(2))) == 0
    assert F.sign(dot(F, row['c'], tuple(lift_K(F, t) for t in h))) > 0
    seed = {str(i): relation(row, ps[f'seed{i}'], pool) for i in range(10)}
    assert seed == record['seed_relations']
    assert all(r['kind'] in ('finite', 'dotted') for r in seed.values())
    assert relation(row, ps[record.get('anchor_name', 'cap2')], pool) == {'kind': 'finite', 'label': record['anchor_label']}
    assert relation(row, ps['old348'], pool) == {'kind': 'finite', 'label': record['polar_label']}
    return row


def run(source):
    data = read(source)
    bank = read(data['source'])
    G, I, h, old, base_source, base_case = configuration(bank.get('configuration', '45'))
    pool = angle_pool(ANGLE_SOURCE)
    selected_indices = sorted(set(i for c in data['cases'] for i in c['selected_caps']))
    rows = {i: verify_cap(bank['caps'][i], f'{bank.get("cap_name_prefix", "gap")}{i}', G, I, h, old, pool)
            for i in selected_indices}
    targets = [tuple(map(parse, e['certificate']['y']))
               for e in read(base_source)['cases'][base_case]['geometry']['escapes']]
    for i, row in rows.items():
        record, F = bank['caps'][i], row['field']
        checks = {n: relation(row, r, pool) for n, r in old.items()}
        assert checks == record['pair_checks']
        dropped = [n for n in old if checks[n]['kind'] not in ('finite', 'dotted')]
        assert dropped == record['dropped'] and len(dropped) <= bank['drop_budget']
        assert set(dropped) <= set(bank['removable'])
        signs = [F.sign(dot(F, row['c'], tuple(lift_K(F, t) for t in y))) for y in targets]
        assert signs == record['target_signs']
        assert [j for j, s in enumerate(signs) if s < 0] == record['cuts']
        assert signs[bank['focus_target']] < 0
    counts = {'configurations': 0, 'all_pair_checks': 0, 'affine_inertia_checks': 0,
              'strict_escape_checks': 0}
    for case in data['cases']:
        ids = case['selected_caps']
        dropped = sorted(set().union(*(set(bank['caps'][i]['dropped']) for i in ids)))
        assert dropped == case['dropped']
        selected = [r for n, r in old.items() if n not in dropped]+[rows[i] for i in ids]
        assert [r['name'] for r in selected] == case['walls']
        assert len(selected) == case['facet_count']
        pairs = {'|'.join(sorted((a['name'], b['name']))): relation(a, b, pool)
                 for a, b in combinations(selected, 2)}
        assert pairs == case['pair_checks']
        assert all(r['kind'] in ('finite', 'dotted') for r in pairs.values())
        counts['all_pair_checks'] += len(pairs)
        counts['affine_inertia_checks'] += census_check(case['walls'], pairs, case['affine_census'])
        if 'geometry' in case:
            assert not case['affine_census']['affine_circuits']
            assert case['status'] == 'EXACT_STRICT_ESCAPES'
            assert case['geometry']['escapes']
            for escape in case['geometry']['escapes']:
                cert = escape['certificate']
                assert strict_escape(tuple(map(parse, cert['y'])), I, h, selected) == cert
                counts['strict_escape_checks'] += 1
        else:
            assert case['status'] == 'REJECTED_POSITIVE_AFFINE_COMPACTNESS'
            assert case['affine_census']['affine_circuits']
        counts['configurations'] += 1
    return {'status': 'PASS_EXACT_GAP_COMPLETION', 'source': source,
            'verified_caps': selected_indices, **counts}


def verify_swap(source):
    data = read(source)
    G, I, h, ps, prior, prior_case = configuration('47')
    _, _, _, original = profiles()
    original.update(ps)
    perm = data['permutation']
    assert perm == [0, 1, 2, 3, 4, 5, 7, 6, 8, 9]
    assert all(G[i][j] == G[perm[i]][perm[j]] for i in range(10) for j in range(10))
    assert all(h[i] == h[perm[i]] for i in range(10))
    pool = angle_pool(ANGLE_SOURCE)
    targets = [tuple(map(parse, e['certificate']['y']))
               for e in read(prior)['cases'][prior_case]['geometry']['escapes']]
    counts = {'configurations': 0, 'all_pair_checks': 0, 'affine_inertia_checks': 0,
              'strict_escape_checks': 0}
    for case in data['cases']:
        source_name = case['transported_source']
        r = original[source_name]
        f = r['field']
        row = dict(r, name='swap_'+source_name,
                   c=tuple(r['c'][i] for i in perm), p=tuple(r['p'][i] for i in perm))
        assert f.sign(f.sub(dot(f, row['c'], row['p']), lift_K(f, row['q']))) == 0
        assert f.sign(dot(f, row['c'], tuple(lift_K(f, t) for t in h))) > 0
        selected = [r for n, r in ps.items() if n not in case['dropped']]+[row]
        assert [r['name'] for r in selected] == case['walls']
        pairs = {'|'.join(sorted((a['name'], b['name']))): relation(a, b, pool)
                 for a, b in combinations(selected, 2)}
        assert pairs == case['pair_checks']
        assert all(r['kind'] in ('finite', 'dotted') for r in pairs.values())
        counts['all_pair_checks'] += len(pairs)
        signs = [f.sign(dot(f, row['c'], tuple(lift_K(f, t) for t in y))) for y in targets]
        assert signs == case['target_signs']
        counts['affine_inertia_checks'] += census_check(case['walls'], pairs, case['affine_census'])
        if 'geometry' in case:
            assert not case['affine_census']['affine_circuits']
            assert case['status'] == 'EXACT_STRICT_ESCAPES'
            for e in case['geometry']['escapes']:
                assert strict_escape(tuple(map(parse, e['certificate']['y'])), I, h, selected) == e['certificate']
                counts['strict_escape_checks'] += 1
        else:
            assert case['status'] == 'REJECTED_POSITIVE_AFFINE_COMPACTNESS'
            assert case['affine_census']['affine_circuits']
        counts['configurations'] += 1
    return {'status': 'PASS_EXACT_SEED_SYMMETRY_REPLACEMENTS', 'source': source, **counts}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', default=PREFIX+'_geometry.json')
    p.add_argument('--verify', action='store_true')
    p.add_argument('--swap', action='store_true')
    a = p.parse_args()
    if a.swap and a.source == PREFIX+'_geometry.json':
        a.source = PREFIX+'_swap.json'
    result = verify_swap(a.source) if a.swap else run(a.source)
    output = ROOT/a.source.replace('.json', '_verify.json')
    if a.verify:
        assert json.loads(output.read_text()) == result
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result, flush=True)
