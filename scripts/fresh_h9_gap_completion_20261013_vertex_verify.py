#!/usr/bin/env python3
"""Pure exact replay of new finite vertices and selected chamber geometry."""
import argparse
import json
from itertools import combinations
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K
from scripts.fresh_h9_replaced_core_20261013_verify import census_check


def row_from_record(record, name):
    f = Quadratic(K, parse(record['radicand']))
    return {'name': name, 'field': f,
            'p': tuple((Z, parse(t)) for t in record['raw_seed_pairings']),
            'c': tuple((Z, parse(t)) for t in record['raw_seed_coefficients']),
            'q': scale(O, 2)}


def run():
    G, I, h, old, source, case = configuration('48')
    bank = read(PREFIX+'_vertex_caps.json')
    assert bank['source'] == source and bank['source_case'] == case
    pool = angle_pool(ANGLE_SOURCE)
    rows, vertices = [], []
    for i, record in enumerate(bank['caps']):
        row = row_from_record(record, f'vertex{i}')
        f, block = row['field'], record['block_names']
        assert len(block) == len(set(block)) == 8
        assert all(old[n]['field'] is K for n in block)
        raw = tuple(map(parse, record['raw_seed_pairings']))
        coeff = tuple(map(parse, record['raw_seed_coefficients']))
        assert matvec(K, I, raw) == coeff
        assert K.sign(pairing(K, raw, I, raw)) > 0
        assert mul(f.d, dot(K, raw, coeff)) == scale(O, 2)
        assert f.sign(dot(f, row['c'], tuple(lift_K(f, t) for t in h))) > 0
        pairs = {n: relation(row, r, pool) for n, r in old.items()}
        assert pairs == record['pair_checks']
        assert all(r == ({'kind': 'finite', 'label': 2} if n in block else {'kind': 'dotted'})
                   for n, r in pairs.items())
        proof = record['finite_vertex']
        gram = [[dot(K, old[a]['c'], old[b]['p']) for b in block] for a in block]
        assert [[enc(t) for t in r] for r in gram] == proof['block_gram']
        assert inertia_K(gram) == (8, 0, 0)
        y = tuple(map(parse, proof['y']))
        norm, cp = pairing(K, y, I, y), pairing(K, y, I, h)
        assert K.sign(norm) < 0 and K.sign(cp) < 0
        assert enc(norm) == proof['norm'] and enc(cp) == proof['center_pairing']
        signs = {n: r['field'].sign(dot(r['field'], r['c'],
                 tuple(lift_K(r['field'], t) for t in y))) for n, r in old.items()}
        assert signs == proof['old_wall_signs']
        assert all(s == (0 if n in block else 1) for n, s in signs.items())
        assert dot(K, coeff, y) == Z and proof['new_cap_sign'] == 0
        rows.append(row)
        vertices.append(y)
    geo_source = PREFIX+'_vertex_geometry.json'
    data = read(geo_source)
    allpairs = {f'{i}|{j}': relation(rows[i], rows[j], pool)
                for i, j in combinations(range(len(rows)), 2)}
    assert allpairs == data['all_cap_pair_relations']
    counts = {'verified_caps': len(rows), 'exact_new_finite_vertices': len(vertices),
              'cap_pair_checks': len(allpairs), 'configurations': 0,
              'all_pair_checks': 0, 'affine_inertia_checks': 0, 'strict_escape_checks': 0,
              'joint_finite_vertex_checks': 0}
    for case in data['cases']:
        ids = case['selected_caps']
        selected = list(old.values())+[rows[i] for i in ids]
        names = [r['name'] for r in selected]
        assert names == case['walls'] and len(names) == case['facet_count']
        pairs = {'|'.join(sorted((a['name'], b['name']))): relation(a, b, pool)
                 for a, b in combinations(selected, 2)}
        assert pairs == case['pair_checks']
        assert all(r['kind'] in ('finite', 'dotted') for r in pairs.values())
        counts['all_pair_checks'] += len(pairs)
        counts['affine_inertia_checks'] += census_check(names, pairs, case['affine_census'])
        for i in ids:
            y = vertices[i]
            for j in ids:
                if i == j:
                    continue
                r, f = rows[j], rows[j]['field']
                assert f.sign(dot(f, r['c'], tuple(lift_K(f, t) for t in y))) > 0
            counts['joint_finite_vertex_checks'] += 1
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
    return {'status': 'PASS_EXACT_FINITE_VERTEX_CAPS_AND_GEOMETRY', 'source': geo_source, **counts}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    result = run()
    output = ROOT/(PREFIX+'_vertex_verify.json')
    if a.verify:
        assert json.loads(output.read_text()) == result
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result, flush=True)
