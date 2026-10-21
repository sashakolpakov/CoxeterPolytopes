#!/usr/bin/env python3
"""Exact coupled-cap and joint-finite-vertex replay, without search imports."""
import argparse
import json
from itertools import combinations
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h9_gap_completion_20261013_vertex_verify import row_from_record
from scripts.fresh_h9_replaced_core_20261013_verify import census_check
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K


def point_pairing(row, field, y):
    common, er = embed(field, row['field'])
    return common, dot(common, tuple(er(x) for x in row['c']),
                       tuple(lift_existing(common, field, x) for x in y))


def run():
    G, I, h, old, source, case = configuration('48')
    vb = read(PREFIX+'_vertex_caps.json')
    cb = read(PREFIX+'_coupled_vertices.json')
    generic = [row_from_record(r, f'vertex{i}') for i, r in enumerate(vb['caps'])]
    new = []
    pool = angle_pool(ANGLE_SOURCE)
    for i, record in enumerate(cb['caps']):
        row = cap_profile(record, f'coupled{i}')
        F, E = row['field'], row['field'].base
        anchor = generic[record['anchor_cap_index']]
        assert key(E) == key(anchor['field'])
        block = record['block_names']
        assert len(block) == len(set(block)) == 8
        assert all(old[n]['field'] is K for n in block)
        gram = [[dot(K, old[a]['c'], old[b]['p']) for b in block] for a in block]
        assert inertia_K(gram) == (8, 0, 0)
        l0, l1 = [tuple(map(E.parse, record[k])) for k in ('line_constant', 'line_direction')]
        for n in block:
            c = tuple(E.lift(t) for t in old[n]['c'])
            assert E.sign(dot(E, c, l0)) == E.sign(dot(E, c, l1)) == 0
        k = {2: Z, 3: O, 5: P}[record['anchor_label']]
        assert E.sign(E.add(dot(E, anchor['c'], l0), E.lift(k))) == 0
        assert E.sign(dot(E, anchor['c'], l1)) == 0
        IE = tuple(tuple(E.lift(t) for t in r) for r in I)
        A = pairing(E, l1, IE, l1)
        B = E.mul(E.rational(2), pairing(E, l0, IE, l1))
        C = E.sub(pairing(E, l0, IE, l0), E.rational(2))
        assert all(E.sign(E.sub(x, y)) == 0 for x, y in
                   zip((A, B, C), map(E.parse, record['norm_polynomial']), strict=True))
        D = E.sub(E.mul(B, B), E.mul(E.rational(4), E.mul(A, C)))
        assert E.sign(E.sub(D, F.d)) == 0 and E.sign(D) > 0 and E.sign(A)
        x = F.parse(record['parameter'])
        assert record['branch'] in (-1, 1)
        expected = (E.div(E.neg(B), E.mul(E.rational(2), A)),
                    E.div(E.rational(record['branch']), E.mul(E.rational(2), A)))
        assert F.sign(F.sub(x, expected)) == 0
        p = tuple(F.add(F.lift(a), F.mul(x, F.lift(b))) for a, b in zip(l0, l1, strict=True))
        assert all(F.sign(F.sub(a, b)) == 0 for a, b in zip(p, row['p'], strict=True))
        GF = tuple(tuple(lift_K(F, t) for t in r) for r in G)
        assert all(F.sign(F.sub(a, b)) == 0 for a, b in
                   zip(matvec(F, GF, row['c']), row['p'], strict=True))
        assert F.sign(F.sub(dot(F, row['c'], row['p']), F.rational(2))) == 0
        assert F.sign(dot(F, row['c'], tuple(lift_K(F, t) for t in h))) > 0
        checks = {n: relation(row, r, pool) for n, r in old.items()}
        assert checks == record['pair_checks']
        assert all(r == ({'kind': 'finite', 'label': 2} if n in block else {'kind': 'dotted'})
                   for n, r in checks.items())
        assert relation(row, anchor, pool) == {'kind': 'finite', 'label': record['anchor_label']}
        new.append(row)
    geometry_source = PREFIX+'_coupled_geometry.json'
    data = read(geometry_source)
    counts = {'verified_coupled_caps': len(new), 'configurations': 0, 'all_pair_checks': 0,
              'affine_inertia_checks': 0, 'joint_finite_vertices': 0, 'strict_escape_checks': 0}
    for case in data['cases']:
        selected = list(old.values())+[generic[i] for i in case['selected_generic']]+[new[i] for i in case['selected_coupled']]
        names = [r['name'] for r in selected]
        assert names == case['walls'] and len(names) == case['facet_count']
        ps = {r['name']: r for r in selected}
        pairs = {'|'.join(sorted((a['name'], b['name']))): relation(a, b, pool)
                 for a, b in combinations(selected, 2)}
        assert pairs == case['pair_checks']
        assert all(r['kind'] in ('finite', 'dotted') for r in pairs.values())
        counts['all_pair_checks'] += len(pairs)
        counts['affine_inertia_checks'] += census_check(names, pairs, case['affine_census'])
        for proof in case['joint_finite_vertices']:
            F = ps[proof['point_field_from_cap']]['field']
            y = tuple(map(F.parse, proof['y']))
            IF = tuple(tuple(lift_K(F, t) for t in r) for r in I)
            norm, cp = pairing(F, y, IF, y), pairing(F, y, IF, tuple(lift_K(F, t) for t in h))
            assert F.sign(norm) < 0 and F.sign(cp) < 0
            assert F.encode(norm) == proof['norm'] and F.encode(cp) == proof['center_pairing']
            signs = {}
            for r in selected:
                field, value = point_pairing(r, F, y)
                signs[r['name']] = field.sign(value)
            assert signs == proof['all_wall_signs']
            zeros = proof['zero_walls']
            assert len(zeros) == len(set(zeros)) == 9
            assert all(s == (0 if n in zeros else 1) for n, s in signs.items())
            # A timelike common point and nine independent incident normals
            # certify an actual finite vertex, independent of the projection recipe.
            cap_index = int(proof['point_field_from_cap'].removeprefix('coupled'))
            record = cb['caps'][cap_index]
            anchor_id = record['anchor_cap_index']
            shared = sorted(set(record['block_names']) & set(vb['caps'][anchor_id]['block_names']))
            order = shared+[f'vertex{anchor_id}', proof['point_field_from_cap']]
            assert len(shared) == 7 and sorted(order) == zeros
            physical = []
            for n in order:
                common, er = embed(F, ps[n]['field'])
                assert key(common) == key(F)
                physical.append((tuple(er(t) for t in ps[n]['c']), tuple(er(t) for t in ps[n]['p'])))
            saved_gram = [[parse(t) for t in r] for r in proof['gram']]
            assert inertia_K(saved_gram) == (9, 0, 0)
            actual_gram = [[dot(F, a[0], b[1]) for b in physical] for a in physical]
            assert all(F.sign(F.sub(actual_gram[i][j], lift_K(F, saved_gram[i][j]))) == 0
                       for i in range(9) for j in range(9))
            counts['joint_finite_vertices'] += 1
        if 'geometry' in case:
            assert not case['affine_census']['affine_circuits']
            assert case['status'] == 'EXACT_STRICT_ESCAPES'
            for escape in case['geometry']['escapes']:
                proof = escape['certificate']
                assert strict_escape(tuple(map(parse, proof['y'])), I, h, selected) == proof
                counts['strict_escape_checks'] += 1
        else:
            assert case['status'] == 'REJECTED_POSITIVE_AFFINE_COMPACTNESS'
            assert case['affine_census']['affine_circuits']
        counts['configurations'] += 1
    return {'status': 'PASS_EXACT_FINITE_ANGLE_COUPLED_GEOMETRY', 'source': geometry_source, **counts}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    result = run()
    output = ROOT/(PREFIX+'_coupled_verify.json')
    if a.verify:
        assert json.loads(output.read_text()) == result
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result, flush=True)
