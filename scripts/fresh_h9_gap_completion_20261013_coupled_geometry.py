#!/usr/bin/env python3
"""Joint caps including exact finite-angle attachments between new facets."""
import json
from itertools import combinations
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h9_gap_completion_20261013_vertex_caps import decode
from scripts.fresh_h4_blocks_independent_verify import inverse
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K
from scripts.fresh_h9_replaced_core_20261013_cores import affine_census
from scripts.fresh_h9_replaced_core_20261013_geometry import volume_gate


def point_pairing(row, field, y):
    common, er = embed(field, row['field'])
    return common, dot(common, tuple(er(x) for x in row['c']),
                       tuple(lift_existing(common, field, x) for x in y))


def joint_vertex(first, new, block, label, I, h, selected):
    F = new['field']
    roots = []
    for row in block+[first, new]:
        common, er = embed(F, row['field'])
        assert key(common) == key(F)
        roots.append((tuple(er(x) for x in row['p']), tuple(er(x) for x in row['c'])))
    gram = [[dot(K, a['c'], b['p']) for b in block] + [Z, Z] for a in block]
    k = {2: Z, 3: O, 5: P}[label]
    gram += [[Z]*7+[scale(O, 2), K.neg(k)], [Z]*7+[K.neg(k), scale(O, 2)]]
    assert inertia_K(gram) == (9, 0, 0)
    gh = tuple(dot(F, c, tuple(lift_K(F, t) for t in h)) for p, c in roots)
    alpha = matvec(F, tuple(tuple(lift_K(F, t) for t in row) for row in inverse(gram)), gh)
    y = tuple(lift_K(F, t) for t in h)
    for (p, _), a in zip(roots, alpha, strict=True):
        y = tuple(F.sub(x, F.mul(a, t)) for x, t in zip(y, p, strict=True))
    IF = tuple(tuple(lift_K(F, t) for t in row) for row in I)
    norm, cp = pairing(F, y, IF, y), pairing(F, y, IF, tuple(lift_K(F, t) for t in h))
    assert F.sign(norm) < 0 and F.sign(cp) < 0
    zero = {r['name'] for r in block+[first, new]}
    signs = {}
    for row in selected:
        field, pair = point_pairing(row, F, y)
        signs[row['name']] = field.sign(pair)
    assert all(s == (0 if n in zero else 1) for n, s in signs.items())
    return {'point_field_from_cap': new['name'], 'y': list(map(F.encode, y)),
            'norm': F.encode(norm), 'center_pairing': F.encode(cp),
            'zero_walls': sorted(zero), 'all_wall_signs': signs,
            'gram': [[enc(t) for t in row] for row in gram]}


def run():
    G, I, h, old, base_source, base_case = configuration('48')
    vertices = read(PREFIX+'_vertex_caps.json')
    coupled_source = PREFIX+'_coupled_vertices.json'
    bank = read(coupled_source)
    generic = [decode(r, f'vertex{i}') for i, r in enumerate(vertices['caps'])]
    new = [cap_profile(r, f'coupled{i}') for i, r in enumerate(bank['caps'])]
    rows = generic+new
    pool = angle_pool(ANGLE_SOURCE)
    graph = [set() for _ in rows]
    relations = {}
    for i, j in combinations(range(len(rows)), 2):
        rel = relation(rows[i], rows[j], pool)
        relations[i, j] = rel
        if rel['kind'] in ('finite', 'dotted'):
            graph[i].add(j); graph[j].add(i)
    n_generic = len(generic)
    cases = []
    for anchor_index in bank['chosen_anchors']:
        attached = [i+n_generic for i, r in enumerate(bank['caps'])
                    if r['anchor_cap_index'] == anchor_index and r['anchor_label'] != 2]
        proposals = []
        for seed in attached:
            ids = [anchor_index, seed]
            available = graph[anchor_index]&graph[seed]
            while available:
                j = max(available, key=lambda j: (j >= n_generic,
                       len(graph[j]&available), -j))
                ids.append(j)
                available &= graph[j]
            ids.sort()
            proposals.append(ids)
        if not proposals:
            continue
        ids = max(proposals, key=lambda ids: (sum(i >= n_generic for i in ids), len(ids), tuple(-i for i in ids)))
        selected = list(old.values())+[rows[i] for i in ids]
        names = [r['name'] for r in selected]
        checks = dict(read(base_source)['cases'][base_case]['pair_checks'])
        for i in ids:
            record = vertices['caps'][i] if i < n_generic else bank['caps'][i-n_generic]
            checks.update({'|'.join(sorted((rows[i]['name'], n))): r for n, r in record['pair_checks'].items()})
        checks.update({'|'.join(sorted((rows[i]['name'], rows[j]['name']))): relations[i, j]
                       for i, j in combinations(ids, 2)})
        assert len(checks) == len(names)*(len(names)-1)//2
        census = affine_census(names, checks)
        case = {'anchor_cap_index': anchor_index, 'selected_generic': [i for i in ids if i < n_generic],
                'selected_coupled': [i-n_generic for i in ids if i >= n_generic],
                'walls': names, 'facet_count': len(names), 'pair_checks': checks, 'affine_census': census}
        joint = []
        for i in case['selected_coupled']:
            record = bank['caps'][i]
            anchor_id = record['anchor_cap_index']
            if anchor_id not in ids:
                continue
            shared = sorted(set(record['block_names'])&set(vertices['caps'][anchor_id]['block_names']))
            assert len(shared) == 7
            joint.append(joint_vertex(generic[anchor_id], new[i], [old[n] for n in shared],
                                      record['anchor_label'], I, h, selected))
        case['joint_finite_vertices'] = joint
        if census['affine_circuits']:
            case['status'] = 'REJECTED_POSITIVE_AFFINE_COMPACTNESS'
        else:
            case['geometry'] = volume_gate(I, h, selected)
            case['status'] = case['geometry']['status']
        cases.append(case)
        out = {'source': coupled_source, 'vertex_cap_source': PREFIX+'_vertex_caps.json',
               'base_source': base_source, 'source_case': base_case,
               'all_new_pair_relations': {f'{i}|{j}': r for (i, j), r in relations.items()},
               'cases': cases, 'scope': 'Coherent chambers with a finite angle between new '
               'vertex facets; actual joint finite-vertex and strict escape certificates. '
               'No exhaustive assembly or cap-family claim.'}
        (ROOT/(PREFIX+'_coupled_geometry.json')).write_text(json.dumps(out, indent=2)+'\n')
        print({'case': len(cases)-1, 'facets': len(names), 'coupled': case['selected_coupled'],
               'joint_vertices': len(joint), 'status': case['status'],
               'max_norm': case.get('geometry', {}).get('maximum_numeric_vertex_norm')}, flush=True)


if __name__ == '__main__':
    run()
