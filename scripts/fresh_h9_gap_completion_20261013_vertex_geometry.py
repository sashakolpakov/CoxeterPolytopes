#!/usr/bin/env python3
"""Joint exact compatibility and volume gates for finite-vertex caps."""
import json
from itertools import combinations
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h9_gap_completion_20261013_vertex_caps import decode
from scripts.fresh_h9_replaced_core_20261013_cores import affine_census
from scripts.fresh_h9_replaced_core_20261013_geometry import volume_gate


def run():
    G, I, h, old, source, case = configuration('48')
    cap_source = PREFIX+'_vertex_caps.json'
    bank = read(cap_source)
    rows = [decode(r, f'vertex{i}') for i, r in enumerate(bank['caps'])]
    pool = angle_pool(ANGLE_SOURCE)
    graph = [set() for _ in rows]
    pairs = {}
    for i, j in combinations(range(len(rows)), 2):
        r = relation(rows[i], rows[j], pool)
        pairs[i, j] = r
        if r['kind'] in ('finite', 'dotted'):
            graph[i].add(j)
            graph[j].add(i)
    maximal = []
    def bk(chosen, possible, excluded):
        if not possible and not excluded:
            maximal.append(tuple(sorted(chosen)))
            return
        pivot = max(possible|excluded, key=lambda i: len(graph[i]&possible)) if possible|excluded else None
        for i in sorted(possible-(graph[pivot] if pivot is not None else set())):
            bk(chosen+[i], possible&graph[i], excluded&graph[i])
            possible.remove(i)
            excluded.add(i)
    bk([], set(range(len(rows))), set())
    maximal.sort(key=lambda ids: (-len(ids), ids))
    print({'nodes': len(rows), 'edges': sum(map(len, graph))//2,
           'maximal_cliques': len(maximal), 'largest': len(maximal[0])}, flush=True)
    cases = []
    for ids in maximal[:4]:
        selected = list(old.values())+[rows[i] for i in ids]
        names = [r['name'] for r in selected]
        checks = dict(read(source)['cases'][case]['pair_checks'])
        for i in ids:
            checks.update({'|'.join(sorted((rows[i]['name'], n))): r
                           for n, r in bank['caps'][i]['pair_checks'].items()})
        checks.update({'|'.join(sorted((rows[i]['name'], rows[j]['name']))): pairs[i, j]
                       for i, j in combinations(ids, 2)})
        assert len(checks) == len(names)*(len(names)-1)//2
        census = affine_census(names, checks)
        result = {'selected_caps': list(ids), 'walls': names, 'facet_count': len(names),
                  'pair_checks': checks, 'affine_census': census}
        if census['affine_circuits']:
            result['status'] = 'REJECTED_POSITIVE_AFFINE_COMPACTNESS'
        else:
            result['geometry'] = volume_gate(I, h, selected)
            result['status'] = result['geometry']['status']
        cases.append(result)
        out = {'source': cap_source, 'base_source': source, 'source_case': case,
               'all_cap_pair_relations': {f'{i}|{j}': r for (i, j), r in pairs.items()},
               'maximal_compatible_subsets': [list(ids) for ids in maximal], 'cases': cases,
               'scope': 'Complete compatibility graph of the displayed30 caps; geometry '
                        'only for saved cases. No exclusion of the continuous cap family.'}
        (ROOT/(PREFIX+'_vertex_geometry.json')).write_text(json.dumps(out, indent=2)+'\n')
        print({'case': len(cases)-1, 'facets': len(names), 'status': result['status'],
               'numeric_vertices': result.get('geometry', {}).get('numeric_vertices'),
               'numeric_spacelike_vertices': result.get('geometry', {}).get('numeric_spacelike_vertices'),
               'max_norm': result.get('geometry', {}).get('maximum_numeric_vertex_norm')}, flush=True)
        if result['status'] == 'UNRESOLVED_VOLUME':
            break


if __name__ == '__main__':
    run()
