#!/usr/bin/env python3
"""Exact compactness gates and volume witnesses for coherent gap completions."""
import argparse
from itertools import combinations
import json
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h9_replaced_core_20261013_cores import affine_census
from scripts.fresh_h9_replaced_core_20261013_geometry import volume_gate


def run(source, output, maximum=12):
    data = read(source)
    G, I, h, original, base_source, base_case = configuration(data.get('configuration', '45'))
    rows = [cap_profile(r, f'{data.get("cap_name_prefix", "gap")}{i}') for i, r in enumerate(data['caps'])]
    pool = angle_pool(ANGLE_SOURCE)
    base_names = list(original)
    base_pairs = read(base_source)['cases'][base_case]['pair_checks']
    cap_pairs = {'|'.join(sorted((row['name'], n))): rel
                 for i, row in enumerate(rows) for n, rel in data['caps'][i]['pair_checks'].items()}
    pair_graph = {}
    for i, j in combinations(range(len(rows)), 2):
        pair_graph[i, j] = relation(rows[i], rows[j], pool)
    proposals = [(i,) for i in range(len(rows))]
    proposals += [pair for pair, rel in pair_graph.items() if rel['kind'] in ('finite', 'dotted')]
    # Every subset of this tiny graph is considered up to cardinality three.
    proposals += [ids for ids in combinations(range(len(rows)), 3)
                  if all(pair_graph[i, j]['kind'] in ('finite', 'dotted') for i, j in combinations(ids, 2))]
    proposals.sort(key=lambda ids: (-len(set().union(*(set(data['caps'][i]['cuts']) for i in ids))), -len(ids), ids))
    cases = []
    for ids in proposals:
        dropped = sorted(set().union(*(set(data['caps'][i]['dropped']) for i in ids)))
        if len(dropped) > data['drop_budget']:
            continue
        names = [n for n in base_names if n not in dropped]+[rows[i]['name'] for i in ids]
        selected = [original[n] for n in base_names if n not in dropped]+[rows[i] for i in ids]
        pairs = dict(base_pairs)|cap_pairs
        for i, j in combinations(ids, 2):
            pairs['|'.join(sorted((rows[i]['name'], rows[j]['name'])))] = pair_graph[i, j]
        pairs = {k: r for k, r in pairs.items() if all(n in names for n in k.split('|'))}
        assert len(pairs) == len(names)*(len(names)-1)//2
        assert all(r['kind'] in ('finite', 'dotted') for r in pairs.values())
        census = affine_census(names, pairs)
        case = {'selected_caps': list(ids), 'dropped': dropped, 'walls': names,
                'facet_count': len(names), 'pair_checks': pairs, 'affine_census': census}
        if census['affine_circuits']:
            case['status'] = 'REJECTED_POSITIVE_AFFINE_COMPACTNESS'
        else:
            case['geometry'] = volume_gate(I, h, selected)
            case['status'] = case['geometry']['status']
        cases.append(case)
        result = {'source': source, 'base_source': base_source, 'source_case': base_case,
                  'candidate_pair_relations': {f'{i}|{j}': r for (i, j), r in pair_graph.items()},
                  'cases': cases, 'scope': 'Exact pair and affine gates for selected coherent '
                  'extensions/replacements. Saved exact strict escapes prove infinite volume; '
                  'absence of a found escape alone has no positive implication.'}
        (ROOT/output).write_text(json.dumps(result, indent=2)+'\n')
        print({'case': len(cases)-1, 'selected': ids, 'dropped': dropped,
               'facets': len(names), 'affine': len(census['affine_circuits']),
               'status': case['status'],
               'maximum_numeric_vertex_norm': case.get('geometry', {}).get('maximum_numeric_vertex_norm')}, flush=True)
        if len(cases) >= maximum or case['status'] == 'UNRESOLVED_VOLUME':
            break


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', default=PREFIX+'_caps_drop0.json')
    p.add_argument('--output', default=PREFIX+'_geometry.json')
    p.add_argument('--maximum', type=int, default=12)
    a = p.parse_args()
    run(a.source, a.output, a.maximum)
