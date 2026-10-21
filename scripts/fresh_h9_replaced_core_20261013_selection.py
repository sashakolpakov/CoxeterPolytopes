#!/usr/bin/env python3
"""Coherent additions to changed cores, with early exact affine and volume gates."""
import argparse
from itertools import combinations
import json
import numpy as np
import mpmath as mp
from scripts.fresh_h9_replaced_core_20261013_base import *
from scripts.fresh_h9_replaced_core_20261013_cores import affine_census
from scripts.fresh_h9_replaced_core_20261013_geometry import volume_gate

SOURCE = 'data/enumerations/fresh_h9_replaced_core_20261013_polar_caps.json'
OUT = ROOT/'data/enumerations/fresh_h9_replaced_core_20261013_selection.json'


def run(limit=6):
    G, I, h, original = profiles()
    data, core = read(SOURCE), read(read(SOURCE)['core_source'])
    required_names = add_required(original, core, data.get('required_caps', []))
    rows = [cap_profile(r, f'cap{i}') for i, r in enumerate(data['caps'])]
    pool = angle_pool(ANGLE_SOURCE)
    core_pairs = dict(core['pair_checks'])
    union_names = list(dict.fromkeys(n for case in core['cases'] for n in case['walls']))
    for a, b in combinations(union_names, 2):
        if a in required_names or b in required_names:
            core_pairs['|'.join(sorted((a, b)))] = relation(original[a], original[b], pool)
    mp.mp.dps = 90
    cc = np.array([[float(value(r['field'], t)) for t in r['c']] for r in rows])
    pn = np.array([[float(value(r['field'], t)) for t in r['p']] for r in rows])
    products = pn @ cc.T
    allowed = np.array([np.cos(np.pi/m)**2 for m in (2, 3, 4, 5, 6, 10)])
    graph = [set() for _ in rows]
    for i, j in combinations(range(len(rows)), 2):
        p = products[i, j]
        ratio = p*p/4
        if p <= 1e-7 and (ratio > 1+1e-7 or min(abs(ratio-allowed)) < 1e-7):
            graph[i].add(j)
            graph[j].add(i)
    cases, seen = [], set()
    seeds = list(range(min(len(rows), 16)))
    seeds += sorted(range(len(rows)), key=lambda i: (-len(graph[i]), i))[:16]
    for seed in seeds:
        selected = [seed]
        eligible = set(graph[seed])
        while eligible:
            j = max(eligible, key=lambda j: (len(graph[j] & eligible), -j))
            selected.append(j)
            eligible &= graph[j]
        selected.sort()
        if tuple(selected) in seen:
            continue
        seen.add(tuple(selected))
        for core_index in (0, 1):
            names = list(core['cases'][core_index]['walls'])
            caps = [rows[i] for i in selected]
            pairchecks = dict(core_pairs)
            valid = True
            for row in caps:
                for name in names:
                    r = relation(row, original[name], pool)
                    pairchecks['|'.join(sorted((row['name'], name)))] = r
                    valid &= r['kind'] in ('finite', 'dotted')
            for a, b in combinations(caps, 2):
                r = relation(a, b, pool)
                pairchecks['|'.join(sorted((a['name'], b['name'])))] = r
                valid &= r['kind'] in ('finite', 'dotted')
            if not valid:
                continue
            allnames = names+[r['name'] for r in caps]
            census = affine_census(allnames, pairchecks)
            case = {'core_index': core_index, 'selected_caps': selected,
                    'walls': allnames, 'facet_count': len(allnames),
                    'pair_checks': {k: v for k, v in pairchecks.items()
                                    if all(t in allnames for t in k.split('|'))},
                    'affine_census': census}
            if census['affine_circuits']:
                case['status'] = 'REJECTED_POSITIVE_AFFINE_COMPACTNESS'
            else:
                case['geometry'] = volume_gate(I, h, [original[n] for n in names]+caps)
                case['status'] = case['geometry']['status']
            cases.append(case)
            result = {'source': SOURCE, 'core_source': data['core_source'],
                      'required_caps': data.get('required_caps', []),
                      'numeric_compatibility_edges': sum(map(len, graph))//2,
                      'cases': cases,
                      'scope': 'Each accepted configuration has exact Coxeter pair relations, '
                               'an exact complete connected affine census and an exact '
                               'strict spacelike volume rejection where a witness is saved. '
                               'Numeric graph proposals are not an exhaustive selection search.'}
            OUT.write_text(json.dumps(result, indent=2)+'\n')
            print({'case': len(cases)-1, 'core': core_index, 'caps': selected,
                   'facets': len(allnames), 'affine': len(census['affine_circuits']),
                   'status': case['status'],
                   'maximum_vertex_norm': case.get('geometry', {}).get('maximum_numeric_vertex_norm')}, flush=True)
            if len(cases) >= limit or case['status'] == 'UNRESOLVED_VOLUME':
                return


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=6)
    p.add_argument('--source')
    p.add_argument('--output')
    a = p.parse_args()
    if a.source:
        SOURCE = a.source
    if a.output:
        OUT = ROOT/a.output
    run(a.limit)
