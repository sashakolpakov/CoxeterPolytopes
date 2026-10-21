#!/usr/bin/env python3
"""Cover-driven selection with exact affine rejection during every addition."""
import argparse
from itertools import combinations
import json
import numpy as np
import mpmath as mp
from scripts.fresh_h9_replaced_core_20261013_base import *
from scripts.fresh_h9_replaced_core_20261013_cores import affine_census
from scripts.fresh_h9_replaced_core_20261013_geometry import volume_gate


def run(source, output, limit):
    G, I, h, old = profiles()
    data = read(source)
    core = read(data['core_source'])
    add_required(old, core, data.get('required_caps', []))
    rows = [cap_profile(r, f'cap{i}') for i, r in enumerate(data['caps'])]
    pool = angle_pool(ANGLE_SOURCE)
    mp.mp.dps = 90
    cc = np.array([[float(value(r['field'], t)) for t in r['c']] for r in rows])
    pp = np.array([[float(value(r['field'], t)) for t in r['p']] for r in rows])
    products = pp@cc.T
    allowed = np.array([np.cos(np.pi/m)**2 for m in (2, 3, 4, 5, 6, 10)])
    graph = [set() for _ in rows]
    exact_pairs = {}
    for i, j in combinations(range(len(rows)), 2):
        p = products[i, j]
        if p > 1e-7 or not (p*p/4 > 1+1e-7 or min(abs(p*p/4-allowed)) < 1e-7):
            continue
        rel = relation(rows[i], rows[j], pool)
        if rel['kind'] in ('finite', 'dotted'):
            graph[i].add(j)
            graph[j].add(i)
            exact_pairs['|'.join(sorted((rows[i]['name'], rows[j]['name'])))] = rel
    cases, attempts = [], []
    for core_index in (1, 0):
        names = core['cases'][core_index]['walls']
        base_pairs = {'|'.join(sorted((a, b))): relation(old[a], old[b], pool)
                      for a, b in combinations(names, 2)}
        eligible, cap_pairs = [], {}
        for i, row in enumerate(rows):
            rs = {'|'.join(sorted((row['name'], n))): relation(row, old[n], pool) for n in names}
            if all(r['kind'] in ('finite', 'dotted') for r in rs.values()):
                eligible.append(i)
                cap_pairs.update(rs)
        allpairs = base_pairs | cap_pairs | exact_pairs
        census_cache = {}
        def selected_census(chosen):
            key = tuple(sorted(chosen))
            if key not in census_cache:
                wanted = set(names+[rows[i]['name'] for i in chosen])
                selected_pairs = {k: r for k, r in allpairs.items()
                                  if all(n in wanted for n in k.split('|'))}
                census_cache[key] = affine_census(sorted(wanted), selected_pairs)
            return census_cache[key]
        seeds = sorted(eligible, key=lambda i: (-len(data['caps'][i]['cuts']), -len(graph[i]), i))
        seen = set()
        for seed in seeds:
            chosen, available, covered = [], set(eligible), set()
            available.discard(seed)
            ordered = [seed]
            while ordered:
                proposed = ordered[0]
                trial = sorted(chosen+[proposed])
                allnames = names+[rows[i]['name'] for i in trial]
                census = selected_census(trial)
                if census['affine_circuits']:
                    attempts.append({'core_index': core_index, 'caps': trial,
                                     'affine_circuits': census['affine_circuits']})
                else:
                    chosen = trial
                    covered.update(data['caps'][proposed]['cuts'])
                    available &= graph[proposed]
                available.discard(proposed)
                ordered = sorted(available, key=lambda i: (
                    -len(set(data['caps'][i]['cuts'])-covered),
                    -len(graph[i]&available), i))
            if not chosen or tuple(chosen) in seen:
                continue
            seen.add(tuple(chosen))
            allnames = names+[rows[i]['name'] for i in chosen]
            census = selected_census(chosen)
            assert not census['affine_circuits']
            selected = [old[n] for n in names]+[rows[i] for i in chosen]
            geo = volume_gate(I, h, selected)
            cases.append({'core_index': core_index, 'selected_caps': chosen,
                          'walls': allnames, 'facet_count': len(allnames),
                          'pair_checks': {k: v for k, v in allpairs.items()
                                          if all(n in allnames for n in k.split('|'))},
                          'affine_census': census, 'status': geo['status'], 'geometry': geo})
            result = {'source': source, 'core_source': data['core_source'],
                      'required_caps': data.get('required_caps', []), 'cases': cases,
                      'affine_rejected_prefixes': attempts,
                      'scope': 'Cover-driven compatible additions, with exact affine-circuit '
                               'rejection at each prefix. Only saved exact escapes decide volume; '
                               'not an exhaustive selection search.'}
            (ROOT/output).write_text(json.dumps(result, indent=2)+'\n')
            print({'case': len(cases)-1, 'core': core_index, 'selected': chosen,
                   'walls': len(allnames), 'status': geo['status'],
                   'maximum_numeric_vertex_norm': geo['maximum_numeric_vertex_norm']}, flush=True)
            if len(cases) >= limit or geo['status'] == 'UNRESOLVED_VOLUME':
                return


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--limit', type=int, default=6)
    a = p.parse_args()
    run(a.source, a.output, a.limit)
