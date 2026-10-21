#!/usr/bin/env python3
"""Exact marked quotient of the new adjacent-discount first extensions.
This reduces literal queue seeds only; it performs no new arithmetic filtering.
"""
import argparse
import json
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path

import fresh_sqrt5_sealed_subcore_closure_20261015_review as I

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_discount_attachment_seeds_20261016'
INPUT = 'fresh_sqrt5_adjacent_discount_frontier_20261016'


def generate():
    source = json.loads((ROOT/f'results/{INPUT}.json').read_text())
    assert source['status'] == 'PASS_COMPLETE_MINIMAL_UNION_ADJACENT_DISCOUNT_SCREEN'
    records, mappings, buckets = [], [], defaultdict(list)
    attempts = 0
    for index, literal in enumerate(source['outgoing_roots']):
        row = literal['row']
        matrix = I.matrix(row)
        key = I.invariant(matrix)
        found = None
        for j, other in buckets[key]:
            attempts += 1
            mp = I.isomorphism(matrix, other)
            if mp is not None:
                I.check_map(matrix, other, mp)
                found = (j, mp)
                break
        if found is None:
            j = len(records)
            buckets[key].append((j, matrix))
            records.append(dict(index=j, row=row, parameter_type=literal['parameter_type'],
                                analysis=literal['analysis'], literal_origins=[]))
            found = (j, list(range(row['nodes'])))
        j, mp = found
        assert records[j]['parameter_type'] == literal['parameter_type']
        I.check_map(matrix, I.matrix(records[j]['row']), mp)
        records[j]['literal_origins'].append(index)
        mappings.append(dict(literal_index=index, parent_class=literal['parent_class'],
                             contacts=literal['contacts'], target=j, map=mp))
    assert len(mappings) == 6848
    assert sorted(i for r in records for i in r['literal_origins']) == list(range(6848))
    zero_touch = [dict(class_index=r['class_index'], row=r['row'], analysis=r['analysis'],
                       sealed=r['inherited_status']=='SEALED_SINGULAR_FINITE_CORE')
                  for r in source['records']
                  if r['singular'] and r['analysis']['affected_if_components_entire']]
    assert len(zero_touch) == 37
    paths = [f'scripts/{TAG}.py', f'scripts/{I.__name__}.py',
             f'scripts/{INPUT}.py', f'results/{INPUT}.json',
             'results/fresh_sqrt5_all_quad_contacts_20261015_review.json']
    return dict(status='PASS_EXACT_DISCOUNT_FIRST_EXTENSION_MARKED_QUOTIENT',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                literal_roots=len(mappings), marked_roots=len(records),
                marked_counts=dict(Counter(r['parameter_type'] for r in records)),
                nodes_counts=dict(Counter(r['row']['nodes'] for r in records)),
                isomorphism_attempts=attempts, literal_maps=mappings, records=records,
                immediate_zero_touch_obligations=zero_touch,
                scope='Exact opposite-pair-preserving marked isomorphism quotient of all6848 still-affected first extensions. Dotted0 is retained symbolically; no parameter is solved or identified with another distance. All37 affected minimal singular zero-touch states are separately exported. No complete attachment closure or payment claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(generate()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
        print('PASS saved marked quotient replay')
    else:
        path.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], data['literal_roots'], data['marked_roots'], data['marked_counts'])
