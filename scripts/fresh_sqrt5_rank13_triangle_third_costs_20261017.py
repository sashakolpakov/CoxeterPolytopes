#!/usr/bin/env python3
"""Literal third vertices and sharper triangle charges for the joint ADE set."""
import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

import fresh_sqrt5_forest_sensitive_weights_20261015_review as C
import fresh_sqrt5_rank13_high_ade_endpoints_20261017 as E

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_rank13_triangle_third_costs_20261017'
TRI = 'fresh_sqrt5_rank13_high_ade_triangles_20261017'


def read(tag):
    data = json.loads((ROOT / f'results/{tag}.json').read_text())
    for path, digest in data['sources'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest, path
    return data


def name(g):
    if all(m == 3 for row in g.values() for m in row.values()):
        out = E.family({i: set(row) for i, row in g.items()})
        assert out is not None
        return out
    return C.family(g)


def run():
    endpoints, triangles = read(E.TAG), read(TRI)
    assert len(endpoints['targets']) == len(triangles['targets']) == 215
    third_types, patterns, outputs, counts = {}, {}, [], Counter()
    for target, tri in zip(endpoints['targets'], triangles['targets']):
        assert target['index'] == tri['index'] and target['forest'] == tri['forest']
        old_edges = {tuple(edge): 3 for edge in target['edges']}
        reduced, h4, old_pairs = set(), set(), set()
        for occurrence_index, occurrence in enumerate(tri['accepted_occurrences']):
            i, j, z = occurrence['port'], occurrence['other'], 13
            attachment = target['ports'][i]['attachments'][occurrence['attachment']]
            assert attachment['destination'] == occurrence['destination']
            edges = dict(old_edges)
            edges.update({tuple(sorted((k, z))): 3 for k in attachment['contacts']})
            edges[(i, z)] = occurrence['label']
            vertices = []
            for omitted, marks in ((z, (i, j)), (i, (j, z)), (j, (i, z))):
                graph = C.graph([k for k in range(14) if k != omitted], edges)
                parts = C.pieces(graph)
                forest = sorted(name(g) for g in parts)
                marked = next((g for g in parts if all(k in g for k in marks)), None)
                corner = F(0) if marked is None else C.base(marked, *marks)
                total = sum((C.base(g, a, b) for g in parts for a, b in combinations(g, 2)), F(0))
                vertices.append(dict(forest=forest, marked_forest_code=sorted(C.code(g, marks) for g in parts),
                                     corner_baseline=str(corner), vertex_baseline=str(total)))
            assert vertices[0]['forest'] == target['forest']
            assert vertices[1]['forest'] == endpoints['targets'][occurrence['destination']]['forest']
            third = vertices[2]
            gthird = C.graph([k for k in range(14) if k != j], edges)
            exceptional = next(g for g in C.pieces(gthird) if i in g)
            assert z in exceptional and name(exceptional) == occurrence['third_exceptional_component']
            assert exceptional[i][z] == occurrence['label']
            w = F(third['corner_baseline'])
            assert w == (2 if name(exceptional) == 'H4' else 1)
            assert sum(F(v['corner_baseline']) for v in vertices) >= 2
            pair = tuple(sorted((i, j)))
            old_pairs.add(pair)
            (h4 if w == 2 else reduced).add(pair)
            counts['occurrences'] += 1
            counts['third_' + name(exceptional)] += 1
            key = '|'.join(third['forest'])
            record = third_types.setdefault(key, dict(forest=third['forest'], baseline=third['vertex_baseline'],
                                                     marked_components=Counter(), occurrences=0))
            assert record['baseline'] == third['vertex_baseline']
            record['marked_components'][name(exceptional)] += 1
            record['occurrences'] += 1
            # Sort the two selected corners only; the third corner is distinguished.
            selected = sorted(vertices[:2], key=lambda v: json.dumps(v['marked_forest_code']))
            pkey = json.dumps(selected + [third], sort_keys=True)
            pat = patterns.setdefault(pkey, dict(selected_corners=selected, third_corner=third,
                                                required_selected_charge=str(2-w), occurrences=0,
                                                witness=dict(target=target['index'], occurrence=occurrence_index)))
            pat['occurrences'] += 1
        assert old_pairs == {tuple(p) for p in tri['possible_two_special_triangle_corners']}
        outputs.append(dict(index=target['index'], forest=target['forest'], baseline=target['baseline'],
                            old_possible_triangle_corners=len(old_pairs),
                            half_charge_corners=[list(p) for p in sorted(reduced)],
                            half_charge_corner_count=len(reduced),
                            h4_only_zero_charge_corners=[list(p) for p in sorted(h4-reduced)],
                            triangle_charge_upper=str(F(len(reduced), 2))))
    assert counts['occurrences'] == 24523
    lookup = {tuple(x['forest']): x for x in outputs}
    for forest, expected in [(('D13',), 5), (('A13',), 6), (('D4', 'D4', 'D5'), 11),
                             (('A1', 'A3', 'A3', 'A3', 'A3'), 24)]:
        assert lookup[forest]['half_charge_corner_count'] == expected
    paths = [f'scripts/{TAG}.py', f'results/{E.TAG}.json', f'results/{TRI}.json',
             f'scripts/{E.TAG}.py', f'scripts/{TRI}.py',
             'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py']
    return dict(status='PASS_COMPLETE_JOINT_ADE_TRIANGLE_THIRD_COSTS',
                sources={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                dimension=13, counts=dict(counts), third_forest_count=len(third_types),
                distinct_marked_triangle_constraints=len(patterns),
                third_forests=sorted(third_types.values(), key=lambda x: x['forest']),
                triangle_patterns=sorted(patterns.values(), key=lambda x: json.dumps(x, sort_keys=True)),
                targets=outputs,
                proof='The third marked pair is the framing bond. Its baseline is2 exactly in H4, and1 in I2(5),H3,I2(10). With third charge2 the selected corners need0 each; otherwise one half each suffices. Count a root pair at the selected vertex only once, even if many local occurrences use it.',
                accounting='This is not compatible with keeping every other corner payment at1: H4 framing5 corners must receive2. Under a generic nonspecial rule of1 at every root pair and one extra at each H4 framing5 bond, a rank13 nonspecial vertex has cost at most binom(13,2)+floor(13/4)=81. More refined actual third-forest costs and marked constraints are saved separately.',
                quad_composition='If J is a complete set of possible all-selected quadrilateral pairs and I is the saved half-charge set, special cost is at most |I|/2+|J minus I|/4. Other corners can be paid as in accounting. Quad completeness is not established by this package.',
                scope='Exact reconstruction and compression of the inherited complete two-selected-triangle census. No joint quadrilateral classification, dimension bound, or simultaneous realization is claimed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], data['third_forest_count'], data['distinct_marked_triangle_constraints'], flush=True)
