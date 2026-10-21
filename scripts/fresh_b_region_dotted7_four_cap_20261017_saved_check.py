#!/usr/bin/env python3
"""Read-only saved-record checks; no producer or numerical geometry imports."""
import json
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_b_region_dotted7_four_cap_20261017'


def main():
    d = json.loads((ROOT / 'results' / (TAG + '.json')).read_text())
    for p, h in d['sources'].items():
        assert sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    primary = json.loads((ROOT / 'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vertices = [tuple(v) for v in primary['retained_old_finite_supports']]
    vertices += [tuple(v['support']) for v in primary['new_finite_vertices']]
    assert len(vertices) == len(set(vertices)) == 199
    edge_ends = defaultdict(set)
    face_vertices = defaultdict(set)
    for i, v in enumerate(vertices):
        for e in combinations(v, 8):
            edge_ends[e].add(i)
        for k in range(8):
            for f in combinations(v, k):
                face_vertices[f].add(i)
    compact = {}
    for f, inds in face_vertices.items():
        if len(inds) < 10 - len(f):
            continue
        if all(len(edge_ends[tuple(j for j in vertices[i] if j != k)]) == 2
               for i in inds for k in vertices[i] if k not in f):
            compact[f] = inds
    for row in d['exact_deletion_comparisons']:
        lost = {i for i, v in enumerate(vertices) if set(v) & set(row['deleted'])}
        assert sorted(vertices[i] for i in lost) == sorted(map(tuple, row['lost_finite_supports']))
        assert len(lost) == row['lost_finite_count']
        assert sum(8 in vertices[i] for i in lost) == row['lost_supports_containing8']
        touched = {f: inds for f, inds in compact.items() if inds & lost}
        assert set(touched) == {tuple(x['support']) for x in row['compact_faces_touched']}
        assert dict(Counter(str(9-len(f)) for f in touched)) == row['compact_faces_by_dimension']
    disk = d['disk']; faces = disk['triangles']; boundary = disk['boundary']
    es = Counter(tuple(sorted(e)) for f in faces for e in combinations(f, 2))
    assert (len(set(sum(faces, []))), len(es), len(faces)) == (8, 17, 10)
    assert {e for e, n in es.items() if n == 1} == {
        tuple(sorted(e)) for e in zip(boundary, boundary[1:] + boundary[:1])}
    assert all(n in (1, 2) for n in es.values())
    for link in disk['links']:
        v = link['vertex']; expected = [tuple(j for j in f if j != v) for f in faces if v in f]
        assert sorted(expected) == sorted(map(tuple, link['edges']))
        degrees = Counter(j for e in expected for j in e)
        seen = {next(iter(degrees))}
        while True:
            new = seen | {j for e in expected if set(e) & seen for j in e}
            if new == seen:
                break
            seen = new
        assert seen == set(degrees)
        assert sorted(degrees.values()) == ([1, 1] + [2] * (len(degrees)-2) if v in boundary else [2] * len(degrees))
    ledgers = d['boundary_ledgers']
    for l in ledgers:
        assert not l['warnings']
        assert dict(Counter(r['kind'] for r in l['rays'])) == l['counts']
        finite = {tuple(sorted(r['active'])) for r in l['rays'] if r['kind'] == 'finite'}
        assert finite == set(map(tuple, l['PD_finite_supports']))
        assert l['open_ends'] == l['PD_open_ends']
        assert l['open_ends'] == sum(1 if e['kind'] == 'one_open' else 2 if e['kind'] == 'two_open' else 0 for e in l['edges'])
    assert len(d['samples']) == 124 and len(d['brackets']) == len(d['events']) == 168
    assert len(d['completions']) == 258
    assert not d['improvements'] and not d['unresolved'] and not d['bank_unresolved']
    assert Counter(e['L_profile'] for e in d['events']) == {0: 142, 1: 16, 4: 10}
    products = []; patterns = Counter()
    for row in d['completions']:
        cap = row['cap']; l = ledgers[cap['score']]
        assert cap['profile'] == 0 and cap['labels'][:2] == [2, 2]
        assert cap['center_product'] < 0 and cap['norm_residual'] < 2e-5
        products.append(dict(cap['native_products'])[69])
        assert products[-1] < -2
        assert [i for i, c in enumerate(row['finite_corners']) if c is None] == [6, 9]
        assert l['counts'] == {'finite': 32, 'exterior': 6} and l['open_ends'] == 10
        assert not row['completed_disk']
        for rec in row['whole_cap_faces']:
            node = rec['node']; got = [r for r in l['rays'] if node in r['active']]
            expected = {tuple(sorted(f)) for f in faces if node in f}
            actual = all(r['kind'] == 'finite' for r in got) and {tuple(sorted(r['active'])) for r in got} == expected
            assert actual == rec['finite_exactly_expected'] == (node == 67)
        patterns[tuple(tuple(r['active']) for r in l['rays'] if r['kind'] == 'exterior')] += 1
    assert sorted(patterns.values()) == [84, 174]
    print('PASS saved bindings, 199-vertex deletion ledger, disk, 8 boundary ledgers, 168 events, 258 completions')
    print('RQ saved product range', min(products), max(products), 'exterior pattern counts', sorted(patterns.values()))


if __name__ == '__main__':
    main()
