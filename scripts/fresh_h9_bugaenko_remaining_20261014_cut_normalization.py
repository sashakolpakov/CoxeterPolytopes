#!/usr/bin/env python3
"""Literal vertex verification of the nested two-truncation interchange."""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014_cut_normalization'


def cut(vertices, face, new):
    out = set()
    for v in vertices:
        if face <= v:
            out.update(frozenset((v - {i}) | {new}) for i in face)
        else:
            out.add(v)
    return out


def run():
    records = []
    for large in range(3, 10):
        S = frozenset(range(large))
        for small in range(2, large):
            T = frozenset(range(small))
            U = (S - T) | {13}
            counts = [0, 0, 0]
            for vertex in combinations(range(12), 9):
                v = frozenset(vertex)
                original = cut(cut({v}, S, 12), T, 13)
                interchanged = cut(cut({v}, T, 13), U, 12)
                assert original == interchanged
                kind = 0 if not T <= v else (1 if not S <= v else 2)
                counts[kind] += 1
                if kind == 2:
                    expected = {frozenset((v - {t}) | {12}) for t in T}
                    expected |= {frozenset((v - {s, t}) | {12, 13}) for s in S - T for t in T}
                    assert original == expected
            records.append(dict(first_codimension=large, second_codimension=small,
                                interchanged_second_codimension=len(U), vertex_cases=counts,
                                total_vertices_checked=sum(counts)))
    assert len(records) == 28 and sum(r['total_vertices_checked'] for r in records) == 6160
    path = str(Path(__file__).resolve().relative_to(ROOT))
    return dict(status='PASS_NESTED_TRUNCATION_VERTEX_IDENTITY',
        sources={path: sha256((ROOT / path).read_bytes()).hexdigest()}, records=records,
        convention='Old facets0..11, first S-cut new facet12=a, second T-cut new facet13=b. '
                   'Equivalent order: T introduces b, then {b} union (S minus T) introduces a.',
        proof='If T is not contained in S, the first missing face S survives. If T is a '
              'proper subset of S, the displayed vertex identity holds. Since S was an old '
              'face, T and {b} union (S minus T) are both minimal nonfaces of the final '
              'complex: no proper subset can contain an old missing face, and the new '
              'facet generators cannot shorten either support. Thus a final compact '
              'Coxeter nerve forces their orders at most5. A second cut of order>5 is '
              'itself a surviving large missing face. All admissible two-cut recipes '
              'therefore normalize to two cuts of codimensions2..5.',
        scope='Every simple H9 polytope on12 facets and two proper face truncations of '
              'codimensions2..9. The finite check confirms the three universal vertex '
              'cases for every cardinality; it does not enumerate arbitrary simple '
              'polytopes. Together with the written minimal-nonface argument this '
              'extends the complete Bugaenko wedge/two-small-cuts exclusion to all '
              'two-face truncations after that wedge.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    a = ap.parse_args()
    result = run()
    target = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert result == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], 6160)
