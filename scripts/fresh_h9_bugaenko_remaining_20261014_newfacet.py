#!/usr/bin/env python3
"""Cut H8, wedge its NEW facet, then cut: a changed 14-facet H9 family."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import time
from fresh_h9_bugaenko_surgery_20261014_incidence import (
    BASE_BLOCKS, PAIRS, masks, tup, minimal, canonical, first_faces, remap_counts)

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014_newfacet'
RAW = 'data/known_examples/bugaenko_compact_h8.json'
OLD = 'results/fresh_h9_bugaenko_surgery_20261014_incidence.json'
GATE = 'results/fresh_h9_bugaenko_remaining_20261014_block_gate.json'


def truncate(M, V, F, new, dimension):
    assert all(m & F != m for m in M)
    out = minimal([F, *[m for m in M if m & F != F],
                   *[(1 << new) | (m & ~F) for m in M]])
    vertices = set()
    for v in V:
        if v & F == F:
            vertices.update((v & ~(1 << i)) | (1 << new) for i in tup(F))
        else:
            vertices.add(v)
    assert all(v.bit_count() == dimension for v in vertices)
    return out, tuple(sorted(vertices))


def wedge(M, V, old=11, new=12):
    missing = tuple(sorted(m | (1 << new) if m >> old & 1 else m for m in M))
    vertices = set()
    for v in V:
        if v >> old & 1:
            vertices.add(v | (1 << new))
        else:
            vertices.update([v | (1 << old), v | (1 << new)])
    assert all(v.bit_count() == 9 for v in vertices)
    return missing, tuple(sorted(vertices))


def second_faces(F, H, M):
    parts = [(tuple(i for i in b if F >> i & 1), tuple(i for i in b if not F >> i & 1))
             for b in BASE_BLOCKS]
    choices = [list(product(range(len(a) + 1), range(len(b) + 1))) for a, b in parts]
    for counts in product(*choices):
        if counts != min(remap_counts(counts, p) for p in H):
            continue
        old_count = sum(sum(c) for c in counts)
        for new_count in range(3):
            if not 2 <= old_count + new_count <= 5:
                continue
            face = sum(1 << i for i in [11, 12][:new_count])
            for (a, b), (ca, cb) in zip(parts, counts):
                face |= sum(1 << i for i in a[:ca] + b[:cb])
            if all(m & face != m for m in M):
                yield face


def block_obstruction(M, V):
    allowed = {(3, 4), (3, 5), (4, 4), (4, 5), (5, 5)}
    for i, a in enumerate(M):
        for j, b in enumerate(M[i + 1:], i + 1):
            if a & b or tuple(sorted([a.bit_count(), b.bit_count()])) not in allowed:
                continue
            union = a | b
            if any(k not in [i, j] and c & union == c for k, c in enumerate(M)):
                continue
            assert all(any(v & (union ^ (1 << x) ^ (1 << y)) == union ^ (1 << x) ^ (1 << y)
                           for v in V) for x in tup(a) for y in tup(b))
            return [i, j]
    return None


def run():
    old = json.loads((ROOT / OLD).read_text())
    gate = json.loads((ROOT / GATE).read_text())
    for d in [old, gate]:
        for p, h in d['sources'].items():
            assert sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    assert gate['status'] == 'PASS_EXACT_TWO_LANNER_BLOCK_OBSTRUCTION'
    old_lookup = {tuple(r['canonical_minimal_nonfaces']): r['case'] for r in old['cases']}
    base = json.loads((ROOT / RAW).read_text())
    V = masks(base['geometry']['vertices'])
    M = masks([BASE_BLOCKS[a] + BASE_BLOCKS[b] for a, b in PAIRS])
    unique = {}
    rows = []
    stats = Counter()
    start = time.monotonic()
    for first_index, (F, counts, H) in enumerate(first_faces(BASE_BLOCKS, M)):
        A, B = truncate(M, V, F, 11, 8)
        A, B = wedge(A, B)
        stats['first_face_orbits'] += 1
        for S in second_faces(F, H, A):
            stats['ordered_recipe_orbits'] += 1
            C, W = truncate(A, B, S, 13, 9)
            if max(map(int.bit_count, C)) > 5:
                stats['large_minimal_nonface'] += 1
                continue
            stats['passes_missing_order'] += 1
            canon, order = canonical(C, 14)
            move = {a: b for b, a in enumerate(order)}
            CW = tuple(sorted(sum(1 << move[i] for i in tup(v)) for v in W))
            if canon in unique:
                r = rows[unique[canon]]
                assert tuple(r['canonical_vertices']) == CW
                r['recipe_multiplicity'] += 1
                continue
            index = len(rows)
            unique[canon] = index
            bad = block_obstruction(canon, CW)
            old_case = old_lookup.get(canon)
            r = dict(case=index, first_face=list(tup(F)), second_face=list(tup(S)),
                     first_face_orbit=first_index, wedge_carrier=11,
                     canonical_order=order, canonical_minimal_nonfaces=canon,
                     canonical_vertices=CW, recipe_multiplicity=1,
                     old_catalogue_case=old_case, block_obstruction=bad)
            rows.append(r)
        if first_index % 10 == 0:
            print('first', first_index, 'recipes', stats['ordered_recipe_orbits'],
                  'types', len(rows), 'seconds', round(time.monotonic() - start, 2), flush=True)
    candidates = [r['case'] for r in rows if r['old_catalogue_case'] is None and r['block_obstruction'] is None]
    for r in rows:
        counts = Counter(v & ~(1 << i) for v in r['canonical_vertices'] for i in tup(v))
        assert set(counts.values()) == {2}
    files = [RAW, OLD, GATE, str(Path(__file__).resolve().relative_to(ROOT)),
             'scripts/fresh_h9_bugaenko_surgery_20261014_incidence.py']
    return dict(status='EXACT_NEW_TRUNCATION_FACET_WEDGE_INCIDENCE',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
        counts=dict(stats), distinct_types=len(rows),
        old_catalogue_types=sum(r['old_catalogue_case'] is not None for r in rows),
        new_types_with_block_obstruction=sum(r['old_catalogue_case'] is None and r['block_obstruction'] is not None for r in rows),
        genuinely_new_unexcluded_types=candidates, cases=rows,
        scope='Original compact H8 is cut at every codimension2..5 face orbit; the newly '
              'created facet11 is wedged into facets11/12; every second face of '
              'codimension2..5 is truncated. Within-block membership stabilizers and '
              'the new twin-facet swap give safe orbit reduction. Full nerves are '
              'deduplicated against the closed old catalogue and screened by the exact '
              'two-Lanner-block gate. Surviving incidence is explicitly polytopal but '
              'has no asserted Coxeter metric.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    a = ap.parse_args()
    result = json.loads(json.dumps(run()))
    target = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert result == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['counts'], 'types', result['distinct_types'],
          'new unexcluded', len(result['genuinely_new_unexcluded_types']), flush=True)
