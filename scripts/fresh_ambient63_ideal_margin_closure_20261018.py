#!/usr/bin/env python3
"""Independent ideal rank62 arithmetic for the specified improved envelopes.

Geometric validity of the new local weights and the finite budget are separate
proof obligations. This script imports no project arithmetic implementation.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_ambient63_ideal_margin_closure_20261018'
FRONT = 'fresh_ambient63_payment_frontier_20261018'
PREVIOUS = 'fresh_ambient64_ideal_margin_closure_20261018'
CONVEX = 'fresh_distance14_unrestricted_dimension156_20261017'
STEMS = (FRONT, PREVIOUS)
RANK = 62


def profiles():
    def visit(rank, left, middle):
        if rank == 18:
            for low in range(9):
                high = left-low
                if high == 0 or high >= 18:
                    yield tuple(([low] if low else [])+middle+([high] if high else []))
            return
        for count in range(left//rank+1):
            yield from visit(rank+1, left-count*rank, middle+[rank]*count)
    rows = sorted(visit(9, RANK, []))
    assert len(rows) == len(set(rows)) == 2204
    assert all(sum(rs) == RANK and list(rs) == sorted(rs) for rs in rows)
    return rows


def envelope(r):
    if r <= 9:
        return F(r, 2)
    return {10: F(54, 11), 11: F(21, 4), 12: F(147, 26), 13: F(6),
            14: F(94, 15), 15: F(107, 16), 16: F(118, 17)}.get(r, F(163, 23))


def scalar(rs):
    h = [1]
    for rank in rs:
        out, window = [], 0
        for j in range(len(h)+rank):
            if j < len(h):
                window += h[j]
            if 0 <= j-rank-1 < len(h):
                window -= h[j-rank-1]
            out.append(window)
        h = out
    V = prod(r+1 for r in rs)
    assert len(h) == RANK+1 and h == h[::-1] and sum(h) == V
    G = sum(h[j]*min(j, RANK-j) for j in range(RANK+1))
    curvature = F(5*G-RANK*V, 2)
    cost = sum(map(envelope, rs), F())
    margin = curvature/V-cost
    assert margin > 0
    return dict(ranks=list(rs), h_polynomial=h, link_vertices=V, G=G,
                normalized_cost=str(cost), curvature=str(curvature),
                normalized_margin=str(margin), absolute_margin=str(V*margin))


def run():
    docs = {stem: json.loads((ROOT/f'results/{stem}.json').read_text()) for stem in STEMS}
    for stem, data in docs.items():
        assert data['status'].startswith('PASS'), stem
        for path, h in data['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == h, (stem, path)
    rows = [scalar(rs) for rs in profiles()]
    front = docs[FRONT]
    baseline = {tuple(r['ranks']): r for r in front['ideal_current_rows']}
    conditional = next(r for r in front['ideal_conditional_scenarios'] if
                       r['total_cost_drops_by_affine_rank'] == {'10': '1/2', '11': '1', '12': '1/2'})
    expected = {tuple(r['ranks']): r for r in conditional['rows']}
    assert set(baseline) == set(expected) == {tuple(r['ranks']) for r in rows}
    for row in rows:
        key = tuple(row['ranks']); old = baseline[key]
        for field in ('h_polynomial', 'link_vertices', 'G', 'curvature'):
            assert row[field] == old[field], (key, field)
        assert row['normalized_margin'] == expected[key]['normalized_margin']
        gain = F(key.count(10), 22)+F(key.count(11), 12)+F(key.count(12), 26)
        assert F(old['normalized_cost'])-F(row['normalized_cost']) == gain
    minimum = min(rows, key=lambda r: (F(r['normalized_margin']), r['ranks']))
    absolute = min(rows, key=lambda r: (F(r['absolute_margin']), r['ranks']))
    assert minimum['normalized_margin'] == '3505/201344' and minimum['ranks'] == [10, 10, 12, 15, 15]
    assert absolute['absolute_margin'] == '139/46' and absolute['ranks'] == [62]
    assert sum(F(r['normalized_margin']) == F(minimum['normalized_margin']) for r in rows) == 1
    assert sum(F(r['absolute_margin']) == F(absolute['absolute_margin']) for r in rows) == 1
    paths = [Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS
                             for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT/f'results/{CONVEX}.md', ROOT/'results/fresh_cusp_link_shape_curvature_20261017.md']
    return dict(status='PASS_AMBIENT63_COMPLETE_IDEAL_MARGIN_ARITHMETIC',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
                ambient_dimension=63, cusp_link_rank=62, complete_terminal_profiles=2204,
                low_concentration_threshold=9, high_concentration_threshold=17,
                envelopes=[dict(rank=r, cost=str(envelope(r))) for r in range(1, 63)],
                required_new_total_cost_envelopes={'11': '54', '12': '63', '13': '147/2'},
                rows=rows, minimum_normalized_margin=minimum, minimum_absolute_margin=absolute,
                all_profiles_positive=True, global_additive_constant='4960',
                finite_budget_proved_here=False, local_geometric_coverage_proved_here=False,
                new_dimension_bound=None,
                scope='Independent multiplicity enumeration and sliding-window arithmetic match every scalar of the earlier conditional2204-profile target. Written concentration gives absolute surplus139/46 for all actual cusp profiles if the specified local envelopes hold. Local geometry, enlarged finite coupon budget and global composition remain separate.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--verify', action='store_true'); args = parser.parse_args()
    out = json.loads(json.dumps(run())); path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        assert not path.exists(), 'Existing artifact; use --verify.'
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], 'profiles', len(out['rows']), 'normalized', out['minimum_normalized_margin']['normalized_margin'],
          'absolute', out['minimum_absolute_margin']['absolute_margin'], flush=True)
