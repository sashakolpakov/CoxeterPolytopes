#!/usr/bin/env python3
"""Read-only independent scalar/domain check; no face-weight or geometry proof."""
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact29_savings_frontier_20261018'
d = json.loads((ROOT / 'results' / (TAG + '.json')).read_text())
assert d['status'] == 'PASS_UNCONDITIONAL_COMPACT29_SAVINGS_DIAGNOSTIC'
for p, h in d['sources'].items():
    assert sha256((ROOT / p).read_bytes()).hexdigest() == h, p
assert len(d['sources']) == 9

# Enumerate multiplicities of every rank, independently of the producer's
# ascending-part recursion. All finite types at nonheavy ranks are irrelevant
# because their contribution cannot exceed the same rank-only old weight.
def rank_weight(r):
    return (r * (r - 1) if r < 5 else r * (r - 1) // 2 if r < 8
            else {8: 24, 9: 27}.get(r, r + 19))

def multiplicities(remaining, largest):
    if largest == 0:
        if remaining == 0:
            yield ()
        return
    for count in range(remaining // largest + 1):
        for rest in multiplicities(remaining - count * largest, largest - 1):
            yield rest + (largest,) * count

partitions = list(multiplicities(29, 29))
assert len(partitions) == len(set(partitions)) == 4565
families = {4: ('A4', 'B4', 'D4', 'F4', 'H4'),
            7: ('A7', 'B7', 'D7', 'E7'), 8: ('A8', 'B8', 'D8', 'E8'),
            9: ('A9', 'B9', 'D9'), 10: ('A10', 'B10', 'D10')}
savings = {4: Fraction(3, 4), 7: Fraction(3, 4), 8: Fraction(3, 4),
           9: Fraction(3, 4), 10: Fraction(1, 8)}
assert d['conditional_rank_savings'] == {str(k): str(v) for k, v in savings.items()}
heavy = {}
typed = {}
for p in partitions:
    old = sum(map(rank_weight, p))
    target = old - sum(savings.get(r, 0) for r in p)
    assert target <= 84
    if old <= 84:
        continue
    assert set(p) <= families.keys()
    counts = Counter(p)
    choices = [list(combinations_with_replacement(families[r], counts[r]))
               for r in sorted(counts)]
    rows = [tuple(t for group in groups for t in group) for groups in product(*choices)]
    assert all(row not in typed for row in rows)
    for row in rows:
        typed[row] = (old, old - 84, old - target, target)
    heavy[p] = (old, old - 84, len(rows), old - target, target)

assert len(heavy) == 9 and len(typed) == 2126
assert Counter(t[0] for t in typed.values()) == {85: 18, 86: 660, 87: 1448}
assert len(d['rank_inequalities']) == 9
assert len(d['above_threshold_typed_forests']) == 2126
assert len(d['exact_maximum_rank_patterns']) == 6
seen = set()
for row in d['rank_inequalities']:
    p = tuple(row['ranks'])
    assert p not in seen
    seen.add(p)
    assert (row['weight'], row['required_saving'], row['typed_count'],
            Fraction(row['conditional_saving']), Fraction(row['conditional_upper'])) == heavy[p]
    assert row['rank_multiplicities'] == {str(k): v for k, v in Counter(p).items()}
assert seen == heavy.keys()
seen = set()
for row in d['above_threshold_typed_forests']:
    p = tuple(row['types'])
    assert p not in seen
    seen.add(p)
    assert (row['old_weight'], row['required_saving'],
            Fraction(row['conditional_saving']), Fraction(row['conditional_upper'])) == typed[p]
assert seen == typed.keys()
assert {tuple(r['ranks']) for r in d['exact_maximum_rank_patterns']} == {
    p for p, row in heavy.items() if row[0] == 87}
assert d['weight_histogram'] == {'85': 18, '86': 660, '87': 1448}

# The two inequalities establish the stated minimax lower bounds. They do
# not certify that a geometric face rule can realize these deductions.
assert heavy[(7, 7, 7, 8)][1] == 3
assert 4 * max(savings.values()) == 3
assert heavy[(9, 10, 10)][1] == 1
assert savings[9] + 2 * savings[10] == 1
assert d['minimax_witness_rank_pattern'] == [7, 7, 7, 8]
assert d['rank10_witness_pattern'] == [9, 10, 10]
assert d['minimal_max_rank_saving'] == '3/4'
assert d['minimal_rank10_saving_at_that_cap'] == '1/8'
assert d['dimension'] == 29 and d['all_rank_partitions'] == 4565
assert d['threshold'] == d['conditional_maximum'] == '84'
assert d['current_maximum'] == '87'
assert d['existing_local_barrier_forest'] == ['E7', 'E7', 'E7', 'E8']
assert d['current_balance_lower_bound'] == '1050'
assert d['necessary_heavy_vertex_count'] == 350 and 3 * 350 == 1050
assert d['arithmetic_assumption'] is False
assert d['scalar_target_is_not_proved_face_rule'] is True
assert d['new_face_theorem'] is False
assert d['new_compact_dimension_bound'] is None
assert d['new_finite_volume_example'] is False
print('PASS saved compact29 scalar diagnostic:4565 partitions,2126 forests,9 inequalities; no geometry replay')
