#!/usr/bin/env python3
"""Produce finite certificates for two continuous compact-filling obstructions.

This is independent of the numerical Gram-mutation scout. Facet indices are
zero based; the source Gram convention is reversed and scaled to diagonal 2.
"""
import hashlib
import json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/enumerations/alexeev_h9_doubling_cubes.json'
OUT = ROOT / 'data/enumerations/fresh_h9_cusp_filling_obstructions_20261013.json'


def run():
    source = json.loads(SOURCE.read_text())
    ids = ['T1_D16_W4_7_12_14', 'T2_D16_W1_5_9_14']
    grams = [[[-F(x) / 2 for x in row] for row in next(
        r for r in source['certified_types'] if r['id'] == rid)['facet_gram']]
        for rid in ids]
    G = grams[0]
    circuits = []
    for nodes in combinations(range(18), 4):
        if all(sum(G[i][j] == -1 for j in nodes) == 2 for i in nodes) and all(
                G[i][j] in (0, -1) for i, j in combinations(nodes, 2)):
            circuits.append({'nodes': list(nodes), 'kernel': [1] * 4,
                             'type': 'affine_A3'})
    for center in range(18):
        neighbors = [i for i in range(18) if G[center][i] == -1]
        for leaves in combinations(neighbors, 4):
            if all(G[i][j] == 0 for i, j in combinations(leaves, 2)):
                circuits.append({'nodes': [center, *leaves],
                                 'kernel': [2, 1, 1, 1, 1],
                                 'type': 'affine_D4'})
    cover = []
    for removed in combinations(range(18), 5):
        chosen = next(k for k, c in enumerate(circuits)
                      if not set(removed).intersection(c['nodes']))
        cover.append({'deleted': list(removed), 'circuit': chosen})
    # Hitting every displayed affine circuit needs at least six deletions.
    six = next(s for s in combinations(range(18), 6)
                if all(set(s).intersection(c['nodes']) for c in circuits))
    H = grams[1]
    centers = [i for i in range(18) if sum(x == -1 for x in H[i]) == 8]
    leaves = [i for i in range(18) if i not in centers]
    pairs = [[i, j] for i, j in combinations(leaves, 2) if H[i][j] != 0]
    counts = Counter(len(set(leaves) - set(s))
                     for s in combinations(range(18), 5))
    result = {
        'schema': 'H9-cusp-filling-continuous-obstructions-v1',
        'source': str(SOURCE.relative_to(ROOT)),
        'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'normalization': 'minus source Gram divided by 2; diagonal 2',
        'global_assumptions': 'All changed and additional outward normal pairings '
                              'are nonpositive, and a strict timelike interior exists.',
        'T1': {'source_id': ids[0], 'circuits': circuits,
               'five_deletion_cover': cover,
               'six_deletion_circuit_transversal': list(six),
               'scope': 'At most five old facet deletions. Change old parallel '
                        'pairings arbitrarily; preserve all other old pairings. '
                        'Any further walls must be acute compatible. No compact '
                        'H9 chamber with a strict interior is possible.'},
        'T2': {'source_id': ids[1], 'centers': centers, 'leaves': leaves,
               'leaf_matching': pairs,
               'five_deletion_leaf_count_histogram': dict(counts),
               'ambient_rank_bound': 10,
               'scope': 'At most five old facet deletions. Change any nonzero '
                        'old pairings, preserving old orthogonalities. Any '
                        'further walls must be acute compatible. No compact '
                        'H9 chamber with a strict interior is possible.'},
        'decision': 'EXACT_SCOPED_COMPACT_OBSTRUCTIONS; finite volume not excluded'}
    OUT.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'output': str(OUT.relative_to(ROOT)),
                      'circuits': dict(Counter(c['type'] for c in circuits)),
                      'covered_deletions': len(cover),
                      'T2_leaf_counts': dict(counts)}))


if __name__ == '__main__':
    run()
