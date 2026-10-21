#!/usr/bin/env python3
"""Complete D5 nonadjacent triangle domain with one whole A1 spectator.

At an A1+D5+E7+E8 compact arithmetic rank21 vertex, E7/E8 spectators
are orthogonal to a D5 triangular third wall. Positive height excludes
an isolated A1 contact of label>=6. The companion report supplies these
geometric arguments. No weight rule, bound20, or chamber is produced.
"""
import argparse
from collections import Counter
import fcntl
from fractions import Fraction as Q
from hashlib import file_digest, sha256
from itertools import combinations, product
import json
from pathlib import Path
import resource
import signal
import sys
import fresh_compact29_small_faces_20261015 as P

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact21_d5_a1_triangle_coupling_20261018'
EDGES = {(0, 1): 3, (1, 2): 3, (2, 3): 3, (2, 4): 3}


def parts(nodes, edges):
    left = set(nodes)
    adjacency = {i: set() for i in nodes}
    for (i, j), label in edges.items():
        if label > 2 and i in adjacency and j in adjacency:
            adjacency[i].add(j)
            adjacency[j].add(i)
    blocks = []
    while left:
        block = [min(left)]
        left.remove(block[0])
        for i in block:
            new = sorted(adjacency[i] & left)
            block.extend(new)
            left.difference_update(new)
        blocks.append(block)
    return blocks, adjacency


def finite(nodes, edges):
    blocks, adjacency = parts(nodes, edges)
    for block in blocks:
        n = len(block)
        if n <= 2:
            continue
        bonds = {pair: m for pair, m in edges.items() if set(pair) <= set(block)}
        if len(bonds) != n-1 or max(map(lambda i: len(adjacency[i]), block)) > 3:
            return False
        forks = [i for i in block if len(adjacency[i]) == 3]
        if not forks:
            current = next(i for i in block if len(adjacency[i]) == 1)
            previous = None
            labels = []
            while adjacency[current]-{previous}:
                nxt, = adjacency[current]-{previous}
                labels.append(bonds[tuple(sorted((current, nxt)))])
                previous, current = current, nxt
            heavy = [(i, m) for i, m in enumerate(labels) if m != 3]
            if not heavy:
                continue
            if len(heavy) != 1:
                return False
            i, m = heavy[0]
            if m == 4 and (i in (0, n-2) or (n == 4 and i == 1)):
                continue
            if m == 5 and n in (3, 4) and i in (0, n-2):
                continue
            return False
        if len(forks) != 1 or any(m != 3 for m in bonds.values()):
            return False
        arms = []
        for neighbor in adjacency[forks[0]]:
            previous, current, length = forks[0], neighbor, 1
            while adjacency[current]-{previous}:
                nxt, = adjacency[current]-{previous}
                previous, current, length = current, nxt, length+1
            arms.append(length)
        arms.sort()
        if arms[:2] != [1, 1] and tuple(arms) not in ((1, 2, 2), (1, 2, 3), (1, 2, 4)):
            return False
    return True


def digest(path):
    with path.open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def run():
    gram = P.base_gram(6, EDGES, 6)
    inv = P.inverse(gram)
    assert P.psd_status(gram) == 1
    cache = {}
    rows = []
    roots = {}
    totals = Counter()
    for a, b in combinations(range(5), 2):
        if (a, b) in EDGES:
            continue
        ledger = sha256()
        eligible = []
        for iteration, labels in enumerate(product(range(2, 6), repeat=6)):
            if iteration % 128 == 0:
                maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                maximum *= 1 if sys.platform == 'darwin' else 1024
                if maximum > 1024**3:
                    raise MemoryError('One-GiB RSS guard exceeded; no certificate written')
            contacts = {i: m for i, m in enumerate(labels) if m > 2}
            edges = EDGES | {(i, 6): m for i, m in contacts.items()}
            corners = [finite([i for i in range(7) if i != omitted], edges)
                       for omitted in (a, b)]
            ledger.update(json.dumps([labels, corners], separators=(',', ':')).encode()+b'\n')
            if not all(corners):
                continue
            h = sum((P.COS[mi]*inv[i][j]*P.COS[mj]
                     for i, mi in contacts.items() for j, mj in contacts.items()), P.K(-2))
            row = {'labels': labels, 'height_squared': h}
            if h.sign() <= 0:
                row['status'] = 'nonpositive_height'
            else:
                cycles = [P.COS[contacts[i]]*P.COS[contacts[j]]
                          for i, j in combinations([i for i in contacts if i < 5], 2)]
                if any(z.v[1] or z.v[3] for z in cycles):
                    row['status'] = 'cycle_outside_Qsqrt5'
                else:
                    assert not (h.v[1] or h.v[3])
                    conjugate = P.K((h.v[0], 0, -h.v[2], 0))
                    row['conjugate_height_squared'] = conjugate
                    if conjugate.sign() > 0:
                        row['status'] = 'conjugate_not_PSD'
                    else:
                        if labels not in cache:
                            g = P.base_gram(6, edges, 7)
                            witness = None
                            for size in range(7):
                                for subset in combinations(range(6), size):
                                    indices = list(subset)+[6]
                                    blocks, _ = parts(indices, edges)
                                    if len(blocks) != 1 or any(
                                        m == 5 and i in indices and j in indices
                                        for (i, j), m in edges.items()):
                                        continue
                                    if P.psd_status(P.sub(g, indices)) == 0:
                                        witness = indices
                                        break
                                if witness is not None:
                                    break
                            cache[labels] = witness
                        row['parabolic_subset'] = cache[labels]
                        row['status'] = ('parabolic_obstruction' if cache[labels] is not None
                                         else 'necessary_arithmetic_triangle_root')
                        if cache[labels] is None:
                            roots[labels] = {'contacts': sorted(contacts.items()), 'height_squared': h}
            totals[row['status']] += 1
            eligible.append(row)
        rows.append({'marked': [a, b], 'all_vectors': 4**6,
                     'finite_domain_digest': ledger.hexdigest(), 'eligible': eligible})
    expected = [((0, 5),), ((3, 5),), ((1, 3), (3, 5)),
                ((4, 5),), ((1, 3), (4, 5)), ((0, 3), (5, 5))]
    assert {tuple(r['contacts']) for r in roots.values()} == set(expected)
    old = json.loads((ROOT/'results/fresh_compact21_d5_two_e8_triangle_coupling_20261018.json').read_text())
    old_roots = [r for r in old['surviving_roots'] if r['component'] == 'D5']
    assert len(old_roots) == 5
    root_rows = list(roots.values())
    pairs = []
    for i, j in combinations(range(len(root_rows)), 2):
        first, second = root_rows[i], root_rows[j]
        aa, bb = dict(first['contacts']), dict(second['contacts'])
        p = sum((P.COS[mi]*inv[x][y]*P.COS[mj]
                 for x, mi in aa.items() for y, mj in bb.items()), P.K())
        margin = p*p-first['height_squared']*second['height_squared']
        assert p.sign() > 0 and margin.sign() > 0
        pairs.append({'roots': [i, j], 'projection_product': p,
                      'positive_square_margin': margin})
    sources = [f'scripts/{TAG}.py', 'scripts/fresh_compact29_small_faces_20261015.py',
               'scripts/fresh_compact21_d5_two_e8_triangle_coupling_20261018.py',
               'results/fresh_compact21_d5_two_e8_triangle_coupling_20261018.json',
               'results/fresh_compact21_single_field_20261016.json',
               'results/fresh_compact21_single_field_20261016.md',
               'results/review_compact21_d5_two_e8_triangle_coupling_20261018.md']
    return P.enc({'status': 'PASS_COMPACT21_D5_A1_COMPLETE_NONADJACENT_TRIANGLE_DOMAIN',
                  'sources': {path: digest(ROOT/path) for path in sources},
                  'vertex_forest': ['A1', 'D5', 'E7', 'E8'], 'ground_field': 'Q(sqrt5)',
                  'old_nodes': {'D5': list(range(5)), 'A1': 5},
                  'high_A1_labels_excluded_by_positive_height': True,
                  'untouched_E7_E8_forced_orthogonal': True,
                  'all_vectors': 6*4**6, 'domains': rows, 'status_counts': dict(totals),
                  'surviving_roots': root_rows, 'acute_excluded_pairs': pairs,
                  'D5_designated_triangle_maximum': 1,
                  'changed_old_face_payments': False, 'new_charge_saving': False,
                  'new_dimension_bound': None, 'Lean_verified': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    # This tiny computation uses one thread and bounded standard-library state.
    # macOS rejects lowering these rlimits in this managed environment; the
    # periodic RSS guard above is explicit rather than claiming an OS cap.
    signal.alarm(59)
    lock_path = ROOT/'formal/.proof-batches/exclusive.lock'
    lock_path.parent.mkdir(exist_ok=True)
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = run()
    target = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert json.loads(target.read_text()) == result
    else:
        assert not target.exists()
        target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'vectors': result['all_vectors'],
                      'roots': len(result['surviving_roots']),
                      'acute_excluded_pairs': len(result['acute_excluded_pairs']),
                      'D5_triangle_maximum': 1, 'new_dimension_bound': None,
                      'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss *
                          (1 if sys.platform == 'darwin' else 1024)}))
