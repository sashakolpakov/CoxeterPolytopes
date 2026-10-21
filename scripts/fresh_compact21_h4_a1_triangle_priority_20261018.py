#!/usr/bin/env python3
"""Complete bounded H4+A1 triangle test for its label-five extra unit.

High bonds not represented here are treated explicitly in the written proof.
No A1 spectator contact is assumed zero in the finite domain.
"""
import argparse
from collections import Counter
import fcntl
from fractions import Fraction as F
from hashlib import file_digest, sha256
from itertools import combinations, product
import json
import os
from pathlib import Path
import resource
import signal
import sys
import fresh_compact21_d5_a1_triangle_coupling_20261018 as A
import run_lean_batches_20261018 as G

P = A.P
ROOT = A.ROOT
TAG = 'fresh_compact21_h4_a1_triangle_priority_20261018'
EDGES = {(0, 1): 5, (1, 2): 3, (2, 3): 3}


def integral(x):
    a, b, c, d = x.v
    return not b and not d and (2*c).denominator == 1 and (a-c).denominator == 1


def payment(nodes, edges, marks):
    blocks, adjacency = A.parts(nodes, edges)
    block = next(b for b in blocks if marks[0] in b)
    if marks[1] not in block:
        return {'ordinary': 0, 'h4_extra': 0, 'total': 0}
    distance = {marks[0]: 0}
    todo = [marks[0]]
    for i in todo:
        for j in sorted(adjacency[i]):
            if j not in distance:
                distance[j] = distance[i]+1
                todo.append(j)
    dd = distance[marks[1]]
    ordinary = int(dd <= 2)
    h4 = len(block) == 4 and any(m == 5 and set(pair) <= set(block)
                                for pair, m in edges.items())
    extra = int(h4 and (dd == 3 or edges.get(tuple(sorted(marks))) == 5))
    removable = int(h4 and edges.get(tuple(sorted(marks))) == 5)
    return {'ordinary': ordinary, 'h4_extra': extra,
            'five_edge_extra': removable, 'total': ordinary+extra,
            'after_deleting_five_edge_extra': ordinary+extra-removable}


def run():
    gram = P.base_gram(5, EDGES, 5)
    inverse = P.inverse(gram)
    assert P.psd_status(gram) == 1
    roots, eligible, counts, ledger = [], [], Counter(), sha256()
    for labels in product(range(2, 6), repeat=5):
        contacts = {i: m for i, m in enumerate(labels) if m > 2}
        edges = EDGES | {(i, 5): m for i, m in contacts.items()}
        finite = all(A.finite([i for i in range(6) if i != omit], edges)
                     for omit in (0, 1))
        ledger.update(json.dumps([labels, finite], separators=(',', ':')).encode()+b'\n')
        if not finite:
            continue
        row = {'contacts': sorted(contacts.items())}
        height = sum((P.COS[mi]*inverse[i][j]*P.COS[mj]
                      for i, mi in contacts.items() for j, mj in contacts.items()), P.K(-2))
        row['height_squared'] = height
        if height.sign() <= 0:
            status = 'nonpositive_height'
        elif any((P.COS[mi]*P.COS[mj]).v[k]
                 for (i, mi), (j, mj) in combinations(contacts.items(), 2)
                 if i < 4 and j < 4 for k in (1, 3)):
            status = 'cycle_outside_Qsqrt5'
        else:
            assert not height.v[1] and not height.v[3]
            conjugate = P.K((height.v[0], 0, -height.v[2], 0))
            row['conjugate_height_squared'] = conjugate
            if conjugate.sign() > 0:
                status = 'conjugate_not_PSD'
            else:
                witness = None
                full = P.base_gram(5, edges, 6)
                for size in range(6):
                    for subset in combinations(range(5), size):
                        nodes = list(subset)+[5]
                        blocks, _ = A.parts(nodes, edges)
                        if len(blocks) != 1 or A.finite(nodes, edges):
                            continue
                        if any(m == 5 and set(pair) <= set(nodes) for pair, m in edges.items()):
                            continue
                        if P.psd_status(P.sub(full, nodes)) == 0:
                            witness = nodes
                            break
                    if witness is not None:
                        break
                row['parabolic_subset'] = witness
                if witness is not None:
                    status = 'parabolic_obstruction'
                else:
                    status = 'necessary_arithmetic_triangle_root'
                    corners = []
                    for omitted in (0, 1, 5):
                        marks = [i for i in (0, 1, 5) if i != omitted]
                        pp = payment([i for i in range(6) if i != omitted], edges, marks)
                        pp.setdefault('five_edge_extra', 0)
                        pp.setdefault('after_deleting_five_edge_extra', pp['total'])
                        corners.append({'omitted': omitted, 'marks': marks, **pp})
                    row['corners'] = corners
                    row['old_triangle_payment'] = sum(r['total'] for r in corners)
                    row['payment_after_all_five_edge_deletions'] = sum(
                        r['after_deleting_five_edge_extra'] for r in corners)
                    row['E8_height_products'] = [height * h for h in
                                                (P.K((4, 0, 2, 0)), P.K((1, 0, 1, 0)))]
                    row['E8_pair_cycle_integral'] = [integral(x) for x in row['E8_height_products']]
                    roots.append(row.copy())
        row['status'] = status
        counts[status] += 1
        eligible.append(row)
    # No optimistic assertion: the tiny certificate records the actual domain.
    dangerous = [i for i, r in enumerate(roots)
                 if r['payment_after_all_five_edge_deletions'] < 2]
    compatible = [i for i in dangerous if any(roots[i]['E8_pair_cycle_integral'])]
    result = P.enc({
        'status': 'PASS_COMPLETE_H4_A1_FIVE_EDGE_TRIANGLE_DOMAIN',
        'forest': ['A1', 'H4', 'E8', 'E8'], 'ground_field': 'Q(sqrt5)',
        'marked_H4_pair': [0, 1], 'A1_node': 4,
        'finite_bounded_vectors': 4**5, 'finite_domain_digest': ledger.hexdigest(),
        'status_counts': dict(counts), 'eligible': eligible, 'surviving_roots': roots,
        'dangerous_triangle_indices': dangerous,
        'dangerous_with_integral_E8_pair_cycle': compatible,
        'conditional_five_edge_extra_deletion_proved_by_this_test': not compatible,
        'new_dimension_bound': None, 'Lean_verified': False,
    })
    result['sources'] = {}
    for path in (f'scripts/{TAG}.py', 'scripts/fresh_compact21_d5_a1_triangle_coupling_20261018.py',
                 'scripts/fresh_compact29_small_faces_20261015.py',
                 'results/fresh_compact21_single_field_20261016.md',
                 'results/fresh_compact21_supported_forest_triangle_coupling_20261018.json'):
        with (ROOT/path).open('rb') as stream:
            result['sources'][path] = file_digest(stream, 'sha256').hexdigest()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, G.interrupt)
    signal.alarm(59)
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        G.resources()
        result = run()
        target = ROOT/f'results/{TAG}.json'
        if args.verify:
            assert json.loads(target.read_text()) == result
        else:
            assert not target.exists()
            temporary = target.with_suffix('.json.tmp')
            try:
                with temporary.open('x') as stream:
                    json.dump(result, stream, indent=2)
                    stream.write('\n')
                    stream.flush()
                    os.fsync(stream.fileno())
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
        G.resources(start=False)
    signal.alarm(0)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
    assert peak < 1024**3
    print(json.dumps({'status': result['status'], 'roots': len(result['surviving_roots']),
                      'dangerous': result['dangerous_triangle_indices'],
                      'dangerous_E8_compatible': result['dangerous_with_integral_E8_pair_cycle'],
                      'status_counts': result['status_counts'], 'peak_rss_bytes': peak}))


if __name__ == '__main__':
    main()
