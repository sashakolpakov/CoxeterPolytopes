#!/usr/bin/env python3
"""Constant-size exact proof certificate: at most one H4 restoration at A5."""
import argparse
import fcntl
from fractions import Fraction as F
from hashlib import file_digest
import json
import os
import signal
import fresh_compact21_d5_a1_triangle_coupling_20261018 as A
import run_lean_batches_20261018 as G

ROOT, P = A.ROOT, A.P
TAG = 'fresh_compact21_a5_h4_source_capacity_20261018'


def run():
    edges = {(i, i+1): 3 for i in range(4)}
    inverse = P.inverse(P.base_gram(5, edges, 5))
    assert inverse[0][0] == inverse[4][4] == P.K(F(5, 6))
    assert inverse[0][4] == P.K(F(1, 6))
    phi2 = P.COS[5]*P.COS[5]
    height = phi2*inverse[0][0]-2
    projection = phi2*inverse[0][4]
    mutual = projection-height
    assert height == P.K((F(-3, 4), 0, F(5, 12), 0)) and height.sign() > 0
    assert mutual == P.K((1, 0, F(-1, 3), 0)) and mutual.sign() > 0
    result = P.enc({
        'status': 'PASS_A5_H4_SOURCE_ANGLE_CAPACITY_ONE',
        'forest': ['A5', 'E8', 'E8'], 'A5_path': list(range(5)),
        'potential_restored_pairs': [[0, 3], [1, 4]],
        'reverse_H4_wall_contacts': [[[0, 5]], [[4, 5]]],
        'inverse_endpoint_diagonal': inverse[0][0],
        'inverse_opposite_endpoints': inverse[0][4],
        'common_height_squared': height,
        'projection_product': projection, 'forced_mutual_product': mutual,
        'forced_product_positive': True, 'old_angle_capacity': 2, 'new_angle_capacity': 1,
        'new_source_capacities_A4_A5_B4_B5': [1, 1, 1, 1],
        'new_dimension_bound': None, 'Lean_verified': False,
    })
    result['sources'] = {}
    for path in (f'scripts/{TAG}.py', 'scripts/fresh_compact29_small_faces_20261015.py',
                 'results/fresh_compact21_h4_triangle_funded_transfer_20261018_freeze.json',
                 'results/fresh_compact21_h4_a1_triangle_priority_20261018.json'):
        with (ROOT/path).open('rb') as stream:
            result['sources'][path] = file_digest(stream, 'sha256').hexdigest()
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, G.interrupt); signal.alarm(59)
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        G.resources()
        data = run()
        target = ROOT/f'results/{TAG}.json'
        if args.verify:
            assert json.loads(target.read_text()) == data
        else:
            assert not target.exists()
            temporary = target.with_suffix('.json.tmp')
            try:
                with temporary.open('x') as stream:
                    json.dump(data, stream, indent=2); stream.write('\n')
                    stream.flush(); os.fsync(stream.fileno())
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
        G.resources(start=False)
    signal.alarm(0)
    print(json.dumps({'status': data['status'], 'capacity': data['new_angle_capacity'],
                      'mutual_product': data['forced_mutual_product'], 'new_dimension_bound': None}))
