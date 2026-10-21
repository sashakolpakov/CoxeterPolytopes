#!/usr/bin/env python3
"""Scoped higher-precision replay; predecessor enclosures remain at 80 digits."""
import copy
import json
import time
from pathlib import Path
import fresh_h9_h5_closure_next_20261013_verify as B

V = B.V
TAG = 'fresh_h9_h5_closure_next_20261013'


def run():
    start = time.monotonic()
    assert V.SCALE == 10**80
    saved = V.read(B.PRE + '_verify.json')
    state = list(B.old_state())
    factor = 10**80

    def lift(interval):
        return tuple(int(x) * factor for x in interval)

    state[3] = [[lift(x) for x in row] for row in state[3]]
    state[4] = [[lift(x) for x in row] for row in state[4]]
    state[6] = [lift(x) for x in state[6]]
    V.SCALE = 10**160
    V.PHI = V.div(V.add(V.rat(1), V.sqrt(V.rat(5))), V.rat(2))
    V.radical.cache_clear()
    B.old_state = lambda: copy.deepcopy(state)
    write = Path.write_text
    blocked = []

    def intercept(path, text, *args, **kwargs):
        assert path == V.ROOT / (B.PRE + '_verify.json'), path
        blocked.append(path)
        return len(text)

    Path.write_text = intercept
    try:
        out = B.run(False)
    finally:
        Path.write_text = write
    assert len(blocked) == 1
    assert out['sources'] == saved['sources']
    for path, digest in saved['sources'].items():
        assert V.sha(path) == digest
    for key in ('facets', 'new_pair_count', 'new_norms', 'positive_affine_census'):
        assert out[key] == saved[key]
    for key in ('fixed_walls', 'facets', 'vertex_census', 'finite_components',
                'finite_edges', 'open_edges', 'edges'):
        assert out['H5_face'][key] == saved['H5_face'][key]
    for key in ('fixed_walls', 'facets', 'sides'):
        assert out['compact_H2_hexagon'][key] == saved['compact_H2_hexagon'][key]
    review = {
        'status': 'PASS_SCOPED_160_DIGIT_TWO_CAP_RECONSTRUCTION',
        'seconds': round(time.monotonic() - start, 2),
        'sources': saved['sources'],
        'review_script_sha256': V.sha('scripts/' + TAG + '_review.py'),
        'certificate_sha256': V.sha(B.PRE + '_verify.json'),
        'facets': out['facets'],
        'new_pair_count': out['new_pair_count'],
        'H5_finite_vertices': len(out['H5_face']['finite_vertices']),
        'H5_finite_edges': out['H5_face']['finite_edges'],
        'H5_open_edges': out['H5_face']['open_edges'],
        'compact_H2_vertices': len(out['compact_H2_hexagon']['finite_vertices']),
        'scope': 'Reconstructs the two new caps and replays all producer gates at '
                 '160-digit interval scale using shared routines. The 123 predecessor '
                 'normals retain exactly their frozen 80-digit rational enclosures. '
                 'No producer output is overwritten. This is not a separate '
                 'arithmetic implementation or a 160-digit predecessor reconstruction.'
    }
    (V.ROOT / 'results' / (TAG + '_review.json')).write_text(
        json.dumps(review, indent=2) + '\n')
    print(review['status'], review['H5_finite_vertices'], review['H5_open_edges'])
    return review


if __name__ == '__main__':
    run()
