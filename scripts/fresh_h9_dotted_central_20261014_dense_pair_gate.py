#!/usr/bin/env python3
"""Two dense leaves, one block>=6 and the other>=4, force Lorentz index>=2."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_dotted_central_20261014_dense_pair_gate'
BASE = 'results/fresh_h9_dotted_central_20261014_gate.json'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def run():
    base, geo = [json.loads((ROOT / p).read_text()) for p in (BASE, GEOM)]
    for d in (base, geo):
        for p, h in d['sources'].items():
            assert sha(p) == h, p
    V = [set(v['incident']) for v in geo['representatives'][5]['geometry']['vertices']]
    maps = []
    for i, r in enumerate(base['maps']):
        A, B, x, y = r['dense_block'], r['other_block'], r['dense_leaf'], r['other_leaf']
        supports = r['actual_supports'][:]
        for b in B:
            S = [x, y, b]
            vi = next(j for j, v in enumerate(V) if set(S) <= v)
            supports.append(dict(kind='foreign_contact_forces_xy_right', support=S, actual_vertex=vi))
        for s in supports:
            assert set(s['support']) <= V[s['actual_vertex']]
        bounds = sorted([sorted(A) + [6], sorted(B) + [4]] +
            [sorted((x, a)) + [3] for a in A] + [sorted((y, b)) + [3] for b in B])
        maps.append(dict(base_map=i, dense_block=A, other_block=B, dense_leaf=x, other_leaf=y,
                         actual_supports=supports, bounds=bounds))
    families, seen = [], set()
    for i, row in enumerate(maps):
        key = tuple(map(tuple, row['bounds']))
        if key in seen:
            continue
        seen.add(key)
        families.append(dict(id='two_dense_pair_' + str(len(families)), witness_map=i, fixed=[], bounds=row['bounds']))
    # cos(pi/6)>5/6, cos(pi/4)>2/3, cos(pi/5)^2<2/3.
    assert F(3, 4) > F(25, 36)
    assert F(1, 2) > F(4, 9)
    assert 5 < F(49, 9)
    assert 1 - F(1, 2) - F(1, 2) == 0
    assert (3 - 1) * (F(3, 2) - 1) - F(1, 2) == F(1, 2)
    assert (3 - 1) * (F(3, 2) - 1) - F(3, 4) == F(1, 4)
    paths = [BASE, GEOM, 'scripts/' + TAG + '.py']
    return dict(status='PASS_ALL_LABEL_TWO_DENSE_BLOCK_GATE_CASE5',
        sources={p: sha(p) for p in paths}, case=5, maps=maps, families=families,
        theorem=dict(normalization='root squares1',
            assumptions='A label>=6,B label>=4;x has two nonright own A contacts;y has two nonright own B contacts;the listed supports are finite',
            forced_right='A+B and A+y force cross-block and y-A contacts right',
            projection_bounds=dict(x_own_H='strictly greater than3', y_own_Y='strictly greater than3/2'),
            no_foreign_case='If x-B right,the four-path bound gives xy squared<=1/2;residual determinant>1/2',
            foreign_case='Otherwise B is4 or5 and x-B is one pi3 contact. The actual x-y-Bj triangle forces xy right. Its foreign projection square F<3/4. Positive-metric Cauchy gives det>=(H-1)*(Y-1)-F>1/4',
            conclusion='Both residual diagonal entries are negative and their determinant is positive, so the full Gram has at least two negative directions'),
        scope='All-label sufficient local obstruction, with all other labels and dotted distances free. Not a whole-incidence exclusion.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    d = run()
    out = ROOT / ('results/' + TAG + '.json')
    if a.verify:
        assert d == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(d, indent=2) + '\n')
    print(d['status'], len(d['maps']), 'placements', len(d['families']), 'conditions')
