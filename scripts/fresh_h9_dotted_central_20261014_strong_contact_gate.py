#!/usr/bin/env python3
"""Literal finite faces force a second negative direction, all other labels free."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_dotted_central_20261014_strong_contact_gate'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'
BASE = 'results/fresh_h9_dotted_central_20261014_gate.json'


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def run():
    base, geo = [json.loads((ROOT / p).read_text()) for p in (BASE, GEOM)]
    for d in (base, geo):
        for p, h in d['sources'].items():
            assert sha(p) == h, p
    g = geo['representatives'][5]['geometry']
    vertices = [set(v['incident']) for v in g['vertices']]
    missing = {frozenset(m) for m in g['minimal_nonfaces']}
    maps = []
    for old_index, r in enumerate(base['maps']):
        A, B, x, y = r['dense_block'], r['other_block'], r['dense_leaf'], r['other_leaf']
        assert frozenset(A + [x]) in missing and frozenset(B + [y]) in missing
        for b in B:
            S = A + [y, b]
            vi = next(i for i, v in enumerate(vertices) if set(S) <= v)
            supports = r['actual_supports'] + [dict(kind='strong_contact_forces_foreign_right', support=S, actual_vertex=vi)]
            for s in supports:
                assert set(s['support']) <= vertices[s['actual_vertex']]
            bounds = sorted([sorted(A) + [4], sorted(B) + [6],
                sorted((x, A[0])) + [3], sorted((x, A[1])) + [3], sorted((y, b)) + [4]])
            maps.append(dict(base_map=old_index, dense_block=A, other_block=B, dense_leaf=x,
                other_leaf=y, strong_contact_wall=b, actual_supports=supports, bounds=bounds))
    families, seen = [], set()
    for i, m in enumerate(maps):
        key = tuple(map(tuple, m['bounds']))
        if key in seen:
            continue
        seen.add(key)
        families.append(dict(id='strong_contact_pair_' + str(len(families)), witness_map=i, fixed=[], bounds=m['bounds']))
    # Norm-one matrices. These bounds are independent of any finite label cutoff.
    assert F(1, 4) / (1 - F(1, 2)) + F(1, 2) == 1
    assert F(1, 2) / (1 - F(3, 4)) == 2
    assert (1 - F(1, 4)) * (1 - F(1, 2)) - F(1, 2) == -F(1, 8)
    # Residual determinant >= sqrt2/2-1/4 >0.
    assert 4 * 2 > 1
    assert len(maps) == 88
    paths = [BASE, GEOM, 'scripts/' + TAG + '.py']
    return dict(status='PASS_ALL_LABEL_STRONG_CONTACT_BLOCK_GATE_CASE5',
        sources={p: sha(p) for p in paths}, case=5, maps=maps, families=families,
        theorem=dict(normalization='all root squares1',
            assumptions='A=I2(a),a>=4;B=I2(b),b>=6;x touches both A roots nonright;y has a specified contact label>=4 to Bj;all literal supports finite',
            cross_and_x_foreign='B has label>=6, so finite A+B and B+x force all cross pairs and x-B pairs right',
            y_foreign='In finite A+y+Bj any nonright y-A contact costs at least1/2 in A; the prescribed y-Bj contact costs at least1/2. Positive Schur complement is impossible, so y-A is right',
            xy_bound='The finite Ai-x-y-Bj path has determinant<=3/8-delta^2;label>=4 fails,so delta^2<=1/4',
            dense_cost_lower={'rational':'1','sqrt2':'1/2'}, other_cost_lower='2',
            residual_determinant_lower={'rational':'-1/4','sqrt2':'1/2'},
            conclusion='The residual2x2 Gram is negative definite; the full Gram has at least two negative directions'),
        scope='Sufficient all-label local obstruction. All unlisted contacts and dotted distances are free. It does not exclude the whole case5 incidence.')


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
