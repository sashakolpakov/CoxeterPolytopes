#!/usr/bin/env python3
"""Exact common gap for current small-radius layers and all central contacts."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from fresh_h20_intrinsic_caps_20261014_quarter_gate import add, scale, mul, sign, pack, Z

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'fresh_h20_intrinsic_caps_20261014'
SOURCE = f'results/{PREFIX}_golden_seed.json'
OUT = ROOT / f'results/{PREFIX}_annular_gate.json'


def run():
    p = json.loads((ROOT/SOURCE).read_text())
    b = [[F(x) for x in row] for row in p['spatial_Gram']]
    sparse = [(i, j, x) for i, row in enumerate(b) for j, x in enumerate(row) if x]
    vv = [[tuple(F(x) for x in c) for c in row[:19]] for row in p['walls'][23:]]

    def dot(v, w):
        out = Z
        for i, j, x in sparse:
            out = add(out, scale(mul(v[i], w[j]), x))
        return out

    norms = [dot(v, v) for v in vv]
    assert set(norms) == {(F(a), F(4)) for a in (302, 326, 446)}
    assert all(sign(add(n, (F(-289), F(0)))) > 0 and
               sign(add(n, (F(-484), F(0)))) < 0 for n in norms)
    x = [scale(v, F(3, 4)) for v in vv[0]]
    xn = dot(x, x)
    assert sign(add(xn, (F(-169), F(0)))) > 0
    assert sign(add(xn, (F(-196), F(0)))) < 0
    gaps = []
    for row in p['walls'][:22]:
        r = [tuple(F(y) for y in c) for c in row]
        gap = add(r[19], scale(dot(r[:19], x), -1))
        assert sign(gap) > 0
        gaps.append(pack(gap))
    # X=(x,(|x|^2-1)/2,1) has square1 and is on the center's future side.
    escape = x+[scale(add(xn, (F(-1), F(0))), F(1, 2)), (F(1), F(0))]
    center_pairing = add(scale(escape[19], -1), (F(-9, 2), F(0)))
    assert sign(center_pairing) < 0
    sources = [SOURCE, f'scripts/{PREFIX}_quarter_gate.py', str(Path(__file__).resolve().relative_to(ROOT))]
    return {'status': 'PASS_EXACT_COMMON_ANNULAR_ESCAPE_FOR_ALL_CURRENT_SMALL_RADIUS_BANKS',
            'sources': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
            'field': 'Q(phi), phi^2=phi+1, phi>1',
            'all252_vertex_norm_squares': sorted([pack(n) for n in set(norms)]),
            'vertex_norm_strict_bounds': ['17', '22'],
            'test_point': '3*z_0/4', 'test_point_norm_square': pack(xn),
            'test_point_norm_strict_bounds': ['13', '14'], 'seed_gaps': gaps,
            'ball_families': [
                {'centers': 'all old product vertices', 'radii': 'at most1',
                 'distance_to_test_point': '>17-14=3', 'ball_power': '>8'},
                {'centers': 'one half of any point of the product alcove', 'radii': 'at most1',
                 'distance_to_test_point': '>13-11=2', 'ball_power': '>3'},
                {'centers': 'one quarter of any point of the product alcove', 'radii': 'at most3/2',
                 'distance_to_test_point': '>13-11/2=15/2', 'ball_power': '>54'},
                {'centers': 'any center whose sphere meets the central radius2 sphere at an acute angle or tangency',
                 'radii': '0<R<=2', 'center_norm': '<=R+2<=4',
                 'distance_to_test_point': '>9', 'ball_power': '>77'},
                {'centers': 'origin', 'radii': 'at most2', 'ball_power': '>165'}],
            'escape_point': [pack(v) for v in escape], 'escape_norm': '1',
            'escape_pairing_with_known_strict_center': pack(center_pairing),
            'conclusion': 'Every finite selection from all displayed families leaves this same strict spacelike escape, independently of pairwise compatibility or wall replacement.',
            'scope': 'Allows arbitrary points of the alcove in the half and quarter families and all real central-contact angles/radii in the displayed ranges. Larger endpoint caps and other placements are outside the theorem.'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    if args.verify:
        assert json.loads(OUT.read_text()) == out
    else:
        OUT.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'])
