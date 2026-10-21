#!/usr/bin/env python3
"""Exact interval incompatibility cover for a literal polar bank.

All old walls are optional. Each witness is strict for the full old
chamber, so any old-wall deletion retains every witness.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import time
import fresh_h9_split_multiple_edge_escape_20261013 as I
from fresh_h9_split_leaf_3034_safe_graph_20261013 import angle_enclosures, AngleGate, PAIR_SCALE

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_makarov_polar_closure_20261013'
SOURCE = ROOT/'results/fresh_h9_makarov_joint_20261013_exact.json'
RECIPES = ROOT/'results'/f'{TAG}_recipes.json'


def k(x):
    return I.add(I.rat(x[0]), I.mul(I.rat(x[1]), I.sqrt(I.rat(5))))


def indexes(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask ^= bit


def run(verify=False):
    start = time.time()
    s = json.loads(SOURCE.read_text())
    recipes = json.loads(RECIPES.read_text())
    assert recipes['sources'][str(SOURCE.relative_to(ROOT))] == sha256(SOURCE.read_bytes()).hexdigest()
    assert len(recipes['polar_recipes']) == 1564 and len(recipes['witnesses']) == 21
    Q = [[k(x) for x in row] for row in s['form']]
    center = [k(x) for x in s['strict_center']]
    core = [[k(x) for x in row] for row in s['seed_normals']]
    for row in s['caps']:
        d = I.sqrt(k(row['radicand']))
        core.append([I.add(k(a), I.mul(k(b), d)) for a, b in row['normal']])
    covectors = [I.mv(Q, r) for r in core]
    assert I.pair(center, Q, center)[1] < 0
    assert all(I.dot(center, a)[1] < 0 for a in covectors)
    polar = []
    details = []
    for row in recipes['polar_recipes']:
        ids = row['nine_walls']
        assert len(ids) == len(set(ids)) == 9
        free, transform, direction = I.solve_line([covectors[i] for i in ids])
        # solve_line proves exact rank nine and defines its exact null line.
        orientation = I.pair(direction, Q, center)
        assert orientation[0] > 0 or orientation[1] < 0
        if orientation[0] > 0:
            direction = [I.neg(a) for a in direction]
        square = I.pair(direction, Q, direction)
        assert square[0] > 0, (row['scout_vertex'], 'not spacelike')
        divisor = I.sqrt(I.div(square, I.rat(2)))
        normal = [I.div(a, divisor) for a in direction]
        assert I.pair(normal, Q, center)[1] < 0
        polar.append(normal)
        details.append({'scout_vertex': row['scout_vertex'], 'nine_walls': ids,
                        'free_coordinate': free, 'direction_square_interval': I.enc(square)})
    n = len(polar)
    polar_cov = [I.mv(Q, r) for r in polar]
    witnesses = []
    cutters = []
    cut_masks = [0]*n
    for wi, row in enumerate(recipes['witnesses']):
        e = [I.rat(a) for a in row['integer_point']]
        square, orientation = I.pair(e, Q, e), I.pair(e, Q, center)
        assert square[0] > 0 and orientation[1] < 0
        old_pairs = [I.dot(a, e) for a in covectors]
        assert all(p[1] < 0 for p in old_pairs), (wi, 'old witness sign')
        values = [I.dot(a, e) for a in polar_cov]
        possible = [i for i, p in enumerate(values) if p[1] >= 0]
        mask = sum(1 << i for i in possible)
        cutters.append(mask)
        for i in possible:
            cut_masks[i] |= 1 << wi
        witnesses.append({'integer_point': row['integer_point'], 'square_interval': I.enc(square),
                          'center_pair_interval': I.enc(orientation),
                          'strict_old_pair_intervals': list(map(I.enc, old_pairs)),
                          'possible_cutters': possible,
                          'all_polar_pair_intervals': list(map(I.enc, values))})
    proof = angle_enclosures()
    gate = AngleGate(proof)
    ratio = I.S//PAIR_SCALE
    assert ratio*PAIR_SCALE == I.S
    all_candidates = 0
    for mask in cutters:
        all_candidates |= mask
    conflicts = {}
    compatibility = {}

    def allowed(i):
        if i in compatibility:
            return compatibility[i]
        good = 0
        for j in indexes(all_candidates):
            if i == j:
                continue
            key = tuple(sorted((i, j)))
            if key in conflicts:
                continue
            pair = I.dot(polar[i], polar_cov[j])
            reason = gate.rejection(pair[0]//ratio, I.ceildiv(pair[1], ratio))
            if reason:
                conflicts[key] = {'polar_indices': list(key), 'reason': reason,
                                  'pair_interval': I.enc(pair)}
            else:
                good |= 1 << j
        compatibility[i] = good
        return good

    nodes = []

    def dfs(available, missing, selected):
        if not missing:
            return selected
        options = None
        smallest = n+1
        chosen = None
        for wi in indexes(missing):
            candidates = available & cutters[wi]
            size = candidates.bit_count()
            if not size:
                nodes.append({'selected': selected, 'witness': wi, 'branches': []})
                return None
            if size < smallest:
                smallest, options, chosen = size, candidates, wi
        order = sorted(indexes(options), key=lambda i: (cut_masks[i] & missing).bit_count(), reverse=True)
        nodes.append({'selected': selected, 'witness': chosen, 'branches': order})
        assert len(nodes) < 100000, 'Stop rather than claim bounded failure exhaustive.'
        for i in order:
            solution = dfs(available & allowed(i), missing & ~cut_masks[i], selected+[i])
            if solution is not None:
                return solution
            available &= ~(1 << i)
        return None

    solution = dfs(all_candidates, (1 << len(witnesses))-1, [])
    assert solution is None, ('Conservative graph has a surviving cover', solution)
    dependencies = [SOURCE, RECIPES, Path(__file__),
                    ROOT/'scripts/fresh_h9_split_multiple_edge_escape_20261013.py',
                    ROOT/'scripts/fresh_h9_split_leaf_3034_safe_graph_20261013.py']
    out = {'status': 'PASS_EXACT_POLAR_BANK_INCOMPATIBILITY_COVER_WITH_OLD_WALL_DELETIONS',
           'sources': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in dependencies},
           'interval_scale': str(I.S), 'polar_count': n, 'polars': details,
           'witness_count': len(witnesses), 'witnesses': witnesses,
           'candidates_cutting_some_witness': all_candidates.bit_count(),
           'angle_intervals': proof, 'exact_conflicts': list(conflicts.values()),
           'exhaustive_cover_tree': nodes,
           'scope': 'No Coxeter-compatible subset of these literal polar recipes, together with any subset of the27old walls, has finite volume. Each recipe is the center-oriented norm-two vector on its exact nine-wall null line. This is not an exhaustive list of all spacelike extreme rays or arbitrary caps; different polar recipes, nonpolar normals and other seeds remain open.'}
    path = ROOT/'results'/f'{TAG}_verify.json'
    if verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], n, 'polars;', len(witnesses), 'witnesses;', len(nodes), 'nodes;',
          len(conflicts), 'conflicts;', round(time.time()-start, 2), 'seconds', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    run(ap.parse_args().verify)
