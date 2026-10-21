#!/usr/bin/env python3
"""Exact comparison with the preceding banks, using normalized seed pairings."""
import argparse
from collections import defaultdict
import json
from scripts.fresh_h9_quartic_completion_20261013_verify import (
    ROOT, read, matrix, load_old, node, angle_pool)
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import (
    K, Quadratic, Z, O, P, parse, add, div, scale, dot, finite_label)

DEFAULT_SOURCES = [
    'data/enumerations/fresh_h9_quartic_completion_20261013_targeted.json',
    'data/enumerations/fresh_h9_quartic_completion_20261013_adaptive.json',
    'data/enumerations/fresh_h9_quartic_completion_20261013_second.json']
OUTPUT = ROOT/'data/enumerations/fresh_h9_quartic_completion_20261013_mirrors.json'


def run(sources, verify=False):
    first = read(sources[0])
    G = matrix(read(first['base_source'])['seed_gram'])
    E = Quadratic(K, parse(first['anchor_radicand']))
    polynomials = angle_pool(first['angle_source'])
    old = load_old(G)
    old += [node(c, G) for c in read('data/enumerations/fresh_h9_joint_completion_20261013.json')['caps']]
    groups = defaultdict(list)
    comparisons, aliases, novel = 0, [], []

    def pattern(field, squares):
        labels = []
        easy = [(2, K.zero), (3, K.one), (4, scale(O, 2)),
                (5, add(O, P)), (6, scale(O, 3))]
        for square in squares:
            cmp = field.sign(field.sub(square, field.rational(4)))
            if cmp >= 0:
                labels.append('dotted' if cmp else 'parallel')
                continue
            label = None
            for candidate, value in easy:
                lift = field.lift(value) if field.base is K else field.lift(E.lift(value))
                if field.sign(field.sub(square, lift)) == 0:
                    label = candidate
                    break
            if label is None:
                label = finite_label(field, square, polynomials)
            labels.append(label)
        return tuple(labels)

    for i, n in enumerate(old):
        f = Quadratic(K, n['d'])
        ps, cs = tuple(zip(n['pa'], n['pb'])), tuple(zip(n['a'], n['b']))
        assert all(f.sign(p) <= 0 for p in ps)
        assert f.sign(f.sub(dot(f, ps, cs), f.lift(n['q']))) == 0
        factor = f.lift(div(scale(O, 2), n['q']))
        squares = tuple(f.mul(f.mul(p, p), factor) for p in ps)
        groups[pattern(f, squares)].append(
            {'kind': 'old', 'd': n['d'], 'squares': squares,
             'source': 'preceding-banks', 'index': i})
    for path in sources:
        data = read(path)
        assert parse(data['anchor_radicand']) == E.d
        for index, c in enumerate(data['caps']):
            f = Quadratic(E, E.parse(c['discriminant_over_E']))
            ps = tuple(map(f.parse, c['seed_pairings']))
            cs = tuple(map(f.parse, c['seed_coefficients']))
            assert all(f.sign(p) <= 0 for p in ps)
            assert f.sign(f.sub(dot(f, ps, cs), f.rational(2))) == 0
            squares = tuple(f.mul(p, p) for p in ps)
            key = pattern(f, squares)
            match = None
            for previous in groups[key]:
                comparisons += 1
                if previous['kind'] == 'old':
                    lift = lambda x: f.lift(E.lift(x))
                else:
                    lift = f.lift
                extension = Quadratic(f, lift(previous['d']))
                if all(extension.sign((f.sub(a, lift(b[0])), f.neg(lift(b[1])))) == 0
                       for a, b in zip(squares, previous['squares'], strict=True)):
                    match = {k: previous[k] for k in ('source', 'index')}
                    break
            origin = {'source': path, 'index': index}
            if match is not None:
                aliases.append({'candidate': origin, 'equals': match})
            else:
                novel.append(origin)
                groups[key].append({'kind': 'new', 'd': f.d, 'squares': squares, **origin})
    result = {'status': 'PASS_EXACT_MIRROR_COMPARISON',
              'sources': sources, 'preceding_bank_entries': len(old),
              'proposal_entries': len(novel)+len(aliases),
              'distinct_new_mirrors': len(novel), 'new_representatives': novel,
              'aliases': aliases, 'exact_same_pattern_comparisons': comparisons,
              'scope': 'Equality of all normalized squared pairings with the independent '
                       'ten-wall seed, using nonpositive orientations. Exact novelty '
                       'relative to the preceding saved banks, not worldwide priority.'}
    if verify:
        assert json.loads(OUTPUT.read_text()) == result
    else:
        OUTPUT.write_text(json.dumps(result, indent=2)+'\n')
    print({k: v for k, v in result.items() if k not in ('new_representatives', 'aliases')}, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--sources', nargs='+', default=DEFAULT_SOURCES)
    p.add_argument('--verify', action='store_true')
    p.add_argument('--output')
    a = p.parse_args()
    if a.output:
        OUTPUT = ROOT/a.output
    run(a.sources, a.verify)
