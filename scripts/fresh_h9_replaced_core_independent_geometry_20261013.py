#!/usr/bin/env python3
"""Second exact geometry replay for the selected 45-wall changed core.

Input extraction reuses source profile loading once. Subsequent replay imports
no construction source: its arithmetic is Golden plus a separate recursive
quadratic tower implementation. No numerical geometry is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction as Q
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden as G

ROOT = Path(__file__).resolve().parents[1]
SOURCE = 'data/enumerations/fresh_h9_replaced_core_20261013_round4_selection.json'
INPUT = ROOT / 'results/fresh_h9_replaced_core_independent_geometry_20261013_input.json'
OUT = ROOT / 'results/fresh_h9_replaced_core_independent_geometry_20261013.json'
Z, O = G(0), G(1)


def golden(raw):
    return G.from_json(raw) if isinstance(raw, dict) else G(*map(Q, raw))


@dataclass(frozen=True)
class Tower:
    base: Tower | None = None
    radicand: object = None

    def parse(self, raw):
        return golden(raw) if self.base is None else tuple(self.base.parse(x) for x in raw)

    def lift(self, x):
        return x if self.base is None else (self.base.lift(x), self.base.lift(Z))

    def add(self, x, y):
        return x+y if self.base is None else tuple(self.base.add(a, b) for a, b in zip(x, y))

    def neg(self, x):
        return -x if self.base is None else tuple(self.base.neg(a) for a in x)

    def sub(self, x, y):
        return self.add(x, self.neg(y))

    def mul(self, x, y):
        if self.base is None:
            return x*y
        a, b = x
        c, d = y
        f = self.base
        return (f.add(f.mul(a, c), f.mul(self.radicand, f.mul(b, d))),
                f.add(f.mul(a, d), f.mul(b, c)))

    def sign(self, x):
        if self.base is None:
            return x.sign()
        a, b = x
        f = self.base
        sa, sb = f.sign(a), f.sign(b)
        if f.sign(self.radicand) == 0:
            return sa
        if not sb:
            return sa
        if not sa or sa == sb:
            return sb
        delta = f.sign(f.sub(f.mul(a, a), f.mul(self.radicand, f.mul(b, b))))
        # Zero is possible for dependent radicals; no independence assumption.
        return 0 if delta == 0 else sa if delta > 0 else sb

    def dot(self, x, y):
        value = self.lift(Z)
        for a, b in zip(x, y, strict=True):
            value = self.add(value, self.mul(a, b))
        return value

    def mv(self, matrix, x):
        return tuple(self.dot(row, x) for row in matrix)


K = Tower()


def field(raw):
    if raw['kind'] == 'golden':
        return K
    base = field(raw['base'])
    rad = base.parse(raw['radicand'])
    # A zero auxiliary radicand is a permitted redundant presentation.
    assert base.sign(rad) >= 0
    return Tower(base, rad)


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def extract():
    # Loading only. Every arithmetic claim below is replayed independently.
    from scripts.fresh_h9_replaced_core_20261013_base import profiles, cap_profile, read
    from scripts.fresh_h9_quartic_completion_20261013_arithmetic import K as SourceK

    def field_encode(f):
        return {'kind': 'golden'} if f is SourceK else {
            'kind': 'quadratic', 'base': field_encode(f.base),
            'radicand': f.base.encode(f.d)}

    data = read(SOURCE)
    case = data['cases'][0]
    cap_data = read(data['source'])
    gram, inv, h, profiles_by_name = profiles()
    sources = {SOURCE, data['source']}
    for i, ref in enumerate(data['required_caps']):
        sources.add(ref['source'])
        profiles_by_name[f'required{i}'] = cap_profile(read(ref['source'])['caps'][ref['index']], f'required{i}')
    for i in case['selected_caps']:
        profiles_by_name[f'cap{i}'] = cap_profile(cap_data['caps'][i], f'cap{i}')
    roots = []
    for name in case['walls']:
        row = profiles_by_name[name]
        f = row['field']
        roots.append({'name': name, 'field': field_encode(f),
                      'coefficients': list(map(f.encode, row['c'])),
                      'seed_pairings': list(map(f.encode, row['p'])),
                      'norm': SourceK.encode(row['q'])})
    bundle = {'source': SOURCE, 'case_index': 0,
              'source_sha256': {s: sha(s) for s in sorted(sources)},
              'gram': [[SourceK.encode(x) for x in row] for row in gram],
              'inverse': [[SourceK.encode(x) for x in row] for row in inv],
              'positive_center_covector': list(map(SourceK.encode, h)),
              'roots': roots,
              'first_escape': case['geometry']['escapes'][0]['certificate']}
    INPUT.write_text(json.dumps(bundle, indent=2)+'\n')


def inverse(matrix):
    n = len(matrix)
    rows = [list(row)+[G(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    pos = neg = 0
    # Exact inertia by symmetric 1x1 elimination. A zero diagonal can require
    # a 2x2 step in general; here a nonzero pivot is always available.
    current = [list(row) for row in matrix]
    while current:
        k = next(i for i in range(len(current)) if current[i][i])
        current[0], current[k] = current[k], current[0]
        for row in current:
            row[0], row[k] = row[k], row[0]
        pivot = current[0][0]
        pos += pivot.sign() > 0
        neg += pivot.sign() < 0
        current = [[current[i][j]-current[i][0]*current[0][j]/pivot
                    for j in range(1, len(current))] for i in range(1, len(current))]
    for j in range(n):
        k = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[k] = rows[k], rows[j]
        p = rows[j][j]
        rows[j] = [v/p for v in rows[j]]
        for i in range(n):
            if i != j:
                p = rows[i][j]
                rows[i] = [x-p*y for x, y in zip(rows[i], rows[j])]
    return tuple(tuple(row[n:]) for row in rows), (pos, neg, 0)


def run():
    # Regressions for the reducible/zero auxiliary radical presentations.
    assert Tower(K, Z).sign((G(-1), G(100))) == -1
    assert Tower(K, O).sign((O, -O)) == 0
    assert Tower(K, O).sign((G(2), -O)) == 1
    data = json.loads(INPUT.read_text())
    assert all(sha(s) == digest for s, digest in data['source_sha256'].items())
    source = json.loads((ROOT/SOURCE).read_text())
    case = source['cases'][data['case_index']]
    assert case['facet_count'] == 45
    assert [r['name'] for r in data['roots']] == case['walls']
    assert data['first_escape'] == case['geometry']['escapes'][0]['certificate']
    gram = tuple(tuple(map(golden, row)) for row in data['gram'])
    inv, signature = inverse(gram)
    assert signature == (9, 1, 0)
    assert inv == tuple(tuple(map(golden, row)) for row in data['inverse'])
    h = tuple(map(golden, data['positive_center_covector']))
    assert all(x.sign() > 0 for x in h)
    center_square = K.dot(h, K.mv(inv, h))
    assert center_square.sign() < 0
    cert = data['first_escape']
    y = tuple(map(golden, cert['y']))
    assert all(x.sign() > 0 for x in y)
    escape_square = K.dot(y, K.mv(inv, y))
    center_pairing = K.dot(y, K.mv(inv, h))
    assert escape_square == golden(cert['norm']) and escape_square.sign() > 0
    assert center_pairing == golden(cert['center_pairing']) and center_pairing.sign() < 0
    signs, center_signs, depths = [], [], []
    for row in data['roots']:
        f = field(row['field'])
        c, p = tuple(map(f.parse, row['coefficients'])), tuple(map(f.parse, row['seed_pairings']))
        lifted_gram = tuple(tuple(f.lift(z) for z in r) for r in gram)
        actual_p = f.mv(lifted_gram, c)
        assert all(f.sign(f.sub(x, z)) == 0 for x, z in zip(actual_p, p, strict=True))
        q = golden(row['norm'])
        assert q.sign() > 0
        assert f.sign(f.sub(f.dot(c, p), f.lift(q))) == 0
        center_signs.append(f.sign(f.dot(c, tuple(f.lift(z) for z in h))))
        signs.append(f.sign(f.dot(c, tuple(f.lift(z) for z in y))))
        depth, current = 0, f
        while current.base is not None:
            depth += 1
            current = current.base
        depths.append(depth)
    assert center_signs == signs == [1]*45 == cert['all_inequality_signs']
    # The source's independently replayed pair census is the additional
    # hypothesis for the facet projection argument; this audit does not
    # repeat compositum-valued pair classification.
    assert len(case['pair_checks']) == 990
    assert all(r['kind'] in ('finite', 'dotted') for r in case['pair_checks'].values())
    assert not case['affine_census']['affine_circuits']
    bad = [Z]*10
    assert K.dot(bad, K.mv(inv, bad)).sign() != 1
    return {'status': 'PASS_INDEPENDENT_45_WALL_STRICT_GEOMETRY',
            'source': SOURCE, 'case_index': 0, 'wall_count': 45,
            'gram_signature': list(signature), 'root_norm_and_gram_checks': 45,
            'common_strict_center_checks': 45, 'strict_escape_checks': 45,
            'maximum_quadratic_tower_depth': max(depths),
            'center_square': center_square.to_json(),
            'escape_square': escape_square.to_json(),
            'center_escape_pairing': center_pairing.to_json(),
            'zero_escape_rejected': True,
            'dependent_and_zero_radical_sign_regressions': True,
            'scope': ('Exact norm, Gram, common interior, and first strict future spacelike '
                      'escape replay with separate arithmetic. Input extraction reuses the '
                      'source loader. The all-pair Coxeter gate and complete affine census '
                      'are independently replayed by the construction agent, not repeated '
                      'here. The displayed 45-wall chamber has infinite volume; arbitrary '
                      'further caps are not excluded.')}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--extract', action='store_true')
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    if args.extract:
        extract()
    result = run()
    if args.verify:
        assert json.loads(OUT.read_text()) == result
    else:
        OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(result)


if __name__ == '__main__':
    main()
