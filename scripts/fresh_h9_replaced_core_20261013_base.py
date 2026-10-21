#!/usr/bin/env python3
"""Exact common-tower helpers and source-anchored profiles for changed cores."""
from functools import lru_cache
from fractions import Fraction as Q
from scripts.fresh_h9_quartic_completion_20261013_opposite import build as source_build
from scripts.fresh_h9_quartic_completion_20261013_polar import (
    ROOT, K, Quadratic, Z, O, P, parse, enc, read, scale, add, sub, mul, div,
    dot, matvec, pairing, lift_K)
from scripts.fresh_h9_quartic_completion_20261013_verify import angle_pool

ANGLE_SOURCE = 'data/enumerations/fresh_h9_quartic_completion_20261013_angles.json'
BASE_SOURCE = 'data/enumerations/fresh_h8_alternate_3_6_20261013.json'


@lru_cache(None)
def key(f):
    return ('K',) if f is K else (key(f.base), f.d)


def embed(target, source):
    """Exact source embedding, reusing structurally identical tower prefixes."""
    if key(target) == key(source):
        return target, lambda x: x
    if target is not K:
        ancestor = target.base
        while True:
            if key(ancestor) == key(source):
                def lift(x):
                    chain = []
                    t = target
                    while key(t) != key(source):
                        chain.append(t)
                        t = t.base
                    for t in reversed(chain):
                        x = t.lift(x)
                    return x
                return target, lift
            if ancestor is K:
                break
            ancestor = ancestor.base
    if source is K:
        return target, lambda x: lift_K(target, x)
    middle, source_base = embed(target, source.base)
    common = Quadratic(middle, source_base(source.d))
    return common, lambda x: (source_base(x[0]), source_base(x[1]))


def lift_existing(target, source, x):
    if target is source:
        return x
    return target.lift(lift_existing(target.base, source, x))


def pair_profiles(a, b):
    f, ca, pa, qa = a['field'], a['c'], a['p'], a['q']
    g = b['field']
    common, eb = embed(f, g)
    value = dot(common, tuple(lift_existing(common, f, x) for x in pa),
                tuple(eb(x) for x in b['c']))
    return common, value


def relation(a, b, pool):
    from scripts.fresh_h9_quartic_completion_20261013_arithmetic import finite_label
    f, p = pair_profiles(a, b)
    if f.sign(p) > 0:
        return {'kind': 'positive'}
    square = f.mul(f.mul(p, p), lift_K(f, div(scale(O, 4), mul(a['q'], b['q']))))
    cmp = f.sign(f.sub(square, f.rational(4)))
    if cmp >= 0:
        return {'kind': 'dotted' if cmp > 0 else 'parallel'}
    label = finite_label(f, square, pool)
    return {'kind': 'finite' if label is not None else 'unclassified', 'label': label}


def profiles():
    record, ctx = source_build()
    G, I, h, old, case, f, p, c, prior = ctx
    out = {}
    for i in range(10):
        out[f'seed{i}'] = {'name': f'seed{i}', 'field': K,
            'c': tuple(O if i == j else Z for j in range(10)),
            'p': tuple(G[j][i] for j in range(10)), 'q': scale(O, 2),
            'reference': {'kind': 'seed', 'index': i}}
    for i in case['retained_old_caps']:
        n = old[i]
        if all(t == Z for t in n['b']):
            ff, cc, pp = K, n['a'], n['pa']
        else:
            ff, cc, pp = Quadratic(K, n['d']), tuple(zip(n['a'], n['b'])), tuple(zip(n['pa'], n['pb']))
        out[f'old{i}'] = {'name': f'old{i}', 'field': ff, 'c': cc, 'p': pp,
                         'q': n['q'], 'reference': {'kind': 'old', 'index': i}}
    for j, (ref, ff, cc, pp) in enumerate(prior):
        out[f'new{j}'] = {'name': f'new{j}', 'field': ff, 'c': cc, 'p': pp,
                         'q': scale(O, 2), 'reference': ref}
    out['opposite'] = {'name': 'opposite', 'field': f, 'c': c, 'p': p,
                      'q': scale(O, 2), 'reference': {'kind': 'opposite'}}
    assert len(out) == 46
    return G, I, h, out


def value(f, x):
    """High-precision conversion; used only by proposal scripts."""
    import mpmath as mp
    if f is K:
        return (mp.mpf(x[0].numerator)/x[0].denominator+
                mp.mpf(x[1].numerator)/x[1].denominator*(1+mp.sqrt(5))/2)
    return value(f.base, x[0])+mp.sqrt(value(f.base, f.d))*value(f.base, x[1])


def strict_escape(y, I, center, selected):
    assert all(K.sign(t) > 0 for t in y)
    norm, cp = pairing(K, y, I, y), pairing(K, y, I, center)
    assert K.sign(norm) > 0 and K.sign(cp) < 0
    signs = []
    for row in selected:
        f = row['field']
        signs.append(f.sign(dot(f, row['c'], tuple(lift_K(f, t) for t in y))))
    assert all(s > 0 for s in signs)
    return {'y': list(map(enc, y)), 'norm': enc(norm),
            'center_pairing': enc(cp), 'all_inequality_signs': signs}


def cap_profile(record, name):
    e = Quadratic(K, parse(record['base_radicand']))
    f = Quadratic(e, e.parse(record['discriminant']))
    return {'name': name, 'field': f, 'c': tuple(map(f.parse, record['seed_coefficients'])),
            'p': tuple(map(f.parse, record['seed_pairings'])), 'q': scale(O, 2)}


def add_required(profiles, core, references):
    names = []
    for i, ref in enumerate(references):
        name = f'required{i}'
        profiles[name] = cap_profile(read(ref['source'])['caps'][ref['index']], name)
        names.append(name)
    for case in core['cases']:
        case['walls'] += names
    return names


def collapsed_anchor(reference):
    """Canonical K(sqrt D) presentation when the saved first radical is sqrt1."""
    record = read(reference['source'])['caps'][reference['index']]
    assert parse(record['base_radicand']) == O
    formal = Quadratic(K, O)
    collapse = lambda x: add(x[0], x[1])
    d = collapse(formal.parse(record['discriminant']))
    e = Quadratic(K, d)
    source = Quadratic(formal, formal.parse(record['discriminant']))
    coefficients = tuple(tuple(collapse(y) for y in source.parse(x))
                         for x in record['seed_coefficients'])
    return e, coefficients
