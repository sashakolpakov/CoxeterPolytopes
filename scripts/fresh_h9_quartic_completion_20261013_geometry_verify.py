#!/usr/bin/env python3
"""Independent exact replay of the coherent tower-field chambers and escapes."""
import argparse
from itertools import combinations
import json
from pathlib import Path
from scripts.fresh_h9_quartic_completion_20261013_verify import (
    read, matrix, load_old, angle_pool)
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import (
    K, Quadratic, Z, O, parse, div, scale, dot, matvec, pairing,
    finite_label, old_pair_gate, new_pair_gate)

ROOT = Path(__file__).resolve().parents[1]


def inertia_K(matrix):
    """Exact symmetric elimination, including zero-diagonal 2x2 pivots."""
    a = [list(row) for row in matrix]
    positive = negative = zero = 0
    while a:
        pivot = next((i for i in range(len(a)) if K.sign(a[i][i])), None)
        if pivot is not None:
            order = [pivot]+[i for i in range(len(a)) if i != pivot]
            a = [[a[i][j] for j in order] for i in order]
            d = a[0][0]
            positive += K.sign(d) > 0
            negative += K.sign(d) < 0
            a = [[K.sub(a[i][j], K.div(K.mul(a[i][0], a[0][j]), d))
                  for j in range(1, len(a))] for i in range(1, len(a))]
            continue
        pair = next(((i, j) for i, j in combinations(range(len(a)), 2)
                     if K.sign(a[i][j])), None)
        if pair is None:
            zero += len(a)
            break
        order = list(pair)+[i for i in range(len(a)) if i not in pair]
        a = [[a[i][j] for j in order] for i in order]
        b = a[0][1]
        positive += 1
        negative += 1
        a = [[K.sub(a[i][j], K.div(K.add(K.mul(a[i][0], a[1][j]),
                                                 K.mul(a[i][1], a[0][j])), b))
              for j in range(2, len(a))] for i in range(2, len(a))]
    return positive, negative, zero


def old_gate(a, b, polynomials):
    f = Quadratic(K, a['d'])
    ext = Quadratic(f, f.lift(b['d']))
    pa = tuple(zip(a['pa'], a['pb']))
    value = (dot(f, pa, tuple(f.lift(x) for x in b['a'])),
             dot(f, pa, tuple(f.lift(x) for x in b['b'])))
    if ext.sign(value) > 0:
        return {'compatible': False, 'kind': 'positive', 'label': None}
    factor = ext.lift(f.lift(div(scale(O, 4), K.mul(a['q'], b['q']))))
    square = ext.mul(ext.mul(value, value), factor)
    cmp = ext.sign(ext.sub(square, ext.rational(4)))
    if cmp > 0:
        return {'compatible': True, 'kind': 'dotted', 'label': None}
    if cmp == 0:
        return {'compatible': False, 'kind': 'parallel', 'label': None}
    label = finite_label(ext, square, polynomials)
    return {'compatible': label is not None,
            'kind': 'finite' if label is not None else 'noncoxeter', 'label': label}


def check(source, data=None):
    data = read(source) if data is None else data
    capdata = read(data['source'])
    base = read(capdata['base_source'])
    G, I = matrix(base['seed_gram']), matrix(base['seed_gram_inverse'])
    center = tuple(map(parse, base['strict_center_y']))
    assert K.sign(pairing(K, center, I, center)) < 0
    assert all(K.sign(t) > 0 for t in center)
    assert all(matvec(K, I, row) == tuple(O if i == j else Z for j in range(10))
               for i, row in enumerate(G))
    assert inertia_K(G) == (9, 1, 0)
    old = load_old(G)
    polynomials = angle_pool(capdata['angle_source'])
    assert all(G[i][i] == scale(O, 2) for i in range(10))
    for i, j in combinations(range(10), 2):
        p = G[i][j]
        assert K.sign(p) <= 0
        square = K.mul(p, p)
        cmp = K.sign(K.sub(square, scale(O, 4)))
        assert cmp > 0 or (cmp < 0 and finite_label(K, square, polynomials) is not None)
    E = Quadratic(K, parse(capdata['anchor_radicand']))
    refs = data.get('cap_references', [{'source': data['source'], 'index': i}
                                      for i in range(len(capdata['caps']))])
    cache = {data['source']: capdata}
    fields, coefficients, pa = [], [], []
    for ref in refs:
        cache.setdefault(ref['source'], read(ref['source']))
        src = cache[ref['source']]
        assert parse(src['anchor_radicand']) == E.d
        c = src['caps'][ref['index']]
        f = Quadratic(E, E.parse(c['discriminant_over_E']))
        cf, pf = tuple(map(f.parse, c['seed_coefficients'])), tuple(map(f.parse, c['seed_pairings']))
        gf = tuple(tuple(f.lift(E.lift(t)) for t in row) for row in G)
        assert all(f.sign(f.sub(x, y)) == 0 for x, y in zip(pf, matvec(f, gf, cf), strict=True))
        assert f.sign(f.sub(dot(f, cf, pf), f.rational(2))) == 0
        assert all(f.sign(x) <= 0 for x in pf)
        assert f.sign(dot(f, cf, tuple(f.lift(E.lift(t)) for t in center))) > 0
        for x in pf:
            square = f.mul(x, x)
            cmp = f.sign(f.sub(square, f.rational(4)))
            assert cmp > 0 or (cmp < 0 and finite_label(f, square, polynomials) is not None)
        fields.append(f)
        coefficients.append(cf)
        pa.append(pf)
    newchecks = {}
    for i, j in combinations(range(len(fields)), 2):
        newchecks[f'{i},{j}'] = new_pair_gate(fields[i], pa[i], fields[j].d,
                                            coefficients[j], polynomials)
    assert newchecks == data['all_new_pair_checks']
    counts = {'new_pair_checks': len(newchecks), 'configurations': 0,
              'configuration_old_new_checks': 0, 'strict_escapes': 0,
              'new_dictionary_witness_signs': 0, 'old_dictionary_witness_signs': 0}
    facets = []
    for case in data['cases']:
        ids, oldids = case['selected_new_caps'], case['retained_old_caps']
        assert len(set(ids)) == len(ids) and len(set(oldids)) == len(oldids)
        for i, j in combinations(ids, 2):
            assert newchecks[f'{min(i,j)},{max(i,j)}']['compatible']
        oldchecks = {f'{i},{j}': old_gate(old[i], old[j], polynomials)
                     for i, j in combinations(oldids, 2)}
        assert oldchecks == case['old_pair_checks']
        assert all(r['compatible'] for r in oldchecks.values())
        for i in ids:
            for j in oldids:
                assert old_pair_gate(fields[i], pa[i], old[j], polynomials)['compatible']
                counts['configuration_old_new_checks'] += 1
        for j in oldids:
            n = old[j]
            f = Quadratic(K, n['d'])
            cc, pp = tuple(zip(n['a'], n['b'])), tuple(zip(n['pa'], n['pb']))
            assert f.sign(f.sub(dot(f, cc, pp), f.lift(n['q']))) == 0
            assert all(f.sign(t) <= 0 for t in pp)
            assert f.sign(dot(f, cc, tuple(f.lift(t) for t in center))) > 0
            factor = f.lift(div(scale(O, 2), n['q']))
            for p in pp:
                square = f.mul(f.mul(p, p), factor)
                cmp = f.sign(f.sub(square, f.rational(4)))
                # Some inherited old caps have parallel seed relations.
                # This is allowed for the finite-volume target. The optional
                # new-cap graph separately rejects new parallel pairs.
                assert cmp >= 0 or finite_label(f, square, polynomials) is not None
        assert case['facet_count'] == 10+len(ids)+len(oldids)
        facets.append(case['facet_count'])
        assert case['escapes'] and case['status'] == 'EXACT_ESCAPES'
        for cert in case['escapes']:
            y = tuple(map(parse, cert['y']))
            assert len(y) == 10 and all(K.sign(t) > 0 for t in y)
            norm, cp = pairing(K, y, I, y), pairing(K, y, I, center)
            assert norm == parse(cert['norm']) and K.sign(norm) > 0
            assert cp == parse(cert['center_pairing']) and K.sign(cp) < 0
            signs = [f.sign(dot(f, cc, tuple(f.lift(E.lift(t)) for t in y)))
                     for f, cc in zip(fields, coefficients)]
            oldsigns = [Quadratic(K, n['d']).sign((dot(K, n['a'], y), dot(K, n['b'], y)))
                        for n in old]
            assert all(signs[i] > 0 for i in ids)
            assert all(oldsigns[j] > 0 for j in oldids)
            assert [i for i, s in enumerate(signs) if s < 0] == cert['new_dictionary_cutters']
            assert [i for i, s in enumerate(oldsigns) if s < 0] == cert['old_dictionary_cutters']
            counts['strict_escapes'] += 1
            counts['new_dictionary_witness_signs'] += len(fields)
            counts['old_dictionary_witness_signs'] += len(old)
        counts['configurations'] += 1
    return {'status': 'PASS_EXACT_COHERENT_CHAMBERS_AND_ESCAPES',
            'source': source, **counts, 'facet_counts': facets,
            'scope': 'All saved pairwise Coxeter relations, strict common centers and '
                     'strict future spacelike escapes. Each specified chamber has infinite '
                     'volume; no global completion or whole-bank obstruction.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', default='data/enumerations/fresh_h9_quartic_completion_20261013_geometry.json')
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    result = check(a.source)
    output = ROOT/a.source.replace('.json', '_verify.json')
    if a.verify:
        assert json.loads(output.read_text()) == result
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result, flush=True)
