#!/usr/bin/env python3
"""Authoritative rebuilt45 source and exact fresh-anchor arithmetic."""
from scripts.fresh_h9_replaced_core_20261013_base import *

SOURCE = 'data/enumerations/fresh_h9_replaced_core_20261013_round4_selection.json'
ANCHOR = {'source': 'data/enumerations/fresh_h9_replaced_core_20261013_round4_caps.json',
          'index': 2}
PREFIX = 'data/enumerations/fresh_h9_gap_completion_20261013'


def current():
    G, I, h, ps = profiles()
    data = read(SOURCE)
    case = data['cases'][0]
    bank = read(data['source'])
    for i, ref in enumerate(data['required_caps']):
        ps[f'required{i}'] = cap_profile(read(ref['source'])['caps'][ref['index']],
                                       f'required{i}')
    for i in case['selected_caps']:
        ps[f'cap{i}'] = cap_profile(bank['caps'][i], f'cap{i}')
    names = case['walls']
    assert len(names) == len(set(names)) == 45
    return G, I, h, {n: ps[n] for n in names}


def configuration(name='45'):
    if name == '48':
        G, I, h, ps, _, _ = configuration('swapped47')
        source = PREFIX+'_left_geometry.json'
        data = read(source)
        case = data['cases'][0]
        bank = read(data['source'])
        for i in case['selected_caps']:
            ps[f'add{i}'] = cap_profile(bank['caps'][i], f'add{i}')
        assert list(ps) == case['walls'] and len(ps) == 48
        return G, I, h, ps, source, 0
    G, I, h, ps = current()
    if name == '45':
        return G, I, h, ps, SOURCE, 0
    assert name in ('47', 'swapped47')
    source = PREFIX+'_geometry.json'
    data = read(source)
    bank = read(data['source'])
    case = data['cases'][0]
    for i in case['selected_caps']:
        ps[f'gap{i}'] = cap_profile(bank['caps'][i], f'gap{i}')
    assert list(ps) == case['walls']
    if name == '47':
        return G, I, h, ps, source, 0
    source = PREFIX+'_swap.json'
    case = read(source)['cases'][1]
    assert case['transported_source'] == 'new6' and case['dropped'] == ['new6']
    original = ps.pop('new6')
    perm = case['permutation']
    ps['swap_new6'] = dict(original, name='swap_new6',
        c=tuple(original['c'][i] for i in perm), p=tuple(original['p'][i] for i in perm))
    assert list(ps) == case['walls']
    return G, I, h, ps, source, 1


def canonical_anchor(reference):
    record = read(reference['source'])['caps'][reference['index']]
    r = parse(record['base_radicand'])
    assert r in (Z, O)
    formal = Quadratic(K, r)
    collapse = (lambda x: x[0]) if r == Z else (lambda x: add(x[0], x[1]))
    d = collapse(formal.parse(record['discriminant']))
    e = Quadratic(K, d)
    source = Quadratic(formal, formal.parse(record['discriminant']))
    coefficients = tuple(tuple(collapse(y) for y in source.parse(x))
                         for x in record['seed_coefficients'])
    return e, coefficients


def selected_coupled_case(index=0):
    """Pure loader for the final73 normals and saved vertex/escape proofs."""
    G, I, h, ps, _, _ = configuration('48')
    geometry = read(PREFIX+'_coupled_geometry.json')
    case = geometry['cases'][index]
    vertices = read(geometry['vertex_cap_source'])
    coupled = read(geometry['source'])
    for i in case['selected_generic']:
        record = vertices['caps'][i]
        f = Quadratic(K, parse(record['radicand']))
        name = f'vertex{i}'
        ps[name] = {'name': name, 'field': f, 'q': scale(O, 2),
                    'p': tuple((Z, parse(t)) for t in record['raw_seed_pairings']),
                    'c': tuple((Z, parse(t)) for t in record['raw_seed_coefficients'])}
    for i in case['selected_coupled']:
        name = f'coupled{i}'
        ps[name] = cap_profile(coupled['caps'][i], name)
    assert list(ps) == case['walls'] and len(ps) == case['facet_count'] == 73
    return G, I, h, ps, case


def exact_cap(spec, G, I, h, ps, pool, anchor_ref=ANCHOR):
    u, v, w, fixed_indices, fixed_labels, anchor_label, polar_label, branch = spec
    E, a = canonical_anchor(anchor_ref)
    b = tuple(E.lift(t) for t in ps['old348']['c'])
    M = E.lift({2: Z, 3: O, 5: P}[anchor_label])
    N = E.lift(Z if polar_label == 2 else O)
    fixed = [E.zero]*10
    for j, label in zip(fixed_indices, fixed_labels, strict=True):
        fixed[j] = E.lift(O if label == 3 else P)
    determinant = E.sub(E.mul(a[v], b[w]), E.mul(a[w], b[v]))
    assert E.sign(determinant)
    r0, r1 = E.sub(M, dot(E, a, fixed)), E.sub(N, dot(E, b, fixed))
    line0, line1 = list(fixed), [E.zero]*10
    line0[v] = E.div(E.sub(E.mul(r0, b[w]), E.mul(a[w], r1)), determinant)
    line0[w] = E.div(E.sub(E.mul(a[v], r1), E.mul(r0, b[v])), determinant)
    line1[u] = E.one
    line1[v] = E.div(E.sub(E.mul(a[w], b[u]), E.mul(a[u], b[w])), determinant)
    line1[w] = E.div(E.sub(E.mul(a[u], b[v]), E.mul(a[v], b[u])), determinant)
    IE = tuple(tuple(E.lift(t) for t in row) for row in I)
    A = pairing(E, line1, IE, line1)
    B = E.mul(E.rational(2), pairing(E, line0, IE, line1))
    C = E.sub(pairing(E, line0, IE, line0), E.rational(2))
    D = E.sub(E.mul(B, B), E.mul(E.rational(4), E.mul(A, C)))
    assert E.sign(D) > 0 and E.sign(A)
    F = Quadratic(E, D)
    x = (E.div(E.neg(B), E.mul(E.rational(2), A)),
         E.div(E.rational(branch), E.mul(E.rational(2), A)))
    pa = tuple(F.neg(F.add(F.lift(t), F.mul(x, F.lift(s))))
               for t, s in zip(line0, line1, strict=True))
    cf = matvec(F, tuple(tuple(F.lift(t) for t in row) for row in IE), pa)
    assert F.sign(F.sub(dot(F, cf, pa), F.rational(2))) == 0
    assert all(F.sign(t) <= 0 for t in pa)
    assert F.sign(dot(F, cf, tuple(lift_K(F, t) for t in h))) > 0
    assert F.sign(F.add(dot(F, pa, tuple(F.lift(t) for t in a)), F.lift(M))) == 0
    assert F.sign(F.add(dot(F, pa, tuple(F.lift(t) for t in b)), F.lift(N))) == 0
    row = {'field': F, 'c': cf, 'p': pa, 'q': scale(O, 2)}
    seedchecks = {str(i): relation(row, ps[f'seed{i}'], pool) for i in range(10)}
    assert all(r['kind'] in ('finite', 'dotted') for r in seedchecks.values())
    return row, {'kind': 'new_quartic_anchor', 'anchor_reference': anchor_ref,
        'anchor_label': anchor_label, 'free_seed_indices': [u, v, w],
        'fixed_seed_indices': list(fixed_indices), 'fixed_seed_labels': list(fixed_labels),
        'polar_label': polar_label, 'branch': branch, 'base_radicand': enc(E.d),
        'line_constant': list(map(E.encode, line0)), 'line_direction': list(map(E.encode, line1)),
        'norm_polynomial': list(map(E.encode, (A, B, C))), 'discriminant': E.encode(D),
        'parameter': F.encode(x), 'seed_pairings': list(map(F.encode, pa)),
        'seed_coefficients': list(map(F.encode, cf)), 'seed_relations': seedchecks}
