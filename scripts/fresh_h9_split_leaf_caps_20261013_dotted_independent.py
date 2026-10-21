#!/usr/bin/env python3
"""Independent Fraction arithmetic and NetworkX audit of the fixed 228 caps.

No repository arithmetic, constructor, graph routine, or search code is imported.
The only nonstdlib dependency is NetworkX for a separate maximal-clique census.
"""
import argparse
import copy
import hashlib
import json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/enumerations/fresh_h9_split_leaf_caps_20261013_exact.json'
CERT = ROOT / 'data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_exact.json'
OUT = ROOT / 'results/fresh_h9_split_leaf_caps_20261013_dotted_independent.json'
K0 = (F(0), F(0))
K1 = (F(1), F(0))
PHI = (F(0), F(1))
Z = (K0, K0)
O = (K1, K0)


def ka(x, y):
    return (x[0] + y[0], x[1] + y[1])


def kn(x):
    return (-x[0], -x[1])


def km(x, y):
    if x == K0 or y == K0:
        return K0
    a, b = x
    c, d = y
    return (a*c + b*d, a*d + b*c + b*d)


def ki(x):
    # phi -> 1-phi, with rational field norm.
    y = (x[0] + x[1], -x[1])
    n = km(x, y)
    assert n[1] == 0 and n[0] != 0
    return (y[0]/n[0], y[1]/n[0])


def sgn(q):
    return (q > 0) - (q < 0)


def ks(x):
    # 2(a+b phi)=(2a+b)+b sqrt(5).
    r, s = 2*x[0] + x[1], x[1]
    sr, ss = sgn(r), sgn(s)
    if not sr:
        return ss
    if not ss or sr == ss:
        return sr
    return sr * sgn(r*r - 5*s*s)


def add(x, y):
    return (ka(x[0], y[0]), ka(x[1], y[1]))


def neg(x):
    return (kn(x[0]), kn(x[1]))


def mul(x, y):
    if x == Z or y == Z:
        return Z
    u, v = x
    p, q = y
    return (ka(km(u, p), km(PHI, km(v, q))), ka(km(u, q), km(v, p)))


def rat(x):
    return ((F(x), F(0)), K0)


def scale(x, q):
    q = F(q)
    return ((q*x[0][0], q*x[0][1]), (q*x[1][0], q*x[1][1]))


def sign(x):
    u, v = x
    su, sv = ks(u), ks(v)
    if not su:
        return sv
    if not sv or su == sv:
        return su
    difference = ka(km(u, u), kn(km(PHI, km(v, v))))
    return su * ks(difference)


def total(xs):
    r = Z
    for x in xs:
        r = add(r, x)
    return r


def dot(x, y):
    assert len(x) == len(y)
    return total(mul(a, b) for a, b in zip(x, y))


def mv(g, x):
    return tuple(dot(row, x) for row in g)


def parse_k(x):
    if isinstance(x, dict):
        return (F(x['a']), F(x['b']))
    return tuple(map(F, x))


def parse(x):
    return tuple(map(parse_k, x))


def enc(x):
    return [{'a': str(a), 'b': str(b)} for a, b in x]


def digest(x):
    return hashlib.sha256(json.dumps(x, separators=(',', ':')).encode()).hexdigest()


def admissible(x):
    finite = x[1] == K0 and x[0] in (K0, (-F(1), F(0)), kn(PHI))
    return finite or sign(add(x, rat(2))) <= 0


def check_norm(g, row):
    assert dot(row, mv(g, row)) == rat(2), 'Cap norm is not two'


def witness(g, rows, raw):
    assert len(raw['coefficients']) == 16
    x = tuple(rat(v) for v in raw['coefficients'])
    gx = mv(g, x)
    square = dot(x, gx)
    center_pairing = total(gx)
    assert sign(square) > 0, 'Escape is not spacelike'
    assert sign(center_pairing) < 0, 'Escape has wrong time orientation'
    assert all(sign(v) < 0 for v in gx), 'Escape fails a seed wall'
    signs = [sign(dot(row, gx)) for row in rows]
    assert all(signs), 'Escape lies on a cap plane'
    cuts = [i for i, s in enumerate(signs) if s > 0]
    assert enc(square) == raw['square']
    assert enc(center_pairing) == raw['center_pairing']
    assert cuts == raw['cut_caps']
    return sum(1 << i for i in cuts)


def audit():
    source = json.loads(SOURCE.read_text())
    saved = json.loads(CERT.read_text())
    meta = saved['metadata']
    a = (F(1, 2), F(1))
    b = kn(km(a, ki(ka((F(1), F(0)), (12*a[0], 12*a[1])))))
    assert parse_k(source['a']) == a and parse_k(source['b']) == b
    assert b == (-F(5, 22), F(1, 11))
    assert ks(a) > 0 and ks(b) < 0 and ks(ka(b, (F(1, 2), F(0)))) > 0
    assert ks(ka(ka(a, b), (-F(2), F(0)))) > 0
    eps = (1,)*4 + (-1,)*4
    gram = []
    for i in range(16):
        row = []
        for j in range(16):
            if i//8 == j//8:
                p = (F(2 if i == j else 0), F(0))
            else:
                q = eps[i % 8]*eps[j % 8]
                p = kn(ka(ka(a, (q*b[0], q*b[1])), (F(2 if i % 8 == j % 8 else 0), F(0))))
            row.append((p, K0))
        gram.append(tuple(row))
    g = tuple(gram)
    assert g == tuple(tuple((parse_k(x), K0) for x in row) for row in source['gram'])
    assert all(admissible(g[i][j]) for i, j in combinations(range(16), 2))
    center = (O,)*16
    assert tuple((parse_k(x), K0) for x in source['strict_center']) == center
    gc = mv(g, center)
    assert gc == (scale((a, K0), -8),)*16
    center_square = dot(center, gc)
    assert sign(center_square) < 0
    # Sum/difference invariant decomposition: one negative, nine positive,
    # six radical modes. Values below are the nonzero signs of the blocks.
    assert ks(ka((F(8), F(0)), (16*b[0], 16*b[1]))) > 0
    assert source['seed_signature'] == [9, 1, 6]
    rows = [tuple((parse_k(x), K0) for x in row['normal']) for row in source['caps']]
    assert len(rows) == 192
    for row, old in zip(rows, source['caps']):
        cov = mv(g, row)
        expected = tuple((kn(parse_k(x)), K0) for x in old['x'] + old['y'])
        assert cov == expected
        assert all(admissible(x) for x in cov)
        assert dot(row, cov) == rat(2)

    d = (PHI, PHI)
    assert sign(add(d, rat(-2))) > 0
    assert parse(meta['dotted_magnitude']) == d
    c = mul(add(d, O), (ki((16*a[0], 16*a[1])), K0))
    w = scale(add(d, neg(O)), F(1, 8))
    patterns = []
    for ii, jj in product(combinations(range(4), 2), combinations(range(4, 8), 2)):
        plus = set(ii+jj)
        eta = [1 if i in plus else -1 for i in range(8)]
        r = tuple(add(c, scale(w, e)) for e in eta) + tuple(add(c, scale(w, -e)) for e in eta)
        gr = mv(g, r)
        expected = tuple(neg(O if e == 1 else d) for e in eta) + tuple(neg(d if e == 1 else O) for e in eta)
        assert gr == expected and dot(r, gr) == rat(2)
        assert all(admissible(x) for x in gr)
        rows.append(r)
        patterns.append(sorted(plus))
    assert patterns == meta['new_cap_patterns']
    assert rows[192:] == [tuple(map(parse, r)) for r in meta['new_cap_normals']]
    covectors = [mv(g, r) for r in rows]
    assert len(set(covectors)) == 228, 'Two cap mirrors coincide'
    assert all(v not in g for v in covectors), 'A cap repeats a seed wall'
    assert all(dot(r, gr) == rat(2) and sign(total(gr)) < 0 for r, gr in zip(rows, covectors))
    print('Independent seed, all 228 norms, leaf pairings, and strict center: PASS', flush=True)

    graph = nx.Graph()
    graph.add_nodes_from(range(228))
    edges = []
    classes = Counter()
    for i, j in combinations(range(228), 2):
        p = dot(rows[i], covectors[j])
        if admissible(p):
            assert sign(p) <= 0
            graph.add_edge(i, j)
            edges.append([i, j, enc(p)])
            if p == Z:
                classes['right'] += 1
            elif p == rat(-1):
                classes['pi_over_3'] += 1
            elif p == (kn(PHI), K0):
                classes['pi_over_5'] += 1
            elif p == rat(-2):
                classes['parallel'] += 1
            else:
                classes['strictly_dotted'] += 1
        else:
            classes['rejected'] += 1
    assert len(edges) == meta['compatibility_edges'] == 8460
    assert digest(edges) == meta['edge_digest']
    print('Independent all-pair classification: 8460 edges, PASS', flush=True)
    cliques = sorted(tuple(sorted(c)) for c in nx.find_cliques(graph))
    assert len(cliques) == meta['maximal_clique_count'] == 156990
    assert max(map(len, cliques)) == meta['maximum_caps'] == 14
    assert digest(cliques) == meta['clique_digest']
    selected = tuple(meta['selected_fourteen_caps'])
    assert len(selected) == 14 and selected in cliques
    assert all(graph.has_edge(i, j) for i, j in combinations(selected, 2))
    print('Independent NetworkX clique census and archive digest: PASS', flush=True)

    cutmasks = [witness(g, rows, raw) for raw in saved['witnesses']]
    assert len(cutmasks) == 12
    covered = Counter()
    individual = [0]*len(cutmasks)
    for clique in cliques:
        mask = sum(1 << i for i in clique)
        indices = [i for i, cut in enumerate(cutmasks) if not mask & cut]
        assert indices, 'A maximal clique has no certified strict escape'
        covered[indices[0]] += 1
        for i in indices:
            individual[i] += 1
    assert sum(covered.values()) == len(cliques)
    assert saved['remaining_maximal_cliques'] == 0
    assert saved['status'] == 'PASS_EXACT_COMBINED228_BANK'
    selected_mask = sum(1 << i for i in selected)
    selected_witnesses = [i for i, cut in enumerate(cutmasks) if not selected_mask & cut]
    assert selected_witnesses
    print('Twelve exact strict escapes cover every maximal clique: PASS', flush=True)

    regressions = []
    forged = list(rows[192]); forged[0] = add(forged[0], rat(2))
    try:
        check_norm(g, forged)
    except AssertionError:
        regressions.append('changed_new_cap_rejected')
    else:
        raise AssertionError('Forged cap was accepted')
    forged_w = copy.deepcopy(saved['witnesses'][0])
    forged_w['coefficients'] = ['0']*16
    try:
        witness(g, rows, forged_w)
    except AssertionError:
        regressions.append('zero_escape_rejected')
    else:
        raise AssertionError('Forged escape was accepted')
    forged_w = copy.deepcopy(saved['witnesses'][0])
    forged_w['cut_caps'] = []
    try:
        witness(g, rows, forged_w)
    except AssertionError:
        regressions.append('forged_coverage_rejected')
    else:
        raise AssertionError('Forged coverage was accepted')
    return {
        'status': 'PASS_INDEPENDENT_COMBINED228_CAP_BANK',
        'source': str(SOURCE.relative_to(ROOT)),
        'certificate': str(CERT.relative_to(ROOT)),
        'arithmetic': 'Independent Fraction implementation of Q(phi,sqrt(phi)); no repository arithmetic imports',
        'graph_implementation': 'networkx.find_cliques',
        'networkx_version': nx.__version__,
        'seed_signature': [9, 1, 6],
        'seed_center_square': enc(center_square),
        'cap_count': len(rows),
        'pair_count': 228*227//2,
        'pair_class_counts': dict(sorted(classes.items())),
        'compatibility_edges': len(edges),
        'edge_digest': digest(edges),
        'maximal_clique_count': len(cliques),
        'maximal_clique_size_histogram': dict(sorted(Counter(map(len, cliques)).items())),
        'clique_digest': digest(cliques),
        'maximum_caps': max(map(len, cliques)),
        'selected_fourteen_caps': list(selected),
        'selected_fourteen_escape_indices': selected_witnesses,
        'exact_witnesses': len(cutmasks),
        'individual_witness_clique_coverage': individual,
        'first_witness_clique_coverage': [covered[i] for i in range(len(cutmasks))],
        'uncovered_maximal_cliques': 0,
        'negative_controls': regressions,
        'scope': 'All compatible subsets of the literal combined 228-cap dictionary have infinite volume. Arbitrary other caps, the 3040-cap numerical scout, and all completions of the seed are not excluded.'
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = audit()
    if args.verify:
        assert json.loads(OUT.read_text()) == json.loads(json.dumps(result))
    else:
        OUT.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'])


if __name__ == '__main__':
    main()
