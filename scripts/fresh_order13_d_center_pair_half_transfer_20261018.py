#!/usr/bin/env python3
"""Universal affine D12 fork-center half transfer; exact local identities.

The accompanying proof supplies all-native and arbitrary-spectator coverage.
No finite forest budget or ambient dimension bound is inferred here.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sympy as s
import fresh_order11_bcd_half_pair_transfer_20261018 as L

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_order13_d_center_pair_half_transfer_20261018'
BASE = 'fresh_order11_order17_dimension63_20261018'
TAIL = 'fresh_deep_finite8_dimension64_20261018'
BCD = 'fresh_order11_bcd_half_pair_transfer_20261018'
STEMS = (BASE, TAIL, BCD, 'fresh_order11_c_endpoint_zero_transfer_20261018',
         'fresh_affine_d12_leaf_pair_transfer_20261017',
         'fresh_spectator_energy_one_boundary_20261017',
         'fresh_spectator_energy_below_one_20261017',
         'fresh_finite7_ideal5_triangle_payments_20261017')


def affine():
    g = L.gram(13, [(i, i+1, 1) for i in range(8)]
               + [(0, 9, 1), (0, 10, 1), (8, 11, 1), (8, 12, 1)])
    nu = s.Matrix([1]*9+[s.Rational(1, 2)]*4)
    assert g*nu == s.zeros(13, 1) and g.rank() == 12
    return g, nu


def reflected(i):
    return 8-i if i <= 8 else {9: 11, 10: 12, 11: 9, 12: 10}[i]


def formulas():
    g, nu = affine()
    W = [i for i in range(13) if i not in (0, 8)]
    J = g.extract(W, W).inv()
    va, vb = [-g.extract([i], W).T for i in (0, 8)]
    old = g.extract([0, 8], [0, 8])-s.Matrix.hstack(va, vb).T*J*s.Matrix.hstack(va, vb)
    assert old == s.Matrix([[1, -1], [-1, 1]])/8
    inv_rows = []
    for side in (0, 1):
        ids = list(range(1, 9))+[11, 12]
        if side: ids = [reflected(i) for i in ids]
        inv = g.extract(ids, ids).inv()
        assert list(inv.diagonal()) == [1, 2, 3, 4, 5, 6, 7, 8, s.Rational(5, 2), s.Rational(5, 2)]
        inv_rows.append(dict(nodes=ids, diagonal=list(map(str, inv.diagonal()))))
    t, u, v = s.symbols('t u v', nonnegative=True)
    profiles = [dict(name='outer', contacts={0: t, 9: u, 10: v}, moment=0,
                     c0=2-(u*u+v*v)/2, height=t+(u+v)/2),
                dict(name='inner3', contacts={0: t, 1: 1, 9: u, 10: v}, moment=1,
                     c0=1-(u*u+v*v)/2, height=t+1+(u+v)/2),
                dict(name='inner4', contacts={0: t, 1: s.sqrt(2)}, moment=s.sqrt(2),
                     c0=0, height=t+s.sqrt(2)),
                dict(name='next3', contacts={0: t, 2: 1}, moment=2,
                     c0=0, height=t+1)]
    rows = []
    for p in profiles:
        x = s.Matrix([p['contacts'].get(i, 0) for i in W])
        assert s.simplify(sum(nu[i]*q for i, q in p['contacts'].items())-p['height']) == 0
        assert s.simplify(8*(vb.T*J*x)[0]-p['moment']) == 0
        assert s.simplify(2-(x.T*J*x)[0]-p['c0']-p['moment']**2/8) == 0
        rows.append(dict(name=p['name'], contacts={str(i): str(q) for i, q in p['contacts'].items()},
                         height=str(p['height']), c0=str(p['c0']), moment=str(p['moment']),
                         cost=str(s.factor(2*p['moment']/p['height']+p['c0']/p['height']**2))))
    crosses = []
    for p, q in product(profiles, repeat=2):
        x = s.Matrix([p['contacts'].get(i, 0) for i in W])
        y = s.Matrix([q['contacts'].get(reflected(i), 0) for i in W])
        value = s.simplify((x.T*J*y)[0])
        assert s.simplify(value-p['moment']*q['moment']/8) == 0
        crosses.append(dict(left=p['name'], right=q['name'], product=str(value)))
    outer = s.factor(profiles[0]['c0']/profiles[0]['height']**2)
    assert outer.subs({t: 0, u: 1, v: 0}) == 6
    assert outer.subs({t: 0, u: s.sqrt(2), v: 0}) == 2
    assert outer.subs({t: 1, u: 1, v: 0}) == s.Rational(2, 3)
    assert outer.subs({t: 0, u: 1, v: 1}) == 1
    inner = 2/profiles[1]['height']+profiles[1]['c0']/profiles[1]['height']**2
    assert inner.subs({t: 0, u: 0, v: 0}) == 3
    h, k, ell, m, c0, d0, p, q, r = s.symbols('h k ell m c0 d0 p q r', nonzero=True)
    a = s.Rational(1, 8)
    T = s.Matrix([[a, -a, -h+ell*a, -m*a], [-a, a, -ell*a, -k+m*a],
                  [-h+ell*a, -ell*a, c0+ell*ell*a-p, -ell*m*a-q],
                  [-m*a, -k+m*a, -ell*m*a-q, d0+m*m*a-r]])
    identity = 1-(2*ell/h+c0/h**2+2*m/k+d0/k**2)/8+(p/h**2-2*q/(h*k)+r/k**2)/8
    assert s.factor(T.det()/(h*h*k*k)-identity) == 0
    assert 1-F(6, 8) == F(1, 4)
    assert F(1, 2)**2*(12-8) == 1
    return dict(resistance=8, inverse_diagonals=inv_rows, profiles=rows, cross_identities=crosses,
                determinant_identity=str(identity), extreme_cost=6, remaining_nonextreme_maximum=3,
                two_nonextreme_strict_gap='1/4', one_extreme_weighted_energy_maximum=1,
                nonzero_extreme_support_energy_lower_bound=2, both_extreme_unscaled_energy=1)


def components(g):
    adj = L.graph(g); todo = set(adj); out = []
    while todo:
        nodes = {min(todo)}; queue = [min(todo)]
        for i in queue:
            for j in sorted(adj[i]-nodes): nodes.add(j); queue.append(j)
        todo -= nodes; out.append(sorted(nodes))
    return out


def immediate_corners(data):
    f, _ = affine(); out = []
    for side, (name, node, magnitude, kind, expected_distance) in product((0, 1),
            [('inner4', 1, s.sqrt(2), 'B', 8), ('next3', 2, 1, 'D', 7)]):
        g = 2*s.eye(14); g[:13, :13] = f
        mark, far = (0, 8) if side == 0 else (8, 0)
        attach = node if side == 0 else reflected(node)
        g[13, attach] = g[attach, 13] = -magnitude
        ids = [i for i in range(14) if i != mark]
        corner = g.extract(ids, ids)
        block = next(b for b in components(corner) if ids.index(13) in b)
        ns = [ids[i] for i in block]; H = g.extract(ns, ns)
        assert len(ns) == 11 and far in ns and H.rank() == 10
        z = H.nullspace()[0]
        if z[0] < 0: z = -z
        assert all(q > 0 for q in z)
        # A positive null vector and acute entries identify the whole affine
        # factor; verify its canonical marked roles and the retained unit.
        adj = L.graph(H); distance = L.distances(adj, ns.index(13))[ns.index(far)]
        assert distance == expected_distance and len(adj[ns.index(far)]) == 3
        assert len(adj[ns.index(13)]) == 1
        canonical, _, _ = L.affine(kind)
        ca = L.graph(canonical)
        center = next(i for i in ca if len(ca[i]) == 3 and L.distances(ca, i)[8 if kind == 'B' else 7] == expected_distance)
        terminal = 8 if kind == 'B' else 7
        pair = tuple(sorted((center, terminal)))
        rule = next(r for r in data[BCD]['factor_rules'] if r['kind'] == kind)
        weights = {tuple(r['pair']): F(r['weight']) for r in rule['weighted_pairs']}
        assert weights[pair] == 1
        out.append(dict(side=side, profile=name, whole_factor='affine_'+kind+'10', nodes=ns,
                        gram=[[str(x) for x in row] for row in H.tolist()], positive_kernel=list(map(str, z)),
                        marks=[far, 13], distance=distance, canonical_pair=list(pair), retained_payment='1'))
    return out


def finite_corner_literals():
    f, _ = affine(); out = []
    for side, leaf_choice, extra in product((0, 1), (0, 1), (0, 1)):
        size = 14+extra; g = 2*s.eye(size); g[:13, :13] = f
        mark, far = (0, 8) if side == 0 else (8, 0)
        leaf = (9+leaf_choice) if side == 0 else (11+leaf_choice)
        g[13, leaf] = g[leaf, 13] = -1
        if extra: g[13, 14] = g[14, 13] = -1
        ids = [i for i in range(size) if i != mark]; H = g.extract(ids, ids)
        assert L.inertia(H) == [len(ids), 0, 0]
        blocks = [[ids[i] for i in b] for b in components(H)]
        small = next(b for b in blocks if 13 in b); large = next(b for b in blocks if far in b)
        assert len(small) == 2+extra and len(large) == 10
        adj = {i: {j for j in ids if j != i and g[i, j]} for i in ids}
        assert len(adj[13]) == 1+extra
        assert len(adj[far]) == 3 and sum(len(adj[i]) == 3 for i in large) == 1
        assert sum(len(adj[i]) == 1 for i in adj[far]) == 2
        out.append(dict(side=side, attached_old_leaf=leaf, spectator_A1=bool(extra), marks=[far, 13],
                        signature=L.inertia(H), components=blocks, whole_small_type='A'+str(2+extra),
                        small_mark_role='middle' if extra else 'endpoint', whole_large_type='D10',
                        large_mark_role='fork_center', coupon_payment='1/2', original_ideal_payment='1/2', total='1'))
    return out


def energy_one_checks():
    rows = []
    # The finite list corroborates the all-rank theorem, not a rank cutoff.
    for rank in range(3, 8):
        if rank == 3: g = L.gram(3, [(0, 1, 1), (0, 2, 1)]); leaves = (1, 2)
        else:
            g = L.gram(rank, [(i, i+1, 1) for i in range(rank-3)]
                       + [(rank-3, rank-2, 1), (rank-3, rank-1, 1)])
            leaves = (rank-2, rank-1)
        x = s.eye(rank)[:, leaves[0]]; y = s.eye(rank)[:, leaves[1]]; J = g.inv()
        assert ((x-y).T*J*(x-y))[0] == 1
        assert (x.T*J*x)[0] == s.Rational(rank, 4) < 2
        assert L.distances(L.graph(g), leaves[0])[leaves[1]]+2 == 4
        rows.append(dict(type='A3' if rank == 3 else 'D'+str(rank), projection=str(s.Rational(rank, 4)),
                         energy='1', opposite_black_path=4, payment='1'))
    assert (s.Matrix([1, -1]).T*(2*s.eye(2)).inv()*s.Matrix([1, -1]))[0] == 1
    return rows


def run():
    data = {tag: json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag, doc in data.items():
        assert doc['status'].startswith('PASS'), tag
        for path, digest in doc['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, (tag, path)
    old = next(r for r in data[TAIL]['short7_small_factor_interfaces'] if (r['order'], r['kind']) == (13, 'D'))
    weights = {tuple(r['pair']): F(r['weight']) for r in old['weighted_pairs']}
    assert sum(weights.values(), F()) == 74 and weights[0, 8] == 1
    f, _ = affine(); ds = L.distances(L.graph(f), 0)
    assert ds[8] == 8
    weights[0, 8] = F(1, 2)
    assert sum(weights.values(), F()) == F(147, 2)
    for pair, weight in weights.items():
        if L.distances(L.graph(f), pair[0])[pair[1]] <= 7: assert weight == 1
    assert data[BASE]['dimension_bound'] == 63
    assert data['fresh_order11_c_endpoint_zero_transfer_20261018']['new_cost'] == '54'
    paths = [Path(__file__), Path(L.__file__)] + [ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS
             for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    return dict(status='PASS_UNIVERSAL_ORDER13_D_CENTER_PAIR_HALF_TRANSFER',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
                affine_type='D12', order=13, changed_pairs=[[0, 8]], old_weight='1', new_weight='1/2',
                results=dict(old_affine_total='74', new_affine_total='147/2', normalized_cost='147/26'),
                weighted_pairs=[dict(pair=list(pair), weight=str(value)) for pair, value in sorted(weights.items())],
                formulas=formulas(), immediate_affine_unit_corners=immediate_corners(data),
                finite_coupon_corners=finite_corner_literals(), mixed_energy_one_cases=energy_one_checks(),
                new_coupons=[dict(small='whole_A2', small_role='either_endpoint', large='whole_D10',
                                  large_role='fork_center', weight='1/2', component_pair_increment='1'),
                             dict(small='whole_A3', small_role='unique_middle', large='whole_D10',
                                  large_role='fork_center', weight='1/2', component_pair_increment='1/2')],
                actual_corner_priority='If a changed D13 center pair occurs, retain its half and obtain a distinct half coupon or a distinct unit by this theorem. Otherwise every actual weight equals the inherited rule, apart from coupon increases.',
                triangle_minimum=2, quadrilateral_minimum=1, finite_budget_proved_here=False, new_dimension_bound=None,
                scope='Universal all-native supporting-component theorem. Original ideal half is retained; one typed whole-component half suffices. Finite budget and global dimension composition are separate.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--verify', action='store_true'); args = parser.parse_args()
    path = ROOT/f'results/{TAG}.json'; output = json.loads(json.dumps(run()))
    if args.verify: assert output == json.loads(path.read_text())
    else: path.write_text(json.dumps(output, indent=2)+'\n')
    print(output['status'], output['results'], flush=True)
