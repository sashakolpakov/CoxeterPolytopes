#!/usr/bin/env python3
"""Universal affine B11 center/terminal transfer; exact local identities.

The accompanying proof handles all native contacts and spectator ranks.
No finite-forest maximum or new ambient dimension bound is asserted here.
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
TAG = 'fresh_order12_b_center_terminal_transfer_20261018'
STEMS = ('fresh_order11_order17_dimension63_20261018',
         'fresh_order12_bc_unit_coupon_transfer_20261017',
         'fresh_affine_d13_long_pair_transfer_20261017',
         'fresh_e6_b10_half_coupon_forest_budget_20261018',
         'fresh_order12_payment_support_20261017')


def affine():
    g = L.gram(12, [(i, i+1, s.sqrt(2) if i == 8 else 1) for i in range(9)]
               + [(0, 10, 1), (0, 11, 1)])
    nu = s.Matrix([1]*9+[1/s.sqrt(2), s.Rational(1, 2), s.Rational(1, 2)])
    assert g*nu == s.zeros(12, 1) and g.rank() == 11
    L.pd(g[:11, :11])
    return g, nu


def profiles():
    t, u, v = s.symbols('t u v', nonnegative=True)
    center = [dict(name='outer', contacts={0: t, 10: u, 11: v},
                   c0=2-(u*u+v*v)/2, moment=0, height=t+(u+v)/2),
              dict(name='inner3', contacts={0: t, 1: 1, 10: u, 11: v},
                   c0=1-(u*u+v*v)/2, moment=1, height=t+1+(u+v)/2),
              dict(name='inner4', contacts={0: t, 1: s.sqrt(2)},
                   c0=0, moment=s.sqrt(2), height=t+s.sqrt(2)),
              dict(name='inner_next3', contacts={0: t, 2: 1},
                   c0=0, moment=2, height=t+1)]
    heavy = [dict(name='isolated', contacts={9: t}, c0=2, moment=0, height=t/s.sqrt(2)),
             dict(name='end3', contacts={9: t, 8: 1}, c0=1, moment=1, height=t/s.sqrt(2)+1),
             dict(name='end4', contacts={9: t, 8: s.sqrt(2)}, c0=0, moment=s.sqrt(2), height=t/s.sqrt(2)+s.sqrt(2)),
             dict(name='next3', contacts={9: t, 7: 1}, c0=0, moment=2, height=t/s.sqrt(2)+1)]
    return center, heavy, (t, u, v)


def formulas():
    g, nu = affine()
    W = [i for i in range(12) if i not in (0, 9)]
    J = g.extract(W, W).inv()
    va, vb = [-nu[i]*g.extract([i], W).T for i in (0, 9)]
    old = s.diag(nu[0], nu[9])*g.extract([0, 9], [0, 9])*s.diag(nu[0], nu[9])
    old -= s.Matrix.hstack(va, vb).T*J*s.Matrix.hstack(va, vb)
    assert old.applyfunc(s.simplify) == s.Matrix([[1, -1], [-1, 1]])/9
    inverse_rows = []
    for ids, endpoint, next_node in [(list(range(1, 10)), 1, 2), (list(range(9))+[10, 11], 8, 7)]:
        inverse = g.extract(ids, ids).inv()
        diagonals = dict(zip(ids, inverse.diagonal()))
        assert diagonals[endpoint] == 1 and diagonals[next_node] == 2
        assert all(value > 2 for i, value in diagonals.items() if i not in (endpoint, next_node))
        inverse_rows.append(dict(nodes=ids, diagonal={str(i): str(v) for i, v in diagonals.items()}))
    center, heavy, variables = profiles()
    rows = []
    for kind, ps, far in [('center', center, vb), ('heavy', heavy, va)]:
        for p in ps:
            w = s.Matrix([p['contacts'].get(i, 0) for i in W])
            assert s.simplify(sum(nu[i]*q for i, q in p['contacts'].items())-p['height']) == 0
            assert s.simplify(9*(far.T*J*w)[0]-p['moment']) == 0
            assert s.simplify(2-(w.T*J*w)[0]-p['c0']-p['moment']**2/9) == 0
            rows.append(dict(end=kind, profile=p['name'], c0=str(p['c0']), moment=str(p['moment']),
                             height=str(p['height']), cost=str(s.factor(2*p['moment']/p['height']+p['c0']/p['height']**2))))
    cross = []
    for p, q in product(center, heavy):
        x = s.Matrix([p['contacts'].get(i, 0) for i in W])
        y = s.Matrix([q['contacts'].get(i, 0) for i in W])
        assert s.simplify((x.T*J*y)[0]-p['moment']*q['moment']/9) == 0
        cross.append(dict(center=p['name'], heavy=q['name'], product=str(s.simplify((x.T*J*y)[0]))))
    h, k, l, m, c0, d0, p, q, r = s.symbols('h k l m c0 d0 p q r', nonzero=True)
    A = s.Rational(1, 9)
    T = s.Matrix([[A, -A, -h+l*A, -m*A], [-A, A, -l*A, -k+m*A],
                  [-h+l*A, -l*A, c0+l*l*A-p, -l*m*A-q], [-m*A, -k+m*A, -l*m*A-q, d0+m*m*A-r]])
    expected = 1-(2*l/h+c0/h**2+2*m/k+d0/k**2)/9+(p/h**2-2*q/(h*k)+r/k**2)/9
    assert s.factor(T.det()/(h*h*k*k)-expected) == 0
    t, u, v = variables
    hc = [s.factor(2*z['moment']/z['height']+z['c0']/z['height']**2) for z in heavy]
    assert hc[0].subs(t, 1) == 4 and hc[0].subs(t, s.sqrt(2)) == 2
    assert hc[1].subs(t, 0) == 3 and hc[2].subs(t, 0) == 2 and hc[3].subs(t, 0) == 4
    assert all(s.diff(z, t).subs(t, 1).is_negative for z in hc)
    assert s.simplify(hc[3].subs(t, 1)-3).is_negative
    cc = [s.factor(2*z['moment']/z['height']+z['c0']/z['height']**2) for z in center]
    assert cc[0].subs({t: 0, u: 1, v: 0}) == 6
    assert cc[0].subs({t: 0, u: s.sqrt(2), v: 0}) == 2
    assert cc[1].subs({t: 0, u: 0, v: 0}) == 3
    assert cc[2].subs(t, 0) == 2 and cc[3].subs(t, 0) == 4
    return dict(resistance=9, finite_inverse_diagonals=inverse_rows, end_rows=rows,
                cross_rows=cross, determinant=str(expected),
                center_extreme_cost=6, center_nonextreme_maximum=4, heavy_maximum=4,
                nonextreme_gap='1/9', spectator_weighted_energy_maximum=1,
                surviving_heavy_profiles=[dict(name='isolated', t='1', cost=4, height='sqrt(2)/2', energy=1),
                    dict(name='end3', t='0', cost=3, height='1', energy=0),
                    dict(name='next3', t='0', cost=4, height='1', energy=1)])


def rational_gram(g, scaled):
    scales = [s.sqrt(2) if i in scaled else 1 for i in range(g.rows)]
    d = s.diag(*scales)
    out = (d*g*d).applyfunc(s.simplify)
    assert all(v.is_Rational for v in out)
    return out, [2 if i in scaled else 1 for i in range(g.rows)]


def finite_components(g, ids, gauge):
    assert L.inertia(g.extract(ids, ids)) == [len(ids), 0, 0]
    adj = {i: {j for j in ids if i != j and g[i, j]} for i in ids}
    todo, out = set(ids), []
    while todo:
        nodes, queue = {min(todo)}, [min(todo)]
        for i in queue:
            for j in sorted(adj[i]-nodes):
                nodes.add(j); queue.append(j)
        todo -= nodes
        ns = sorted(nodes)
        degrees = {i: len(adj[i]) for i in ns}
        assert sum(degrees.values()) == 2*(len(ns)-1)
        double = [(i, j) for i in ns for j in adj[i] if i < j and g[i, j]**2 == 2*gauge[i]*gauge[j]]
        assert all(g[i, j]**2 in (gauge[i]*gauge[j], 2*gauge[i]*gauge[j]) for i in ns for j in adj[i])
        branches = [i for i in ns if degrees[i] == 3]
        if double:
            assert len(double) == 1 and not branches and max(degrees.values()) <= 2
            terminal = [i for i in double[0] if degrees[i] == 1]
            assert len(terminal) == 1
            kind, roles = 'B'+str(len(ns)), {str(terminal[0]): 'doubled_terminal'}
        elif branches:
            assert len(branches) == 1
            leaves = [i for i in adj[branches[0]] if degrees[i] == 1]
            assert len(leaves) == 2
            kind = 'D'+str(len(ns))
            roles = {str(branches[0]): 'fork_center', **{str(i): 'short_leaf' for i in leaves}}
        else:
            assert max(degrees.values()) <= 2
            kind, roles = 'A'+str(len(ns)), {str(i): 'endpoint' for i in ns if degrees[i] <= 1}
        out.append(dict(type=kind, nodes=ns, roles=roles))
    return out


def literals():
    f, _ = affine()
    rows = []
    for leaf, name in product((10, 11), ('isolated', 'end3', 'next3')):
        size = 15 if name == 'isolated' else 14
        g = 2*s.eye(size); g[:12, :12] = f
        g[12, leaf] = g[leaf, 12] = -1
        target = {'isolated': 9, 'end3': 8, 'next3': 7}[name]
        g[13, target] = g[target, 13] = -1
        if name == 'isolated':
            g[13, 14] = g[14, 13] = -1
        scaled = {9, 13, 14} if name == 'isolated' else {9}
        rat, gauge = rational_gram(g, scaled)
        W = [i for i in range(size) if i not in (0, 9, 12, 13)]
        # next3 has a whole affine corner; its forced spectator block is not
        # specified by this zero-spectator template. Only that corner is used.
        if name != 'next3':
            assert L.inertia(rat) == [size-2, 1, 1]
            assert L.inertia(rat.extract(W, W)) == [size-4, 0, 0]
            for black in (0, 9, 12, 13):
                assert L.inertia(rat.extract(W+[black], W+[black])) == [size-3, 0, 0]
        corners = []
        total = F()
        for marks in ((0, 9), (9, 12), (12, 13), (13, 0)):
            ids = sorted(set(W)|set(marks))
            sig = L.inertia(rat.extract(ids, ids))
            if name == 'next3' and marks != (13, 0):
                continue
            assert sig[1] == 0
            if sig[2]:
                if marks == (0, 9):
                    corners.append(dict(marks=list(marks), signature=sig, payment='0', kind='removed_B11'))
                    continue
                assert name == 'next3' and marks == (13, 0)
                block = g.extract(ids, ids)
                kernel = block.nullspace()[0]
                if kernel[0] < 0: kernel = -kernel
                assert all(v > 0 for v in kernel)
                adj = L.graph(block); distance = L.distances(adj, ids.index(0))[ids.index(13)]
                assert distance == 8 and len(ids) == 12 and all(block[i, j] in (0, -1) for i in range(12) for j in range(i+1, 12))
                assert sorted(len(adj[i]) for i in adj).count(3) == 2
                assert len(adj[ids.index(0)]) == 3 and len(adj[ids.index(13)]) == 1
                corners.append(dict(marks=list(marks), signature=sig, kind='whole_affine_D11',
                                    positive_kernel=list(map(str, kernel)), marked_distance=8, payment='1'))
                total += 1
                continue
            parts = finite_components(rat, ids, gauge)
            marked = [next(p for p in parts if m in p['nodes']) for m in marks]
            payment, reason = F(), 'none'
            if marked[0] is marked[1]:
                assert marked[0]['type'] == 'D12'
                distance = L.distances({i: {j for j in ids if i != j and rat[i, j]} for i in ids}, marks[0])[marks[1]]
                assert distance == 9
                payment, reason = F(1), 'whole_D12_center_long_endpoint_unit'
            else:
                for i, j in ((0, 1), (1, 0)):
                    small, large = marked[i], marked[j]
                    role = large['roles'].get(str(marks[j]))
                    if small['type'] == 'A2' and large['type'] == 'D11' and role == 'fork_center':
                        payment, reason = F(1, 2), 'existing_A2_D11_fork_center_half'
                    if small['type'] == 'A2' and large['type'] == 'B9' and role == 'doubled_terminal':
                        payment, reason = F(1, 2), 'new_A2_B9_doubled_terminal_half'
            total += payment
            corners.append(dict(marks=list(marks), signature=sig, components=parts, payment=str(payment), reason=reason))
        assert total >= 1
        rows.append(dict(center_leaf=leaf, heavy_profile=name, gram=[[str(x) for x in row] for row in rat.tolist()],
                         gauge_squared=gauge, supporting_white=W, corners=corners, total=str(total),
                         arbitrary_support_scope='Whole affine corner sealed by its positive kernel; other spectator contacts immaterial.' if name == 'next3'
                         else 'All unlisted supporting roots are orthogonal by the zero-energy or isolated-A1 equality proof.'))
    assert len(rows) == 6
    return rows


def run():
    docs = {tag: json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag, d in docs.items():
        assert d['status'].startswith('PASS')
        for p, h in d['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, (tag, p)
    target = next(r for r in docs['fresh_order12_payment_support_20261017']['records'] if r['kind'] == 'B')
    paid = set(map(tuple, target['basic_pairs']))
    assert len(paid) == 64 and (0, 9) in paid
    g, _ = affine(); adj = L.graph(g)
    assert L.distances(adj, 0)[9] == 9
    assert all(L.distances(adj, a)[b] <= 8 for a, b in target['short8_pairs'])
    assert (0, 9) not in set(map(tuple, target['short8_pairs']))
    paths = [Path(__file__), Path(L.__file__)] + [ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS
            for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    return dict(status='PASS_UNIVERSAL_ORDER12_B_CENTER_TERMINAL_TRANSFER',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
                affine_type='B11', order=12, removed_pair=[0, 9], old_cost=64, new_cost=63,
                formulas=formulas(), literal_cases=literals(),
                new_paid_pairs=[list(p) for p in sorted(paid-{(0, 9)})],
                new_finite_coupon=dict(small='whole A2', small_marks='either endpoint',
                    large='whole B9', large_mark='doubled terminal', pair_weight='1/2', component_pair_cost='1'),
                finite_budget_recomputed=False, new_dimension_bound=None,
                scope='Universal actual-quadrilateral proof with all native profiles and arbitrary positive-definite spectators. Exact identities, inverse diagonals and six literal boundary alternatives are checked. Finite-forest budget is separate. No chamber construction or new global bound.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--verify', action='store_true'); args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        assert not path.exists(), 'Existing artifact; use --verify.'
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['new_cost'], 'six boundary cases', flush=True)
