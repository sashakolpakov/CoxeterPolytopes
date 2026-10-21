#!/usr/bin/env python3
"""Coupled universal B/C/D order11 half transfers, with one typed coupon.

The accompanying proof establishes arbitrary supporting-root scope.
Exact identities and literal boundary corners are checked here.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_order11_bcd_half_pair_transfer_20261018'
STEMS = ('fresh_cycle11_half_pair_transfer_20261018',
         'fresh_spectator_five_quarters_transfer_20261017',
         'fresh_spectator_energy_below_one_20261017',
         'fresh_affine_d11_leaf_pair_transfer_20261017',
         'fresh_large_finite_d_leaf_half_transfer_20261017',
         'fresh_order13_affine_b_singleton_transfer_20261017',
         'fresh_deep_finite8_dimension64_20261018')


def gram(n, edges):
    g = 2*s.eye(n)
    for i, j, q in edges:
        g[i, j] = g[j, i] = -q
    return g


def graph(g):
    return {i: {j for j in range(g.rows) if j != i and g[i, j] != 0} for i in range(g.rows)}


def distances(adj, start):
    ds = {start: 0}
    queue = [start]
    for i in queue:
        for j in sorted(adj[i]):
            if j not in ds:
                ds[j] = ds[i]+1
                queue.append(j)
    return ds


def pd(g):
    g = g.copy()
    pivots = []
    while g.rows:
        p = s.simplify(g[0, 0])
        assert p.is_positive is True, p
        pivots.append(str(p))
        g = (g[1:, 1:]-g[1:, :1]*g[:1, 1:]/p).applyfunc(s.simplify)
    return pivots


def inertia(matrix):
    g = [[F(x) for x in row] for row in matrix.tolist()]
    pos = neg = zero = 0
    while g:
        n = len(g)
        i = next((i for i in range(n) if g[i][i]), None)
        if i is not None:
            ids = [i]+[j for j in range(n) if j != i]
            g = [[g[a][b] for b in ids] for a in ids]
            p = g[0][0]
            pos += p > 0
            neg += p < 0
            g = [[g[a][b]-g[a][0]*g[0][b]/p for b in range(1, n)] for a in range(1, n)]
        else:
            pair = next(((i, j) for i in range(n) for j in range(i+1, n) if g[i][j]), None)
            if pair is None:
                zero += n
                break
            i, j = pair
            ids = [i, j]+[k for k in range(n) if k not in pair]
            g = [[g[a][b] for b in ids] for a in ids]
            p = g[0][1]
            pos += 1
            neg += 1
            g = [[g[a][b]-(g[a][0]*g[1][b]+g[a][1]*g[0][b])/p for b in range(2, n)] for a in range(2, n)]
    return [pos, neg, zero]


def affine(kind):
    if kind == 'C':
        edges = [(i, i+1, s.sqrt(2) if i in (0, 9) else 1) for i in range(10)]
        nu = s.Matrix([1/s.sqrt(2)]+[1]*9+[1/s.sqrt(2)])
        marks = (0, 10)
    elif kind == 'B':
        edges = [(i, i+1, s.sqrt(2) if i == 7 else 1) for i in range(8)]+[(0, 9, 1), (0, 10, 1)]
        nu = s.Matrix([1]*8+[1/s.sqrt(2), s.Rational(1, 2), s.Rational(1, 2)])
        marks = (8, 9)
    else:
        assert kind == 'D'
        edges = [(i, i+1, 1) for i in range(6)]+[(0, 7, 1), (0, 8, 1), (6, 9, 1), (6, 10, 1)]
        nu = s.Matrix([1]*7+[s.Rational(1, 2)]*4)
        marks = (7, 9)
    g = gram(11, edges)
    assert g*nu == s.zeros(11, 1) and g.rank() == 10
    pd(g[:10, :10])
    return g, nu, marks


def profiles(kind, mark, nu):
    g, _, _ = affine(kind)
    adj = graph(g)
    neighbor = next(iter(adj[mark]))
    t = s.Symbol('t', nonnegative=True)
    if nu[mark] == s.Rational(1, 2):
        mate = next(i for i in adj[neighbor] if i != mark and len(adj[i]) == 1)
        positions = (mate, neighbor)
        end = 'fork'
        heights = (t/2, (t+1)/2, (t+s.sqrt(2))/2, 1+t/2)
    else:
        next_node = next(i for i in adj[neighbor] if i != mark)
        positions = (neighbor, next_node)
        end = 'heavy'
        heights = (t/s.sqrt(2), t/s.sqrt(2)+1, t/s.sqrt(2)+s.sqrt(2), t/s.sqrt(2)+1)
    out = []
    for name, extra, c0, moment, h in zip(
            ('isolated', 'end3', 'end4', 'next3'),
            ({}, {positions[0]: 1}, {positions[0]: s.sqrt(2)}, {positions[1]: 1}),
            (2, 1, 0, 0), (0, 1, s.sqrt(2), 2), heights):
        out.append(dict(name=name, contacts={mark: t, **extra}, c0=c0, moment=moment,
                        height=h, variable=t, end=end))
    return out


def marked_half_type(g, marks):
    adj = graph(g)
    heavy = [(i, j) for i in range(g.rows) for j in range(i+1, g.rows) if g[i, j]**2 == 2]
    branches = [i for i in adj if len(adj[i]) == 3]
    distance = distances(adj, marks[0])[marks[1]]
    assert all(len(adj[i]) == 1 for i in marks)
    if len(heavy) == 2 and not branches:
        assert distance == 10
        assert all(any(i in e for e in heavy) for i in marks)
        return 'C'
    if len(heavy) == len(branches) == 1:
        assert distance == 9
        assert sum(any(i in e for e in heavy) for i in marks) == 1
        return 'B'
    assert not heavy and len(branches) == 2 and distance == 8
    return 'D'


def formulas():
    rows, cross, immediate = [], [], []
    R = 10
    for kind in 'BCD':
        f, nu, (a, b) = affine(kind)
        W = [i for i in range(11) if i not in (a, b)]
        J = f.extract(W, W).inv().applyfunc(s.simplify)
        va, vb = [-nu[i]*f.extract([i], W).T for i in (a, b)]
        old = s.diag(nu[a], nu[b])*f.extract([a, b], [a, b])*s.diag(nu[a], nu[b])
        old -= s.Matrix.hstack(va, vb).T*J*s.Matrix.hstack(va, vb)
        assert old.applyfunc(s.simplify) == s.Matrix([[1, -1], [-1, 1]])/R
        for mark, other, far in ((a, b, vb), (b, a, va)):
            for p in profiles(kind, mark, nu):
                w = s.Matrix([p['contacts'].get(i, 0) for i in W])
                assert s.simplify(sum(nu[i]*q for i, q in p['contacts'].items())-p['height']) == 0
                assert s.simplify((far.T*J*w)[0]-p['moment']/R) == 0
                assert s.simplify(2-(w.T*J*w)[0]-p['c0']-p['moment']**2/R) == 0
                cost = s.factor(2*p['moment']/p['height']+p['c0']/p['height']**2)
                t = p['variable']
                assert s.diff(cost, t).subs(t, 1).is_negative is True
                if p['name'] == 'isolated':
                    assert s.simplify(cost.subs(t, 1)) == (8 if p['end'] == 'fork' else 4)
                    assert s.simplify(cost.subs(t, s.sqrt(2))) == (4 if p['end'] == 'fork' else 2)
                if p['name'] == 'end3':
                    assert s.simplify(cost.subs(t, 0)) == (8 if p['end'] == 'fork' else 3)
                    if p['end'] == 'fork':
                        assert s.simplify(cost.subs(t, 1)) == 3
                rows.append(dict(kind=kind, mark=mark, end=p['end'], profile=p['name'],
                                 c0=str(p['c0']), moment=str(p['moment']), height=str(p['height']), cost=str(cost)))
                if p['name'] in ('end4', 'next3'):
                    ids = [i for i in range(11) if i != mark]
                    g = 2*s.eye(11)
                    g[:10, :10] = f.extract(ids, ids)
                    for i, q in p['contacts'].items():
                        if i != mark:
                            g[ids.index(i), 10] = g[10, ids.index(i)] = -q
                    null = g.nullspace()
                    assert len(null) == 1
                    kernel = null[0] if null[0][0].is_positive else -null[0]
                    assert all(v.is_positive for v in kernel)
                    target = marked_half_type(g, (10, ids.index(other)))
                    immediate.append(dict(kind=kind, mark=mark, profile=p['name'], target=target,
                                          positive_kernel=list(map(str, kernel)), marked_distance=distances(graph(g), 10)[ids.index(other)]))
        for p, q in product(profiles(kind, a, nu), profiles(kind, b, nu)):
            x = s.Matrix([p['contacts'].get(i, 0) for i in W])
            y = s.Matrix([q['contacts'].get(i, 0) for i in W])
            value = s.simplify((x.T*J*y)[0])
            assert s.simplify(value-p['moment']*q['moment']/R) == 0
            cross.append(dict(kind=kind, left=p['name'], right=q['name'], value=str(value)))
    h, k, l, m, c0, d0, p, q, r = s.symbols('h k l m c0 d0 p q r', nonzero=True)
    A = s.Rational(1, 10)
    T = s.Matrix([[A, -A, -h+l*A, -m*A], [-A, A, -l*A, -k+m*A],
                  [-h+l*A, -l*A, c0+l*l*A-p, -l*m*A-q], [-m*A, -k+m*A, -l*m*A-q, d0+m*m*A-r]])
    identity = 1-(2*l/h+c0/h**2+2*m/k+d0/k**2)/10+(p/h**2-2*q/(h*k)+r/k**2)/10
    assert s.factor(T.det()/(h*h*k*k)-identity) == 0
    assert len(rows) == 24 and len(cross) == 48 and len(immediate) == 12
    return dict(resistance=10, end_profiles=rows, cross_identities=cross,
                whole_affine_half_alternatives=immediate, determinant=str(identity),
                remaining_heavy_maximum='4', remaining_heavy_end3_maximum='3',
                fork_extreme_height='1/2', fork_extreme_cost='8', fork_nonextreme_maximum='4',
                extreme_nonextreme_energy_maximum='2', two_extreme_unscaled_energy='3/2')


def components(g, ids):
    unseen = set(ids)
    out = []
    while unseen:
        start = min(unseen)
        seen, queue = {start}, [start]
        for i in queue:
            for j in sorted(unseen-seen):
                if g[i, j] != 0:
                    seen.add(j)
                    queue.append(j)
        unseen -= seen
        out.append(sorted(seen))
    return out


def selected_small_corners():
    rows = []
    # Canonical automorphisms cover every marked pair and both D ends.
    for kind, max_rank in (('B', 1), ('D', 3)):
        f, _, (a, b) = affine(kind)
        old = b if kind == 'B' else a
        far = a if kind == 'B' else b
        adj = graph(f)
        center = next(iter(adj[old]))
        mate = next(i for i in adj[center] if i != old and len(adj[i]) == 1)
        for rank, attach in product(range(max_rank+1), (old, mate)):
            g = 2*s.eye(12+rank)
            g[:11, :11] = f
            g[11, attach] = g[attach, 11] = -1
            if rank:
                g[11, 12] = g[12, 11] = -1
                for j in range(rank-1):
                    g[12+j, 13+j] = g[13+j, 12+j] = -1
            ids = [i for i in range(g.rows) if i != old]
            pivots = pd(g.extract(ids, ids))
            parts = components(g, ids)
            small = next(p for p in parts if 11 in p)
            large = next(p for p in parts if far in p)
            if attach == old:
                assert small != large and len(small) == rank+1 and len(large) == 10
                old_adj = graph(g.extract(large, large))
                j = large.index(far)
                assert len(old_adj[j]) == 1
                if kind == 'B':
                    neighbor = next(iter(old_adj[j]))
                    assert g[large[j], large[neighbor]]**2 == 2
                    payment = 'new wholeA1/B10 doubled-terminal half' if rank == 0 else 'existing wholeA2/B10 doubled-terminal half'
                else:
                    assert len(old_adj[next(iter(old_adj[j]))]) == 3
                    payment = f'existing wholeA{rank+1}/D10 short-leaf half'
                weight = '1/2'
            else:
                assert small == large and len(large) == 11+rank <= 14
                payment, weight = 'whole finite component<=14 unit', '1'
            rows.append(dict(kind=kind, attached_spectator_A_rank=rank, old_attachment=attach,
                             selected_corner_marks=[11, far], components=parts,
                             positive_pivots=pivots, payment=payment, weight=weight))
    assert len(rows) == 12
    return rows


def mixed_spectrum():
    # Far components at total3/2: arbitrary-rank A inequalities force only
    # endpoint pairs in A5,A6,A7; D6 and E6 complete the inherited list.
    cases = []
    for kind, a, b, c, energy in (
            ('A5', F(5, 6), F(5, 6), F(1, 6), F(4, 3)),
            ('A6', F(6, 7), F(6, 7), F(1, 7), F(10, 7)),
            ('A7', F(7, 8), F(7, 8), F(1, 8), F(3, 2)),
            ('D6', F(1), F(3, 2), F(1, 2), F(3, 2)),
            ('E6', F(4, 3), F(4, 3), F(2, 3), F(4, 3))):
        assert a+b-2*c == energy and a-c >= 0 and b-c >= 0
        # Either contact above1 is at leastsqrt2. The positive-remainder
        # argument gives the minimum on that constrained quadrant here.
        nonunit = [s.Rational(2*a+b)-2*s.sqrt(2)*s.Rational(c),
                   s.Rational(a+2*b)-2*s.sqrt(2)*s.Rational(c)]
        assert all(s.simplify(z).is_positive for z in (
            s.Rational(a)*s.sqrt(2)-s.Rational(c),
            s.Rational(b)-s.Rational(c)*s.sqrt(2),
            s.Rational(a)-s.Rational(c)*s.sqrt(2),
            s.Rational(b)*s.sqrt(2)-s.Rational(c)))
        assert all(s.simplify(z-s.Rational(3, 2)).is_positive for z in nonunit)
        cases.append(dict(kind=kind, inverse_block=list(map(str, (a, b, c))),
                          unit_energy=str(energy), nonunit_minima=list(map(str, nonunit))))
    # A far distance d requires N>=d+2 and E>=2d/(d+2).
    assert F(14, 9) > F(3, 2)
    bounds = [dict(distance=d, N_at_least=d+2, N_at_most=str(F(2*d*d, 2*d-3))) for d in (4, 5, 6)]
    assert [r['N_at_most'] for r in bounds] == ['32/5', '50/7', '8']
    assert all(F(r['unit_energy'])+F(1, 2) > F(3, 2) for r in cases)
    return dict(A_far_bounds=bounds, cases=cases, exact_far_components=['A7 endpoints', 'D6 long/short leaves'],
                minimum_any_active_component='1/2', no_mixed_resolution='At most3 active components. For3 all are wholeA1; for2 some side is a single wholeA_r withr<=3. Select that root, or a spectator-free root, for a retained half.')


def mixed_corners():
    f, _, (a, b) = affine('D')
    rows = []
    for kind, rank, edges, contacts in (
            ('A7', 7, [(i, i+1) for i in range(6)], (0, 6)),
            ('D6', 6, [(i, i+1) for i in range(3)]+[(3, 4), (3, 5)], (0, 4))):
        for orientation in range(1 if kind == 'A7' else 2):
            for left, right in product((7, 8), (9, 10)):
                g = 2*s.eye(13+rank)
                g[:11, :11] = f
                for i, j in [(11, left), (12, right)]+[(13+i, 13+j) for i, j in edges]:
                    g[i, j] = g[j, i] = -1
                cp = contacts if not orientation else contacts[::-1]
                for root, node in zip((11, 12), cp):
                    g[root, 13+node] = g[13+node, root] = -1
                W = [i for i in range(g.rows) if i not in (a, b, 11, 12)]
                edge_sigs = [inertia(g.extract(W+[root], W+[root])) for root in (a, b, 11, 12)]
                corners = []
                for marks in ((a, b), (b, 11), (11, 12), (a, 12)):
                    ids = W+list(marks)
                    sig = inertia(g.extract(ids, ids))
                    parts = components(g, ids)
                    part = next(p for p in parts if marks[0] in p)
                    paid = False
                    distance = None
                    if marks[1] in part and sig[1] == sig[2] == 0:
                        sub = g.extract(part, part)
                        adj = graph(sub)
                        distance = distances(adj, part.index(marks[0]))[part.index(marks[1])]
                        paid = (distance <= 8 or (distance == 9 and any(len(adj[part.index(v)]) == 1 for v in marks)))
                        # Every such component has order<=18. Thus a leaf9
                        # payment is not among the changed largeD>=20 halves,
                        # and no distance8 payment can have feature depth11.
                        assert len(part) <= 18
                    corners.append(dict(marks=list(marks), signature=sig, components=parts,
                                        marked_distance=distance, retained_finite_unit=paid))
                rejected = any(sig[1] or sig[2] for sig in edge_sigs) or any(r['signature'][1] for r in corners)
                if not rejected:
                    assert any(r['retained_finite_unit'] for r in corners)
                    assert inertia(g) == [g.rows-2, 1, 1]
                rows.append(dict(spectator=kind, orientation=orientation, old_attachments=[left, right],
                                 edge_signatures=edge_sigs, corners=corners,
                                 outcome='excluded by actual edge/corner' if rejected else 'retained finite unit',
                                 gram=[[str(x) for x in row] for row in g.tolist()]))
    assert len(rows) == 12
    return rows


def run():
    packets = {stem: json.loads((ROOT/f'results/{stem}.json').read_text()) for stem in STEMS}
    for stem, packet in packets.items():
        assert packet['status'].startswith('PASS'), stem
        for path, digest in packet['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, (stem, path)
    records = []
    for kind in 'BCD':
        g, _, _ = affine(kind)
        adj = graph(g)
        target_distance = {'B': 9, 'C': 10, 'D': 8}[kind]
        changed = [p for p in combinations(range(11), 2)
                   if distances(adj, p[0])[p[1]] == target_distance]
        assert len(changed) == {'B': 2, 'C': 1, 'D': 4}[kind]
        assert target_distance > 7
        weights = [dict(pair=list(p), weight='1/2' if p in changed else '1') for p in combinations(range(11), 2)]
        total = sum(F(r['weight']) for r in weights)
        records.append(dict(kind=kind, changed_pairs=list(map(list, changed)), changed_distance=target_distance,
                            weighted_pairs=weights, total_cost=str(total)))
    assert [r['total_cost'] for r in records] == ['54', '109/2', '53']
    form, small, spectrum, mixed = formulas(), selected_small_corners(), mixed_spectrum(), mixed_corners()
    paths = [Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS
             for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md',
              ROOT/'results/fresh_finite8_d12_b10_coupon_budget_20261017.md',
              ROOT/'results/fresh_finite7_ideal5_triangle_payments_20261017.md']
    return dict(status='PASS_COUPLED_ORDER11_BCD_HALF_PAIR_TRANSFER',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
                factor_rules=records, formulas=form, selected_whole_small_corners=small,
                three_halves_spectator=spectrum, sealed_far_spectator_corners=mixed,
                new_finite_coupon='Weight1/2 between wholeA1 and the doubled terminal of separate wholeB10; exactly1/2 added per such component pair.',
                new_order11_costs=dict(A='99/2', B='54', C='109/2', D='53'),
                new_order11_envelope='109/22', finite_budget_proved_here=False,
                theorem='Halve the2 terminal/fork pairs of affineB10, the endpoint pair of affineC10 and4 opposite-fork pairs of affineD10. With the stated finite coupon, every affected quadrilateral retains total>=1. Short7 and triangle units stay fixed; the separate cycle11 half theorem composes.',
                scope='Universal local theorem with all native labels and arbitrary supporting spectators. Finite coupon budget and new global scalar composition remain separate; no ambient bound or chamber asserted.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['new_order11_costs'], flush=True)
