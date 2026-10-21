#!/usr/bin/env python3
"""Universal four cross-fork-leaf transfers for affine D12 of order13."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sympy as s
import fresh_affine_d14_long_pair_transfer_20261017 as M

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_affine_d12_leaf_pair_transfer_20261017'
C = M.C
MARKS = list(product((9, 10), (11, 12)))


def affine():
    g = 2*s.eye(13)
    for i, j in [(i, i+1) for i in range(8)]+[(0, 9), (0, 10), (8, 11), (8, 12)]:
        g[i, j] = g[j, i] = -1
    nu = s.Matrix([1]*9+[s.Rational(1, 2)]*4)
    assert g*nu == s.zeros(13, 1) and C.inertia(g.tolist()) == [12, 0, 1]
    return g, nu


def profiles(mark):
    center = 0 if mark < 11 else 8
    mate = 19-mark if mark < 11 else 23-mark
    x = s.symbols('x' if mark < 11 else 'y', nonnegative=True)
    return [dict(name='isolated', contacts={mark: x}, c0=2, l=0, h=x/2),
        dict(name='mate3', contacts={mark: x, mate: 1}, c0=1, l=1, h=(x+1)/2),
        dict(name='mate4', contacts={mark: x, mate: s.sqrt(2)}, c0=0,
             l=s.sqrt(2), h=(x+s.sqrt(2))/2),
        dict(name='center3', contacts={mark: x, center: 1}, c0=0, l=2, h=1+x/2)]


def formulas():
    g, nu = affine()
    rows, cross = [], []
    R = 12
    for a, b in MARKS:
        W = sorted(set(range(13))-{a, b})
        J = g.extract(W, W).inv()
        va = -nu[a]*g.extract([a], W).T
        vb = -nu[b]*g.extract([b], W).T
        old = s.diag(nu[a], nu[b])*g.extract([a, b], [a, b])*s.diag(nu[a], nu[b])
        old -= s.Matrix.hstack(va, vb).T*J*s.Matrix.hstack(va, vb)
        assert old == s.Matrix([[1, -1], [-1, 1]])/R
        for side, mark, far in [('left', a, vb), ('right', b, va)]:
            for p in profiles(mark):
                w = s.Matrix([p['contacts'].get(i, 0) for i in W])
                h = sum(nu[i]*q for i, q in p['contacts'].items())
                contact = (far.T*J*w)[0]
                square = 2-(w.T*J*w)[0]
                assert s.simplify(h-p['h']) == 0
                assert s.simplify(contact-p['l']/R) == 0
                assert s.simplify(square-p['c0']-p['l']**2/R) == 0
                rows.append(dict(marks=[a, b], side=side, profile=p['name'],
                    contacts={str(i): str(q) for i, q in p['contacts'].items()},
                    c0=str(p['c0']), moment=str(p['l']), height=str(p['h']),
                    far_contact=str(contact), residual_square=str(s.factor(square))))
        for p, q in product(profiles(a), profiles(b)):
            w = s.Matrix([p['contacts'].get(i, 0) for i in W])
            z = s.Matrix([q['contacts'].get(i, 0) for i in W])
            value = (w.T*J*z)[0]
            assert s.simplify(value-p['l']*q['l']/R) == 0
            cross.append(dict(marks=[a, b], left=p['name'], right=q['name'], value=str(value)))
    assert len(rows) == 32 and len(cross) == 64
    h, k, l, m, c0, d0, p, q, r = s.symbols('h k l m c0 d0 p q r', nonzero=True)
    A = s.Rational(1, R)
    T = s.Matrix([[A, -A, -h+l*A, -m*A], [-A, A, -l*A, -k+m*A],
        [-h+l*A, -l*A, c0+l*l*A-p, -l*m*A-q],
        [-m*A, -k+m*A, -l*m*A-q, d0+m*m*A-r]])
    expected = 1-(2*l/h+c0/h**2+2*m/k+d0/k**2)/R+(p/h**2-2*q/(h*k)+r/k**2)/R
    assert s.factor(T.det()/(h*h*k*k)-expected) == 0
    costs = [x for x in M.cost_checks()['rows'] if x['end_type'] == 'leaf']
    assert len(costs) == 6
    assert 1-F(8+4, R) == 0 and 1-F(8+8, R) == -F(1, 3)
    assert F(16-R, 4) == 1
    return dict(end_rows=rows, cross_rows=cross, native_cost_rows=costs,
        resistance=R, determinant=str(expected), nonextreme_gap='0',
        extreme_height='1/2', extreme_required_unscaled_spectator_energy='1')


def simple_finite_components(g, ids):
    assert all(g[i, i] == 2 for i in ids)
    assert all(g[i, j] in (0, -1) for i in ids for j in ids if i != j)
    assert C.inertia(g.extract(ids, ids).tolist()) == [len(ids), 0, 0]
    parts = M.components(g, ids)
    for p in parts:
        assert p['type'].startswith(('A', 'D')) and not p['type'].startswith('affine')
        adj = {i: [j for j in p['nodes'] if i != j and g[i, j]] for i in p['nodes']}
        assert sum(map(len, adj.values())) == 2*(len(adj)-1)
        if p['type'].startswith('D'):
            assert len(p['branch']) == 1 and len(p['short_leaves']) == 2
    return parts


def payment(g, W, p, q):
    ids = sorted(set(W)|{p, q})
    parts = simple_finite_components(g, ids)
    for part in parts:
        part['marks'] = [v for v in (p, q) if v in part['nodes']]
    A = next(part for part in parts if p in part['nodes'])
    B = next(part for part in parts if q in part['nodes'])
    unit, half, upgraded = 0, F(0), False
    if A is B:
        assert len(A['nodes']) <= 17
        unit = 1
    else:
        for small, large in ((A, B), (B, A)):
            if large['type'] == 'D12' and any(m in large['short_leaves'] for m in large['marks']):
                if small['type'] == 'A1':
                    unit, upgraded = 1, True
                elif len(small['nodes']) <= 5:
                    half = F(1, 2)
    return dict(marks=[p, q], signature=[len(ids), 0, 0], components=parts,
        unit=unit, old_half=str(half), full_A1_D12_upgrade=upgraded)


def put(g, root, mark, name):
    center = 0 if mark < 11 else 8
    mate = 19-mark if mark < 11 else 23-mark
    j = mark if name in ('isolated', 'isolated4') else mate if name in ('mate3', 'mate4') else center
    assert name in ('isolated', 'mate3', 'isolated4', 'mate4', 'center3')
    if name in ('isolated4', 'mate4'):
        # Congruence: scale only the new norm-two root by sqrt2.
        g[root, root] = 4
        value = 2
    else:
        value = 1
    g[root, j] = g[j, root] = -value


def zero_cases():
    f, _ = affine()
    out = []
    for (a, b), side, ext, ne in product(MARKS, ('left', 'right'),
            ('isolated', 'mate3'), ('isolated4', 'mate4', 'center3')):
        g = 2*s.eye(15)
        g[:13, :13] = f
        put(g, 13, a, ne if side == 'left' else ext)
        put(g, 14, b, ext if side == 'left' else ne)
        W = sorted(set(range(13))-{a, b})
        assert C.inertia(g.tolist()) == [13, 1, 1]
        assert C.inertia(g.extract(W, W).tolist()) == [11, 0, 0]
        for root in (a, b, 13, 14):
            edge = W+[root]
            assert C.inertia(g.extract(edge, edge).tolist()) == [12, 0, 0]
        corners = []
        for p, q in ((a, b), (b, 13), (13, 14), (14, a)):
            ids = sorted(set(W)|{p, q})
            sig = C.inertia(g.extract(ids, ids).tolist())
            assert sig[1] == 0
            corners.append(dict(marks=[p, q], signature=sig))
        chosen = payment(g, W, 14, a) if side == 'left' else payment(g, W, b, 13)
        assert chosen['unit'] == 1
        assert chosen['full_A1_D12_upgrade'] == (ext == 'isolated')
        out.append(dict(marks=[a, b], nonextreme_side=side, nonextreme_profile=ne,
            extreme_profile=ext, gram=[[str(v) for v in row] for row in g.tolist()],
            signature=[13, 1, 1], white=W, corner_signatures=corners, chosen=chosen))
    assert len(out) == 48
    return out


def one_sided_shortcut():
    f, _ = affine()
    out = []
    for (a, b), L, R, contacting in product(MARKS, (9, 10), (11, 12), (13, 14)):
        g = 2*s.eye(15)
        g[:13, :13] = f
        for i, j in ((13, L), (14, R)):
            g[i, j] = g[j, i] = -1
        W = sorted(set(range(13))-{a, b})
        p, q = (14, a) if contacting == 13 else (b, 13)
        chosen = payment(g, W, p, q)
        assert chosen['unit'] == 1
        assert chosen['full_A1_D12_upgrade'] == ((R == b) if contacting == 13 else (L == a))
        out.append(dict(marks=[a, b], leaf_attachments=[L, R],
            only_spectator_contacting_black=contacting, chosen=chosen,
            justification='All spectators are orthogonal to F and the other new black. The chosen actual corner excludes the contacting black, so its displayed finite components are whole regardless of the spectator diagram.'))
    assert len(out) == 32
    return out


def two_isolated_spectators():
    f, _ = affine()
    out, tally = [], Counter()
    for (a, b), L, R in product(MARKS, (9, 10), (11, 12)):
        g = 2*s.eye(17)
        g[:13, :13] = f
        for i, j in ((13, L), (14, R), (13, 15), (14, 16)):
            g[i, j] = g[j, i] = -1
        W = sorted(set(range(13))-{a, b})+[15, 16]
        assert C.inertia(g.tolist()) == [15, 1, 1]
        assert C.inertia(g.extract(W, W).tolist()) == [13, 0, 0]
        for root in (a, b, 13, 14):
            edge = W+[root]
            assert C.inertia(g.extract(edge, edge).tolist()) == [14, 0, 0]
        orig = sorted(set(W)|{a, b})
        assert C.inertia(g.extract(orig, orig).tolist()) == [14, 0, 1]
        corners = [payment(g, W, p, q) for p, q in ((b, 13), (13, 14), (14, a))]
        unit = sum(r['unit'] for r in corners)
        half = sum(F(r['old_half']) for r in corners)
        assert unit+half >= 1
        assert not any(r['full_A1_D12_upgrade'] for r in corners)
        tally[(unit, str(half))] += 1
        out.append(dict(marks=[a, b], leaf_attachments=[L, R],
            gram=[[str(v) for v in row] for row in g.tolist()], signature=[15, 1, 1],
            white=W, corners=corners, unit_total=unit, half_total=str(half)))
    assert len(out) == 16
    return dict(rows=out, payment_classes=[dict(units=k[0], halves=k[1], cases=v)
        for k, v in sorted(tally.items())],
        other_two_sided_cases='Whole A3 endpoints or D4..D7 short twins give a path of length4 between the new black marks, hence a retained finite8/ideal5 payment.')


def run():
    stems = ('fresh_affine_d14_long_pair_transfer_20261017',
        'fresh_spectator_energy_one_boundary_20261017', 'fresh_finite8_d12_b10_coupon_budget_20261017',
        'fresh_order13_payment_support_20261017')
    inputs = {}
    for stem in stems:
        d = json.loads((ROOT / f'results/{stem}.json').read_text())
        assert d['status'].startswith('PASS')
        for path, h in d['sources'].items():
            assert sha256((ROOT / path).read_bytes()).hexdigest() == h, (stem, path)
        inputs[stem] = d
    diag = inputs['fresh_order13_payment_support_20261017']
    rec = next(r for r in diag['records'] if r['kind'] == 'D')
    g, _ = affine()
    assert rec['target'] == [[3 if i != j and g[i, j] else 0 for j in range(13)] for i in range(13)]
    paid = set(map(tuple, rec['required_union']))
    basic = set(map(tuple, rec['basic_pairs']))
    assert len(paid) == 78 and len(basic) == 74 and paid-basic == set(MARKS)
    out = dict(formulas=formulas(), nonextreme_zero_cases=zero_cases(),
        arbitrary_one_sided_cases=one_sided_shortcut(), two_A1_spectators=two_isolated_spectators())
    paths = [Path(__file__), Path(C.__file__)] + [ROOT / f'{folder}/{stem}.{ext}'
        for stem in stems for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT / 'results/fresh_cusp_factor_payment_framework_20261017.md']
    out.update(status='PASS_UNIVERSAL_AFFINE_D12_CROSS_FORK_LEAF_TRANSFER',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        removed_pairs=[list(p) for p in MARKS], old_cost=78, new_cost=74,
        new_paid_pairs=[list(p) for p in sorted(basic)],
        counts=dict(marked_pairs=4, end_formulas=32, cross_formulas=64,
            nonextreme_zero_configurations=48, arbitrary_one_sided_corners=32,
            two_A1_configurations=16, two_A1_finite_corners=48),
        theorem='Each actual quadrilateral with marks at opposite short fork leaves of an affine D12 factor transfers to a retained short path, a whole finite component of order<=17, a full wholeA1/D12-shortleaf coupon, or two inherited wholeA2/D12-shortleaf halves.',
        scope='Universal face-local proof for arbitrary supporting spectators and all native labels. Both selected-core and immediate-factor uses of the four pairs are covered. No global dimension improvement or new chamber construction is asserted.')
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    d = json.loads(json.dumps(run()))
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert d == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], d['counts'], flush=True)
