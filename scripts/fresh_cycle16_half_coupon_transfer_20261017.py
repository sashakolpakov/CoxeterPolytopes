#!/usr/bin/env python3
"""Universal affine-cycle16 antipodal transfer; no ambient bound."""
import argparse
import json
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_cycle16_half_coupon_transfer_20261017'


def inertia(a):
    a = [list(map(F, r)) for r in a]
    out = [0, 0, 0]
    while a:
        n = len(a)
        i = next((i for i in range(n) if a[i][i]), None)
        if i is not None:
            q = a[i][i]
            out[0 if q > 0 else 1] += 1
            rest = [j for j in range(n) if j != i]
            a = [[a[j][k] - a[j][i]*a[i][k]/q for k in rest] for j in rest]
        else:
            e = next(((i,j) for i in range(n) for j in range(i+1,n) if a[i][j]), None)
            if e is None:
                out[2] += n
                break
            i,j = e
            q = a[i][j]
            out[0] += 1
            out[1] += 1
            rest = [k for k in range(n) if k not in e]
            a = [[a[k][l] - (a[k][i]*a[j][l] + a[k][j]*a[i][l])/q
                  for l in rest] for k in rest]
    return out


def distances(adj, nodes, a):
    nodes = set(nodes)
    out = {a: 0}
    q = deque([a])
    while q:
        i = q.popleft()
        for j in sorted((adj[i] & nodes) - out.keys()):
            out[j] = out[i] + 1
            q.append(j)
    return out


def tree_pd(nodes, edges):
    # Leaf Schur elimination needs only squared bond magnitudes.
    adj = {i: {} for i in nodes}
    for i,j,w2 in edges:
        adj[i][j] = adj[j][i] = F(w2)
    diag = {i: F(2) for i in nodes}
    pivots = []
    while adj:
        i = next(i for i in sorted(adj) if len(adj[i]) <= 1)
        q = diag[i]
        assert q > 0
        pivots.append(str(q))
        if adj[i]:
            j,w2 = next(iter(adj[i].items()))
            diag[j] -= w2/q
            del adj[j][i]
        del adj[i]
    return pivots


def single_profiles():
    old = list(range(1,16))
    path = [(i,i+1,1) for i in range(1,15)]
    rows = []
    for side in (0,1):
        for name,pos,w2,kind in [('end3',1,1,'A16'), ('end4',1,2,'B16'),
                                 ('next3',2,1,'D16')]:
            contact = pos if side == 0 else 16-pos
            edges = path + [(16,contact,w2)]
            adj = {i:set() for i in old+[16]}
            for i,j,_ in edges:
                adj[i].add(j)
                adj[j].add(i)
            fixed_leaf = 15 if side == 0 else 1
            neighbor = next(iter(adj[fixed_leaf]))
            assert len(adj[fixed_leaf]) == 1 and len(adj[neighbor]) == 2
            dist = distances(adj,old+[16],16)[8]
            assert dist == 9-pos and dist <= 8
            rows.append(dict(profile=f'{name}_{side}',contact=contact,bond_square=w2,
                finite_core_type=kind,order=16,positive_pivots=tree_pd(old+[16],edges),
                marked_distance=dist,fixed_simple_leaf=fixed_leaf,
                fixed_degree_two_neighbor=neighbor))
    assert len(rows) == 6
    return rows


def profile_rectangles():
    profiles = [('isolated',F(2),F(0),F(1)), ('both_ends3',F(0),F(1),F(1,2))]
    rows = []
    for c,d in product(profiles,repeat=2):
        vals = []
        for ub,vb in product((0,1),repeat=2):
            u,v = ub*c[3],vb*d[3]
            # Equal moments on the two equal arcs cancel the square term.
            val = 1-(c[2]*u+d[2]*v)/2-(c[1]*u*u+d[1]*v*v)/4
            vals.append(dict(upper_reciprocals=[ub,vb],value=str(val)))
        minimum = min(F(x['value']) for x in vals)
        expected = F(0) if c[0] == d[0] == 'isolated' else F(1,2) if c[0] == d[0] else F(1,4)
        assert minimum == expected
        rows.append(dict(c=c[0],d=d[0],rectangle_values=vals,minimum=str(minimum)))
    # Isolated equality is exact, not a numerical boundary observation.
    u,v,h,k,p,q,r = s.symbols('u v h k p q r', real=True)
    isolated = 1-(u*u+v*v)/2
    assert s.expand(isolated-((1-u*u)+(1-v*v))/2) == 0
    S = s.Matrix([[s.Rational(1,4),-s.Rational(1,4),-h,0],
                  [-s.Rational(1,4),s.Rational(1,4),0,-k],
                  [-h,0,2,0],[0,-k,0,2]])
    corrected = S.copy()
    corrected[2,2] -= p
    corrected[2,3] -= q
    corrected[3,2] -= q
    corrected[3,3] -= r
    assert s.expand(S.det()-(h*h*k*k-(h*h+k*k)/2)) == 0
    correction = s.factor(corrected.det()-S.det())
    assert s.expand(correction-(p*k*k-2*q*h*k+r*h*h)/4) == 0
    assert s.expand(correction.subs({h:1,k:1})-(p-2*q+r)/4) == 0
    return dict(rows=rows,rectangle_values=16,
        isolated_formula='1-(u^2+v^2)/2 = ((1-u^2)+(1-v^2))/2',
        equality={'h':'1','k':'1','native_a_c':'-1','native_b_d':'-1'},
        spectator_correction='(p/h^2-2*q/(h*k)+r/k^2)/4',
        equality_spectator_correction='(p-2*q+r)/4')


def equality_core():
    nodes = set(range(18))
    edges = [(i,(i+1)%16) for i in range(16)] + [(0,16),(8,17)]
    adj = {i:set() for i in nodes}
    for i,j in edges:
        adj[i].add(j)
        adj[j].add(i)
    g = [[F(2 if i==j else -int(j in adj[i])) for j in range(18)] for i in range(18)]
    assert inertia(g) == [16,1,1]
    radical = [-4,-3,-2,-1,0,1,2,3,4,3,2,1,0,-1,-2,-3,-2,2]
    assert all(sum(g[i][j]*radical[j] for j in nodes) == 0 for i in nodes)
    W = [i for i in range(16) if i not in (0,8)]
    def sub(ids): return [[g[i][j] for j in ids] for i in ids]
    assert inertia(sub(W)) == [14,0,0]
    corners = []
    for a,b in ((0,8),(8,16),(16,17),(17,0)):
        ids = sorted(set(W)|{a,b})
        sig = inertia(sub(ids))
        da = distances(adj,ids,a)
        if (a,b) == (0,8):
            assert sig == [15,0,1] and da[b] == 8
            corners.append(dict(marks=[a,b],kind='ideal_cycle16',distance=8,coupon='0'))
            continue
        assert sig == [16,0,0] and b not in da
        db = distances(adj,ids,b)
        sizes = sorted((len(da),len(db)))
        coupon = F(0)
        parts = []
        for ds,mark in ((da,a),(db,b)):
            leaves = sorted(i for i in ds if len(adj[i]&ds.keys()) <= 1)
            assert max(len(adj[i]&ds.keys()) for i in ds) <= 2
            central = len(ds) == 15 and len(leaves) == 2 and all(distances(adj,ds,mark)[j] == 7 for j in leaves)
            parts.append(dict(type=f'A{len(ds)}',nodes=sorted(ds),mark=mark,
                              unique_middle=central))
        if sizes == [1,15]:
            assert any(p['unique_middle'] for p in parts)
            coupon = F(1,2)
        corners.append(dict(marks=[a,b],kind='finite',signature=sig,
            marked_components=parts,finite9_or_small17_payment='0',coupon=str(coupon)))
    assert sum(F(c['coupon']) for c in corners) == 1
    paid = [[i,j] for i in range(16) for j in range(i+1,16) if min(j-i,16-j+i) <= 7]
    dropped = [[i,i+8] for i in range(8)]
    assert len(paid) == 112 and len(dropped) == 8
    return dict(gram=[[str(x) for x in row] for row in g],signature=[16,1,1],
        radical=radical,white=W,corners=corners,coupon_total='1',
        ideal_short7_pairs=paid,dropped_antipodal_pairs=dropped,
        ideal_factor_cost=112,old_all_pair_cost=120)


def run():
    out = dict(single_profile_checks=single_profiles(),residual_profiles=profile_rectangles(),
               equality_core=equality_core())
    paths = [Path(__file__)]
    for stem in ('fresh_ideal_distance7_mixed_replacement_20261017',
                 'fresh_distance8_ideal_corner_replacement_20261017'):
        paths += [ROOT/f'scripts/{stem}.py', ROOT/f'results/{stem}.json', ROOT/f'results/{stem}.md']
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    out.update(status='PASS_UNIVERSAL_CYCLE16_ANTIPODAL_HALF_COUPON_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        counts=dict(single_profiles=6,residual_profile_pairs=4,rectangle_vertices=16,
                    finite_equality_corners=3,half_coupon_corners=2,removed_antipodal_pairs=8),
        theorem='An actual quadrilateral having an antipodal marked pair in a 16-node varying affine A factor is paid by finite distance9/ideal distance7, or has two sealed finite A1 x A15 corners with the A15 mark its unique middle. A half-payment at each such whole-component pair pays the quadrilateral.',
        scope='Universal local transfer including arbitrary supporting spectators and nonsimple cusps. No global finite-vertex budget, ambient dimension bound, finite-volume construction, or proof-assistant claim.')
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    args = ap.parse_args()
    data = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],flush=True)
