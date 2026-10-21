#!/usr/bin/env python3
"""Independent bounded beta131 census in project Golden arithmetic.

No producer or producer arithmetic imports. Natural covector order, literal
E8 matrix, independently enlarged coefficient boxes and rational inverses.
"""
import argparse
import json
import sys
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import ceil, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from coxeter_polytopes.exact.golden import Golden as K

TAG = 'fresh_h9_beta131_pilot_20261016_review'
PRIMARY = 'results/fresh_h9_beta131_pilot_20261016.json'


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def enc(value):
    if isinstance(value, K): return value.to_json()
    if isinstance(value, F): return str(value)
    if isinstance(value, dict): return {k: enc(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [enc(v) for v in value]
    return value


def inverse(matrix):
    n = len(matrix)
    rows = [[F(v) for v in row] + [F(i == j) for j in range(n)]
            for i, row in enumerate(matrix)]
    for j in range(n):
        k = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[k] = rows[k], rows[j]
        pivot = rows[j][j]
        rows[j] = [v / pivot for v in rows[j]]
        for i in range(n):
            if i != j:
                factor = rows[i][j]
                rows[i] = [a - factor*b for a, b in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


def dot(a, b): return sum((x*y for x, y in zip(a, b)), K())
def mv(matrix, vector): return [dot(row, vector) for row in matrix]
def pair(a, b, matrix): return dot(a, mv(matrix, b))
def key(vector): return tuple((x.a, x.b) for x in vector)
def decode(vector): return [K.from_json(v) for v in vector]


def sqrt_ceiling(value):
    n = isqrt(value.numerator // value.denominator)
    return n + int(n*n < value)


def box(physical, conjugate):
    # z=a+b*phi, z'=a+b*(1-phi), sqrt(5)>11/5,
    # phi<17/10 and phi-1<7/10. Both coefficient bounds are strict.
    aa = ceil((F(17, 10)*conjugate + F(7, 10)*physical)/F(11, 5)) + 1
    bb = ceil((physical + conjugate)/F(11, 5)) + 1
    return aa, bb


def run():
    data = json.loads((ROOT / PRIMARY).read_text())
    for path, digest in data['sources'].items(): assert sha(path) == digest, path
    cap = F(data['height_bound'])
    assert cap > 0
    phi = K(0, 1)
    delta, beta = 1 + 12*phi, 12 + 13*phi
    assert phi*delta == beta and delta.norm() == -131 and beta.norm() == 131
    assert all(131 % d for d in range(2, 12))
    assert beta > 0 and beta.conjugate() > 0
    assert delta > 20 and -delta.conjugate() > 6
    assert phi < K(F(17, 10)) and phi-1 < K(F(7, 10)) and 2*phi-1 > K(F(11, 5))
    norms = [K(2), 2*beta]
    assert all(q < 67 and q.conjugate() < 8 for q in norms)
    edges = {(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)}
    cartan = [[F(2) if i == j else F(-1) if (min(i,j), max(i,j)) in edges else F(0)
               for j in range(8)] for i in range(8)]
    inv = inverse(cartan)
    assert all(v >= 0 and v.denominator == 1 for row in inv for v in row)
    # Exact positive pivots also check the literal positive block.
    block = [row[:] for row in cartan]
    determinant = F(1)
    while block:
        pivot = block[0][0]
        assert pivot > 0
        determinant *= pivot
        block = [[block[i][j]-block[i][0]*block[0][j]/pivot
                  for j in range(1,len(block))] for i in range(1,len(block))]
    assert determinant == 1
    form = [[K() for _ in range(10)] for _ in range(10)]
    form[0][0], form[0][9], form[9][0], form[9][9] = -2*phi, K(-1), K(-1), K(6)
    for i in range(8):
        for j in range(8): form[i+1][j+1] = K(cartan[i][j])
    assert [decode(row) for row in data['form']] == form
    assert decode(data['norm_classes']) == norms
    base = [K(6)] + [K()]*8 + [K(1)]
    assert decode(data['basepoint']) == base
    assert mv(form, base) == [-delta] + [K()]*9
    assert pair(base, base, form) == -6*delta
    assert (-2*phi).conjugate() > 0 and (-delta).conjugate() > 0

    # Separate bounds: x^2<6*67*T/20, (x')^2<8;
    # l^2<402*(1+T), (l')^2<48. They also bound every lambda coordinate.
    x_bound = sqrt_ceiling(F(6*67, 20)*cap)
    l_bound = sqrt_ceiling(F(402)*(1+cap))
    xa, xb = box(x_bound, 3)
    za, zb = box(l_bound, 7)
    controls = []
    for q in norms:
        for a in range(-xa, xa+1):
            for b in range(-xb, xb+1):
                x = K(a, b)
                height = delta*x*x/(6*q)
                if x > 0 and height <= cap and (6*q+delta*x*x).sign(1) >= 0:
                    controls.append((height, q, x))
    controls.sort(key=lambda row: row[0])
    assert len({h for h,q,x in controls}) == len(controls)
    assert controls == [(K.from_json(r['height']), K.from_json(r['norm']),
                         K.from_json(r['control'])) for r in data['levels']]
    inverse_energy = [[F(1)]+[F(0)]*8] + [[F(0)]+[v/6 for v in row] for row in cartan]
    forms = [inverse([row[:n] for row in inverse_energy[:n]]) for n in range(1,10)]
    pool = [K(a,b) for a in range(-za,za+1) for b in range(-zb,zb+1)]
    squares = {z:z*z for z in pool}
    normals = [[K(-int(i==j)) for i in range(10)] for j in range(1,9)]
    assert normals == [decode(r) for r in data['zero_simple_walls']]
    zero = [decode(r) for r in data['ordinary_zero_roots']]
    assert len({key(r) for r in zero}) == len(zero) == 240
    assert all(len(r)==8 and all(v.b==0 and v.is_integral for v in r)
               and pair(r,r,cartan)==2 for r in zero)
    total = Counter()
    records = []
    for index, (height, q, x) in enumerate(controls):
        target = 6*q+delta*x*x
        divisor = q/2
        divisible = [z for z in pool if (z/divisor).is_integral]
        signed = [z for z in divisible if (target-squares[z]).sign()>=0
                  and (target-squares[z]).sign(1)>=0 and ((x-z)/6).is_integral]
        positive = [z for z in divisible if z>=0 and (target-2*squares[z]).sign()>=0
                    and (target-2*squares[z]).sign(1)>=0]
        saved = data['levels'][index]
        assert set(signed) == set(decode(saved['signed_l_alphabet']))
        # Recurse over the weaker factor-two bound; separately check the
        # producer's valid factor-three bound (C_ii=2 gives 6/Jinv_ii=3).
        producer_positive = [z for z in positive if (target-3*squares[z]).sign()>=0
                             and (target-3*squares[z]).sign(1)>=0]
        assert set(producer_positive) == set(decode(saved['lambda_alphabet']))
        rows = [[r[0]-6*r[9]]+[6*r[i] for i in range(1,9)] for r in normals[8:]]
        bounds = [delta*x*r[0] for r in normals[8:]]
        assert all(z>=0 for row in rows for z in row[1:])
        found = []
        counts = Counter()
        def visit(prefix, partials):
            n = len(prefix)
            counts['nodes'] += 1
            if n:
                matrix = forms[n-1]
                energy = sum((prefix[i]*prefix[j]*matrix[i][j]
                              for i in range(n) for j in range(n)), K())
                if (target-energy).sign()<0 or (target-energy).sign(1)<0:
                    counts['energy_prunes'] += 1
                    return
                partials = [v+prefix[-1]*row[n-1] for v,row in zip(partials,rows)]
                if any(v>bound for v,bound in zip(partials,bounds)):
                    counts['wall_prunes'] += 1
                    return
            if n == 9:
                if energy != target: return
                vector = [x,*mv(inv,prefix[1:]),(x-prefix[0])/6]
                assert all(v.is_integral for v in vector)
                assert pair(vector,vector,form) == q
                assert all((v/divisor).is_integral for v in mv(form,vector))
                assert all(pair(vector,r,form)<=0 for r in normals)
                found.append(vector)
                return
            for z in (signed if n==0 else positive): visit(prefix+[z],partials)
        visit([], [K() for _ in rows])
        assert {key(r) for r in found} == {key(decode(r)) for r in saved['roots']}
        assert all(pair(r,s,form)<=0 for i,r in enumerate(found) for s in found[i+1:])
        normals.extend(found)
        total.update(counts)
        records.append(dict(index=index+1,height=height,norm=q,control=x,counts=dict(counts),roots=found))
        print('independent beta131 level',index+1,'/',len(controls),'walls',len(normals),flush=True)
    assert {key(r) for r in normals} == {key(decode(r)) for r in data['roots']}
    # Reorder only after independently reconstructing the entire prefix.
    normals = [decode(r) for r in data['roots']]
    gram = [[pair(r,s,form) for s in normals] for r in normals]
    assert gram == [decode(r) for r in data['gram']]
    center = decode(data['strict_center'])
    assert pair(center,center,form)<0 and pair(center,base,form)<0
    assert all(pair(center,r,form)<0 for r in normals)
    assert len(data['actual_facet_points']) == len(normals)
    for i,saved in enumerate(data['actual_facet_points']):
        point = decode(saved)
        assert pair(point,point,form)<0 and pair(point,base,form)<0
        assert all(pair(point,r,form).sign()==(0 if i==j else -1) for j,r in enumerate(normals))
    allowed = {K():2,K(1):3,K(2):4,1+phi:5,K(3):6,2+phi:10,K(4):'parallel'}
    relations = Counter()
    for i in range(len(normals)):
        for j in range(i):
            assert gram[i][j]<=0
            value=4*gram[i][j]*gram[i][j]/(gram[i][i]*gram[j][j])
            if value>4: kind='dotted'
            else:
                label=allowed[value]
                kind='parallel' if label=='parallel' else 'right' if label==2 else f'pi/{label}'
            relations[kind]+=1
    assert dict(relations)==data['pair_counts']
    paths=[PRIMARY,f'scripts/{TAG}.py','src/coxeter_polytopes/exact/golden.py']
    return enc(dict(status='PASS_INDEPENDENT_COMPLETE_BETA131_PREFIX',sources={p:sha(p) for p in paths},
                    height_bound=cap,norm_classes=norms,determinant_norm=-131,
                    coordinate_order=['l']+[f'lambda{i}' for i in range(1,9)],
                    independent_boxes=dict(control=[xa,xb],coordinates=[za,zb]),
                    controls=len(controls),walls=len(normals),counts=dict(total),levels=records,
                    pair_counts=dict(relations),
                    scope='Complete bounded prefix and actual-facet witnesses. No complete cone census, finite-volume or nonreflectivity conclusion. Ring, unit, norm-alphabet and zero-sector completeness require the accompanying mathematical argument.'))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    path=ROOT/'results'/f'{TAG}.json'
    if args.verify: assert result==json.loads(path.read_text())
    else: path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['controls'],result['walls'],result['counts'])
