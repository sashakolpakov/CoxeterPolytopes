#!/usr/bin/env python3
"""Independent Fraction-only E8 prefix and chamber-symmetry verification.

Dominant covectors are enumerated with a physical energy bound only; no
producer import, leading-Cartan Schur bounds, or producer recursion is used.
"""
import argparse
import json
from fractions import Fraction as F
from functools import cmp_to_key, lru_cache
from hashlib import sha256
from itertools import combinations, groupby
from pathlib import Path
import fresh_silver_f4_independent_verify as I

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_silver_e8_full_chamber_20261016_review'
PREFIX = 'fresh_silver_e8_h9_continuation_20261016'
SYMMETRY = PREFIX + '_symmetry'
Z, O = I.Z, I.O
EDGES = {(i, i+1) for i in range(6)} | {(2, 7)}
E8 = tuple(tuple(I.scale(O, 2) if i == j else I.neg(O) if (min(i,j),max(i,j)) in EDGES
                 else Z for j in range(8)) for i in range(8))
EI = I.inverse(E8)
G = tuple(tuple(I.neg(I.UNIT) if i == j == 0 else O if i == j == 9 else E8[i-1][j-1]
                if 1 <= i <= 8 and 1 <= j <= 8 else Z for j in range(10)) for i in range(10))


def positive(matrix):
    a = [list(row) for row in matrix];pivots = []
    while a:
        p = a[0][0];assert I.sign(p) > 0;pivots.append(p)
        a = [[I.sub(a[i][j], I.div(I.mul(a[i][0], a[0][j]), p))
              for j in range(1, len(a))] for i in range(1, len(a))]
    return pivots


@lru_cache(None)
def blocks(target, norm):
    alphabet = tuple(w for w in I.bounded((I.scale(target, 2), I.scale(target, 2)))
                     if I.sign(w) >= 0 and I.integral(I.div(I.scale(w, 2), norm)))
    assert all(I.sign(x) > 0 for row in EI for x in row)
    found = [];visits = 0
    def visit(prefix, energy):
        nonlocal visits
        visits += 1
        if I.sign(I.sub(target, energy)) < 0:
            return
        k = len(prefix)
        if k == 8:
            if I.sign(I.conj(I.sub(target, energy))) < 0:
                return
            vector = I.mv(EI, prefix)
            assert all(map(I.integral, vector)) and I.pair(vector, E8, vector) == energy
            found.append((energy, vector))
            return
        cross = I.total(I.mul(EI[k][j], prefix[j]) for j in range(k))
        for w in alphabet:
            increment = I.add(I.mul(EI[k][k], I.mul(w,w)), I.scale(I.mul(w,cross),2))
            visit(prefix+(w,), I.add(energy,increment))
    visit((), Z)
    return tuple(found), visits, len(alphabet)


def candidates(x, norm):
    target = I.add(norm, I.mul(I.UNIT, I.mul(x,x)))
    choices, visits, alphabet = blocks(target, norm)
    scalars = [z for z in I.bounded((target,target)) if I.sign(z) >= 0
               and I.integral(I.div(I.scale(z,2),norm))]
    squares = {I.mul(z,z):z for z in scalars};assert len(squares) == len(scalars)
    out = []
    for energy, spatial in choices:
        remaining = I.sub(target,energy)
        if remaining not in squares:
            continue
        vector = (x,*spatial,squares[remaining]);cov = I.mv(G,vector)
        assert I.dot(vector,cov) == norm and all(map(I.integral,vector))
        assert all(I.integral(I.div(I.scale(c,2),norm)) for c in cov)
        out.append(vector)
    return out, visits, alphabet


def compare(left,right):
    for a,b in zip(left,right):
        s = I.sign(I.sub(a,b))
        if s:return s
    return 0


def prefix(data):
    assert I.matrix(data['gram']) == G and I.determinant(E8) == O
    assert I.determinant(G) == I.neg(I.UNIT) and I.integral(I.inv(I.UNIT))
    positive(E8);positive(tuple(tuple(I.conj(x) for x in row) for row in G))
    pool = (O,I.ETA,I.scale(O,2));assert tuple(map(I.parse,data['norm_classes'])) == pool
    zero = []
    for q, saved in zip(pool,data['zero_height_dominant_census']):
        roots,visits,alphabet = candidates(Z,q)
        assert set(roots) == set(tuple(map(I.parse,r)) for r in saved['roots'])
        zero.append(dict(norm=I.enc(q),roots=len(roots),visits=visits,alphabet=alphabet))
    cap = I.parse(data['height']);assert cap == I.scale(O,100)
    controls = []
    for q in pool:
        for x in I.bounded((I.div(I.mul(cap,q),I.UNIT), I.neg(I.div(q,I.UNIT)))):
            if I.sign(x)>0 and I.integral(I.div(I.scale(I.mul(I.UNIT,x),2),q)):
                controls.append((I.div(I.mul(I.UNIT,I.mul(x,x)),q),q,x))
    controls.sort(key=cmp_to_key(compare))
    assert len(controls) == data['time_controls'] == 29
    assert controls == [(I.parse(r['height']),I.parse(r['norm']),I.parse(r['time_coordinate'])) for r in data['levels']]
    walls = [tuple(I.neg(O) if i==j else Z for i in range(10)) for j in range(1,10)]
    keys = {I.oriented_key(v) for v in walls};records = [];visits_total = 0
    saved_walls = I.matrix(data['all_walls'])
    for h, packet in groupby(controls,key=lambda r:r[0]):
        added = []
        for _,q,x in packet:
            roots,visits,alphabet = candidates(x,q);visits_total += visits
            selected = duplicates = 0
            for r in roots:
                if any(I.sign(I.pair(r,G,w))>0 for w in walls):continue
                key = I.oriented_key(r)
                if key in keys:duplicates += 1;continue
                keys.add(key);added.append(r);selected += 1
            saved = data['levels'][len(records)]
            assert (len(roots),selected,duplicates) == (saved['candidates'],saved['roots'],saved['parallel_duplicates'])
            records.append(dict(height=I.enc(h),norm=I.enc(q),control=I.enc(x),
                                candidates=len(roots),selected=selected,duplicates=duplicates,
                                physical_bound_visits=visits,covector_alphabet=alphabet))
        assert all(I.sign(I.pair(r,G,s))<=0 for r,s in combinations(added,2))
        walls.extend(added)
        expected = [r for r in saved_walls if I.sign(I.sub(
            I.div(I.mul(I.UNIT,I.mul(r[0],r[0])),I.pair(r,G,r)),h))<=0]
        assert {I.oriented_key(v) for v in walls} == {I.oriented_key(v) for v in expected}
        print('independent E8 level',len(records),'walls',len(walls),flush=True)
    assert len(walls) == data['wall_count'] == 48
    return dict(zero_height=zero,time_controls=29,walls=48,physical_bound_visits=visits_total,levels=records)


def witness(data, proof):
    normals = I.matrix(data['all_walls']);cov = [I.mv(G,r) for r in normals]
    gram = tuple(tuple(I.dot(r,v) for v in cov) for r in normals)
    assert gram == I.matrix(proof['all_pair_gram'])
    center = tuple(map(I.parse,proof['strict_center']));base = (O,)+(Z,)*9
    assert I.sign(I.pair(center,G,center))<0 and I.sign(I.pair(center,G,base))<0
    assert all(I.sign(I.dot(center,v))<0 for v in cov)
    for i,r in enumerate(normals):
        assert I.sign(gram[i][i])>0 and all(I.sign(gram[i][j])<=0 for j in range(i))
        assert all(I.integral(I.div(I.scale(v,2),gram[i][i])) for v in cov[i])
        point = tuple(map(I.parse,proof['actual_facet_points'][i]))
        assert I.sign(I.pair(point,G,point))<0 and I.sign(I.pair(point,G,base))<0
        assert all(I.sign(I.dot(point,v)) == (0 if i==j else -1) for j,v in enumerate(cov))
    sg = I.matrix(proof['spherical_gram']);positive(sg);inverse = I.inverse(sg)
    assert all(I.sign(v)>=0 for row in inverse for v in row)
    assert sg == tuple(tuple(G[i][j] for j in range(1,10)) for i in range(1,10))
    for row in proof['corners']:
        ids = row['ordered_indices'];point = tuple(map(I.parse,row['point']))
        assert tuple(tuple(gram[i][j] for j in ids) for i in ids) == sg
        assert I.sign(I.pair(point,G,point))<0 and I.sign(I.pair(point,G,base))<0
        pairings = [I.dot(point,v) for v in cov]
        assert [i for i,v in enumerate(pairings) if v==Z] == sorted(ids)
        assert all(I.sign(v)<=0 for v in pairings)
        a = I.div(I.pair(base,G,point),I.pair(point,G,point))
        b = I.mv(inverse,[I.dot(base,cov[i]) for i in ids])
        assert I.sign(a)>0 and all(I.sign(v)<=0 for v in b)
        assert base == tuple(I.add(I.mul(a,point[k]), I.total(I.mul(b[j],normals[i][k])
                           for j,i in enumerate(ids))) for k in range(10))
        assert a == I.parse(row['basepoint_decomposition']['vertex_coefficient'])
        assert b == tuple(map(I.parse,row['basepoint_decomposition']['normal_coefficients']))
    maps = []
    for row in proof['generators']:
        matrix = I.matrix(row['matrix']);inverse = I.inverse(matrix)
        assert inverse == I.matrix(row['inverse']) and I.mm(matrix,inverse) == I.eye(10)
        assert all(I.integral(v) for r in matrix+inverse for v in r)
        assert I.mm(I.tr(matrix),I.mm(G,matrix)) == G
        assert all(I.mv(matrix,normals[i]) == normals[j] for i,j in zip(row['source_indices'],row['target_indices']))
        target = next(c for c in proof['corners'] if c['ordered_indices']==row['target_indices'])
        scale = I.parse(row['positive_corner_scale']);assert I.sign(scale)>0
        assert I.mv(matrix,base) == tuple(I.mul(scale,I.parse(v)) for v in target['point'])
        maps.append(matrix)
    product = I.mm(maps[1],maps[0]);matrix = I.mm(product,product)
    assert product == I.matrix(proof['product_matrix']) and matrix == I.matrix(proof['matrix'])
    assert I.inverse(matrix) == I.matrix(proof['inverse']) and I.determinant(matrix) == O
    trace = I.total(matrix[i][i] for i in range(10))
    assert trace == I.parse(proof['trace']) == (F(10),F(4))
    assert I.sign(I.sub(trace,I.scale(O,10)))>0
    escape = tuple(map(I.parse,proof['strict_escape']['point']))
    assert I.sign(I.pair(escape,G,escape))>0 and I.sign(I.pair(escape,G,base))<0
    assert all(I.sign(I.dot(escape,v))<0 for v in cov)
    return dict(integral_reflections=48,actual_facets=48,persistent_outward_corners=3,
                integral_future_generators=2,stabilizer_word=[1,0,1,0],determinant=I.enc(O),
                trace=I.enc(trace),trace_minus_dimension=I.enc(I.sub(trace,I.scale(O,10))),
                strict_prefix_escape=True)


def run():
    data = I.read(f'results/{PREFIX}.json');proof = I.read(f'results/{SYMMETRY}.json')
    enumeration = prefix(data);checked_witness = witness(data,proof)
    paths = [f'scripts/{TAG}.py','scripts/fresh_silver_f4_independent_verify.py',
             f'results/{PREFIX}.json',f'results/{SYMMETRY}.json',
             'scripts/fresh_silver_e8_h9_20261015.py']
    visited = set()
    def check(path):
        if path in visited:return
        visited.add(path)
        if path.endswith('.json'):
            for source,digest in I.read(path).get('sources',{}).items():
                assert sha256((ROOT/source).read_bytes()).hexdigest()==digest,(path,source)
                check(source)
    for path in paths:check(path)
    return dict(status='PASS_INDEPENDENT_COMPLETE_E8_PREFIX_AND_INFINITE_ORDER_CHAMBER_SYMMETRY',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                independent_prefix=enumeration,independent_witness=checked_witness,
                recursive_manifest_files=len(visited),
                scope='Independent Fraction field arithmetic, all29 height controls, full three-norm dominant-covector prefix,48 actual facets,three persistent complete E8+A1 corners,two integral future sector maps and infinite-order squared product. The unimodular primitive-norm argument and corner-persistence theorem are explicitly audited mathematical inputs. The full finite-corner/edge census is not independently replayed. The literal lattice full reflection group is nonreflective; no other host or construction dimension is excluded.')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['independent_prefix']['physical_bound_visits'])
