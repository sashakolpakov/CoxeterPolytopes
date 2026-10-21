#!/usr/bin/env python3
"""Small exact certificate for a forced noncrystallographic ideal link.

The four-wall Gram forces a cusp and an I2(5) subgroup on its horosphere.
The accompanying report supplies the Bieberbach finite-volume argument.
No search/producer or numerical module is imported.
"""
import argparse
import copy
import json
from pathlib import Path
from scripts.fresh_h9_quartic_completion_20261013_verify import read, matrix, load_old
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import (
    K, O, Z, P, enc, parse, dot, matvec, pairing)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = 'data/enumerations/fresh_h9_quartic_completion_20261013_second_geometry.json'
OUT = ROOT/'data/enumerations/fresh_h9_fixed_core_cusp_obstruction_20261013.json'


def construct():
    geometry = read(SOURCE)
    assert 31 in geometry['cases'][0]['retained_old_caps']
    capsource = read(geometry['source'])
    base = read(capsource['base_source'])
    G,I = matrix(base['seed_gram']),matrix(base['seed_gram_inverse'])
    cy = tuple(map(parse,base['strict_center_y']))
    assert inertia_K(G) == (9,1,0)
    for i,row in enumerate(G):
        assert matvec(K,I,row) == tuple(O if i==j else Z for j in range(10))
    center = tuple(K.neg(x) for x in matvec(K,I,cy))
    assert K.sign(pairing(K,center,G,center)) < 0
    assert all(K.sign(x)>0 for x in cy)
    old = load_old(G)[31]
    assert old['d'] == Z and old['q'] == K.rational(2)
    assert all(x == Z for x in old['b'])
    unit = lambda i: tuple(O if i==j else Z for j in range(10))
    normals = (old['a'],unit(8),unit(6),unit(7))
    gram = tuple(tuple(pairing(K,x,G,y) for y in normals) for x in normals)
    two,minus_two,minus_phi = K.rational(2),K.rational(-2),K.neg(P)
    expected = ((two,minus_two,Z,Z),(minus_two,two,Z,Z),
                (Z,Z,two,minus_phi),(Z,Z,minus_phi,two))
    assert gram == expected
    center_pairs = tuple(pairing(K,center,G,x) for x in normals)
    assert all(K.sign(x)<0 for x in center_pairs)
    v = tuple(K.add(x,y) for x,y in zip(normals[0],normals[1],strict=True))
    assert any(x != Z for x in v)
    assert pairing(K,v,G,v) == Z
    assert K.sign(pairing(K,center,G,v)) < 0
    assert all(pairing(K,v,G,x)==Z for x in normals)

    # Reflection matrices in the ten-dimensional Lorentz space, norm(r)=2.
    reflections = []
    for normal in normals[2:]:
        covector = matvec(K,G,normal)
        reflections.append(tuple(tuple(K.sub(O if i==j else Z,K.mul(normal[i],covector[j]))
                                       for j in range(10)) for i in range(10)))
    product = tuple(tuple(dot(K,row,tuple(reflections[1][k][j] for k in range(10)))
                          for j in range(10)) for row in reflections[0])
    assert matvec(K,product,v) == v
    trace = K.zero
    for i in range(10): trace = K.add(trace,product[i][i])
    assert trace == K.add(K.rational(7),P)
    quotient_trace = K.sub(trace,K.rational(2))
    assert quotient_trace == K.add(K.rational(5),P)
    # phi is irrational: t²-t-1 has nonsquare rational discriminant5.
    assert quotient_trace[1] == 1
    return {'status':'PASS_EXACT_FOUR_WALL_CUSP_OBSTRUCTION',
            'base_source':capsource['base_source'],'geometry_source':SOURCE,
            'wall_names':['old31','seed8','seed6','seed7'],
            'normal_coefficients':[[enc(x) for x in row] for row in normals],
            'gram':[[enc(x) for x in row] for row in gram],
            'seed_inertia':[9,1,0],
            'strict_center_coefficients':list(map(enc,center)),
            'strict_center_square':enc(pairing(K,center,G,center)),
            'center_wall_pairings':list(map(enc,center_pairs)),
            'positive_null_sum_coefficients':[1,1],
            'null_vector_coefficients':list(map(enc,v)),
            'null_square':enc(Z),'null_center_pairing':enc(pairing(K,center,G,v)),
            'lorentz_reflection_product_trace':enc(trace),
            'horosphere_linear_dimension':8,'horosphere_linear_trace':enc(quotient_trace),
            'scope':'No finite-volume Coxeter completion by adding acute-compatible walls '
                    'while retaining these four facet mirrors, hence none for the fixed '
                    '44-,45-,46-wall cores. Does not obstruct replacing inherited walls '
                    'or finite-volume completion of the ten-wall seed alone.'}


def verify(saved):
    expected = construct()
    assert saved == expected
    return expected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args = parser.parse_args()
    if args.verify:
        result = verify(json.loads(OUT.read_text()))
    else:
        result = construct()
        OUT.write_text(json.dumps(result,indent=2)+'\n')
    for kind in ('gram','trace'):
        bad = copy.deepcopy(result)
        if kind == 'gram': bad['gram'][0][1]['a']='-1'
        else: bad['horosphere_linear_trace']['b']='0'
        try: verify(bad)
        except AssertionError: pass
        else: raise AssertionError('Forged certificate accepted')
    print({'status':result['status'],'horosphere_trace':'5+phi',
           'forged_gram_rejected':True,'forged_trace_rejected':True})


if __name__ == '__main__':
    main()
