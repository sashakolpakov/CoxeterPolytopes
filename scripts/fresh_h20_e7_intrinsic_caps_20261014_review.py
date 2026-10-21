#!/usr/bin/env python3
"""Independent sparse-field replay of the85wall replacement and source transfer.

Uses sparse square-class multiplication/interval signs, unlike the producer's
recursive quadratic arithmetic. Source vertex completeness is inherited from
the separate sparse-support source review, not recounted as a new census.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_h20_e8_interior_caps_20261014_exact as K

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_intrinsic_caps_20261014_review'
PRIMARY = 'results/fresh_h20_e7_intrinsic_caps_20261014.json'
BASE = 'results/fresh_h20_e8_interior_caps_20261014_exact.json'
FACE = 'results/fresh_h20_e7_face_20261014.json'
SOURCE = 'results/fresh_h13_a11a2_completion_20261014.json'
SOURCE_REVIEW = 'results/fresh_h13_a11a2_completion_20261014_review.json'


def field(x): return tuple(map(F, x))
def vector(x): return tuple(map(field, x))
def total(xs):
    out = K.Z
    for x in xs: out = K.add(out, x)
    return out


def run():
    d, old, face, source, review = [json.loads((ROOT/p).read_text()) for p in [PRIMARY, BASE, FACE, SOURCE, SOURCE_REVIEW]]
    for doc in [d, old, face, source, review]:
        for p, h in doc['sources'].items(): assert sha256((ROOT/p).read_bytes()).hexdigest() == h
    Q = [[F(x) for x in r] for r in d['form']]
    assert Q == [[F(x) for x in r] for r in old['form']]
    sparse = [[(i,q) for i,q in enumerate(row) if q] for row in Q]
    cache = {}
    def dot(v, w):
        if w not in cache:
            cache[w] = [tuple(sum(q*w[i][k] for i,q in row) for k in range(8)) for row in sparse]
        return total(K.mul(a,b) for a,b in zip(v,cache[w]))
    def linear(vectors, coefficients):
        return tuple(total(K.mul(v[i],a) for v,a in zip(vectors,coefficients)) for i in range(21))
    N0, N = list(map(vector, old['normals'])), list(map(vector, d['normals']))
    retained = [i for i in range(87) if i not in [17,26,75,76]]
    assert d['retained_old87_walls'] == retained and d['removed_old87_walls'] == [17,26,75,76]
    assert N[:83] == [N0[i] for i in retained]
    T = vector(face['intrinsic_normalization']['timelike_basis'])
    A = face['face_walls']
    J = K.rat.inv([[F(x) for x in row] for row in face['face_gram']])
    weights = [linear([N0[i] for i in A],list(map(K.c,row))) for row in J]
    sqrt2 = (F(0),F(1))+(F(0),)*6
    capR = linear([weights[10],weights[11],T],[K.c(-1),K.scale(sqrt2,-1),K.scale(sqrt2,F(1,4))])
    capS = linear([weights[2],T],[K.c(-1),K.scale(sqrt2,F(1,4))])
    assert N[83:] == [capR,capS]
    norms = [dot(v,v) for v in N]
    assert norms == [K.c(F(x)) for x in d['norms']]
    kinds = Counter()
    new_kinds = Counter()
    for i, v in enumerate(N):
        for j, w in enumerate(N[:i]):
            pair = dot(v,w)
            assert K.sign(pair) <= 0
            sq, prod = K.mul(pair,pair), K.mul(norms[i],norms[j])
            if pair == K.Z: kind='right'
            elif sq == K.scale(prod,F(1,4)): kind='pi3'
            elif sq == K.scale(prod,F(1,2)): kind='pi4'
            elif sq == K.scale(prod,F(3,4)): kind='pi6'
            elif sq == prod: kind='parallel'
            else:
                assert K.sign(K.add(sq,K.scale(prod,-1))) > 0
                kind='dotted'
            kinds[kind] += 1
            if i>=83: new_kinds[kind] += 1
    assert dict(new_kinds) == d['new_pair_counts']
    H = vector(d['strict_center'])
    assert K.sign(dot(H,H)) < 0 and all(K.sign(dot(H,r)) < 0 for r in N)
    for i,row in enumerate(d['facet_points']):
        p = vector(row)
        assert K.sign(dot(p,p)) < 0 and K.sign(dot(p,H)) < 0
        assert all(dot(p,r)==K.Z if i==j else K.sign(dot(p,r))<0 for j,r in enumerate(N))
    # Explicit source-frame isometry: source V -> sqrt2*T, its first11
    # weights unchanged, its A2 weights -> swapped target weights/sqrt2.
    transferred_weights = weights[:11]+[
        linear([weights[12]],[K.scale(sqrt2,F(1,2))]),
        linear([weights[11]],[K.scale(sqrt2,F(1,2))])]
    source_time = linear([T],[sqrt2])
    saved = {vector(row['ambient_point']):row for row in d['face_finite_vertices']+d['face_ideal_vertices']}
    seen = set()
    for vertex in source['projective_vertices']:
        z = list(map(F,vertex['gaps']))
        p = linear([source_time]+transferred_weights,[K.c(1)]+[K.c(-q) for q in z])
        assert p in saved and dot(p,p)==K.c(F(vertex['norm']))
        rec=saved[p]; seen.add(p)
        incident=[i for i,r in enumerate(N) if dot(p,r)==K.Z]
        assert incident==rec['ambient_incident']
        assert all(K.sign(dot(p,r))<=0 for r in N)
        if F(vertex['norm'])==0:
            recession=vector(rec['ambient_spherical_recession'])
            assert dot(recession,p)==K.Z and K.sign(dot(recession,recession))>0
            assert all(K.sign(dot(recession,N[i]))<=0 for i in incident)
        else:
            assert len(incident)==20
    assert seen==set(saved) and len(seen)==46
    # Other retained walls restrict to future timelike normals, and thus
    # add no boundaries on the complete transferred hyperbolic face.
    other=[i for i in retained if i not in A+face['E7_support']]
    assert other==d['causal_nonfacet_old87_walls'] and len(other)==63
    for i in other:
        p=vector(face['projected_walls'][i])
        assert K.sign(dot(p,p))<0 and K.sign(dot(p,T))<0
    E=vector(d['strict_escape']['point'])
    assert E==vector(old['strict_escape']['point'])
    assert K.sign(dot(E,E))>0 and K.sign(dot(E,H))<0
    assert all(K.sign(dot(E,r))<0 for r in N)
    paths=[PRIMARY,BASE,FACE,SOURCE,SOURCE_REVIEW,f'scripts/{TAG}.py',
           'scripts/fresh_h20_e8_interior_caps_20261014_exact.py',
           'scripts/fresh_h20_coupled_lift_20261014.py']
    return dict(status='PASS_INDEPENDENT_85_WALL_REPLACEMENT_AND_COMPLETE_FACE_TRANSFER',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        all_pairs_checked=sum(kinds.values()),pair_counts=dict(kinds),new_pair_counts=dict(new_kinds),
        actual_facets=85,transferred_finite_vertices=44,transferred_ideal_rays=2,
        redundant_timelike_face_projections=63,strict_escape_survives=True,
        scope='Independent full85 pair/facet/strict-center replay and exact isometry of every complete-source vertex. Intrinsic completeness inherits the separate sparse-support source review. Both ambient cusp links remain unbounded. Does not repeat all16 rejected transfer attempts or claim a full ambient vertex census.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify: assert json.loads(path.read_text())==result
    else: path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['pair_counts'])
