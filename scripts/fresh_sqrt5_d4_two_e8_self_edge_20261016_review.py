#!/usr/bin/env python3
"""Independent project-field review of the V-arithmetic local edge witness."""
import argparse
import json
from fractions import Fraction as Q
import fresh_h9_beta131_h300_symmetry_20261016_review as I

K, ROOT = I.K, I.ROOT
TAG = 'fresh_sqrt5_d4_two_e8_self_edge_20261016_review'
PRIMARY = 'results/fresh_sqrt5_d4_two_e8_self_edge_20261016.json'


def dec(v):
    return [K(Q(a),Q(b)) for a,b in v]


def ldl(a):
    a = [row[:] for row in a]
    out = []
    while a:
        p = a[0][0]
        assert p
        out.append(p)
        a = [[a[i][j]-a[i][0]*a[0][j]/p for j in range(1,len(a))]
             for i in range(1,len(a))]
    return out


def run():
    saved = I.read(PRIMARY)
    phi = K(0,1)
    g = [[K(2*int(i==j)) for j in range(22)] for i in range(22)]
    for i,j in ((0,1),(0,2),(0,3),(1,4),(1,5)):
        g[i][j] = g[j][i] = K(-1)
    g[4][5] = g[5][4] = -phi
    for offset in (6,14):
        for i,j in ((0,1),(0,2),(2,3),(0,4),(4,5),(5,6),(6,7)):
            g[offset+i][offset+j] = g[offset+j][offset+i] = K(-1)
    assert g == [dec(row) for row in saved['root_gram']]
    inv,det = I.inverse(g)
    assert det == -4*(1+3*phi)
    physical = ldl(g)
    conjugate = ldl([[x.conjugate() for x in row] for row in g])
    assert sum(x>0 for x in physical) == 21 and sum(x<0 for x in physical) == 1
    assert all(x>0 for x in conjugate)
    assert all(x.is_integral for row in g for x in row)
    assert g[1][4]*g[4][5]*g[5][1] == -phi
    def pair(a,b):
        return I.dot(a,I.mv(g,b))
    endpoints = []
    for omitted,raw in zip((4,5),saved['endpoints']):
        v = [-row[omitted] for row in inv]
        vals = I.mv(g,v)
        assert vals == [-K(int(i==omitted)) for i in range(22)]
        assert v == dec(raw['coefficients'])
        assert pair(v,v) == (4-3*phi)/5 and pair(v,v) < 0
        active = [i for i in range(22) if i!=omitted]
        assert all(x>0 for x in ldl([[g[i][j] for j in active] for i in active]))
        endpoints.append(v)
    cross = pair(*endpoints)
    assert cross == (1-2*phi)/5 and cross < 0
    assert cross/pair(endpoints[0],endpoints[0]) == 1+phi
    u = [-x for x in I.mv(inv,[K(Q(1,10))]*4+[K(1)]*2+[K()]*16)]
    h = [-x for x in I.mv(inv,[K()]*6+[K(1)]*16)]
    uu,hh = pair(u,u),pair(h,h)
    assert uu == (176-169*phi)/50 and hh == 1240 and pair(u,h) == 0
    assert uu < 0 and hh > 0
    t2 = -uu/hh
    assert t2 > 0 and t2 == (169*phi-176)/62000
    assert uu+t2*hh == 0
    gu,gh = I.mv(g,u),I.mv(g,h)
    assert all(x<0 for x in gu[:6]) and all(x==0 for x in gu[6:])
    assert gh[:6] == [K()]*6 and gh[6:] == [K(-1)]*16
    escape = saved['strict_null_escape']
    assert u == dec(escape['u']) and h == dec(escape['h'])
    assert t2 == dec([escape['kappa']])[0]
    w = [a+Q(1,1000)*b for a,b in zip(u,h)]
    assert pair(w,w) < 0 and pair(w,u) < 0 and all(x<0 for x in I.mv(g,w))
    assert w == dec(saved['strict_interior_timelike_point']['coefficients'])
    paths = [f'scripts/{TAG}.py',PRIMARY,'src/coxeter_polytopes/exact/golden.py',
             'scripts/fresh_h9_beta131_h300_symmetry_20261016_review.py']
    return I.enc(dict(status='PASS_INDEPENDENT_ADMISSIBLE_SELF_EDGE_AND_STRICT_REAL_NULL_ESCAPE',
                      sources={p:I.sha(p) for p in paths},determinant=det,
                      physical_pivots=physical,conjugate_pivots=conjugate,
                      endpoints=endpoints,endpoint_pairing=cross,cosh_length=1+phi,
                      escape_u=u,escape_h=h,escape_t_squared=t2,
                      scope='Independent exact local V-arithmetic edge and infinite-volume certificate. The null vector uses a positive real radical, not a K-rational isotropic vector. No finite-volume or arithmetic-lattice claim, and no dimension-bound improvement.'))


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert result==json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
