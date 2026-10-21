#!/usr/bin/env python3
"""Exact checks for the accompanying all-rank narrow-cusp separation proof.

No producer arithmetic is imported. The all-rank statements are proved by the
closed formulas in the report; bounded matrix/row checks are regression checks,
not substitutes for those arguments. No full Nikulin <300 claim is made.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nikulin_narrow_20261015'


def gram(kind, r):
    g = [[F(0) for _ in range(r)] for _ in range(r)]
    for i in range(r):
        g[i][i] = 2
    if kind in ('A', 'B_scaled'):
        edges = [(i, i+1) for i in range(r-1)]
        if kind == 'B_scaled':
            # The short end root has norm 1, so all entries are rational.
            g[0][0] = 1
    else:
        assert kind == 'D' and r >= 4
        edges = [(0, 2), (1, 2)] + [(i, i+1) for i in range(2, r-1)]
    for i, j in edges:
        g[i][j] = g[j][i] = -1
    return g


def inverse_formula(kind, r):
    j = [[F(0) for _ in range(r)] for _ in range(r)]
    for a in range(r):
        for b in range(r):
            if kind == 'A':
                j[a][b] = F((min(a,b)+1)*(r-max(a,b)), r+1)
            elif kind == 'B_scaled':
                j[a][b] = F(r-max(a,b))
            elif max(a,b) < 2:
                j[a][b] = F(r if a == b else r-2, 4)
            elif min(a,b) < 2:
                j[a][b] = F(r-max(a,b), 2)
            else:
                j[a][b] = F(r-max(a,b))
    return j


def check_inverse(g, j):
    r = len(g)
    for a in range(r):
        row = [(b, x) for b, x in enumerate(g[a]) if x]
        for c in range(r):
            assert sum(x*j[b][c] for b, x in row) == (a == c)
    assert all(x > 0 for row in j for x in row)


def lorentz(x, y):
    return -x[0]*y[1]-x[1]*y[0]+sum(a*b for a,b in zip(x[2:],y[2:]))


def norm_two(sigma, z):
    assert sigma > 0
    return [(sum(a*a for a in z)-2)/(2*sigma), sigma] + z


def cycle_row(p):
    # Norm-two affine A_{n-1}: alpha0=c+en-e1, alphaj=ej-e(j+1).
    n = len(p)
    tails = [sum(p[j:]) for j in range(1,n)] + [F(0)]
    average = sum(tails)/n
    z = [average-a for a in tails]
    return norm_two(sum(p), z)


def path_row(p, delta):
    # All contact positions are internal. Delta is their common null weight:
    # 1 for affine C, 2 for affine B and D. End contacts are zero here.
    m = len(p)-1
    z = [-sum(p[j:]) for j in range(1,m+1)]
    if delta == 2:
        assert p[0] == p[1] == p[-1] == p[-2] == 0
        assert z[0] == z[1]
    return norm_two(delta*sum(p), z)


def check_row_pair(x,y):
    assert lorentz(x,x) == lorentz(y,y) == 2
    s,t = x[1], y[1]
    energy = sum((a/s-b/t)**2 for a,b in zip(x[2:],y[2:]))
    assert -lorentz(x,y) == s*t*energy/2-s/t-t/s
    return energy


def row_at(size, indices, numerators):
    p = [F(0)]*size
    for i,a in zip(indices,numerators):
        p[i] += F(a)
    return p


def run():
    dependencies = [
        'scripts/dimension_bound_prokhorov_20261015_proof_review.py',
        'results/dimension_bound_prokhorov_20261015_proof_review.json',
    ]
    audit = json.loads((ROOT/dependencies[1]).read_text())
    assert audit['status'] == 'PASS_PRIMARY_GEOMETRY_AND_SCOPE_AUDIT'
    for p,h in audit['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest() == h

    inverse_checks = {}
    for kind, lower in [('A',1),('B_scaled',2),('D',4)]:
        for r in range(lower,81):
            g,j = gram(kind,r), inverse_formula(kind,r)
            check_inverse(g,j)
            if kind == 'A' and r >= 9:
                n = r+1
                assert min(j[k][k] for k in range(2,r-2)) > 2
                assert min(j[k][k] for k in range(r)) > F(1,2)
            if kind == 'B_scaled' and r >= 5:
                # Restoring short-root norm2 divides inverse J00 by2.
                assert j[0][0]/2 > 2
                assert all(j[k][k] > 2 for k in range(1,r-2))
                assert j[-2][-2] == 2 and j[-1][-1] == 1
            if kind == 'D' and r >= 9:
                assert j[0][0] == j[1][1] == F(r,4) > 2
                assert all(j[k][k] > 2 for k in range(2,r-2))
                assert j[-2][-2] == 2 and j[-1][-1] == 1
        inverse_checks[kind] = 81-lower

    row_checks = 0
    for n in (16,31,64,135,136,137,200):
        i,j = 2, 2+n//2
        l1,l2 = j-i-4,n-(j-i)-4
        for a in range(-2,3):
            for b in range(-2,3):
                p = row_at(n,[(i+a)%n],[1])
                q = row_at(n,[(j+b)%n],[1])
                e = check_row_pair(cycle_row(p),cycle_row(q))
                assert e >= F(l1*l2,l1+l2)
                row_checks += 1
        p = row_at(n,[(i-2)%n,i,(i+2)%n],[1,3,2])
        q = row_at(n,[(j-2)%n,j,(j+2)%n],[2,1,4])
        assert check_row_pair(cycle_row(p),cycle_row(q)) >= F(l1*l2,l1+l2)
        row_checks += 1
    for m in (47,51,55,80):
        for delta in (1,2):
            i,j = 10,m-10
            for a in range(-2,3):
                for b in range(-2,3):
                    p = row_at(m+1,[i+a],[1])
                    q = row_at(m+1,[j+b],[1])
                    e = check_row_pair(path_row(p,delta),path_row(q,delta))
                    assert e >= F(j-i-4,delta*delta)
                    row_checks += 1

    # Exact scalar square identity for the two long cycle arcs.
    for l1,l2 in ((64,64),(64,65),(2,9)):
        for a in (F(-2),F(0),F(2,3),F(1),F(3)):
            assert l1*a*a+l2*(a-1)**2-F(l1*l2,l1+l2) == ((l1+l2)*a-l2)**2/(l1+l2)

    # Choosing deleted nodes while avoiding the narrow-facet node. These checks
    # include every narrow-node position at boundary and several larger ranks.
    selection_checks = 0
    thresholds = [('C',47,5,5,1),('B',51,9,5,2),('D',55,9,9,2)]
    records = []
    for family, m0,left,right,delta in thresholds:
        for m in (m0,m0+1,m0+13,120):
            for s in range(m+1):
                i,j = left,m-right
                if s == i: i += 1
                elif s == j: j -= 1
                assert s not in (i,j)
                assert i >= left and m-j >= right
                assert j-i-4 >= m-left-right-5 >= 32
                assert F(j-i-4,2)-2 >= 14
                selection_checks += 1
        records.append(dict(affine_family=family, rank_parameter='m=n-1',
                            first_excluded_m=m0, upper_ambient_dimension=m0,
                            interior_null_weight=delta,
                            guaranteed_energy=f'(m-{left+right+5})/{delta*delta}',
                            guaranteed_negative_pairing=f'(m-{left+right+9})/2'))
    for n in (136,137,200):
        for s in range(n):
            i,j=0,n//2
            if s in (i,j): i,j=1,1+n//2
            assert s not in (i,j)
            l1,l2=n//2-4,n-n//2-4
            assert l1 >= 64 and l2 >=64
            assert F(l1*l2,l1+l2) >=32
            selection_checks += 1

    sources = [f'scripts/{TAG}.py'] + dependencies
    return dict(
        status='PASS_EXPLICIT_CLASSICAL_IRREDUCIBLE_NARROW_CUSP_BOUNDS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        hypothesis='Finite-sided finite-volume real hyperbolic Coxeter polytope; Nikulin narrow facet contains an ideal vertex whose full cusp diagram is one irreducible affine classical diagram.',
        bounds=dict(A=135,B=51,C=47,D=55),
        exact_checks=dict(inverse_formulas=inverse_checks, norm_two_row_pairs=row_checks,
                          narrow_node_selections=selection_checks, path_bounds=records,
                          affine_A_boundary=dict(n=136,L1=64,L2=64,energy=32,negative_pairing=14)),
        proof_scope='Closed-form inverse, localization, Lorentz pairing and energy arguments in the report prove all ranks. Finite exact checks only validate their implementation and boundary indexing.',
        unresolved=['Products of affine cusp components','The case where the chosen narrow facet has only finite vertices','A complete proof of Nikulin dimension<300','Novelty or strongest-known-bound status'],
        primary_source=dict(url='https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf',
                            printed_pages=[655,657,658,659],pdf_pages_one_based=[757,759,760,761],
                            pdf_sha256='0375d371245a73ed40ebd95ffc66846b2990d54fd25a447f75ecd358f4d1fcff'),
    )


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--verify',action='store_true')
    args=parser.parse_args(); data=run(); out=ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert data == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
