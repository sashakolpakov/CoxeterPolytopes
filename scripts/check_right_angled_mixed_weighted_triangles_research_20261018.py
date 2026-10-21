#!/usr/bin/env python3
"""Tiny exact arithmetic certificate; the universal geometric proof is in the report.

No Lean, optimizer, numerical search, external data or large files are used.
Without --write, replay the saved compact certificate without modifying it.
"""
import argparse
from fractions import Fraction as Q
from itertools import combinations
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/right_angled_mixed_weighted_triangles_research_20261018.json'


def solve_rows(rows, unknowns):
    rows = [[Q(x) for x in row] for row in rows]
    for col in range(unknowns):
        pivot = next(i for i in range(col, len(rows)) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        divisor = rows[col][col]
        rows[col] = [x/divisor for x in rows[col]]
        for i in range(len(rows)):
            if i != col:
                factor = rows[i][col]
                rows[i] = [a-factor*b for a,b in zip(rows[i],rows[col])]
    return [[Q(0)]*unknowns+[-x for x in row[unknowns:]] for row in rows]


def ds_rows(n):
    width = n+2  # f,a1,...,a(n-1),c,constant.
    b = [[Q(0)]*width for _ in range(n+1)]
    b[0][0], b[0][n] = 1, 2**(n-1)
    for k in range(1,n-1):
        b[k][k], b[k][n] = 1, 2**(n-1-k)*comb(n-1,k)
    b[n-1][n-1], b[n-1][n], b[n][n+1] = 1, 1, 1
    h = [[sum(Q(comb(k,i)*(-1)**(k-i))*b[k][j] for k in range(i,n+1))
          for j in range(width)] for i in range(n+1)]
    count = (n+1)//2
    return solve_rows([[x-y for x,y in zip(h[i],h[n-i])] for i in range(count)],count)


def expression(row,n):
    names = ['f']+[f'a{k}' for k in range(1,n)]+['c','constant']
    return {name:str(value) for name,value in zip(names,row) if value}


def vertex_certificate(d):
    constraints=[]
    for i in range(d):constraints.append(([Q(int(i==j)) for j in range(d)],Q(0)))
    for i,j in combinations(range(d),2):
        constraints.append(([Q(int(k in (i,j))) for k in range(d)],Q(1)))
    vertices=set()
    for selected in combinations(constraints,d):
        # solve_rows returns the solution to Mx+B=0, hence append -rhs.
        try:rows=solve_rows([a+[-b] for a,b in selected],d)
        except StopIteration:continue
        point=tuple(row[-1] for row in rows)
        if min(point)>=0 and all(point[i]+point[j]<=1 for i,j in combinations(range(d),2)):
            vertices.add(point)
    expected={(Q(0),)*d}
    expected.update(tuple(Q(int(i==j)) for j in range(d)) for i in range(d))
    expected.update(tuple(Q(int(i in subset),2) for i in range(d))
                    for k in range(3,d+1) for subset in combinations(range(d),k))
    assert vertices==expected
    norms=[sum(x*x for x in point) for point in vertices]
    assert all(x<=1 for x in norms)
    return dict(dimension=d,walls=len(constraints),vertices=len(vertices),
        finite_vertices=sum(x<1 for x in norms),ideal_vertices=sum(x==1 for x in norms),
        active_constraint_choices_examined=comb(len(constraints),d),
        coordinates=[[str(x) for x in point] for point in sorted(vertices)])


def generate():
    records={}
    expected_residuals={
        11:{'a7':'-44/5','a8':'44','a9':'-132','a10':'308','constant':'-616'},
        12:{'a7':'-22/5','a8':'88/5','a9':'-198/5','a10':'308/5','a11':'-308/5'}}
    for n in (11,12):
        solved=ds_rows(n);f,a2=solved[0],solved[2]
        alpha,beta,U=comb(n,2),(n-1)*2**(n-2),Q(22,5)
        weight=(U*a2[n]-alpha*f[n])/beta
        assert weight==Q(2123,1280)
        residual=[alpha*x-U*y+(weight*beta if i==n else 0)
                  for i,(x,y) in enumerate(zip(f,a2))]
        assert expression(residual,n)==expected_residuals[n]
        low=(3*weight-U)/(U-3)*beta/alpha
        high=(2*U/5-weight)/(1-U/5)*beta/alpha
        A,B=3+weight-U,U-1-2*weight
        assert A>0 and B>0 and weight<2
        for a,b in ((0,3),(2,2),(3,1),(5,0)):
            assert Q(a)+weight*b-U>=A*a/3
        face_density=A/(2-weight)
        finite_density=alpha*A/(A+3*B)
        cube2=comb(n-1,2)*2**(n-3)
        b3_f=Q(5*(n-2),2)*alpha/U-Q(7*comb(n,3),4)
        b3_c=Q(5*(n-2),2)*weight*beta/U-Q(7*cube2,2)
        assert b3_c>0
        b3_finite=b3_f+b3_c/high
        k4_free=Q(n**3,27)
        assert 2*b3_finite>k4_free
        k4_fraction=(2*b3_finite-k4_free)/(comb(n,3)-k4_free)
        k4_free_integer=k4_free.numerator//k4_free.denominator
        # Delete at most one edge per original K4; each edge loses <=n-2 triangles.
        p4_finite=(2*b3_finite-k4_free_integer)/(5*(n-2))
        records[str(n)]=dict(finite_corner_incidence=alpha,ideal_corner_incidence=beta,
            ideal_weight=str(weight),strict_weighted_average_upper=str(U),
            ds_finite=expression(f,n),ds_two_faces=expression(a2,n),
            strict_residual=expression(residual,n),
            strict_finite_per_ideal_lower=str(low),strict_finite_per_ideal_upper=str(high),
            strict_triangle_face_fraction_lower=str(face_density),
            strict_triangles_per_finite_vertex_lower=str(finite_density),
            b3_lower_coefficients={'f':str(b3_f),'c':str(b3_c)},
            strict_b3_per_finite_vertex_lower=str(b3_finite),
            k4_free_triangle_upper=str(k4_free),
            strict_k4_corner_fraction_lower=str(k4_fraction),
            strict_p4_per_finite_vertex_lower=str(p4_finite))
    assert records['12']['strict_finite_per_ideal_lower']=='1474/21'
    assert records['12']['strict_finite_per_ideal_upper']=='1298/9'
    assert records['12']['strict_b3_per_finite_vertex_lower']=='2515/59'
    assert records['12']['strict_p4_per_finite_vertex_lower']=='627/1475'
    return dict(status='EXACT_ARITHMETIC_AND_LOCAL_VERTEX_CERTIFICATE',
        dimensions=records,local_polytopes=[vertex_certificate(3),vertex_certificate(4)],
        scope=dict(universal_geometry_in_written_proof=True,
            independent_geometric_review_complete=False,lean_compilation_performed=False,
            dimension_bound_improved=False,existence_claim=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert record==saved,'Saved compact certificate differs from exact replay'
    print('PASS exact DS11/12, fractions, polygon minima, complete P3/P4 vertices; geometric proof remains written; no Lean invoked')


if __name__=='__main__':main()
