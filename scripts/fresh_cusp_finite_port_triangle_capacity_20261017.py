#!/usr/bin/env python3
"""Exact counting for the cusp-fiber T1+T2 capacity; no optimizer."""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_cusp_finite_port_triangle_capacity_20261017'
SHARP = 'fresh_finite_edge_cusp_tradeoff_28over3_20261017'
SEVEN = 'fresh_finite_edge_cusp_tradeoff10_20261017'
EIGHT = 'fresh_finite_edge_cusp_tradeoff12_20261017'
COUNTER = 'fresh_finite_vertex_t2_port_counterexample_20261017'


def sha(p):
    return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def read(tag):
    j = json.loads((ROOT/f'results/{tag}.json').read_text())
    for p,h in j['sources'].items():
        assert sha(p) == h, p
    return j


def inverse(a):
    n = len(a)
    b = [list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for j in range(n):
        k = next(k for k in range(j,n) if b[k][j])
        b[j],b[k] = b[k],b[j]
        p = b[j][j]
        b[j] = [x/p for x in b[j]]
        for i in range(n):
            if i != j:
                p = b[i][j]
                b[i] = [x-p*y for x,y in zip(b[i],b[j])]
    return [row[n:] for row in b]


def rational(x):
    assert F(x[1]) == 0
    return F(x[0])


def run():
    sharp, seven, eight, counter = map(read, (SHARP,SEVEN,EIGHT,COUNTER))
    assert seven['seven_simplex_ideal_count_envelope'] == [0,1,2,3,8]
    assert eight['eight_simplex_case_count'] == 34
    assert sharp['uniform_finite_edge_coefficient'] == '28/3'
    # f counts finite vertices/ports; i counts ideal ports OTHER than the
    # distinguished cusp. The actual simplex has dimension k=f+i and i+1
    # ideal vertices. Unproved realizations remain in the safe envelope.
    cases = []
    for k in range(2,10):
        for f in range(1,k+1):
            i = k-f
            if k == 6 and f >= 2 and i+1 > 4:
                continue
            if k >= 7 and i+1 > 3:
                continue
            value = F(comb(f,2)+f*i, f)
            assert value <= 5 and value <= k-1
            cases.append(dict(dimension=k,finite_ports=f,other_ideal_ports=i,
                              T1=comb(f,2),T2=f*i,per_finite_port=str(value)))
    by_rank = []
    for r in range(1,65):
        relevant = [x for x in cases if x['dimension'] <= r+1]
        best = max(F(x['per_finite_port']) for x in relevant)
        assert best == min(r,5)
        by_rank.append(dict(factor_rank=r,maximum=str(best),
                            maximizing_types=[x for x in relevant if F(x['per_finite_port'])==best]))
    # Reconstruct the classical L24 simplex's rational inverse. Its causal
    # cone and complete affine cusp links were independently checked in the
    # bound28/3 package; this run checks the literal input and sharp counts.
    old = sharp['sharpness']
    edges = [(0,2),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(6,8),(8,9)]
    g = [[F(2 if i==j else 0) for j in range(10)] for i in range(10)]
    for i,j in edges:
        g[i][j]=g[j][i]=F(-1)
    assert g == [[rational(x) for x in row] for row in old['gram']]
    gi = inverse(g)
    assert gi == [[rational(x) for x in row] for row in old['inverse']]
    ideals = [i for i in range(10) if gi[i][i] == 0]
    finite = [i for i in range(10) if gi[i][i] < 0]
    assert ideals == [0,1,9] and len(finite)==7
    assert all(x <= 0 for row in gi for x in row)
    for link in old['cusp_links']:
        ids = link['incident_indices']
        kernel = [rational(x) for x in link['positive_affine_kernel']]
        assert all(x>0 for x in kernel)
        assert all(sum(g[i][j]*x for j,x in zip(ids,kernel))==0 for i in ids)
    triangles = [0]*4
    local = []
    for triangle in combinations(range(10),3):
        triangles[sum(v in ideals for v in triangle)] += 1
    for v in ideals:
        t1 = sum(v in u and sum(x in ideals for x in u)==1
                 for u in combinations(range(10),3))
        t2 = sum(v in u and sum(x in ideals for x in u)==2
                 for u in combinations(range(10),3))
        assert (t1,t2)==(21,14)
        assert t1+t2==7*5
        local.append(dict(cusp=v,factor_ranks=[8],finite_ended_ports=7,
                          T1=t1,T2=t2,capacity=35))
    assert triangles==[35,63,21,1]
    assert triangles[1]+2*triangles[2]==5*7*3
    # Read the literal T2 counterexample's boundary. Its two triangles at
    # the nonsimple cusp occupy DIFFERENT simplex factor fibers.
    ce = {tuple(sorted(e)) for e in counter['edges']}
    ideal_set = {x['vertex'] for x in counter['cusps']}
    cusp_rows=[]
    for v in sorted(ideal_set):
        ports=sorted(w for e in ce if v in e for w in e if w!=v)
        link_edges=[]; triangular=[]
        for label,face in counter['faces'].items():
            if v not in face:
                continue
            adjacent=sorted(w for w in ports if w in face)
            assert len(adjacent)==2
            link_edges.append(adjacent)
            if len(face)==3:
                triangular.append(dict(face=label,ports=adjacent,
                    ideal_vertices=sum(w in ideal_set for w in face)))
        if v == 'Iinf':
            assert len(ports)==4 and len(link_edges)==4
            assert {tuple(x) for x in link_edges}=={('A','I1'),('A','I2'),('F','I1'),('F','I2')}
            fibers=link_edges
            ranks=[1,1]
            assert len(triangular)==2 and triangular[0]['ports']!=triangular[1]['ports']
        else:
            assert len(ports)==3 and len(link_edges)==3
            fibers=[ports];ranks=[2]
        a=sum(w not in ideal_set for w in ports)
        t1=sum(t['ideal_vertices']==1 for t in triangular)
        t2=sum(t['ideal_vertices']==2 for t in triangular)
        assert t1+t2<=a*sum(min(r,5) for r in ranks)
        cusp_rows.append(dict(cusp=v,ports=ports,link_edges=link_edges,
            factor_fibers=fibers,factor_ranks=ranks,finite_ended_ports=a,
            triangular_link_edges=triangular,T1=t1,T2=t2))
    paths=[f'scripts/{TAG}.py']
    for t in (SHARP,SEVEN,EIGHT,COUNTER):
        paths += [f'scripts/{t}.py',f'results/{t}.json',f'results/{t}.md']
    paths += ['results/review_finite_edge_28over3_20261017.json',
              'results/review_finite_edge_28over3_20261017.md',
              'results/fresh_multi_ideal_triangle_scope_20261017.md',
              'results/fresh_finite_volume_simplex_truncation_20261017.md',
              'results/fresh_cusp_one_ideal_simplex_capacity_20261017.md']
    return dict(status='PASS_EXACT_CUSP_FINITE_PORT_MIXED_TRIANGLE_CAPACITY',
        sources={p:sha(p) for p in paths},block_types=cases,
        factor_rank_checks=by_rank,all_rank_argument='Every actual simplex has dimension<=9; each allowed block has (f-1)/2+i<=5 and <=k-1<=r. No rank extrapolation is used.',
        cusp_bound='T1(v)+T2(v)<=a_v*sum_i min(r_vi,5)',
        global_bound='T1+2*T2<=sum_v a_v*sum_i min(r_vi,5)',
        simple_cusp_corollary='T1+2*T2<=5*E_fi if every ideal vertex is simple',
        residual_budget='2*T3+Q >= (b_n-14*n/3)*v_f + sum_v a_v*(14/3-R5(v)) + alpha_n*sum_v E(L_v) + T1',
        sharpness=dict(diagram=24,dimension=9,finite_vertices=7,ideal_vertices=ideals,
            triangle_counts=triangles,finite_ideal_edges=21,cusps=local,
            global_left=105,global_capacity=105),
        counterexample_cusp_fibers=cusp_rows,
        scope='Finite-volume acute Coxeter polytopes. Maximal simplex blocks are within one cusp factor face, and only their finite-ended port sets are disjoint. No transitivity at an arbitrary finite vertex is assumed. Arbitrary other ideal vertices may be nonsimple. The ideal-count lemmas concern every actual simplex face, not only ambient-maximal ones. The coefficient5 is sharp; factor multiplicity cannot be deleted by this proof. No T3/quadrilateral upper bound, all-ideal exclusion, optimizer or ambient dimension improvement is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert result==json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],len(result['block_types']),'block types; sharp capacity',result['sharpness']['global_capacity'])
