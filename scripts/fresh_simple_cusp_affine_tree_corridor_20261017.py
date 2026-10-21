#!/usr/bin/env python3
"""Exact central-corridor branch for long affine B/C/D quadrilateral corners."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_simple_cusp_affine_tree_corridor_20261017'
BASE='fresh_simple_cusp_affine_a_quadrilateral_20261017'
REVIEW='review_simple_cusp_affine_a_quadrilateral_20261017'
spec=importlib.util.spec_from_file_location('affine_a_exact',ROOT/f'scripts/{BASE}.py')
M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)
Q,SQ,Z=M.Q,M.SQ,M.Z

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def endblock(kind,r):
    if kind=='B':
        edges=[(i,i+1,Q(1) if i<r-2 else SQ) for i in range(r-1)]
        column=[Q(1)]*(r-1)+[SQ/2]
    else:
        edges=[(i,i+1,Q(1)) for i in range(r-3)]+[(r-3,r-2,Q(1)),(r-3,r-1,Q(1))]
        column=[Q(1)]*(r-2)+[Q(F(1,2)),Q(F(1,2))]
    A=[[Q(2) if i==j else Z for j in range(r)] for i in range(r)]
    for i,j,w in edges:A[i][j]=A[j][i]=-w
    result=[sum((A[i][j]*column[j] for j in range(r)),Z) for i in range(r)]
    assert result==[Q(1)]+[Z]*(r-1)
    assert column[0]==1 and all(t.sign()>0 for t in column)
    return dict(type=kind,order=r,inverse_port_diagonal='1',
                checked_column=[t.wire() for t in column])

def run():
    b=json.loads((ROOT/f'results/{BASE}.json').read_text())
    assert b['status']=='PASS_ALL_RANK_AFFINE_A_QUADRILATERAL_PAYMENT'
    rv=json.loads((ROOT/f'results/{REVIEW}.json').read_text())
    assert rv['status']=='PASS_INDEPENDENT_AFFINE_A_RECTANGLE_AND_KERNEL_REVIEW'
    for data in (b,rv):
        for p,h in data['sources'].items():assert sha(p)==h,p
    rows=[]
    for name,c0,l,u in [('isolated',Q(2),Z,SQ/2),
                        ('outer3',Q(1),Z,SQ-1),
                        ('inner3',Q(1),Q(1),SQ-1),
                        ('both3',Z,Q(1),1-SQ/2)]:
        g=2*l*u+c0*u*u
        assert g.sign()>=0 and (1-g).sign()>=0
        rows.append(dict(name=name,c0=c0.wire(),moment=l.wire(),umax=u.wire(),
                         maximum_single_root_cost=g.wire(),upper_bound='1'))
    checks=[endblock(kind,r) for r in range(9,101) for kind in ('B','D')]
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',
           f'results/{BASE}.md',f'scripts/{REVIEW}.py',f'results/{REVIEW}.json',
           f'results/{REVIEW}.md',
           'results/fresh_simple_cusp_vinberg_weight_frontier_20261017.md',
           'results/fresh_simple_cusp_affine_exclusion_20261017.json',
           'results/fresh_simple_cusp_affine_exclusion_20261017.md']
    return dict(status='PASS_ALL_RANK_AFFINE_TREE_CENTRAL_CORRIDOR_PAYMENT',
        sources={p:sha(p) for p in paths},profiles=rows,endblock_checks=checks,
        assumptions=dict(ambient='All-simple finite-volume Coxeter polytope',
            ideal_corner='affine B,C,or D',
            shape='B/D endblock of rank r, simple path a--...--b, B/D endblock of rank s',
            r_min=9,s_min=9,black_distance_min=4,ambient_dimension='r+s+P+1>=23'),
        all_rank_identities=dict(B_inverse_column='1 on the path except last node 1/sqrt(2)',
            D_inverse_column='1 on the main path and1/2 on both fork leaves',
            port_inverse_diagonal='1',schur_alpha='1/P',
            determinant='detS/(h^2*k^2)=1-(2*(l*u+m*v)+c0*u^2+d0*v^2)/P',
            strict_lower_bound='1-2/P>=1/2'),
        theorem='The opposite corner has adjacent marked roots; hence the quadrilateral is paid.',
        scope='All-rank central-corridor branch only. The finite inverse-column checks confirm explicit all-rank identities. Cases with a marked root in a bounded decorated-end region and general quadrilateral payment remain open. No new dimension bound.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],len(result['profiles']),'profiles;',len(result['endblock_checks']),'inverse columns')
