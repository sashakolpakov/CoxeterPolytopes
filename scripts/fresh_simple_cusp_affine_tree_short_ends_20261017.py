#!/usr/bin/env python3
"""All-rank affine-tree payment with short outer end blocks; no dimension claim."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import fresh_simple_cusp_affine_a_quadrilateral_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_simple_cusp_affine_tree_short_ends_20261017'
Q,Z,SQ=A.Q,A.Z,A.SQ

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def end(kind,r):
    assert (kind=='B' and r>=1) or (kind=='D' and r>=3)
    G=[[Q(2) if i==j else Z for j in range(r)] for i in range(r)]
    if kind=='B':
        edges=[(i,i+1,SQ if i==r-2 else Q(1)) for i in range(r-1)]
        if r==1:J=[[Q(F(1,2))]]
        else:
            J=[[Q(min(i,j)+1) if max(i,j)<r-1 else
                Q(F(r,2)) if i==j==r-1 else SQ*F(min(i,j)+1,2)
                for j in range(r)] for i in range(r)]
        v=[SQ if r==1 else Q(1)]+[Z]*(r-1)
    else:
        edges=[(i,i+1,Q(1)) for i in range(r-3)]+[(r-3,r-2,Q(1)),(r-3,r-1,Q(1))]
        J=[[Q(min(i,j)+1) if max(i,j)<r-2 else
            Q(F(r if i==j else r-2,4)) if min(i,j)>=r-2 else
            Q(F(min(i,j)+1,2)) for j in range(r)] for i in range(r)]
        v=[Q(1)]+[Z]*(r-1)
    for i,j,w in edges:G[i][j]=G[j][i]=-w
    for i in range(r):
        # Sparse multiplication checks the displayed closed inverse formula.
        for j in range(r):
            assert sum((G[i][k]*J[k][j] for k in range(r) if G[i][k]!=0),Z)==Q(i==j)
    t=[sum((J[i][j]*v[j] for j in range(r)),Z) for i in range(r)]
    assert sum((v[i]*t[i] for i in range(r)),Z)==1
    assert all((x-Q(F(1,2))).sign()>=0 for x in t)
    low=[i for i in range(r) if (J[i][i]+t[i]**2-2).sign()<0]
    predicted=([0] if r==1 else [1] if r==2 else []) if kind=='B' else ([r-2,r-1] if r<=6 else [])
    assert low==predicted
    assert all((2*(J[i][i]+t[i]**2)-2).sign()>=0 for i in range(r))
    return dict(type=kind,rank=r,low_simple_contact_nodes=low,
                diagonal=[J[i][i].wire() for i in range(r)],
                cusp_contact_coefficients=[x.wire() for x in t])

def run():
    previous='fresh_simple_cusp_affine_tree_corridor_20261017'
    old=json.loads((ROOT/f'results/{previous}.json').read_text())
    assert old['status']=='PASS_ALL_RANK_AFFINE_TREE_CENTRAL_CORRIDOR_PAYMENT'
    for p,h in old['sources'].items():assert sha(p)==h,p
    checks=[end('B',r) for r in range(1,41)]+[end('D',r) for r in range(3,41)]
    joins=[]
    for name,q,t in [('B1_finite_end3',Q(F(1,2)),SQ/2),
                     ('ordinary_affine_port',Q(1),Q(1)),
                     ('B2_swapped_port',Q(1),SQ/2),
                     ('D3_finite_leaf',Q(F(3,4)),Q(F(1,2))),
                     ('D4_other_leaf',Q(1),Q(F(1,2)))]:
        u=Q(1)/(1+t);value=2*u+(1-q)*u*u
        assert (1-q).sign()>=0 and (Q(F(13,9))-value).sign()>=0
        joins.append(dict(name=name,outer_cost=q.wire(),outer_cusp_contact=t.wire(),
                          maximum_reciprocal_height=u.wire(),contribution=value.wire()))
    assert joins[3]['contribution']==Q(F(13,9)).wire()
    # q>=t^2 and t>=1/2 imply (2-q)/(1+t)^2<=1.
    assert 2*F(1,2)**2+2*F(1,2)-1==F(1,2)>0
    assert 2*(SQ-1)+(SQ-1)**2==1
    assert 2*(SQ/2)**2==1
    f4=A.affine_witness('affine_F4_low_outer_leaf',5,
        [(0,1,Q(1)),(1,2,SQ),(2,3,Q(1)),(3,4,Q(1))],
        [SQ,2*SQ,Q(3),Q(2),Q(1)])
    edges=[];kernel=[Q(6)];nextnode=1
    for length in (1,2,5):
        last=0
        for j in range(1,length+1):
            edges.append((last,nextnode,Q(1)));last=nextnode;nextnode+=1
            kernel.append(Q(F(6*(length+1-j),length+1)))
    e8=A.affine_witness('affine_E8_low_D_fork_leaf',9,edges,kernel)
    remaining=[]
    for P in range(4,8):
        for kind,lo in [('B',1),('D',3)]:
            for r in range(lo,9-P):remaining.append(dict(short_end_type=kind,short_end_rank=r,black_distance=P,inner_rank=r+P))
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',
           f'results/{A.TAG}.md',f'scripts/{previous}.py',f'results/{previous}.json',
           f'results/{previous}.md',
           'results/fresh_simple_cusp_vinberg_weight_frontier_20261017.md',
           'results/fresh_simple_cusp_affine_exclusion_20261017.json',
           'results/fresh_simple_cusp_affine_exclusion_20261017.md']
    return dict(status='PASS_ALL_RANK_AFFINE_TREE_SHORT_OUTER_END_PAYMENT',
        sources={p:sha(p) for p in paths},end_formula_checks=checks,
        joining_profile_costs=joins,small_affine_witnesses=[f4,e8],
        assumptions=dict(ambient='all-simple finite-volume Coxeter polytope of dimension n>=14',
            ideal_corner='affine B,C,D tree',shape='F_r -- a -- central path -- b -- F_s',
            allowed_endblocks='B_r,r>=1 (B1=A1 with attaching bond4), or D_r,r>=3 (D3=A3 attached at its middle)',
            central_black_distance='P>=4',inner_rank_conditions=['r+P>=9','s+P>=9'],
            node_count='n=r+s+P+1'),
        complete_profile_classes=['isolated','outer-only, arbitrary finite extension','inner-only long-port3','join both components, full classical finite/affine deletion list'],
        all_rank_inverse='Closed B/D inverse formulas, not a rank cutoff; checked in78 examples.',
        cost_bounds=dict(isolated='1',outer_only='1',inner_only='1',join='13/9',
            determinant_ratio='1-(g_c+g_d)/P >= 1-26/(9P) >= 5/18'),
        remaining_central_positions_up_to_left_right_symmetry=remaining,
        theorem='The opposite marked roots are adjacent; therefore every quadrilateral in this branch is paid.',
        scope='All-rank short-outer-end reduction. Remaining inner-rank<=8 central positions and fork-leaf/center positions stay open; no dimension bound or full quadrilateral theorem.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],len(result['end_formula_checks']),'inverse formulas;',
          len(result['joining_profile_costs']),'join costs;',
          len(result['remaining_central_positions_up_to_left_right_symmetry']),'remaining central positions')
