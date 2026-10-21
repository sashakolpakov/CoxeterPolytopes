#!/usr/bin/env python3
"""Exact seven-simplex partition and nine finite-black hard-case reductions.

Uses the frozen standalone Fraction Qsqrt2 matrix engine; no producer census
arithmetic imports. Four other hard cases are explicitly outside this package.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_edge_cusp_tradeoff12_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_finite_edge_seven_simplex_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
Q=A.Q

def partition(diagrams):
    rows=[]
    for L in diagrams:
        n=L['order']
        if n<8:continue
        G=A.gram(n,L['edges']);assert A.inertia(G)==[n-1,1,0]
        invs={};matrices={};aff=[]
        for j in range(n):
            ids=[k for k in range(n) if k!=j];H=A.sub(G,ids);I=A.inertia(H)
            assert I in ([n-1,0,0],[n-2,0,1]);matrices[j]=(ids,H)
            if I[-1]:aff.append(j);A.kernel1(H)
            else:invs[j]=A.inverse(H)
        assert aff==L['ideal_deletions']
        for W in combinations(range(n),n-8):
            black=[j for j in range(n) if j not in W];I=[j for j in aff if j in black]
            row=dict(diagram=L['index'],order=n,internal_whites=list(W),black=black,
                     black_affine_deletions=I)
            if len(I)>=2:
                row.update(method='external_whites_orthogonal_to_L',possible_ideal_counts=[len(I)])
            elif len(I)==1:
                a=I[0];ts=[]
                for j in black:
                    if j==a:continue
                    ids,H=matrices[j];t=(1/invs[j][ids.index(a)][ids.index(a)]).rational()
                    assert t>0;ts.append([j,str(t)])
                minimum=min(F(t) for _,t in ts);ties=[j for j,t in ts if F(t)==minimum]
                for j in black:
                    if j==a:continue
                    ids,H=matrices[j];K=[r[:] for r in H]
                    K[ids.index(a)][ids.index(a)]-=Q(minimum)
                    assert A.inertia(K)==([n-2,0,1] if j in ties else [n-1,0,0])
                row.update(method='rank_one_diagonal_decrease',decreased_black=a,
                           thresholds=ts,minimum=str(minimum),ties=ties,
                           possible_ideal_counts=sorted(set([1,1+len(ties)])))
            else:
                row.update(method='no_black_affine_deletion',
                           scope='finite_black_E8_reduction' if L['index']==25 else 'separate_root_lane')
            if I:
                allowed=[i for i in row['possible_ideal_counts'] if 8-i>=2]
                assert all(i<=3 for i in allowed)
                row['ideal_counts_with_at_least_two_finite_vertices']=allowed
            rows.append(row)
    counts=Counter(r['method'] for r in rows)
    assert len(rows)==175
    assert counts==Counter(external_whites_orthogonal_to_L=77,
                           rank_one_diagonal_decrease=85,no_black_affine_deletion=13)
    hard=[(r['diagram'],r['internal_whites']) for r in rows if not r['black_affine_deletions']]
    assert hard==[(20,[8]),(21,[8]),(22,[8]),(23,[8,9])]+[(25,[a,9]) for a in range(9)]
    allideal=[r for r in rows if 8 in r.get('possible_ideal_counts',[])]
    assert len(allideal)==1 and allideal[0]['diagram']==19 and allideal[0]['internal_whites']==[8]
    return rows,dict(counts)

def finite_black_reduction():
    edges=[[0,1,3],[0,2,3],[0,4,3],[2,3,3],[4,5,3],[5,6,3],[6,7,3],[7,8,3]]
    G=A.gram(9,edges);assert A.inertia(G)==[8,0,1];A.kernel1(G)
    inv={};cost=[];inverse_records=[]
    for j in range(9):
        ids=[i for i in range(9) if i!=j];H=A.sub(G,ids)
        assert A.inertia(H)==[8,0,0];J=A.inverse(H);inv[j]=(ids,J)
        assert all(x.sign()>=0 for r in J for x in r)
        c=F(0) if j==8 else J[ids.index(8)][ids.index(8)].rational();cost.append(c)
        inverse_records.append(dict(deleted=j,indices=ids,inverse=A.matrix_data(J)))
    assert cost==list(map(F,['5/6','8/9','7/8','1','4/5','3/4','2/3','1/2','0']))
    ids8,J8=inv[8]
    assert [J8[i][i].rational() for i in range(8)]==list(map(F,[30,8,14,4,20,12,6,2]))
    ids3,J3=inv[3]
    assert [[J3[ids3.index(i)][ids3.index(j)] for j in (7,8)] for i in (7,8)]==[[Q(2),Q(1)],[Q(1),Q(1)]]
    profiles=[('A',8,Q(1)),('B',8,Q(0,1)),('C',7,Q(1))]
    pair_checks=[]
    for vertex in (1,3):
        ids,J=inv[vertex];C=[]
        for name,node,m in profiles:C.append([m if i==node else Q() for i in ids])
        costs=[[sum((C[a][i]*J[i][j]*C[b][j] for i in range(8) for j in range(8)),Q())
                for b in range(3)] for a in range(3)]
        dets=[]
        for k in (1,2):
            det=(Q(2)-costs[0][0])*(Q(2)-costs[k][k])-costs[0][k]*costs[0][k]
            assert det.sign()<0
            dets.append(dict(other_profile=profiles[k][0],right_native_contact_residual_determinant=det.data()))
        expected=([Q(F(-4,3)).data(),Q(F(-1,9)).data()] if vertex==1
                  else [Q(-2).data(),Q(-1).data()])
        assert [x['right_native_contact_residual_determinant'] for x in dets]==expected
        pair_checks.append(dict(actual_black_vertex=vertex,cost_gram=A.matrix_data(costs),
                               incompatible_with_existing_A=dets))
    cases=[]
    for a in range(9):
        black=[i for i in range(9) if i!=a];maximum=max(cost[j] for j in black)
        where=[j for j in black if cost[j]==maximum];assert where==([1] if a==3 else [3])
        cases.append(dict(diagram=25,internal_whites=[a,9],finite_black=black,
            existing_affine_E8_internal_white=a,existing_A_white=9,
            A_B_C_pair_test_vertex=1 if a==3 else 3,
            black_vertex_costs=[[j,str(cost[j])] for j in black],
            unique_maximum_cost=str(maximum),unique_possible_ideal_vertex=where[0],
            all_individual_black_plus_white_PSD_branch_maximum_ideals=1,
            indefinite_individual_extension_branch='Transfers to an order9 minimal diagram with a black affine deletion, already in the162easycases; finiteB excludes its allideal exception.',
            overall_ideal_upper_bound_with_at_least_two_finite_vertices=3))
    return dict(affine_E8_gram=A.matrix_data(G),black_deletion_inverses=inverse_records,
                direct_profile_pair_exclusions=pair_checks,
                all_A_projection_costs=list(map(str,cost)),cases=cases,
                unrestricted_white_tail_argument='After eliminating whites Z orthogonal to the entire affineE8 block, K=G_U-R^tG_Z^-1R is positive definite because U+Z is a subset of the spherical whiteblock. For blackdeletionj the residual isK-c_j11^t. PSD at the unique largest c_j forcesPD at every smaller c_j, independently of the cardinality or graph ofZ andU.')

def run():
    d=json.loads((ROOT/f'results/{CENSUS}.json').read_text())
    assert d['status']=='PASS_EXACT_HIGH_RANK_QUASILANNER_ATTACHMENT_CENSUS'
    cases,counts=partition(d['diagrams']);hard=finite_black_reduction()
    paths=[f'scripts/{TAG}.py',f'scripts/{CENSUS}.py',f'results/{CENSUS}.json',
           'scripts/fresh_finite_edge_cusp_tradeoff12_20261017.py',
           'results/fresh_finite_edge_cusp_tradeoff12_20261017.json',
           'results/fresh_finite_edge_cusp_tradeoff12_20261017.md',
           'results/fresh_finite_ended_simplex_ports_20261017.md',
           'results/fresh_multi_ideal_triangle_scope_20261017.md']
    return dict(status='PASS_EXACT_SEVEN_SIMPLEX_162_EASY_AND_NINE_FINITE_BLACK_REDUCTIONS',
        scope='All175minimaldiagram choices classified;162easycases and9finiteblackhardcases controlled forf>=2. Fourhardcases remain a separate dependency; no fullseven-simplextheorem or newdimensionclaim here.',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=counts,total_cases=len(cases),cases=cases,finite_black_reduction=hard,
        outstanding_cases=[[20,[8]],[21,[8]],[22,[8]],[23,[8,9]]],
        proved_cases=171,maximum_ideals_in_proved_cases_with_at_least_two_finite_vertices=3,
        all_ideal_boundary_exception=dict(diagram=19,internal_whites=[8],
                                         possible_ideal_counts=[1,8],excluded_from_finite_edge_blocks=True))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==data
    else:p.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print(data['status'],'saved_replay=',args.verify)
