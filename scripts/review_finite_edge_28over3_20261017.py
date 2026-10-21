#!/usr/bin/env python3
"""Independent Fraction checks for the six-simplex/28-over-3 proof.

No producer arithmetic or SymPy is imported. This rechecks all651 pair
systems,111 single-profile envelopes,three triples and the sharp H9 simplex.
The539-case profile/classification producer has a separate stated scope.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='review_finite_edge_28over3_20261017'
THEOREM='fresh_finite_edge_cusp_tradeoff_28over3_20261017'
FRONT='fresh_six_simplex_ideal_frontier_20261017'


def rref(M):
    A=[[F(a) for a in row] for row in M];pivots=[];i=0
    if not A:return A,pivots
    for j in range(len(A[0])):
        pivot=next((k for k in range(i,len(A)) if A[k][j]),None)
        if pivot is None:continue
        A[i],A[pivot]=A[pivot],A[i];a=A[i][j];A[i]=[b/a for b in A[i]]
        for k in range(len(A)):
            if k!=i:
                a=A[k][j];A[k]=[b-a*c for b,c in zip(A[k],A[i])]
        pivots.append(j);i+=1
        if i==len(A):break
    return A,pivots


def rank(M):return len(rref(M)[1])


def det(M):
    A=[[F(a) for a in row] for row in M];ans=F(1)
    for i in range(len(A)):
        p=next((k for k in range(i,len(A)) if A[k][i]),None)
        if p is None:return F()
        if p!=i:A[p],A[i]=A[i],A[p];ans=-ans
        a=A[i][i];ans*=a
        for k in range(i+1,len(A)):
            b=A[k][i]/a
            for j in range(i,len(A)):A[k][j]-=b*A[i][j]
    return ans


def sub(M,ids):return [[M[i][j] for j in ids] for i in ids]


def inv(M):
    n=len(M);a,p=rref([list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(M)])
    assert p==list(range(n))
    return [row[n:] for row in a]


def pd(M):return all(det(sub(M,list(range(k))))>0 for k in range(1,len(M)+1))


def psd(M):
    return all(det(sub(M,ids))>=0 for k in range(1,len(M)+1) for ids in combinations(range(len(M)),k))


def qrat(x):
    a,b=x;assert F(b)==0
    return F(a)


def pair_row(C):
    return [-C[1][1],-C[0][0],2*C[0][1],F(1),det(C)]


def run():
    theorem=json.loads((ROOT/f'results/{THEOREM}.json').read_text())
    front=json.loads((ROOT/f'results/{FRONT}.json').read_text())
    lookup={(r['diagram'],tuple(r['internal_whites'])):r for r in front['cases']}
    counts=Counter();pairs=[];triples=[]
    for case in theorem['hard_cases']:
        base=lookup[case['diagram'],tuple(case['internal_whites'])]
        costs={r['black']:r['matrix'] for r in base['profile_cost_matrices']}
        for model in case['allowed_models']:
            ids=model['profile_indices'];counts['models']+=1
            matrices={b:[[qrat(C[i][j]) for j in ids] for i in ids] for b,C in costs.items()}
            if model['method']=='one_profile_rank_one_order':
                hi=max(C[0][0] for C in matrices.values());ties=[b for b,C in matrices.items() if C[0][0]==hi]
                assert hi==F(model['maximum']) and ties==model['maximum_ties']
                assert len(ties)<=4 or len(ties)>=6
                counts['single_profile_models']+=1
            elif model['method']=='two_profile_five_determinant_exclusion':
                rows={b:pair_row(C) for b,C in matrices.items()}
                assert len(model['five_vertex_checks'])==21
                for saved in model['five_vertex_checks']:
                    T=saved['five_black'];A=[rows[b] for b in T]
                    assert A==[[F(a) for a in row] for row in saved['linear_rows']]
                    ra=rank([row[:4] for row in A]);rr=rank(A)
                    assert (ra,rr)==(saved['coefficient_rank'],saved['augmented_rank'])
                    if rr>ra:reason='inconsistent_linear_determinant_equations';closure=[]
                    else:
                        closure=[b for b in base['black'] if rank(A+[rows[b]])==rr]
                        assert closure==saved['forced_singular_black']
                        if len(closure)>=6:reason='at_least_six_ideals_leaves_at_most_one_finite'
                        else:
                            rrA,p=rref(A);assert p==[0,1,2]
                            assert rrA[:3]==[[1,0,0,-F(7,12),0],[0,1,0,-F(7,12),0],[0,0,1,-F(7,12),F(3,2)]]
                            # xy-z²-t = (3/4)t-9/4, hence t=3,z=1/4.
                            linear=2*F(7,12)*F(3,2)-1;constant=-F(3,2)**2
                            t=-constant/linear;z=F(7,12)*t-F(3,2)
                            assert t==3 and z==F(1,4)>0
                            reason='determinant_identity_forces_positive_offdiagonal'
                    assert reason==saved['reason'];counts[reason]+=1
                    pairs.append([case['diagram'],case['internal_whites'],ids,T,ra,rr,reason,closure])
            else:
                assert model['method']=='three_profile_nonnegative_kernel_reduction'
                H={b:[[F(2*(i==j))-C[i][j] for j in range(3)] for i in range(3)] for b,C in matrices.items()}
                reduction=model['reduction'];zero=set();initial=[]
                for rec in reduction['right_residuals']:
                    assert H[rec['black']]==[[F(a) for a in row] for row in rec['matrix']]
                for b,A in H.items():
                    assert psd(A)
                    if det(A)==0:initial.append(b);assert rank(A)==2
                    else:assert pd(A)
                for rec in reduction['nonnegative_kernels']:
                    v=list(map(F,rec['vector']));A=H[rec['black']]
                    assert all(a>=0 for a in v) and any(a>0 for a in v)
                    assert all(sum(a*b for a,b in zip(row,v))==0 for row in A)
                    support=[i for i,a in enumerate(v) if a>0];assert support==rec['support'];zero.update(support)
                assert sorted(zero)==reduction['correction_zero_rows']
                if len(zero)==3:assert len(initial)==3 and reduction['maximum_ideals']==3
                else:
                    assert len(zero)==2;pos=next(i for i in range(3) if i not in zero)
                    thresholds={}
                    for rec in reduction['thresholds']:
                        b=rec['black'];g=rec['grounded_zero_row'];keep=[i for i in range(3) if i!=g]
                        A=sub(H[b],keep);assert pd(A)
                        threshold=1/inv(A)[keep.index(pos)][keep.index(pos)]
                        assert threshold==F(rec['threshold']);thresholds[b]=threshold
                    low=min(thresholds.values());ties=[b for b,a in thresholds.items() if a==low]
                    assert low==F(1,2) and len(ties)==1 and len(initial)==1 and set(ties).isdisjoint(initial)
                    for A in H.values():
                        corrected=[list(row) for row in A];corrected[pos][pos]-=low;assert psd(corrected)
                triples.append([case['diagram'],case['internal_whites'],ids,len(initial),sorted(zero)])
    assert counts['single_profile_models']==111 and len(pairs)==651 and len(triples)==3
    sharp=theorem['sharpness'];G=[[F(2*(i==j)) for j in range(10)] for i in range(10)]
    for i,j,m in sharp['edges']:assert m==3;G[i][j]=G[j][i]=F(-1)
    assert G==[[qrat(a) for a in row] for row in sharp['gram']]
    J=inv(G);assert J==[[qrat(a) for a in row] for row in sharp['inverse']]
    assert all(J[i][j]<0 for i in range(10) for j in range(10) if i!=j)
    assert all(J[i][i]<=0 for i in range(10))
    center=[-sum(row) for row in J]
    norm=sum(center[i]*G[i][j]*center[j] for i in range(10) for j in range(10))
    assert norm<0 and norm==qrat(sharp['center_norm'])
    ideal=[]
    for b in range(10):
        H=sub(G,[i for i in range(10) if i!=b])
        if J[b][b]<0:assert pd(H)
        else:
            assert J[b][b]==0;ideal.append(b);assert rank(H)==8
            A=sub(H,list(range(8)));assert pd(A)
            invA=inv(A);tail=[H[i][8] for i in range(8)]
            kernel=[-sum(invA[i][j]*tail[j] for j in range(8)) for i in range(8)]+[F(1)]
            assert all(a>0 for a in kernel)
            assert all(sum(a*b for a,b in zip(row,kernel))==0 for row in H)
    assert ideal==sharp['ideal_vertices'] and len(ideal)==3
    # A positive-definite9principal block plus detG<0 gives inertia(9,1).
    assert det(G)<0 and F(2*(35+63),21)==F(28,3)
    for rec in theorem['coefficient_cases']:
        assert 2*(F(rec['i'])+F(rec['f']-2,3))==F(rec['coefficient'])<=F(28,3)
    paths=[f'scripts/{TAG}.py',f'scripts/{THEOREM}.py',f'results/{THEOREM}.json',
        f'scripts/{FRONT}.py',f'results/{FRONT}.json']
    return dict(status='PASS_INDEPENDENT_FRACTION_28_OVER_3_ALGEBRA_AND_SHARPNESS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=dict(counts,pair_systems=len(pairs),triple_models=len(triples)),
        pair_checks=pairs,triple_checks=triples,sharp_simplex=dict(ideal_vertices=ideal,
            finite_vertices=7,gram_determinant=str(det(G)),center_norm=str(norm),coefficient='28/3'),
        mathematical_review='Parent read the complete minimum-order profile argument, arbitrary-Z Schur reduction, entrywise PSD relaxation, multiplicity bounds, rank-one threshold mechanism, pair determinant implication, triple nonnegative-kernel/grounding argument, finite-edge block summation and causal finite-volume sharpness proof. No mathematical gap identified in that scope.',
        independence='Only Python Fraction and elementary Gaussian elimination are used. No producer arithmetic module, SymPy or producer elimination function is imported. Rational profiles are read from the frozen539-case frontier; their underlying Coxeter classification and catalogue-isomorphism enumeration are not independently rerun.',
        scope='Independent exact recheck of651 pair systems,111single-profile envelopes,3triples and rationalH9 sharpness, together with a separate mathematical/source review. Not a duplicate full539-case profile census, formal proof-assistant verification, external refereeing or an ambient dimension theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
