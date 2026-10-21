#!/usr/bin/env python3
"""Six-simplex finite-edge obstruction, all-size white Schur reduction.

Finite exact profile multisets are consequences of minimal order and PSD,
not a cutoff on the supporting white block. Indirect whites are arbitrary.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations,product
import json
from pathlib import Path
import sympy as s
import fresh_finite_edge_cusp_tradeoff12_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_finite_edge_cusp_tradeoff_28over3_20261017'
FRONT='fresh_six_simplex_ideal_frontier_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
PREV='fresh_finite_edge_cusp_tradeoff10_20261017'
Q=A.Q

def sym(M):return s.Matrix([[s.Rational(a)+s.Rational(b)*s.sqrt(2) for a,b in row] for row in M])
def qmat(M):return [[Q(F(a),F(b)) for a,b in row] for row in M]
def data(M):return [[str(s.simplify(M[i,j])) for j in range(M.cols)] for i in range(M.rows)]
def pd(M):return all(s.simplify(M[:k,:k].det())>0 for k in range(1,M.rows+1))

def pair_five_ideal_exclusion(matrices,black):
    rows={b:[-C[1,1],-C[0,0],2*C[0,1],s.Integer(1),s.factor(C.det())] for b,C in matrices.items()}
    assert all(v.is_Rational for row in rows.values() for v in row)
    out=[];x,y,z,t=s.symbols('x y z t')
    for T in combinations(black,5):
        M=s.Matrix([rows[b] for b in T]);rank=M.rank();rankA=M[:,:4].rank()
        rec=dict(five_black=list(T),linear_rows=data(M),coefficient_rank=rankA,augmented_rank=rank)
        if rank>rankA:rec['reason']='inconsistent_linear_determinant_equations'
        else:
            closure=[b for b in black if M.col_join(s.Matrix([rows[b]])).rank()==rank]
            rec['forced_singular_black']=closure
            if len(closure)>=6:rec['reason']='at_least_six_ideals_leaves_at_most_one_finite'
            else:
                sol=list(s.linsolve((M[:,:4],-M[:,4]),(x,y,z,t)))
                assert len(sol)==1 and sol[0]==(7*t/12,7*t/12,7*t/12-s.Rational(3,2),t)
                relation=s.factor(sol[0][0]*sol[0][1]-sol[0][2]**2-sol[0][3])
                roots=s.solve(relation,t);assert roots==[3]
                K=s.Matrix([[sol[0][0],sol[0][2]],[sol[0][2],sol[0][1]]]).subs(t,3)
                assert K==s.Matrix([[s.Rational(7,4),s.Rational(1,4)],[s.Rational(1,4),s.Rational(7,4)]])
                rec.update(reason='determinant_identity_forces_positive_offdiagonal',
                    linear_solution=list(map(str,sol[0])),determinant_relation=str(relation),
                    sole_root='3',forced_K=data(K))
        out.append(rec)
    assert len(out)==21
    return out

def triple_reduction(matrices):
    H={b:2*s.eye(3)-C for b,C in matrices.items()};nulls=[];zero_rows=set();initial=[]
    for b,M in H.items():
        if M.det()!=0:assert pd(M);continue
        initial.append(b);ns=M.nullspace();assert len(ns)==1
        v=ns[0]
        if all(a<=0 for a in v):v=-v
        assert all(a>=0 for a in v) and any(a>0 for a in v)
        support=[i for i,a in enumerate(v) if a>0];zero_rows.update(support)
        nulls.append(dict(black=b,vector=list(map(str,v)),support=support))
    complement=[i for i in range(3) if i not in zero_rows]
    assert len(complement)<=1
    out=dict(right_residuals=[dict(black=b,matrix=data(M)) for b,M in H.items()],
        nonnegative_kernels=nulls,correction_zero_rows=sorted(zero_rows),
        native_contact_matrix_forced_zero=True)
    if not complement:
        assert len(initial)==3
        out.update(indirect_correction_forced_zero=True,maximum_ideals=3)
    else:
        pos=complement[0];thresholds=[]
        for b,M in H.items():
            if M.det()!=0:threshold=1/M.inv()[pos,pos];ground=None
            else:
                ground=next(i for i in sorted(zero_rows) if pd(M.extract([j for j in range(3) if j!=i],[j for j in range(3) if j!=i])))
                keep=[i for i in range(3) if i!=ground];threshold=1/M.extract(keep,keep).inv()[keep.index(pos),keep.index(pos)]
            thresholds.append(dict(black=b,threshold=str(threshold),grounded_zero_row=ground))
        low=min(s.Rational(v['threshold']) for v in thresholds)
        ties=[v['black'] for v in thresholds if s.Rational(v['threshold'])==low]
        assert low==s.Rational(1,2) and len(ties)==1 and len(initial)==1 and not set(ties)&set(initial)
        out.update(indirect_correction_coordinate=pos,thresholds=thresholds,
            minimum='1/2',minimum_ties=ties,existing_ideal_blacks=initial,maximum_ideals=2)
    return out

def sharp_simplex(census):
    L=census['diagrams'][24];assert L['order']==10 and len(L['ideal_deletions'])==3
    G=A.gram(10,L['edges']);assert A.inertia(G)==[9,1,0]
    J=A.inverse(G);assert all(v.sign()<0 for i,row in enumerate(J) for j,v in enumerate(row) if i!=j)
    ideal=[i for i in range(10) if J[i][i].sign()==0]
    assert ideal==L['ideal_deletions'] and all(J[i][i].sign()<=0 for i in range(10))
    center=[-sum(row,Q()) for row in J]
    norm=sum((sum((G[i][j]*center[j] for j in range(10)),Q())*center[i] for i in range(10)),Q())
    assert norm.sign()<0
    cusps=[]
    for b in ideal:
        ids=[i for i in range(10) if i!=b];H=A.sub(G,ids)
        assert A.inertia(H)==[8,0,1]
        kernel=A.kernel1(H)
        cusps.append(dict(ideal_vertex=b,incident_indices=ids,positive_affine_kernel=[v.data() for v in kernel]))
    assert 2*(s.binomial(7,3)+3*s.binomial(7,2))/s.binomial(7,2)==s.Rational(28,3)
    return dict(diagram=24,edges=L['edges'],gram=A.matrix_data(G),inverse=A.matrix_data(J),
        center_coefficients=[v.data() for v in center],center_norm=norm.data(),
        ideal_vertices=ideal,finite_vertices=7,cusp_links=cusps,
        triangle_T0=35,triangle_T1=63,finite_edges=21,twice_triangle_per_finite_edge='28/3',
        geometric_argument='The invertibleLorentz Gram and strict center define a simplex cone generated by the10rays -G^-1 e_i; allarefuturecausal because inverseentries are nonpositive and centerpairingsnegative. All9wall ideal links are connected affine rank8simplex diagrams, hencebounded. Therefore finitevolumeH9 Coxetersimplex; this is a classical census simplex, not a newexample.')

def run():
    front=json.loads((ROOT/f'results/{FRONT}.json').read_text());census=json.loads((ROOT/f'results/{CENSUS}.json').read_text())
    assert front['proved_easy_cases']==465 and front['hard_cases']==74
    hard=[];model_counts=Counter();five_counts=Counter();triples=[]
    for r in front['cases']:
        if r['method']!='no_black_affine_deletion':continue
        Cq={z['black']:qmat(z['matrix']) for z in r['profile_cost_matrices']}
        Cs={z['black']:sym(z['matrix']) for z in r['profile_cost_matrices']}
        p=len(r['profiles']);m=r['mandatory_profile'];bounds=[]
        for i in range(p):
            hi=max(C[i][i].rational() for C in Cq.values());assert hi>0
            bounds.append(int(F(2)/hi))
        trials=[];models=[]
        for counts in product(*(range(b+1) for b in bounds)):
            if not counts[m]:continue
            ids=[i for i,k in enumerate(counts) for _ in range(k)];bad=None
            for b,C in Cq.items():
                H=[[Q(2 if i==j else 0)-C[u][v] for j,v in enumerate(ids)] for i,u in enumerate(ids)]
                if A.inertia(H)[1]:bad=b;break
            trials.append(dict(multiplicities=list(counts),first_indefinite_black=bad))
            if bad is not None:continue
            matrices={b:C.extract(ids,ids) for b,C in Cs.items()}
            rec=dict(profile_indices=ids)
            if len(set(ids))==1:
                costs={b:C[0,0] for b,C in matrices.items()};hi=max(costs.values());ties=[b for b,c in costs.items() if c==hi]
                assert len(ties)<=4 or len(ties)>=6
                rec.update(method='one_profile_rank_one_order',costs=[[b,str(c)] for b,c in costs.items()],
                    maximum=str(hi),maximum_ties=ties,ideal_count_envelope=sorted(set([0,len(ties)])))
            elif len(ids)==2:
                checks=pair_five_ideal_exclusion(matrices,r['black'])
                five_counts.update(z['reason'] for z in checks)
                rec.update(method='two_profile_five_determinant_exclusion',five_vertex_checks=checks)
            else:
                assert len(ids)==3
                reduction=triple_reduction(matrices)
                rec.update(method='three_profile_nonnegative_kernel_reduction',reduction=reduction)
                triples.append([r['diagram'],r['internal_whites'],ids])
            model_counts[rec['method']]+=1;models.append(rec)
        assert models
        hard.append(dict(diagram=r['diagram'],internal_whites=r['internal_whites'],
            per_profile_multiplicity_upper_bounds=bounds,complete_multiset_trials=trials,allowed_models=models,
            conclusion='At least2finite vertices implies at most4ideal vertices.'))
    assert len(hard)==74
    assert triples==[[18,[7],[0,1,2]],[22,[4,8],[0,0,1]],[22,[7,8],[0,1,1]]]
    assert model_counts['two_profile_five_determinant_exclusion']==31
    assert five_counts['determinant_identity_forces_positive_offdiagonal']==2
    coefficient=[]
    for vertices in range(2,11):
        for f in range(2,vertices+1):
            i=vertices-f
            if vertices==7 and i>4:continue
            if vertices in (8,9) and i>3:continue
            if vertices==10 and (f,i) not in [(7,3),(8,2),(9,1)]:continue
            c=2*(F(i)+F(f-2,3));assert c<=F(28,3)
            coefficient.append(dict(dimension=vertices-1,f=f,i=i,coefficient=str(c)))
    assert [r for r in coefficient if r['coefficient']=='28/3']==[dict(dimension=9,f=7,i=3,coefficient='28/3')]
    paths=[f'scripts/{TAG}.py',f'scripts/{FRONT}.py',f'results/{FRONT}.json',f'scripts/{CENSUS}.py',f'results/{CENSUS}.json',
        'scripts/fresh_finite_edge_cusp_tradeoff12_20261017.py',f'scripts/{PREV}.py',f'results/{PREV}.json',f'results/{PREV}.md',
        'results/fresh_finite_ended_simplex_ports_20261017.md','results/fresh_multi_ideal_triangle_scope_20261017.md']
    return dict(status='PASS_EXACT_SIX_SIMPLEX_NO_TWO_FINITE_FIVE_IDEAL_AND_SHARP_FINITE_EDGE_28_OVER_3',
        scope='Proof candidate: actual6simplex with at least2finite vertices has at most4ideals, allowing arbitrary supportingwhites. Sharp uniform finiteedge triangle coefficient28/3. Noambientdimensionbound.',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        easy_case_count=465,hard_case_count=74,all_case_count=539,hard_cases=hard,
        direct_white_model_counts=dict(model_counts),five_determinant_check_counts=dict(five_counts),
        coefficient_cases=coefficient,uniform_finite_edge_coefficient='28/3',sharpness=sharp_simplex(census),
        inequality='2(T0+T1) <= (28/3) E_ff = (14/3)n v_f -(14/3)E_fi',
        residual_budget='2(T2+T3)+Q >= (b_n-(14/3)n)v_f +(14/3)E_fi+alpha_n sum E(L_v)',
        outside_white_argument='WriteF=L-h affine with hwhite for minimum-orderL. DirectU uses the complete catalogueprofiles. IndirectZ is rightF and arbitrary. K=G_U-R^tG_Z^-1R is positive definite, has nonpositiveoffdiagonal, and diagonal<=2. Actual blackvertex matrices aftercongruence areK-C_b. Entrywise monotonicity ofPSD Zmatrices allows the necessary relaxation2I-C_b>=0. Repetition bounds and all multisets follow withoutwhite-rankcutoff.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(out.read_text())==d
    else:out.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],d['direct_white_model_counts'],d['five_determinant_check_counts'],'saved_replay=',a.verify,flush=True)
