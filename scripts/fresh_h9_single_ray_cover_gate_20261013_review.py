#!/usr/bin/env python3
"""Separate sparse-radical replay of support bounds and covering margins."""
import argparse,json,sys
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_makarov_frontier_20261013 as E
A=E.A;R=E.R
TAG='fresh_h9_single_ray_cover_gate_20261013'
SOURCE='results/'+TAG+'.json';OUT='results/'+TAG+'_review.json'
def run(verify=False):
    d=R.read(SOURCE);s=R.read(A.SOURCE)
    assert d['status']=='PASS_EXACT_FOUR_SEPARATE_SINGLE_RAY_NO_COMPLETION_GATES'
    assert d['script_sha256']==R.sha('scripts/'+TAG+'.py') and d['source_sha256']==R.sha(A.SOURCE)
    G=[[A.k(z)for z in row[:6]]for row in s['form'][:6]]
    W=[[A.k(z)for z in row[:6]]for row in s['seed_normals'][:7]]
    V=[[A.k(z)for z in row]for row in s['compact_H5_vertices']]
    B=[[A.k(z)for z in row]for row in s['H4_inverse']]
    rows=[A.mv(G,w)for w in W];low=[R.C(1)]*6+[R.C(F(2,5))]
    assert R.sign(A.sub(R.scale(A.pair(W[6],G,W[6]),F(1,2)),R.mul(low[6],low[6])))>0
    inverses={i:E.matrix_inverse([r for j,r in enumerate(rows)if j!=i])for i in range(7)}
    assert not d['singular_six_constraint_omissions']
    target=A.sub(B[3][3],R.C(2));assert target=={1:F(1),5:F(1)}
    minima={};empty=0;numvertices=0;admissible=[]
    assert len(d['support_vertex_census'])==127
    for bits,saved in enumerate(d['support_vertex_census'],1):
        nz=[i for i in range(7)if bits>>i&1];zz=[i for i in range(7)if i not in nz]
        assert saved['nonzero_indices']==nz and saved['zero_indices']==zz
        rhs=[R.scale(low[i],-1)if i in nz else {}for i in range(7)];candidates=[]
        for omitted in nz:
            ids=[i for i in range(7)if i!=omitted];u=A.mv(inverses[omitted],[rhs[i]for i in ids]);values=A.mv(rows,u)
            if any(R.sign(A.sub(a,b))>0 for a,b in zip(values,rhs,strict=True)):continue
            assert all(not values[i]for i in zz)
            q=R.scale(A.pair(u,G,u),-1);assert R.sign(q)>0
            candidates.append((omitted,u,q));numvertices+=1
            n=len(zz)
            if n not in minima or R.sign(A.sub(q,minima[n]))<0:minima[n]=q
            if n<=4:assert R.sign(A.sub(q,target))>0
        assert len(candidates)==len(saved['relaxed_polyhedron_vertices'])
        for (omitted,u,q),v in zip(candidates,saved['relaxed_polyhedron_vertices'],strict=True):
            assert omitted==v['omitted'] and u==list(map(A.k,v['normal_old_part'])) and q==A.k(v['positive_lorentz_norm_square'])
        if not candidates:empty+=1
        elif any(R.sign(A.sub(q,target))<=0 for _,_,q in candidates):
            assert len(zz)==5 and zz in s['vertex_active_sets'];admissible.append(zz)
    assert empty==d['empty_relaxed_polyhedra']==34
    assert admissible==d['j3_m3_possible_zero_sets_after_relaxation']
    assert {n:A.k(z)for n,z in d['minimum_lorentz_norm_square_by_zero_count'].items()}=={str(n):z for n,z in minima.items()}

    x=list(map(A.k,d['strict_old_interior_gap_point']));xx=A.pair(x,G,x)
    assert R.sign(xx)<0 and all(R.sign(z)<0 for z in A.mv(rows,x))
    C3=E.div(R.mul(B[3][0],B[3][0]),R.mul(B[0][0],target));assert C3==A.k(d['j3_m3_ball_cosh_radius_square'])
    gaps=[]
    for v,saved in zip(V,d['ten_strict_vertex_ball_gaps'],strict=True):
        px=A.pair(v,G,x);assert R.sign(px)<0
        gap=A.sub(R.mul(px,px),R.mul(C3,R.mul(A.pair(v,G,v),xx)))
        assert gap==A.k(saved) and R.sign(gap)>0;gaps.append(R.enc(gap))
    i,j=d['diameter_vertex_indices'];vv=A.pair(V[i],G,V[i]);ww=A.pair(V[j],G,V[j]);vw=A.pair(V[i],G,V[j])
    diameter=E.div(R.mul(vw,vw),R.mul(vv,ww));assert diameter==A.k(d['diameter_lower_cosh_square'])
    records=[]
    for row in d['four_separate_ray_gates']:
        j=row['ray'];i=row['worst_axis'];P=E.div(R.mul(B[j][i],B[j][i]),B[i][i]);D=A.sub(B[j][j],P)
        assert D==A.k(row['D'])
        for k in range(4):assert R.sign(A.sub(R.mul(B[i][i],R.mul(B[j][k],B[j][k])),R.mul(B[k][k],R.mul(B[j][i],B[j][i]))))>=0
        if j<2:assert R.sign(A.sub(D,R.C(2)))>0
        else:
            t2=row['minimum_isolated_t_square'];margin=A.sub(R.scale(D,t2+1),R.C(2))
            assert margin==A.k(row['strict_ball_separation_margin']) and R.sign(margin)>0
            C=E.div(P,A.sub(B[j][j],R.C(F(2,t2))));assert C==A.k(row['largest_isolated_ball_cosh_radius_square']) and R.sign(A.sub(C,R.C(1)))>0
            twice_radius=R.mul(A.sub(R.scale(C,2),R.C(1)),A.sub(R.scale(C,2),R.C(1)))
            gap=A.sub(diameter,twice_radius);assert gap==A.k(row['diameter_margin'])and R.sign(gap)>0
            ineffective=2 if j==2 else 3;assert R.sign(A.sub(D,R.C(F(2,ineffective))))>0
        records.append({'ray':j,'axis':i,'D':R.enc(D)})
    out={'status':'PASS_INDEPENDENT_SINGLE_RAY_SUPPORT_AND_BALL_MARGIN_REVIEW',
         'sources':{p:R.sha(p)for p in [SOURCE,'scripts/'+TAG+'.py',A.SOURCE,'scripts/'+Path(A.__file__).name,'scripts/'+Path(E.__file__).name,'scripts/fresh_h22_multiquadratic_20261013_review.py','scripts/'+Path(__file__).name]},
         'support_count':127,'empty_relaxed_polyhedra':empty,'enumerated_polyhedron_vertex_records':numvertices,
         'minimum_square_by_zero_count':{str(n):R.enc(q)for n,q in sorted(minima.items())},
         'possible_zero_sets':admissible,'ten_positive_gap_values':gaps,'four_ray_margins':records,
         'scope':'Independent squarefree-radical arithmetic replays the entire 127-support relaxation and all saved covering margins. Separate mathematical review supplies the Minkowski-Weyl/reverse-triangle reduction, finite separated-ball covering argument, and strict ideal-patch infinite-volume conclusion. Four separate fixed-ray families only; no mixed-ray or retained-twelfth-cap obstruction is inferred.'}
    if verify:assert out==R.read(OUT)
    else:(R.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
