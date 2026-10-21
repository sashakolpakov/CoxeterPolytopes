#!/usr/bin/env python3
"""Explicit integral Q(sqrt5) compact local triangle missed by the saved generic class rule."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as A
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_discriminant_gate_20261017 as G
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_a13_triangle_failure_20261017'
LP='fresh_sqrt5_rank13_discriminant_simplex_lp_20261017'
SEP='fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017'
K,Z,O=A.K,A.Z,A.O

def run():
    phi=K(F(1,2),F(1,2));M=[[K(2*(i==j)) for j in range(14)] for i in range(14)]
    edges=[[i,i+1,3] for i in range(12)]+[[1,13,5]]
    graph={i:{} for i in range(14)}
    for i,j,m in edges:graph[i][j]=graph[j][i]=m;M[i][j]=M[j][i]=-phi if m==5 else -O
    assert A.inertia(M)==[13,1,0] and A.inertia([[x.conj() for x in row] for row in M])==[14,0,0]
    assert all(x.integral() for row in M for x in row)
    old=A.sub(M,list(range(13)));inverse_old=A.inverse(old)
    assert inverse_old[1][1]==K(F(12,7))
    schur=K(2)-phi*phi*inverse_old[1][1]
    determinant=14*schur;assert determinant==K(-8,-12)
    black=[1,2,13];white=sorted(set(range(14))-set(black));inverse=A.inverse(M);rays=[];deletions=[]
    finite={B.C.code(g):name for name,g in B.C.models() if len(g)<=13}
    lp=B.read(LP);gate=B.read(G.TAG)
    assert all(G.ksqrt(determinant/K(*c['representative'])) is None for c in gate['classes'])
    case=next(r for r in lp['results'] if r['discriminant_class'] is None)
    primal=next(c for c in case['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL');weights={i:F(a) for i,a in primal['weights']}
    lookup={(tuple(v['forest']),B.dump(v['marked_forest_code'])):v['index'] for v in lp['variables']};corners=[]
    for removed in range(14):
        ids=[i for i in range(14) if i!=removed];inertia=A.inertia(A.sub(M,ids));deletions.append(dict(removed=removed,physical_inertia=inertia))
        if removed not in black:assert inertia!=[13,0,0];continue
        assert inertia==[13,0,0]
        end={i:{j:m for j,m in row.items() if j!=removed} for i,row in graph.items() if i!=removed}
        pieces=B.C.pieces(end);forest=tuple(sorted(finite[B.C.code(p)] for p in pieces));marks=[i for i in black if i!=removed]
        code=B.dump(sorted(B.C.code(p,marks) for p in pieces));index=lookup[forest,code]
        ray=[-inverse[i][removed] for i in range(14)]
        assert A.mv(M,ray)==[-O if i==removed else Z for i in range(14)]
        assert A.pair(M,ray,ray).sign()<0;rays.append(ray)
        corners.append(dict(removed=removed,forest=list(forest),marked_pair=marks,variable=index,payment=str(weights.get(index,F(0))),
            ray=[x.enc() for x in ray],norm=A.pair(M,ray,ray).enc()))
    assert sorted(c['variable'] for c in corners)==[3069,12403,12481]
    assert all(A.pair(M,v,w).sign()<0 for v in rays for w in rays)
    h=[sum((r[i] for r in rays),Z) for i in range(14)]
    push=A.mv(inverse,[-O if i in white else Z for i in range(14)]);epsilon=F(1)
    while True:
        center=[h[i]+epsilon*push[i] for i in range(14)]
        if A.pair(M,center,center).sign()<0 and all(A.pair(M,center,r).sign()<0 for r in rays):break
        epsilon/=2
    assert all(v.sign()<0 for v in A.mv(M,center))
    total=sum((F(c['payment']) for c in corners),F(0));assert total==0
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{G.TAG}.py',f'results/{LP}.json',f'results/{G.TAG}.json']
    return dict(status='PASS_EXACT_A13_COMPACT_LOCAL_TRIANGLE_PAYMENT_FAILURE',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        gram=A.enc(M),edges=edges,black=black,white=white,physical_inertia=[13,1,0],conjugate_inertia=[14,0,0],
        old_A13_inverse_diagonal=str(F(12,7)),determinant=determinant.enc(),last_schur_pivot=schur.enc(),
        absent_from_saved_fixed_discriminant_classes=True,integral_Cartan_entries=196,deletions=deletions,corners=corners,
        center_epsilon=str(epsilon),strict_center=[x.enc() for x in center],strict_center_norm=A.pair(M,center,center).enc(),
        aggregate_payment=str(total),required_payment='2',
        scope='A genuine compact local triangle with a coherent integral reflection chamber and an admissible two-place quadratic form. It disproves coverage of the saved unrepresented-discriminant proposal on this additional star. It does not invalidate its restricted certificate, prove compact ambient realization, or exclude a repaired rule.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()))
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['aggregate_payment'],d['determinant'],flush=True)
