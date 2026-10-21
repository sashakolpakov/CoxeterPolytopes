#!/usr/bin/env python3
"""Exact common-escape gate for all Coxeter walls bounding the retained S75 H2.

The numerical scout only proposes a21-column basis. This checker solves and
verifies its identity exactly; all real contact magnitudes are covered.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from flint import fmpq_mat,fmpq
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_s75_completion_20261014_gate'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
PRE='results/fresh_h20_e7_bridge_face_20261014_exact.json'
CURRENT='results/fresh_h20_bridge_edge_closure_20261014.json'
FACE='results/fresh_h20_bridge_faces_20261014_review.json'
RAD=[1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]
FORCED=[3,4,5,6,7,8,13,14,15,75]
POSITIVE=[0,2,18,31,65,79,84,85,88,89,92]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(tuple,list)):return[enc(v)for v in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x
def run():
    old,pre,cur,face=[json.loads((ROOT/p).read_text())for p in[BASE,PRE,CURRENT,FACE]]
    for doc in[old,pre,cur,face]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert old['field']['base_square_root_basis']==RAD
    f=next(row for row in face['faces']if 75 in row['support']);S=f['support']
    G=[[F(q)for q in r]for r in f['rational_gram']];J=[[F(q)for q in r]for r in f['inverse']]
    assert len(S)==18 and all(G[i][i]==2 for i in range(18))
    assert all(q>=0 for row in J for q in row)
    assert all(sum(G[i][k]*J[k][j]for k in range(18))==int(i==j)for i in range(18)for j in range(18))
    assert [S[i]for i in range(18)if J[i][i]>=2]==FORCED
    ids=FORCED+POSITIVE;assert len(ids)==len(set(ids))==21
    vectors=[]
    for i in ids:
        v=[]
        for a,b in old['normals'][i]:
            assert all(F(t)==0 for t in b)
            v.append(tuple(map(F,a)))
        vectors.append(v)
    assert cur['strict_escape']['point']==pre['strict_escape']['point']
    E=[]
    for ab,cd in cur['strict_escape']['point']:
        a,b=ab;c,d=cd
        assert all(F(x)==0 for x in b+c+d+a[1:])
        E.append(F(a[0]))
    Q=[[F(x)for x in row]for row in old['form']]
    qE=sum(E[i]*Q[i][j]*E[j]for i in range(21)for j in range(21))
    assert qE==F(20181174579467302741,10**18)>0
    # Regular representation of the exact16-element field. Multiplication
    # sqrt(RAD[k])*sqrt(RAD[l]) = RAD[k&l]*sqrt(RAD[k xor l]).
    size=21*16;M=fmpq_mat(size,size);rhs=fmpq_mat(size,1)
    for i in range(21):
        rhs[16*i,0]=fmpq(str(E[i]))
        for j,v in enumerate(vectors):
            for k,a in enumerate(v[i]):
                if a:
                    for ell in range(16):M[16*i+(k^ell),16*j+ell]=fmpq(str(a*RAD[k&ell]))
    solution=M.solve(rhs,algorithm='dixon');assert M*solution==rhs
    assert M.rank()==size
    coefficients=[tuple(F(str(solution[16*j+k,0]))for k in range(16))for j in range(21)]
    signs=[B.sign(v)for v in coefficients]
    assert all(s>0 for s in signs[len(FORCED):])
    # Independently recompute the physical coefficient identity using the
    # inherited Fraction radical product, rather than only FLINT matrix code.
    reconstructed=[]
    for k in range(21):
        total=B.Z
        for c,v in zip(coefficients,vectors):total=B.add(total,B.mul(c,v[k]))
        target=(E[k],)+(F(0),)*15;assert total==target
        reconstructed.append(total)
    paths=[BASE,PRE,CURRENT,FACE,f'scripts/{TAG}.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py']
    return enc(dict(status='PASS_ALL_S75_BOUNDARY_CAPS_STRICTLY_MISS_CURRENT_ESCAPE',sources={p:sha(p)for p in paths},field_square_root_basis=RAD,support=S,support_gram=G,support_inverse=J,mandatory_right_walls=FORCED,mandatory_right_inverse_diagonals={i:J[S.index(i)][S.index(i)]for i in FORCED},positive_wall_indices=POSITIVE,identity_coefficients=[dict(wall=i,coefficient=c,sign=s,role='free_forced_right'if j<len(FORCED)else'strictly_positive_acute_wall')for j,(i,c,s)in enumerate(zip(ids,coefficients,signs))],physical_basis_rank=21,regular_representation_rank=size,escape_point=E,escape_norm=qE,identity='E = sum(forced μ_i N_i) + sum(positive λ_j N_j); all λ_j>0 and the21 normals span the ambient space',geodesic_gate='For r²=q>0, old support norms2 and acute Coxeter pairings, any nonright support contact has (r.N_i)²>=q/2. Entrywise nonnegative inverse makes projection cost>=J_ii*q/2. IfJ_ii>=2, projected square<=0, so a boundary geodesic requires r.N_i=0.',strict_pairing_gate='If r is nonzero, acute to the11 positive walls and right to all10 forced walls, then r.E<0. Equality would force orthogonality to all21 independent normals.',scope='Every acute Coxeter cap inducing a geodesic boundary on retainedS75 strictly misses the current exactescape, for arbitrary real normal/radius/field and every label. Any finite extension using only such caps retains a strictspacelikeescape and has infinitevolume. Doesnot exclude caps not boundingS75, mixtures with them, or changed retainedwalls.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'positive=',out['positive_wall_indices'],flush=True)
