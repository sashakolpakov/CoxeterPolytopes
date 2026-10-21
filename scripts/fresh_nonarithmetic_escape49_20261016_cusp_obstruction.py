#!/usr/bin/env python3
"""Positive affine G2 kernel and an orthogonal retained pi/5 ridge."""
import argparse,json
from pathlib import Path
from hashlib import sha256
import fresh_nonarithmetic_escape49_20261016_cluster as U
B,D,A,R,C,K,Z,O=U.B,U.D,U.A,U.R,U.C,U.K,U.Z,U.O
I=U.I;ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_escape49_20261016_cusp_obstruction'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 old,cone,Q,N,H,E,J,a,qa,k,qk,records,caps,mutual=U.build()
 cov=[D.mv(Q,n)for n in N];QI=[[U.iv(D.lift(x))for x in row]for row in Q];HI=[U.iv(x)for x in H]
 phi=(O+R.BASIS[1])/2
 assert all(D.pair(N[i],N[i],Q)==D.lift(2*O)for i in[0,1,30,33])
 assert D.pair(N[0],N[1],Q)==D.lift(-phi)
 assert D.pair(N[30],N[33],Q)==D.ZERO
 assert all(D.pair(N[i],N[j],Q)==D.ZERO for i in[0,1]for j in[30,33])
 assert {j:dict(records[0]['profile']).get(j,2)for j in[0,1,30,33]}=={0:2,1:2,30:3,33:6}
 # G*[1,sqrt(3),2]=0, all coefficients positive, and H.Nu<0.
 s3=U.sqrt(I(3));nu=[U.rnd(U.iv(N[30][j])+s3*U.iv(N[33][j])+2*caps[0][j])for j in range(10)]
 hp=U.idot(nu,U.imv(QI,HI));assert hp.h<0
 all_normals=[[U.iv(x)for x in n]for n in N]+caps
 pairings=[];zeros=[]
 for i,n in enumerate(all_normals):
  if i<49:
   g30=D.pair(N[30],N[i],Q);g33=D.pair(N[33],N[i],Q)
   # Diagonal entries of affine generators are handled by exact kernel identity.
   if i in[30,33]:kind='zero'
   else:
    assert D.sign(g30)<=0 and D.sign(g33)<=0
    kind='zero'if g30==g33==D.ZERO and records[0]['relations'][i]['kind']=='right'else'strict'
  elif i==49:kind='zero'
  else:kind='strict' # Its pairing with r49 is strictly negative; all other summands are nonpositive.
  v=U.idot(nu,U.imv(QI,n))
  if kind=='strict':assert v.h<0,(i,v.enc())
  else:zeros.append(i);assert v.l<=0<=v.h
  pairings.append(dict(wall=i,kind=kind,interval=v.enc()))
 assert zeros==[0,1,2,3,5,6,30,33,49]
 # Literal actual pi/5 ridge witness in the52 chamber, separate from the ideal chart.
 J2=C.inverse([[2*O,-phi],[-phi,2*O]])
 h=[D.pair(H,N[j],Q)for j in[0,1]]
 coeff=[D.sumfield(D.scale(x,c)for x,c in zip(h,row))for row in J2]
 assert all(D.sign(x)<0 for x in coeff)
 p=B.vecsub(H,tuple(D.sumfield(D.mul(coeff[i],N[i][j])for i in range(2))for j in range(10)))
 pn=D.pair(p,p,Q);ph=D.pair(p,H,Q);assert D.sign(pn)<0 and D.sign(ph)<0
 pv=[D.pair(p,n,Q)for n in N];assert pv[0]==pv[1]==D.ZERO and all(D.sign(x)<0 for x in pv[2:])
 pI=[U.iv(x)for x in p];newpv=[U.idot(pI,U.imv(QI,n))for n in caps];assert all(x.h<0 for x in newpv)
 paths=[f'scripts/{TAG}.py',f'scripts/{U.TAG}.py',f'results/{U.TAG}.json',U.CAP,
 'results/fresh_h9_fixed_core_cusp_obstruction_20261013.md',
 'results/fresh_h20_weighted_binary_20261014_dotted_completion_review.md']
 return dict(status='PASS_PERMANENT_NONCRYSTALLOGRAPHIC_CUSP_IN_FIRST52_CLUSTER',sources={p:sha(p)for p in paths},
 retained_wall_indices=[0,1,30,33,49],affine_indices=[30,33,49],affine_normalized_gram=[['2','0','-1'],['0','2','-sqrt(3)'],['-1','-sqrt(3)','2']],positive_kernel=['1','sqrt(3)','2'],null_recipe='nu=n30+sqrt(3)*n33+2*r49',null_norm_exact=0,null_interval=[x.enc()for x in nu],future_center_pairing=hp.enc(),null_pairings=pairings,incident_walls=zeros,
 forbidden_edge=dict(indices=[0,1],norms=[2,2],pairing=B.encode(D.lift(-phi)),label=5,euclidean_link_dimension=8,reflection_product_linear_trace='5+phi'),
 actual_ridge=dict(point=B.encode(p),norm=B.encode(pn),future_center_pairing=B.encode(ph),old_pairings=B.encode(pv),new_pairings=[x.enc()for x in newpv]),
 theorem='Every finite full-dimensional acute-compatible Coxeter extension retaining these five walls has infinite volume. Positive affine coefficients preserve nu in the closed chamber; a strict center of the proposed extension makes nu+epsilon*center feasible timelike. A finite-volume end there would have a bounded Euclidean Coxeter link. Orthogonal projection preserves the inherited pi/5 facets/ridge in that link. Their reflection product has irrational linear trace5+phi, impossible for a group preserving a full translation lattice.',
 scope='The exact five-wall anchors, positive future null vector and actual52-wall pi/5 ridge are checked here. The general completion implication uses the separately cited frozen cusp proof; it is independent of the complete52 cone. This does not obstruct deleting/replacing one of the five anchors or the untouched49 chamber.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved permanent cusp replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
