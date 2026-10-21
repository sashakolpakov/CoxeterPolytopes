#!/usr/bin/env python3
"""Exact conic Farkas certificate: every acute wall perpendicular to33 missesE."""
import argparse,json
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_affine_filtered_20261016_cluster as U
B,D,A,R,C,K,Z,O=U.B,U.D,U.A,U.R,U.C,U.K,U.Z,U.O
I=U.I;ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_pi10_ridge_20261016_gate';CL=f'results/{U.TAG}.json';CO='results/fresh_nonarithmetic_affine_filtered_20261016_cone.json'
BASIS=[0,10,11,13,19,27,30,33,49,51]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 old=B.read(U.CAP);cl=B.read(CL);co=B.read(CO)
 N=[[U.iv(U.dec(x))for x in n]for n in old['normals']]+[[I(*x)for x in r['normal_interval']]for r in cl['caps']]
 E=[U.iv(U.dec(x))for x in co['strict_escape']['point']]
 # Every entry denotes an exact algebraic quantity via bound construction recipes.
 # Nonzero interval pivots certify unique exact solution M*c=E; no rounded equality is used.
 M=[[N[i][j]for i in BASIS]+[E[j]]for j in range(10)];pivots=[];swaps=[]
 for k in range(10):
  choices=[r for r in range(k,10)if M[r][k].l>0 or M[r][k].h<0];assert choices,(k,'unresolved pivot')
  p=max(choices,key=lambda r:min(abs(M[r][k].l),abs(M[r][k].h)));M[k],M[p]=M[p],M[k];swaps.append([k,p]);q=M[k][k];pivots.append(q.enc())
  M[k]=[U.rnd(x/q)for x in M[k]]
  for r in range(10):
   if r==k:continue
   f=M[r][k];M[r]=[U.rnd(x-f*y)for x,y in zip(M[r],M[k])]
 coeff=[M[i][-1]for i in range(10)]
 for i,c in zip(BASIS,coeff):
  if i==33:assert c.h<0
  else:assert c.l>0,(i,c.enc())
 # An independent interval residual enclosure is a sanity check only; equality is the inverse definition.
 residual=[U.rnd(sum((c*N[i][j]for i,c in zip(BASIS,coeff)),I())-E[j])for j in range(10)];assert all(x.l<=0<=x.h for x in residual)
 paths=[f'scripts/{TAG}.py',f'scripts/{U.TAG}.py',U.CAP,CL,CO]
 return dict(status='PASS_ALL_ACUTE_WALLS_RIGHT_TO33_MISS_RENEWED52_ESCAPE',sources={p:sha(p)for p in paths},retained_walls=52,escape=co['strict_escape']['point'],basis=BASIS,coefficient_intervals=[x.enc()for x in coeff],nonzero_pivot_intervals=pivots,row_swaps=swaps,residual_enclosures=[x.enc()for x in residual],
  exact_identity='Let M have columns the exact normals with indices BASIS. The nonzero-pivot certificates prove det(M)!=0. Define c=M^-1 E. Outward rational interval Gaussian elimination proves c33<0 and every other coefficient strictly positive. Thus E=sum_i c_i n_i is an exact algebraic identity.',
  theorem='For every real normal r acute-compatible with all52 retained walls and with <r,n33>=0, <r,E>=sum_(i!=33)c_i<r,n_i><=0. In fact the inequality is strict for nonzero r: equality forces r orthogonal to all10 independent basis normals. No norm normalization, Coxeter alphabet, field, new-interior-point premise or finite-corner premise is needed.',
  scope='This excludes escape-cutting additions perpendicular to old33, including every cap used to close any of the three new pi/10 ridges. It does not exclude useful noncutting ridge-closing walls, or caps with a nonright/dotted contact to33. It is not a no-completion theorem for the52 chamber.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved exact all-real perpendicular33 gate')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
