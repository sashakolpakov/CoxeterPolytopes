#!/usr/bin/env python3
"""Exact positive finite-vertex decomposition of the renewed53 escape projection."""
import argparse,json
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_escape52_20261016_cap as C53
U=C53.U;B,D,A,I=U.B,U.D,U.A,U.I
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_h3_escape_projection_20261016'
CAP='results/fresh_nonarithmetic_escape52_20261016_cap.json';CONE='results/fresh_nonarithmetic_escape52_20261016_cone.json'
S=[0,1,2,3,5,6];TEST=[7,8,30,33]
SUPPORTS=[S+[7,8,11],S+[7,10,33],S+[8,10,44],S+[8,12,34]]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def positive_definite(M):
 a=[list(row)for row in M];pivots=[]
 for k in range(len(a)):
  q=a[k][k];assert q.l>0,('not certified PD',k,q.enc());pivots.append(q.enc())
  for i in range(k+1,len(a)):
   for j in range(i,len(a)):
    a[i][j]=a[j][i]=U.rnd(a[i][j]-a[i][k]*a[j][k]/q)
 return pivots
def run():
 old=B.read(U.CAP);cl=B.read(C53.CL);cap=B.read(CAP);co=B.read(CONE)
 Q=[[A.B.db(x)for x in row]for row in old['form']];QI=[[U.iv(D.lift(x))for x in row]for row in Q]
 N0=[tuple(U.dec(x)for x in row)for row in old['normals']];N=[[U.iv(x)for x in n]for n in N0]+[[I(*x)for x in row['normal_interval']]for row in cl['caps']]+[[I(*x)for x in cap['normal_interval']]]
 assert len(N)==53
 H=[U.iv(U.dec(x))for x in old['center']];E=[U.iv(U.dec(x))for x in co['strict_escape']['point']];cov=[U.imv(QI,n)for n in N];hc=U.imv(QI,H);hh=U.idot(H,hc);assert hh.h<0
 pair=lambda x,y:U.idot(x,U.imv(QI,y))
 SG=[[U.idot(N[i],cov[j])for j in S]for i in S];spd=positive_definite(SG);SJ,sp=C53.inverse(SG);se=[U.idot(E,cov[i])for i in S];d=U.imv(SJ,se);spoint=C53.combine(d,[N[i]for i in S]);p=[U.rnd(x-y)for x,y in zip(E,spoint)];L=-U.idot(p,hc);assert L.l>0;p=[U.rnd(x/L)for x in p]
 vertices=[];V=[]
 for support in SUPPORTS:
  assert support in co['complete_census']['finite_vertices']
  G=[[U.idot(N[i],cov[j])for j in support]for i in support];pd=positive_definite(G);J,pivots=C53.inverse(G);hs=[U.idot(H,cov[i])for i in support];cs=U.imv(J,hs);v=[U.rnd(x-y)for x,y in zip(H,C53.combine(cs,[N[i]for i in support]))];qv=U.rnd(hh-U.idot(hs,cs));vh=U.idot(v,hc);assert qv.h<0 and vh.h<0;scale=-vh;v=[U.rnd(x/scale)for x in v];strict=[]
  for j,c in enumerate(cov):
   if j in support:continue
   z=U.idot(v,c);assert z.h<0,('nonincident sign',support,j,z.enc());strict.append(dict(wall=j,pairing=z.enc()))
  V.append(v);vertices.append(dict(support=support,positive_gram_pivots=pd,inverse_pivots=pivots,point_interval=[x.enc()for x in v],unscaled_norm=qv.enc(),unscaled_future_chart=vh.enc(),normalization=scale.enc(),nonincident_pairings=strict))
 # The four tests are restricted covectors on S-perp. Their invertible
 # pairing matrix with the four independent finite vertices proves that
 # they form a full dual basis. Thus the resulting four equations prove
 # the full vector identity, not merely four approximate coordinates.
 M=[[U.idot(v,cov[i])for v in V]for i in TEST];MI,pivots=C53.inverse(M);b=[U.idot(p,cov[i])for i in TEST];weights=U.imv(MI,b);assert all(w.l>0 for w in weights)
 total=sum(weights,I());assert total.l<=1<=total.h
 qproj=pair(p,p);assert qproj.h<0
 paths=[f'scripts/{TAG}.py',f'scripts/{C53.TAG}.py',f'scripts/{U.TAG}.py',U.CAP,C53.CL,CAP,CONE,'results/fresh_nonarithmetic_escape52_20261016_cap_review.json','results/fresh_nonarithmetic_escape52_20261016_cone_review.json']
 for path in paths:
  if path.endswith('.json'):B.read(path)
 return dict(status='PASS_EXACT_RENEWED53_ESCAPE_POSITIVE_FINITE_VERTEX_PROJECTION',sources={p:sha(p)for p in paths},retained_wall_count=53,right_support=S,support_positive_pivots=spd,support_inverse_pivots=sp,escape=co['strict_escape']['point'],support_projection_coefficients=[x.enc()for x in d],projected_future_scale=L.enc(),normalized_projected_escape=[x.enc()for x in p],projected_norm=qproj.enc(),vertices=vertices,dual_test_walls=TEST,pairing_matrix=[[x.enc()for x in row]for row in M],inverse_pivots=pivots,positive_weights=[x.enc()for x in weights],weight_sum=total.enc(),
 exact_identity='Let v_j be the future finite corner projections of H, each normalized by -<H,v_j>=1. Let p=(E-N_S G_S^-1 <N_S,E>)/L, L=-<H,E-N_S G_S^-1<N_S,E>>0. The invertible4x4 matrix M_ij=<n_TEST[i],v_j> defines lambda=M^-1(<n_TEST[i],p>). Since S has rank6 in Lorentz dimension10 and M is invertible, the restricted test covectors are a basis of(S-perp)*. Consequently p=sum lambda_j v_j exactly. All lambda_j>0, and normalization gives sum lambda_j=1.',
 persistence_proof='Let r be a spacelike acute new normal with <r,n_i><=0 for all retained walls, and assume there exists a future timelike X strictly satisfying all retained inequalities and <X,r><0. For each saved finite vertex v with positive-definite acute Gram G, write X=a v+sum b_i n_i. Here a=<X,v>/<v,v>>0, and b=G^-1(<X,n_i>) is entrywise nonpositive because G is a positive-definite Z-matrix. Therefore <X,r>=a<v,r>+sum b_i<n_i,r>>=a<v,r>. Since <X,r><0, one has <v,r><0. This proof does not require that the old chosen center H is retained by r.',
 conclusion='Every spacelike acute new wall admitting a common strict future interior with the53 retained walls and perpendicular to all six roots S has <r,E>=L*sum lambda_j<r,v_j><0. Hence no such wall cuts or touches the renewed53 escape. This covers all real coefficients, normalization choices and contact labels. Any finite family of such additions retains that strict spacelike escape.',
 scope='This is a fixed-support escape obstruction, not a no-completion theorem for the53-wall chamber. A successful global completion must include a wall nonright to at least one root of S. Closing the H3 face by right-S walls may still be useful in combination with separate nonright cutters. No unsupported sign assertion is made for acute normals whose halfspaces have no common future interior.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved finite-vertex projection gate')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['positive_weights'])
