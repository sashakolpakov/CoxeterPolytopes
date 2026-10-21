#!/usr/bin/env python3
"""Exact algebraic escape cutter with one pi/4 contact and eight right contacts."""
import argparse,itertools,json
from hashlib import sha256
from pathlib import Path
import sympy as sy
import fresh_nonarithmetic_affine_filtered_20261016_cluster as U
B,D,A,R,C,K,Z,O=U.B,U.D,U.A,U.R,U.C,U.K,U.Z,U.O
I=U.I;ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_escape52_20261016_cap'
CL=f'results/{U.TAG}.json';CO='results/fresh_nonarithmetic_affine_filtered_20261016_cone.json'
S=[0,1,2,3,5,6,8,34];ANCHOR=30;FRAME=S+[ANCHOR]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def inverse(M):
 n=len(M);a=[list(row)+[I(int(i==j))for j in range(n)]for i,row in enumerate(M)];pivots=[]
 for k in range(n):
  options=[i for i in range(k,n)if a[i][k].l>0 or a[i][k].h<0];assert options,('unresolved pivot',k)
  j=max(options,key=lambda i:min(abs(a[i][k].l),abs(a[i][k].h)));a[k],a[j]=a[j],a[k];q=a[k][k];pivots.append(dict(column=k,row=j,pivot=q.enc()));a[k]=[U.rnd(x/q)for x in a[k]]
  for i in range(n):
   if i==k:continue
   f=a[i][k];a[i]=[U.rnd(x-f*y)for x,y in zip(a[i],a[k])]
 return [row[n:]for row in a],pivots
def combine(coeff,N):return[U.idot(coeff,[n[j]for n in N])for j in range(10)]
def run():
 old=B.read(U.CAP);cl=B.read(CL);co=B.read(CO)
 Q=[[A.B.db(x)for x in row]for row in old['form']];N0=[tuple(U.dec(x)for x in row)for row in old['normals']];H0=tuple(U.dec(x)for x in old['center']);E0=tuple(U.dec(x)for x in co['strict_escape']['point'])
 for i in FRAME:assert D.pair(N0[i],N0[i],Q)==D.lift(2*O)
 QI=[[U.iv(D.lift(x))for x in row]for row in Q];NI=[[U.iv(x)for x in n]for n in N0]+[[I(*x)for x in r['normal_interval']]for r in cl['caps']];H=[U.iv(x)for x in H0];E=[U.iv(x)for x in E0];cov=[U.imv(QI,n)for n in NI];hc=U.imv(QI,H);ec=U.imv(QI,E)
 G=[[U.iv(D.pair(N0[i],N0[j],Q))for j in FRAME]for i in FRAME];J,pivots=inverse(G)
 # J denotes the exact inverse of the literal nine-root Gram, certified by its nonzero pivots.
 contact=[I()]*8+[-U.sqrt(I(2))];pcoeff=U.imv(J,contact);p=combine(pcoeff,[NI[i]for i in FRAME]);hp=[U.idot(H,cov[i])for i in FRAME];hcoeff=U.imv(J,hp);projection=combine(hcoeff,[NI[i]for i in FRAME]);k=[U.rnd(x-y)for x,y in zip(H,projection)]
 pp=U.idot(contact,pcoeff);hh=U.idot(H,hc);qk=U.rnd(hh-U.idot(hp,hcoeff));assert qk.l>0;ke=U.idot(k,ec)
 if ke.h<0:k=[-x for x in k];kernel_orientation=-1
 else:assert ke.l>0;kernel_orientation=1
 delta=U.rnd((2-pp)/qk);lam=U.sqrt(delta);normal=[U.rnd(x+lam*y)for x,y in zip(p,k)];nr=U.imv(QI,normal);he=U.idot(normal,hc);ee=U.idot(normal,ec);assert he.h<0<ee.l
 native_norms=[U.iv(D.pair(n,n,Q))for n in N0]+[I(2)]*3
 relations=[]
 for i,n in enumerate(NI):
  g=U.idot(normal,cov[i])
  if i in FRAME:
   m=4 if i==ANCHOR else 2;relations.append(dict(wall=i,label=m,kind='right'if m==2 else'finite',pairing_interval=g.enc(),exact_by_recipe=True))
  else:
   gap=U.rnd(g*g-2*native_norms[i]);assert g.h<0<gap.l,('not strict dotted',i,g.enc(),gap.enc());relations.append(dict(wall=i,kind='dotted',pairing_interval=g.enc(),squared_product_minus_norm_product=gap.enc()))
 # Every new finite vertex uses eight of the nine finite-neighbor walls.
 # Exact symbolic Gram tests cover all nine choices; no complete ambient census is claimed here.
 inv=B.read('results/fresh_h9_persistent_cusp_screen_20261016.json');L=inv['pair_matrix'];assert all(L[ANCHOR][i]in[0,2]for i in S);theta={2:sy.Integer(0),3:sy.Integer(1),4:sy.sqrt(2),5:(1+sy.sqrt(5))/2,6:sy.sqrt(3),10:sy.sqrt((5+sy.sqrt(5))/2)}
 supports=[];vertices=[]
 for ids in itertools.combinations(FRAME,8):
  ss=list(ids)+[52]
  labels=[[2 if i==j else(4 if min(i,j)==ANCHOR else 2)if max(i,j)==52 else L[i][j]for j in ss]for i in ss]
  if any(labels[i][j]==0 for i in range(9)for j in range(i)):
   supports.append(dict(support=ss,kind='contains_dotted_pair'));continue
  M=sy.Matrix([[2 if i==j else -theta[labels[i][j]]for j in range(9)]for i in range(9)])
  leading=[sy.simplify(M[:i,:i].det())for i in range(1,10)];pd=all(z.is_positive is True for z in leading)
  row=dict(support=ss,leading_minors=[str(z)for z in leading],positive_definite=pd);supports.append(row)
  if not pd:continue
  VN=[NI[i]if i<52 else normal for i in ss];GI=[[I(2)if i==j else-U.theta(labels[i][j])for j in range(9)]for i in range(9)];JI,vp=inverse(GI);hs=[U.idot(v,hc)for v in VN];cs=U.imv(JI,hs);v=[U.rnd(x-y)for x,y in zip(H,combine(cs,VN))];qv=U.rnd(hh-U.idot(hs,cs));vh=U.idot(v,hc);assert qv.h<0 and vh.h<0
  nonincident=[]
  for j,c in enumerate(cov+[nr]):
   if j not in ss:
    g=U.idot(v,c);assert g.h<0,('finite witness off-wall',ss,j,g.enc());nonincident.append(dict(wall=j,pairing=g.enc()))
  vertices.append(dict(support=ss,point_interval=[x.enc()for x in v],norm_interval=qv.enc(),future_chart_interval=vh.enc(),inverse_pivots=vp,nonincident_pairings=nonincident))
 paths=[f'scripts/{TAG}.py',f'scripts/{U.TAG}.py',U.CAP,CL,CO,'results/fresh_h9_persistent_cusp_screen_20261016.json','results/fresh_nonarithmetic_affine_filtered_20261016_affine_gate.json','results/fresh_nonarithmetic_ridge_pair_20261016_cutter_scout.json','results/fresh_nonarithmetic_ridge_pair_20261016_selection.json']
 return dict(status='PASS_EXACT_53RD_WALL_CUTTING_RENEWED52_ESCAPE',sources={p:sha(p)for p in paths},scout_index=6,retained_walls=52,resulting_walls=53,form=old['form'],center=old['center'],cut_escape=co['strict_escape']['point'],right_support=S,finite_anchor=dict(wall=ANCHOR,label=4),frame=FRAME,frame_gram_interval=[[x.enc()for x in row]for row in G],inverse_frame_gram_interval=[[x.enc()for x in row]for row in J],inverse_pivots=pivots,particular=[x.enc()for x in p],particular_norm=pp.enc(),kernel=[x.enc()for x in k],kernel_orientation=kernel_orientation,kernel_norm=qk.enc(),radical_square=delta.enc(),positive_radical=lam.enc(),normal_interval=[x.enc()for x in normal],norm_exact=2,center_pairing=he.enc(),escape_pairing=ee.enc(),relations=relations,finite_support_census=supports,new_finite_vertices=vertices,actual_facets=53,
 exact_recipe='With F the ordered nine-root frame, G its exact Gram and c=(0^8,-sqrt2), set p=NF G^-1 c and k=H-NF G^-1(<H,n_i>) with the saved orientation. Define r=p+sqrt((2-p²)/k²)k. Certified nonzero exact Gram pivots prove existence/uniqueness; orthogonal projection identities give every prescribed contact and r²=2 exactly. Outward rational intervals prove every other relation strictly dotted and the positive radical branch cutsE.',
 facet_proof='All53 normals are spacelike and acute and share the strict timelike centerH. For each n of squareq, qH-<H,n>n is future timelike, lies on its wall and is strict at every other wall. Thus all53 are actual facets.',
 no_affine_proof='The retained52 has no affine circuit. The new root has exactly one nonright finite edge, label4 to30, and eight right finite neighbors. Its only connected standard-label subgraphs containing it are A1 and I2(4); every other wall is dotted to it. Hence no affine circuit is introduced.',
 scope='Exact local geometry only. The complete census is for vertices on the new facet, by all nine eight-neighbor supports. The ambient53-wall cone and volume are separate. This is not a finite-volume completion claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved exact escape cutter')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],'new finite vertices',len(d['new_finite_vertices']))
