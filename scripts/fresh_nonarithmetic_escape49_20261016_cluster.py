#!/usr/bin/env python3
"""Three exact algebraic caps via prescribed contacts and one positive square root."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
import fresh_nonarithmetic_joint_regions_20261016_batch as B
D,A,R,C,K,Z,O=B.D,B.A,B.R,B.C,B.K,B.Z,B.O
I=R.L.I;ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_escape49_20261016_cluster'
CAP='results/fresh_nonarithmetic_common_escape_20261016_cap.json';CONE='results/fresh_nonarithmetic_common_escape_20261016_cone.json'
S=[0,1,2,3,5,6,7,30];ANCHOR=33
PROFILES=[{30:3,7:6,33:6},{3:3,7:5,33:6},{5:6,6:4,33:6}]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def dec(x):return tuple(A.B.db(a)for a in x)
def rnd(x):return B.rnd(x,256)
def iv(a):return rnd(sum((D.biv(x)*b for x,b in zip(a,D.BIV)),I()))
def sqrt(x):
 assert x.l>0
 def lo(a):return F(isqrt(a.numerator*(1<<512)//a.denominator),1<<256)
 return I(lo(x.l),lo(x.h)+F(1,1<<256))
def idot(a,b):return rnd(sum((x*y for x,y in zip(a,b)),I()))
def imv(M,v):return[idot(row,v)for row in M]
def theta(m):
 if m==2:return I()
 if m==3:return I(1)
 if m==4:return sqrt(I(2))
 if m==5:return iv(D.lift((O+R.BASIS[1])/2))
 if m==6:return sqrt(I(3))
 raise AssertionError(m)
def inverse01(x):
 assert not any(x[2:]);den=x[0]*x[0]-D.D[0]*x[1]*x[1];assert den
 out=(x[0]/den,-x[1]/den)+(Z,)*6;assert D.mul(x,out)==D.ONE;return out

def build():
 old=B.read(CAP);cone=B.read(CONE)
 Q=[[A.B.db(x)for x in row]for row in old['form']];N=[tuple(dec(x)for x in n)for n in old['normals']]
 H=tuple(dec(x)for x in old['center']);E=tuple(dec(x)for x in cone['strict_escape']['point']);cov=[D.mv(Q,n)for n in N]
 G=[[D.dot(N[i],cov[j])for j in S]for i in S];assert all(not any(x[1:])for row in G for x in row)
 GG=[[x[0]for x in row]for row in G];assert R.pd(GG)[0];J=C.inverse(GG)
 assert all(D.dot(N[i],cov[i])==D.lift(2*O)for i in S+[ANCHOR])
 def proj(v):
  coeff=[D.sumfield(D.scale(D.dot(v,cov[j]),a)for a,j in zip(row,S))for row in J]
  return tuple(D.sumfield(D.mul(c,N[j][t])for c,j in zip(coeff,S))for t in range(10))
 a=B.vecsub(N[ANCHOR],proj(N[ANCHOR]));assert all(D.pair(a,N[i],Q)==D.ZERO for i in S);qa=D.pair(a,a,Q);assert D.sign(qa)<0
 t=B.vecsub(H,proj(H));k=B.vecsub(t,B.vecscale(a,D.mul(D.pair(t,a,Q),inverse01(qa))))
 assert all(D.pair(k,N[i],Q)==D.ZERO for i in S+[ANCHOR]);qk=D.pair(k,k,Q);assert D.sign(qk)>0
 if D.sign(D.pair(k,E,Q))<0:k=tuple(D.scale(x,-O)for x in k)
 assert D.sign(D.pair(k,E,Q))>0
 QI=[[iv(D.lift(x))for x in row]for row in Q];NI=[[iv(x)for x in n]for n in N];HI=[iv(x)for x in H];EI=[iv(x)for x in E];covI=[imv(QI,n)for n in NI]
 ai=[iv(x)for x in a];ki=[iv(x)for x in k];Ji=[[iv(D.lift(x))for x in row]for row in J];qai,qki=iv(qa),iv(qk);assert qai.h<0<qki.l
 caps=[];records=[]
 for number,profile in enumerate(PROFILES):
  c=[-theta(profile.get(i,2))for i in S];coeff=imv(Ji,c);pS=[idot(coeff,[NI[j][t]for j in S])for t in range(10)]
  alpha=rnd((-theta(profile[ANCHOR])-idot(pS,covI[ANCHOR]))/qai)
  p=[rnd(x+alpha*y)for x,y in zip(pS,ai)];pp=rnd(idot(c,coeff)+alpha*alpha*qai);delta=rnd((2-pp)/qki);lam=sqrt(delta)
  r=[rnd(x+lam*y)for x,y in zip(p,ki)];vals=[idot(r,x)for x in covI];hp=idot(r,imv(QI,HI));ep=idot(r,imv(QI,EI));assert hp.h<0<ep.l
  rel=[]
  for i,g in enumerate(vals):
   if i in S+[ANCHOR]:rel.append(dict(wall=i,label=profile.get(i,2),kind='right'if profile.get(i,2)==2 else'finite',pairing_interval=g.enc(),exact_pairing=f'-2*cos(pi/{profile.get(i,2)})'if profile.get(i,2)!=2 else'0'))
   else:
    qi=iv(D.dot(N[i],cov[i]));gap=rnd(g*g-2*qi);assert g.h<0<gap.l,(number,i,g.enc(),gap.enc())
    rel.append(dict(wall=i,kind='dotted',pairing_interval=g.enc(),squared_product_minus_norm_product=gap.enc()))
  caps.append(r);records.append(dict(index=49+number,scout_index=[372,1815,1782][number],profile=[[i,m]for i,m in sorted(profile.items())],support_contacts=[c.enc()for c in c],projection_coefficients=[x.enc()for x in coeff],anchor_coefficient=alpha.enc(),particular=[x.enc()for x in p],particular_norm=pp.enc(),radical_square=delta.enc(),positive_radical=lam.enc(),normal_interval=[x.enc()for x in r],norm_exact=2,relations=rel,center_pairing=hp.enc(),cut_escape_pairing=ep.enc()))
 mutual=[]
 for i,j in combinations(range(3),2):
  g=idot(caps[i],imv(QI,caps[j]));assert g.h<-2;mutual.append(dict(caps=[49+i,49+j],pairing_interval=g.enc(),kind='dotted'))
 return old,cone,Q,N,H,E,J,a,qa,k,qk,records,caps,mutual

def run():
 old,cone,Q,N,H,E,J,a,qa,k,qk,records,caps,mutual=build()
 paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',CAP,CONE,'results/fresh_nonarithmetic_common_escape_20261016_cap_review.json','results/fresh_nonarithmetic_common_escape_20261016_cone_review.json','results/fresh_nonarithmetic_escape49_20261016_target_scout.json']
 return dict(status='PASS_EXACT_THREE_COMPATIBLE_ESCAPE_CUTTING_CAPS',sources={p:sha(p)for p in paths},ambient_dimension=9,retained_walls=49,resulting_walls=52,form=B.encode(Q),center=B.encode(H),cut_escape=B.encode(E),spherical_support=S,anchor=ANCHOR,inverse_support_gram=B.encode(J),projected_anchor=B.encode(a),projected_anchor_norm=B.encode(qa),positive_kernel=B.encode(k),kernel_norm=B.encode(qk),caps=records,mutual=mutual,actual_facets=52,
  exact_recipe='Let c_i=-2cos(pi/m_i) on S, with omitted labels2. pS=NS*J*c; alpha=(-2cos(pi/m33)-<pS,n33>)/qa; p=pS+alpha*a; r=p+sqrt((2-<p,p>)/qk)*k. The positive square root is used. Exact field identities certify a perpendicular to S, k perpendicular to S and33, and qa<0<qk. Thus every prescribed contact and r²=2 holds identically. All unspecified contacts and mutual pairs have strict rational dotted certificates.',
  facet_proof='The inherited center is strict for all52 normals. Every pair is acute. For each normal n of positive squareq, qH-<H,n>n has negative norm, future orientation, exactly one zero wall pairing and all other wall pairings negative. Thus every wall remains an actual facet.',
  scope='Exact compatible three-cap collection cutting the renewed49-wall escape. This uses algebraic construction recipes and outward rational intervals, without expansion into a large primitive field. Complete52-wall cone and volume are separate; no finite-volume claim follows from this local certificate.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved exact three-cap replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
