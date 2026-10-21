#!/usr/bin/env python3
"""Exact812-facet incidence chamber and23 mutually compatible pi/7 caps.

The new normals are p*u+sqrt(rad)*v with rational u,v, p proportional to
T=2cos(pi/7), and rad in Q(T). Polynomial identities prove finite angles;
outward-rounded integer intervals prove all strict inequalities.
"""
import argparse,hashlib,json,math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from fresh_h22_incidence_core_20261013 import inertia
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_incidence_prism_20261013';OUT='results/'+PRE+'_verify.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';BANK='data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json';REC='data/enumerations/'+PRE+'_recipes.json';CORE='results/fresh_h22_incidence_core_20261013.json'
S=10**70
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def rat(x):
 x=F(x)*S;return x.numerator//x.denominator,-((-x.numerator)//x.denominator)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def mul(a,b):
 z=[x*y for x in a for y in b];return min(z)//S,-((-max(z))//S)
def sqrt(a):
 assert a[0]>=0;a0=math.isqrt(a[0]*S);a1=math.isqrt(a[1]*S);return a0,a1+(a1*a1<a[1]*S)
def isum(v):
 out=(0,0)
 for a in v:out=add(out,a)
 return out
def scale(a,c):return mul(a,rat(c))
def tau():
 lo=18*S//10;hi=19*S//10
 def poly(x):return x**3-x*x*S-2*x*S*S+S**3
 assert poly(lo)<0<poly(hi)
 while hi-lo>1:
  mid=(lo+hi)//2
  if poly(mid)<0:lo=mid
  else:hi=mid
 # f'=3x²-2x-2>0 on[1.8,1.9]; this isolates its largest real root.
 return lo,hi
T=tau();ZERO=(F(0),)*3;ONE=(F(1),F(0),F(0));TWO=(F(0),F(0),F(1))
def ka(a,b):return tuple(x+y for x,y in zip(a,b))
def ks(a,c):return tuple(x*c for x in a)
def km(a,b):
 z=[F(0)]*5
 for i in range(3):
  for j in range(3):z[i+j]+=a[i]*b[j]
 for i in(4,3):z[i-1]+=z[i];z[i-2]+=2*z[i];z[i-3]-=z[i]
 return tuple(z[:3])
def ki(a):return isum([rat(a[0]),scale(T,a[1]),scale(mul(T,T),a[2])])
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mv(M,v):return[dot(r,v)for r in M]
def pair(M,u,v):return dot(u,mv(M,v))
def angle_identity(L,B,rad,q,n,m,p):
 if m==2:E=ZERO
 elif m in(3,4,6):E=ks(ONE,q*n*{3:F(1,4),4:F(1,2),6:F(3,4)}[m])
 else:assert m==7;E=ks(TWO,q*n/4)
 Y=ka(ka(E,ks(TWO,-L*L)),ks(rad,-B*B))
 if not L*B:assert Y==ZERO
 else:
  assert km(Y,Y)==ks(km(TWO,rad),4*L*L*B*B)
  yi=ki(Y);assert yi[0]>0 if L*B>0 else yi[1]<0
 if m==2:assert p[0]<=0<=p[1]
 else:assert p[1]<0
def dotted(p,q,n):assert p[1]<0 and p[1]*p[1]>q*n*S*S
def setup():
 raw=read(RAW)['coxeter_input'];Q=[[-int(x)for x in row]for row in raw['basis_gram']];M=[r+[0]for r in Q]+[[0]*22+[2]];W=[[int(x)for x in row]+[0]for row in raw['wall_coordinates_as_columns_transposed']];core=read(CORE);points,lines=core['incidence_parts'];vp=core['old_projections'][1]
 chosen=[c for c in read(BANK)['candidates']if c['twice_height_square']==4];C=[]
 for c in chosen:
  a=c['A'];mask=c['mask'];u=[mask>>i&1 for i in range(21)];assert sum(u)==c['weight'] and a*a-sum(u)+4==4
  twice=[a*vp[j]-sum(u[i]*W[points[i]][j]for i in range(21))for j in range(22)];assert twice==c['twice_old_projection']and all(x%2==0 for x in twice);C.append([x//2 for x in twice])
 assert len(C)==301 and len(set(map(tuple,C)))==301
 R=W+[r+[1]for r in C]+[r+[-1]for r in C];h=core['center']+[0];return M,R,h
def run():
 M,R,h=setup();assert inertia(M)==[22,1,0]and pair(M,h,h)==-14
 m=np.array(M,dtype=np.int64);r=np.array(R,dtype=np.int64);assert 23**2*int(abs(m).max())*int(abs(r).max())**2<2**63;G=r@m@r.T;norms=np.diag(G).tolist();assert Counter(norms)=={2:644,4:168}
 assert all(x<0 for x in r@(m@np.array(h,dtype=np.int64)))
 counts=Counter()
 for i in range(812):
  for j in range(i):
   p=int(G[i,j]);assert p<=0;c=F(p*p,norms[i]*norms[j]);assert c in(0,F(1,4),F(1,2),F(3,4))or c>=1;counts[str(c)if c<1 else'parallel'if c==1 else'dotted']+=1
 cap_intervals=[];cap_norms=[];records=[]
 for rec in read(REC)['recipes']:
  u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm']);qu=pair(M,u,u);qv=pair(M,v,v);assert qu==F(rec['u_square'])and qv==F(rec['v_square'])and qv>0 and pair(M,u,v)==0
  assert pair(M,u,R[rec['anchor']])==1 and pair(M,v,R[rec['anchor']])==0 and q==norms[rec['anchor']]
  assert all(pair(M,u,R[i])==pair(M,v,R[i])==0 for i in rec['fixed_right'])
  branch=rec['radical_branch'];assert branch in(-1,1)
  rad=(q/qv,F(0),-q*q*qu/(4*qv));rad_i=ki(rad);assert rad_i[0]>0;root=scale(sqrt(rad_i),branch);p=scale(T,-q/2);coords=[add(scale(p,a),scale(root,b))for a,b in zip(u,v)];center=isum(scale(x,y)for x,y in zip(coords,mv(M,h)));assert center[1]<0
  finite={x['wall']:x['label']for x in rec['finite_contacts']};pairs=[]
  for i,w in enumerate(R):
   L=-q*pair(M,u,w)/2;B=branch*pair(M,v,w);value=add(scale(T,L),scale(sqrt(rad_i),B))
   if i in finite:angle_identity(L,B,rad,q,F(norms[i]),finite[i],value)
   else:dotted(value,q,norms[i])
   pairs.append(value)
  target=rec['target_integer_ray'];tq=pair(M,target,target);assert tq>0 and pair(M,target,h)<0 and max(pair(M,target,w)for w in R)<=0
  cut=isum(scale(x,y)for x,y in zip(coords,mv(M,target)));assert cut[0]>0
  newpairs=[]
  for other,n in zip(cap_intervals,cap_norms):
   co=[isum(scale(x,a)for x,a in zip(other,row))for row in M];value=isum(mul(a,b)for a,b in zip(coords,co));dotted(value,q,n);newpairs.append(value)
  cap_intervals.append(coords);cap_norms.append(q);records.append({'candidate':rec['candidate_index'],'norm':str(q),'radicand_coefficients':list(map(str,rad)),'center_pairing':list(map(str,center)),'finite_contacts':{str(i):m for i,m in finite.items()},'old_pair_intervals':[list(map(str,x))for x in pairs],'previous_cap_pair_intervals':[list(map(str,x))for x in newpairs],'target_norm':str(tq),'target_cut':list(map(str,cut))});print('verified',rec['candidate_index'],flush=True)
 out={'status':'PASS_EXACT_835_FACET_PI7_INCIDENCE_BRANCH','sources':{p:sha(p)for p in(RAW,BANK,REC,CORE)},'coordinate_convention':'last coordinate multiplies sqrt2; form diag(Q,2)','interval_denominator':str(S),'T_polynomial':[1,-2,-1,1],'T_interval':list(map(str,T)),'base_facets':812,'new_facets':len(records),'facets':812+len(records),'base_pair_counts':dict(counts),'new_caps':records,'scope':'Complete exact Coxeter-pair gate and strict center prove835 actual facets. Twenty-three prescribed positive-square parent rays are cut. No compactness, finite-volume, arithmeticity or literature-novelty conclusion is claimed.'}
 assert out['facets']==835
 return out,(M,R,h,cap_intervals,cap_norms)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out,_=run()
 if a.verify:assert out==read(OUT)
 else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
