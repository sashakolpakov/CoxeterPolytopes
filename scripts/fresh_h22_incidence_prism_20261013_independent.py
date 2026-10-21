#!/usr/bin/env python3
"""Separate Fraction/polynomial audit of two representative pi/7 caps.

No producer or project arithmetic imports. Reconstructs the812 parent roots
from raw old walls and binary profiles; independently isolates T and square
roots by rational bisection, then checks all1624 cap-parent relations.
"""
import argparse, hashlib, json
from collections import Counter, deque
from fractions import Fraction as F
from pathlib import Path
P=Path(__file__).resolve().parents[1]
PRE='fresh_h22_incidence_prism_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
BANK='data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json'
REC='data/enumerations/'+PRE+'_recipes.json'
OUT='results/'+PRE+'_independent.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),F(0))
def mv(A,v):return[dot(r,v)for r in A]
def pair(A,a,b):return dot(a,mv(A,b))
def pa(a,b):return tuple(x+y for x,y in zip(a,b))
def ps(a,c):return tuple(x*c for x in a)
def pm(a,b):
 c=[F(0)]*5
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 for i in range(4,2,-1):
  z=c[i];c[i-1]+=z;c[i-2]+=2*z;c[i-3]-=z
 return tuple(c[:3])
ZERO=(F(0),)*3;ONE=(F(1),F(0),F(0));TT=(F(0),F(0),F(1))
def ia(a,b):return(a[0]+b[0],a[1]+b[1])
def im(a,b):
 c=[x*y for x in a for y in b];return min(c),max(c)
def isc(a,s):return im(a,(s,s))
def isum(v):
 out=(F(0),F(0))
 for a in v:out=ia(out,a)
 return out
def tau():
 lo,hi=F(9,5),F(19,10)
 def f(t):return t**3-t*t-2*t+1
 assert f(lo)<0<f(hi) and 3*lo*lo-2*hi-2>0
 for _ in range(230):
  mid=(lo+hi)/2
  if f(mid)<0:lo=mid
  else:hi=mid
 return lo,hi
def root(a):
 assert a[0]>0;lo=F(0);hi=F(1)
 while hi*hi<a[1]:hi*=2
 for _ in range(240):
  mid=(lo+hi)/2
  if mid*mid<a[0]:lo=mid
  elif mid*mid>a[1]:hi=mid
  else:break
 assert lo*lo<=a[0] and hi*hi>=a[1]
 return lo,hi
def run():
 raw=read(RAW)['coxeter_input'];Q=[[-F(x)for x in r]for r in raw['basis_gram']];W=[[F(x)for x in r]for r in raw['wall_coordinates_as_columns_transposed']]
 cov=[mv(Q,w)for w in W];G=[[dot(w,c)for c in cov]for w in W];short=[i for i in range(210)if G[i][i]==2]
 color={short[0]:0};todo=deque(color)
 while todo:
  i=todo.popleft()
  for j in short:
   if i!=j and G[i][j]==-1:
    if j not in color:color[j]=1-color[i];todo.append(j)
    assert color[j]!=color[i]
 parts=[[i for i in short if color[i]==c]for c in(0,1)];assert list(map(len,parts))==[21,21]
 sums=[[sum(W[i][j]for i in part)for j in range(22)]for part in parts]
 vp=[(5*a+2*b)/21 for a,b in zip(*sums)];h=[(a+b)/3 for a,b in zip(*sums)]+[F(0)]
 M=[r+[F(0)]for r in Q]+[[F(0)]*22+[F(2)]];R=[w+[F(0)]for w in W]
 C=[]
 for rec in read(BANK)['candidates']:
  if rec['twice_height_square']!=4:continue
  a=rec['A'];m=rec['mask'];ids=[parts[0][i]for i in range(21)if(m>>i)&1]
  v=[(a*vp[j]-sum(W[i][j]for i in ids))/2 for j in range(22)]
  assert [2*x for x in v]==rec['twice_old_projection'] and pair(Q,v,v)==0
  C.append(v)
 assert len(C)==301
 R+=[v+[F(1)]for v in C]+[v+[F(-1)]for v in C];RC=[mv(M,w)for w in R];norm=[dot(w,c)for w,c in zip(R,RC)];assert len(R)==812
 assert pair(M,h,h)==-14 and all(dot(h,c)<0 for c in RC)
 T=tau();T2=im(T,T)
 def pi(a):return isum(((a[0],a[0]),isc(T,a[1]),isc(T2,a[2])))
 outputs=[];new=[]
 for index in(0,22):
  rec=read(REC)['recipes'][index];u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm']);qu=pair(M,u,u);qv=pair(M,v,v)
  assert qv>0 and pair(M,u,v)==0 and q>0
  assert qu==F(rec['u_square']) and qv==F(rec['v_square'])
  rad=(q/qv,F(0),-q*q*qu/(4*qv));assert pa(ps(TT,q*q*qu/4),ps(rad,qv))==ps(ONE,q)
  branch=rec['radical_branch'];assert branch in(-1,1);ri=pi(rad);rt=root(ri)
  coords=[ia(isc(T,-q*a/2),isc(rt,branch*b))for a,b in zip(u,v)]
  cp=isum(isc(x,y)for x,y in zip(coords,mv(M,h)));assert cp[1]<0
  finite={x['wall']:x['label']for x in rec['finite_contacts']};counts=Counter();twice=0;zero_nontrivial=0
  for i,covw in enumerate(RC):
   L=-q*dot(u,covw)/2;B=branch*dot(v,covw);p=ia(isc(T,L),isc(rt,B))
   if i not in finite:assert p[1]<0 and p[1]*p[1]>q*norm[i];counts['dotted']+=1;continue
   m=finite[i];E=ZERO if m==2 else ps(TT,q*norm[i]/4)if m==7 else ps(ONE,q*norm[i]*{3:F(1,4),4:F(1,2),6:F(3,4)}[m])
   Y=pa(pa(E,ps(TT,-L*L)),ps(rad,-B*B))
   if not L*B:assert Y==ZERO
   else:
    assert pm(Y,Y)==ps(pm(TT,rad),4*L*L*B*B)
    yr=pi(Y);assert (yr[0]>0)if L*B>0 else(yr[1]<0)
    twice+=1;zero_nontrivial+=m==2
   if m==2:assert p[0]<=0<=p[1]
   else:assert p[1]<0
   counts['pi/'+str(m)]+=1
  assert dot(u,RC[rec['anchor']])==1 and dot(v,RC[rec['anchor']])==0
  assert all(dot(u,RC[i])==dot(v,RC[i])==0 for i in rec['fixed_right'])
  target=list(map(F,rec['target_integer_ray']));assert pair(M,target,target)>0 and pair(M,target,h)<0
  assert max(dot(target,c)for c in RC)<=0
  cut=isum(isc(a,b)for a,b in zip(coords,mv(M,target)));assert cut[0]>0
  outputs.append({'recipe_index':index,'candidate':rec['candidate_index'],'norm':str(q),'relations':dict(counts),'twice_squaring_nonzero_cross_terms':twice,'zero_angle_nonzero_cross_terms':zero_nontrivial,'strict_center':True,'strict_target_cut':True,'norm_identity':True})
  new.append((coords,q))
 cross=isum(im(x,isum(isc(y,z)for y,z in zip(new[1][0],row)))for x,row in zip(new[0][0],M))
 assert cross[1]<0 and cross[1]*cross[1]>new[0][1]*new[1][1]
 return {'status':'PASS_SEPARATE_FRACTION_TWO_REPRESENTATIVE_PI7_CAPS','sources':{p:sha(p)for p in(RAW,BANK,REC)},'arithmetic':'Independent Fraction polynomial quotient and rational bisection; no producer/project imports. T is the unique root in(1.8,1.9) of T^3-T^2-2T+1.','records':outputs,'new_new_pair':'strictly dotted','scope':'Reconstructs812 parent normals; verifies two caps,1624 new-parent relations,mutual dotted relation,exact norm identities,centers,target feasibility/cuts. It does not independently replay the other21 caps or the whole835 geometry.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
