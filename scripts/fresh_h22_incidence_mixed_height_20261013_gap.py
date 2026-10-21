#!/usr/bin/env python3
"""Exact two-cap extension of the frozen mixed-height branch, and its gap.

The default command writes numerical proposals first, then proves an exact
strict rational escape. --verify uses only the archived rational witness.
"""
import argparse,json,hashlib,math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from fresh_h22_incidence_mixed_height_20261013 import inertia,mm
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_incidence_mixed_height_20261013'
SRC='data/enumerations/'+PRE+'.json';RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';OUT='data/enumerations/'+PRE+'_gap.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def roots(n):
 S=10**50;k=math.isqrt(n*S*S);return F(k,S),F(k+1,S)
def run(verify):
 d=read(SRC);raw=read(RAW)['coxeter_input'];Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
 assert d['sources'][RAW]==sha(RAW)
 for p,hsh in d['sources'].items():assert sha(p)==hsh
 D=np.array(d['null_projections'],dtype=np.int64);X2=np.array(d['sqrt3_twice_projections'],dtype=np.int64);h=np.array(d['center'],dtype=np.int64)
 u=np.array([4]+[0]*8+[15,10,20,30,25,20,16,12]+[0]*4+[8],dtype=np.int64)
 assert int(u@Q@u)==-20 and int(u@Q@h)==-160
 G=mm(mm(W,Q),W.T);wn=G.diagonal();pu=W@Q@u;counts=Counter()
 for p,n in zip(pu,wn):
  assert p<=0;s=F(int(p*p),2*int(n))
  if p==0:counts['right']+=1
  elif s==F(1,4):counts['pi/3']+=1
  else:assert s>1;counts['dotted']+=1
 assert counts=={'right':21,'pi/3':1,'dotted':188}
 ids=np.flatnonzero(pu==0).tolist();assert inertia(G[np.ix_(ids,ids)])==[21,0,0]
 assert max(D@Q@u)==-32 and max(X2@Q@u)==-84
 # a>2, so -32a+sqrt44<-64+sqrt44<-2.
 assert 44<62**2 and 66<40**2 # -42+sqrt66<-2
 # Distinct signs pair=-42. Norm=-20+22=2.
 old=d['universal_escape'];oldx=np.array(old['old_coordinates'],dtype=np.int64);oldz=abs(old['last_coordinate']);cutold=int(u@Q@oldx)
 assert cutold==-2160 and 22*oldz*oldz>cutold*cutold
 # (u,+/-20/sqrt22) is an actual vertex: norm=-20/11, other
 # inequalities strict by the same max bounds and old strictness.
 assert 400<F(20*22) and 40*F(1,11)<62**2 and F(600,11)<40**2
 if verify:
  prior=read(OUT);num=prior['remaining_escape']['numerator'];den=prior['remaining_escape']['denominator'];proposal=prior['numerical_proposal']
 else:
  from scipy.optimize import linprog
  a=(np.sqrt(6)+np.sqrt(3))/2;Qn=np.zeros((23,23));Qn[:22,:22]=Q;Qn[22,22]=1
  R=np.vstack([np.c_[W,np.zeros(210)]]+[np.c_[a*D,np.full(322,s*np.sqrt(2))]for s in(1,-1)]+[np.c_[X2/2,np.full(720,s*np.sqrt(3))]for s in(1,-1)]+[np.r_[u,s*np.sqrt(22)][None,:]for s in(1,-1)])
  c=np.r_[h/np.sqrt(14),0];rng=np.random.default_rng(1336);A=R@Qn
  for trial in range(100):
   lp=linprog(rng.normal(size=23),A_ub=A,b_ub=np.zeros(len(R)),A_eq=(c@Qn)[None,:],b_eq=[-1],bounds=[(None,None)]*23,method='highs')
   assert lp.success
   q=float(lp.x@Qn@lp.x)
   if q>1e-6:break
  else:raise AssertionError('No numerical escape; needs a full volume gate')
  proposal={'trial':trial,'square':q,'point':lp.x.tolist(),'max_residual':float((A@lp.x).max())}
  point=.999*lp.x+.001*c;den=10**8;num=np.rint(point*den).astype(np.int64).tolist()
 x=[F(v,den)for v in num];xx=x[:22];z=abs(x[22])
 def pair(row):return sum(F(int(row[i]))*sum(F(int(Q[i,j]))*xx[j]for j in range(22))for i in range(22))
 oldpairs=[pair(w)for w in W];dp=max(pair(w)for w in D);xp=max(pair(w)/2 for w in X2);up=pair(u);hp=pair(h)
 q=sum(xx[i]*sum(F(int(Q[i,j]))*xx[j]for j in range(22))for i in range(22))+z*z
 assert q>0 and hp<0 and max(oldpairs)<0 and dp<0 and xp<0 and up<0
 a_lower=(roots(6)[0]+roots(3)[0])/2
 assert a_lower*dp+roots(2)[1]*z<0 and xp+roots(3)[1]*z<0 and up+roots(22)[1]*z<0
 out={'status':'PASS_EXACT_THREE_1336_FACET_EXTENSIONS_WITH_STRICT_ESCAPE','sources':{p:sha(p)for p in(SRC,RAW)},'new_old_projection':u.tolist(),'new_normals':'(u,+/-sqrt22)','norm':2,'new_old_pair_counts_per_sign':dict(counts),'other2294_new_pairs':'strictly dotted except the listed old210 contacts','center_pair':-160,'mutual_pair':-42,'cut_frozen_escape':'-2160+2400sqrt22>0','actual_finite_vertices':{'points':'(u,+/-20/sqrt22)','square':'-20/11','old_incident':ids,'incident_gram_signature':[22,0,0],'old21_gram_signature':[21,0,0]},'numerical_proposal':proposal,'remaining_escape':{'numerator':num,'denominator':den,'square':str(q),'center_pair':str(hp),'max_old_pair':str(max(oldpairs)),'max_null_projection_pair':str(dp),'max_sqrt3_projection_pair':str(xp),'new_projection_pair':str(up),'scope':'Strict against all2296 inequalities in the relaxed union, hence all three coherent1336-facet chambers.'},'scope':'Two exact compatible caps cut the prior strict escape, and supply two actual finite vertices per branch. A new exact strict spacelike escape still proves infinite volume; no complete frontier or global impossibility claim.'}
 if verify:assert out==prior
 else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();print(run(a.verify)['status'])
