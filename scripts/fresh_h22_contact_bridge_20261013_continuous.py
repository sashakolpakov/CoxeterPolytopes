#!/usr/bin/env python3
"""Numerical continuous acuteness and all-finite-anchor interval diagnostic."""
import json,math
from collections import Counter
import numpy as np
import fresh_h22_contact_bridge_20261013_scout as S
TAG=S.TAG+'_continuous';OUT=S.V.P/'results'/f'{TAG}.json'
def add_bound(lo,hi,a,b,tag,active):
 # a*T-b <= 0.
 if abs(a)<1e-9:
  if b< -1e-8:return lo,-math.inf,tag,tag
  return lo,hi,*active
 t=b/a;l,u=active
 if a<0 and t>lo:lo=t;l=tag
 if a>0 and t<hi:hi=t;u=tag
 return lo,hi,l,u
def finite_interval(a,b,n):
 # Existing acuteness gives a*T-b<=0; a finite or parallel pair also
 # satisfies (a*T-b)^2 <= 4*n*T.
 aa=a*a;bb=-2*a*b-4*n;cc=b*b
 if aa<1e-18:
  if abs(bb)<1e-12:return (0.,math.inf)if cc<=1e-9 else None
  return(-cc/bb,math.inf)if bb<0 else(0.,-cc/bb)
 D=bb*bb-4*aa*cc
 if D< -1e-8:return None
 D=max(D,0.);large=(-bb+math.sqrt(D))/(2*aa)
 small=cc/(aa*large)if large else(-bb-math.sqrt(D))/(2*aa)
 return small,large
def values(a,b,n,T):return(-a*math.sqrt(T)+b/math.sqrt(T))/(2*np.sqrt(n))
def ranges(a,b,n,lo,hi):
 c0=values(a,b,n,lo)
 if math.isinf(hi):c1=np.where(a< -1e-10,math.inf,np.where(a>1e-10,-math.inf,0.))
 else:c1=values(a,b,n,hi)
 mn=np.minimum(c0,c1);mx=np.maximum(c0,c1);t=np.divide(-b,a,out=np.full_like(a,-1.),where=abs(a)>1e-10);ok=(t>lo)&(t<hi)
 if np.any(ok):c=values(a[ok],b[ok],n[ok],t[ok])if False else(-a[ok]*np.sqrt(t[ok])+b[ok]/np.sqrt(t[ok]))/(2*np.sqrt(n[ok]));mn[ok]=np.minimum(mn[ok],c);mx[ok]=np.maximum(mx[ok],c)
 return mn,mx
def run():
 M,W,h,G,n,meta,doc,vertices=S.load();WM=W@M;hc=M@h;rows=[];counts=Counter();gaps=[(m,math.cos(math.pi/m),math.cos(math.pi/(m+1)))for m in range(2,200)]
 for vi,src in enumerate(doc['compact_H2_vertices']):
  inc=src['incident'];AI=np.linalg.inv(G[np.ix_(inc,inc)]);v=vertices[vi]/math.sqrt(-src['norm'])
  for k,j in enumerate(inc):
   edge=[i for i in inc if i!=j];e=(AI[k]@W[inc])/math.sqrt(AI[k,k]);a=WM@(v-e);b=WM@(v+e);a[edge]=0.;b[edge]=0.;lo,hi=1.,math.inf;l,u='outgoing',None
   for i in range(860):lo,hi,l,u=add_bound(lo,hi,a[i],b[i],i,(l,u))
   lo,hi,l,u=add_bound(lo,hi,(v-e)@hc,(v+e)@hc,'center',(l,u))
   lo,hi,l,u=add_bound(lo,hi,-(v-e)[-1],-(v+e)[-1],'upper_height',(l,u))
   row={'vertex':vi,'released':j,'edge':edge,'T_interval':[lo,hi],'active_lower':l,'active_upper':u,'anchor_intervals':[]};rows.append(row)
   if not lo<hi:counts['empty_acute_edges']+=1;continue
   counts['nonempty_acute_edges']+=1
   for anchor in range(210,511):
    if anchor in edge:continue
    ii=finite_interval(a[anchor],b[anchor],n[anchor]);counts['anchor_tests']+=1
    if ii is None:counts['finite_anchor_no_real']+=1;continue
    al,ah=max(lo,ii[0]),min(hi,ii[1])
    if not al<ah-1e-8:counts['finite_anchor_empty']+=1;continue
    mn,mx=ranges(a,b,n,al,ah);bad=[]
    for m,c0,c1 in gaps:
     inds=np.where((mn>c0+1e-7)&(mx<c1-1e-7))[0]
     for z in inds:bad.append({'wall':int(z),'gap_labels':[m,m+1],'cosine_range':[float(mn[z]),float(mx[z])]})
    rr={'anchor':anchor,'T_interval':[al,ah],'forbidden_range_witness':bad[0]if bad else None};row['anchor_intervals'].append(rr);counts['nonempty_finite_anchor_intervals']+=1
    counts['whole_interval_forbidden_angle'if bad else'remaining_intervals']+=1
   print('edge',vi,j,'acute',lo,hi,'finiteanchors',len(row['anchor_intervals']),'remaining',sum(r['forbidden_range_witness']is None for r in row['anchor_intervals']),flush=True)
 out={'status':'NUMERICAL_CONTINUOUS_ALL_RIGHT_EDGE_DIAGNOSTIC','sources':{p:S.V.sha(p)for p in(S.P.OUT,'scripts/'+S.TAG+'_scout.py','scripts/'+TAG+'.py')},'counts':dict(counts),'rows':rows,'parameter':'r=(t(v-e)-t^-1(v+e))/sqrt2, T=t^2>1, v^2=-1,e^2=1,v.e=0','scope':'Numerical bounds on every all-right edge family at four retained vertices. Acuteness, upper height and common-center signs are continuous constraints. Finite-anchor intervals include every finite angle and the parallel limit, independently of anchor labels. A whole-interval forbidden-angle witness is only numerical until exactified. No assertion covers nonright edge profiles or wall replacement.'}
 OUT.write_text(json.dumps(out,indent=2)+'\n');print(dict(counts),OUT,flush=True)
if __name__=='__main__':run()
