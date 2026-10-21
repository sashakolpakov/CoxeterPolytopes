import sys,numpy as np,json
sys.path.insert(0,'scripts')
from fresh_h20_coupled_completion_20261013 import *
Qe,he,hc,roots=setup();Q=np.array([[num(B,x)for x in r]for r in Qe]);h=np.array([num(B,x)for x in he]);R=np.array([[num(f,x)for x in r]for f,r in roots]);p=np.array([num(B,x)for x in vec(read(SOURCE)['final_hyperideal_point']['vector'])])*15/4
D,theta,x,y=frames(Qe,he,hc,roots);de=np.array([num(F10,t)for t in D]);dirs=[tuple(-x for x in t)for t in theta]+[x,tuple(-t for t in x),y,tuple(-t for t in y)];R=np.vstack([R,*[np.array([num(B,t)for t in v])+5*de for v in dirs]])
wn=np.einsum('ij,jk,ik->i',R,Q,R);S=(R[26]+R[28]+R[42])/3;C=R[49]*np.sqrt(6)/3;C0=1+2*np.sqrt(6)/3;target=np.array([num(B,t)for t in vec(read(SOURCE)['strict_escape']['vector'])]);qs=p@Q@p;hc0=p@Q@C
allowed=np.array([0]+[2*np.cos(np.pi/m)for m in (3,4,5,6,8,10,12)]+[2])
out=[];best=[]
for ai,a in enumerate(allowed):
 for bi,v in enumerate(allowed):
  b=v/C0;r0=a*S+b*C;bb=2*p@Q@r0;cc=r0@Q@r0-2;disc=bb*bb-4*qs*cc
  if disc<-1e-9:continue
  for k in ((-bb+np.sqrt(max(0,disc)))/(2*qs),(-bb-np.sqrt(max(0,disc)))/(2*qs)):
   r=k*p+r0
   if r@Q@h>=0 or r@Q@target<=0:continue
   rat=-R@Q@r/np.sqrt(2*wn);g=(rat>=1-1e-7)|(np.min(abs(rat[:,None]-allowed[None,:]/2),axis=1)<1e-7)
   bad=np.flatnonzero(~g).tolist();rec={'a_label':ai,'b_label':bi,'a':a,'b':b,'k':k,'bad':bad,'bad_ratios':[float(rat[i])for i in bad],'normal':r.tolist(),'cut':float(r@Q@target)};best.append(rec)
   if not bad:out.append(rec)
best.sort(key=lambda x:len(x['bad']));print('good',len(out),'best',[(x['a_label'],x['b_label'],x['k'],x['bad'],x['bad_ratios'])for x in best[:12]])
json.dump({'good':out,'best':best[:30]},open('data/enumerations/fresh_h20_coupled_completion_20261013_next_polar_scout.json','w'),indent=2)
