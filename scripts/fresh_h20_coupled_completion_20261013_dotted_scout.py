exec(open('scripts/fresh_h20_coupled_completion_20261013_next_polar_scout.py').read().split('out=[];best=[]')[0])
out=[];best=[]
for ai,a in enumerate(allowed[:-1]):
 for li,L in enumerate(allowed):
  u=(2*a+L)*p+a*S;v=C-p
  aa=v@Q@v;bb=2*v@Q@u;cc=u@Q@u-2;disc=bb*bb-4*aa*cc
  if disc<-1e-8:continue
  for b in ((-bb+np.sqrt(max(0,disc)))/(2*aa),(-bb-np.sqrt(max(0,disc)))/(2*aa)):
   r=u+b*v;k=2*a+L-b
   if b*C0<2-1e-8 or r@Q@h>=0 or r@Q@target<=0:continue
   rat=-R@Q@r/np.sqrt(2*wn);g=(rat>=1-1e-7)|(np.min(abs(rat[:,None]-allowed[None,:]/2),axis=1)<1e-7)
   bad=np.flatnonzero(~g).tolist();rec={'a_label':ai,'leaf_label':li,'a':a,'b':b,'k':k,'bad':bad,'bad_ratios':[float(rat[i])for i in bad],'normal':r.tolist(),'cut':float(r@Q@target)};best.append(rec)
   if not bad:out.append(rec)
best.sort(key=lambda x:len(x['bad']));print('good',len(out),'best',[(x['a_label'],x['leaf_label'],x['b'],x['k'],x['bad'],x['bad_ratios'])for x in best[:16]])
json.dump({'good':out,'best':best[:30]},open('data/enumerations/fresh_h20_coupled_completion_20261013_dotted_scout.json','w'),indent=2)
