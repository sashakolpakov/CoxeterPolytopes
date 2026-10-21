import sys,json,numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from fresh_h22_new_cusp_completion_20261013_verify import setup,families,DEN,PRE,P,read
p,Q,W,h,meta,gens,D,V,B,A,ids,sizes=setup();cert=read('results/'+PRE+'_verify.json');seed=np.array(cert['base_old_numerators'])/DEN;norm=np.array([float(__import__('fractions').Fraction(x))for x in cert['new_cap_norms']]);shifts=np.array([3,3,4,4,4,4,4,10,5,.05,.05]);data=[]
for factor in(1,1.1,1.25,1.5,2):
 cs=seed+factor*shifts[:,None]*D[0]/6;worst=-1e100;count=0;where=None
 for index,(name,base,ix)in enumerate(families(D,V,B,A,ids)):
  if index==0:
   t=np.where(sizes>=4,3.,4.);n=np.einsum('ij,jk,ik->i',A,Q,A)
  else:t=np.full(len(ix),[10,5,.05,.05][index-1]);n=np.full(len(ix),norm[6+index])
  rows=base/DEN+factor*t[:,None]*D[ix]/6
  pair=(cs@Q)@rows.T+44*factor**2*shifts[:,None]*t
  measure=pair+np.sqrt(norm[:,None]*n[None,:]);measure[:,ix==0]=-np.inf
  maximum=float(np.max(measure));bad=int(np.sum(measure>1e-6));count+=bad
  if maximum>worst:worst=maximum;ij=np.unravel_index(np.argmax(measure),measure.shape);where=(name,int(ij[0]),int(ij[1]),int(ix[ij[1]]))
 print('factor',factor,'worst_dotted_excess',worst,'bad',count,'where',where,flush=True)
 data.append({'factor':factor,'worst_dotted_excess':worst,'bad':count,'where':where})
 if count==0:break
(P/('data/enumerations/'+PRE+'_tight_scout.json')).write_text(json.dumps(data,indent=2)+'\n')
