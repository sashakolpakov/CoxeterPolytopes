import sys,json,numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from fresh_h22_new_cusp_completion_20261013_verify import setup,families,DEN,PRE,P,read
from fractions import Fraction as F
p,Q,W,h,meta,gens,D,V,B,A,ids,sizes=setup();cert=read('results/'+PRE+'_verify.json');seed=np.array(cert['base_old_numerators'])/DEN;norm=np.array([float(F(x))for x in cert['new_cap_norms']]);shifts=np.array([3,3,4,4,1,1,2,10,5,.05,.05]);cs=seed+shifts[:,None]*D[0]/6;an=np.einsum('ij,ij->i',A,A@Q);labels=np.arange(2,101);cos=np.cos(np.pi/labels);records=[]
for index,(name,base,ix)in enumerate(families(D,V,B,A,ids)):
 if index==0:t=np.where(sizes>=4,3.,np.where(sizes==3,4.,np.where(an==2,1.,2.)));n=an
 else:t=np.full(len(ix),[10,5,.05,.05][index-1]);n=np.full(len(ix),norm[6+index])
 rows=base/DEN+t[:,None]*D[ix]/6;pair=(cs@Q)@rows.T+44*shifts[:,None]*t
 ratios=-pair/np.sqrt(norm[:,None]*n[None,:]);ratios[:,ix==0]=2
 bad=[];finite=[]
 for j in range(11):
  ii=np.flatnonzero(ratios[j]<1-1e-6);vals=ratios[j,ii];wrong=np.min(abs(vals[:,None]-cos),axis=1)>1e-5
  bad.extend((j,int(i),float(ratios[j,i]))for i in ii[wrong]);finite.append(len(ii))
 print(name,'finite',finite,'bad',len(bad),'first',bad[:10],flush=True)
 records.append({'family':name,'finite_counts':finite,'bad_count':len(bad),'first_bad':bad[:20]})
(P/('data/enumerations/'+PRE+'_lower_orbit_scout.json')).write_text(json.dumps(records,indent=2)+'\n')
