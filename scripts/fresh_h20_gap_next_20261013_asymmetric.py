import sys,json,numpy as np,itertools
sys.path[:0]=['scripts','src']
from fresh_h20_gap_next_20261013 import setup,read
from fresh_h20_edge_family_completion_20261013_bank import gate
from scipy.optimize import linprog
Q,R,h,w=setup();d1=read('data/enumerations/fresh_h20_gap_next_20261013_coupled.json');d2=read('data/enumerations/fresh_h20_gap_next_20261013_second_coupled.json');P=np.array(d2['polar'])/np.sqrt(2);U=np.vstack([R,np.array(d1['new_normal'])]);A=np.array([R[5],R[14],R[36]/3]);target=np.array(d2['escape']['vector']);ns=np.einsum('ij,jk,ik->i',U,Q,U)
labels=[3,4,5,6,8,10,12,0];values=[2*np.cos(np.pi/m)if m else 2 for m in labels];caps=[]
for p in itertools.product(range(len(labels)),repeat=3):
 co=np.array([values[i]for i in p]);s=sum(co);D=s*s-2*sum(co*co)+2
 if D<0:continue
 k=s+np.sqrt(D);r=k*P+co@A;ratios=-(U@Q@r)/np.sqrt(2*ns)
 if not np.all(gate(ratios))or r@Q@target<=1e-7:continue
 caps.append({'labels':[labels[i]for i in p],'coefficients':co.tolist(),'kappa':k,'normal':r.tolist(),'target_pairing':float(r@Q@target),'center_pairing':float(r@Q@h),'old97_pairing':float(r@Q@np.array(d2['new_normal']))})
print('candidates',len(caps),flush=True)
for i,c in sorted(enumerate(caps),key=lambda z:-z[1]['target_pairing'])[:8]:print(i,c['labels'],c['target_pairing'],c['old97_pairing'],flush=True)
for i,c in enumerate(caps):
 if c['labels'] in [[3,3,4],[3,3,6]]:print('SIMPLE',i,c,flush=True)
json.dump({'status':'NUMERICAL_ALTERNATIVE_ASYMMETRIC_LOCAL_METRIC_CAPS','source':'data/enumerations/fresh_h20_gap_next_20261013_second_coupled.json','retained_count':97,'target':target.tolist(),'labels':labels,'caps':caps},open('data/enumerations/fresh_h20_gap_next_20261013_asymmetric.json','w'),indent=2)
