#!/usr/bin/env python3
"""Numerical facet-complement diagnosis for the exact176-wall H20 chamber."""
import json,numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from fresh_h20_coupled_completion_20261013_second_exact import setup as setup175,OUT as SOURCE
from fresh_h20_coupled_completion_20261013 import num,read,ROOT
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension

def setup():
 Q,h,hc,roots,f=setup175();d=read(SOURCE);g=Extension(f,f.parse(d['field']['radicand']));r=tuple(map(g.parse,d['new_cap']['normal']));roots.append((g,r))
 Qn=np.array([[num(B,x)for x in rr]for rr in Q]);Rn=np.array([[num(ff,x)for x in rr]for ff,rr in roots]);hn=np.array([num(B,x)for x in h]);ids=[p['old']for p in d['new_cap']['old_pair_checks']if p['kind']!='dotted']
 return Qn,Rn,hn,ids,Q,h,hc,roots,f,g
if __name__=='__main__':
 Q,R,h,ids,*_=setup();N=null_space(R[ids+[175]]@Q);vals,V=np.linalg.eigh(N.T@Q@N);N=N@V@np.diag(1/np.sqrt(vals));P=R@Q@N;hp=h@Q@N
 print('Gram eigen',np.linalg.eigvalsh(R[ids+[175]]@Q@R[ids+[175]].T).tolist());print('orthogonal plane norms',N.T@Q@N,'center',hp)
 objectives=[np.array([1.,0.]),np.array([-1.,0.]),np.array([0.,1.]),np.array([0.,-1.])]
 for v in objectives:
  r=linprog(v,A_ub=P,b_ub=np.zeros(len(P)),bounds=[(-1,1)]*2,method='highs');print('cone',v.tolist(),r.success,None if not r.success else r.x.tolist())
 print('oldshort projections',[(i,np.round(P[i],8).tolist())for i in range(24)if np.linalg.norm(P[i])>1e-7])
 print('small projected covectors',sorted((float(np.linalg.norm(v)/np.sqrt(R[i]@Q@R[i])),i)for i,v in enumerate(P)if np.linalg.norm(v)>1e-7)[:20])
 (ROOT/'data/enumerations/fresh_h20_facet_bridge_20261013_diagnosis.json').write_text(json.dumps({'status':'NUMERICAL_FACET_COMPLEMENT_DIAGNOSIS','source':SOURCE,'finite_neighbors':ids,'complement_basis':N.tolist(),'wall_projections':P.tolist(),'center_projection':hp.tolist()},indent=2)+'\n')
