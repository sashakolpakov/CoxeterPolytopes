#!/usr/bin/env python3
"""Produce rational two-plane recipes for selected pi/7 truncations."""
import json,math
from functools import reduce
from pathlib import Path
import numpy as np
import sympy as s
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_incidence_prism_20261013'
def primitive(v):
 den=s.ilcm(*[x.q for x in v]);v=[int(x*den)for x in v];d=reduce(math.gcd,v);return s.Matrix([x//d for x in v])
def run():
 raw=json.loads((P/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json').read_text())['coxeter_input'];Q=-s.Matrix([[int(x)for x in row]for row in raw['basis_gram']]);Q2=s.diag(Q,2)
 W=[[int(x)for x in row]+[0]for row in raw['wall_coordinates_as_columns_transposed']];bank=json.loads((P/'data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json').read_text())['candidates'];C=[[s.Rational(x,2)for x in c['twice_old_projection']]for c in bank if c['twice_height_square']==4];R=s.Matrix(W+[r+[1]for r in C]+[r+[-1]for r in C]);assert all(x.q==1 for x in R)
 old=json.loads((P/('data/enumerations/'+PRE+'_planes.json')).read_text())['candidates'];selected=json.loads((P/('data/enumerations/'+PRE+'_selected.json')).read_text())['candidate_indices'];vertices={r['trial']:r for r in json.loads((P/('data/enumerations/'+PRE+'_scout.json')).read_text())['runs']};out=[]
 h=s.Matrix(json.loads((P/'results/fresh_h22_incidence_core_20261013.json').read_text())['center']+[0])
 for ci in selected:
  c=old[ci];fixed=c['fixed_right'];anchor=c['anchor'];M=R[fixed,:]*Q2;E=M.nullspace();assert len(E)==2
  w=R[anchor,:].T;aa=[(e.T*Q2*w)[0]for e in E];k=next(i for i,a in enumerate(aa)if a);u=E[k]/aa[k]
  v=primitive((R[fixed+[anchor],:]*Q2).nullspace()[0]);qv=(v.T*Q2*v)[0];assert qv>0;u-=((u.T*Q2*v)[0]/qv)*v;qu=(u.T*Q2*u)[0];q=(w.T*Q2*w)[0]
  wanted=np.array(c['physical_normal']);wanted[-1]/=np.sqrt(2);wanted*=np.sqrt(float(q)/2);scalar=float((s.Matrix(wanted).T*Q2*v)[0]/qv);branch=1 if scalar>0 else-1
  target=primitive((R[vertices[c['vertex_trial']]['active'],:]*Q2).nullspace()[0]);target=target if(target.T*Q2*h)[0]<0 else-target
  out.append({'candidate_index':ci,'fixed_right':fixed,'anchor':anchor,'label':7,'norm':str(q),'u':list(map(str,u)),'v':list(map(str,v)),'u_square':str(qu),'v_square':str(qv),'radical_branch':branch,'finite_contacts':c['finite_contacts'],'target_vertex_trial':c['vertex_trial'],'target_integer_ray':list(map(int,target))})
  print(ci,'q',q,'qu',qu,'qv',qv,'branch',branch,flush=True)
 (P/('data/enumerations/'+PRE+'_recipes.json')).write_text(json.dumps({'status':'RATIONAL_RECIPE_PROPOSALS','coordinate_convention':'last coefficient of sqrt2; form diag(Q,2)','recipes':out},indent=2)+'\n')
if __name__=='__main__':run()
