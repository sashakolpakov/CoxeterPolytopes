import json,numpy as np,sympy as sp
from fractions import Fraction
from pathlib import Path
P=Path(__file__).resolve().parents[1];d=json.loads((P/'results/fresh_h22_gap_face_family_20261013.json').read_text());raw=json.loads((P/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json').read_text())['coxeter_input'];Q=np.array(raw['basis_gram'],dtype=np.int64)*-1;W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64);h=W.sum(0);h//=np.gcd.reduce(h)
seed=np.array(d['six_old_projection']);partner=np.array(d['new_cusp']['six_partner']);delta=(seed+partner)/6
ids=d['new_cusp']['common_old_walls'];G=W[ids]@Q@W[ids].T
unseen=set(range(len(ids)));cs=[]
while unseen:
 c={unseen.pop()};todo=list(c)
 for i in todo:
  for j in list(unseen):
   if G[i,j]:c.add(j);todo.append(j);unseen.remove(j)
 cs.append([ids[i]for i in sorted(c)])
ns=sp.Matrix(np.vstack([W[ids]@Q,seed@Q,partner@Q]).tolist()).nullspace();SM=sp.Matrix(Q.tolist());v=ns[0];w=ns[1]-v*((v.T*SM*ns[1])[0]/(v.T*SM*v)[0]);vs=[-sp.Matrix(W[c].sum(0))for c in cs]+[v,-v,w,-w]
U=np.load(P/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz')['twice_old_projections']/2
z=np.load(P/'data/enumerations/fresh_h22_cusp_orbits_20261013.npz');D=z['old_null_projections'];ht=z['four_heights']/4
R=np.load(P/'data/enumerations/fresh_h22_gap_face_orbit_20261013.npz')['six_old_projections']/6
families=[('old',W,np.zeros(len(W)),np.diag(W@Q@W.T)),('U',U,np.full(len(U),np.sqrt(2.5)),np.ones(len(U))),('D',D,ht,ht**2),('R',R,np.full(len(R),np.sqrt(11)),np.full(len(R),2/3)),('center',h[None,:],np.array([4.]),np.array([2.]))]
x=np.array([float(Fraction(s)) for s in d['uniform_strict_escape']['old_coordinates']]);xz=2*np.sqrt(11)
out=[]
for i,sv in enumerate(vs):
 a=np.array(sv,dtype=float).ravel();nv=float((sv.T*SM*sv)[0]);bound=-1.;which=None
 for name,old,z,norm in families:
  dp=old@(Q@delta)+2*np.sqrt(11)*z;vp=old@(Q@a);mask=dp<-1e-7
  ratio=(vp[mask]+np.sqrt(nv*norm[mask]))/(-dp[mask]);m=float(max(ratio))
  if m>bound:bound=m;which=name
 t=int(np.ceil(bound))+1
 cap=a+t*delta;pair=cap@Q@x+2*t*np.sqrt(11)*xz
 row={'i':i,'component':cs[i]if i<7 else 'tangent','base':list(map(str,sv)),'norm':str((sv.T*SM*sv)[0]),'numeric_min_shift':bound,'binding_family':which,'shift':t,'old_escape_pair':float(pair)}
 out.append(row);print({k:v for k,v in row.items()if k!='base'},flush=True)
(P/'data/enumerations/fresh_h22_new_cusp_completion_20261013_scout.json').write_text(json.dumps({'caps':out,'components':cs,'scope':'Numerical sufficient shifts only; exact replay needed.'},indent=2)+'\n')
