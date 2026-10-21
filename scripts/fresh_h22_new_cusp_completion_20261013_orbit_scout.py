#!/usr/bin/env python3
"""Numerical guide to a complete orbit of new-cusp link repairs."""
import json,numpy as np
from pathlib import Path
from fractions import Fraction as F
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_new_cusp_completion_20261013'
def read(n):return json.loads((P/n).read_text())
def run():
 d=read('results/fresh_h22_gap_face_family_20261013.json');s=read('data/enumerations/'+PRE+'_scout.json');raw=read('data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json')['coxeter_input'];Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
 z=np.load(P/('data/enumerations/'+PRE+'_frames.npz'));D=z['six_cusp_old'];V=z['first_tangent_numerator'];B=z['second_tangent_numerator'];N=len(D);den=5208
 G=W@Q@W.T;adj=[sum(1<<int(j) for j in np.flatnonzero(row))for row in G]
 # Seven A-type spherical components at every transported cusp.
 theta=[];cusp=[];sizes=[]
 for start in range(0,N,2048):
  pairs=D[start:start+2048]@Q@W.T
  for off,row in enumerate(pairs):
   ids=np.flatnonzero(row==0).tolist();assert len(ids)==18
   unseen=sum(1<<j for j in ids);local=[]
   while unseen:
    bit=unseen&-unseen;unseen^=bit;todo=bit;component=0
    while todo:
     bit=todo&-todo;todo^=bit;j=bit.bit_length()-1;component|=bit
     added=unseen&adj[j];unseen^=added;todo|=added
    js=[]
    while component:
     bit=component&-component;component^=bit;js.append(bit.bit_length()-1)
    local.append(js);theta.append(-W[js].sum(0));cusp.append(start+off);sizes.append(len(js))
   assert sorted(map(len,local))==[1,1,1,3,3,4,5]
  if start%20480==0:print('components',start,len(theta),flush=True)
 A=np.array(theta,dtype=np.int64);idx=np.array(cusp,dtype=np.int64)
 np.savez_compressed(P/('data/enumerations/'+PRE+'_highest_roots.npz'),old_bases=A,cusp_indices=idx,component_sizes=np.array(sizes,dtype=np.int8))
 families=[('theta',A*den,idx),('v+',V*(den//28),np.arange(N)),('v-',-V*(den//28),np.arange(N)),('w+',B,np.arange(N)),('w-',-B,np.arange(N))]
 seeds=np.array([[int(F(x)*den)for x in r['base']]for r in s['caps']],dtype=np.int64)
 records=[]
 for shift in(10,12,16,24,32,48,64,96,128):
  worst=-1e100;loc=None
  seedcaps=seeds/den+shift*D[0]/6
  for name,base,ids in families:
   roots=base/den+shift*D[ids]/6
   pairs=(seedcaps@Q)@roots.T+44*shift**2
   pairs[:,ids==0]=-np.inf
   ij=np.unravel_index(np.argmax(pairs),pairs.shape);m=float(pairs[ij])
   if m>worst:worst=m;loc=[name,int(ij[0]),int(ij[1]),int(ids[ij[1]])]
  records.append({'shift':shift,'maximum_nonlocal_pair':worst,'where':loc});print('shift',records[-1],flush=True)
  if worst<-20:break
 out={'status':'NUMERICAL_FULL_CUSP_ORBIT_REPAIR_SHIFT_SCOUT','positive_cusps':N,'new_positive_caps':11*N,'signed_new_caps':22*N,'records':records,'scope':'Numerical maximum pairs; exact verification required.'}
 (P/('data/enumerations/'+PRE+'_orbit_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':run()
