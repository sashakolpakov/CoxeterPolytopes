#!/usr/bin/env python3
"""Separate exact interval check of the same missed cusp in the pi7 branch.

Imports only the separate Fraction audit arithmetic, never the producer.
The reflected mates of two caps have identical old projection and therefore
the same strict cusp pairing; their insertion cannot change this diagnosis.
"""
import argparse,json,hashlib
from fractions import Fraction as F
from pathlib import Path
from fresh_h22_incidence_prism_20261013_independent import tau,root,ia,isc,im,isum,pair,mv,dot
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_incidence_prism_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';REC='data/enumerations/'+PRE+'_recipes.json';CUSP='data/enumerations/fresh_h22_incidence_mixed_height_20261013_covering.json';BASE='data/enumerations/fresh_h22_incidence_mixed_height_20261013.json';OUT='results/'+PRE+'_uncapped.json'
BANK='data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def run():
 raw=read(RAW)['coxeter_input'];Q=[[-F(x)for x in r]for r in raw['basis_gram']];M=[r+[F(0)]for r in Q]+[[F(0)]*22+[F(2)]]
 c=read(CUSP);delta=list(map(F,c['old_null_ray']));b=read(BASE);h=list(map(F,b['center']));assert pair(Q,delta,delta)==0 and pair(Q,delta,h)==-8
 W=[[F(x)for x in r]for r in raw['wall_coordinates_as_columns_transposed']];assert max(pair(Q,delta,w)for w in W)==0
 D=[[F(x)for x in r]for r in b['null_projections']];assert max(pair(Q,delta,w)for w in D)==-2
 Dset=set(map(tuple,D));parents=[tuple(F(x,2)for x in r['twice_old_projection'])for r in read(BANK)['candidates']if r['twice_height_square']==4]
 assert len(parents)==len(set(parents))==301 and set(parents)<=Dset
 T=tau();T2=im(T,T)
 def pi(a):return isum(((a[0],a[0]),isc(T,a[1]),isc(T2,a[2])))
 records=[];coords=[]
 for i,rec in enumerate(read(REC)['recipes']):
  u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm']);qu=pair(M,u,u);qv=pair(M,v,v)
  assert qv>0 and pair(M,u,v)==0;rad=(q/qv,F(0),-q*q*qu/(4*qv));rt=root(pi(rad));branch=rec['radical_branch'];assert branch in(-1,1)
  r=[ia(isc(T,-q*a/2),isc(rt,branch*z))for a,z in zip(u,v)]
  p=isum(isc(x,y)for x,y in zip(r,mv(M,delta+[F(0)])));assert p[1]<0
  records.append({'recipe_index':i,'candidate':rec['candidate_index'],'cusp_pair_interval':list(map(str,p))});coords.append(r)
 # A rational strict escape for all835, including both reflections of any
 # saved cap: using abs(height) handles both signs uniformly.
 for k in range(1,7):
  n=100**k;z=4*10**k;x=[n*a+c for a,c in zip(delta,h)];point=x+[F(z)];cp=mv(M,point)
  oldok=max(pair(Q,x,w)for w in W)<0
  nullok=max(pair(Q,x,w)for w in D)+2*z<0
  vals=[]
  for r in coords:
   po=isum(isc(r[j],cp[j])for j in range(22));height=max(abs(r[22][0]),abs(r[22][1]));vals.append(po[1]+2*z*height)
  if oldok and nullok and max(vals)<0:break
 else:raise AssertionError('No strict witness in the bounded exact recipe')
 square=pair(M,point,point);assert square==16*n-14>0 and pair(M,point,h+[F(0)])<0
 return {'status':'PASS_EXACT_UNCAPPED_OLD_CUSP_IN_835_AND_REFLECTED_MATE_EXTENSIONS','sources':{p:sha(p)for p in(RAW,REC,CUSP,BASE,BANK)},'literal_parent_null_entries_covered':301,'old_null_ray':list(map(int,delta)),'strict_pi7_old_projection_pairs':records,'strict_escape':{'old_coordinates':list(map(str,x)),'last_coordinate_in_sqrt2_units':z,'square':str(square),'center_pair':str(pair(M,point,h+[F(0)])),'max_cap_pair_upper_bound':str(max(vals))},'scope':'Checks the23 cap recipes and all322 null projections (a superset of the812 parent301, explicitly checked against the literal source bank). It also allows reflection of every saved cap in the added-coordinate plane. Thus the same cusp is untouched in the835 branch and its837 two-reflected-mate extension. This does not inspect unrelated later caps.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
