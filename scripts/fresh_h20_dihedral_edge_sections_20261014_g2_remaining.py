#!/usr/bin/env python3
"""Exact boundary audit of every numerically unresolved G2 closure recipe.

Only the discrete recipe is imported. All edge-plane quadratics and signs
are rebuilt with Fraction arithmetic and an exact square-root ring.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import fresh_h20_cusp_sections_20261013 as P
import fresh_h22_multiquadratic_20261013_review as R
ROOT=P.ROOT; TAG='fresh_h20_dihedral_edge_sections_20261014_g2'
SUMMARY=f'results/{TAG}_summary.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def run():
 d=json.loads((ROOT/SUMMARY).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 assert sha(f'results/{TAG}_scout.json.gz')==d['archive_sha256']
 raw=json.loads((ROOT/P.RAW).read_text())['coxeter_input']
 Q=[[-int(x)for x in r]for r in raw['basis_gram']]
 W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 C=[P.mv(Q,w)for w in W];G=[[P.dot(w,c)for c in C]for w in W]
 h=[sum(w[j]for w in W)for j in range(22)];g=math.gcd(*h);h=[F(x,g)for x in h];hc=P.mv(Q,h)
 walk=json.loads((ROOT/P.WALK).read_text())['corner_orbits']; rows=[]
 assert len(d['unresolved_cases'])==66
 for r in d['unresolved_cases']:
  fam=d['families'][r['family']];corner=walk[fam['corner']];edge=corner['edge_orbits'][fam['edge_index']]
  E=fam['edge_walls'];assert E==edge['edge_walls'];a=fam['anchor'];cc=fam['component'];l=r['closure_wall']
  assert len(cc)==2 and a in cc and all(G[i][i]==2 for i in cc)and G[cc[0]][cc[1]]==-1
  assert l not in E and all(G[l][i]==0 for i in cc if i!=a)
  left=list(map(F,corner['primitive_ray']));right=list(map(F,edge['endpoint_primitive_ray']))
  left=P.scale(-1/P.dot(left,hc),left);right=P.scale(-1/P.dot(right,hc),right)
  B=[[P.pair(x,y,Q)for y in[left,right]]for x in[left,right]]
  p=[P.dot(x,C[l])for x in[left,right]];pp=P.dot(p,p);assert pp>0
  gamma=F(1,12);k=gamma*G[l][a];z0=P.scale(k/pp,p);y=[p[1],-p[0]]
  aa=P.dot(y,P.mv(B,y));bb=2*P.dot(z0,P.mv(B,y));dd=P.dot(z0,P.mv(B,z0))-gamma
  disc=bb*bb-4*aa*dd
  # The scalar quadratic has at most two roots. Reconstruct the named sign,
  # even when the floating-point scout split a repeated rational root.
  if aa:
   assert disc>=0
   root=R.add(R.C(-bb/(2*aa)),R.scale(R.sqrtq(disc)if disc else{},F(-1 if r['quadratic_branch']==0 else 1,2)/aa))
  else:
   assert bb;root=R.C(-dd/bb)
  zc=[R.add(R.C(z0[i]),R.scale(root,y[i]))for i in range(2)]
  vals=R.mv(B,zc);norm=R.pair(zc,R.mv(B,zc));assert norm==R.C(gamma)
  assert R.pair(zc,R.ratvec(p))==R.C(k)
  signs=[R.sign(x)for x in vals]
  assert 0 in signs and disc==0
  endpoint=signs.index(0)
  rows.append({'scout_candidate':r['index'],'family':r['family'],'closure_wall':l,
   'quadratic_coefficients':[aa,bb,dd],'exact_discriminant':disc,'edge_plane_coefficients':zc,
   'endpoint_pairings':vals,'endpoint_signs':signs,'zero_endpoint':endpoint,
   'status':'EXACT_DOUBLE_ROOT_AT_OLD_ENDPOINT_NOT_TRANSVERSE_INTERIOR'})
 sources=[SUMMARY,f'results/{TAG}_scout.json.gz',P.RAW,P.WALK,
  f'scripts/{TAG}_remaining.py','scripts/fresh_h20_cusp_sections_20261013.py',
  'scripts/fresh_h22_multiquadratic_20261013_review.py']
 return enc({'status':'PASS_ALL_66_NUMERICAL_REMAINDERS_ARE_EXACT_ENDPOINT_TANGENCIES',
  'sources':{p:sha(p)for p in sources},'cases':rows,
  'scope':'All66 unresolved recipes of this numerical census are exactly boundary double roots, not transverse interior sections. The481016 other numerical ridge rejections are not exactified by this certificate.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();z=run();out=ROOT/f'results/{TAG}_remaining.json'
 if args.verify:assert z==json.loads(out.read_text())
 else:out.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status'])
