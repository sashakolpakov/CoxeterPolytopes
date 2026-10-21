#!/usr/bin/env python3
"""Exact obstruction to one cap closing all three old partner edges.

Allows every finite or parabolic endpoint profile and all old contact types.
The obstruction is already in the literal21-root core.
"""
import argparse,hashlib,itertools,json
import sympy as s
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_opposite_prism_gate_20261017'
END='results/fresh_h20_e8_partner_endpoints_20261017.json';AFF='results/fresh_h20_parabolic_endpoint_profiles_20261017.json';CORE='results/fresh_h20_e8_cusp_completion_20261017.json'
def run():
 D=json.loads((ROOT/END).read_text());A=json.loads((ROOT/AFF).read_text());C=json.loads((ROOT/CORE).read_text())
 Ws=[r['support']for r in D['endpoint_censuses']]
 rows=[[dict(r,kind='finite')for r in d['profiles']]+[dict(r,kind='parabolic')for r in a['profiles']]for d,a in zip(D['endpoint_censuses'],A['edges'])]
 third={tuple(tuple(z)for z in r['contacts']):r for r in rows[2]};matches=[]
 for a,b in itertools.product(*rows[:2]):
  da,db=dict(a['contacts']),dict(b['contacts'])
  if any(da.get(j,2)!=db.get(j,2)for j in set(Ws[0])&set(Ws[1])):continue
  d={**da,**db};key=tuple(sorted((j,m)for j,m in d.items()if j in Ws[2]))
  if key not in third:continue
  c=third[key];matches.append(dict(id=len(matches),endpoint_kinds=[x['kind']for x in(a,b,c)],endpoint_profiles=[x['id']for x in(a,b,c)],contacts=sorted([j,m]for j,m in d.items())))
 assert len(matches)==45
 ids=C['core_ids'];G=s.Matrix(C['core_gram']);GI=G.inv();jj=[i for i,j in enumerate(ids)if j!=158];SB=s.eye(21)[:,jj];H=G.extract(jj,jj);J=H.inv();ee=s.eye(21)[:,ids.index(158)]
 u=(ee-SB*J*SB.T*G*ee).applyfunc(s.simplify);sigma=s.simplify((u.T*G*u)[0]);assert sigma==s.Rational(1,30)
 assert (SB.T*G*u).applyfunc(s.simplify)==s.zeros(20,1)
 V=[-GI[:,ids.index(i)]for i in(159,160,161)]
 corners=[]
 for k,v in zip((159,160,161),V):
  assert s.simplify((v.T*G*v)[0]).is_negative
  ii=[i for i,j in enumerate(ids)if j!=k];K=G.extract(ii,ii);L,diag=K.LDLdecomposition(hermitian=False)
  assert all(s.simplify(diag[i,i]).is_positive for i in range(20))
  assert all(s.simplify(x).is_nonnegative for x in K.inv())
  assert s.simplify((u.T*G*v)[0])==0
  corners.append(dict(omitted_partner=k,square=str(s.simplify((v.T*G*v)[0])),positive_LDL_pivots=[str(s.simplify(diag[i,i]))for i in range(20)]))
 output=[];special=[]
 for row in matches:
  f=s.zeros(20,1)
  for j,m in row['contacts']:f[jj.index(ids.index(j))]=-2 if m=='parallel'else -2*s.cos(s.pi/m)
  p=(SB*J*f).applyfunc(s.simplify);q=s.simplify((p.T*G*p)[0]);vprod=[s.simplify((p.T*G*v)[0])for v in V]
  rec=dict(row,particular_square=str(q),old_vertex_products=list(map(str,vprod)),norm_remainder=str(s.simplify(2-q)))
  if all(x==0 for x in vprod):rec['exclusion']='Every norm completion vanishes on the old finite vertices; incompatible with a common strict center and acute old contacts.'
  else:
   assert row['contacts']==[[160,3]]
   assert q==(-5+s.sqrt(2))/23
   assert s.simplify((p.T*G*ee)[0])==0
   sq=s.simplify((2-q)*sigma);assert s.simplify(sq-(51-s.sqrt(2))/690)==0
   assert s.Rational(49,690)>0 and s.Rational(50,690)<1
   rec.update(exclusion='The two contacts with158 are opposite nonzero signs with absolute value less than1: one obtuse, one strictly between right andpi/3.',
    pair158_squared=str(sq),pair158_squared_rational_bounds=['49/690','50/690'],
    two_branch_mutual_pair=str(s.simplify(2*q-2)),particular_coefficients=list(map(str,p)))
   special.append(rec)
  output.append(rec)
 assert len(special)==1
 # The entire unbounded singleton151 family, without choosing its angle.
 p151=(-SB*J*SB.T*s.eye(21)[:,ids.index(151)]).applyfunc(s.simplify)
 assert all(s.simplify((p151.T*G*v)[0])==0 for v in V)
 assert p151==-s.eye(21)[:,ids.index(151)]/2
 counts={k:sum(all(x==k for x in r['endpoint_kinds'])for r in output)for k in('finite','parabolic')}
 assert counts==dict(finite=40,parabolic=5)
 paths=[f'scripts/{TAG}.py',END,AFF,CORE,'results/fresh_h20_e8_exact_cluster_20261017.json','results/fresh_h20_e8_f4_cusps_20261017.json']
 return dict(status='PASS_EXACT_COMPLETE_ONE_OPPOSITE_PRISM_CAP_OBSTRUCTION',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
  fixed_profile_counts=counts,old_finite_corner_checks=corners,perpendicular158_residual_coefficients=list(map(str,u)),perpendicular158_residual_square=str(sigma),
  profiles=output,special_profile=special[0],singleton151_all_label_proof=dict(particular_coefficients=list(map(str,p151)),old_vertex_products=['0','0','0'],domain='Every real magnitude d for which the spherical/affine endpoint Gram is permitted; in particular every integer finite label m>=3 and the parallel endpoint.'),
  scope='No one new coherent acute Coxeter wall can make all three W159160+r,W159161+r,W160161+r endpoint Grams positive definite or PSD corank1 while retaining the old21-root core. No dotted-other restriction and no label cutoff. This does not exclude replacing158 or adding several intermediate facets.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['fixed_profile_counts'],flush=True)
