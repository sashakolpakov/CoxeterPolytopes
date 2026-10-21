#!/usr/bin/env python3
"""Exact future-ray and linear-height gates for all parabolic endpoint profiles.

The five positive projected-norm profiles are tested on relaxed continuous
height domains. No finite-volume construction is claimed.
"""
import argparse, hashlib, json
from fractions import Fraction as F
import sympy as s
import fresh_h20_e8_exact_cluster_20261017 as P
import fresh_h20_e8_endpoint_belt_20261017_exact_gate as E
ROOT=P.ROOT;TAG='fresh_h20_parabolic_endpoint_gate_20261017'
CENSUS='results/fresh_h20_parabolic_endpoint_profiles_20261017.json'
END='results/fresh_h20_e8_partner_endpoints_20261017.json'
pt,add,sub,neg,mul,div,sqrt=P.pt,P.add,P.sub,P.neg,P.mul,P.div,P.sqrt
scale,plus,sx=E.scale,E.plus,E.sx

def build():
 C=json.loads((ROOT/P.CUSP).read_text());D=json.loads((ROOT/END).read_text())
 A=json.loads((ROOT/CENSUS).read_text());Z=json.loads((ROOT/'results/fresh_h20_e8_exact_cluster_20261017.json').read_text())
 H=json.loads((ROOT/P.GEOM).read_text());K=json.loads((ROOT/'results/fresh_h20_e8_f4_cusps_20261017.json').read_text())
 ids=C['core_ids'];G=s.Matrix(C['core_gram']);GI=G.inv();IG=[[P.alg(x)for x in row]for row in G.tolist()]
 def gp(v):return P.mv(IG,v)
 def dot(v,w):return P.inner(v,gp(w))
 Q=[list(map(F,row))for row in H['form']];nz=[(i,j,q)for i,row in enumerate(Q)for j,q in enumerate(row)if q]
 N=[[tuple(map(F,z))for z in v]for v in H['normals']]+[[tuple(map(F,z))for z in r['physical_normal']]for r in Z['steps']]
 qs=list(map(F,H['norms']))+[F(2)]*4
 def physical(v,w):return P.total(mul(pt(q),mul(v[i],w[j]))for i,j,q in nz)
 BP=[[physical(N[i],n)for n in N]for i in ids]
 for cusp in C['cusps']:
  nn=s.Matrix(cusp['null_coefficients'])
  for b in cusp['base_new_root_coefficients']:
   v=s.Matrix(b)+K['selected_common_height']*nn;g=sx(G*v)
   for i in range(21):BP[i].append(g[i])
   qs.append(F(2))
 center=[tuple(map(F,z))for z in K['center']]
 hp=[physical(N[i],center)for i in ids]
 thresholds=[sqrt(pt(2*q))for q in qs]
 bounds={(0,2):(F(29,20),145),(1,0):(F(0),None),(1,1):(F(31,100),160),(1,7):(F(31,100),160),(2,0):(F(81,100),159)}
 rays=[];families=[];edge_data=[]
 for ci,cen in enumerate(D['endpoint_censuses']):
  wi=[ids.index(i)for i in cen['support']];SB=s.eye(21)[:,wi];WJ=s.Matrix(cen['inverse']);ee=s.eye(21)[:,ids.index(158)]
  u=(ee-SB*WJ*SB.T*G*ee).applyfunc(s.simplify);sigma=s.simplify((u.T*G*u)[0]);assert sigma.is_positive
  missing=next(i for i in(159,160,161)if i not in cen['edge']);mi=ids.index(missing)
  V=(-GI[:,mi]/s.sqrt(-GI[mi,mi])).applyfunc(s.simplify)
  unit=(u/s.sqrt(sigma)).applyfunc(s.simplify)
  nu=((V-unit)/s.sqrt(sigma)).applyfunc(s.simplify)
  other=((V+unit)/s.sqrt(sigma)).applyfunc(s.simplify)
  assert s.simplify((nu.T*G*nu)[0])==0 and (SB.T*G*nu).applyfunc(s.simplify)==s.zeros(19,1)
  assert s.simplify((nu.T*G*ee)[0])==-1 and s.simplify((other.T*G*ee)[0])==1
  ni=sx(nu);nc=P.inner(ni,hp);assert nc[1]<0
  products=[P.inner(ni,[BP[k][j]for k in range(21)])for j in range(168)]
  strict=[j for j in range(168)if j not in cen['support']]
  assert all(products[j][1]<0 for j in strict)
  rays.append(dict(edge=ci,coefficients=list(map(str,nu)),future_center_pair=P.save(nc),
    old_pairs={str(j):P.save(products[j])for j in strict},other_future_ray_pair158='1'))
  edge_data.append(dict(ci=ci,cen=cen,SB=SB,WJ=WJ,V=sx(V),unit=sx(unit),sigma=sigma,nu=ni))
  for row in A['edges'][ci]['profiles']:
   if not row['can_meet_target_H3']:continue
   pp=SB*s.Matrix(row['projection_coefficients']);off=s.simplify((pp.T*G*ee)[0]);base=(pp+off*nu).applyfunc(s.simplify)
   assert s.simplify((base.T*G*base)[0])==2 and s.simplify((base.T*G*nu)[0])==0
   assert s.simplify((base.T*G*ee)[0])==0
   L,block=bounds[ci,row['id']];v=sx(base+L*nu)
   witness=None
   if block is not None:
    slope=products[block];test=add(P.inner(v,[BP[k][block]for k in range(21)]),thresholds[block])
    assert slope[1]<0 and test[0]>0
    witness=dict(wall=block,slope=P.save(slope),dotted_test_at_lower_bound=P.save(test))
   f=dict(ci=ci,profile=row['id'],contacts=row['contacts'],base=sx(base),nu=ni,lo=v,L=L,off=off,pp=pp,nu_expr=nu,
          record=dict(edge=ci,profile=row['id'],contacts=row['contacts'],T_residual_norm=row['T_residual_norm'],
           base_coefficients=list(map(str,base)),height_parameter='x=-r.158',lower_bound=str(L),lower_bound_witness=witness))
   families.append(f)
  print('PARABOLIC_EDGE',ci,'families',sum(f['ci']==ci for f in families),flush=True)
 return dict(C=C,D=D,A=A,ids=ids,G=G,GI=GI,IG=IG,gp=gp,dot=dot,BP=BP,hp=hp,qs=qs,thresholds=thresholds,
             rays=rays,families=families,edge_data=edge_data)

def run():
 B=build();fs=B['families'];dot=B['dot'];pairs=[];same=[]
 for i,f in enumerate(fs):
  for g in fs[:i]:
   if f['ci']==g['ci']:
    val=s.simplify((f['pp'].T*B['G']*g['pp'])[0]);assert val.is_positive
    same.append(dict(first=[f['ci'],f['profile']],second=[g['ci'],g['profile']],constant_pair=str(val)))
   else:
    co=[add(dot(f['lo'],g['lo']),pt(2)),dot(f['nu'],g['lo']),dot(f['lo'],g['nu']),dot(f['nu'],g['nu'])]
    assert all(v[1]<0 for v in co)
    pairs.append(dict(first=[f['ci'],f['profile']],second=[g['ci'],g['profile']],
      polynomial='(r(x).s(y)+2), x=L_first+t,y=L_second+u; coefficients 1,t,u,tu',coefficients=[P.save(v)for v in co]))
 paths=[f'scripts/{TAG}.py',CENSUS,END,P.CUSP,P.GEOM,'results/fresh_h20_e8_exact_cluster_20261017.json',
  'results/fresh_h20_e8_f4_cusps_20261017.json',f'scripts/{P.TAG}.py',f'scripts/{E.TAG}.py']
 return dict(status='PASS_EXACT_PARABOLIC_ENDPOINT_CONTINUOUS_HEIGHT_SEPARATION',
  sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},future_rays=B['rays'],
  meeting_families=[f['record']for f in fs],cross_edge_dotted_pairs=pairs,same_edge_obtuse_pairs=same,
  scope='All five positive-T-norm parabolic endpoint profiles, with all unprescribed old168 contacts other than158 strictly dotted and acute158. Continuous height relaxation; no label cutoff. Distinct such meeting profiles cannot be actual neighbor facets. This does not yet classify mixed finite/parabolic pairs or provide a finite-volume construction.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
