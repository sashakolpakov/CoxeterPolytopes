#!/usr/bin/env python3
"""Exact all-radius incenter bound and half-point cutters for balanced endpoints."""
import argparse,json,hashlib
from fractions import Fraction as F
from pathlib import Path
import fresh_h20_balanced_endpoint_20261014 as E
ROOT=E.ROOT;TAG='fresh_h20_balanced_endpoint_20261014_bounds'
SRC='results/fresh_h20_balanced_endpoint_20261014.json'
sub=lambda a,b:E.add(a,E.scale(b,-1))
def dec(x):return tuple(map(F,x))
def inv(a):
 cols=[E.mul(a,tuple(F(i==j)for i in range(8)))for j in range(8)]
 M=[[cols[j][i]for j in range(8)]+[F(i==0)]for i in range(8)]
 for j in range(8):
  k=next(i for i in range(j,8)if M[i][j]);M[j],M[k]=M[k],M[j]
  q=M[j][j];M[j]=[x/q for x in M[j]]
  for i in range(8):
   if i!=j:
    q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[j],strict=True)]
 out=tuple(r[-1]for r in M);assert E.mul(a,out)==E.c(1)
 return out
def div(a,b):return E.mul(a,inv(b))
def pmul(a,b):
 out=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out
def peval(p,x):
 v=F(0)
 for a in p[::-1]:v=v*x+a
 return v
def run():
 d=json.loads((ROOT/SRC).read_text())
 for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 typ=d['types'][0];I=[[F(x)for x in r]for r in typ['inverse_Gram']];ids=typ['edge_support']
 c4=(F(0),F(1,2),F(0),F(0),F(0),F(0),F(0),F(0))
 phi=(F(1,2),F(0),F(0),F(0),F(1,2),F(0),F(0),F(0))
 assert E.mul(c4,c4)==E.c(F(1,2))
 assert E.mul(phi,phi)==E.add(E.c(1),phi)
 scales=[]
 for name,s in [('unit',E.c(1)),('golden_third_factor',phi)]:
  A=[E.c(7),E.c(7),E.scale(s,8)]
  N=E.add(E.c(112),E.scale(E.mul(s,s),84));records=[]
  for profile in typ['profiles']:
   Bs=[dec(r['subtract'])for r in profile['omitted_contact_ratios']]
   cos=[dec(x)for x in profile['cosine_vector']];cost=dec(profile['normalized_spatial_cost'])
   bounds=[div(a,E.add(b,c4))for a,b in zip(A,Bs)]
   base=bounds[0]
   for r in bounds[1:]:
    if E.sign(sub(r,base))<0:base=r
   candidates=[base];exceptions=[]
   for factor,(a,b)in enumerate(zip(A,Bs)):
    for target in [E.ZERO,E.c(F(1,2))]:
     den=E.add(b,target)
     if den==E.ZERO:continue
     r=div(a,den);ratios=[sub(div(aa,r),bb)for aa,bb in zip(A,Bs)]
     good=all(v in(E.ZERO,E.c(F(1,2)))or E.sign(sub(v,c4))>=0 for v in ratios)
     exceptions.append(dict(factor=factor,target_cosine=target,radius=r,ratios=ratios,retained=good))
     if good:candidates.append(r)
   R=candidates[0]
   for r in candidates[1:]:
    if E.sign(sub(r,R))>0:R=r
   zy=E.ZERO
   for i,w in enumerate(ids):
    scale=s if w>=14 else E.c(1)
    for j in range(19):zy=E.add(zy,E.scale(E.mul(scale,cos[j]),2*I[i][j]))
   assert E.sign(zy)>=0 and E.sign(sub(E.c(1),cost))>=0
   power=E.add(sub(N,E.scale(E.mul(R,zy),2)),E.mul(sub(cost,E.c(1)),E.mul(R,R)))
   assert E.sign(sub(power,E.c(68)))>=0
   if name!='unit':assert E.sign(sub(power,E.c(75)))>0
   assert E.sign(sub(E.c(14),R))>=0
   records.append(dict(profile_id=profile['id'],continuous_radius_bound=base,
       exceptional_radii=exceptions,all_radius_upper_bound=R,vertex_dot_direction=zy,
       normalized_spatial_cost=cost,incenter_ball_power_lower_bound=power))
  scales.append(dict(name=name,factor_scales=[E.c(1),E.c(1),s],all_vertex_square=N,profiles=records))
 # All-right unit-scale cutter: c=cos(pi/7), 7/8<c<1, R=7/c.
 # T7(c)+1=(c+1)*(8c^3-4c^2-4c+1)^2 identifies the largest root.
 poly=[F(1),F(-4),F(-4),F(8)]
 assert pmul([F(1),F(1)],pmul(poly,poly))==[1,-7,0,56,0,-112,0,64]
 assert peval(poly,F(7,8))<0<peval(poly,F(1))
 assert 24*F(7,8)**2-8*F(7,8)-4>0 and 48*F(7,8)-8>0
 # Golden all-right cutter: R^2=98 and both omitted A6 contacts pi/4.
 Ng=E.add(E.c(112),E.scale(E.mul(phi,phi),84));R2=E.c(98)
 cutmargin=sub(R2,E.scale(Ng,F(1,4)))
 assert E.sign(cutmargin)>0
 assert E.sign(sub(E.scale(E.mul(phi,phi),64),R2))>0
 paths=[SRC,'scripts/fresh_h20_balanced_endpoint_20261014.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return E.enc(dict(status='PASS_ALL_RADIUS_BALANCED_ENDPOINT_INCENTER_OBSTRUCTION_AND_HALF_POINT_CUTTERS',
  sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
  scales=scales,uniform_incenter_ball_power_lower_bound=68,uniform_radius_upper_bound=14,
  common_strict_timelike_point=dict(spatial_zero=True,D=128,U=1,norm=-256),
  common_strict_escape=dict(spatial_zero=True,D='-1/2',U=1,norm=1,center_pairing='-255/2'),
  unit_half_point_cutter=dict(profile_id=0,cosine_minimal_polynomial=poly,
    cosine_root_interval=['7/8','1'],radius='7/cos(pi/7)',omitted_labels=[7,7,'dotted'],
    half_vertex_distance_squared=49,radius_squared_strictly_greater_than=49),
  golden_half_point_cutter=dict(profile_id=0,radius_squared=98,omitted_labels=[4,4,'dotted'],
    half_vertex_cut_margin=cutmargin),
  scope='All real radii and all omitted-wall labels for the complete41 profiles, transported to392 edges, at unit or(1,1,phi)scales. Incenter is a strict escape for endpoint caps alone. This does not obstruct a central/interior cap, and the old half-vertex witness is explicitly cut.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
