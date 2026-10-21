#!/usr/bin/env python3
"""Highest-root construction inputs for the four newH20 cusps."""
import argparse,json
from collections import Counter,deque
from fractions import Fraction as F
from fresh_h20_changed_cusp_repair_20261013_second_exact import (
 construct,R,OUT as SOURCE,ROOT,read,sha,B,E,Z,O,vec,dot,mv,encoded,
 pairing,norm,e,label,strict_point,bind_dependencies,inverse,inertia)
from fresh_dimension_bridge_20261013_verify import rank
OUT='data/enumerations/fresh_h20_four_cusp_completion_20261013.json'
POLAR=[48,-16,-12,6,6,12,28,32,36,-10,-14,-8,-6,-6,-12,14,16,24,-2,-4,-7]

def state():
 Q,h,hc,roots,f,s,u,previous=construct();roots.append((f,previous));roots.append((B,tuple(E(x)for x in R)))
 assert len(roots)==215
 return Q,h,hc,roots

def highest(G):
 n=len(G);seen=set();queue=deque()
 for i in range(n):
  v=tuple(int(j==i)for j in range(n));seen.add(v);queue.append(v)
 while queue:
  v=queue.popleft()
  for i in range(n):
   p=sum(v[j]*G[j][i]for j in range(n));w=list(v);w[i]-=p;w=tuple(w)
   if w not in seen:seen.add(w);queue.append(w)
  assert len(seen)<1000
 positive=[v for v in seen if all(x>=0 for x in v)];height=max(map(sum,positive));top=[v for v in positive if sum(v)==height];assert len(top)==1
 assert len(seen)==(90 if n==9 else 60),(n,len(seen),G)
 for v in seen:assert sum(v[i]*G[i][j]*v[j]for i in range(n)for j in range(n))==2
 return list(top[0]),len(seen)

def frames(Q,roots):
 out=[]
 for ci,cusp in enumerate(read(SOURCE)['four_new_parallel_cusps']):
  delta=vec(cusp['null_ray']);dc=mv(Q,delta)
  for component in cusp['components']:
   if component['kind']!='spherical':continue
   ids=component['walls'];simple=[];scales=[]
   for i in ids:
    f,r=roots[i];assert f is B;q=norm(f,r,Q);assert q in(E(2),E(18));scale=O if q==E(2)else E(F(1,3));scales.append(scale);simple.append(tuple(scale*x for x in r))
   G=[[dot(r,mv(Q,s))for s in simple]for r in simple];assert all(x.b==0 and x.a.denominator==1 for row in G for x in row)
   Gi=[[int(x.a)for x in row]for row in G];assert all(Gi[i][j]in((2,)if i==j else(0,-1))for i in range(len(ids))for j in range(len(ids)))
   coefficients,count=highest(Gi);theta=tuple(sum((E(w)*r[j]for w,r in zip(coefficients,simple,strict=True)),Z)for j in range(21));b=tuple(-x for x in theta)
   assert dot(theta,mv(Q,theta))==E(2)and dot(theta,dc)==Z
   assert all(dot(b,mv(Q,roots[i][1])).sign()<=0 for i in cusp['incident_walls'])
   out.append({'cusp':ci,'parallel_wall':cusp['parallel_old_wall'],'delta':delta,'component':ids,'scales':scales,'gram':Gi,'coefficients':coefficients,'root_count':count,'theta':theta,'base':b})
 assert len(out)==9
 return out

def run():
 dependencies=bind_dependencies(SOURCE);Q,h,hc,roots=state();old=roots[:];fs=frames(Q,roots);data=read(SOURCE)
 target=vec(data['strict_escape']['vector']);tc=mv(Q,target);additions=[];lower_bounds=[]
 for k,row in enumerate(fs):
  b=row['base'];delta=row['delta'];bc=mv(Q,b);dc=mv(Q,delta)
  choices=[]
  for j,(f,r)in enumerate(old[:71]):
   assert f is B
   a=-dot(r,dc)
   if a.sign()>0:choices.append((dot(r,bc)/a,j,a,dot(r,bc)))
  assert all(v[0].b==0 for v in choices)
  lower,j,slope,constant=max(choices,key=lambda v:v[0].a);assert lower==E(2)
  lower_bounds.append({'family':k,'old_wall':j,'positive_slope':slope.enc(),'base_pairing':constant.enc(),'necessary_shift_lower_bound':['2','0']})
  r=tuple(x+E(2)*y for x,y in zip(b,delta,strict=True));rc=mv(Q,r);assert dot(r,rc)==E(2)and dot(r,hc).sign()<0
  pairs=[]
  for j,(f,v)in enumerate(roots):
   H,p=pairing(B,rc,f,v);kind=label(H,p,e(H,E(2)),norm(f,v,Q));pairs.append({'old':j,'pairing':H.encode(p),'kind':kind})
  additions.append({'wall':len(roots),'family':k,'cusp':row['cusp'],'parallel_old_wall':row['parallel_wall'],'component':row['component'],'simple_scales':list(map(E.enc,row['scales'])),'simple_gram':row['gram'],'highest_root_coefficients':row['coefficients'],'finite_root_orbit_size':row['root_count'],'theta':encoded(row['theta']),'null_ray':encoded(delta),'shift':['2','0'],'normal':encoded(r),'norm':['2','0'],'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(p['kind']for p in pairs).items())),'previous_escape_pairing':dot(r,tc).enc()})
  roots.append((B,r))
 assert len(roots)==224 and sum(dot(vec(c['normal']),tc).sign()>0 for c in additions)==1
 links=[]
 for ci,cusp in enumerate(data['four_new_parallel_cusps']):
  delta=vec(cusp['null_ray']);dc=mv(Q,delta);weights={};own=[c for c in additions if c['cusp']==ci]
  for component in cusp['components']:
   if component['kind']=='affine':
    for i,w in zip(component['walls'],component['positive_kernel'],strict=True):weights[i]=F(w)
  for c in own:
   for i,w,scale in zip(c['component'],c['highest_root_coefficients'],map(lambda z:E(F(z[0]),F(z[1])),c['simple_scales']),strict=True):
    assert scale.b==0;weights[i]=F(w)*scale.a
   weights[c['wall']]=F(1)
  expected=set(cusp['incident_walls'])|{c['wall']for c in own};assert set(weights)==expected and all(w>0 for w in weights.values())
  actual=[];signs=[]
  for j,(f,r)in enumerate(roots):
   H,p=pairing(B,dc,f,r);sg=H.sign(p);assert sg<=0;signs.append(sg)
   if not sg:actual.append(j)
  assert set(actual)==expected
  rs=[roots[i][1]for i in actual];assert all(roots[i][0]is B for i in actual)and all(x.b==0 for r in rs for x in r)
  assert rank([[x.a for x in r]for r in rs])==20
  total=tuple(sum((E(w)*roots[i][1][j]for i,w in weights.items()),Z)for j in range(21));nz=next(j for j,x in enumerate(delta)if x!=Z);mu=total[nz]/delta[nz];assert mu.sign()>0 and total==tuple(mu*x for x in delta)
  links.append({'cusp':ci,'null_ray':encoded(delta),'incident_walls':actual,'all224_pairing_signs':signs,'positive_weights':[[i,str(weights[i])]for i in actual],'positive_delta_multiple':mu.enc(),'normal_span_rank':20,'positive_quotient_span_rank':19,'status':'bounded'})
 p=tuple(E(x)for x in POLAR);pq=dot(p,mv(Q,p));assert pq.sign()>0
 for f,r in roots:
  H,z=pairing(B,mv(Q,p),f,r);assert H.sign(z)<=0
 escape=None
 for k in range(12,65):
  x=tuple(a+E(F(1,2**k))*b for a,b in zip(p,h,strict=True))
  try:escape=strict_point(x,Q,hc,roots)
  except AssertionError:continue
  escape['center_perturbation_power']=k;break
 assert escape is not None
 xc=mv(Q,vec(escape['vector']));dominance=[]
 for k,row in enumerate(fs):
  delta_pair=dot(row['delta'],xc);minimum_pair=dot(vec(additions[k]['normal']),xc)
  assert delta_pair.sign()<0 and minimum_pair.sign()<0
  dominance.append({'family':k,'null_ray_pairing':delta_pair.enc(),'minimal_cap_pairing':minimum_pair.enc(),'identity':'(x,-theta+t*delta)=(x,-theta+2*delta)+(t-2)*(x,delta)<0 for every t>=2'})
 return {'status':'PASS_EXACT_224_WALL_FOUR_CUSP_COMPLETION_AND_ALL_SHIFT_ESCAPE','source':SOURCE,'source_sha256':sha(SOURCE),'bound_prior_dependencies':dependencies,'wall_count':224,'new_caps':additions,'four_completed_links':links,'sharp_acuteness_lower_bounds':lower_bounds,'strict_escape':escape,'all_real_shift_domination':dominance,
         'scope':'Nine mutually compatible norm-two negative-highest-root caps at the sharp minimal shift2 close all four specified cusp links and cut the215 escape. A new strict escape survives224. Every acute-compatible cap in these nine affine-ray families has shift>=2 and is weaker on the retained215 cone, so no finite collection from the whole all-real-shift family can complete this fixed215 chamber. Other profiles and retained-wall changes are unrestricted. This is not a full cusp census of224 walls.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print('newpairs',Counter(p['kind']for c in d['new_caps']for p in c['old_pair_checks']));print('escape',d['strict_escape']['norm'])
