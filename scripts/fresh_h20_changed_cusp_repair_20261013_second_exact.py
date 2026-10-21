#!/usr/bin/env python3
"""A rational215th wall, its escape, and four unbounded new cusp links."""
import argparse,json
from fractions import Fraction as F
from collections import Counter
from fresh_h20_changed_cusp_repair_20261013_gap_exact import (
 construct,OUT as SOURCE,ROOT,read,sha,B,E,Z,O,vec,dot,mv,encoded,
 pairing,norm,e,label,strict_point,bind_dependencies,inverse,inertia)
from fresh_h20_changed_assembly_20261013_cusps import kernel1
from fresh_dimension_bridge_20261013_verify import rank
OUT='data/enumerations/fresh_h20_changed_cusp_repair_20261013_second_exact.json'
SCOUT='data/enumerations/fresh_h20_changed_cusp_repair_20261013_second_volume_scout.json'
R=[16,-6,-4,3,2,4,9,11,12,-3,-6,-3,-2,-2,-4,6,6,8,0,-2,1]

def run():
 dependencies=bind_dependencies(SOURCE);Q,h,hc,roots,f,s,u,previous=construct();roots.append((f,previous));assert len(roots)==214
 r=tuple(E(x)for x in R);rc=mv(Q,r);assert dot(r,rc)==O and dot(r,hc).sign()<0
 pairs=[];parallel=[]
 for i,(g,v)in enumerate(roots):
  H,p=pairing(B,rc,g,v);kind=label(H,p,e(H,O),norm(g,v,Q));pairs.append({'old':i,'pairing':H.encode(p),'kind':kind})
  if kind=='parallel':parallel.append(i)
 assert parallel==[30,34,46,50]
 target=vec(read(SOURCE)['strict_escape']['vector']);cut=dot(r,mv(Q,target));assert cut.sign()>0
 roots.append((B,r));source_vertex=tuple(map(f.parse,read(SOURCE)['actual_finite_vertex']['point']));vp=pairing(B,rc,f,source_vertex);assert vp[0].sign(vp[1])<0
 cusps=[]
 for j in parallel:
  assert roots[j][0]is B and norm(B,roots[j][1],Q)==E(9)
  delta=tuple(a+b/E(3)for a,b in zip(r,roots[j][1],strict=True));dc=mv(Q,delta);assert dot(delta,dc)==Z and dot(delta,hc).sign()<0
  inc=[];signs=[]
  for i,(g,v)in enumerate(roots):
   H,p=pairing(B,dc,g,v);sg=H.sign(p);assert sg<=0;signs.append(sg)
   if not sg:inc.append(i)
  assert len(inc)==20 and all(roots[i][0]is B for i in inc)
  normals=[roots[i][1]for i in inc];assert all(x.b==0 for v in normals for x in v)
  assert rank([[x.a for x in v]for v in normals])==20
  G=[[dot(a,mv(Q,b)).a for b in normals]for a in normals];assert inertia(G)==(19,0,1)
  todo=set(range(len(inc)));components=[];affine_rank=spherical_rank=0;recession=None
  while todo:
   seed=min(todo);todo.remove(seed);cc=[seed]
   for i in cc:
    for k in sorted(todo.copy()):
     if G[i][k]<0:cc.append(k);todo.remove(k)
   M=[[G[i][k]for k in cc]for i in cc];sig=inertia(M);row={'walls':[inc[i]for i in cc],'inertia':list(sig)}
   if sig==(len(cc)-1,0,1):
    weights=kernel1(M);v=tuple(sum((E(w)*normals[i][k]for w,i in zip(weights,cc,strict=True)),Z)for k in range(21));nz=next(k for k,x in enumerate(delta)if x!=Z);mu=v[nz]/delta[nz]
    assert mu.sign()>0 and v==tuple(mu*x for x in delta)
    row.update(kind='affine',positive_kernel=list(map(str,weights)),positive_delta_multiple=mu.enc());affine_rank+=len(cc)-1
   else:
    assert sig==(len(cc),0,0);row['kind']='spherical';spherical_rank+=len(cc)
    if recession is None:
     I=inverse(M);recession=tuple(-sum((E(I[k][0])*normals[i][c]for k,i in enumerate(cc)),Z)for c in range(21))
   components.append(row)
  assert affine_rank==1 and spherical_rank==18 and recession is not None
  wc=mv(Q,recession);wq=dot(recession,wc);assert wq.sign()>0 and dot(delta,wc)==Z
  wp=[dot(v,wc)for v in normals];assert all(p.sign()<=0 for p in wp)and any(p.sign()<0 for p in wp)
  cusps.append({'parallel_old_wall':j,'null_ray':encoded(delta),'incident_walls':inc,'all215_pairing_signs':signs,'components':components,'normal_span_rank':20,'total_quotient_rank':19,'positive_affine_rank':1,'spherical_rank':18,'recession_vector':encoded(recession),'recession_square':wq.enc(),'incident_recession_pairings':list(map(E.enc,wp)),'status':'unbounded despite full quotient rank'})
 proposal=read(SCOUT)['vector'];escape=None
 for power in range(5,13):
  x=tuple(E(F(round(float(z)*10**10),10**10))+E(F(1,10**power))*c for z,c in zip(proposal,h,strict=True))
  try:escape=strict_point(x,Q,hc,roots)
  except AssertionError:continue
  escape['center_perturbation_power']=power;escape['rounding_denominator']=10**10;break
 assert escape is not None
 return {'status':'PASS_EXACT_215_WALL_BRANCH_WITH_FOUR_DEFICIENT_NEW_CUSPS','source':SOURCE,'source_sha256':sha(SOURCE),'scout':SCOUT,'scout_sha256':sha(SCOUT),'bound_prior_dependencies':dependencies,'wall_count':215,'new_cap':{'wall':214,'normal':encoded(r),'norm':['1','0'],'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(x['kind']for x in pairs).items()))},'previous_escape_cut':cut.enc(),'previous_finite_vertex_strict_pairing':vp[0].encode(vp[1]),'four_new_parallel_cusps':cusps,'strict_escape':escape,'scope':'The rational norm-one wall cuts the214 escape and preserves every214 facet and its archived finite vertex. It introduces four feasible new parallel cusps whose total quotient rank19 does not imply boundedness: positive affine rank1 plus spherical rank18 leaves explicit recession. A new strict escape proves215 infinite volume. This is an explored branch with explicit completion obligations, not a completed cusp construction or an unrestricted obstruction.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print(d['new_cap']['pair_counts']);print('escape',d['strict_escape']['norm'])
