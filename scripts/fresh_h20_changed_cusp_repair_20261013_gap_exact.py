#!/usr/bin/env python3
"""Exact remote-gap cap outside the seventeen closed cusp planes."""
import argparse,json
from fractions import Fraction as F
from collections import Counter
from fresh_h20_changed_cusp_repair_20261013 import (
 state,ROOT,read,sha,OUT as SOURCE,CUSPS,B,E,Z,O,vec,dot,mv,encoded,
 pairing,norm,e,label,strict_point,bind_dependencies)
from fresh_h20_coupled_completion_20261013_second_exact import emb
from fresh_dimension_bridge_two_angle_20261013_independent import Extension,fdot,fmv
from fresh_dimension_bridge_20261013_verify import inverse,inertia
OUT='data/enumerations/fresh_h20_changed_cusp_repair_20261013_gap_exact.json'
S=[0,1,2,3,5,8,9,13,14,15,18,19,20,21,22,23,24,11,25]
U=[16,-8,-2,1,1,4,9,11,12,-7,-5,-4,-1,-1,-2,5,6,8,-3,-2,0]
POLAR=[72,-24,-18,15,9,18,41,49,54,-11,-25,-12,-9,-9,-18,28,26,36,2,-8,10]

def construct():
 Q,h,hc,roots=state();d=read(SOURCE)
 roots.extend((B,vec(c['normal']))for c in d['new_caps']);assert len(roots)==213
 f=Extension(B,E(34));s=(E(F(144,11)),E(F(32,11)));u=tuple(E(F(x,16))for x in U)
 r=tuple(f.mul(s,e(f,x))for x in u[:-1])+(f.sub(f.rational(F(1,2)),f.mul(f.rational(F(3,32)),s)),)
 return Q,h,hc,roots,f,s,u,r

def run():
 dependencies=bind_dependencies(SOURCE);Q,h,hc,roots,f,s,u,r=construct();q=norm(f,r,Q);assert f.sign(f.sub(q,f.rational(2)))==0
 assert dot(u,mv(Q,u))==E(F(-1,64))and u[-1]==Z
 assert all(dot(roots[i][1],mv(Q,u))==Z for i in S)
 G=[[dot(roots[i][1],mv(Q,roots[j][1])).a for j in S]for i in S];assert inertia(G)==(19,0,0);I=inverse(G);assert all(x>=0 for row in I for x in row)
 assert dot(roots[55][1],mv(Q,u))==E(F(-9,16))and roots[55][1][-1]==E(-2)
 assert f.sign(f.sub(f.add(f.mul(f.rational(11),f.mul(s,s)),f.mul(f.rational(-288),s)),f.rational(1280)))==0
 rc=fmv(f,tuple(tuple(e(f,x)for x in row)for row in Q),r);cp=fdot(f,r,tuple(e(f,x)for x in hc));assert f.sign(cp)<0
 pairs=[]
 for i,(g,v)in enumerate(roots):
  H,p=pairing(f,rc,g,v);qo=norm(g,v,Q)
  if H is f:qold=emb(f,g,qo);qnew=q
  else:
   assert H.base is f and g.base is B
   qold=(e(f,qo[0]),e(f,qo[1]));qnew=H.lift(q)
  kind=label(H,p,qnew,qold);pairs.append({'old':i,'pairing':H.encode(p),'kind':kind})
 assert all(pairs[i]['kind']=='2'for i in S)and pairs[55]['kind']=='3'
 target=vec(read(SOURCE)['strict_escape']['vector']);cut=fdot(f,r,tuple(e(f,x)for x in mv(Q,target)));assert f.sign(cut)>0
 cusp_pairs=[]
 for row in read(CUSPS)['rays']:
  delta=tuple(E(x)for x in row['ray']);p=fdot(f,r,tuple(e(f,x)for x in mv(Q,delta)));assert f.sign(p)<0;cusp_pairs.append(f.encode(p))
 roots.append((f,r))
 # Project the center onto the PD19 short block and the orthogonal new root.
 hs=[dot(roots[i][1],hc)for i in S];coeff=[sum((E(a)*b for a,b in zip(row,hs,strict=True)),Z)for row in I]
 p0=tuple(h[k]-sum((a*roots[i][1][k]for a,i in zip(coeff,S,strict=True)),Z)for k in range(21))
 point=tuple(f.sub(e(f,x),f.mul(f.mul(f.rational(F(1,2)),cp),v))for x,v in zip(p0,r,strict=True));pc=fmv(f,tuple(tuple(e(f,x)for x in row)for row in Q),point);pq=fdot(f,point,pc);ph=fdot(f,point,tuple(e(f,x)for x in hc));assert f.sign(pq)<0 and f.sign(ph)<0
 vertexpairs=[];incident=S+[213]
 for i,(g,v)in enumerate(roots):
  H,p=pairing(f,pc,g,v);assert H.sign(p)==(0 if i in incident else-1);vertexpairs.append({'wall':i,'pairing':H.encode(p)})
 p=tuple(E(x)for x in POLAR);assert dot(p,mv(Q,p))==E(180)
 assert all(pairing(B,mv(Q,p),g,v)[0].sign(pairing(B,mv(Q,p),g,v)[1])<=0 for g,v in roots)
 escape=None
 for k in range(12,65):
  x=tuple(a+E(F(1,2**k))*b for a,b in zip(p,h,strict=True))
  try:escape=strict_point(x,Q,hc,roots)
  except AssertionError:continue
  escape['center_perturbation_power']=k;break
 assert escape is not None
 return {'status':'PASS_EXACT_214_WALL_GAP_CUTTER_VERTEX_AND_STRICT_ESCAPE','source':SOURCE,'source_sha256':sha(SOURCE),'bound_prior_dependencies':dependencies,'wall_count':214,'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':['34','0']},
         'formula':{'u':encoded(u),'u_square':['-1/64','0'],'s':f.encode(s),'s_equation':[11,-288,-1280],'root':'s*u+(1/2-3*s/32)*z','short_orthogonality_indices':S,'anchor':55,'anchor_pairing':'-3'},
         'new_cap':{'wall':213,'normal':list(map(f.encode,r)),'norm':f.encode(q),'strict_center_pairing':f.encode(cp),'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(x['kind']for x in pairs).items()))},
         'previous_escape_cut':f.encode(cut),'all_204_cusp_strict_pairings':cusp_pairs,'actual_finite_vertex':{'incident':incident,'short_gram':[[str(x)for x in row]for row in G],'short_inertia':[19,0,0],'point':list(map(f.encode,point)),'norm':f.encode(pq),'center_pairing':f.encode(ph),'all_wall_pairings':vertexpairs},'next_feasible_polar':POLAR,'strict_escape':escape,
         'scope':'One quadratic cap retains all213 facets, cuts their saved strict escape, is strict at all204 inherited cusps and creates an actual finite vertex. It is outside the seventeen cusp planes. A new exact strict escape proves the214 chamber still has infinite volume; no exhaustive completion claim.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print(d['new_cap']['pair_counts']);print('escape',d['strict_escape']['norm'])
