#!/usr/bin/env python3
"""Changed19short-block cap outside the closed fixed17block family."""
import argparse,json
from fractions import Fraction as F
from collections import Counter
from fresh_h20_four_cusp_completion_20261013_gap_exact import (
 construct as previous,ROOT,read,sha,OUT as SOURCE,SOURCE as CUSPS,B,E,Z,O,vec,dot,mv,encoded,
 pairing,norm,e,label,strict_point,bind_dependencies,POLAR)
from fresh_h20_coupled_completion_20261013_second_exact import emb
from fresh_dimension_bridge_two_angle_20261013_independent import Extension,fdot,fmv
from fresh_dimension_bridge_20261013_verify import inverse,inertia
OUT='data/enumerations/fresh_h20_four_cusp_completion_20261013_next_exact.json'
S=[0,1,2,3,4,5,8,9,10,13,15,16,18,19,20,21,22,23,24]
RAW=[3562453098,4477180785,0,3562453098,3562453098,0,3562453098,3562453098,0,3562453098,3562453098,2238590392,0,0,0,3562453098,3562453098,0,3562453098,3562453098,1781226549]
def construct():
 Q,h,hc,roots,r=previous();roots.append((B,r));assert len(roots)==225
 f=Extension(B,E(7));u=tuple(E(F(x,3))for x in POLAR[:-1])+(Z,)
 r=tuple(e(f,x)for x in u[:-1])+((Z,E(F(-1,3))),)
 return Q,h,hc,roots,f,u,r

def run():
 dependencies=bind_dependencies(SOURCE);Q,h,hc,roots,f,u,r=construct();q=norm(f,r,Q);assert f.sign(f.sub(q,f.rational(1)))==0
 assert dot(u,mv(Q,u))==E(F(-4,3))and u[-1]==Z
 assert all(dot(roots[i][1],mv(Q,u))==Z for i in S)
 G=[[dot(roots[i][1],mv(Q,roots[j][1])).a for j in S]for i in S];assert inertia(G)==(19,0,0);I=inverse(G);assert all(x>=0 for row in I for x in row)
 rc=fmv(f,tuple(tuple(e(f,x)for x in row)for row in Q),r);cp=fdot(f,r,tuple(e(f,x)for x in hc));assert f.sign(cp)<0
 pairs=[]
 for i,(g,v)in enumerate(roots):
  H,p=pairing(f,rc,g,v);qo=norm(g,v,Q)
  if H is f:qold=emb(f,g,qo);qnew=q
  else:
   assert H.base is f and g.base is B
   qold=(e(f,qo[0]),e(f,qo[1]));qnew=H.lift(q)
  kind=label(H,p,qnew,qold);pairs.append({'old':i,'pairing':H.encode(p),'kind':kind})
 assert all(pairs[i]['kind']=='2'for i in S)and pairs[17]['kind']=='4'and pairs[25]['kind']=='4'
 target=vec(read(SOURCE)['strict_escape']['vector']);cut=fdot(f,r,tuple(e(f,x)for x in mv(Q,target)));assert f.sign(cut)>0
 cusp_pairs=[]
 for row in read(CUSPS)['four_completed_links']:
  delta=vec(row['null_ray']);p=fdot(f,r,tuple(e(f,x)for x in mv(Q,delta)));assert f.sign(p)<0;cusp_pairs.append(f.encode(p))
 roots.append((f,r))
 # Project the center onto the PD19 short block and the orthogonal new root.
 hs=[dot(roots[i][1],hc)for i in S];coeff=[sum((E(a)*b for a,b in zip(row,hs,strict=True)),Z)for row in I]
 p0=tuple(h[k]-sum((a*roots[i][1][k]for a,i in zip(coeff,S,strict=True)),Z)for k in range(21))
 point=tuple(f.sub(e(f,x),f.mul(cp,v))for x,v in zip(p0,r,strict=True));pc=fmv(f,tuple(tuple(e(f,x)for x in row)for row in Q),point);pq=fdot(f,point,pc);ph=fdot(f,point,tuple(e(f,x)for x in hc));assert f.sign(pq)<0 and f.sign(ph)<0
 vertexpairs=[];incident=S+[225]
 for i,(g,v)in enumerate(roots):
  H,p=pairing(f,pc,g,v);assert H.sign(p)==(0 if i in incident else-1);vertexpairs.append({'wall':i,'pairing':H.encode(p)})
 p=tuple(E(F(x,10**10))for x in RAW)
 escape=None
 for k in range(3,13):
  x=tuple(a+E(F(1,10**k))*b for a,b in zip(p,h,strict=True))
  try:escape=strict_point(x,Q,hc,roots)
  except AssertionError:continue
  escape['center_perturbation_power10']=k;break
 assert escape is not None
 return {'status':'PASS_EXACT_226_WALL_CHANGED_BLOCK_CAP_FINITE_VERTEX_AND_ESCAPE','source':SOURCE,'source_sha256':sha(SOURCE),'bound_prior_dependencies':dependencies,'wall_count':226,'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':['7','0']},
         'formula':{'u':encoded(u),'u_square':['-4/3','0'],'root':'u-sqrt(7)*z/3','short_orthogonality_indices':S,'anchors':[17,25],'anchor_pairing':'-1'},
         'new_cap':{'wall':225,'normal':list(map(f.encode,r)),'norm':f.encode(q),'strict_center_pairing':f.encode(cp),'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(x['kind']for x in pairs).items()))},
         'previous_escape_cut':f.encode(cut),'four_completed_cusp_strict_pairings':cusp_pairs,'actual_finite_vertex':{'incident':incident,'short_gram':[[str(x)for x in row]for row in G],'short_inertia':[19,0,0],'point':list(map(f.encode,point)),'norm':f.encode(pq),'center_pairing':f.encode(ph),'all_wall_pairings':vertexpairs},'rounded_escape_numerators':RAW,'strict_escape':escape,
         'scope':'One quadratic cap retains all225 facets, cuts their saved strict escape, preserves all four completed cusp links and creates an actual finite vertex. It changes the previous fixed17block. A new exact strict escape proves the226 chamber still has infinite volume; no exhaustive completion or global cusp-census claim.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print(d['new_cap']['pair_counts']);print('escape',d['strict_escape']['norm'])
