#!/usr/bin/env python3
"""Exact affineD8^2 x affineA1^3 completion of the deficient18/96 cusp."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from fresh_h20_coupled_completion_20261013 import setup,frames,F10,D8IDS,COEFF,SOURCE,DIAG,ROOT,read,sha
from fresh_h20_gap_next_20261013_exact import root_checks,cusp_signs
from fresh_h20_gap_next_20261013_cusp_diagnostic import pair
from fresh_dimension_bridge_20261013_verify import E,Z,O,vec,parse,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv,same,rank
from fresh_h20_edge_family_completion_20261013_exact import positive_solve
from fresh_h20_targeted_escape_20261013_volume import check
OUT='data/enumerations/fresh_h20_coupled_completion_20261013_exact.json'

def plus(a,b):return tuple(F10.add(x,y)for x,y in zip(a,b,strict=True))
def scaled(a,c):return tuple(F10.mul(F10.lift(c),x)for x in a)
def rr(r):return tuple(map(F10.lift,r))
def run():
 Q,h,hc,roots=setup();delta,theta,p,q=frames(Q,h,hc,roots);Qf=tuple(tuple(map(F10.lift,row))for row in Q)
 dc=fmv(F10,Qf,delta);assert F10.sign(fdot(F10,delta,dc))==0 and F10.sign(fdot(F10,delta,rr(hc)))<0
 bases=[tuple(-x for x in t)for t in theta]+[p,tuple(-x for x in p),q,tuple(-x for x in q)]
 out=[]
 for i,v in enumerate(bases,98):
  assert F10.sign(fdot(F10,rr(v),dc))==0
  r=plus(rr(v),scaled(delta,E(F(5))));norm,center,checks=root_checks(F10,r,Q,hc,roots)
  assert same(F10,norm,F10.lift(dot(v,mv(Q,v))))
  cusps=cusp_signs(F10,r,Q)
  out.append({'wall':i,'formula':'base+5*delta','base_normal':[x.enc()for x in v], 'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':['10','0']},'normal':list(map(F10.encode,r)),'norm':F10.encode(norm),'strict_center_pairing':F10.encode(center),'old_pair_checks':checks,'pair_counts':dict(sorted(Counter(t['kind']for t in checks).items())),'inherited_cusp_strict_pairings':cusps})
  roots.append((F10,r))
 incidence=[];pairs=[]
 for i,(g,r)in enumerate(roots):
  H,z=pair(F10,dc,g,r);sg=H.sign(z);assert sg<=0
  if not sg:incidence.append(i)
  pairs.append({'wall':i,'pairing':H.encode(z),'sign':sg,'outer_radicand':None if H is F10 else g.d.enc()})
 components=[]
 for ids,co in [(D8IDS[0]+[98],[F(c,3)if i>=24 else F(c)for i,c in zip(D8IDS[0],COEFF)]+[F(1)]),(D8IDS[1]+[99],[F(c,3)if i>=24 else F(c)for i,c in zip(D8IDS[1],COEFF)]+[F(1)]),([18,96],[F(1),F(1)]),([100,101],[F(1),F(1)]),([102,103],[F(1),F(1)])]:
  rs=[rr(roots[i][1])if roots[i][0]is B else roots[i][1]for i in ids];M=[]
  for r in rs:
   row=[]
   for s in rs:
    z=fdot(F10,r,fmv(F10,Qf,s));assert F10.sign((Z,z[1]))==0;row.append(z[0])
   M.append(row)
  assert all(M[i][j].sign()<=0 for i in range(len(M))for j in range(i))
  assert positive_solve([row[:-1]for row in M[:-1]],[Z]*(len(M)-1))is not None
  assert rank(M)==len(M)-1 and all(sum((x*E(c)for x,c in zip(row,co,strict=True)),Z)==Z for row in M)
  total=F10.zero;sv=tuple(F10.zero for _ in range(21))
  for c,r in zip(co,rs,strict=True):sv=plus(sv,scaled(r,E(c)))
  mu=E(F(5)if len(ids)==9 else F(1)if ids==[18,96]else F(10))
  assert all(same(F10,x,y)for x,y in zip(sv,scaled(delta,mu),strict=True))
  components.append({'walls':ids,'gram':[[x.enc()for x in row]for row in M],'rank':len(M)-1,'positive_kernel':[str(x)for x in co],'kernel_sum_multiple':mu.enc()})
 assert sorted(sum((x['walls']for x in components),[]))==incidence and sum(c['rank']for c in components)==19
 escape=check(vec(read(SOURCE)['strict_escape']['vector']),Q,hc,roots)
 return {'status':'EXACT_104_WALL_H20_FULL_AFFINE_CUSP_COMPLETION_WITH_STRICT_ESCAPE','source':SOURCE,'source_sha256':sha(SOURCE),'diagnostic_source':DIAG,'diagnostic_sha256':sha(DIAG),'retained_walls':98,'wall_count':104,'new_caps':out,'cusp':{'ray':list(map(F10.encode,delta)),'ray_norm':F10.encode(fdot(F10,delta,dc)),'ray_center_pairing':F10.encode(fdot(F10,delta,rr(hc))),'all_wall_pairings':pairs,'incident_walls':incidence,'quotient_rank':19,'components':components,'type':'affineD8 x affineD8 x affineA1^3','bounded_link':True},'strict_escape':escape,'scope':'Six exact compatible walls complete the previously deficient18/96 cusp. Five positive affine components exhaust its19-dimensional Euclidean quotient, proving the cusp link bounded. All603 new pairs and642 inherited-cusp signs are exact. The old strict spacelike escape survives; volume is infinite. No complete census of other cusps.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print([x['pair_counts']for x in d['new_caps']]);print('escape',d['strict_escape']['norm'])
