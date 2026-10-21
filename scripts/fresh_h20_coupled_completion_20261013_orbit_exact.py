#!/usr/bin/env python3
"""Exact coherent72-diagram-symmetry cap orbit and its remaining strict escape."""
import argparse,json
from fractions import Fraction as F
from collections import Counter
from math import isqrt
from fresh_h20_coupled_completion_20261013_dotted_exact import setup as setup104,OUT as SOURCE
from fresh_h20_coupled_completion_20261013 import ROOT,read,sha,num
from fresh_h20_gap_next_20261013_exact import root_checks
from fresh_h20_gap_next_20261013_cusp_diagnostic import pair
from fresh_dimension_bridge_20261013_verify import E,Z,O,vec,parse,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv,same
AUTOS='data/enumerations/fresh_h20_coupled_completion_20261013_automorphisms.json'
SCOUT='data/enumerations/fresh_h20_coupled_completion_20261013_orbit_scout.json'
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
OUT='data/enumerations/fresh_h20_coupled_completion_20261013_orbit_exact.json'
def inv(f,x):
 assert f.sign(x)
 if f is B:return O/x
 a,b=x;den=f.base.sub(f.base.mul(a,a),f.base.mul(f.d,f.base.mul(b,b)));di=inv(f.base,den)
 return f.base.mul(a,di),f.base.neg(f.base.mul(b,di))
def div(f,x,y):return f.mul(x,inv(f,y))
def kernel(f,rows):
 a=[list(r)for r in rows];n=len(a[0]);piv=[];k=0
 for j in range(n):
  v=next((i for i in range(k,len(a))if f.sign(a[i][j])),None)
  if v is None:continue
  a[k],a[v]=a[v],a[k];z=a[k][j];a[k]=[div(f,x,z)for x in a[k]]
  for i in range(len(a)):
   if i==k:continue
   z=a[i][j]
   if f.sign(z):a[i]=[f.sub(x,f.mul(z,y))for x,y in zip(a[i],a[k],strict=True)]
  piv.append(j);k+=1
  if k==len(a):break
 assert len(piv)==n-1
 free=next(j for j in range(n)if j not in piv);v=[f.zero]*n;v[free]=f.one
 for i,j in enumerate(piv):v[j]=f.neg(a[i][free])
 assert all(f.sign(fdot(f,row,v))==0 for row in rows)
 return tuple(v),piv

def point_check(f,v,Q,hc,roots,strict):
 Qf=tuple(tuple(map(f.lift,row))for row in Q);vc=fmv(f,Qf,v);n=fdot(f,v,vc);ch=fdot(f,v,tuple(map(f.lift,hc)));assert f.sign(n)>0 and f.sign(ch)<0
 ps=[];zeros=[]
 for i,(g,r)in enumerate(roots):
  H,z=pair(f,vc,g,r);sg=H.sign(z);assert sg<0 if strict else sg<=0
  if not sg:zeros.append(i)
  ps.append({'wall':i,'pairing':H.encode(z),'sign':sg,'outer_radicand':None if H is f else g.d.enc()})
 return {'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':f.d.enc()},'vector':list(map(f.encode,v)),'norm':f.encode(n),'center_pairing':f.encode(ch),'all_wall_pairings':ps,'zero_indices':zeros}

def setup():
 Q,h,hc,roots=setup104();c=read(SOURCE)['new_cap'];f=Extension(B,parse(c['field']['radicand']));r=tuple(map(f.parse,c['normal']));roots.append((f,r))
 return Q,h,hc,roots,f,r

def run():
 Q,h,hc,roots,f,r=setup()
 # N(D) is not a rational square, so D cannot be a square in Qsqrt6.
 # Thus the two-coefficient field encodings used for distinctness are injective.
 nd=f.d.a*f.d.a-6*f.d.b*f.d.b;assert nd==F(41700835,48)
 assert isqrt(nd.numerator)**2!=nd.numerator or isqrt(nd.denominator)**2!=nd.denominator
 base=read(BASE);W=list(map(vec,base['old_normals']));perms=read(AUTOS)['automorphisms'];Pset={tuple(p)for p in perms};assert len(Pset)==72
 assert all(sorted(p)==list(range(36))for p in perms)
 assert all(tuple(p[q[i]]for i in range(36))in Pset for p in perms for q in perms)
 G=[[dot(a,mv(Q,b))for b in W]for a in W];eps=dict(base['height_signs']);full=[];allmats=[]
 for pi,p in enumerate(perms):
  assert all(G[p[i]][p[j]]==G[i][j]for i in range(36)for j in range(36))
  sg=eps[p[24]]//eps[24];assert all(eps[p[i]]==sg*eps[i]for i in eps)
  M=[list(W[p[i]])for i in range(20)]+[[Z]*20+[E(sg)]]
  MC=[mv(Q,row)for row in M];assert all(dot(M[i],MC[j])==Q[i][j]for i in range(21)for j in range(21))
  v=tuple(fdot(f,r,tuple(f.lift(M[k][j])for k in range(21)))for j in range(21));full.append(v);allmats.append(M)
 assert len({tuple(f.encode(x)[0][0]+'|'+str(f.encode(x))for x in v)for v in full})==72
 ids=[c['permutation_index']for c in read(SCOUT)['accepted']];assert len(ids)==70 and len(set(ids))==70
 identity=next(i for i,v in enumerate(full)if all(same(f,a,b)for a,b in zip(v,r,strict=True)));assert identity not in ids
 other=set(range(72))-set(ids)-{identity};assert len(other)==1
 records=[]
 for j,pi in enumerate(ids,105):
  v=full[pi];q,ch,ps=root_checks(f,v,Q,hc,roots);assert same(f,q,f.rational(2))
  records.append({'wall':j,'permutation_index':pi,'normal':list(map(f.encode,v)),'norm':f.encode(q),'strict_center_pairing':f.encode(ch),'old_pair_checks':ps,'pair_counts':dict(sorted(Counter(x['kind']for x in ps).items()))});roots.append((f,v))
 active=read(SCOUT)['all_accepted_escape']['active_indices'];Qf=tuple(tuple(map(f.lift,row))for row in Q);active_roots=[]
 for i in active:
  g,v=roots[i]
  if g is B:active_roots.append(tuple(map(f.lift,v)))
  else:assert g.d==f.d;active_roots.append(v)
 v,piv=kernel(f,[fmv(f,Qf,r)for r in active_roots]);ch=fdot(f,v,tuple(map(f.lift,hc)))
 if f.sign(ch)>0:v=tuple(map(f.neg,v))
 hyper=point_check(f,v,Q,hc,roots,False);strict=None
 for n in range(12,81):
  w=tuple(f.add(x,f.lift(E(F(1,2**n))*y))for x,y in zip(v,h,strict=True))
  try:strict=point_check(f,w,Q,hc,roots,True)
  except AssertionError:continue
  strict['perturbation_power']=n;break
 assert strict is not None
 return {'status':'EXACT_175_WALL_H20_COHERENT_DIAGRAM_ORBIT_WITH_STRICT_ESCAPE','source':SOURCE,'source_sha256':sha(SOURCE),'automorphism_source':AUTOS,'automorphism_sha256':sha(AUTOS),'active_source':SCOUT,'active_source_sha256':sha(SCOUT),'base_source':BASE,'base_sha256':sha(BASE),'group_order':72,'verified_group_closure':True,'orbit_size':72,'identity_permutation_index':identity,'unselected_orbit_indices':sorted(other),'retained_walls':105,'wall_count':175,'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':f.d.enc()},'new_caps':records,'final_hyperideal_point':dict(hyper,selected_active_indices=active,kernel_rank=len(piv)),'strict_escape':strict,'scope':'A72-element exact diagram-isometry group transports the new cap to72 distinct mirrors. Retain the original and70 further mutually compatible mirrors, checking all9765 added pairs exactly. One orbit entry is omitted; no maximality claim. The175-wall chamber remains infinite volume by an exact strict spacelike escape. This is not a proof about arbitrary later caps or other symmetry orbits.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print('escape',d['strict_escape']['norm'])
