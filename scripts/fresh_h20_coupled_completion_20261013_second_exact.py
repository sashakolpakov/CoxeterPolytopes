#!/usr/bin/env python3
"""Exact second coupled cap, after a coherent symmetry orbit closes the old gap."""
import argparse,json
from fractions import Fraction as F
from collections import Counter
from fresh_h20_coupled_completion_20261013_orbit_exact import setup as setup105,OUT as SOURCE,inv,div
from fresh_h20_coupled_completion_20261013 import ROOT,read,sha
from fresh_h20_gap_next_20261013_exact import kernel as base_kernel
from fresh_h20_edge_family_completion_20261013_cusp_repairs import base_root,sumvec
from fresh_dimension_bridge_20261013_verify import E,Z,O,vec,parse,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv,same,label
SCOUT='data/enumerations/fresh_h20_coupled_completion_20261013_second_orbit_scout.json'
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
OUT='data/enumerations/fresh_h20_coupled_completion_20261013_second_exact.json'
def e(f,x):return x if f is B else f.lift(e(f.base,x))
def contains(f,g):
 if f is g:return True
 if f is B:return False
 return contains(f.base,g)
def emb(f,g,x):
 if f is g:return x
 assert f is not B
 return f.lift(emb(f.base,g,x))
def pairing(f,fc,g,r):
 if contains(f,g):return f,fdot(f,fc,tuple(emb(f,g,x)for x in r))
 if contains(g,f):return g,fdot(g,tuple(emb(g,f,x)for x in fc),r)
 assert g.base is B
 H=Extension(f,e(f,g.d));rr=tuple((e(f,a),e(f,b))for a,b in r)
 return H,fdot(H,tuple(map(H.lift,fc)),rr)
def norm(f,r,Q):return fdot(f,r,fmv(f,tuple(tuple(e(f,x)for x in row)for row in Q),r))
def point(g,x,Q,hc,roots,strict):
 xc=fmv(g,tuple(tuple(e(g,z)for z in row)for row in Q),x);q=fdot(g,x,xc);hp=fdot(g,x,tuple(e(g,z)for z in hc));assert g.sign(q)>0 and g.sign(hp)<0
 pairs=[];zero=[]
 for i,(f,r)in enumerate(roots):
  H,p=pairing(g,xc,f,r);sg=H.sign(p);assert sg<0 if strict else sg<=0
  if not sg:zero.append(i)
  pairs.append({'wall':i,'pairing':H.encode(p),'sign':sg})
 return {'vector':list(map(g.encode,x)),'norm':g.encode(q),'center_pairing':g.encode(hp),'all_wall_pairings':pairs,'zero_indices':zero}
def setup():
 Q,h,hc,roots,f,r=setup105();d=read(SOURCE)
 for c in d['new_caps']:roots.append((f,tuple(map(f.parse,c['normal']))))
 assert len(roots)==175
 return Q,h,hc,roots,f

def run():
 Q,h,hc,roots,f=setup();d=read(SOURCE);p=tuple(map(f.parse,d['final_hyperideal_point']['vector']));r18=base_root(*roots[18]);z=fdot(f,p,tuple(map(f.lift,mv(Q,r18))));p=tuple(div(f,f.neg(x),z)for x in p)
 leaves=[r18,tuple(x/E(3)for x in base_root(*roots[30])),tuple(x/E(3)for x in base_root(*roots[46]))]
 assert [[dot(a,mv(Q,b))for b in leaves]for a in leaves]==[[E(2)if i==j else Z for j in range(3)]for i in range(3)]
 ps=[fdot(f,p,tuple(map(f.lift,mv(Q,r))))for r in leaves];assert all(same(f,a,f.rational(b))for a,b in zip(ps,[-1,F(-3,2),-1],strict=True))
 u=sumvec([(O,leaves[0]),(E(0,F(1,2)),leaves[1]),(O,leaves[2])]);assert dot(u,mv(Q,u))==E(7)
 A=norm(f,p,Q);BB=f.mul(f.rational(2),fdot(f,p,tuple(map(f.lift,mv(Q,u)))));CC=f.rational(6);D=f.sub(f.mul(BB,BB),f.mul(f.rational(4),f.mul(A,CC)));assert f.sign(D)>0
 g=Extension(f,D);twice=f.mul(f.rational(2),A);t=(div(f,f.neg(BB),twice),inv(f,twice));assert g.sign(g.sub(t,g.rational(F(353,100))))>0 and g.sign(g.sub(g.rational(F(177,50)),t))>0
 r=tuple(g.add(e(g,x),g.mul(t,g.lift(y)))for x,y in zip(u,p,strict=True));q=norm(g,r,Q);assert same(g,q,g.rational(1));rc=fmv(g,tuple(tuple(e(g,x)for x in row)for row in Q),r);center=fdot(g,r,tuple(e(g,x)for x in hc));assert g.sign(center)<0
 pairs=[]
 for i,(h0,v)in enumerate(roots):
  H,pr=pairing(g,rc,h0,v);q0=norm(h0,v,Q)
  if H is g:qold=emb(g,h0,q0);qnew=q
  else:
   assert H.base is g and h0.base is B
   qold=(e(g,q0[0]),e(g,q0[1]));qnew=H.lift(q)
  kind=label(H,pr,qnew,qold);pairs.append({'old':i,'pairing':H.encode(pr),'kind':kind})
 target=tuple(map(f.parse,d['strict_escape']['vector']));old=point(f,target,Q,hc,roots,True);tc=fmv(f,tuple(tuple(map(f.lift,row))for row in Q),target);cut=fdot(g,r,tuple(map(g.lift,tc)));assert g.sign(cut)>0
 cusps=[]
 for c in read(BASE)['cusp_rows']:
  cv=mv(Q,vec(c['ray']));cp=fdot(g,r,tuple(e(g,x)for x in cv));assert g.sign(cp)<0;cusps.append(g.encode(cp))
 roots.append((g,r));ids=read(SCOUT)['escape']['active_indices'];rows=[mv(Q,base_root(*roots[i]))for i in ids];v,piv=base_kernel(rows)
 if dot(v,hc).sign()>0:v=tuple(-x for x in v)
 hyper=point(B,v,Q,hc,roots,False);escape=None
 for k in range(12,81):
  x=tuple(a+E(F(1,2**k))*b for a,b in zip(v,h,strict=True))
  try:escape=point(B,x,Q,hc,roots,True)
  except AssertionError:continue
  escape['perturbation_power']=k;break
 assert escape is not None
 return {'status':'EXACT_176_WALL_H20_SECOND_COUPLED_CUTTER_WITH_STRICT_ESCAPE','source':SOURCE,'source_sha256':sha(SOURCE),'active_source':SCOUT,'active_source_sha256':sha(SCOUT),'retained_walls':175,'wall_count':176,'field':{'kind':'quadratic','base':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':f.d.enc()},'radicand':f.encode(D)},'formula':{'polar':list(map(f.encode,p)),'polar_norm':f.encode(A),'leaf_walls':[18,30,46],'leaf_scales':['1','1/3','1/3'],'coefficients':[['1','0'],['0','1/2'],['1','0']],'u':'r18+(sqrt6/2)*(r30/3)+r46/3','root':'u+t*p','t_quadratic':[f.encode(A),f.encode(BB),f.encode(CC)],'t':g.encode(t),'t_isolation':['353/100','177/50']},'new_cap':{'wall':175,'normal':list(map(g.encode,r)),'norm':g.encode(q),'strict_center_pairing':g.encode(center),'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(p['kind']for p in pairs).items())),'inherited_cusp_strict_pairings':cusps},'strict_target':dict(old,new_pairing=g.encode(cut)),'final_hyperideal_point':dict(hyper,selected_active_indices=ids,kernel_rank=len(piv)),'strict_escape':escape,'scope':'One exact norm-one cap cuts the strict175-wall escape and has only valid Coxeter or dotted relations to all175 retained walls. Its new finite nonright labels are4 and6. All107 inherited cusp rays are strict. The176 chamber still has an exact strict spacelike escape. The wider72-copy second orbit remains a numerical scout, not part of this exact checkpoint.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print(d['new_cap']['pair_counts']);print('escape',d['strict_escape']['norm'])
