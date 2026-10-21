#!/usr/bin/env python3
"""Separate88-wall exact escape gate, preserving the frozen86/87-wall files."""
import argparse,json
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from fresh_h20_targeted_escape_20261013_volume import check,read,sha,ROOT,BASE,POLAR,EXT,CAP,E,parse,vec,Extension,B
from fresh_h20_targeted_escape_20261013_second_volume import extended_setup as prior_setup,SECOND
from fresh_h20_targeted_escape_20261013 import ev
THIRD='data/enumerations/fresh_h20_targeted_escape_20261013_third.json'
OUT='data/enumerations/fresh_h20_targeted_escape_20261013_third_volume.json'

def extended_setup():
 Qn,Rn,hn,Q,h,hc,roots=prior_setup();c=read(THIRD);d=ev(c['field']['radicand'])
 r=np.array([ev(x[0])+np.sqrt(d)*ev(x[1])for x in c['normal']]);Rn=np.vstack([Rn,r])
 f=Extension(B,parse(c['field']['radicand']));roots.append((f,tuple(map(f.parse,c['normal']))))
 return Qn,Rn,hn,Q,h,hc,roots

def run(attempts=200):
 Qn,Rn,hn,Q,h,hc,roots=extended_setup();norms=np.einsum('ij,jk,ik->i',Rn,Qn,Rn)
 cov=Rn@Qn/np.sqrt(norms[:,None]);sh=hn/np.sqrt(-hn@Qn@hn);rng=np.random.default_rng(20262088);tries=[];cert=None
 for i in range(attempts):
  lp=linprog(rng.normal(size=21),A_ub=cov,b_ub=np.zeros(len(Rn)),A_eq=(sh@Qn)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
  if not lp.success:continue
  q=lp.x@Qn@lp.x;ids=np.flatnonzero(abs(cov@lp.x)<1e-7).tolist();tries.append({'attempt':i,'numeric_norm':float(q),'active_indices':ids})
  if q<=1e-7:continue
  scale=max(abs(lp.x));base=tuple(E(F(round(x/scale*10**9),10**9))for x in lp.x)
  for power in range(12,61):
   v=tuple(x+E(F(1,2**power))*z for x,z in zip(base,h,strict=True))
   try:cert=check(v,Q,hc,roots)
   except AssertionError:continue
   cert.update({'source_attempt':i,'perturbation_power':power,'numerical_active_indices':ids});break
  if cert:break
 return {'status':'EXACT_STRICT_ESCAPE_FOR_88_WALL_H20'if cert else'BOUNDED_88_WALL_VOLUME_SCOUT_UNDECIDED',
  'source_hashes':{p:sha(p)for p in(BASE,POLAR,EXT,CAP,SECOND,THIRD)},'attempt_limit':attempts,'attempts':tries,'certificate':cert,
  'scope':'Exact strict future spacelike witness gate against all88 roots. LP onlyproposespoints; noexactpositivevolumeclaim fromnumericalfailure.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 if a.verify:
  out=read(OUT);assert all(sha(p)==h for p,h in out['source_hashes'].items());_,_,_,Q,h,hc,roots=extended_setup();c=out['certificate'];assert c
  actual=check(vec(c['vector']),Q,hc,roots);assert all(actual[k]==c[k]for k in actual)
 else:out=run();(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':out['status'],'attempts':len(out['attempts']),'norm':out['certificate']['norm']if out['certificate']else None}))
