#!/usr/bin/env python3
"""Exact doubled-H19 source lift, simultaneous H20 core replacement and polars.

Old-versus-old geometry is a dependency of the frozen176 certificate. Every
changed/new normal and every new pair, target cut, center and escape is checked
here in exact rational/radical arithmetic. Scouts choose indices only.
"""
import argparse,json
from collections import Counter,deque
from fractions import Fraction as F
from functools import reduce
from math import gcd,lcm

from fresh_h20_coupled_completion_20261013_second_exact import setup,OUT as SOURCE,e,pairing,norm
from fresh_h20_coupled_completion_20261013 import ROOT,read,sha
from fresh_dimension_bridge_20261013_verify import E,Z,O,vec,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,fdot,fmv,label

BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
H19='data/enumerations/vinberg_x3_h19.json'
CROSS='data/enumerations/fresh_h20_changed_assembly_20261013_double_cross_scout.json'
POLARS='data/enumerations/fresh_h20_changed_assembly_20261013_double_polars_scout.json'
OUT='data/enumerations/fresh_h20_changed_assembly_20261013_exact.json'

def encoded(v):return [x.enc()for x in v]
def key(v):return tuple(tuple(x.enc())for x in v)
def vp(x,y,c=O):return tuple(a+c*b for a,b in zip(x,y,strict=True))
def classify(f,p,qa,qb):return label(f,p,qa,qb)

def bind_dependencies(path):
 """Bind the old-old proof to every hash-declared reconstruction input."""
 checked={};pending=[path]
 def add(p,s):
  assert isinstance(p,str) and isinstance(s,str) and len(s)==64
  assert sha(p)==s,(p,'dependency hash mismatch')
  if p not in checked:checked[p]=s;pending.append(p)
 def scan(obj):
  if isinstance(obj,list):
   for v in obj:scan(v)
  elif isinstance(obj,dict):
   for k,v in obj.items():
    if k.endswith('_sha256') and isinstance(v,str):
     p=obj.get(k[:-7])
     if isinstance(p,str):add(p,v)
    if isinstance(v,dict):
     for p,s in v.items():
      if isinstance(p,str) and '/'in p and isinstance(s,str) and len(s)==64 and all(c in'0123456789abcdef'for c in s):add(p,s)
    scan(v)
 checked[path]=sha(path)
 while pending:
  p=pending.pop()
  if p.endswith('.json'):scan(read(p))
 return dict(sorted(checked.items()))

def strict_point(x,Q,hcov,roots):
 xc=mv(Q,x);q=dot(x,xc);hp=dot(x,hcov)
 assert q.sign()>0 and hp.sign()<0
 ps=[]
 for i,(f,r)in enumerate(roots):
  H,p=pairing(B,xc,f,r);assert H.sign(p)<0
  ps.append({'wall':i,'pairing':H.encode(p)})
 return {'vector':encoded(x),'norm':q.enc(),'center_pairing':hp.enc(),'all_wall_pairings':ps}

def run():
 dependencies=bind_dependencies(SOURCE)
 Q,h,hcov,oldroots,_=setup();frozen=read(SOURCE);base=read(BASE)
 assert frozen['wall_count']==176 and len(oldroots)==175
 assert base['source']==H19 and base['source_sha256']==sha(H19)
 original=list(map(vec,base['old_normals']));a=original[24]
 assert dot(a,mv(Q,a))==E(6)
 source_gram=read(H19)['wall_gram']
 assert [[dot(x,mv(Q,y))for y in original]for x in original]==[[E(-int(v))for v in row]for row in source_gram]
 assert all(x.b==0 for row in Q for x in row)
 # All angles at this H19 wall are right or pi/6; the finite-volume double is
 # therefore justified by the even-wall doubling lemma.
 source_angles=[]
 for i,r in enumerate(original):
  if i==24:continue
  p=dot(a,mv(Q,r));q=dot(r,mv(Q,r));kind=label(B,p,E(6),q)
  if kind not in('dotted','parallel'):assert kind in('2','6')
  source_angles.append({'old':i,'kind':kind})
 doubled=[];provenance=[];seen={}
 for i,r in enumerate(original):
  if i==24:continue
  reflected=vp(r,a,-dot(r,mv(Q,a))/E(3))
  for branch,v in enumerate((r,reflected)):
   k=key(v)
   if k not in seen:seen[k]=len(doubled);doubled.append(v);provenance.append({'old':i,'reflected':bool(branch)})
 assert len(doubled)==49
 dn=[dot(r,mv(Q,r))for r in doubled];short=[i for i,q in enumerate(dn)if q==E(2)];long=[i for i,q in enumerate(dn)if q==E(6)]
 assert len(short)==27 and len(long)==22
 doubled_pairs=[]
 for i in range(49):
  for j in range(i):doubled_pairs.append([i,j,label(B,dot(doubled[i],mv(Q,doubled[j])),dn[i],dn[j])])
 adj={i:[j for j in long if j!=i and dot(doubled[i],mv(Q,doubled[j]))==E(-6)]for i in long}
 signs={long[0]:1};queue=deque([long[0]])
 while queue:
  i=queue.popleft()
  for j in adj[i]:
   if j not in signs:signs[j]=-signs[i];queue.append(j)
   assert signs[j]==-signs[i]
 assert len(signs)==22
 U=[doubled[i]for i in short];names=[{'double_index':i,'kind':'short'}for i in short]
 for i in long:
  for tag,height in(('A',2*signs[i]),('B',-signs[i])):
   U.append(doubled[i][:-1]+(E(height),));names.append({'double_index':i,'kind':tag,'height':height})
 assert len(U)==71
 center=vp(h,a,-dot(h,mv(Q,a))/E(6));hc=mv(Q,center)
 assert dot(center,hc).sign()<0 and dot(center,mv(Q,a))==Z
 un=[dot(r,mv(Q,r))for r in U];uc=[mv(Q,r)for r in U]
 assert all(dot(r,hc).sign()<0 for r in U)
 core_pairs=[]
 for i in range(71):
  for j in range(i):core_pairs.append([i,j,label(B,dot(U[i],uc[j]),un[i],un[j])])
 kept=[i for i in read(CROSS)['retained_old']if i>=48]
 assert len(kept)==102 and 175 not in kept
 roots=[(B,r)for r in U];cross=[]
 for oldid in kept:
  f,r=oldroots[oldid];qn=norm(f,r,Q);cp=fdot(f,r,tuple(e(f,x)for x in hc));assert f.sign(cp)<0
  for j,(v,q,cov)in enumerate(zip(U,un,uc,strict=True)):
   H,p=pairing(B,cov,f,r);assert H is f
   cross.append({'new_core':j,'old':oldid,'pairing':f.encode(p),'kind':label(f,p,e(f,q),qn)})
  roots.append((f,r))
 assert len(roots)==173
 oldtarget=vec(frozen['strict_escape']['vector']);cuts=[dot(r,mv(Q,oldtarget))for r in U]
 assert cuts[59].sign()>0
 additions=[];scout=read(POLARS)
 assert scout['retained_old']==kept and len(scout['added_polars'])==6
 for entry in scout['added_polars']:
  r=tuple(E(x)for x in entry['integer_normal']);q=dot(r,mv(Q,r));rc=mv(Q,r)
  assert q.sign()>0 and dot(r,hc).sign()<0
  pairs=[];zeros=[]
  for i,(f,v)in enumerate(roots):
   H,p=pairing(B,rc,f,v);kind=label(f,p,e(f,q),norm(f,v,Q));assert f.sign(p)<=0
   if not f.sign(p):zeros.append(i)
   pairs.append({'old':i,'pairing':f.encode(p),'kind':kind})
  assert sum(p['kind']not in('dotted','parallel')for p in pairs)>=20
  # The polar itself is feasible, and a small move toward the strict center
  # gives an exact strict old-chamber target that it still cuts.
  target=None
  for power in range(12,65):
   x=vp(r,center,E(F(1,2**power)))
   try:proof=strict_point(x,Q,hc,roots)
   except AssertionError:continue
   cut=dot(r,mv(Q,x))
   if cut.sign()>0:target=dict(proof,new_pairing=cut.enc(),perturbation_power=power);break
  assert target is not None
  additions.append({'wall':len(roots),'normal':encoded(r),'norm':q.enc(),'old_pair_checks':pairs,
                    'pair_counts':dict(sorted(Counter(p['kind']for p in pairs).items())),'zero_indices':zeros,'strict_target':target})
  roots.append((B,r))
 assert len(roots)==179
 proposal=scout['stages'][-1]['escape']['vector'];escape=None
 for power in range(5,13):
  x=tuple(E(F(round(float(z)*10**10),10**10))+E(F(1,10**power))*c for z,c in zip(proposal,center,strict=True))
  try:escape=strict_point(x,Q,hc,roots)
  except AssertionError:continue
  escape['rational_rounding_denominator']=10**10;escape['center_perturbation_power']=power;break
 assert escape is not None
 return {'status':'PASS_EXACT_179_WALL_H20_CHANGED_CORE_WITH_STRICT_ESCAPE',
         'sources':{p:sha(p)for p in(SOURCE,BASE,H19,CROSS,POLARS)},'bound_prior_dependencies':dependencies,'source_operation':'Double VinbergX3 H19 across original long wall24; paired lift in the opposite original long-wall phase; retain102 later compatible normals; add6 exact rational gap polars.',
         'doubling_wall':24,'source_doubling_pair_labels':source_angles,'doubled_H19_normals':list(map(encoded,doubled)),
         'doubled_provenance':provenance,'doubled_H19_pair_counts':dict(sorted(Counter(x[2]for x in doubled_pairs).items())),
         'doubled_short_indices':short,'doubled_long_indices':long,'long_height_signs':[[i,signs[i]]for i in long],
         'new_core_names':names,'new_core_normals':list(map(encoded,U)),'new_core_norms':list(map(E.enc,un)),
         'new_core_pair_counts':dict(sorted(Counter(x[2]for x in core_pairs).items())),
         'retained_old_indices':kept,'new_core_retained_old_pairs':cross,'new_core_retained_pair_counts':dict(sorted(Counter(x['kind']for x in cross).items())),
         'strict_center':encoded(center),'strict_center_norm':dot(center,hc).enc(),
         'previous176_escape_cut':{'new_core_wall':59,'pairing':cuts[59].enc(),'source':SOURCE},
         'new_polars':additions,'wall_count':179,'strict_escape':escape,
         'scope':'Exact changed-core179-facet acute Coxeter chamber in H20. Old-old pairs depend on the frozen176 certificate. All new pairs, centers, sequential cuts and final strict escape are checked here. The last old facet175 is removed; this is not an additive extension of the closed family. Final volume is infinite. No complete new-cusp census or finite-volume conclusion is inferred from the doubled H19 source.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print('core',d['new_core_pair_counts'],'retained',len(d['retained_old_indices']));print('escape',d['strict_escape']['norm'])
