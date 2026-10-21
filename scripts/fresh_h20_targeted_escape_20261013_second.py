#!/usr/bin/env python3
"""Exact87th cap: preserve the next polar's old projection and set angle pi/3."""
import argparse,json,hashlib
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
from fresh_dimension_bridge_20261013_verify import E,Z,O,parse,vec,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv,label,same

ROOT=Path(__file__).resolve().parents[1]
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
POLAR='data/enumerations/fresh_dimension_bridge_polar_20261013.json'
EXT='data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json'
FIRST='data/enumerations/fresh_h20_targeted_escape_20261013_exact.json'
VOLUME='data/enumerations/fresh_h20_targeted_escape_20261013_volume.json'
OUT='data/enumerations/fresh_h20_targeted_escape_20261013_second.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def solve(A,b):
 n=len(A);a=[list(row)+[v]for row,v in zip(A,b,strict=True)]
 for j in range(n):
  k=next(i for i in range(j,n)if a[i][j].sign());a[j],a[k]=a[k],a[j]
  s=a[j][j];a[j]=[x/s for x in a[j]]
  for i in range(n):
   if i==j:continue
   s=a[i][j]
   if s.sign():a[i]=[x-s*y for x,y in zip(a[i],a[j],strict=True)]
 return tuple(row[-1]for row in a)

def run():
 base,polar,ext,first,volume=map(read,(BASE,POLAR,EXT,FIRST,VOLUME))
 Q=tuple(map(vec,base['form']));h=vec(base['strict_center']);hc=mv(Q,h)
 roots=list(map(vec,base['normals']))+[vec(c['vector'])for c in polar['cuts']]
 cov=[mv(Q,r)for r in roots];ns=[dot(r,c)for r,c in zip(roots,cov,strict=True)]
 ids=volume['certificate']['numerical_active_indices'];assert len(ids)==20 and all(i<82 for i in ids)
 spatial=solve([cov[i][:20]for i in ids],[cov[i][-1]for i in ids])+(Z,)
 polar_ray=spatial[:-1]+(-O,)
 assert all(dot(polar_ray,cov[i])==Z for i in ids)and dot(polar_ray,mv(Q,polar_ray)).sign()>0
 anchor=21;p=dot(spatial,cov[anchor]);assert p.sign()<0
 scale=-O/p;spatial=tuple(scale*x for x in spatial)
 assert dot(spatial,cov[anchor])==-O
 d=(E(F(2))-dot(spatial,mv(Q,spatial)))/E(F(3));assert d.sign()>0
 f=Extension(B,d);r=tuple(f.lift(x)for x in spatial[:-1])+((Z,-O),)
 Qf=tuple(tuple(f.lift(x)for x in row)for row in Q);rc=fmv(f,Qf,r);q=fdot(f,r,rc)
 assert same(f,q,f.rational(2))
 center=fdot(f,r,tuple(map(f.lift,hc)));assert f.sign(center)<0
 w=vec(volume['certificate']['vector']);wc=mv(Q,w)
 cut=fdot(f,r,tuple(map(f.lift,wc)));assert f.sign(cut)>0
 pairs=[]
 for i,(cv,n)in enumerate(zip(cov,ns,strict=True)):
  p=fdot(f,r,tuple(map(f.lift,cv)));kind=label(f,p,q,f.lift(n))
  pairs.append({'old':i,'pairing':f.encode(p),'kind':kind})
 for i,c in enumerate(ext['caps']+[{'field':first['field'],'vector':first['normal'],'norm':None}],82):
  g=Extension(f,f.lift(parse(c['field']['radicand'])))
  other=tuple((f.lift(parse(x[0])),f.lift(parse(x[1])))for x in c['vector'])
  pp=fdot(g,tuple(map(g.lift,rc)),other)
  if i==85:
   # Norm of the firstcap is76−24sqrt10, expressed in the outer radical.
   oq=(f.rational(76),f.rational(-24))
  else:oq=g.lift(f.lift(parse(c['norm'])))
  kind=label(g,pp,g.lift(q),oq);assert kind=='dotted'
  pairs.append({'old':i,'pairing':g.encode(pp),'kind':kind})
 assert len(pairs)==86 and pairs[21]['kind']=='3'
 cusp_pairs=[]
 for c in base['cusp_rows']:
  pp=fdot(f,r,tuple(map(f.lift,mv(Q,vec(c['ray'])))))
  assert f.sign(pp)<0;cusp_pairs.append(f.encode(pp))
 return {'status':'EXACT_87_WALL_H20_TARGETED_SECOND_CUT_VOLUME_PENDING',
  'source_sha256':{p:sha(p)for p in(BASE,POLAR,EXT,FIRST,VOLUME)},
  'construction':'Solve the next polar20equations with lastcoordinate−1; scale its oldprojection to pairing−1 with shortwall21; choose negative height so norm2.',
  'projection_active_indices':ids,'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':d.enc()},
  'normal':list(map(f.encode,r)),'norm':[['2','0'],['0','0']],
  'strict_center_pairing':f.encode(center),'previous_escape_pairing':f.encode(cut),
  'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(a['kind']for a in pairs).items())),
  'inherited_cusp_strict_pairings':cusp_pairs,
  'scope':'Exact norm2, all86 relations, strict center, positive cut of the frozen86strictescape and107inheritedcusps preserved. All87walls actual byacuteprojection. Finitevolume notcertified.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k]for k in ['status','field','pair_counts','previous_escape_pairing']}))
