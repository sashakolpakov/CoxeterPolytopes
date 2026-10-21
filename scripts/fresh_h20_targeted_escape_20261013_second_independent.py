#!/usr/bin/env python3
"""Independent Fraction field replay of the 87th H20 wall and strict escape."""
import argparse
from collections import Counter
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
POLAR='data/enumerations/fresh_dimension_bridge_polar_20261013.json'
CAPS='data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json'
FIRST='data/enumerations/fresh_h20_targeted_escape_20261013_exact.json'
OLDVOL='data/enumerations/fresh_h20_targeted_escape_20261013_volume.json'
SOURCE='data/enumerations/fresh_h20_targeted_escape_20261013_second.json'
VOLUME='data/enumerations/fresh_h20_targeted_escape_20261013_second_volume.json'
OUT='results/fresh_h20_targeted_escape_20261013_second_independent.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def add(x,y):return tuple(a+b for a,b in zip(x,y,strict=True))
def neg(x):return tuple(-a for a in x)
def sub(x,y):return add(x,neg(y))
def qm(x,y):return(x[0]*y[0]+6*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def qs(x):
 a,b=x;sa=(a>0)-(a<0);sb=(b>0)-(b<0)
 if not sb:return sa
 if not sa or sa==sb:return sb
 d=a*a-6*b*b;return sa*((d>0)-(d<0))

class Field:
 def __init__(self,d):self.d=tuple(map(Q,d));assert qs(self.d)>0
 def rational(self,x):return(Q(x),Q(0),Q(0),Q(0))
 def base(self,x):return(Q(x[0]),Q(x[1]),Q(0),Q(0))
 def parse(self,x):return tuple(map(Q,x[0]+x[1]))
 def mul(self,x,y):return add(qm(x[:2],y[:2]),qm(self.d,qm(x[2:],y[2:])))+add(qm(x[:2],y[2:]),qm(x[2:],y[:2]))
 def sign(self,x):
  a,b=x[:2],x[2:];sa,sb=qs(a),qs(b)
  if not sb:return sa
  if not sa or sa==sb:return sb
  return sa*qs(sub(qm(a,a),qm(self.d,qm(b,b))))
 def dot(self,x,y):
  out=self.rational(0)
  for a,b in zip(x,y,strict=True):out=add(out,self.mul(a,b))
  return out
 def mv(self,M,x):return tuple(self.dot(r,x)for r in M)
 def extmul(self,x,y,d):return(add(self.mul(x[0],y[0]),self.mul(d,self.mul(x[1],y[1]))),add(self.mul(x[0],y[1]),self.mul(x[1],y[0])))
 def extsign(self,x,d):
  sa,sb=self.sign(x[0]),self.sign(x[1])
  if not sb:return sa
  if not sa or sa==sb:return sb
  return sa*self.sign(sub(self.mul(x[0],x[0]),self.mul(d,self.mul(x[1],x[1]))))

def enc(x):return list(map(str,x))
def check(data):
 deps=[BASE,POLAR,CAPS,FIRST,OLDVOL]
 assert data['source_sha256']=={p:sha(p)for p in deps}
 b,p,c,first,oldvol=map(read,deps);f=Field(data['field']['radicand']);zero=f.rational(0)
 assert f.d==(Q(268,3),Q(36))
 M=tuple(tuple(f.base(x)for x in row)for row in b['form'])
 roots=[tuple(map(f.base,x))for x in b['normals']]+[tuple(map(f.base,x['vector']))for x in p['cuts']]
 cov=[f.mv(M,r)for r in roots];norms=[f.dot(r,v)for r,v in zip(roots,cov,strict=True)]
 h=tuple(map(f.base,b['strict_center']));hc=f.mv(M,h);assert f.sign(f.dot(h,hc))<0
 w=tuple(map(f.base,oldvol['certificate']['vector']));wc=f.mv(M,w)
 assert f.sign(f.dot(w,wc))>0 and f.sign(f.dot(w,hc))<0
 r=tuple(map(f.parse,data['normal']));assert len(r)==21
 assert all(x[2:]==(0,0)for x in r[:-1])and r[-1]==(0,0,-1,0)
 rc=f.mv(M,r);q=f.dot(r,rc);assert q==f.rational(2)==f.parse(data['norm'])
 center=f.dot(r,hc);cut=f.dot(r,wc)
 assert f.sign(center)<0 and center==f.parse(data['strict_center_pairing'])
 assert f.sign(cut)>0 and cut==f.parse(data['previous_escape_pairing'])
 records=[]
 for i,(v,n)in enumerate(zip(cov,norms,strict=True)):
  pp=f.dot(r,v);assert f.sign(pp)<=0
  if pp==zero:kind='2'
  elif f.mul(pp,pp)==f.mul(f.rational(Q(1,4)),f.mul(q,n)):kind='3'
  else:assert f.sign(sub(f.mul(pp,pp),f.mul(q,n)))>0;kind='dotted'
  saved=data['old_pair_checks'][i]
  assert saved['old']==i and saved['kind']==kind and pp==f.parse(saved['pairing'])
  records.append({'old':i,'kind':kind,'pairing':enc(pp)})
 assert records[21]['kind']=='3'
 extras=[]
 for i,cap in enumerate(c['caps']+[{'field':first['field'],'vector':first['normal']}],82):
  d=f.base(cap['field']['radicand']);assert f.sign(d)>0
  other=tuple((f.base(x[0]),f.base(x[1]))for x in cap['vector'])
  pp=(f.dot(rc,tuple(x[0]for x in other)),f.dot(rc,tuple(x[1]for x in other)))
  norm=(f.base(cap['norm']),zero)if i<85 else(f.rational(76),f.rational(-24))
  square=f.extmul(pp,pp,d);delta=tuple(sub(x,f.mul(q,y))for x,y in zip(square,norm,strict=True))
  assert f.extsign(pp,d)<0 and f.extsign(delta,d)>0
  saved=data['old_pair_checks'][i]
  assert saved['old']==i and saved['kind']=='dotted'and pp==tuple(map(f.parse,saved['pairing']))
  records.append({'old':i,'kind':'dotted','pairing':[enc(x)for x in pp]});extras.append((d,other))
 assert len(records)==len(data['old_pair_checks'])==86
 counts=dict(sorted(Counter(x['kind']for x in records).items()))
 assert counts==data['pair_counts']=={'2':18,'3':1,'dotted':67}
 cusp_pairs=[f.dot(r,f.mv(M,tuple(map(f.base,c['ray']))))for c in b['cusp_rows']]
 assert len(cusp_pairs)==107 and all(f.sign(x)<0 for x in cusp_pairs)
 assert cusp_pairs==list(map(f.parse,data['inherited_cusp_strict_pairings']))
 volume=read(VOLUME);assert all(sha(path)==digest for path,digest in volume['source_hashes'].items())
 assert volume['source_hashes'][SOURCE]==sha(SOURCE)
 escape=tuple(map(f.base,volume['certificate']['vector']));ec=f.mv(M,escape)
 eq=f.dot(escape,ec);eh=f.dot(escape,hc)
 assert eq==f.base(volume['certificate']['norm'])and eq==f.rational(Q(73993,8388608))and f.sign(eq)>0
 assert eh==f.base(volume['certificate']['center_pairing'])and eh==f.rational(Q(-28919,2048))
 assert all(f.sign(f.dot(escape,v))<0 for v in cov)and f.sign(f.dot(escape,rc))<0
 for d,other in extras:
  pp=(f.dot(ec,tuple(x[0]for x in other)),f.dot(ec,tuple(x[1]for x in other)))
  assert f.extsign(pp,d)<0
 return {'status':'PASS_INDEPENDENT_87_WALL_H20_CONSTRUCTION_AND_ESCAPE','source':SOURCE,'source_sha256':sha(SOURCE),
         'dependencies':data['source_sha256'],'new_pair_counts':counts,'new_pairs':records,
         'new_norm':enc(q),'strict_center_pairing':enc(center),'previous_escape_cut_pairing':enc(cut),
         'inherited_cusps_strictly_preserved':107,
         'remaining_escape':{'source':VOLUME,'source_sha256':sha(VOLUME),'norm':enc(eq),'center_pairing':enc(eh),'strict_inequalities':87},
         'scope':'Separate Fraction Qsqrt6 and parameterized quadratic arithmetic checks all86 new pairs,norm2,center,previous escape cut and107cusp signs. A new strict future spacelike escape proves infinite volume. Existing86-wall certificates are hash-anchored dependencies.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=read(SOURCE);out=check(data)
 for field in('normal','norm','previous_escape_pairing'):
  bad=deepcopy(data)
  if field=='normal':bad[field][-1][1][0]='1'
  else:bad[field][0][0]=str(Q(bad[field][0][0])+1)
  try:check(bad)
  except AssertionError:pass
  else:raise AssertionError('forged '+field+' accepted')
 out['forged_normal_norm_cut_rejected']=True
 if a.verify:assert read(OUT)==out
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['new_pair_counts'])
