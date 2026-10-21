#!/usr/bin/env python3
"""Independent recursive-radical replay of the three coupled H20 caps.

Uses the separate Fraction Q(sqrt6) implementation, never either producer.
The cross-pair computations reverse the producer's radical-adjunction order.
"""
import argparse
from collections import Counter
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from fresh_dimension_bridge_20261013_verify import E,Z,O,parse,vec,dot,mv,inertia

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json'
OUT='results/fresh_dimension_bridge_two_angle_20261013_independent.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

class Base:
 zero=Z;one=O
 @staticmethod
 def parse(x):return parse(x)
 @staticmethod
 def encode(x):return x.enc()
 @staticmethod
 def add(x,y):return x+y
 @staticmethod
 def neg(x):return -x
 @staticmethod
 def sub(x,y):return x-y
 @staticmethod
 def mul(x,y):return x*y
 @staticmethod
 def sign(x):return x.sign()
 @staticmethod
 def rational(x):return E(Q(x))
B=Base()

class Extension:
 def __init__(self,base,d):
  assert base.sign(d)>0
  self.base=base;self.d=d;self.zero=(base.zero,base.zero);self.one=(base.one,base.zero)
 def lift(self,x):return x,self.base.zero
 def rational(self,x):return self.lift(self.base.rational(x))
 def parse(self,x):return tuple(self.base.parse(t)for t in x)
 def encode(self,x):return list(map(self.base.encode,x))
 def add(self,x,y):return self.base.add(x[0],y[0]),self.base.add(x[1],y[1])
 def neg(self,x):return self.base.neg(x[0]),self.base.neg(x[1])
 def sub(self,x,y):return self.add(x,self.neg(y))
 def mul(self,x,y):
  b=self.base
  return b.add(b.mul(x[0],y[0]),b.mul(self.d,b.mul(x[1],y[1]))),b.add(b.mul(x[0],y[1]),b.mul(x[1],y[0]))
 def sign(self,x):
  b=self.base;sa,sb=b.sign(x[0]),b.sign(x[1])
  if not sb:return sa
  if not sa:return sb
  if sa==sb:return sa
  gap=b.sub(b.mul(x[0],x[0]),b.mul(self.d,b.mul(x[1],x[1])))
  return sa*b.sign(gap)

def fdot(f,a,b):
 out=f.zero
 for x,y in zip(a,b,strict=True):out=f.add(out,f.mul(x,y))
 return out
def fmv(f,M,x):return tuple(fdot(f,row,x)for row in M)
def same(f,a,b):return f.sign(f.sub(a,b))==0
def label(f,p,q,r):
 assert f.sign(p)<=0
 if f.sign(p)==0:return'2'
 z=f.mul(p,p);qr=f.mul(q,r);delta=f.sign(f.sub(z,qr))
 if delta>0:return'dotted'
 if delta==0:return'parallel'
 for m,c in((3,Q(1,4)),(4,Q(1,2)),(6,Q(3,4))):
  if same(f,z,f.mul(f.rational(c),qr)):return str(m)
 raise AssertionError('forbidden angle')

def rank(matrix):
 a=[list(r)for r in matrix];count=0
 for j in range(len(a[0])):
  k=next((i for i in range(count,len(a))if a[i][j].sign()),None)
  if k is None:continue
  a[count],a[k]=a[k],a[count];p=a[count][j];a[count]=[x/p for x in a[count]]
  for i in range(count+1,len(a)):
   p=a[i][j];a[i]=[x-p*y for x,y in zip(a[i],a[count],strict=True)]
  count+=1
  if count==len(a):break
 return count

def check(data):
 assert data['source']=='data/enumerations/fresh_dimension_bridge_two_angle_20261013_scout.json'
 assert sha(data['source'])==data['source_sha256']
 assert data['base_source']=='data/enumerations/fresh_dimension_bridge_polar_20261013.json'
 assert sha(data['base_source'])==data['base_sha256']
 scout=read(data['source']);polar=read(data['base_source']);assert scout['base_source']==data['base_source']
 assert sha(polar['source'])==polar['source_sha256'];base=read(polar['source'])
 assert sha(base['source'])==base['source_sha256']
 form=tuple(map(vec,base['form']));assert len(form)==21 and all(len(r)==21 for r in form)
 assert all(x.b==0 for row in form for x in row)and inertia([[x.a for x in row]for row in form])==(20,1,0)
 roots=list(map(vec,base['normals']))+[vec(r['vector'])for r in polar['cuts']];assert len(roots)==82
 cov=[mv(form,r)for r in roots];norms=[dot(r,c)for r,c in zip(roots,cov,strict=True)]
 assert all(q.sign()>0 for q in norms)
 h=vec(base['strict_center']);hc=mv(form,h);assert dot(h,hc).sign()<0 and all(dot(r,hc).sign()<0 for r in roots)
 escape=vec(polar['strict_escape']['vector']);ec=mv(form,escape);eq=dot(escape,ec);eh=dot(escape,hc)
 assert eq.sign()>0 and eh.sign()<0 and all(dot(r,ec).sign()<0 for r in roots)
 assert data['strict_escape_from_base']==polar['strict_escape']
 assert eq==parse(polar['strict_escape']['norm'])and eh==parse(polar['strict_escape']['center_pairing'])
 assert len(data['caps'])==3 and data['wall_count']==85
 rows=[];live=[];counts=Counter()
 for cap,expected in zip(data['caps'],(0,580,914),strict=True):
  assert cap['scout_index']==expected
  proposal=scout['candidates'][expected];assert proposal['labels']==[6,6]
  assert cap['target_index']==proposal['target_index']
  target=polar['rejected_polar_proposals'][cap['target_index']]
  ids=[i for i in target['active_indices']if i not in proposal['released']]+proposal['anchors']
  assert ids==cap['active_indices']and len(ids)==len(set(ids))==20 and all(0<=i<82 for i in ids)
  assert rank([cov[i]for i in ids])==20
  field=cap['field'];assert field['kind']=='quadratic'and field['base']=={'kind':'sqrt6'}
  f=Extension(B,parse(field['radicand']));r=tuple(map(f.parse,cap['vector']));assert len(r)==21
  Fform=tuple(tuple(f.lift(x)for x in row)for row in form);rc=fmv(f,Fform,r);q=fdot(f,r,rc)
  assert same(f,q,f.lift(parse(cap['norm'])))and same(f,q,f.rational(Q(2,27)))
  pairs=[fdot(f,r,tuple(map(f.lift,c)))for c in cov]
  kinds=[label(f,p,q,f.lift(n))for p,n in zip(pairs,norms,strict=True)];assert kinds==cap['old_pair_kinds']
  prescribed=list(map(parse,cap['prescribed_pairings']));assert len(prescribed)==20 and prescribed[:18]==[Z]*18
  assert all(same(f,pairs[j],f.lift(v))for j,v in zip(ids,prescribed,strict=True))
  assert all(kinds[i]=='6'for i in ids[-2:])and all(kinds[i]=='2'for i in ids[:18])
  center=fdot(f,r,tuple(map(f.lift,hc)));survives=fdot(f,r,tuple(map(f.lift,ec)))
  assert f.sign(center)<0 and f.sign(survives)<0
  assert same(f,center,f.parse(cap['center_pairing']))and same(f,survives,f.parse(cap['escape_pairing']))
  t=vec(target['vector']);tc=mv(form,t);tq=dot(t,tc);th=dot(t,hc)
  assert tq.sign()>0 and th.sign()<0 and all(dot(t,c).sign()<=0 for c in cov)
  assert len(target['active_indices'])==20 and rank([cov[i]for i in target['active_indices']])==20
  assert tq==parse(target['norm'])and all(dot(t,cov[i])==Z for i in target['active_indices'])
  cut=fdot(f,r,tuple(map(f.lift,tc)));assert f.sign(cut)>0 and same(f,cut,f.parse(cap['target_pairing']))
  counts.update(kinds);live.append((f,r,rc,q))
  rows.append({'scout_index':expected,'norm':cap['norm'],'old_pair_count':len(kinds),'old_pair_types':dict(sorted(Counter(kinds).items())),
               'prescribed_rank':20,'prescribed_right_angles':18,'prescribed_pi6_angles':2,'center_pairing':f.encode(center),
               'cut_target_index':cap['target_index'],'target_norm':tq.enc(),'target_old_wall_feasible':True,
               'target_pairing':f.encode(cut),'escape_pairing':f.encode(survives)})
 cross=[]
 for i,(fi,ri,rci,qi)in enumerate(live):
  for j,(fj,rj,rcj,qj)in enumerate(live[:i]):
   # Reverse the producer: inner sqrt(D_j), outer sqrt(D_i).
   f=Extension(fj,fj.lift(fi.d))
   ri2=tuple((fj.lift(x[0]),fj.lift(x[1]))for x in ri);cj2=tuple(f.lift(x)for x in rcj)
   pp=fdot(f,ri2,cj2);qi2=f.lift(fj.lift(parse(data['caps'][i]['norm'])));qj2=f.lift(qj)
   assert label(f,pp,qi2,qj2)=='dotted'
   saved=next(row for row in data['new_cross_pairs']if row['indices']==[j,i]);assert saved['kind']=='dotted'
   s=saved['pairing'];swapped=((parse(s[0][0]),parse(s[1][0])),(parse(s[0][1]),parse(s[1][1])))
   assert same(f,pp,swapped)
   cross.append({'indices':[j,i],'kind':'dotted','pairing_reversed_tower':f.encode(pp)})
 assert len(cross)==len(data['new_cross_pairs'])==3
 return {'status':'PASS_INDEPENDENT_85_WALL_H20_TWO_ANGLE_CAPS','source':SOURCE,'source_sha256':sha(SOURCE),
         'source_chain_sha256':{p:sha(p)for p in(data['source'],data['base_source'],polar['source'],base['source'])},
         'wall_count':85,'new_caps':3,'new_old_pairs':246,'new_cross_pairs':3,'new_old_pair_type_counts':dict(sorted(counts.items())),
         'caps':rows,'cross_pairs':cross,'strict_escape_norm':eq.enc(),'strict_escape_center_pairing':eh.enc(),
         'scope':'Separate Fraction Qsqrt6 and recursive-radical arithmetic checks all3 new normals,249 new pairs,20 prescribed equations per cap, ranks, strict common center, feasible spacelike targets and surviving strict escape. Cross pairs use reversed radical order. The inherited82-wall angle/cusp certificate is an anchored dependency, not rerun here. All85 facets are actual by acute projection; volume remains infinite.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=read(SOURCE);out=check(data)
 for label0,field in(('forged normal','vector'),('forged norm','norm')):
  bad=deepcopy(data)
  if field=='vector':bad['caps'][0][field][0][0][0]=str(Q(bad['caps'][0][field][0][0][0])+1)
  else:bad['caps'][0][field]=['0','0']
  try:check(bad)
  except AssertionError:pass
  else:raise AssertionError(label0+' accepted')
 bad=deepcopy(data);bad['caps'][0]['target_pairing']=[['0','0'],['0','0']]
 try:check(bad)
 except AssertionError:pass
 else:raise AssertionError('forged target pairing accepted')
 out['forged_normal_norm_and_target_pairing_rejected']=True
 if a.verify:assert read(OUT)==out
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['new_old_pair_type_counts'])
