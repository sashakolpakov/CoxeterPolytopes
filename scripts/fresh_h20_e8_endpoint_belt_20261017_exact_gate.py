#!/usr/bin/env python3
"""Exact fixed-profile endpoint separation, with every ambient158 Coxeter angle.

The tail parameter is the future conic parameter k, at ambient magnitude9/5, below2cos(pi/7).
Discrete labels2..6 and the full tail x>=9/5 exhaust Coxeter values.
"""
import argparse,hashlib,json,math,itertools
from fractions import Fraction as F
import sympy as s
import fresh_h20_e8_exact_cluster_20261017 as P
ROOT=P.ROOT;TAG='fresh_h20_e8_endpoint_belt_20261017_exact_gate';END='results/fresh_h20_e8_partner_endpoints_20261017.json'
pt,add,sub,neg,mul,div,sqrt=P.pt,P.add,P.sub,P.neg,P.mul,P.div,P.sqrt
def scale(v,c):return[mul(c,x)for x in v]
def plus(v,w):return[add(a,b)for a,b in zip(v,w)]
def lin(vs,cs):return[P.total(mul(c,v[i])for c,v in zip(cs,vs))for i in range(len(vs[0]))]
def sx(v):return[P.alg(x)for x in v]
def arctan(q):
 z=F(0)
 for k in range(180):z+=(-1)**k*q**(2*k+1)/F(2*k+1)
 e=q**361/F(361);return z,z+e
PI=sub(scale([arctan(F(1,5))],pt(16))[0],scale([arctan(F(1,239))],pt(4))[0])
def cosine(m):
 x=div(PI,pt(m));xx=mul(x,x);term=pt(1);ans=pt(1)
 for k in range(1,80):
  term=div(mul(neg(xx),term),pt((2*k-1)*(2*k)));ans=add(ans,term)
 er=xx[1]*max(abs(v)for v in term)/F(159*160)
 return ans[0]-er,ans[1]+er
def run():
 C=json.loads((ROOT/P.CUSP).read_text());E=json.loads((ROOT/END).read_text());D=json.loads((ROOT/'results/fresh_h20_e8_exact_cluster_20261017.json').read_text());H=json.loads((ROOT/P.GEOM).read_text())
 ids=C['core_ids'];SG=s.Matrix([[s.sympify(x)for x in row]for row in C['core_gram']]);SI=SG.inv();IG=[[P.alg(x)for x in row]for row in SG.tolist()]
 def gp(v):return P.mv(IG,v)
 def dot(v,w):return P.inner(v,gp(w))
 # Independently evaluate the literal coefficient vectors against the saved
 # exact physical normals. These gates only remove inadmissible discrete roots.
 Q=[list(map(F,row))for row in H['form']];nz=[(i,j,q)for i,row in enumerate(Q)for j,q in enumerate(row)if q]
 N=[[tuple(map(F,z))for z in v]for v in H['normals']]+[[tuple(map(F,z))for z in r['physical_normal']]for r in D['steps']]
 qs=list(map(F,H['norms']))+[F(2)]*4;thresholds=[sqrt(pt(2*q))for q in qs]
 def physical(v,w):return P.total(mul(pt(q),mul(v[i],w[j]))for i,j,q in nz)
 BP=[[physical(N[i],n)for n in N]for i in ids]
 hp=[physical(N[i],[tuple(map(F,z))for z in D['strict_center']])for i in ids]
 families=[[],[],[]];discrete=[];rejections=[];xlabels=[(2,s.Integer(0)),(3,s.Integer(1)),(4,s.sqrt(2)),(5,(1+s.sqrt(5))/2),(6,s.sqrt(3))];phi=s.Rational(9,5);assert cosine(7)[0]>F(9,10)
 for ci,cen in enumerate(E['endpoint_censuses']):
  wi=[ids.index(i)for i in cen['support']];SB=s.eye(21)[:,wi];WJ=s.Matrix([[s.sympify(x)for x in row]for row in cen['inverse']]);ix=ids.index(158);ee=s.eye(21)[:,ix];u=(ee-SB*WJ*SB.T*SG*ee).applyfunc(s.simplify);sigma=s.simplify((u.T*SG*u)[0]);assert sigma.is_positive
  missing=next(i for i in(159,160,161)if i not in cen['edge']);mi=ids.index(missing);V=sx((-SI[:,mi]/s.sqrt(-SI[mi,mi])).applyfunc(s.simplify));unit=sx((u/s.sqrt(sigma)).applyfunc(s.simplify));assert P.inner(V,hp)[1]<0
  for row in cen['profiles']:
   dd=s.zeros(19,1)
   for j,label in row['contacts']:dd[cen['support'].index(j)]=2*s.cos(s.pi/label)
   pp=(-SB*WJ*dd).applyfunc(s.simplify);p=sx(pp);off=s.simplify((pp.T*SG*ee)[0]);h=s.sympify(row['residual_norm']);assert h.is_positive
   root_h=sqrt(P.alg(h));a=scale(plus(V,scale(unit,pt(-1))),div(root_h,pt(2)));b=scale(plus(scale(V,pt(-1)),scale(unit,pt(-1))),div(root_h,pt(2)))
   t=div(P.alg(phi+off),sqrt(P.alg(sigma*h)));assert t[0]>1;k=add(t,sqrt(sub(mul(t,t),pt(1))));aa=scale(a,k);bb=scale(b,div(pt(1),k));nums=[plus(plus(p,aa),bb),plus(p,scale(aa,pt(2))),aa];gn=[gp(v)for v in nums]
   fam=dict(ci=ci,profile=row['id'],nums=nums,gnums=gn,tail_k=k,off=off,h=h,sigma=sigma,a=a,b=b,p=p,discrete=[])
   families[ci].append(fam)
   for label,x in xlabels:
    rad=s.simplify((x+off)**2/sigma-h);ival=P.alg(rad)
    if ival[1]<=0 or rad==0:rejections.append(dict(edge=ci,profile=row['id'],label158=label,reason='no coherent norm completion',radicand=str(rad)));continue
    assert ival[0]>0
    r=plus(plus(p,scale(unit,neg(div(P.alg(x+off),sqrt(P.alg(sigma)))))),scale(V,sqrt(ival)))
    # Coherence must use the future branch; the opposite branch pairs
    # positively with a strict center because all projection coefficients
    # on a spherical support are nonpositive.
    center=P.inner(r,hp) # Recorded orientation aid only; a new center may move.
    fail=None;margin=F(-10**9)
    for j in range(162):
     if j in cen['support']or j==158:continue
     g=P.inner(r,[BP[t][j]for t in range(21)]);test=add(g,thresholds[j])
     if test[0]>0:fail=dict(wall=j,pair=P.save(g),dotted_test=P.save(test));break
     assert test[1]<0,(ci,row['id'],label,j,P.save(test))
     margin=max(margin,test[1])
    if fail:rejections.append(dict(edge=ci,profile=row['id'],label158=label,reason='unprescribed old contact not dotted',witness=fail));continue
    rec=dict(ci=ci,profile=row['id'],label=label,r=r,gr=gp(r));fam['discrete'].append(rec);discrete.append(rec)
   if row['id']%10==0:print('PROFILE',ci,row['id'],flush=True)
  print('EXACT_EDGE',ci,'profiles',len(families[ci]),'discretes',sum(len(f['discrete'])for f in families[ci]),flush=True)
 paircounts={};coeffmax={};gaps=[];exceptions=[];digest=hashlib.sha256()
 def coeffcheck(name,co,record):
  upper=max(z[1]for z in co);paircounts[name]=paircounts.get(name,0)+1;coeffmax[name]=max(coeffmax.get(name,F(-10**99)),upper)
  digest.update(json.dumps([record,[P.save(z)for z in co]],sort_keys=True).encode())
  if upper>=0:exceptions.append(dict(kind=name,record=record,coefficients=[P.save(z)for z in co]))
 for ci,cj in[(0,1),(0,2),(1,2)]:
  for f in families[ci]:
   for g in families[cj]:
    base=[ci,cj,f['profile'],g['profile']]
    co=[add(P.inner(v,w),pt(2 if a<2 and b<2 else 0))for a,v in enumerate(f['nums'])for b,w in enumerate(g['gnums'])];coeffcheck('tail_tail',co,base)
    for r in f['discrete']:
     co=[add(P.inner(r['r'],w),pt(2 if b<2 else 0))for b,w in enumerate(g['gnums'])];coeffcheck('discrete_tail',co,base+[r['label']])
    for r in g['discrete']:
     co=[add(P.inner(v,r['gr']),pt(2 if a<2 else 0))for a,v in enumerate(f['nums'])];coeffcheck('tail_discrete',co,base+[r['label']])
    for r in f['discrete']:
     for z in g['discrete']:
      dotrz=P.inner(r['r'],z['gr']);record=base+[r['label'],z['label']]
      paircounts['discrete_discrete']=paircounts.get('discrete_discrete',0)+1
      if dotrz[1]<-2:continue
      if dotrz[0]>0:gaps.append(dict(record=record,type='obtuse',pair=P.save(dotrz)));continue
      co=scale([dotrz],pt(F(-1,2)))[0]
      assert 0<co[0]<co[1]<1,(record,P.save(co))
      guess=max(2,int(math.pi/math.acos(float(sum(co)/2))))
      low=cosine(guess);high=cosine(guess+1)
      assert low[1]<co[0]and co[1]<high[0],(record,guess,P.save(co),P.save(low),P.save(high))
      gaps.append(dict(record=record,type='strict consecutive finite-label gap',labels=[guess,guess+1],cosine=P.save(co),lower_cosine=P.save(low),upper_cosine=P.save(high)))
  print('PAIR_BLOCK',ci,cj,'exceptions',len(exceptions),'gaps',len(gaps),flush=True)
 assert not exceptions,exceptions[:1]
 paths=[f'scripts/{TAG}.py',END,P.CUSP,P.GEOM,'results/fresh_h20_e8_exact_cluster_20261017.json',f'scripts/{P.TAG}.py']
 return dict(status='PASS_EXACT_ALL158LABEL_FIXED_ENDPOINT_SEPARATION',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},fixed_profile_counts=[len(z)for z in families],discrete_counts=[sum(len(f['discrete'])for f in z)for z in families],rejected_discretes=rejections,tail_threshold='ambient158 magnitude9/5<2cos(pi/7); includes every label>=7, parallel and dotted contacts',pair_counts=paircounts,strict_coefficient_maxima={k:str(v)for k,v in coeffmax.items()},coefficient_stream_sha256=digest.hexdigest(),finite_gaps=gaps,exceptions=exceptions,scope='All159 fixed endpoint profiles. Every unprescribed old162 contact other than158 is required strictly dotted. Any158 Coxeter contact is allowed; labels2..6 plus x>=9/5 exhaust all possibilities. Every cross-edge pair of such endpoint caps is dotted or non-Coxeter. The three unbounded singleton151 endpoint profiles, additional finite old contacts, and different closure combinatorics remain open. No general H3 nonclosure theorem.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
