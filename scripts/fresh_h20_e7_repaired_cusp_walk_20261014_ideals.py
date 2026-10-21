#!/usr/bin/env python3
"""Exact full-link diagnostic for the16 ideal first contacts of the cusp walk."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from collections import Counter
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
import fresh_h20_e7_finite_edgewalk_20261014_frontier as G
from fresh_h20_e8_interior_caps_20261014_exact_review import pd
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_repaired_cusp_walk_20261014_ideals'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json';WALK='results/fresh_h20_e7_repaired_cusp_walk_20261014.json'
Z=B.Z;add,sub,mul,scale,sg=B.add,B.sub,B.mul,B.scale,B.sign
EZ=(Z,Z)
def sc(q):return(F(q),)+(F(0),)*15
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return[enc(y)for y in x]
 if isinstance(x,dict):return{str(k):enc(y)for k,y in x.items()}
 return x
def run():
 d,w=[json.loads((ROOT/p).read_text())for p in[BASE,WALK]]
 for doc in[d,w]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 D=tuple(map(F,d['field']['discriminant']))
 def ea(a,b):return(add(a[0],b[0]),add(a[1],b[1]))
 def es(a,q):return(scale(a[0],q),scale(a[1],q))
 def em(a,b):
  if a[1]==Z:return(mul(a[0],b[0]),mul(a[0],b[1]))
  if b[1]==Z:return(mul(a[0],b[0]),mul(a[1],b[0]))
  return(add(mul(a[0],b[0]),mul(D,mul(a[1],b[1]))),add(mul(a[0],b[1]),mul(a[1],b[0])))
 def sign(a):
  x,y=a;sx,sy=sg(x),sg(y)
  if not sy:return sx
  if not sx or sx==sy:return sy
  return sx*sg(sub(mul(x,x),mul(D,mul(y,y))))
 def et(xs):
  a=EZ
  for x in xs:a=ea(a,x)
  return a
 def parse(v):return tuple(tuple(tuple(map(F,a))for a in x)for x in v)
 N=list(map(parse,d['normals']));H=parse(d['strict_center']);Q=[[F(x)for x in row]for row in d['form']];sp=[[(j,q)for j,q in enumerate(row)if q]for row in Q];cache={}
 def dot(v,z):
  if z not in cache:cache[z]=tuple(et(es(z[j],q)for j,q in row)for row in sp)
  return et(em(a,b)for a,b in zip(v,cache[z]))
 def lin(vs,cs):return tuple(et(em(c,v[k])for c,v in zip(cs,vs))for k in range(21))
 def rat(a):assert a[1]==Z and not any(a[0][1:]);return a[0][0]
 paircache={}
 def pair(i,j):
  k=tuple(sorted([i,j]))
  if k not in paircache:paircache[k]=dot(N[i],N[j])
  return paircache[k]
 records=[]
 for number,row in enumerate(w['ideal_contact_points']):
  X=parse(row['point']);ids=row['incident'];assert len(ids)==20 and dot(X,X)==EZ and sign(dot(X,H))<0
  assert all(sign(dot(X,n))==0 if i in ids else sign(dot(X,n))<0 for i,n in enumerate(N))
  rem=set(ids);components=[];spherical=[];rank=0
  while rem:
   seed=min(rem);seen={seed};todo=[seed];factors={seed:(sc(1),Z)}
   for i in todo:
    for j in sorted(rem-seen):
     g=pair(i,j)
     if sign(g):
      assert sign(g)<0;seen.add(j);todo.append(j);factors[j]=es(em(factors[i],g),-1)
   rem-=seen;cc=sorted(seen);RR=[lin([N[i]],[factors[i]])for i in cc]
   M=[[rat(em(em(factors[i],factors[j]),pair(i,j)))for j in cc]for i in cc]
   rr=len(G.rref(M)[1]);rank+=rr;record=dict(walls=cc,scales=[factors[i]for i in cc],gram=M,rank=rr)
   if rr==len(cc):
    assert pd(M);record['type']='spherical';spherical.append((cc,RR,M))
   else:
    assert rr==len(cc)-1;kk=G.kernel(M);assert len(kk)==1;k=kk[0]
    if all(a<0 for a in k):k=[-a for a in k]
    assert all(a>0 for a in k)and pd([r[:-1]for r in M[:-1]])
    s=lin(RR,[(sc(a),Z)for a in k]);coord=next(i for i,a in enumerate(X)if sign(a))
    # Cross multiplication proves the positive null sum is on this ray.
    assert all(em(a,X[coord])==em(b,s[coord])for a,b in zip(s,X))and sign(dot(s,H))<0
    record.update(type='affine',positive_kernel=k)
   components.append(record)
  assert rank==19 and spherical
  cc,RR,M=spherical[0];mi=B.inv(M);v=lin(RR,[(sc(-a),Z)for a in mi[0]])
  assert dot(v,v)==(sc(mi[0][0]),Z)and mi[0][0]>0 and dot(v,X)==EZ
  vg=[dot(v,N[i])for i in ids];assert all(sign(a)<=0 for a in vg)and any(sign(a)<0 for a in vg)
  records.append(dict(ideal_point=number,incident=ids,quotient_rank=rank,components=components,bounded_link=False,recession=dict(point=v,norm=dot(v,v),incident_pairings=vg,negative_weight_component=cc)))
 assert len(records)==16
 paths=[BASE,WALK,f'scripts/{TAG}.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py','scripts/fresh_h20_e7_finite_edgewalk_20261014_frontier.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
 return enc(dict(status='PASS_16_IDEAL_CONTACTS_ALL_UNBOUNDED_LINKS',sources={p:sha(p)for p in paths},records=records,scope='Full incident sets for all16 ideal first-contact points from the1728-germ census. Each quotient spans19 but has a spherical component, with an explicit positive-norm tangent recession. None is a proper bounded-link cusp. No new repair walls or whole-chamber census.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
