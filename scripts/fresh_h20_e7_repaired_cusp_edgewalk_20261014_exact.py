#!/usr/bin/env python3
"""All-real escape-cut exclusion for complete finite extensions on four supports."""
import argparse,json,itertools
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from collections import Counter,defaultdict
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_repaired_cusp_edgewalk_20261014_exact'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json';WALK='results/fresh_h20_e7_repaired_cusp_walk_20261014.json';SCOUT='results/fresh_h20_e7_repaired_cusp_edgewalk_20261014_scout.json'
P=[2,79,11,3,5];Z=(F(0),)*32;EZ=(Z,Z)
def sc(q):return(F(q),)+(F(0),)*31
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,q):return tuple(x*q for x in a)
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
 if not any(a)or not any(b):return(F(0),)*len(a)
 if not any(a[1:]):return scale(b,a[0])
 if not any(b[1:]):return scale(a,b[0])
 h=len(a)//2
 if not h:return(a[0]*b[0],)
 x,y,u,v=a[:h],a[h:],b[:h],b[h:];p=P[h.bit_length()-1]
 return add(mul(x,u),scale(mul(y,v),p))+add(mul(x,v),mul(y,u))
def sg(a):
 if not any(a):return 0
 if len(a)==1:return 1 if a[0]>0 else-1
 h=len(a)//2;x,y=a[:h],a[h:];sx,sy=sg(x),sg(y)
 if not sy:return sx
 if not sx or sx==sy:return sy
 t=sg(sub(mul(x,x),scale(mul(y,y),P[h.bit_length()-1])));assert t
 return sx*t

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return[enc(y)for y in x]
 if isinstance(x,dict):return{str(k):enc(y)for k,y in x.items()}
 return x

def run():
 d,w,s=[json.loads((ROOT/p).read_text())for p in[BASE,WALK,SCOUT]]
 for doc in[d,w,s]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 def pad(a):return tuple(map(F,a))+(F(0),)*(32-len(a))
 D=pad(d['field']['discriminant']);assert sg(D)>0
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
 def parse(v):return tuple(tuple(pad(a)for a in x)for x in v)
 N=list(map(parse,d['normals']));H=parse(d['strict_center']);E=parse(d['strict_escape']['point']);nu=parse(w['cusp_point']);Q=[[F(x)for x in row]for row in d['form']];sp=[[(j,q)for j,q in enumerate(row)if q]for row in Q];cache={}
 def dot(v,z):
  if z not in cache:cache[z]=tuple(et(es(z[j],q)for j,q in row)for row in sp)
  return et(em(a,b)for a,b in zip(v,cache[z]))
 def lin(vs,cs):return tuple(et(em(c,v[k])for c,v in zip(cs,vs))for k in range(21))
 one=(sc(1),Z);two=(sc(2),Z)
 sq2=(tuple(F(i==1)for i in range(32)),Z);sq3=(tuple(F(i==8)for i in range(32)),Z);sq5=(tuple(F(i==16)for i in range(32)),Z)
 mag={3:one,4:sq2,5:es(ea(one,sq5),F(1,2)),6:sq3}
 def rat(a):assert a[1]==Z and not any(a[0][1:]);return a[0][0]
 h=dot(H,nu);en=dot(E,nu);assert sign(h)<0 and sign(en)<0
 nn=[dot(nu,n)for n in N];assert all(sign(a)<=0 for a in nn)
 lookup={}
 for row in s['continuous_rejections']:
  if not row['profile'].get('I2_tail_sample'):lookup[row['edge'],tuple(map(tuple,row['profile']['contacts']))]=row
 assert len(lookup)==418
 choices=[{r['omitted']:r for r in cc}for cc in w['component_projections']]
 finite=[];tails=[];edges=[]
 for erec in s['edge_templates']:
  ei=erec['cusp_germ'];edge=w['edge_germs'][ei];ids=edge['support'];assert edge['type']=='open_opposite_end'and ids==erec['support'];parts=[choices[i][j]for i,j in enumerate(edge['omitted'])]
  W=lin([H]+[parse(r['projection'])for r in parts],[one]+[es(one,-1)]*6);w2=dot(W,W);assert sign(w2)<0 and dot(W,nu)==h
  Wn=[dot(W,n)for n in N];WE=dot(W,E);ZE=ea(es(em(h,WE),-2),em(w2,en));assert sign(ZE)<0
  assert all(sign(Wn[i])==0 for i in ids)
  nonshared=[i for i in range(95)if i not in ids]
  Zn={i:ea(es(em(h,Wn[i]),-2),em(w2,nn[i]))for i in nonshared};assert all(sign(a)<0 for a in Zn.values())
  rr=[lin([N[i]],[sq2 if i in[17,79]else one])for i in ids];M=[[rat(dot(a,b))for b in rr]for a in rr];MI=B.inv(M);si=[sq2 if i in[17,79]else one for i in ids]
  dual=[lin(rr,[es(si[i],v)for v in MI[i]])for i in range(19)]
  for i,u in enumerate(dual):
   assert dot(u,nu)==EZ and dot(u,W)==EZ
   assert all(dot(u,N[j])==(one if i==k else EZ)for k,j in enumerate(ids))
  diag=[MI[i][i]*(2 if ids[i]in[17,79]else 1)for i in range(19)];assert min(diag)>=F(1,2)
  comps=[];remaining=set(range(19))
  while remaining:
   cc={min(remaining)};todo=list(cc)
   for i in todo:
    for j in sorted(remaining-cc):
     if M[i][j]:cc.add(j);todo.append(j)
   remaining-=cc;comps.append(sorted(cc))
  which={i:k for k,c in enumerate(comps)for i in c};tail_positions=[c[0]for c in comps if len(c)==1]
  assert [ids[i]for i in tail_positions]==[r['wall']for r in erec['symbolic_I2_tails']]
  profiles=[[]]
  for size in[1,2,3]:
   for inds in itertools.combinations(range(19),size):
    if len({which[i]for i in inds})!=size or sum(diag[i]for i in inds)>=2:continue
    for labels in itertools.product([3,4,5,6],repeat=size)if size<3 else[(3,3,3)]:
     cost=et(es(em(mag[m],mag[m]),diag[i])for i,m in zip(inds,labels))
     if sign(ea(two,es(cost,-1)))>0:profiles.append([[ids[i],m]for i,m in zip(inds,labels)])
  assert profiles==[r['contacts']for r in erec['finite_pool_profiles']]
  # Negative scale2ah converts positive wall and negative escape pairings
  # into respectively negative and positive exact polynomial numerators.
  def numerator(pn,q,a,j):
   vn,wn=(en,WE)if j=='E'else(nn[j],Wn[j])
   return et([es(em(h,pn),2*a),es(em(h,wn),2*a*a),em(ea(q,es(w2,-a*a)),vn)])
  for contacts in profiles:
   row=lookup[ei,tuple(map(tuple,contacts))];a=F(format((row['lower']+row['upper'])/2,'.12g'));assert a>0
   pp=lin([dual[ids.index(i)]for i,m in contacts],[es(mag[m],-1)for i,m in contacts])if contacts else tuple(EZ for _ in range(21))
   cost=dot(pp,pp);q=ea(two,es(cost,-1));assert sign(q)>0 and dot(pp,W)==dot(pp,nu)==EZ
   pe=dot(pp,E);ke=numerator(pe,q,a,'E');assert sign(ke)>0
   j=row['blocking_wall'];kn=numerator(dot(pp,N[j]),q,a,j);assert j in nonshared and sign(kn)<0
   finite.append(dict(edge=ei,contacts=contacts,residual_norm=q,separating_parameter=a,blocking_wall=j,wall_numerator=kn,escape_numerator=ke))
  # A single isolated A1 contact has arbitrary c=2cos(pi/m), m>=7.
  # At fixed a each numerator is a quadratic in c on[√3,2].
  def polynomial(u,a,j):
   vn,wn,un=(en,WE,dot(u,E))if j=='E'else(nn[j],Wn[j],dot(u,N[j]))
   A=es(vn,F(-1,2));BB=es(em(h,un),-2*a);CC=ea(es(em(h,wn),2*a*a),em(ea(two,es(w2,-a*a)),vn))
   return A,BB,CC
  def evalpoly(pol,c):return et([em(pol[0],em(c,c)),em(pol[1],c),pol[2]])
  for pos in tail_positions:
   wall=ids[pos];u=dual[pos];assert dot(u,u)==(sc(F(1,2)),Z)
   sampled=[r for r in s['continuous_rejections']if r['edge']==ei and r['profile'].get('I2_tail_sample')and r['profile']['contacts'][0][0]==wall]
   lo=max(r['upper']for r in sampled);hi=min(r['lower']for r in sampled);assert lo<hi;a=F(format((lo+hi)/2,'.12g'))
   e_pol=polynomial(u,a,'E');A,BB,CC=e_pol;assert sign(A)>0
   left,right=evalpoly(e_pol,sq3),evalpoly(e_pol,two);dl=ea(es(em(A,sq3),2),BB);dr=ea(es(A,4),BB)
   assert sign(left)>0 and sign(right)>0
   if sign(dl)>=0:minimum='left';mincert=left
   elif sign(dr)<=0:minimum='right';mincert=right
   else:minimum='interior';mincert=ea(es(em(A,CC),4),es(em(BB,BB),-1));assert sign(mincert)>0
   found=None
   for j in [r['blocking_wall']for r in sampled]+nonshared:
    if j not in nonshared:continue
    pol=polynomial(u,a,j);assert sign(pol[0])>=0
    l,r=evalpoly(pol,sq3),evalpoly(pol,two)
    if sign(l)<0 and sign(r)<0:found=(j,pol,l,r);break
   assert found is not None,(ei,wall)
   j,pol,l,r=found
   tails.append(dict(edge=ei,contact_wall=wall,labels='all integers m>=7',c_interval=['sqrt3','2'],separating_parameter=a,blocking_wall=j,wall_polynomial=pol,wall_endpoint_values=[l,r],escape_polynomial=e_pol,escape_endpoint_values=[left,right],escape_derivatives=[dl,dr],escape_minimum_case=minimum,escape_minimum_certificate=mincert))
  edges.append(dict(edge=ei,support=ids,omitted=edge['omitted'],finite_profiles=len(profiles),I2_tails=len(tail_positions),h=h,W_norm=w2,escape_null_pairing=ZE))
  print('EXACT_EDGE',ei,len(profiles),len(tail_positions),flush=True)
 assert len(finite)==418 and len(tails)==15
 paths=[BASE,WALK,SCOUT,f'scripts/{TAG}.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py']
 return enc(dict(status='PASS_FOUR_SUPPORT_ALL_FINITE_ENDPOINT_PROFILES_MISS_ESCAPE',sources={p:sha(p)for p in paths},field=dict(base_primes=P,extra_quadratic_discriminant=D),edges=edges,finite_profile_certificates=finite,I2_tail_certificates=tails,counts=dict(edges=4,finite_profiles=418,all_label_I2_tails=15),scope='All norm2 outward caps making a finite endpoint on one of these four cusp edges, with a spherical extension of its19 retained walls and acute pairing with all95, cannot cut the saved95 strict escape. Finite-profile census is complete, including15 all-label I2 tails. Other cusp edges, ideal endpoint profiles, changed walls, and non-cutting finite endpoint additions remain outside this statement. The numerical538-row scout is merely separator provenance; every scalar/sign and completeness claim used here is rechecked exactly.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'])
