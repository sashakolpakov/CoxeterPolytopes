#!/usr/bin/env python3
"""One exact bridge-face cutter; renewed strict escape, no finite-volume claim."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from collections import Counter
from decimal import Decimal, localcontext
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
from fresh_h20_e7_edge_cusp_repair_20261014_exact import inv as binv
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_bridge_face_20261014_exact'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json';WALK='results/fresh_h20_e7_repaired_cusp_walk_20261014.json'
Z=B.Z;O=(Z,Z);add,sub,mul,scale,sg=B.add,B.sub,B.mul,B.scale,B.sign
S=list(range(1,11))+[12,13,14,15,17,31,61,87]
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
 D=tuple(map(F,d['field']['discriminant']));one=(sc(1),Z);two=(sc(2),Z)
 def a(x,y):return(add(x[0],y[0]),add(x[1],y[1]))
 def s(x,q):return(scale(x[0],q),scale(x[1],q))
 def m(x,y):
  if x[1]==Z:return(mul(x[0],y[0]),mul(x[0],y[1]))
  if y[1]==Z:return(mul(x[0],y[0]),mul(x[1],y[0]))
  return(add(mul(x[0],y[0]),mul(D,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0])))
 def sign(x):
  u,v=x;su,sv=sg(u),sg(v)
  if not sv:return su
  if not su or su==sv:return sv
  return su*sg(sub(mul(u,u),mul(D,mul(v,v))))
 def inverse(x):
  u,v=x;den=sub(mul(u,u),mul(D,mul(v,v)));assert sg(den)!=0
  y=binv(den);r=(mul(u,y),scale(mul(v,y),-1));assert m(x,r)==one;return r
 def total(xs):
  y=O
  for x in xs:y=a(y,x)
  return y
 def parse(v):return tuple(tuple(tuple(map(F,b))for b in x)for x in v)
 N=list(map(parse,d['normals']));H=parse(d['strict_center']);E=parse(d['strict_escape']['point']);Q=[[F(x)for x in row]for row in d['form']];sp=[[(j,q)for j,q in enumerate(row)if q]for row in Q];cache={}
 def dot(v,z):
  if z not in cache:cache[z]=tuple(total(s(z[j],q)for j,q in row)for row in sp)
  return total(m(x,y)for x,y in zip(v,cache[z]))
 def lin(vs,cs):return tuple(total(m(c,v[k])for c,v in zip(cs,vs))for k in range(21))
 def rat(x):assert x[1]==Z and not any(x[0][1:]);return x[0][0]
 G=[[rat(dot(N[i],N[j]))for j in S]for i in S];from fresh_h20_e8_interior_caps_20261014_exact_review import pd
 assert pd(G);GI=B.inv(G)
 def project(v):
  hh=[dot(v,N[i])for i in S];cs=[total(s(b,q)for b,q in zip(hh,row))for row in GI]
  return lin([v]+[N[i]for i in S],[one]+[s(x,-1)for x in cs])
 F0=project(H);A0,A1=[project(N[i])for i in[0,28]];G2=[[dot(x,y)for y in[A0,A1]]for x in[A0,A1]];det=a(m(G2[0][0],G2[1][1]),s(m(G2[0][1],G2[0][1]),-1));assert sign(det)<0
 invdet=inverse(det);sq3=(tuple(F(i==8)for i in range(16)),Z);rhs=[s(one,-1),s(sq3,-1)]
 cs=[m(a(m(G2[1][1],rhs[0]),s(m(G2[0][1],rhs[1]),-1)),invdet),m(a(m(G2[0][0],rhs[1]),s(m(G2[0][1],rhs[0]),-1)),invdet)]
 part=lin([A0,A1],cs);assert dot(part,N[0])==rhs[0]and dot(part,N[28])==rhs[1]
 def residual(v):
  z=project(v);t0,t1=dot(z,A0),dot(z,A1);c0=m(a(m(G2[1][1],t0),s(m(G2[0][1],t1),-1)),invdet);c1=m(a(m(G2[0][0],t1),s(m(G2[0][1],t0),-1)),invdet)
  return lin([z,A0,A1],[one,s(c0,-1),s(c1,-1)])
 aux='center';z=residual(H)
 if sign(dot(z,z))==0:
  for k in range(21):
   z=residual(tuple(one if i==k else O for i in range(21)))
   if sign(dot(z,z))>0:aux=k;break
 z2=dot(z,z);p2=dot(part,part);assert sign(z2)>0 and dot(z,part)==O and all(dot(z,N[i])==O for i in S+[0,28])
 K=m(a(two,s(p2,-1)),inverse(z2));assert sign(K)>0
 # New field (old32) plus positive sqrt(K).
 ZZ=(O,O)
 def na(x,y):return(a(x[0],y[0]),a(x[1],y[1]))
 def ns(x,q):return(s(x[0],q),s(x[1],q))
 def nm(x,y):return(a(m(x[0],y[0]),m(K,m(x[1],y[1]))),a(m(x[0],y[1]),m(x[1],y[0])))
 def sn(x):
  u,v=x;su,sv=sign(u),sign(v)
  if not sv:return su
  if not su or su==sv:return sv
  return su*sign(a(m(u,u),s(m(K,m(v,v)),-1)))
 def nt(xs):
  y=ZZ
  for x in xs:y=na(y,x)
  return y
 def lift(v):return tuple((x,O)for x in v)
 nc={}
 def ndot(v,t):
  if t not in nc:nc[t]=tuple(nt(ns(t[j],q)for j,q in row)for row in sp)
  return nt(nm(x,y)for x,y in zip(v,nc[t]))
 def nlin(vs,cs):return tuple(nt(nm(c,v[k])for c,v in zip(cs,vs))for k in range(21))
 NN=list(map(lift,N));HH=lift(H);EE=lift(E);chosen=[]
 for branch in[-1,1]:
  r=tuple((p,s(zz,branch))for p,zz in zip(part,z));pairs=[ndot(r,n)for n in NN]
  if all(sn(x)<=0 for x in pairs)and sn(ndot(r,EE))>0:chosen.append((branch,r,pairs))
 assert len(chosen)==1;branch,r,pairs=chosen[0];assert ndot(r,r)==(two,O)
 counts=Counter();pairrows=[]
 for i,g in enumerate(pairs):
  prod=2*F(d['norms'][i]);gg=nm(g,g)
  if sn(g)==0:kind='right'
  elif gg==((sc(prod/4),Z),O):kind='pi3'
  elif gg==((sc(prod/2),Z),O):kind='pi4'
  elif gg==((sc(3*prod/4),Z),O):kind='pi6'
  elif gg==((sc(prod),Z),O):kind='parallel'
  else:assert sn(na(gg,((sc(-prod),Z),O)))>0;kind='dotted'
  counts[kind]+=1;pairrows.append(dict(wall=i,pairing=g,relation=kind))
 assert [i for i,g in enumerate(pairs)if sn(g)==0]==S
 assert pairs[0]==((sc(-1),Z),O)and pairs[28]==(s(sq3,-1),O)
 NN.append(r)
 # Whole-face point projected to the new geodesic is an actual rank19 edge point.
 facev=lift(F0);fr=ndot(facev,r);edgepoint=nlin([facev,r],[(one,O),ns(fr,F(-1,2))])
 assert sn(ndot(edgepoint,edgepoint))<0 and sn(ndot(edgepoint,HH))<0
 edgepairs=[ndot(edgepoint,n)for n in NN];assert all(sn(g)==0 if i in S+[95]else sn(g)<0 for i,g in enumerate(edgepairs))
 # A strict center follows by an explicit small inward perturbation.
 hr=ndot(HH,r);kk=1
 while sn(na(hr,((sc(-2*kk),Z),O)))>=0:kk*=2
 center=None
 for power in range(1,50):
  eps=F(1,2**power);C=nlin([edgepoint,HH,r],[(one,O),((sc(eps),Z),O),((sc(-eps*kk),Z),O)])
  cp=[ndot(C,n)for n in NN]
  if sn(ndot(C,C))<0 and sn(ndot(C,HH))<0 and all(sn(t)<0 for t in cp):center=C;break
 assert center is not None
 facet=nlin([center,r],[(one,O),ns(ndot(center,r),F(-1,2))]);assert sn(ndot(facet,facet))<0
 assert all(sn(ndot(facet,n))<0 for n in NN[:-1])and ndot(facet,r)==ZZ
 # Renewed strict spacelike escape near one of the previously certified open null ends.
 nu=parse(w['cusp_point']);h=dot(H,nu);components=[{row['omitted']:row for row in cc}for cc in w['component_projections']];escape=None
 sq2=(tuple(F(i==1)for i in range(16)),Z)
 for ei,rec in enumerate(w['edge_germs']):
  if rec['type']!='open_opposite_end':continue
  pp=[parse(components[k][j]['projection'])for k,j in enumerate(rec['omitted'])];W=lin([H]+pp,[one]+[s(one,-1)]*6);ww=dot(W,W);Znull=lin([W,nu],[s(h,-2),ww]);NZ=lift(Znull)
  if sn(ndot(NZ,r))>=0:continue
  ids=rec['support'];RR=[lin([N[i]],[sq2 if i in[17,79]else one])for i in ids];MG=[[rat(dot(x,y))for y in RR]for x in RR];MI=B.inv(MG);scales=[sq2 if i in[17,79]else one for i in ids]
  coef=[total(s(scales[j],-MI[i][j])for j in range(19))for i in range(19)];U=lin(RR,coef)
  assert dot(U,Znull)==O and all(dot(U,N[i])==s(one,-1)for i in ids)
  u2=dot(U,U);zh=dot(Znull,H);assert sign(u2)>0 and sign(zh)<0;t=m(u2,inverse(s(zh,-4)))
  for power2 in range(1,70):
   eps2=F(1,2**power2);X=lin([Znull,U,H],[one,s(one,eps2),s(t,eps2*eps2)]);XX=lift(X);xx=ndot(XX,XX)
   if sn(xx)>0 and sn(ndot(XX,HH))<0 and all(sn(ndot(XX,n))<0 for n in NN):escape=XX;break
  if escape is not None:break
 assert escape is not None
 # Round only as a proposal; exact rational signs/norm certify the saved escape.
 geometric_escape=escape
 with localcontext() as ctx:
  ctx.prec=100
  radicals=[Decimal(q).sqrt() for q in d['field']['base_square_root_basis']]
  def value16(v):return sum((Decimal(q.numerator)/Decimal(q.denominator))*u for q,u in zip(v,radicals))
  sqrtD=value16(D).sqrt()
  def value32(v):return value16(v[0])+sqrtD*value16(v[1])
  sqrtK=value32(K).sqrt()
  vals=[value32(v[0])+sqrtK*value32(v[1]) for v in geometric_escape]
  rational_escape=None
  for denominator in [10**6,10**9,10**12,10**15,10**18]:
   rv=[F(int((v*denominator).to_integral_value()),denominator) for v in vals]
   XX=tuple(((sc(q),Z),O) for q in rv)
   if sn(ndot(XX,XX))>0 and sn(ndot(XX,HH))<0 and all(sn(ndot(XX,n))<0 for n in NN):rational_escape=XX;break
 assert rational_escape is not None
 escape=rational_escape
 cusp_points=[tuple((tuple(map(F,a)),Z) for a in d['completed_cusp']['point'])]+[tuple((tuple(map(F,a)),Z) for a in c['point']) for c in d['preserved_old_bounded_cusps']]
 cusp_pairs=[ndot(lift(v),r) for v in cusp_points];assert all(sn(g)<0 for g in cusp_pairs)
 paths=[BASE,WALK,f'scripts/{TAG}.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py','scripts/fresh_h20_e7_edge_cusp_repair_20261014_exact.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
 return enc(dict(status='PASS_96_WALL_BRIDGE_FACE_CUTTER_STILL_INFINITE',sources={p:sha(p)for p in paths},form=Q,field=dict(base=d['field'],new_positive_sqrt_radical=K,coordinate_encoding='[[base16 a,base16 b],[base16 c,base16 d]]=(a+b sqrt(D))+(c+d sqrt(D))sqrt(K)'),recipe=dict(right_support=S,support_gram=G,anchors=[0,28],anchor_values=rhs,residual_anchor_gram=G2,particular=part,particular_norm=p2,auxiliary=aux,complement=z,complement_norm=z2,radical=K,branch=branch),new_normal=r,new_pairs=pairrows,new_pair_counts=dict(counts),old_escape_pairing=ndot(EE,r),preserved_bounded_cusp_pairings=cusp_pairs,finite_vertex_preservation=dict(known_count_at_least=116,proof='Common strict center and all pair acuteness imply strict retention by projection to each old PD20 support.'),strict_center=center,strict_center_norm=ndot(center,center),center_epsilon=F(1,2**power),center_inward_coefficient=kk,new_facet_point=facet,actual_edge=dict(support=S+[95],point=edgepoint,norm=ndot(edgepoint,edgepoint),all96_pairings=edgepairs),strict_escape=dict(point=escape,norm=ndot(escape,escape),all96_pairings=[ndot(escape,n)for n in NN],source_cusp_germ=ei,epsilon=F(1,2**power2),timelike_coefficient=t,rational_rounding_denominator=denominator),scope='Retains all95 and adds one norm2 cap, right to the18-root bridge support and pi3/pi6 to old0/28. All95 new pair gates, common strict center, facet actuality, actual new rank19 edge, old escape cut and renewed strict escape are exact. No finite vertex or complete finite-area H2 polygon is claimed.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');arg=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if arg.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'])
