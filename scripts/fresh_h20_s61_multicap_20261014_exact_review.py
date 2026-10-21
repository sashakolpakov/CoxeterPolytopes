#!/usr/bin/env python3
"""Independent physical-field and complex-power replay of the S61 right chain."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from collections import Counter
from functools import lru_cache
import sympy as sp
from fresh_h20_bridge_edge_closure_20261014_review import Field,parse
import fresh_h20_bridge_faces_20261014_review as B
ROOT=Path(__file__).resolve().parents[1];BASE='fresh_h20_s61_multicap_20261014_exact';TAG=BASE+'_review'
BITS=512;DEN=1<<BITS
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json','results/fresh_h20_e7_bridge_face_20261014_exact.json','results/fresh_h20_bridge_edge_closure_20261014.json','results/fresh_h20_polar_pair_20261014.json','results/fresh_h20_coupled_roles_20261014_exact.json',f'results/{BASE}.json']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
def out(a):return F(floor(a[0]*DEN),DEN),F(ceil(a[1]*DEN),DEN)
def pt(x):return F(x),F(x)
def plus(a,b):return out((a[0]+b[0],a[1]+b[1]))
def neg(a):return -a[1],-a[0]
def minus(a,b):return plus(a,neg(b))
def times(a,b):
 v=[x*y for x in a for y in b];return out((min(v),max(v)))
def div(a,b):
 assert b[0]>0 or b[1]<0
 return times(a,(1/b[1],1/b[0]))
def root(a):
 assert a[0]>0
 p=isqrt(floor(a[0]*DEN*DEN));q=isqrt(floor(a[1]*DEN*DEN))+1
 r=F(p,DEN),F(q,DEN);assert r[0]**2<=a[0] and r[1]**2>a[1];return r
def total(xs):
 z=pt(0)
 for x in xs:z=plus(z,x)
 return z
def inner(a,b):return total(times(x,y)for x,y in zip(a,b))
def lor(a,b):return plus(neg(times(a[0],b[0])),inner(a[1:],b[1:]))
def overlap(a,b):return a[0]<=b[1]and b[0]<=a[1]
def cmul(x,y):return minus(times(x[0],y[0]),times(x[1],y[1])),plus(times(x[0],y[1]),times(x[1],y[0]))
def cadd(x,y):return plus(x[0],y[0]),plus(x[1],y[1])
def cpow(x,n):
 a=(pt(1),pt(0))
 while n:
  if n&1:a=cmul(a,x)
  x=cmul(x,x);n//=2
 return a
def chain(a):
 aa=times(a,a);z=root(plus(aa,pt(1)));den=plus(aa,pt(2))
 C=div(aa,den);S=div(times(pt(2),z),den);rot=C,neg(S);initial=neg(z),pt(1)
 half=[(a,*cmul(initial,cpow(rot,j)))for j in range(13)]
 return list(reversed(half))+[(x,y,neg(z))for x,y,z in half],(C,S,z)
def endpoint(a,k,rw,r2):
 co,(C,S,z)=chain(a);aa,b,c=co[0]
 f=minus(plus(times(aa,k),div(b,r2)),times(c,rw))
 den=plus(times(a,a),pt(2));den2=times(den,den)
 Cp=div(times(pt(4),a),den2)
 Sp=neg(div(times(pt(2),times(times(a,a),a)),times(z,den2)))
 rot=C,neg(S);init=neg(z),pt(1);initp=neg(div(a,z)),pt(0)
 term1=cmul(initp,cpow(rot,12));term2=cmul(cmul(init,cpow(rot,11)),(times(pt(12),Cp),neg(times(pt(12),Sp))))
 bp,cp=cadd(term1,term2)
 fp=minus(plus(k,div(bp,r2)),times(cp,rw))
 return f,fp

def run():
 old,bridge,edge,polar,cur,d=[json.loads((ROOT/p).read_text())for p in FILES]
 for doc in [old,bridge,edge,polar,cur,d]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 assert d['status']=='PASS_133_WALL_CHAMBER_WITH_COMPACT_S61_29_GON'
 ff=Field([parse(old['field']['discriminant'],0),parse(bridge['field']['new_positive_sqrt_radical'],1),parse(edge['field']['new_positive_sqrt_radical'],2)])
 add,mul,sc,ft,c,sg=ff.add,ff.mul,ff.scale,ff.total,ff.const,ff.sign;zero=ff.zeros[3]
 Q=[list(map(F,row))for row in old['form']];rows=[[(j,q)for j,q in enumerate(row)if q]for row in Q];cache={}
 def vp(raw,k):return tuple(ff.lift(parse(x,k),k)for x in raw)
 def dot(v,w):
  if id(w)not in cache:cache[id(w)]=(w,tuple(ft(sc(w[j],q)for j,q in row)for row in rows))
  return ft(mul(x,y)for x,y in zip(v,cache[id(w)][1]))
 def lin(vs,cs):return tuple(ft(mul(a,v[j])for a,v in zip(cs,vs))for j in range(21))
 def rat(x):
  for k in [3,2,1]:assert x[1]==ff.zeros[k-1];x=x[0]
  assert not any(x[1:]);return x[0]
 @lru_cache(None)
 def iv(x):return out(ff.interval(x,3,BITS))
 N=[vp(v,1)for v in old['normals']]+[vp(bridge['new_normal'],2)]+[vp(v,3)for v in edge['new_normals']]+[tuple(c(x)for x in v)for v in polar['new_normals']]+[vp(cur['cap_R']['normal'],1)]
 nq=[c(x)for x in old['norms']]+[c(2)]+[parse(edge['new_normal_norm'],3)]*2+[c(x)for x in polar['new_norms']]+[ff.lift(parse(cur['cap_R']['norm'],1),1)]
 assert len(N)==len(nq)==101
 rec=edge['recipe'];S=rec['support'];assert S==d['support']==cur['S61_limitation']['support']
 roots=[N[i]for i in S];G=[[rat(dot(v,w))for w in roots]for v in roots];B.pd(G)
 assert G==[list(map(F,row))for row in d['support_gram']]
 frame=[vp(rec['X'],2),N[95],vp(rec['W'],2)];X,rr,W=frame
 Av=ff.lift(parse(rec['A'],2),2);wv=ff.lift(parse(rec['W_norm'],2),2)
 assert [[dot(v,w)for w in frame]for v in frame]==[[sc(Av,-1),zero,zero],[zero,c(2),zero],[zero,zero,wv]]
 assert sg(Av)>0 and sg(add(wv,c(-2)))>0 and all(dot(v,n)==zero for v in frame for n in roots)
 om=vp(rec['dual_weight'],2);assert dot(om,om)==c(F(2,3)) and all(dot(om,v)==zero for v in frame)
 part=lin([om,rr],[c(-1),c(F(-1,2))]);assert part==vp(rec['particular'],2)
 delta=ff.lift(parse(edge['field']['new_positive_sqrt_radical'],2),2)
 assert delta==mul(Av,add(wv,c(F(-5,6))))
 # Coefficientwise tie to the raw old side roots, before new normalization.
 for raw,sgn in zip(edge['new_normals'],[1,-1]):
  vv=[parse(x,3)for x in raw]
  low=lin([part,W],[Av,sc(Av,sgn)])
  assert tuple(x[0]for x in vv)==tuple(x[0]for x in low) and tuple(x[1]for x in vv)==tuple(x[0]for x in X)
 rootsq=[root(iv(v))for v in [Av,c(2),wv]]
 P=[]
 for i,(n,q)in enumerate(zip(N,nq)):
  vals=[dot(v,n)for v in frame]
  row=[div(iv(x),times(a,root(div(iv(q),pt(2)))))for x,a in zip(vals,rootsq)]
  if i in S:assert vals==[zero]*3;row=[pt(0)]*3
  P.append(row)
 hb=[vp(cur['hex_basis'][key],1)for key in ['C','U','V']];hq=[ff.lift(parse(x,1),1)for x in cur['hex_basis']['squared_lengths']]
 for co in cur['hexagon']['frame_coefficients']:
  co=[ff.lift(parse(x,0),0)for x in co]
  P.append([total(div(iv(mul(b,dot(v,u))),times(den,root(iv(q))))for b,u,q in zip(co,hb,hq))for v,den in zip(frame,rootsq)])
 assert len(P)==107
 rw=rootsq[2];r2=rootsq[1];kk=root(minus(iv(wv),pt(F(5,6))))
 expected={95:[pt(0),r2,pt(0)],96:[neg(kk),neg(div(pt(1),r2)),rw],97:[neg(kk),neg(div(pt(1),r2)),neg(rw)]}
 for i,p in expected.items():assert all(overlap(a,b)for a,b in zip(P[i],p));P[i]=p
 assert all(overlap(a,tuple(map(F,b)))for p,row in zip(P,d['frame']['old107_unit_normal_pair_bounds'])for a,b in zip(p,row))
 print('PASS_INDEPENDENT_PHYSICAL_FRAME_AND_107_COVECTORS',flush=True)
 lo,hi=map(F,d['parameter']['isolating_interval']);assert hi>lo>0
 fl,_=endpoint(pt(lo),kk,rw,r2);fh,_=endpoint(pt(hi),kk,rw,r2);_,fp=endpoint((lo,hi),kk,rw,r2)
 assert fl[1]<0<fh[0] and fp[0]>0
 coeff,(Crot,Srot,z)=chain((lo,hi));assert len(coeff)==26
 for v,raw in zip(coeff,d['chain']['frame_coefficient_bounds']):assert all(overlap(x,tuple(map(F,y)))for x,y in zip(v,raw))
 # Exact symbolic identities justify every equality, rather than interval zero overlap.
 a,zz=sp.symbols('a z');CC=a*a/(a*a+2);SS=2*zz/(a*a+2)
 assert sp.simplify((CC*CC+SS*SS-1)*(a*a+2)**2-4*(-a*a+zz*zz-1))==0
 assert sp.expand(-a*a+zz*zz-1).subs(zz**2,a*a+1)==0
 assert sp.simplify((a*a+2)*CC-a*a)==0
 counts=Counter();neighbors=[]
 for j,v in enumerate(coeff):
  fn=[]
  for i,p in enumerate(P):
   g=inner(v,p);right=i in S or(j==0 and i==96)or(j==25 and i==97)
   if right:assert g[0]<=0<=g[1];lab='right';fn.append(i)
   else:assert g[1]<-2;lab='dotted'
   assert d['chain']['old107_pairs'][j][i]['relation']==lab;counts[lab]+=1
  for k,u in enumerate(coeff):
   if k==j:continue
   g=lor(v,u);right=abs(j-k)==1
   if right:assert g[0]<=0<=g[1];fn.append(107+k)
   else:assert g[1]<-2
   if k<j:counts['right'if right else'strict_dotted']+=1
  neighbors.append(sorted(fn))
 counts['dotted']+=counts.pop('strict_dotted',0)
 assert dict(counts)==d['new_pair_counts'] and sum(counts.values())==3107
 assert neighbors==[sorted(x)for x in d['chain']['finite_neighbors']]
 print('PASS_ALL3107_NEW_PAIRS_AND_ROOT_ISOLATION',flush=True)
 face=[[pt(0),r2,pt(0)],[kk,neg(div(pt(1),r2)),rw]]+coeff+[[kk,neg(div(pt(1),r2)),neg(rw)]]
 ids=[95,96]+list(range(107,133))+[97]
 vertices=[];bad=0
 for i in range(29):
  for j in range(i+1,29):
   u,v=face[i],face[j];det=minus(times(u[1],v[2]),times(u[2],v[1]));assert det[0]>0 or det[1]<0
   x=div(minus(times(u[0],v[2]),times(u[2],v[0])),det);y=div(minus(times(u[1],v[0]),times(u[0],v[1])),det);point=[pt(1),x,y]
   vals=[lor(point,n)for n in face];consecutive=j==i+1 or(i==0 and j==28)
   if not consecutive:assert any(g[0]>0 for g in vals);bad+=1;continue
   assert lor(point,point)[1]<0 and all(g[1]<0 for k,g in enumerate(vals)if k not in[i,j])
   own={ids[i],ids[j]}
   for k,p in enumerate(P):
    g=inner(point,p)
    if k in S or k in own:assert g[0]<=0<=g[1]
    else:assert g[1]<0
   vertices.append(sorted(own))
 assert len(vertices)==29 and bad==377 and sorted(vertices)==sorted(v['walls']for v in d['face']['vertices'])
 assert sum(bool(set(v)&set(range(107,133)))for v in vertices)==27
 # Spatial normals of central pair and r95 span and have a positive zero-sum.
 assert coeff[12][1][1]<0 and coeff[12][2][0]>0 and coeff[13][2][1]<0
 for axis in [1,2]:assert total([coeff[12][axis],coeff[13][axis],times(times(r2,z),face[0][axis])])[0]<=0<=total([coeff[12][axis],coeff[13][axis],times(times(r2,z),face[0][axis])])[1]
 Hold=vp(cur['strict_center']['point'],1);qH=dot(Hold,Hold);assert sg(qH)<0;hden=root(iv(sc(qH,-1)))
 hf=[div(iv(dot(v,Hold)),times(den,hden))for v,den in zip(frame,rootsq)]
 center=[pt(1),pt(F(-1,2)),pt(0)];eps=F(d['strict_center']['epsilon']);assert eps>0 and inner(center,hf)[1]<0
 assert plus(plus(pt(F(-3,4)),times(pt(2*eps),inner(center,hf))),pt(-eps*eps))[1]<0
 for v in coeff:assert plus(lor(v,center),times(pt(eps),inner(v,hf)))[1]<0
 for i,p in enumerate(P):
  if i in S:assert p==[pt(0)]*3
  else:assert inner(center,p)[1]<0
 # Inherited strict Hold supplies the remaining inward support directions.
 assert all(sg(dot(Hold,n))<0 for n in N)
 for raw in [old['completed_cusp']['point']]+[x['point']for x in old['preserved_old_bounded_cusps']]:
  nu=tuple(ff.lift(tuple(map(F,x)),0)for x in raw);assert dot(nu,nu)==zero
  pp=[div(iv(dot(v,nu)),den)for v,den in zip(frame,rootsq)]
  assert all(inner(v,pp)[1]<0 for v in coeff)
 E=tuple(c(x)for x in cur['strict_escape']['rational_point']);qE=rat(dot(E,E));assert qE==F(38017800175279,500000000000)
 ep=[div(iv(dot(v,E)),den)for v,den in zip(frame,rootsq)]
 assert all(inner(v,ep)[1]<0 for v in coeff)
 assert plus(inner(center,ep),times(pt(eps),div(iv(dot(E,Hold)),hden)))[1]<0
 S75=cur['hex_basis']['support'];assert set(S75)-set(S)=={75} and nq[75]==c(2)
 J75=[list(map(F,row))for row in cur['hex_basis']['support_inverse']];assert J75[S75.index(75)][S75.index(75)]==5 and all(x>=0 for row in J75 for x in row)
 assert all(minus(pt(2),times(pt(5),times(inner(v,P[75]),inner(v,P[75]))))[1]<0 for v in coeff)
 print('PASS_29_GON_27_NEW_VERTICES_CENTER_CUSPS_ESCAPE_AND_S75',flush=True)
 paths=FILES+[f'scripts/{BASE}.py','scripts/fresh_h20_coupled_roles_20261014_field.py',f'scripts/{TAG}.py','scripts/fresh_h20_bridge_edge_closure_20261014_review.py','scripts/fresh_h20_bridge_faces_20261014_review.py','results/fresh_h20_coupled_roles_20261014_exact_review.json']
 return dict(status='PASS_INDEPENDENT_133_WALL_RIGHT_CHAIN_AND_COMPACT_S61_FACE',sources={p:sha(p)for p in paths},interval_bits=BITS,parameter_interval=list(map(str,[lo,hi])),endpoint_signs=[[str(x)for x in fl],[str(x)for x in fh]],derivative_bound=list(map(str,fp)),new_pairs=dict(counts),all_pair_count=3107,polygon_vertices=29,nonconsecutive_pair_exclusions=377,complete_new_finite_vertices=27,preserved_bounded_cusps=3,strict_escape_norm=str(qE),S75_redundancy=True,implementation='Prior independent sparse-radical tower and new dyadic intervals; scalar chain evaluated using complex powers by repeated squaring, derivative by a closed power rule. No producer arithmetic imports. Old107 geometry inherited through its frozen independent certificate; all3107 new pairs and406 polygon intersections replayed.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
