#!/usr/bin/env python3
"""Exact coupled PD18 repairs of the frozen 60-wall H20 lift.

Rational form; each new root lies in one explicitly positive quadratic
extension of Q. Strict cross-field inequalities use integer outward bounds.
The numerical scout supplies dyadic parameters and a proposed escape only.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from math import isqrt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_coupled_lift_20261014'
SOURCE='results/fresh_dimension_construction_routes_20261014.json';SCOUT=f'results/{TAG}_scout.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def rref(A):
 A=[list(r)for r in A];p=[];row=0
 for j in range(len(A[0])):
  k=next((i for i in range(row,len(A))if A[i][j]),None)
  if k is None:continue
  A[row],A[k]=A[k],A[row];a=A[row][j];A[row]=[x/a for x in A[row]]
  for i in range(len(A)):
   if i!=row and A[i][j]:
    a=A[i][j];A[i]=[x-a*y for x,y in zip(A[i],A[row])]
  p.append(j);row+=1
 return A,p

def inv(A):
 n=len(A);R,p=rref([r+[F(i==j)for j in range(n)]for i,r in enumerate(A)]);assert p[:n]==list(range(n));return[r[n:]for r in R]
def mv(A,x):return[sum(a*b for a,b in zip(r,x))for r in A]
def dot(x,y,Q):return sum(a*b for a,b in zip(x,mv(Q,y)))
def lin(x,y,c):return[a+c*b for a,b in zip(x,y)]
def pd(A):
 A=[r[:]for r in A];ps=[]
 for i in range(len(A)):
  p=A[i][i];assert p>0;ps.append(p)
  for j in range(i+1,len(A)):
   for k in range(i+1,len(A)):A[j][k]-=A[j][i]*A[i][k]/p
 return ps

def sign1(a,b,d):
 assert d>0
 if not b:return (a>0)-(a<0)
 if not a:return (b>0)-(b<0)
 if a*b>0:return (a>0)-(a<0)
 z=a*a-b*b*d;return ((a>0)-(a<0))*((z>0)-(z<0))
def decode(v):return [F(x['a'])for x in v],[F(x['b'])for x in v]
def root(a,b=None,d=F(2)):
 return dict(a=a,b=b if b is not None else[F(0)]*len(a),d=d)
def qp(r,p,Q):return dot(r['a'],p,Q),dot(r['b'],p,Q)
def val_sign(r,p,Q):a,b=qp(r,p,Q);return sign1(a,b,r['d'])
def norm(r,Q):
 assert dot(r['a'],r['b'],Q)==0
 return dot(r['a'],r['a'],Q)+r['d']*dot(r['b'],r['b'],Q)
def zero_pair(r,s,Q):
 aa=dot(r['a'],s['a'],Q);ba=dot(r['b'],s['a'],Q);ab=dot(r['a'],s['b'],Q);bb=dot(r['b'],s['b'],Q)
 if r['d']==s['d']:return aa+r['d']*bb==0 and ba+ab==0
 return aa==ba==ab==bb==0
SCALE=1<<100
def root_bound(d):
 n=isqrt(d.numerator*SCALE*SCALE//d.denominator)
 return n,n+(n*n*d.denominator!=d.numerator*SCALE*SCALE)
def coeff_bounds(r):
 lo,hi=root_bound(r['d']);res=[]
 for a,b in zip(r['a'],r['b']):
  v,w=(lo,hi)if b>=0 else(hi,lo)
  l=a*SCALE+b*v;u=a*SCALE+b*w
  res.append((l.numerator//l.denominator,-((-u.numerator)//u.denominator)))
 return res

def ipair(X,Y):
 lo=hi=0
 for (a,b),(c,d)in zip(X,Y):
  z=(a*c,a*d,b*c,b*d);lo+=min(z);hi+=max(z)
 return lo,hi

def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return{str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return[enc(x)for x in v]
 return v

def run():
 src=json.loads((ROOT/SOURCE).read_text());sc=json.loads((ROOT/SCOUT).read_text())
 for p,h in src['sources'].items():assert sha(p)==h,p
 assert src['status']=='PASS_TWO_AXIS_H18_LIFT_THREE_BOUNDED_LINKS_AND_STRICT_H20_ESCAPE'
 assert all(F(x['b'])==0 for row in src['form']for x in row)
 Q=[[F(x['a'])for x in row]for row in src['form']]
 old=[root(*decode(v))for v in src['normals']];assert len(old)==60
 ha,hb=decode(src['strict_center']);H=[a+b for a,b in zip(ha,hb)]
 v,vb=decode(src['strict_escape']['boundary_point']);assert not any(vb)
 assert dot(H,H,Q)==-4296 and all(val_sign(r,H,Q)<0 for r in old)
 assert sign1(dot(H,ha,Q),dot(H,hb,Q),F(2))<0
 assert dot(v,v,Q)==25
 Eold=root(*decode(src['strict_escape']['point']))
 # Every old source root direction in the triangle has a retained long -1 contact.
 raw=json.loads((ROOT/'data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json').read_text())
 labels=raw['labels'];gram=raw['gram'];tri=[labels.index(x)for x in ['L18','L25','L26']]
 short=[i for i in range(40)if gram[i][i]=={'a':'1','b':'0'}]
 long=[i for i in range(40)if gram[i][i]=={'a':'2','b':'0'}]
 assert len(short)==14 and len(long)==26 and set(tri)<=set(short)
 def g(i,j):return F(gram[i][j]['a']),F(gram[i][j]['b'])
 rawchecks=[]
 for i in tri:
  witness=next(j for j in long if g(i,j)==(-1,0));rawchecks.append([labels[i],labels[witness]])
 for i,j in combinations(tri,2):
  assert g(i,j)==(-1,0)
  for k in range(40):
   a,b=g(i,k);c,d=g(j,k);s=sign1(a+c,b+d,F(2));assert s<=0
   if k in short and k not in (i,j):assert s<0
 for i in short:
  p=F(-1)if i in tri else F(-2)
  assert any(g(i,j)==(p,0)for j in long)
  # New norm3 gives rational 2cos(2theta)=-4/3 or2/3, nonintegral.
  assert (4*p*p/F(6)-2).denominator!=1
 templates=[];base=[];zero=[F(0)]*21
 for rec in sc['recipes']:
  S=rec['support'];assert len(S)==18 and len(set(S))==18
  assert all(not any(old[i]['b'])and qp(old[i],v,Q)==(0,0)for i in S)
  B=[[dot(old[i]['a'],old[j]['a'],Q)for j in S]for i in S];pivots=pd(B)
  coeff=mv(inv(B),[dot(old[i]['a'],H,Q)for i in S]);h0=H[:]
  for c,i in zip(coeff,S):h0=lin(h0,old[i]['a'],-c)
  assert dot(h0,h0,Q)<0 and all(dot(h0,old[i]['a'],Q)==0 for i in S)
  R,p=rref([mv(Q,old[i]['a'])for i in S]+[mv(Q,v),mv(Q,h0)])
  assert len(p)==20;free=next(i for i in range(21)if i not in p);w=zero[:];w[free]=F(1)
  for row,i in zip(R,p):w[i]=-row[free]
  assert dot(w,w,Q)>0 and dot(w,v,Q)==dot(w,h0,Q)==0
  assert all(dot(w,old[i]['a'],Q)==0 for i in S)
  t=F(rec['dyadic_t']);assert t>0;m=lin(v,h0,t);q=dot(m,m,Q);assert q>0;d=q/dot(w,w,Q)
  pair=[root(m,w,d),root(m,[-x for x in w],d)]
  assert norm(pair[0],Q)==norm(pair[1],Q)==2*q and zero_pair(*pair,Q)
  assert dot(m,H,Q)<0
  point=lin(h0,m,-dot(h0,m,Q)/q);assert dot(point,point,Q)<0 and dot(point,H,Q)<0
  for r in pair:assert qp(r,point,Q)==(0,0)
  templates.append(dict(support=S,t=t,h0=h0,w=w,m=m,radicand=d,cap_norm=2*q,finite_vertex=point,pd18_pivots=pivots))
  base.extend(pair)
 normals=old[:];groups=[];vertices=[];target_cuts=[];permutations=[]
 oldkeys={(tuple(r['a']),tuple(r['b'])):i for i,r in enumerate(old)};assert len(oldkeys)==60
 for signs in product((1,-1),repeat=2):
  def flip(x):
   x=x[:];x[19]*=signs[0];x[20]*=signs[1];return x
  perm=[oldkeys[(tuple(flip(r['a'])),tuple(flip(r['b'])))]for r in old];assert sorted(perm)==list(range(60));permutations.append(perm)
  first=len(normals)
  for r in base:normals.append(root(flip(r['a']),flip(r['b']),r['d']))
  for j,rec in enumerate(templates):
   point=flip(rec['finite_vertex']);ids=[perm[i]for i in rec['support']]+[first+2*j,first+2*j+1]
   vertices.append(dict(incident=ids,point=point,norm=dot(point,point,Q),signs=list(signs),template=j))
  # Exact contraction with old irrational E; evaluate by integer bounds below.
  target_cuts.append(dict(a=flip(Eold['a']),b=flip(Eold['b']),d=F(2),owners=list(range(first,first+8))))
  groups.append(dict(signs=list(signs),indices=list(range(first,first+8))))
 assert len(normals)==92 and len(vertices)==16
 norms=[norm(r,Q)for r in normals];assert all(q>0 for q in norms)
 boxes=[coeff_bounds(r)for r in normals]
 cov=[coeff_bounds(root(mv(Q,r['a']),mv(Q,r['b']),r['d']))for r in normals]
 counts=Counter();allright={i:[]for i in range(60,92)}
 for i in range(60,92):
  for j in range(i):
   if zero_pair(normals[i],normals[j],Q):
    counts['right']+=1;allright[i].append(j)
    if j>=60:allright[j].append(i)
   else:
    lo,hi=ipair(boxes[i],cov[j]);prod=norms[i]*norms[j]
    assert hi<0 and hi*hi*prod.denominator>prod.numerator*SCALE**4,(i,j)
    counts['dotted']+=1
 assert sum(counts.values())==32*60+32*31//2 and counts=={'right':592,'dotted':1824}
 assert all(len(ids)==19 for ids in allright.values())
 for row in vertices:
  inc=row['incident'];p=row['point'];assert len(inc)==20 and len(set(inc))==20
  zeroids=[i for i,r in enumerate(normals)if val_sign(r,p,Q)==0];assert sorted(inc)==zeroids
  assert all(val_sign(r,p,Q)<0 for i,r in enumerate(normals)if i not in inc)
  assert dot(p,p,Q)<0 and dot(p,H,Q)<0
 assert len({tuple(sorted(r['incident']))for r in vertices})==16
 assert all(val_sign(r,H,Q)<0 for r in normals)
 cuts=[]
 for rec in target_cuts:
  covtarget=coeff_bounds(root(mv(Q,rec['a']),mv(Q,rec['b']),rec['d']))
  vals=[]
  for i in rec['owners']:
   lo,hi=ipair(boxes[i],covtarget);assert lo>0;vals.append([i,str(F(lo,SCALE*SCALE)),str(F(hi,SCALE*SCALE))])
  cuts.append(vals)
 cusps=[]
 for row in src['repaired_cusps']:
  da,db=decode(row['delta']);assert not any(db)and dot(da,da,Q)==0
  signs=[val_sign(r,da,Q)for r in normals];assert all(s<=0 for s in signs)
  cusps.append(dict(pair=row['pair'],new_signs=signs[60:],previous_bounded_link_retained=True))
 # The LP point is only a proposal, rationalized and checked exactly here.
 pp=sc['relaxed_full_bank_escapes'][0]['point'];escape=[F(round(x*10**9),10**9)+c/1000 for x,c in zip(pp,H)]
 eq=dot(escape,escape,Q);assert eq>0 and dot(escape,H,Q)<0
 assert all(val_sign(r,escape,Q)<0 for r in normals)
 # Explicit open continuation along each of the20 edges of each new vertex.
 # Every outside facet is dotted to a retained new cap on each such edge;
 # hence none can contain a second finite or ideal endpoint.
 for row in vertices:
  inc=set(row['incident']);caps=inc.intersection(range(60,92));assert len(caps)==2
  for released in inc:
   remaining=inc-{released};kept=caps-{released};assert kept
   for j in set(range(92))-inc:
    assert any(j not in allright[c]for c in kept)
 paths=[f'scripts/{TAG}.py',SOURCE,SCOUT,f'scripts/{TAG}_scout.py','data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json']
 return enc(dict(status='PASS_92_WALL_COUPLED_LIFT_16_FINITE_VERTICES_WITH_STRICT_ESCAPE',sources={p:sha(p)for p in paths},
  form=Q,strict_center=H,center_norm=dot(H,H,Q),raw_height_direction_gate=dict(triangle_long_minus_one_witnesses=rawchecks,
   scope='All finite planar height lifts of each original short direction with origin in its height convex hull; all26 original long walls unchanged; pairwise acute Coxeter gates. No new projected directions allowed.',
   conclusion='Every such raw lift has a strict spacelike escape; no finite-volume model. Singletonzero or antipodalunit special heights and equal/perpendicular triangle lines force a missing cusp direction.'),
  base_recipes=templates,height_sign_permutations=permutations,orbit_groups=groups,normals=normals,norms=norms,
  new_pair_counts=dict(counts),strict_pair_interval_scale=SCALE,finite_vertices=vertices,
  old_escape_orbit_cuts=cuts,preserved_three_cusp_links=cusps,
  strict_escape=dict(point=escape,norm=eq,center_pairing=dot(escape,H,Q),all92_strict=True),
  new_vertex_frontier=dict(vertices=16,outgoing_germs=320,all_open=True,reason='Every outgoing edge retains a new cap; all92 facets outside its starting20 are dotted to at least one retained cap.'),
  conditional_nonquasi=src['conditional_nonquasi'],
  scope='Exact92facet Coxeter chamber with16 specified actual new finite vertices and allfour height-sign images of previous strictescape cut. Old60certificate/signature/boundedlink proofs inherited by hash. Every new pair right or strictlydotted, so no new noncrystallographic angle. Global exactstrictescape remains: no finitevolume. Newcorner edges remain open; future completion needs normals meeting these caps with changed supports, not only further dotted repairs. No whole corner/cusp census.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'],d['strict_escape']['norm'])
