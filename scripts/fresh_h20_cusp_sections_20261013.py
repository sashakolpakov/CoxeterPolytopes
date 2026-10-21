#!/usr/bin/env python3
"""Exact all-real fork-leaf cusp-section certificates.

Only Python's standard library is used.  The numerical scout is not imported.
The matching report supplies the elementary geometric implications.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from collections import Counter
import argparse, hashlib, json, math

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_cusp_sections_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
OUT=ROOT/f'results/{TAG}.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mv(A,x):return [dot(r,x)for r in A]
def pair(x,y,Q):return dot(x,mv(Q,y))
def add(a,b):return [x+y for x,y in zip(a,b,strict=True)]
def scale(c,a):return [c*x for x in a]
def enc(x):
 if isinstance(x,list):return [enc(y)for y in x]
 if isinstance(x,F):return str(x)
 return x
def inv(A):
 n=len(A);M=[[F(x)for x in r]+[F(i==j)for j in range(n)]for i,r in enumerate(A)]
 for j in range(n):
  k=next(i for i in range(j,n)if M[i][j]);M[j],M[k]=M[k],M[j]
  q=M[j][j];M[j]=[x/q for x in M[j]]
  for i in range(n):
   if i!=j:
    q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[j],strict=True)]
 return [r[n:]for r in M]
def pd(A):
 A=[[F(x)for x in r]for r in A];p=[]
 while A:
  q=A[0][0];assert q>0;p.append(q)
  A=[[A[i][j]-A[i][0]*A[0][j]/q for j in range(1,len(A))]for i in range(1,len(A))]
 return p
def comps(A):
 left=set(range(len(A)));out=[]
 while left:
  q=[min(left)];left.remove(q[0])
  for i in q:
   ns=[j for j in sorted(left)if A[i][j]];q+=ns;left.difference_update(ns)
  out.append(sorted(q))
 return out
def project(h,ids,W,G,Q,C):
 A=[[G[i][j]for j in ids]for i in ids];piv=pd(A);B=inv(A)
 assert all(x>=0 for r in B for x in r)
 cs=mv(B,[dot(h,C[i])for i in ids])
 x=[h[k]-sum(c*W[i][k]for c,i in zip(cs,ids,strict=True))for k in range(22)]
 ps=[dot(x,r)for r in C]
 assert pair(x,x,Q)<0 and pair(x,h,Q)<0
 assert [i for i,p in enumerate(ps)if p==0]==sorted(ids)
 assert all(ps[i]<0 for i in range(210)if i not in ids)
 return x,{'walls':ids,'point':enc(x),'positive_LDL_pivots':enc(piv),'norm':enc(pair(x,x,Q)),
           'center_pairing':enc(pair(x,h,Q)),'maximum_other_pairing':enc(max(ps[i]for i in range(210)if i not in ids))}
def plane_key(a,b):
 M=[[F(x)for x in a],[F(x)for x in b]];row=0
 for j in range(len(a)):
  k=next((i for i in range(row,2)if M[i][j]),None)
  if k is None:continue
  M[row],M[k]=M[k],M[row];q=M[row][j];M[row]=[x/q for x in M[row]]
  for i in range(2):
   if i!=row:
    q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[row],strict=True)]
  row+=1
  if row==2:break
 assert row==2
 return tuple(tuple(r)for r in M)

def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
 Q=[[-int(x)for x in r]for r in raw['basis_gram']]
 W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 assert len(W)==210 and all(len(w)==22 for w in W)
 C=[mv(Q,w)for w in W];G=[[dot(w,c)for c in C]for w in W]
 assert all(G[i][j]<=0 for i,j in combinations(range(210),2))
 total=[sum(w[k]for w in W)for k in range(22)];g=math.gcd(*total);h=[x//g for x in total]
 assert pair(h,h,Q)==-14 and set(map(lambda c:dot(h,c),C))=={-1,-8}
 walk=json.loads((ROOT/WALK).read_text())
 assert walk['status']=='PASS_COMPLETE_CORNER_ORBIT_CENSUS' and not walk['pending_orbit_ids']
 assert walk['sources'][RAW]==sha(RAW)
 cusps=[r for r in walk['corner_orbits']if r['type']=='ideal'];assert len(cusps)==18
 rows=[];planes=set();by_cusp=Counter();types=Counter();face_cache={}
 def face(ids):
  key=tuple(sorted(ids))
  if key not in face_cache:face_cache[key]=project(h,list(key),W,G,Q,C)
  return face_cache[key]
 for cu in cusps:
  u=cu['primitive_ray'];ids=cu['active_walls'];up=[dot(u,c)for c in C]
  assert pair(u,u,Q)==0 and pair(u,h,Q)<0 and max(up)==0
  assert [i for i,p in enumerate(up)if not p]==ids
  components=cu['affine_components'];assert sorted(sum([c['walls']for c in components],[]))==ids
  assert sum(c['rank']for c in components)==20
  for c in components:
   marks=list(map(F,c['positive_kernel']));cs=c['walls'];lam=F(c['ray_multiple']);A=[[G[i][j]for j in cs]for i in cs]
   assert min(marks)>0 and lam>0 and len(comps(A))==1
   assert all(dot(r,marks)==0 for r in A)
   assert [sum(k*W[i][d]for k,i in zip(marks,cs,strict=True))for d in range(22)]==scale(lam,u)
   assert len(cs)-1==c['rank']
  for aa,bb in combinations(ids,2):
   if G[aa][aa]!=G[bb][bb] or G[aa][bb]!=0:continue
   if any(G[aa][i]!=G[bb][i]for i in ids if i not in[aa,bb]):continue
   # Every equal-length orthogonal twin in this literal census is short.
   assert G[aa][aa]==2
   comp=next(c for c in components if aa in c['walls']);cs=comp['walls'];assert bb in cs
   marks=dict(zip(cs,map(F,comp['positive_kernel'])));m=marks[aa]
   assert marks[bb]==m and F(comp['ray_multiple'])==m
   for a,b in[(aa,bb),(bb,aa)]:
    e0=scale(F(1,2),add(W[a],scale(-1,W[b])))
    assert pair(e0,e0,Q)==1 and pair(e0,u,Q)==pair(e0,h,Q)==0
    v=[dot(e0,c)for c in C]
    assert all(v[i]==(1 if i==a else -1 if i==b else 0)for i in ids)
    # e(1/2) is a strictly positive combination of all component walls except b.
    ew={i:(F(1)if i==a else marks[i]/(2*m))for i in cs if i!=b}
    assert min(ew.values())>0
    assert add(e0,scale(F(1,2),u))==[sum(k*W[i][d]for i,k in ew.items())for d in range(22)]
    keep=[i for i in ids if i!=b];H=[[G[i][j]-v[i]*v[j]for j in keep]for i in keep]
    assert all(H[i][i]>0 for i in range(len(keep)))
    ratios=Counter()
    for i,j in combinations(range(len(keep)),2):
     assert H[i][j]<=0
     ratio=F(H[i][j]**2,H[i][i]*H[j][j]);assert ratio in[0,F(1,4),F(1,2),F(3,4),1]
     ratios[str(ratio)]+=1
    link=[]
    for js in comps(H):
     ws=[keep[i]for i in js];cc=next(c for c in components if ws[0]in c['walls'])
     cm=dict(zip(cc['walls'],map(F,cc['positive_kernel'])));ks=[2*cm[i]if i==a else cm[i]for i in ws]
     A=[[H[i][j]for j in js]for i in js];assert all(dot(r,ks)==0 for r in A)
     assert min(ks)>0 and len(comps(A))==1
     # Constant and linear coefficients of sum k_i Pi_t(r_i).
     p0=[sum(k*(W[i][d]-v[i]*e0[d])for k,i in zip(ks,ws,strict=True))for d in range(22)]
     p1=scale(-sum(k*v[i]for k,i in zip(ks,ws,strict=True)),u)
     lam=F(cc['ray_multiple']);slope=-2*m if a in ws else F(0)
     assert p0==scale(lam,u)and p1==scale(slope,u)
     assert lam>0 and lam+slope/2>=0
     link.append({'walls':ws,'positive_kernel':enc(ks),'null_multiple_coefficients':enc([lam,slope]),'rank':len(ws)-1})
    assert sum(z['rank']for z in link)==19
    # Prefer an auxiliary wall attached only to the omitted leaf.
    candidates=[]
    for attached in[b,a]:
     for j in range(210):
      if G[j][j]==2 and G[j][attached]==-1 and all(G[j][i]==0 for i in cs if i!=attached):candidates.append((attached,j))
    assert candidates,(cu['id'],a,b)
    selected=None
    for attached,j in candidates:
     try:zT,certT=face([i for i in cs if i!=b]+[j])
     except AssertionError:continue
     selected=attached,j,zT,certT;break
    assert selected is not None
    attached,j,zT,certT=selected;z0,cert0=face([a,b,j])
    assert dot(z0,mv(Q,e0))==0 and dot(zT,mv(Q,add(e0,scale(F(1,2),u))))==0
    A=-pair(z0,u,Q);B=-pair(zT,u,Q);assert min(A,B)>0
    assert set(cert0['walls'])&set(certT['walls'])=={a,j}
    assert pair(z0,zT,Q)<0
    # e_j=+1/2-t for attached=b, and -1/2-t for attached=a.
    sign=1 if attached==b else -1
    assert v[j]==F(sign,2) and up[j]==-1 and G[a][j]==(0 if sign==1 else -1)
    # Both have off-diagonal -(1/2-t), projected a norm one.
    # Four times the numerator subtracted from denominator is
    # 3/4 + (4+sign)t - 5t^2, concave and positive at 0 and 1/2.
    gap=[F(3,4),F(4+sign),F(-5)]
    assert gap[0]>0 and sum(gap[i]*F(1,2)**i for i in range(3))>0
    assert pair(e0,e0,Q)==1
    rows.append({'cusp_orbit':cu['id'],'retained_leaf':a,'omitted_leaf':b,'cusp':u,'base_normal':enc(e0),
      'parameter_domain':'0<t<1/2; reversed leaf order covers negative t; t=0 has a positive duplicate',
      'endpoint_positive_combination':[[i,str(k)]for i,k in ew.items()],
      'local_link':link,'local_squared_cosine_counts':dict(ratios),'auxiliary_wall':j,'auxiliary_attached_leaf':attached,
      'ridge_walls':[a,j],'left_face':cert0,'right_face':certT,
      'interpolation_weights':{'A':str(A),'B':str(B),'formula':'(1/2-t)*B*z0+t*A*zT'},
      'endpoint_inner_product':enc(pair(z0,zT,Q)),
      'cosine_squared_numerator_coefficients':['1/4','-1','1'],
      'cosine_squared_denominator_coefficients':[str(F(7,4)),str(sign),'-1'],
      'denominator_minus_four_numerator':enc(gap)})
    by_cusp[str(cu['id'])]+=1;types['attached_to_omitted'if sign==1 else'attached_to_retained']+=1
    planes.add(plane_key(e0,u))
 assert len(rows)==118 and len(planes)==59
 return {'status':'PASS_ALL_REAL_FORK_LEAF_CUSP_SECTION_PENCILS_EXCLUDED',
  'sources':{p:sha(p)for p in[RAW,WALK,'scripts/'+Path(__file__).name]},
  'ambient_dimension':21,'section_dimension':20,'ambient_ideal_orbits':18,'cusp_orbits_with_twins':len(by_cusp),
  'oriented_pencil_count':len(rows),'distinct_normal_planes':len(planes),'pencils_by_cusp':dict(by_cusp),
  'auxiliary_attachment_counts':dict(types),'strict_center':h,'pencils':rows,
  'scope':'Every genuine interior H20 section whose normal lies in any of these59 literal two-dimensional normal planes has finite volume and a full Coxeter cusp link, but an actual ridge has cosine squared strictly between0 and1/4. This excludes all real parameters, not other cusp-adapted normal directions, other cusp representatives outside symmetry transport, or arbitrary H20 polytopes.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['oriented_pencil_count'],out['distinct_normal_planes'],out['auxiliary_attachment_counts'])
