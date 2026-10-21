#!/usr/bin/env python3
"""Independent cusp-star review: exact geometry and outward interval signs.

Reconstructs all25 projections. For open germs, a separate interval proof
checks the opposite null ray directly, without repeating the producer's
first-contact maximization. Imports no producer arithmetic.
"""
import argparse,json,itertools
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_e7_replacement_escape_20261014_review as R
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_repaired_cusp_walk_20261014'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
PREV='results/fresh_h20_e7_finite_edgewalk_20261014_exact.json'
PRIMARY=f'results/{TAG}.json';IDEALS=f'results/{TAG}_ideals.json'
RAD=[1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for s,h in d.get('sources',{}).items():assert sha(s)==h,s
 return d
def run():
 d,old,census,ideals=map(read,[BASE,PREV,PRIMARY,IDEALS])
 def base(a):return{r:F(x)for r,x in zip(RAD,a)if F(x)}
 D=base(d['field']['discriminant']);assert R.sign(D)>0
 zero=({},{});one=(R.c(1),{})
 def num(q):return(R.c(q),{})
 def add(a,b):return(R.add(a[0],b[0]),R.add(a[1],b[1]))
 def scale(a,q):return(R.scale(a[0],q),R.scale(a[1],q))
 def mul(a,b):return(R.add(R.mul(a[0],b[0]),R.mul(D,R.mul(a[1],b[1]))),R.add(R.mul(a[0],b[1]),R.mul(a[1],b[0])))
 def sign(a):
  sa,sb=R.sign(a[0]),R.sign(a[1])
  if not sb:return sa
  if not sa or sa==sb:return sb
  return sa*R.sign(R.sub(R.mul(a[0],a[0]),R.mul(D,R.mul(a[1],a[1]))))
 def total(xs):
  z=zero
  for x in xs:z=add(z,x)
  return z
 def parse(a):return(base(a[0]),base(a[1]))
 def vec(a):return tuple(map(parse,a))
 def bvec(a):return tuple((base(x),{})for x in a)
 Q=[list(map(F,row))for row in d['form']]
 def cov(v):return tuple(total(scale(v[j],q)for j,q in enumerate(row)if q)for row in Q)
 def rawdot(v,w):return total(mul(a,b)for a,b in zip(v,w))
 def dot(v,w):return rawdot(v,cov(w))
 def lin(vs,cs):return tuple(total(mul(v[k],c)for v,c in zip(vs,cs))for k in range(21))
 def rat(x):assert not x[1];return R.rational(x[0])
 N=list(map(vec,d['normals']));CN=list(map(cov,N));H=vec(d['strict_center']);CH=cov(H)
 nu=bvec(d['completed_cusp']['point']);assert nu==vec(census['cusp_point'])
 Hn=[rawdot(H,n)for n in CN];nun=[rawdot(nu,n)for n in CN];H2=dot(H,H);Hnu=dot(H,nu)
 assert dot(nu,nu)==zero and sign(H2)<0 and sign(Hnu)<0
 active=[i for i,g in enumerate(nun)if g==zero];assert active==d['completed_cusp']['incident']
 assert all(sign(g)<=0 for g in nun)and all(sign(g)<0 for g in Hn)
 sq2=({2:F(1)},{})
 roots={i:lin([N[i]],[sq2 if i in[17,79]else one])for i in active}
 components=d['completed_cusp']['components'];choices=[]
 assert len(components)==len(census['component_projections'])==6
 assert sorted(sum((c['walls']for c in components),[]))==active
 for component,saved_choices in zip(components,census['component_projections']):
  ids=component['walls'];assert [r['omitted']for r in saved_choices]==ids
  full=[[rat(dot(roots[i],roots[j]))for j in ids]for i in ids]
  kernel=list(map(F,component['positive_kernel']))
  assert len(kernel)==len(ids)and all(x>0 for x in kernel)
  assert all(sum(a*b for a,b in zip(row,kernel))==0 for row in full)
  rows=[]
  for rec in saved_choices:
   retained=[i for i in ids if i!=rec['omitted']];assert retained==rec['retained']
   M=[[rat(dot(roots[i],roots[j]))for j in retained]for i in retained]
   R.pd(M);I=R.inverse(M);h=[dot(H,roots[i])for i in retained]
   coefficients=[total(scale(h[j],q)for j,q in enumerate(row))for row in I]
   P=lin([roots[i]for i in retained],coefficients)
   assert P==vec(rec['projection'])and dot(P,nu)==zero
   assert all(dot(P,roots[i])==dot(H,roots[i])for i in retained)
   p2=dot(P,P);pairings=[rawdot(P,n)for n in CN]
   assert p2==parse(rec['norm'])and pairings==list(map(parse,rec['pairings']))
   rows.append(dict(omitted=rec['omitted'],retained=retained,projection=P,norm=p2,pairings=pairings))
  choices.append(rows)
 assert [len(a)for a in choices]==census['component_sizes']and sum(map(len,choices))==25
 assert len(choices)==6
 for i,j in itertools.combinations(range(6),2):
  assert all(dot(roots[a],roots[b])==zero for a in components[i]['walls']for b in components[j]['walls'])
 print('review:25 physical projections rebuilt',flush=True)
 # Exact 128-bit outward enclosures; no binary floating point is used.
 def ia(a,b):return(a[0]+b[0],a[1]+b[1])
 def ine(a):return(-a[1],-a[0])
 def im(a,b):
  v=[x*y for x in a for y in b];return(min(v),max(v))
 def isum(xs):
  z=(F(0),F(0))
  for x in xs:z=ia(z,x)
  return z
 dl,du=R.bounds(D,128);assert dl>0;db=R.sqrt_interval(dl,du,128)
 def interval(a):return ia(R.bounds(a[0],128),im(R.bounds(a[1],128),db))
 hnb=list(map(interval,Hn));nub=list(map(interval,nun));h2b=interval(H2);hub=interval(Hnu)
 for rows in choices:
  for r in rows:r['norm_bounds']=interval(r['norm']);r['pair_bounds']=list(map(interval,r['pairings']))
 pointmap={};known=old['preserved_old_vertices']+old['new_cap_finite_vertices']+d['new_edge_endpoints']
 known_sets={tuple(sorted(r['incident'])):i for i,r in enumerate(known)}
 for kind,key in [('finite','finite_vertices'),('ideal_contact','ideal_contact_points')]:
  for k,rec in enumerate(census[key]):
   X=vec(rec['point']);ids=tuple(rec['incident']);q=dot(X,X)
   assert q==parse(rec['norm'])and sign(q)==(-1 if kind=='finite'else 0)
   assert sign(rawdot(X,CH))<0 and sign(dot(X,nu))<0
   assert all(sign(rawdot(X,n))==0 if i in ids else sign(rawdot(X,n))<0 for i,n in enumerate(CN))
   assert len(ids)==20 and len(set(ids)&set(active))==19
   if kind=='finite':assert rec['known_vertex']==known_sets.get(ids)
   assert ids not in pointmap
   pointmap[ids]=(kind,k,X,rec['germs'])
 counts=Counter();seen=set();fallbacks=0;interval_signs=0
 for idx,(selected,rec)in enumerate(zip(itertools.product(*choices),census['edge_germs'])):
  support=sorted(sum((r['retained']for r in selected),[]));omitted=[r['omitted']for r in selected]
  assert len(support)==19 and rec['support']==support and rec['omitted']==omitted
  assert tuple(omitted)not in seen;seen.add(tuple(omitted))
  kind=rec['type'];counts[kind]+=1
  if kind!='open_opposite_end':
   ids=tuple(rec['endpoint_key']);saved_kind,k,X,germs=pointmap[ids]
   assert saved_kind==kind and idx in germs and set(ids)&set(active)==set(support)
   assert rec['first_contact']in ids and rec['first_contact']not in support
   assert rec['contact_walls']==[i for i in ids if i not in support]
   assert all(sign(rawdot(X,CN[i]))<0 or sign(nun[i])<0 for i in range(95)if i not in support)
   continue
  assert rec['opposite_null_incident']==support
  w2b=ia(h2b,ine(isum(r['norm_bounds']for r in selected)));assert w2b[1]<0
  # Z=(-2<H,nu>)W+<W,W>nu is future null, orthogonal to support.
  # Its strict non-support signs prove this edge meets no wall before
  # the opposite light-cone point, independently of the stored winner.
  minus2hub=im((F(-2),F(-2)),hub)
  for j in range(95):
   if j in support:continue
   wnb=ia(hnb[j],ine(isum(r['pair_bounds'][j]for r in selected)))
   zb=ia(im(minus2hub,wnb),im(w2b,nub[j]));interval_signs+=1
   if zb[1]>=0:
    fallbacks+=1
    wn=add(Hn[j],scale(total(r['pairings'][j]for r in selected),-1))
    w2=add(H2,scale(total(r['norm']for r in selected),-1))
    assert sign(add(mul(scale(Hnu,-2),wn),mul(w2,nun[j])))<0
  if (idx+1)%256==0:print('review',idx+1,dict(counts),flush=True)
 assert len(seen)==1728==len(census['edge_germs'])and dict(counts)==census['counts']
 assert counts==dict(finite=26,ideal_contact=16,open_opposite_end=1686)
 assert len(census['finite_vertices'])==26 and sum(r['known_vertex']is None for r in census['finite_vertices'])==20
 # Every finite/ideal record's germ list must be exactly the inverse map.
 for ids,(kind,k,X,germs)in pointmap.items():
  assert germs==[i for i,r in enumerate(census['edge_germs'])if tuple(r.get('endpoint_key',[]))==ids]
 # Rank19 alone does not imply a bounded cusp. Check explicit recession
 # directions at all16 ideal contacts, independently of their classification.
 for k,rec in enumerate(ideals['records']):
  assert rec['ideal_point']==k and rec['bounded_link']is False
  saved=census['ideal_contact_points'][k];X=vec(saved['point']);ids=saved['incident']
  assert ids==rec['incident'];recession=rec['recession'];z=vec(recession['point'])
  assert dot(z,X)==zero and dot(z,z)==parse(recession['norm'])and sign(dot(z,z))>0
  g=[rawdot(z,CN[i])for i in ids]
  assert g==list(map(parse,recession['incident_pairings']))and all(sign(a)<=0 for a in g)
 assert len(ideals['records'])==16
 paths=[BASE,PREV,PRIMARY,IDEALS,'scripts/fresh_h20_e7_replacement_escape_20261014_review.py',f'scripts/{TAG}_review.py']
 return dict(status='PASS_INDEPENDENT_COMPLETE_CUSP_STAR_AND_IDEAL_RECESSIONS',sources={p:sha(p)for p in paths},
  checked_component_projections=25,complete_germs=1728,counts=dict(counts),distinct_finite_vertices=26,
  finite_vertices_new_to_known96=20,verified_finite_total_at_least=116,unbounded_ideal_contact_links=16,
  open_null_strict_interval_signs=interval_signs,interval_bits=128,exact_sign_fallbacks=fallbacks,
  scope='Separate sparse arithmetic reconstructs25component projections and all26finite/16ideal endpoint geometries. All1728supports independently checked. For1686open ends, outward rational interval bounds independently certify all76non-support inequalities at the opposite future null ray, with exact fallback if needed; stored first-contact winners beyond the light cone are not independently maximized. All16ideal contacts have independently checked positive-norm horospherical recession vectors. This is one complete cusp star, not a whole95-wall vertex/edge census or finite-volume example.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'],d['exact_sign_fallbacks'],flush=True)
