#!/usr/bin/env python3
"""Independent exact construction of the five-root cusp repair.

The formulas below use only frozen90 input; no numerical scout, repair
producer, or producer arithmetic is imported.
"""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from collections import Counter
from flint import fmpq,fmpq_mat
import fresh_h20_e7_replacement_escape_20261014_review as R
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_edge_cusp_repair_20261014'
BASE='results/fresh_h20_e7_finite_edgewalk_20261014_exact.json'
FRONT='results/fresh_h20_e7_finite_edgewalk_20261014_frontier.json'
RAD=[1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for s,h in d.get('sources',{}).items():assert sha(s)==h,s
 return d
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(tuple,list)):return[enc(v)for v in x]
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 return x

def run(endpoints=False):
 d,f=map(read,[BASE,FRONT]);Q=[list(map(F,row))for row in d['form']]
 def base(a):return{r:F(x)for r,x in zip(RAD,a)if F(x)}
 D=base(d['field']['discriminant']);assert R.sign(D)>0
 zero=({},{});one=(R.c(1),{});sq2=({2:F(1)},{});sq3=({3:F(1)},{})
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
 def cov(v):return tuple(total(scale(v[j],q)for j,q in enumerate(row)if q)for row in Q)
 def rawdot(v,w):return total(mul(a,b)for a,b in zip(v,w))
 def dot(v,w):return rawdot(v,cov(w))
 def lin(vs,cs):return tuple(total(mul(v[k],c)for v,c in zip(vs,cs))for k in range(21))
 def rat(x):assert not x[1];return R.rational(x[0])
 def encode_field(x):return[[str(a.get(r,F(0)))for r in RAD]for a in x]
 def encode_vector(v):return list(map(encode_field,v))
 N=list(map(vec,d['normals']));nu=bvec(f['new_parallel_cusp']['point']);z=bvec(f['new_parallel_cusp']['lineality_vector'])
 assert dot(nu,nu)==zero and dot(z,nu)==zero and dot(z,z)==num(F(98,27))
 quotient=[lin([N[i]for i in[12,13,18,31,65]],[num(-1)]*5),
           lin([N[15]],[num(-1)]),
           lin([N[17],N[79],N[83]],[scale(sq2,-1),scale(sq2,-1),num(-1)]),
           lin([z],[scale(sq3,F(3,7))]),lin([z],[scale(sq3,F(-3,7))])]
 shifts=[num(1),num(2),num(1),add(one,scale(sq3,F(4,7))),scale(sq3,F(10,7))]
 caps=[lin([a,nu],[one,k])for a,k in zip(quotient,shifts)]
 assert all(dot(a,a)==num(2)and dot(a,nu)==zero for a in quotient)
 assert all(dot(r,r)==num(2)and dot(r,nu)==zero for r in caps)
 N+=caps;CN=list(map(cov,N));norms=list(map(F,d['norms']))+[F(2)]*5
 H=vec(d['strict_center']);CH=cov(H);E=vec(d['strict_escape']['point'])
 assert sign(dot(H,H))<0 and all(sign(rawdot(H,n))<0 for n in CN)
 counts=Counter();records=[]
 for i in range(90,95):
  for j in range(i):
   g=rawdot(N[i],CN[j]);p=norms[i]*norms[j];assert sign(g)<=0
   sq=mul(g,g)
   if g==zero:kind='right'
   elif sq==num(p/4):kind='pi3'
   elif sq==num(p/2):kind='pi4'
   elif sq==num(3*p/4):kind='pi6'
   elif sq==num(p):kind='parallel'
   else:assert sign(add(sq,num(-p)))>0;kind='dotted'
   counts[kind]+=1;records.append(dict(walls=[i,j],relation=kind,pairing=encode_field(g)))
 assert len(records)==460
 for i,n in enumerate(N):
  c=scale(rawdot(H,CN[i]),-1/norms[i]);p=lin([H,n],[one,c])
  assert sign(c)>0 and sign(dot(p,p))<0 and sign(rawdot(p,CH))<0 and rawdot(p,CN[i])==zero
 # Acuteness and strict center imply strictness at every other wall.
 old_vertices=d['preserved_old_vertices']+d['new_cap_finite_vertices']
 for row in old_vertices:
  V=bvec(row['point']);assert all(sign(rawdot(V,CN[i]))<0 for i in range(90,95))
 for row in d['preserved_cusps']:
  V=bvec(row['point']);assert dot(V,V)==zero and all(sign(rawdot(V,CN[i]))<0 for i in range(90,95))
 ep=[rawdot(E,n)for n in CN]
 assert all(sign(x)<0 for x in ep)and sign(dot(E,E))>0 and sign(rawdot(E,CH))<0
 # Complete cusp link: all95 inequalities and every connected Gram factor.
 np=[rawdot(nu,n)for n in CN];assert all(sign(x)<=0 for x in np)
 active=[i for i,x in enumerate(np)if x==zero];assert len(active)==25
 scaled=[lin([N[i]],[sq2 if i in[17,79]else one])for i in active]
 G=[[rat(dot(a,b))for b in scaled]for a in scaled]
 def fm(M):return fmpq_mat([[fmpq(str(x))for x in row]for row in M])
 assert fm(G).rank()==19
 remaining=set(active);components=[]
 while remaining:
  ids={min(remaining)};todo=list(ids)
  for i in todo:
   for j in sorted(remaining-ids):
    if G[active.index(i)][active.index(j)]:ids.add(j);todo.append(j)
  ids=sorted(ids);remaining-=set(ids)
  M=[[G[active.index(i)][active.index(j)]for j in ids]for i in ids];size=len(M)
  assert fm(M).rank()==size-1
  # All but one root form a positive-definite matrix; the full block is
  # then PSD of corank one. Its positive kernel makes the Euclidean
  # intersection bounded in the block's span.
  assert all(fm([row[:k]for row in M[:k]]).det()>0 for k in range(1,size))
  GI=R.inverse([row[:-1]for row in M[:-1]])
  kernel=[-sum(GI[i][j]*M[j][-1]for j in range(size-1))for i in range(size-1)]+[F(1)]
  assert all(x>0 for x in kernel)and all(sum(a*b for a,b in zip(row,kernel))==0 for row in M)
  components.append(dict(walls=ids,rank=size-1,gram=M,positive_kernel=kernel))
 assert sorted(x['rank']for x in components)==[1,1,1,3,5,8]
 assert sum(x['rank']for x in components)==19
 paths=[BASE,FRONT,'scripts/fresh_h20_e7_replacement_escape_20261014_review.py',f'scripts/{TAG}_review.py']
 continuation={}
 if endpoints:
  primary=f'results/{TAG}_exact.json';p=read(primary);paths.append(primary)
  assert list(map(vec,p['normals']))==N
  assert p['form']==d['form'] and p['field']==d['field'] and p['new_pairs']==records
  assert p['new_pair_counts']==dict(counts)
  vertices=d['preserved_old_vertices']+d['new_cap_finite_vertices']
  cache={};incidences=[];wall_counts=Counter();ideal_edges=[]
  for endpoint_index,rec in enumerate(p['new_edge_endpoints']):
   vi,j,k=rec['old_vertex'],rec['released'],rec['new_wall']
   assert 44<=vi<84 and 90<=k<95
   if vi not in cache:
    old=vertices[vi];ids=old['incident'];V=bvec(old['point'])
    rr=[lin([N[i]],[sq2 if i in[17,79]else one])for i in ids]
    GM=[[rat(dot(a,b))for b in rr]for a in rr];GI=R.inverse(GM)
    cache[vi]=(ids,V,rr,GI)
   ids,V,rr,GI=cache[vi];pos=ids.index(j);s=sq2 if j in[17,79]else one
   dual=lin(rr,[scale(s,q)for q in GI[pos]])
   assert dot(V,dual)==zero and dot(dual,N[j])==one
   assert all(dot(dual,N[i])==zero for i in ids if i!=j)
   t=(base(rec['parameter']),{});point=bvec(rec['point']);endpoint_active=tuple(sorted(rec['incident']))
   assert sign(t)>0 and point==lin([V,dual],[one,scale(t,-1)])
   assert set(endpoint_active)==(set(ids)-{j})|{k} and len(endpoint_active)==20
   assert sign(dot(point,point))<0 and sign(rawdot(point,CH))<0
   assert dot(point,point)==(base(rec['norm']),{})
   assert all(sign(rawdot(point,n))==0 if i in endpoint_active else sign(rawdot(point,n))<0 for i,n in enumerate(CN))
   incidences.append(endpoint_active);wall_counts[k]+=1
   common=sorted(set(endpoint_active)&set(active))
   if len(common)==19:
    assert k==93
    assert all(rawdot(point,CN[i])==rawdot(nu,CN[i])==zero for i in common)
    assert all(sign(rawdot(point,CN[i]))<0 or sign(rawdot(nu,CN[i]))<0 for i in range(95)if i not in common)
    ideal_edges.append(dict(finite_endpoint_index=endpoint_index,cusp_support=common))
  assert len(set(incidences))==len(incidences)==p['edge_continuation_scope']['actual_new_endpoints']
  assert len(ideal_edges)==6
  for rec in p['actual_nonincident_ridges']:
   point=bvec(rec['point']);ids=rec['walls'];assert len(set(ids))==2
   assert sign(dot(point,point))<0 and sign(rawdot(point,CH))<0
   assert all(sign(rawdot(point,n))==0 if i in ids else sign(rawdot(point,n))<0 for i,n in enumerate(CN))
  continuation=dict(primary_normals_match=True,actual_new_endpoints=len(incidences),
   endpoint_wall_counts=dict(wall_counts),actual_nonincident_ridges=len(p['actual_nonincident_ridges']),
   finite_to_proper_ideal_edges=ideal_edges,
   scope='All saved accepted endpoints independently verified as distinct future timelike vertices and first intersections along their stated edges. All saved ridge witnesses checked. Six endpoints share19independent facets with the bounded cusp and their open segments are strict at every other wall, certifying six finite-to-proper-ideal edges. The full bounded541-germ search is not independently re-enumerated; no complete new-cap vertex census.')
 return enc(dict(status='PASS_INDEPENDENT_FIVE_ROOT_95_WALL_CUSP_REPAIR_STILL_INFINITE',sources={p:sha(p)for p in paths},
  new_normals=[encode_vector(v)for v in caps],new_pairs=records,new_pair_counts=dict(counts),
  shifts=[encode_field(k)for k in shifts],actual_facets=95,preserved_finite_vertices=len(old_vertices),preserved_bounded_cusps=2,
  repaired_cusp=dict(point=encode_vector(nu),incident=active,rank=19,components=components),
  continuation=continuation,
  strict_escape_norm=str(rat(dot(E,E))),new_escape_pairings=[encode_field(g)for g in ep[90:]],
  scope='Independent exact formulas using only frozen90 data, with separate sparse radical arithmetic. Checks460newpairs,95actualfacets,preservation of84finitevertices and two prior boundedcusps. The new nullray has25incidentfacets and six PSD affine components with positive kernels spanning rank19, hence a bounded full horospherical link. The saved strict spacelike escape survives. No enumeration of newly created finite vertices or global finite-volume claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--endpoints',action='store_true');a=ap.parse_args();d=run(a.endpoints);suffix='_endpoints_review'if a.endpoints else'_review';p=ROOT/f'results/{TAG}{suffix}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'],flush=True)
