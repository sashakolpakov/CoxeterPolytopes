#!/usr/bin/env python3
"""Independent saved-certificate review; no producer arithmetic imports.

Reads the frozen high-core catalogue and the new certificate. Reconstructs
port and subset coverage; checks every rational witness and explicit inertia.
An independent Q(sqrt2) class verifies all26 input core inertias/deletions.
"""
import argparse,json
from collections import Counter,defaultdict,deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from math import isqrt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_distance8_high_disjoint_20261017'
TAG=BASE+'_review'
CAT='fresh_high_rank_quasilanner_census_20261017'

class K:
 def __init__(self,a=0,b=0):self.a,self.b=(a.a,a.b)if isinstance(a,K)else(F(a),F(b))
 def __add__(self,o):o=K(o);return K(self.a+o.a,self.b+o.b)
 __radd__=__add__
 def __neg__(self):return K(-self.a,-self.b)
 def __sub__(self,o):return self+-K(o)
 def __rsub__(self,o):return K(o)+-self
 def __mul__(self,o):o=K(o);return K(self.a*o.a+2*self.b*o.b,self.a*o.b+self.b*o.a)
 __rmul__=__mul__
 def __truediv__(self,o):o=K(o);d=o.a*o.a-2*o.b*o.b;assert d;return self*K(o.a/d,-o.b/d)
 def __bool__(self):return bool(self.a or self.b)
 def sign(self):
  if not self.b:return (self.a>0)-(self.a<0)
  if not self.a:return (self.b>0)-(self.b<0)
  sa=(self.a>0)-(self.a<0);sb=(self.b>0)-(self.b<0)
  if sa==sb:return sa
  d=self.a*self.a-2*self.b*self.b
  return sa*((d>0)-(d<0))

def sign(x):return x.sign()if isinstance(x,K)else (x>0)-(x<0)
def sub(a,ids):return [[a[i][j]for j in ids]for i in ids]
def inertia(a):
 a=[r[:]for r in a];out=[0,0,0]
 while a:
  k=next((i for i in range(len(a))if a[i][i]),None)
  if k is not None:
   ids=[k]+[j for j in range(len(a))if j!=k];a=sub(a,ids);q=a[0][0];out[sign(q)<0]+=1
   a=[[a[i][j]-a[i][0]*a[0][j]/q for j in range(1,len(a))]for i in range(1,len(a))]
  else:
   ij=next(((i,j)for i in range(len(a))for j in range(i+1,len(a))if a[i][j]),None)
   if ij is None:out[2]+=len(a);break
   i,j=ij;a=sub(a,[i,j]+[k for k in range(len(a))if k not in ij]);q=a[0][1];out[0]+=1;out[1]+=1
   a=[[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/q for j in range(2,len(a))]for i in range(2,len(a))]
 return out

def inverse(a):
 n=len(a);x=[[F(v)for v in row]+[F(i==j)for j in range(n)]for i,row in enumerate(a)]
 for j in range(n):
  p=next(k for k in range(j,n)if x[k][j]);x[j],x[p]=x[p],x[j];q=x[j][j];x[j]=[v/q for v in x[j]]
  for i in range(n):
   if i!=j:
    q=x[i][j];x[i]=[v-q*w for v,w in zip(x[i],x[j])]
 return [r[n:]for r in x]

def rroot(q):
 q=F(q);a=isqrt(q.numerator);b=isqrt(q.denominator);assert a*a==q.numerator and b*b==q.denominator
 return F(a,b)

def gauge(g):
 bits={}
 for start in range(len(g)):
  if start in bits:continue
  bits[start]=0;todo=[start]
  for i in todo:
   for j,m in enumerate(g[i]):
    if i==j or not m:continue
    want=bits[i]^(m==4)
    if j in bits:assert bits[j]==want
    else:bits[j]=want;todo.append(j)
 sc=[2**bits[i]for i in range(len(g))]
 return scaled(g,sc),sc

def scaled(g,sc,unit=False):
 assert len(sc)==len(g)and all(v>0 for v in sc)
 return [[F(2*sc[i])if i==j else (-F(1)if unit else-rroot(F(sc[i]*sc[j])*(g[i][j]-2)))if g[i][j]else F(0)for j in range(len(g))]for i in range(len(g))]

def check_matrix(saved,g,sc):
 a=scaled(g,sc);assert [[F(x)for x in row]for row in saved]==a
 return a

def witness(row,g,strict=False):
 ids=row['nodes'];assert len(ids)==len(set(ids));p=sub(g,ids)
 mode=row.get('comparison','positive_diagonal_congruence')
 if mode=='unit_bond_relaxation':
  assert row['gauge_squared']==[1]*len(ids)
  a=[[F(2 if i==j else-int(bool(p[i][j])))for j in range(len(ids))]for i in range(len(ids))]
 else:
  assert mode=='positive_diagonal_congruence';a=scaled(p,row['gauge_squared'])
 v=list(map(F,row['vector']));assert len(v)==len(ids)and any(v)and all(x>=0 for x in v)
 q=sum((v[i]*a[i][j]*v[j]for i in range(len(v))for j in range(len(v))),F())
 assert q==F(row['square'])and (q<0 if strict else q<=0)
 return q

def run():
 ds={r['index']:r for r in json.loads((ROOT/'results'/f'{CAT}.json').read_text())['diagrams']}
 data=json.loads((ROOT/'results'/f'{BASE}.json').read_text());graphs={};input_checks=Counter()
 for di,d in ds.items():
  n=d['order'];g=[[0]*n for _ in range(n)];a=[[K(2 if i==j else 0)for j in range(n)]for i in range(n)]
  for i,j,m in d['edges']:
   assert m in [3,4];g[i][j]=g[j][i]=m;a[i][j]=a[j][i]=K(-1)if m==3 else K(0,-1)
  assert inertia(a)==[n-1,1,0];input_checks['lorentz_cores']+=1
  for b in d['finite_deletions']:
   assert inertia(sub(a,[i for i in range(n)if i!=b]))==[n-1,0,0];input_checks['finite_deletions']+=1
  for b in d['ideal_deletions']:
   sig=inertia(sub(a,[i for i in range(n)if i!=b]));assert sig[1]==0 and sig[2]>0;input_checks['ideal_deletions']+=1
  graphs[di]=g
 @lru_cache(None)
 def inv(di,deleted):
  ids=[i for i in range(len(graphs[di]))if i not in deleted];a,sc=gauge(sub(graphs[di],ids));assert inertia(a)==[len(a),0,0]
  j=inverse(a);return {i:sc[k]*j[k][k]for k,i in enumerate(ids)}
 contexts=[];ports=[]
 for di,d in ds.items():
  g=graphs[di]
  for black in combinations(d['finite_deletions'],2):
   dist={i:0 for i in black};todo=deque(black)
   while todo:
    i=todo.popleft()
    for j,m in enumerate(g[i]):
     if m and j not in dist:dist[j]=dist[i]+1;todo.append(j)
   assert len(dist)==len(g);contexts.append((di,black,dist))
   for i in range(len(g)):
    if dist[i]<1:continue
    vals=[inv(di,(b,))[i]for b in black]
    ports.append(dict(diagram=di,black=list(black),node=i,distance=dist[i],inverse_diagonals=list(map(str,vals)),allowed=max(vals)<2))
 assert len(contexts)==384 and max(max(d.values())for _,_,d in contexts)==7
 assert ports==data['ports'];allowed=[p for p in ports if p['allowed']]
 assert len(ports)==2404 and len(allowed)==data['allowed_ports']==754
 def united(di,dj,edges):
  a,b=graphs[di],graphs[dj];n=len(a);g=[[0]*(n+len(b))for _ in range(n+len(b))]
  for i in range(len(a)):
   for j in range(len(a)):g[i][j]=a[i][j]
  for i in range(len(b)):
   for j in range(len(b)):g[n+i][n+j]=b[i][j]
  for i,j,m in edges:g[i][n+j]=g[n+j][i]=m
  return g,n
 expected={};good=defaultdict(list)
 for x,y in product(allowed,repeat=2):
  if x['distance']+y['distance']<8:continue
  di,dj=x['diagram'],y['diagram'];u,v=tuple(x['black']),tuple(y['black']);i,j=x['node'],y['node']
  wp=inv(di,u)[i]*inv(dj,v)[j];key=(di,u,i,dj,v,j);assert key not in expected
  mx=[max(map(F,z['inverse_diagonals']))for z in [x,y]]
  labels=[]if wp>=1 else[3]if max(mx)>=1 or 2*wp>=1 else[3,4]
  if labels==[3,4]:assert F(5,2)*max(mx)>=2
  expected[key]=(wp,mx,labels)
  if labels:good[di,u,dj,v].append((i,j,tuple(labels)))
 counts=Counter();type_keys={}
 for row in data['single_types']:
  pp=tuple(tuple(x)for x in row['ports']);key=pp+(row['label'],);assert key not in type_keys
  (di,i),(dj,j)=pp;g,n=united(di,dj,[(i,j,row['label'])])
  if row['method']=='rational_inertia':
   # Type keys are unordered; the saved matrix may retain either block order.
   a=None
   for dd,ii,ee,jj in [(di,i,dj,j),(dj,j,di,i)]:
    gg,nn=united(dd,ee,[(ii,jj,row['label'])])
    try:a=check_matrix(row['gram'],gg,row['gauge_squared']);break
    except AssertionError:pass
   assert a is not None
   sig=inertia(a);assert sig==row['signature']and sig[1]==row['negative_index']
   counts['explicit_single_inertias']+=1
  else:
   assert row['method']=='zero_inverse_rank_one_schur';zi,zn=row['zero_endpoint'];assert (zi,zn)in pp
   ids=[k for k in range(len(graphs[zi]))if k!=zn];z=row['zero_cofactor'];a=check_matrix(z['deletion_gram'],sub(graphs[zi],ids),z['gauge_squared'])
   sig=inertia(a);assert sig==z['deletion_signature']and sig[1]==0 and sig[2]>0
   assert row['negative_index']==2 and row['zero_inverse_diagonal']==z['inverse_diagonal']=='0'
   counts['zero_cofactor_schur_types']+=1
  type_keys[key]=row
 seen=set()
 for row in data['single_joins']:
  di,u,i=row['left'];dj,v,j=row['right'];u=tuple(u);v=tuple(v);key=(di,u,i,dj,v,j);assert key not in seen;seen.add(key)
  wp,mx,labs=expected[key];assert F(row['white_inverse_product'])==wp
  g,n=united(di,dj,[(i,j,3)]);black=set(u)|{n+k for k in v};white=sorted(set(range(len(g)))-black)
  if not labs:
   ob=row['white_rejection'];assert ob['nodes']==white;witness(ob,g);counts['white_witnesses']+=1;continue
  assert row['labels']==labs;gate=row['higher_label_gate'];assert list(map(F,gate['endpoint_maxima']))==mx
  assert gate['first_excluded_label']==(4 if labs==[3] else 5)
  if labs==[3]:assert F(gate['twice_white_product'])==2*wp
  else:assert F(gate['five_halves_maximum'])==F(5,2)*max(mx)
  assert [r['label']for r in row['label_rejections']]==labs
  for r in row['label_rejections']:
   m=r['label'];gg,nn=united(di,dj,[(i,j,m)])
   if r['kind']=='two_negative_directions':
    tk=tuple(sorted(((di,i),(dj,j))))+(m,);assert type_keys[tk]['negative_index']==2;counts['single_index2_cases']+=1
   else:
    assert r['kind']=='actual_corner_not_psd';p,q=r['corner'];assert p in u and q in {n+k for k in v};assert r['nodes']==sorted(white+[p,q])
    witness(r,gg,True);counts['strict_corner_witnesses']+=1
 assert seen==set(expected)and len(seen)==12135
 multi={}
 for key,edges in good.items():
  for size in range(2,len(edges)+1):
   for es in combinations(edges,size):multi[key+((tuple((i,j)for i,j,l in es)),)]=es
 seen=set()
 for row in data['multiple_edge_subsets']:
  di,u=row['left'];dj,v=row['right'];u=tuple(u);v=tuple(v);edges=tuple(tuple(e)for e in row['edges']);key=(di,u,dj,v,edges)
  assert key in multi and key not in seen;seen.add(key);es=multi[key];g,n=united(di,dj,[(i,j,3)for i,j in edges]);black=set(u)|{n+k for k in v};white=set(range(len(g)))-black
  def edge_witness(ob,gg):
   assert ob['mark']in black and ob['nodes']==sorted(white|{ob['mark']});witness(ob,gg)
  if 'all_labels_actual_edge_rejection'in row:
   edge_witness(row['all_labels_actual_edge_rejection'],g);counts['monotone_multiple_witnesses']+=1
  else:
   labs=list(product(*(e[2]for e in es)));assert [tuple(r['labels'])for r in row['label_cases']]==labs;counts['residual_multiple_subsets']+=1
   for r in row['label_cases']:
    gg,nn=united(di,dj,[(i,j,m)for (i,j),m in zip(edges,r['labels'])])
    if r['kind']=='actual_edge_not_pd':edge_witness(r,gg);counts['residual_edge_witnesses']+=1
    else:
     assert r['kind']=='two_negative_directions';a=check_matrix(r['gram'],gg,r['gauge_squared']);sig=inertia(a);assert sig==r['signature']and sig[1]==2;counts['residual_index2_inertias']+=1
 assert seen==set(multi)and len(seen)==10990
 assert counts['strict_corner_witnesses']==20 and counts['residual_multiple_subsets']==3 and counts['residual_edge_witnesses']==3 and counts['residual_index2_inertias']==3
 paths=[Path(__file__),ROOT/'scripts'/f'{BASE}.py',ROOT/'results'/f'{BASE}.json',ROOT/'results'/f'{BASE}.md',ROOT/'results'/f'{CAT}.json']
 return dict(status='PASS_INDEPENDENT_DISTANCE8_HIGH_DISJOINT_CERTIFICATE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},input_core_checks=dict(input_checks),marked_contexts=len(contexts),ports=len(ports),allowed_ports=len(allowed),single_joins=len(expected),multiple_subsets=len(multi),checks=dict(counts),scope='Independent Fraction/Qsqrt2 saved-certificate coverage and witness/inertia checks. No producer/helper arithmetic imports. Complete high-core input catalogue is inherited, but its26 listed core inertias and all listed finite/ideal deletions are independently verified. This does not rederive the catalogue completeness or prove any overlap/extension/ambient dimension theorem.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],out['checks'],flush=True)
