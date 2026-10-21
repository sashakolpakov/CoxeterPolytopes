#!/usr/bin/env python3
"""Independent Fraction/Q(sqrt5) replay of the compressed binary orbit.

No producer/project imports and no million-by-million Gram construction.
"""
from fractions import Fraction as F
from itertools import combinations, product
from collections import Counter
from hashlib import sha256
from math import comb
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_weighted_binary_20261014'
ZERO=(F(0),F(0));ONE=(F(1),F(0));PHI=(F(1,2),F(1,2))
def rat(x):return(F(x),F(0))
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def sub(x,y):return add(x,neg(y))
def scale(x,t):return(x[0]*t,x[1]*t)
def mul(x,y):return(x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def sign(x):
 a,b=x
 if not a:return(b>0)-(b<0)
 if not b:return(a>0)-(a<0)
 if a*b>0:return(a>0)-(a<0)
 d=a*a-5*b*b;assert d
 return((a>0)-(a<0))*((d>0)-(d<0))
def decode(x):
 a,b=map(F,x);return(a+b/2,b/2)
WEIGHTS=[ONE]*10+[rat(2)]+[add(ONE,PHI)]*8
def dot(x,y):
 s=ZERO
 for w,a,b in zip(WEIGHTS,x,y,strict=True):s=add(s,scale(w,a*b))
 return s
def distance(x,y):
 d=[a-b for a,b in zip(x,y,strict=True)];return dot(d,d)
def center(mask):
 parity=sum((mask>>j)&1 for j in list(range(10))+list(range(11,19)))%2
 full={0,1,2}if parity==0 else{0,3,4}
 return[F((-1 if(mask>>j)&1 else 1)*(2 if j in full else 1),2)for j in range(19)]
def relation(d):
 c=scale(sub(d,rat(2)),F(1,2))
 if c==ZERO:return'2'
 if c==rat(F(1,2)):return'3'
 if c==scale(PHI,F(1,2)):return'5'
 assert sign(sub(c,ONE))>=0
 return'0'
def masks_at_counts(row):
 if row[0]=='same':
  ha,hn,hh,hg=row[1:5]
  flip=list(range(3))[:ha]+list(range(3,10))[:hn]
  multiple=comb(3,ha)*comb(7,hn)*comb(8,hg)
  assert(ha+hn+hg)%2==0
 else:
  hi,hd,hn,hh,hg=row[1:6]
  flip=[0][:hi]+[1,2,3,4][:hd]+[5,6,7,8,9][:hn]
  multiple=comb(4,hd)*comb(5,hn)*comb(8,hg)
  assert(hi+hd+hn+hg)%2==1
 return sum(1<<j for j in flip+[10]*hh+list(range(11,19))[:hg]),multiple

def run():
 path=ROOT/f'results/{TAG}.json';data=json.loads(path.read_text())
 assert data['status']=='PASS_EXACT_PARITY_TRIPLE_CAP_ORBIT_WITH_GLOBAL_ESCAPE'
 for p,h in data['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
 reference=center(0);CN=dot(reference,reference)
 assert CN==decode(data['center_square'])
 central=scale(sub(CN,rat(5)),F(1,4))
 assert central==decode(data['central_cosine'])and sign(sub(central,ONE))>0
 assert sign(add(rat(100),scale(sub(CN,ONE),F(1,2))))>0
 seedcounts=Counter()
 for w,x in zip(WEIGHTS,reference,strict=True):
  for e in[-1,1]:
   gap=1-e*x;assert gap>=0
   csq=scale(w,gap*gap)
   pool={ZERO:'right',rat(F(1,4)):'pi3',rat(F(1,2)):'pi4',scale(add(ONE,PHI),F(1,4)):'pi5'}
   if csq in pool:seedcounts[pool[csq]]+=1
   else:assert sign(sub(csq,ONE))>0;seedcounts['dotted']+=1
 assert dict(seedcounts)==data['seed_relations']
 # The retained central radius-two wall has seed labels3/4/5.
 assert all(scale(w,F(1,4))in[rat(F(1,4)),rat(F(1,2)),scale(add(ONE,PHI),F(1,4))]
            for w in WEIGHTS)

 counts={'same':Counter(),'different':Counter()};row_count=0
 for row in data['compressed_pair_cases']:
  mask,multiple=masks_at_counts(row);c=center(mask)
  assert dot(c,c)==CN
  d=distance(reference,c);assert d==decode(row[-3])
  typ='identity'if mask==0 else relation(d)
  assert typ==str(row[-2])and multiple==row[-1]
  counts[row[0]][typ]+=multiple;row_count+=1
 assert row_count==828
 assert dict(counts['same'])==data['reference_same_type_counts']
 assert dict(counts['different'])==data['reference_other_type_counts']
 assert all(sum(c.values())==2**18 for c in counts.values())

 neighbors=[]
 for k in range(4):
  for uu in combinations(range(10),k):
   for hh in[0,1]:
    for gg in[()]+[(j,)for j in range(11,19)]:
     mask=sum(1<<j for j in uu+gg)+(hh<<10)
     if not mask:continue
     typ=relation(distance(reference,center(mask)))
     if typ!='0':neighbors.append([mask,int(typ)])
 assert neighbors==data['finite_neighbor_masks']and len(neighbors)==95
 assert len({r[0]for r in neighbors})==95
 assert sum(counts[t][label]for t in counts for label in['2','3','5'])==95

 # Matched half-probe and all endpoint/cap cross pairs use exact minimum
 # distances; signs are transferred because every coordinate is nonzero.
 assert sub(distance(reference,[F(1,2)]*19),ONE)==rat(F(-1,4))
 N=dot([F(1)]*19,[F(1)]*19)
 assert sign(sub(N,rat(16)))>0
 min_cross=distance(reference,[F(1)]*19)
 assert min_cross==decode(data['binary_endpoint_min_distance_square'])
 assert sign(sub(min_cross,rat(6)))>0 # 3+2sqrt2 <6
 endpoint_counts=Counter()
 for h1,hh,hg in product(range(11),range(2),range(9)):
  H=add(rat(h1+2*hh),scale(add(ONE,PHI),hg))
  if H==ZERO:continue
  cosine=sub(H,ONE)
  assert cosine==ZERO or sign(sub(cosine,ONE))>=0
  endpoint_counts['right'if cosine==ZERO else'parallel_or_dotted']+=comb(10,h1)*comb(8,hg)
 assert dict(endpoint_counts)==data['endpoint_reference_pair_counts']

 escape=data['strict_escape'];x=list(map(F,escape['spatial']));D=decode(escape['D'])
 assert sub(dot(x,x),scale(D,2))==ONE and sign(add(rat(100),D))>0
 assert min(1-e*v for v in x for e in[-1,1])==F(1,100)
 bp=sub(distance(x,reference),ONE)
 ep=sub(distance(x,[F(1)]*19),rat(2));cp=sub(dot(x,x),rat(4))
 for key,val in [('binary_power',bp),('endpoint_power',ep),('central_power',cp)]:
  assert val==decode(escape[key])and sign(sub(val,ONE))>0
 # At the escape every unit coordinate is zero, so type changes do not
 # change their19/4 contribution. Matching nonunit signs minimizes power.
 assert sum(c*c for c in reference[:10])==F(19,4)
 assert data['coherent_actual_facets']==39+2**20
 assert data['preserved_finite_endpoints']==2**19
 sources=[f'scripts/{TAG}.py',f'results/{TAG}.json',f'scripts/{TAG}_review.py']
 return dict(status='PASS_INDEPENDENT_BINARY_COMPRESSED_GEOMETRY_REVIEW',
  sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sources},
  arithmetic='stdlib Fraction in Q(sqrt5); no producer/project imports',
  compressed_pair_cases=row_count,finite_neighbors=95,
  exact_seed_counts=dict(seedcounts),reference_pair_counts={k:dict(v)for k,v in counts.items()},
  actual_facets=39+2**20,preserved_finite_endpoints=2**19,
  scope='Independent compressed representative reconstruction, pair multiplicities, explicit seed relations, norm/center/escape checks, and proof review of all-pair transfer and actual facets/endpoints. No full million-wall matrix or finite-vertex census. Local frontier files are outside this review.')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=run();out=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert json.loads(out.read_text())==data
 else:out.write_text(json.dumps(data,indent=2)+'\n')
 print(data['status'])
