#!/usr/bin/env python3
"""Independent graph-distance and algebra audit of a Prokhorov-weight refinement."""
import argparse,json
from collections import deque
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1];TAG='dimension_bound_prokhorov_20261015_review'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def path(n):return [(i,i+1)for i in range(n-1)]
def fork(n):return [(0,2),(1,2)]+[(i,i+1)for i in range(2,n-1)]
def exceptional(n):
 edges=[];last=1
 for length in [1,2,n-4]:
  old=0
  for _ in range(length):edges.append((old,last));old=last;last+=1
 assert last==n
 return edges
def histogram(n,edges):
 adj=[[]for _ in range(n)]
 for a,b in edges:adj[a].append(b);adj[b].append(a)
 hist=[0]*n
 for a in range(n):
  dist=[-1]*n;dist[a]=0;q=deque([a])
  while q:
   u=q.popleft()
   for v in adj[u]:
    if dist[v]<0:dist[v]=dist[u]+1;q.append(v)
  assert min(dist)>=0
  for b in range(a+1,n):hist[dist[b]]+=1
 return hist
def run():
 rows=[]
 for n in range(1,81):
  families=[('path',path(n))]+([('D',fork(n))]if n>=4 else[])+([('E',exceptional(n))]if n in[6,7,8]else[])
  for family,edges in families:
   hist=histogram(n,edges);w3=sum((3 if 1<=k<=7 else 1 if 8<=k<=15 else 0)*h for k,h in enumerate(hist));deficit=29*n-w3
   if n>=17:
    assert family in['path','D']
    assert hist[1:16]==[n-k+(family=='D'and k==2)for k in range(1,16)]
    assert deficit==(176 if family=='path' else 173)
   else:assert 17*deficit>=173*n
   rows.append(dict(family=family,rank=n,distance_histogram=hist,three_times_weight=w3,deficit=deficit))
 # No graph census is needed for the all-small-component case at the target dimensions:
 # w3<=3*binom(r,2), so deficit>=r*(61-3r)/2>=13r/2 for1<=r<=16.
 for r in range(1,17):assert 29*r-F(3*r*(r-1),2)==F(r*(61-3*r),2)>=F(13*r,2)
 assert F(13*27,2)>173
 # Independent direct application of published Lemma2, C=(29N-173)/(3N).
 # Necessary: n<96C+68, equivalently n^2-997n+6532<0.
 for n in range(991,1992):
  N=n-1;c=F(29*N-173,3*N);assert c>0
  q=n*n-997*n+6532;t=n-991
  assert q==t*t+985*t+586>0
  assert F(n)-(96*c+68)==F(q,N)>0
 # Independent check against the stronger published parity incidence inequalities.
 parity=[]
 for n in[990,991,992,993,994,995,996,1991,1992]:
  N=n-1;W=F(29*N-173,3)
  if n%2:
   a31=F(12*n,n-3);a32=F(6*(n-1),n-3);bad=F(3*n+9,4*(n-2))
   poly=n**3-997*n*n+6679*n-147;den=4*(n-3)*(n-2)*(n-1)
  else:
   a31=F(12*(n-1),n-4);a32=F(6*(n-2),n-4);bad=F(3*n+6,4*(n-3))
   poly=n**3-1000*n*n+9604*n-19824;den=4*(n-4)*(n-3)*(n-2)
  gap=7-a32-bad-2*W*a31/F((n-1)*(n-2))
  assert gap==F(poly,den)
  if n>=991:assert gap>0
  else:assert gap<0
  parity.append(dict(dimension=n,lower_minus_upper=str(gap)))
 paths=[f'scripts/{TAG}.py']
 return dict(status='PASS_INDEPENDENT_PROKHOROV_WEIGHT_REFINEMENT_ARITHMETIC',sources={p:sha(p)for p in paths},
  arithmetic='Independent breadth-first graph distances and Fraction arithmetic, no primary or producer imports.',
  connected_graph_checks=len(rows),graphs=rows,large_component_deficits=dict(path=176,D=173),
  disconnected_argument='For totalN>=27: anyrank>=17 component hasdeficit>=173; else everyrank<=16component hasdeficit>=13r/2, givingtotaldeficit>173. Allotherdeficits are nonnegative.',
  graph_bound='W<= (29N-173)/3 for every finite spherical Coxeter forest of totalrankN>=27.',
  lemma_substitution='C_n=29/3-173/(3(n-1)); publishedLemma2 gives n<996-5536/(n-1).',
  contradiction_polynomial='n^2-997n+6532=(n-991)^2+985(n-991)+586>0 for n>=991',
  parity_checks=parity,
  scope='Arithmetic verification of a proposed refinement of the published995argument, with spherical-graph classification and publishedweight/facelemma as explicit geometric inputs. Literature status and the older stronger announcement require separate audit; this certificate alone is not a novelty/current-best claim.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['connected_graph_checks'])
