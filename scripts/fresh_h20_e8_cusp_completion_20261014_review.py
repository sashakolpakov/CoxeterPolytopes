#!/usr/bin/env python3
"""Separate small Fraction/interval replay of the two E8 cusp fillers.
No producer imports; inherited77 geometry is an explicit dependency.
"""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from math import isqrt
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e8_cusp_completion_20261014_review'
PRIMARY='results/fresh_h20_e8_cusp_completion_20261014.json';BASE='results/fresh_h20_contact_closure_20261014_h12_replacement.json'
Z=(F(0),)*4
def scalar(x):return(F(x),F(0),F(0),F(0))
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def scale(x,s):return tuple(a*s for a in x)
def mul(x,y):
 out=[F(0)]*4
 for i,a in enumerate(x):
  for j,b in enumerate(y):out[i^j]+=a*b*(2 if i&j&1 else 1)*(79 if i&j&2 else 1)
 return tuple(out)
DEN=1<<160;BOUNDS=[]
for rad in [1,2,79,158]:
 n=isqrt(rad*DEN*DEN);BOUNDS.append((F(n,DEN),F(n+int(n*n!=rad*DEN*DEN),DEN)))
def sign(x):
 if x==Z:return 0
 lo=sum(a*(v[0]if a>=0 else v[1])for a,v in zip(x,BOUNDS))
 hi=sum(a*(v[1]if a>=0 else v[0])for a,v in zip(x,BOUNDS))
 assert lo>0 or hi<0
 return 1 if lo>0 else -1
def matrix_inverse(M):
 n=len(M);A=[row[:]+[F(i==j)for j in range(n)]for i,row in enumerate(M)]
 for i in range(n):
  k=next(j for j in range(i,n)if A[j][i]);A[i],A[k]=A[k],A[i];a=A[i][i];A[i]=[x/a for x in A[i]]
  for j in range(n):
   if j!=i:
    a=A[j][i];A[j]=[x-a*y for x,y in zip(A[j],A[i])]
 return[row[n:]for row in A]
def pd(M):
 A=[r[:]for r in M];pivots=[]
 for i in range(len(A)):
  v=A[i][i];assert v>0;pivots.append(v)
  for j in range(i+1,len(A)):
   for k in range(i+1,len(A)):A[j][k]-=A[j][i]*A[i][k]/v
 return pivots
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(tuple,list)):return[enc(v)for v in x]
 return x
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/PRIMARY).read_text());base=json.loads((ROOT/BASE).read_text())
 for data in [d,base]:
  for p,h in data['sources'].items():assert sha(p)==h,p
 Q=[[F(x)for x in row]for row in d['form']]
 def vec(row):return[tuple(map(F,x))for x in row]
 N=list(map(vec,d['normals']));old=list(map(vec,base['normals']));assert N[:77]==old and len(N)==79
 H=vec(d['strict_center']);escape=vec(d['strict_escape']['point']);cache={}
 def dot(x,y):
  key=tuple(y)
  if key not in cache:cache[key]=[tuple(sum(Q[i][j]*y[j][k]for j in range(21))for k in range(4))for i in range(21)]
  out=Z
  for a,b in zip(x,cache[key]):out=add(out,mul(a,b))
  return out
 def linear(vs,cs):return[tuple(sum(c*v[j][k]for v,c in zip(vs,cs))for k in range(4))for j in range(21)]
 E=base['E8_support'];B=[]
 for i in E:
  row=[]
  for j in E:
   g=dot(old[i],old[j]);assert g[1:]==(0,0,0);row.append(g[0])
  B.append(row)
 pd(B);I=matrix_inverse(B);assert I==[[F(x)for x in r]for r in base['E8_inverse']]
 assert I[0][0]==2 and all(x>0 for row in I for x in row)
 w=linear([old[i]for i in E],[-x for x in I[0]])
 rays=[vec(row['point'])for row in base['H12_face_ideal_vertices']]
 assert w==vec(d['negative_fundamental_weight'])and dot(w,w)==scalar(2)
 assert sign(dot(*rays))<0
 for k,nu in enumerate(rays):
  assert dot(nu,nu)==Z and dot(nu,w)==Z and sign(dot(nu,H))<0
  assert N[77+k]==[add(a,scale(v,8))for a,v in zip(w,nu)]
 norms=[]
 for r,q in zip(N,d['norms']):assert dot(r,r)==scalar(F(q))and F(q)>0;norms.append(F(q))
 counts=Counter()
 for i in [77,78]:
  for j in range(i):
   g=dot(N[i],N[j]);assert sign(g)<=0;sq=mul(g,g);prod=norms[i]*norms[j]
   if g==Z:kind='right'
   elif sq==scalar(prod/4):kind='pi3'
   else:assert sign(add(sq,scalar(-prod)))>0;kind='dotted'
   counts[kind]+=1
 assert counts==Counter(right=39,pi3=2,dotted=114)
 assert sign(dot(H,H))<0 and all(sign(dot(H,r))<0 for r in N)
 assert sign(dot(escape,escape))>0 and sign(dot(escape,H))<0 and all(sign(dot(escape,r))<0 for r in N)
 links=[]
 for nu,row in zip(rays,d['completed_cusp_links']):
  vals=[dot(nu,r)for r in N];inc=[i for i,g in enumerate(vals)if g==Z]
  assert inc==row['incident']and all(sign(g)<=0 for g in vals)
  seen=set();rank=0;blocks=[]
  for block in row['affine_blocks']:
   C=block['incident'];assert not seen.intersection(C);seen.update(C)
   A=[]
   for i in C:
    rr=[]
    for j in C:
     g=dot(N[i],N[j]);assert g[1:]==(0,0,0);rr.append(g[0])
    A.append(rr)
   weights=list(map(F,block['positive_kernel']));assert all(x>0 for x in weights)
   assert all(sum(a*x for a,x in zip(rr,weights))==0 for rr in A)
   pd([rr[:-1]for rr in A[:-1]])
   mu=F(block['null_multiple']);assert mu>0 and linear([N[i]for i in C],weights)==[scale(x,mu)for x in nu]
   assert all(dot(N[i],N[j])==Z for i in C for j in inc if j not in C)
   rank+=len(C)-1;blocks.append(len(C))
  assert seen==set(inc)and rank==19
  links.append(dict(block_sizes=sorted(blocks),quotient_rank=rank,all_incident_normals_in_positive_kernels=True))
 assert [r['block_sizes']for r in links]==[[4,9,9],[9,12]]
 for row in base['H12_face_finite_vertices']:
  p=vec(row['point']);assert sign(dot(p,p))<0 and all(sign(dot(p,r))<0 for r in N[77:])
 paths=[PRIMARY,BASE,'scripts/fresh_h20_e8_cusp_completion_20261014.py',f'scripts/{TAG}.py']
 return enc(dict(status='PASS_SEPARATE_FRACTION_79_CUSP_COMPLETION_REVIEW',sources={p:sha(p)for p in paths},new_pair_counts=dict(counts),completed_links=links,all79_norms_center_escape=True,preserved_finite_vertex_new_signs=70,escape_norm=dot(escape,escape),new_escape_pairings=[dot(escape,r)for r in N[77:]],scope='No producer imports. Separate sparse radical multiplication and160bit rational intervals recheck caps,155newpairs,all79norms/center/escape,positive full-link kernels/rank19,and70vertex-preservation signs. Actual facet geometry and inherited77/source completeness are reviewed certified dependencies, not independently re-enumerated.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
