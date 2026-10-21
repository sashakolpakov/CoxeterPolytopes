#!/usr/bin/env python3
"""Separate sparse-radical arithmetic replay of the all-E7-face escape gate.

No producer or project imports. Fraction elimination and rational192bit
enclosures replace the primary recursive quadratic sign implementation.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_global_escape_gate_20261014'
P=f'results/{TAG}.json';CH='results/fresh_h20_e8_interior_caps_20261014_exact.json';FA='results/fresh_h20_e7_face_20261014.json'
RAD=[1,2,79,158,11,22,869,1738];Z=(F(0),)*8;DEN=1<<192
BOUND=[]
for q in RAD:
 n=isqrt(q*DEN*DEN);BOUND.append((F(n,DEN),F(n+int(n*n!=q*DEN*DEN),DEN)))
def field(v):return tuple(map(F,v))
def vector(v):return list(map(field,v))
def c(x):return(F(x),)+(F(0),)*7
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def scale(x,s):return tuple(a*s for a in x)
def sub(x,y):return add(x,scale(y,-1))
def mul(x,y):
 z=[F(0)]*8
 for i,a in enumerate(x):
  for j,b in enumerate(y):
   if a and b:z[i^j]+=a*b*RAD[i&j]
 return tuple(z)
def total(values):
 v=Z
 for x in values:v=add(v,x)
 return v
def sign(x):
 if x==Z:return 0
 lo=sum(a*(u if a>=0 else v)for a,(u,v)in zip(x,BOUND));hi=sum(a*(v if a>=0 else u)for a,(u,v)in zip(x,BOUND))
 assert lo>0 or hi<0,('unresolved192bit',x,lo,hi)
 return 1 if lo>0 else -1
def inv(A):
 n=len(A);B=[list(map(F,row))+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
 for k in range(n):
  p=next(i for i in range(k,n)if B[i][k]);B[k],B[p]=B[p],B[k];a=B[k][k];B[k]=[x/a for x in B[k]]
  for i in range(n):
   if i!=k:
    a=B[i][k];B[i]=[x-a*y for x,y in zip(B[i],B[k])]
 return[r[n:]for r in B]
def pd(A):
 B=[list(map(F,r))for r in A]
 for k in range(len(B)):
  p=B[k][k];assert p>0
  for i in range(k+1,len(B)):
   for j in range(k+1,len(B)):B[i][j]-=B[i][k]*B[j][k]/p
def rational(x):assert x[1:]==(0,)*7;return x[0]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[enc(v)for v in x]
 return x
def run():
 p,d,f=map(read,[P,CH,FA]);Q=[list(map(F,r))for r in d['form']];cache={};dc={}
 def dot(v,w):
  key=(tuple(v),tuple(w))
  if key in dc:return dc[key]
  if key[1]not in cache:cache[key[1]]=[tuple(sum(Q[i][j]*w[j][k]for j in range(21))for k in range(8))for i in range(21)]
  out=total(mul(a,b)for a,b in zip(v,cache[key[1]]));dc[key]=out;dc[(key[1],key[0])]=out;return out
 N=list(map(vector,d['normals']));E=vector(d['strict_escape']['point']);H=vector(d['strict_center']);T=vector(f['intrinsic_normalization']['timelike_basis']);omega=vector(f['unique_nonright_meeting_profile']['fundamental_weight']);A=f['face_walls'];S=f['E7_support']
 assert sign(dot(E,E))>0 and sign(dot(E,H))<0 and all(sign(dot(E,r))<0 for r in N)
 assert sign(dot(T,H))<0 and dot(T,T)==c(-2)and dot(omega,omega)==c(F(3,2))and dot(omega,T)==Z
 assert all(dot(omega,N[i])==c(int(i==33))for i in S)and all(dot(omega,N[i])==Z for i in A)
 G=[[rational(dot(N[i],N[j]))for j in A]for i in A];assert enc(G)==f['face_gram'];pd(G);J=inv(G)
 CJ=S+A;CG=[[rational(dot(N[i],N[j]))for j in CJ]for i in CJ];pd(CG);CI=inv(CG)
 assert min(x for row in CI for x in row)>=0 and all(dot(T,N[i])==Z for i in CJ)
 assert CJ==p['finite_corner']['walls']and enc(CI)==p['finite_corner']['inverse']
 W=[[total(scale(N[j][k],q)for j,q in zip(A,row))for k in range(21)]for row in J]
 assert all(dot(W[i],N[j])==c(int(i==k))for i in range(13)for k,j in enumerate(A))
 a=[scale(dot(w,E),-1)for w in W];B=scale(dot(T,E),-1);w=scale(dot(omega,E),-1)
 assert all(sign(x)>0 for x in a+[B,w]);assert enc(dict(a=a,B=B,w=w))==p['target_coefficients']
 co=[[dot(lam,n)for lam in W]+[scale(dot(T,n),F(-1,2))]for n in N]
 assert enc(co)==f['intrinsic_normalization']['projected_wall_coordinates']
 ow=[dot(omega,n)for n in N]
 mu={30:F('0.002038242873333'),75:F('0.002461273481333'),76:F('1.639210640401')}
 coeff=[list(map(lambda x:scale(x,-1),row[:13]))+[scale(row[13],-2)]for row in co]
 residual=[sub(x,total(scale(coeff[j][i],s)for j,s in mu.items()))for i,x in enumerate(a+[scale(B,-1)])]
 assert all(sign(x)<0 for x in residual)
 assert enc(dict(multipliers=mu,residual_coefficients=residual))==p['right_lift']
 rows={};rhs={}
 for j in [16,76]:
  assert sign(co[j][13])>0
  rows[j]=[sub(scale(mul(B,co[j][i]),-1),scale(mul(co[j][13],a[i]),2))for i in range(13)]
  rhs[j]=add(mul(B,ow[j]),scale(mul(co[j][13],w),2))
 mixed=[add(scale(rows[16][i],F(194,1000)),scale(rows[76][i],F(843,1000)))for i in range(13)]
 mixed_rhs=add(scale(rhs[16],F(194,1000)),scale(rhs[76],F(843,1000)))
 assert all(sign(x)>0 for x in mixed)and sign(mixed_rhs)>0 and sign(sub(mixed[12],mixed_rhs))>0
 assert sign(rows[16][0])<0 and sign(sub(rows[16][0],rhs[16]))>0 and all(sign(rows[16][i])>0 for i in range(1,12))
 B2=mul(B,B);rankgap=sub(B2,scale(total(scale(mul(a[i],a[j]),G[i][j])for i in range(13)for j in range(13)),2));assert sign(rankgap)>0
 sqrt2=(F(0),F(1))+(F(0),)*6;C=add(mul(a[0],sqrt2),w);D=sub(B2,scale(mul(a[11],a[11]),3));assert sign(C)>0 and sign(D)>0
 assert J[0][0]==F(11,12)and J[11][11]==F(2,3)and all(J[i][11]==0 for i in range(11))
 grads=[sub(scale(mul(D,sqrt2),J[i][0]),scale(mul(a[i],C),2))for i in range(11)];assert all(sign(x)>0 for x in grads)
 margin=sub(scale(D,F(4,3)),scale(mul(C,C),2));assert sign(margin)>0
 result=dict(cut_boundary_rows=rows,cut_boundary_rhs=rhs,mixed_coefficients=mixed,mixed_rhs=mixed_rhs,rank_one_positive_gap=rankgap,denominator=D,C=C,positive_gradient_numerators=grads,minimum_norm_minus_half_times_D=margin)
 assert enc(result)==p['nonright_lift']and enc(dot(E,E))==p['escape_norm']
 paths=[P,CH,FA,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py']
 return enc(dict(status='PASS_SEPARATE_SPARSE_FIELD_ALL_E7_FACE_ESCAPE_GATE',sources={x:sha(x)for x in paths},field='Q(sqrt2,sqrt79,sqrt11)',sign_method='outward192bit rational radical enclosures',
  independently_reconstructed=dict(old_escape_signs=87,finite_corner_rank=20,corner_inverse_nonnegative=True,reduced_wall_coordinate_vectors=87,right_dual_residuals=14,nonright_linear_rows=2,positive_KKT_gradients=11),
  nonright_norm_gap_times_D=margin,proof_scope='All actual additional E7-face facets compatible with retained87 and a common strict interior. Both allowed lift types and all finite labels/parallel/dotted pairs included. Changing retainedwalls or adding walls that are not facets on thisface is outside scope.',
  independence='No project/producer imports. Physical field coefficients, inverses, duals, elimination rows and every positivity gate reconstructed independently. No repeat of old87 complete pair/corner census.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();out=ROOT/f'results/{TAG}_review.json'
 if args.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
