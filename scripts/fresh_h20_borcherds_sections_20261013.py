#!/usr/bin/env python3
"""Exact high-support section gates for the frozen Borcherds H21 chamber.

Standard-library integer/Fraction replay. The geometric implications and
continuous parameter coverage are stated in the matching report.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
from collections import Counter
import argparse, hashlib, json, math

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_borcherds_sections_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
OUT=ROOT/'results'/f'{TAG}.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(x,y):return sum(a*b for a,b in zip(x,y,strict=True))
def mv(A,x):return [dot(r,x)for r in A]
def pair(x,y,Q):return dot(x,mv(Q,y))
def add(x,y):return [a+b for a,b in zip(x,y,strict=True)]
def scale(a,x):return [a*b for b in x]
def inv(A):
 n=len(A);M=[[F(v)for v in row]+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
 for j in range(n):
  k=next(i for i in range(j,n)if M[i][j]);M[j],M[k]=M[k],M[j]
  z=M[j][j];M[j]=[v/z for v in M[j]]
  for i in range(n):
   if i!=j:
    z=M[i][j];M[i]=[x-z*y for x,y in zip(M[i],M[j],strict=True)]
 return [r[n:]for r in M]
def projection(h,ids,W,G):
 A=[[G[i][j]for j in ids]for i in ids];B=inv(A)
 assert all(x>=0 for r in B for x in r)
 cc=mv(B,[dot(h,COV[i])for i in ids])
 return [h[j]-sum(a*W[i][j]for a,i in zip(cc,ids,strict=True))for j in range(22)]
def enc(x):return str(x)if isinstance(x,F)else x
def vec(x):return [enc(v)for v in x]
def p_add(a,b):
 return [(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b)))]
def p_mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def vpoly_pair(A,B,Q):
 ans=[0]*(len(A)+len(B)-1)
 for i,x in enumerate(A):
  for j,y in enumerate(B):ans[i+j]+=pair(x,y,Q)
 return ans

def run():
 global COV
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
 Q=[[-int(x)for x in row]for row in raw['basis_gram']]
 W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 assert len(W)==210 and all(len(r)==22 for r in W)
 COV=[mv(Q,r)for r in W];G=[[dot(x,c)for c in COV]for x in W]
 assert all(G[i][j]<=0 for i,j in combinations(range(210),2))
 short=[i for i in range(210)if G[i][i]==2];assert len(short)==42
 assert Counter(G[i][i]for i in range(210))=={2:42,4:168}
 s=[sum(r[j]for r in W)for j in range(22)];gg=math.gcd(*s);h=[x//gg for x in s]
 assert pair(h,h,Q)==-14 and set(dot(h,c)for c in COV)=={-1,-8}
 color={short[0]:0};queue=[short[0]]
 for i in queue:
  for j in short:
   if G[i][j]==-1:
    if j not in color:color[j]=1-color[i];queue.append(j)
    assert color[j]!=color[i]
 assert len(color)==42
 parts=[[i for i in short if color[i]==k]for k in range(2)]
 assert list(map(len,parts))==[21,21]
 assert all(G[i][j]==2*(i==j)for part in parts for i in part for j in part)
 assert all(sum(G[i][j]==-1 for j in parts[1-k])==5 for k in range(2)for i in parts[k])
 assert all(sum(G[i][l]==G[j][l]==-1 for l in parts[1-k])==1
            for k in range(2)for i,j in combinations(parts[k],2))
 sums=[[sum(W[i][j]for i in part)for j in range(22)]for part in parts]
 a,b=sums
 assert pair(a,a,Q)==pair(b,b,Q)==42 and pair(a,b,Q)==-105
 assert add(a,b)==scale(3,h)
 assert all(dot(a,COV[i])==dot(b,COV[i])==-12 for i in range(210)if i not in short)
 walk=json.loads((ROOT/WALK).read_text())
 assert walk['status']=='PASS_COMPLETE_CORNER_ORBIT_CENSUS' and walk['sources'][RAW]==sha(RAW)
 assert not walk['pending_orbit_ids']
 vertices=[];availability=[];examples=[]
 for side,part in enumerate(parts):
  v=add(h,scale(F(1,2),sums[side]));vv=pair(v,v,Q)
  ps=[dot(v,c)for c in COV]
  assert vv==F(-49,2) and pair(v,h,Q)<0
  assert [i for i,p in enumerate(ps)if p==0]==part and all(ps[i]<0 for i in range(210)if i not in part)
  # Twenty-one orthogonal positive vectors and this negative vector give
  # the full ambient signature (21,1), without a numerical eigenvalue test.
  assert all(dot(v,COV[i])==0 for i in part)
  vertices.append({'walls':part,'point':vec(v),'norm':str(vv),'maximum_other_pairing':str(max(p for p in ps if p))})
  counts=Counter()
  for p,q,r in permutations(part,3):
   choices=[l for l in parts[1-side]if G[l][q]==-1 and G[l][p]==G[l][r]==0]
   assert len(choices)in(3,4);counts[len(choices)]+=1
  availability.append({str(k):v for k,v in sorted(counts.items())})
  p,q,r=part[:3]
  l=next(l for l in parts[1-side]if G[l][q]==-1 and G[l][p]==G[l][r]==0)
  faces=[]
  for ids in ([p,q,l],[r,p,l],[r,q,l]):
   z=projection(h,ids,W,G);pp=[dot(z,c)for c in COV]
   assert pair(z,z,Q)<0 and pair(z,h,Q)<0
   assert [i for i,x in enumerate(pp)if x==0]==sorted(ids) and max(pp)==0
   faces.append({'walls':ids,'point':vec(z),'norm':str(pair(z,z,Q)),'maximum_other_pairing':str(max(x for x in pp if x))})
  examples.append({'p':p,'q':q,'r':r,'auxiliary':l,'ambient_faces':faces})
 # The full plane span(a,b) has exactly e=a-t*b, t>0, as its interior
 # sections, up to rescaling: same-sign coefficients pair strictly with P.
 pencil=[]
 for side in range(2):
  aa,bb=sums[side],sums[1-side];i,j=parts[side][:2]
  e=[aa,scale(-1,bb)]
  center=[add(scale(5,aa),scale(2,bb)),add(scale(2,aa),scale(5,bb))]
  assert vpoly_pair(center,e,Q)==[0,0,0]
  cq=vpoly_pair(center,center,Q);assert all(x<0 for x in cq)
  assert all(all(x<=0 for x in [dot(z,COV[k])for z in center]) and any(dot(z,COV[k])<0 for z in center)for k in range(210))
  A=add(scale(2,aa),scale(-2,add(W[i],W[j])))
  B=add(scale(2,bb),scale(5,add(W[i],W[j])))
  point=[add(scale(95,A),scale(38,B)),add(scale(17,A),scale(95,B))]
  assert vpoly_pair(point,e,Q)==[0,0,0]
  nq=vpoly_pair(point,point,Q);np=vpoly_pair(point,[h],Q)
  assert all(x<0 for x in nq+np)
  pp=[[dot(z,c)for z in point]for c in COV]
  assert [k for k,p in enumerate(pp)if p==[0,0]]==[i,j]
  # Strict for every 0<t<=1, including all nonincident old inequalities.
  assert all(p[0]<=0 and sum(p)<0 for k,p in enumerate(pp)if k not in[i,j])
  qe=vpoly_pair(e,e,Q);assert qe==[42,210,42]
  v=[dot(z,COV[i])for z in e];assert v==[2,5]
  assert [dot(z,COV[j])for z in e]==v
  vv=p_mul(v,v);diagonal=p_add([2*x for x in qe],[-x for x in vv])
  gap=p_add(diagonal,[-2*x for x in vv]);assert gap==[72,360,9]
  # Projected diagonal=diagonal/q, offdiagonal=-vv/q.
  # Hence 0<cos<1/2, and determinant>0, for every positive t.
  assert all(x>0 for x in vv+gap)
  pencil.append({'side':side,'parameter_domain':'0<t<=1','normal_coefficients':e,'strict_center_coefficients':center,
                 'ridge_walls':[i,j],'ridge_point_coefficients':point,'ridge_point_norm_coefficients':nq,
                 'ridge_point_center_coefficients':np,'all_wall_pairing_coefficients':pp,
                 'normal_norm_coefficients':qe,'scaled_diagonal':diagonal,'negative_scaled_offdiagonal':vv,
                 'diagonal_minus_twice_absolute_offdiagonal':gap})
 return {'status':'PASS_EXACT_ORTHOGONAL_VERTEX_AND_BIPARTITION_SECTION_OBSTRUCTIONS',
  'sources':{p:sha(p)for p in(RAW,WALK,'scripts/'+Path(__file__).name)},'base_dimension':21,'section_dimension':20,
  'ambient_finite_volume_dependency':WALK,'bipartitions':parts,'bipartition_sums':sums,'strict_center':h,
  'orthogonal_vertices':vertices,'ordered_triple_auxiliary_line_count_distributions':availability,
  'representative_actual_face_points':examples,'bipartition_pencil':pencil,
  'scope':'Every interior codimension-one section through either recorded A1^21 vertex is non-Coxeter, for arbitrary real normal coefficients. Every interior section with normal in span(sum_points,sum_lines) is non-Coxeter, for every real coefficient ratio. All such genuine interior sections have finite volume inherited from the frozen H21 polytope. Other section normals and other ambient chambers remain open; no new Coxeter polytope is claimed.'}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');aa=ap.parse_args();out=run()
 if aa.verify:assert json.loads(OUT.read_text())==out
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
