#!/usr/bin/env python3
"""Narrow independent E8-face/projection/isometry replay; no producer imports."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_contact_closure_20261014_e8_face_review'
BASE='results/fresh_h20_contact_closure_20261014_replacement.json'
FACE='results/fresh_h20_contact_closure_20261014_e8_face.json'
NEXT='results/fresh_h20_contact_closure_20261014_h12_replacement.json'
TRANSFER='results/fresh_h20_contact_closure_20261014_h12_transfer.json'
SOURCE='results/fresh_h12_a12_source_20261014.json'
ZERO=(F(0),)*4
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def c(x):return(F(x),F(0),F(0),F(0))
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def scale(x,s):return tuple(a*s for a in x)
def mul(x,y):
 z=[F(0)]*4
 for i,a in enumerate(x):
  for j,b in enumerate(y):z[i^j]+=a*b*(2 if i&j&1 else 1)*(79 if i&j&2 else 1)
 return tuple(z)
def sgn(x):return int(x>0)-int(x<0)
def sign2(x):
 a,b=x;sa,sb=sgn(a),sgn(b)
 if not sa:return sb
 if not sb or sa==sb:return sa
 return sa*sgn(a*a-2*b*b)
def sq2(x):
 a,b=x;return(a*a+2*b*b,2*a*b)
def sign(x):
 A=x[:2];B=x[2:];sa,sb=sign2(A),sign2(B)
 if not sa:return sb
 if not sb or sa==sb:return sa
 aa,bb=sq2(A),sq2(B)
 return sa*sign2(tuple(a-79*b for a,b in zip(aa,bb)))
def vec(x):return[tuple(map(F,y))for y in x]
def ratmat(x):return[[F(str(v))for v in row]for row in x.tolist()]
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[enc(v)for v in x]
 return x
def run():
 d,f,n,tr,src=[load(p)for p in(BASE,FACE,NEXT,TRANSFER,SOURCE)]
 Q=[[F(x)for x in row]for row in d['form']];terms=[(i,j,Q[i][j])for i in range(21)for j in range(21)if Q[i][j]]
 def dot(x,y):
  z=ZERO
  for i,j,a in terms:z=add(z,scale(mul(x[i],y[j]),a))
  return z
 def linear(rows,co):
  return[tuple(sum(a*r[i][k]for a,r in zip(co,rows))for k in range(4))for i in range(21)]
 N=list(map(vec,d['normals']));H=list(map(c,d['strict_center']));E=f['E8_walls'];A=f['orthogonal_face_walls']
 B=sp.Matrix([[dot(N[i],N[j])[0]for j in E]for i in E])
 assert all(dot(N[i],N[j])==c(B[a,b])for a,i in enumerate(E)for b,j in enumerate(E))
 assert B.tolist()==[[F(x)for x in row]for row in f['gram']]
 _,D=B.LDLdecomposition(hermitian=False);assert all(D[i,i]>0 for i in range(8))and B.det()==1
 BI=ratmat(B.inv());assert BI==[[F(x)for x in row]for row in f['inverse']]
 assert all(x>0 for row in BI for x in row)
 assert [BI[i][i]for i in range(8)]==list(map(F,[2,6,12,20,8,30,14,4]))
 adjacency={i:[j for j in range(8)if j!=i and B[i,j]!=0]for i in range(8)}
 fork=[i for i in adjacency if len(adjacency[i])==3];assert len(fork)==1
 lengths=[]
 for start in adjacency[fork[0]]:
  prev,cur=fork[0],start;length=1
  while len(adjacency[cur])==2:
   nxt=next(x for x in adjacency[cur]if x!=prev);prev,cur=cur,nxt;length+=1
  lengths.append(length)
 assert sorted(lengths)==[1,2,4]
 P=vec(f['face_strict_interior']);assert sign(dot(P,P))<0 and sign(dot(P,H))<0
 eh=[dot(N[i],H)for i in E]
 coeff=[tuple(sum(BI[i][j]*eh[j][k]for j in range(8))for k in range(4))for i in range(8)]
 assert all(sign(x)<0 for x in coeff)
 projected=[add(h,scale(sumv,-1))for h,sumv in zip(H,[tuple(sum(mul(a,r[i])[k]for a,r in zip(coeff,[N[j]for j in E]))for k in range(4))for i in range(21)])]
 assert projected==P
 assert all(dot(P,r)==ZERO if i in E else sign(dot(P,r))<0 for i,r in enumerate(N))
 orth=[i for i,r in enumerate(N)if all(dot(r,N[j])==ZERO for j in E)];assert orth==sorted(A)
 for row in f['nonorthogonal_projected_walls']:
  i=row['wall'];pairs=[dot(N[j],N[i])for j in E];cost=ZERO
  for a in range(8):
   for b in range(8):cost=add(cost,scale(mul(pairs[a],pairs[b]),BI[a][b]))
  q=add(dot(N[i],N[i]),scale(cost,-1));orientation=dot(P,N[i])
  assert q==tuple(map(F,row['projection_norm']))and sign(q)<0
  assert orientation==tuple(map(F,row['orientation']))and sign(orientation)<0
 assert len(f['nonorthogonal_projected_walls'])==56
 C=sp.Matrix([[dot(N[i],N[j])[0]for j in A]for i in A])
 assert all(dot(N[i],N[j])==c(2 if a==b else -1 if abs(a-b)==1 else 0)for a,i in enumerate(A)for b,j in enumerate(A))
 assert C.det()==13
 T=vec(f['normalized_timelike_basis']);assert dot(T,T)==c(F(-1,13))
 assert all(dot(T,N[i])==ZERO for i in E+A)and sign(dot(T,H))<0
 # Reconstruct both reverse-chain completion normals independently.
 CI=ratmat(C.inv());caps=[]
 for j,k,q in((10,2,2),(2,3,1)):
  spatial=linear([N[i]for i in A],[-CI[j-1][i]for i in range(12)])
  r=[add(x,scale(t,F(k)))for x,t in zip(spatial,T)]
  assert dot(r,r)==c(q)and all(dot(r,N[e])==ZERO for e in E);caps.append(r)
 assert dot(*caps)==ZERO
 NN=list(map(vec,n['normals']));assert n['deleted_predecessor_wall']==0 and NN==N[1:]+caps
 assert 0 not in E+A and next(r for r in f['nonorthogonal_projected_walls']if r['wall']==0)
 # The only excluded predecessor for the reversed pair is old wall0.
 for cap,row in zip(caps,tr['orientations'][1]['caps']):assert cap==vec(row['normal'])
 assert [x['forbidden_walls']for x in tr['orientations'][1]['caps']]==[[0],[]]
 mapping={int(k):v for k,v in n['source_wall_map'].items()}
 for i in range(14):
  for j in range(14):assert dot(NN[mapping[i]],NN[mapping[j]])==c(F(src['gram'][i][j]))
 # Any further acute intrinsic normal must lie in the already complete
 # source's future chamber cone. Verify the literal inequalities used.
 SR=[[F(v)for v in row]for row in src['normals']]
 for i in range(11):assert SR[i]==[F(0)]+[F(-1 if j==i else 1 if j==i+1 else 0)for j in range(12)]
 assert SR[11]==[F(0)]*12+[F(-1)]
 assert SR[12]==[F(1)]*4+[F(0)]*9
 for row in src['projective_vertices']:
  x=list(map(F,row['point']));q=-x[0]*x[0]+sum(a*a for a in x[1:])
  assert q==F(row['norm'])and q<=0
 assert len(src['projective_vertices'])==37
 assert len(n['H12_face_finite_vertices'])==35 and len(n['H12_face_ideal_vertices'])==2
 # A common escape checks the scope, not all ambient vertex/pair certificates.
 X=vec(n['strict_escape']['point']);HH=vec(n['strict_center'])
 assert sign(dot(HH,HH))<0 and all(sign(dot(HH,r))<0 for r in NN)
 assert sign(dot(X,X))>0 and sign(dot(X,HH))<0 and all(sign(dot(X,r))<0 for r in NN)
 paths=[BASE,FACE,NEXT,TRANSFER,SOURCE,f'scripts/{TAG}.py',
  'scripts/fresh_h20_contact_closure_20261014_e8_face.py',
  'scripts/fresh_h20_contact_closure_20261014_h12_replacement.py']
 return enc(dict(status='PASS_INDEPENDENT_E8_FACE_MATRIX_PROJECTION_AND_COMPLETION_SCOPE_REVIEW',sources={p:sha(p)for p in paths},
  E8_arm_lengths=sorted(lengths),E8_determinant=1,E8_inverse=BI,minimum_inverse_diagonal=2,
  strict_timelike_redundant_projections=56,retained_after_wall0_deletion=55,A12_determinant=13,face_form_determinant=-1,
  reconstructed_reverse_cap_norms=[2,1],source_gram_entries_checked=196,strict_escape_inequalities_checked=77,
  intrinsic_source_ray_norms_rechecked=37,
  further_orthogonal_cap_gate='The literal B12 inequalities force x1>=...>=x12>=0; the norm2 root at source index12 forces t>=x1+x2+x3. A nonzero acute vector has t>0, so the inherited complete t1 slice forces nonpositive norm. No additional spacelike E8-orthogonal acute normal exists after the14 source walls are retained.',
  scope='Separate Fraction multiquadratic arithmetic and SymPy matrix inversion; no producer arithmetic imports. Reconstructs E8 Gram/inverse and strict face point, all56 nonorthogonal projection norms, the A12 form, the two reverse caps and the14-root H12 Gram isometry. Rechecks the77-wall strict escape. It does not duplicate2926 ambient pair classifications, all35 finite-vertex checks or the separate H12 source volume census.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
