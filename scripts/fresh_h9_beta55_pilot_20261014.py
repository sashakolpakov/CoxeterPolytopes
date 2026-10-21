#!/usr/bin/env python3
"""Complete bounded Vinberg pilot on E8 + [[-2phi,-1],[-1,4]]."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import cmp_to_key
from hashlib import sha256
from math import lcm,isqrt,ceil
from pathlib import Path
import fresh_h9_beta11_vinberg_h300_20261014_review as A
K=A.K;PHI=K(0,1);ALPHA=2+PHI;BETA=4-PHI;DELTA=1+8*PHI
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_beta55_pilot_20261014'
SOURCE='results/fresh_h9_beta11_vinberg_h300_20261014.json'
ZERO='results/fresh_h9_beta11_zero_sector_20261014.json'
HELPER='scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py'
RING='results/fresh_h9_h8_cap_height_20261014_norm_pool_review.json'
HEIGHT=F(20)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d['sources'].items():assert sha(f)==h,f
 return d
def run(height_bound=HEIGHT):
 height_bound=F(height_bound);assert height_bound>0
 old=read(SOURCE);zero=read(ZERO);ring=read(RING)
 C=[[A.dec(x).a for x in row[1:9]]for row in old['form'][1:9]]
 assert all(A.dec(x).b==0 for row in old['form'][1:9]for x in row[1:9])
 J=A.inverse(C)
 assert all(v.denominator==1 and v>=0 for row in J for v in row)
 G=[[K()for _ in range(10)]for _ in range(10)]
 G[0][0]=-2*PHI;G[0][9]=G[9][0]=K(-1);G[9][9]=K(4)
 for i in range(8):
  for j in range(8):G[i+1][j+1]=K(C[i][j])
 H=[K(4)]+[K()]*8+[K(1)]
 assert A.pair(H,H,G)==-4*DELTA and A.mv(G,H)==[-DELTA]+[K()]*9
 assert DELTA==PHI*ALPHA*BETA and ALPHA*ALPHA.conjugate()==5 and BETA*BETA.conjugate()==11
 assert all(b.sign()>0 and b.sign(True)>0 for b in (ALPHA,BETA))
 assert (DELTA-13).sign()>0 and (-DELTA.conjugate()-3).sign()>0
 norms=[2*K(1),2*ALPHA,2*BETA,2*ALPHA*BETA]
 assert all((18-q).sign()>0 and (13-q.conjugate()).sign()>0 for q in norms)
 def sqrt_ceiling(f):
  n=isqrt(f.numerator//f.denominator)
  while n*n<f:n+=1
  return n
 # Uniform rational bounds, with exact filtering at both embeddings below.
 X=sqrt_ceiling(4*18*height_bound/F(13));Z=sqrt_ceiling(4*18*(1+height_bound))
 assert 4*18*height_bound/F(13)<=X*X and F(4*13,3)<25
 assert 4*18*(1+height_bound)<=Z*Z and 4*13<64
 assert (PHI-F(8,5)).sign()>0 and (F(17,10)-PHI).sign()>0
 xb=ceil(F(5*(X+5),11));xa=ceil(F(17*5+7*X,22))
 zb=ceil(F(5*(Z+8),11));za=ceil(F(17*8+7*Z,22))
 controls=[]
 for q in norms:
  for a in range(-xa,xa+1):
   for b in range(-xb,xb+1):
    x=K(a,b);target=4*q+DELTA*x*x
    if x.sign()>0 and (height_bound-DELTA*x*x/(4*q)).sign()>=0 and target.sign(True)>=0:
     controls.append((DELTA*x*x/(4*q),q,x))
 controls.sort(key=lambda row:row[0]);assert all(a[0]!=b[0]for a,b in zip(controls,controls[1:]))
 order=list(range(7,-1,-1))
 Finv=[[F(1)]+[F(0)]*8]+[[F(0)]+[C[i][j]/4 for j in order]for i in order]
 forms=[]
 for n in range(1,10):
  M=A.inverse([r[:n]for r in Finv[:n]]);den=lcm(*(z.denominator for r in M for z in r))
  forms.append((den,[[int(z*den)for z in r]for r in M]))
 simple=[[K(-int(i==j))for i in range(10)]for j in range(1,9)]
 assert simple==[list(map(A.old,r))for r in zero['simple_zero_walls'][:8]]
 # x=0: all norm2b reflector coordinates y,t are divisible by b. A nonunit
 # b then contradicts evenness of norm(r/b)=2/b. The q2 shell is ordinary E8.
 zero_e8=[list(map(K,r[:-1]))for r in zero['complete_integer_norm2_roots']if r[-1]==0]
 assert len(zero_e8)==240 and all(A.pair(r,r,[[K(x)for x in row]for row in C])==2 for r in zero_e8)
 previous=[];levels=[];total=Counter()
 for index,(height,q,x)in enumerate(controls):
  target=4*q+DELTA*x*x;ti=target.ri();xi=x.ri();div=q/2
  alphabet=[];lalpha=[]
  for a in range(-za,za+1):
   for b in range(-zb,zb+1):
    z=K(a,b);zi=(a,b)
    if not (z/div).integral():continue
    rem=target-z*z
    if rem.sign()<0 or rem.sign(True)<0:continue
    if (x-z).a.denominator==1 and (x-z).b.denominator==1 and all(t%4==0 for t in (a-xi[0],b-xi[1])):lalpha.append(zi)
    if z.sign()>=0 and (target-2*z*z).sign()>=0 and (target-2*z*z).sign(True)>=0:alphabet.append(zi)
  key=cmp_to_key(lambda a,b:A.sign(A.sub(a,b)))
  alphabet.sort(key=key);lalpha.sort(key=key)
  co=[];bounds=[]
  for r in previous:
   L=r[0]-4*r[9];c=[L]+[4*r[i+1]for i in order]
   assert all(z.sign()>=0 for z in c[1:])
   co.append([z.ri()for z in c]);bounds.append((DELTA*x*r[0]).ri())
  visited=[0]*10;schur=[0]*10;walls=[0]*10;roots=[]
  def recurse(prefix,minimum,partial):
   depth=len(prefix);visited[depth]+=1
   if depth:
    num,den=minimum;rem=A.sub(A.scale(ti,den),num)
    if A.sign(rem)<0 or A.sign(rem,True)<0:schur[depth]+=1;return
    z=prefix[-1];partial=[A.add(p,A.mul(z,c[depth-1]))for p,c in zip(partial,co)]
    # l has already been fixed; every still-unassigned lambda contribution
    # is nonnegative. No incorrect l>=0 restriction is imposed.
    if any(A.sign(A.sub(bound,p))<0 for bound,p in zip(bounds,partial)):walls[depth]+=1;return
   if depth==9:
    num,den=minimum
    if num!=A.scale(ti,den):return
    l=K(*prefix[0]);lam=[K(*z)for z in reversed(prefix[1:])];y=A.mv(J,lam);t=(x-l)/4;r=[x,*y,t]
    assert all(z.integral()for z in r)and A.pair(r,r,G)==q
    assert all((2*z/q).integral()for z in A.mv(G,r))
    assert all(A.pair(r,n,G).sign()<=0 for n in simple+previous)
    roots.append(r);return
   den,M=forms[depth];base=(0,0);linear=(0,0)
   for i in range(depth):
    linear=A.add(linear,A.scale(prefix[i],2*M[i][depth]))
    for j in range(depth):base=A.add(base,A.scale(A.mul(prefix[i],prefix[j]),M[i][j]))
   for z in (lalpha if depth==0 else alphabet):
    num=A.add(A.add(base,A.mul(linear,z)),A.scale(A.mul(z,z),M[depth][depth]))
    recurse(prefix+[z],(num,den),partial)
  recurse([],None,[(0,0)]*len(previous))
  assert all(A.pair(a,b,G).sign()<=0 for i,a in enumerate(roots)for b in roots[i+1:])
  previous+=roots;total['nodes']+=sum(visited);total['levels']+=1
  levels.append(dict(height=height,norm=q,control=x,target=target,lambda_alphabet=[K(*z)for z in alphabet],
   signed_l_alphabet=[K(*z)for z in lalpha],nodes=visited,schur_prunes=schur,wall_prunes=walls,roots=roots))
  print('level',index+1,'of',len(controls),'roots',len(roots),'walls',8+len(previous),'nodes',sum(visited),flush=True)
 N=simple+previous;gram=[[A.pair(r,s,G)for s in N]for r in N]
 assert all(gram[i][j].sign()<=0 for i in range(len(N))for j in range(i))
 angle_values={K():2,K(1):3,K(2):4,1+PHI:5,K(3):6,2+PHI:10,K(4):'parallel'}
 pair_rows=[];pair_counts=Counter()
 for i in range(len(N)):
  for j in range(i):
   value=4*gram[i][j]*gram[i][j]/(gram[i][i]*gram[j][j])
   if (value-4).sign()>0:kind='dotted'
   else:
    assert value in angle_values
    label=angle_values[value];kind='parallel'if label=='parallel'else'right'if label==2 else f'pi/{label}'
   pair_counts[kind]+=1;pair_rows.append(dict(pair=[j,i],cartan_cycle=value,kind=kind))
 Z=[K()]+A.mv(J,[K(1)]*8)+[K()];epsilon=F(1)
 while True:
  center=[a+epsilon*b for a,b in zip(H,Z)]
  if A.pair(center,center,G).sign()<0 and all(A.pair(center,n,G).sign()<0 for n in N):break
  epsilon/=2
 facets=[]
 for i,n in enumerate(N):
  point=[a-A.pair(center,n,G)/gram[i][i]*b for a,b in zip(center,n)]
  assert A.pair(point,point,G).sign()<0
  assert all(A.pair(point,r,G).sign()==(0 if j==i else -1)for j,r in enumerate(N))
  facets.append(point)
 sources={p:sha(p)for p in [f'scripts/{TAG}.py',HELPER,SOURCE,ZERO,RING]}
 return A.enc(dict(status='PASS_COMPLETE_BOUNDED_BETA55_PREFIX',sources=sources,height_bound=height_bound,
  arithmetic='Exact Fraction Q(phi) and integer-pair search; Schur forms reconstructed from literal E8 Cartan.',
  form=G,basepoint=H,determinant=-DELTA,norm_classes=norms,zero_simple_walls=simple,
  ordinary_zero_roots=zero_e8,
  inverse_E8=J,positive_coordinate_order=['l']+[f'lambda{i+1}'for i in order],
  uniform_bounds=dict(control_physical=X,control_conjugate=5,control_coefficient_box=[xa,xb],
   coordinate_physical=Z,coordinate_conjugate=8,coordinate_coefficient_box=[za,zb]),
  levels=levels,roots=N,gram=gram,strict_center=center,center_epsilon=epsilon,actual_facet_points=facets,
  pair_relations=pair_rows,pair_counts=dict(pair_counts),
  counts=dict(total,positive_walls=len(previous),walls=len(N)),
  completeness=dict(alphabet='Even primitive reflector norm2b implies b divides delta; squarefree primes alpha and beta give four totallypositive squareclasses.',
   zero_sector='E8 only: norm4 coordinate vector is not a reflector; nonunit b forces all x0 coordinates divisible by b, impossible by evenness.',
   controls='x²<72T/13,abs(x_conjugate)<5; explicit rational coefficient boxes, then exact both-place filtering.',
   coordinates='l²<72(1+T),abs(l_conjugate)<8; explicit rational coefficient boxes; lambda bound twice as strong.',
   congruences='lambda and l divisible by b; l=x mod4O; gcd(4,b)=1 gives remaining binary reflection filter.',
   previous_wall_pruning='Only after signed l is fixed; remaining 4*lambda·Y coefficients are nonnegative.'),
  scope='A complete bounded-height prefix on a distinct discriminant55 lattice. No reflectivity or finite-volume conclusion from this bounded run.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text());print('PASS saved beta55 prefix replay')
 else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['counts'])
