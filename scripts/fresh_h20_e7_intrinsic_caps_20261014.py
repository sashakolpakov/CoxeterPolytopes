#!/usr/bin/env python3
"""Exact H13-source transplant after four-wall backtracking in the87 chamber.

All arithmetic is Fraction in Q(sqrt2,sqrt79,sqrt11), reusing the earlier
independent nested-quadratic arithmetic.  The face census is rational.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h20_e8_interior_caps_20261014_exact_review as K

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_intrinsic_caps_20261014'
BASE='results/fresh_h20_e8_interior_caps_20261014_exact.json'
FACE='results/fresh_h20_e7_face_20261014.json'
SOURCE='results/fresh_h13_a11a2_completion_20261014.json'
Z=K.Z
SQ2=(F(0),F(1))+(F(0),)*6
SQ22=(F(0),)*5+(F(1),F(0),F(0))

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(tuple,list)):return[enc(v)for v in x]
 return x

def rank_null(a):
 """Rational RREF and a complete nullspace basis."""
 a=[list(row)for row in a];m=len(a);n=len(a[0]);r=0;piv=[]
 for j in range(n):
  hits=[i for i in range(r,m)if a[i][j]]
  if not hits:continue
  i=hits[0];a[r],a[i]=a[i],a[r];q=a[r][j];a[r]=[v/q for v in a[r]]
  for i in range(m):
   if i!=r:
    q=a[i][j];a[i]=[v-q*w for v,w in zip(a[i],a[r])]
  piv.append(j);r+=1
  if r==m:break
 out=[]
 for j in range(n):
  if j in piv:continue
  v=[F(0)]*n;v[j]=F(1)
  for i,k in enumerate(piv):v[k]=-a[i][j]
  out.append(v)
 return r,out

def run():
 d,f,source=[json.loads((ROOT/p).read_text())for p in [BASE,FACE,SOURCE]]
 for doc in [d,f,source]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 Q=[[F(x)for x in row]for row in d['form']];old=list(map(K.vector,d['normals']))
 H=K.vector(d['strict_center']);escape=K.vector(d['strict_escape']['point'])
 sparse=[[(j,x)for j,x in enumerate(row)if x]for row in Q];cache={};dots={}
 def dot(x,y):
  key=(x,y)
  if key in dots:return dots[key]
  if y not in cache:cache[y]=tuple(tuple(sum(q*y[j][k]for j,q in row)for k in range(8))for row in sparse)
  out=Z
  for a,b in zip(x,cache[y]):
   if any(a)and any(b):out=K.add(out,K.multiply(a,b))
  dots[x,y]=dots[y,x]=out;return out
 def lin(vectors,coefficients):
  out=(Z,)*21
  for v,a in zip(vectors,coefficients):out=K.vadd(out,K.vscale(v,a))
  return out
 A=f['face_walls'];S=f['E7_support'];assert A==list(range(11))+[18,83]
 J=K.inverse([[F(x)for x in row]for row in f['face_gram']])
 W=[lin([old[i]for i in A],list(map(K.scalar,row)))for row in J]
 T=K.vector(f['intrinsic_normalization']['timelike_basis'])
 omega=K.vector(f['unique_nonright_meeting_profile']['fundamental_weight'])
 assert dot(old[26],old[26])==K.scalar(2)
 assert dot(T,T)==K.scalar(-2)and dot(omega,omega)==K.scalar(F(3,2))
 assert all(dot(T,old[i])==Z for i in A+S)
 assert all(dot(W[i],old[j])==K.scalar(int(i==k))for i in range(13)for k,j in enumerate(A))
 # All four oriented corner isometries; each facet has both necessary lifts.
 transfer=[];selected=None
 for reverse in [False,True]:
  for swap in [False,True]:
   rpos=0 if reverse else 10;spos=8 if reverse else 2;apos=11 if swap else 12
   R=lin([W[rpos],W[apos],T],[K.scalar(-1),K.neg(SQ2),K.scale(SQ2,F(1,4))])
   U=lin([W[spos],T],[K.scalar(-1),K.scale(SQ2,F(1,4))])
   assert dot(R,R)==dot(U,U)==K.scalar(2)and dot(R,U)==Z
   assert dot(R,old[26])==K.sub(SQ2,K.scale(SQ22,F(1,4)))
   assert K.sign(dot(R,old[26]))>0
   assert dot(U,old[26])==K.scale(SQ22,F(-1,4))
   assert K.multiply(dot(U,old[26]),dot(U,old[26]))==K.scalar(F(11,8))
   assert F(1,4)<F(11,32)<F(1,2)
   for name,y,pos in [('R',R,rpos),('S',U,spos)]:
    for lift in [0,1]:
     r=K.vadd(y,K.vscale(omega,K.scalar(-2*lift)));q=dot(r,r)
     assert q==K.scalar(2 if not lift else 8)
     if lift:
      assert dot(r,old[A[pos]])==K.scalar(-1)
      blocker=dict(wall=A[pos],pairing=K.scalar(-1),squared_cosine=F(1,16),reason='strictly_between_right_and_pi3')
     else:
      blocker=dict(wall=26,pairing=dot(r,old[26]),reason='positive_pairing'if name=='R'else'squared_cosine_between_pi3_and_pi4')
     transfer.append(dict(A11_reversed=reverse,A2_swapped=swap,facet=name,nonright_lift=bool(lift),norm=q,blocker=blocker))
   if not reverse and swap:selected=[R,U]
 assert len(transfer)==16 and selected is not None
 removed=[17,26,75,76];retained=[i for i in range(87)if i not in removed]
 old_to_new={i:j for j,i in enumerate(retained)}
 N=[old[i]for i in retained]+selected;assert len(N)==85
 norms=[K.rational(dot(n,n))for n in N]
 counts=Counter();newpairs=[]
 for i in [83,84]:
  assert norms[i]==2
  for j in range(i):
   g=dot(N[i],N[j]);assert K.sign(g)<=0,(i,j,g)
   square=K.multiply(g,g);prod=norms[i]*norms[j]
   if g==Z:kind='right'
   elif square==K.scalar(prod/4):kind='pi3'
   elif square==K.scalar(prod/2):kind='pi4'
   else:assert K.sign(K.sub(square,K.scalar(prod)))>0,(i,j,g);kind='dotted'
   counts[kind]+=1;newpairs.append(dict(walls=[i,j],pairing=g,relation=kind))
 assert len(newpairs)==167
 # The four deletions are explicit. Every retained old pair was certified
 # in the immutable87 package. Check the full Gram's acuteness anew.
 Gram=[[dot(n,m)for m in N]for n in N]
 assert all(K.sign(Gram[i][j])<=0 for i in range(85)for j in range(i))
 h=[dot(H,n)for n in N];assert K.sign(dot(H,H))<0 and all(K.sign(v)<0 for v in h)
 facets=[]
 for i,r in enumerate(N):
  p=K.vadd(H,K.vscale(r,K.scale(h[i],-1/norms[i])))
  assert K.sign(dot(p,p))<0 and K.sign(dot(p,H))<0
  assert all(dot(p,w)==Z if j==i else K.sign(dot(p,w))<0 for j,w in enumerate(N))
  facets.append(p)
 # A rational face basis: A11, sqrt2*A2, sqrt2*T.  Both source caps
 # have rational coordinates in this basis, whose last square is -4.
 basis=[old[i]for i in A[:11]]+[K.vscale(old[i],SQ2)for i in A[11:]]+[K.vscale(T,SQ2)]
 C=[[K.rational(dot(v,w))for w in basis]for v in basis]
 assert C==[[F(2 if i<11 else 4)if i==j and i<13 else F(-4)if i==j==13 else F(-1)if i<11 and j<11 and abs(i-j)==1 else F(-2)if {i,j}=={11,12}else F(0)for j in range(14)]for i in range(14)]
 spatial=C[:13];spatial=[row[:13]for row in spatial];JI=K.inverse(spatial);assert K.pd(spatial)
 facevectors=basis[:13]+selected
 source_permutation=list(range(11))+[12,11,13,14]
 for row,i in zip(source['normalized_transfer']['normal_coordinates'],source_permutation):
  coefficients=[(F(v[0]),F(v[1]))+(F(0),)*6 for v in row]
  coefficients[11],coefficients[12]=coefficients[12],coefficients[11]
  assert lin([old[j]for j in A]+[T],coefficients)==facevectors[i]
 # Every other retained projection is future timelike, hence redundant
 # throughout the hyperbolic E7 face, not merely at the saved vertices.
 causal_nonfacets=[]
 for i in retained:
  if i in A+S:continue
  projection=K.vector(f['projected_walls'][i])
  assert K.sign(dot(projection,projection))<0 and K.sign(dot(projection,T))<0
  causal_nonfacets.append(i)
 assert len(causal_nonfacets)==63
 reduced=[]
 for v in facevectors:
  pairs=[K.rational(dot(v,b))for b in basis]
  c=[sum(JI[i][j]*pairs[j]for j in range(13))for i in range(13)]+[-pairs[13]/4]
  assert lin(basis,list(map(K.scalar,c)))==v;reduced.append(c)
 assert reduced[:13]==[[F(int(i==j))for j in range(14)]for i in range(13)]
 fcov=[[sum(v[i]*C[i][j]for i in range(14))for j in range(14)]for v in reduced]
 FG=[[sum(a*b for a,b in zip(row,v))for v in reduced]for row in fcov]
 assert all(FG[i][i]>0 for i in range(15))and all(FG[i][j]<=0 for i in range(15)for j in range(i))
 # Positive combination of every15 face normals equals Tc/2: a future
 # timelike dual, giving a bounded projective section with Tc-coordinate1.
 dual=[-reduced[13][i]-reduced[14][i]for i in range(13)]+[F(1),F(1)]
 assert all(v>0 for v in dual)
 assert [sum(dual[i]*reduced[i][j]for i in range(15))for j in range(14)]==[F(0)]*13+[F(1,2)]
 facephysical=[old_to_new[i]for i in A]+[83,84]
 E7=[old_to_new[i]for i in S]
 vertices={};census=[]
 for subset in combinations(range(15),13):
  rank,null=rank_null([fcov[i]for i in subset]);assert rank<=13
  rec=dict(support=list(subset),rank=rank)
  if rank<13:rec['status']='rank_deficient';census.append(rec);continue
  v=null[0]
  if not v[-1]:
   signs=[sum(a*b for a,b in zip(row,v))for row in fcov]
   assert min(signs)<0<max(signs);rec['status']='outside_cone_at_zero_time';census.append(rec);continue
  v=[q/v[-1]for q in v];pair=[sum(a*b for a,b in zip(row,v))for row in fcov]
  if max(pair)>0:rec['status']='infeasible';census.append(rec);continue
  q=sum(v[i]*C[i][j]*v[j]for i in range(14)for j in range(14));assert q<=0
  incident=[i for i,g in enumerate(pair)if not g]
  rec['status']='finite'if q<0 else'ideal';rec['point']=v;census.append(rec)
  vertices[tuple(v)]=dict(point=v,norm=q,incident=incident)
 assert len(census)==105
 finite=[];ideal=[]
 for row in vertices.values():
  v=row['point'];point=lin(basis,list(map(K.scalar,v)));assert dot(point,point)==K.scalar(row['norm'])and K.sign(dot(point,H))<0, (row, dot(point,point), dot(point,H))
  incident=sorted(E7+[facephysical[i]for i in row['incident']])
  assert all(dot(point,n)==Z if i in incident else K.sign(dot(point,n))<0 for i,n in enumerate(N))
  row.update(ambient_point=point,ambient_incident=incident)
  if row['norm']<0:
   assert len(row['incident'])==13 and K.pd([[FG[i][j]for j in row['incident']]for i in row['incident']])
   assert len(incident)==20;finite.append(row)
  else:
   unseen=set(row['incident']);blocks=[]
   while unseen:
    todo=[min(unseen)];block=set()
    while todo:
     i=todo.pop()
     if i in block:continue
     block.add(i);todo.extend(j for j in row['incident']if i!=j and FG[i][j])
    unseen-=block;B=sorted(block);G=[[FG[i][j]for j in B]for i in B]
    rr,nn=rank_null(G);assert rr==len(B)-1 and len(nn)==1
    kernel=nn[0]
    if kernel[0]<0:kernel=[-x for x in kernel]
    assert min(kernel)>0 and K.pd([r[:-1]for r in G[:-1]])
    ray=[sum(a*reduced[i][j]for a,i in zip(kernel,B))for j in range(14)]
    mu=ray[-1];assert mu>0 and ray==[mu*x for x in v]
    blocks.append(dict(incident=B,gram=G,positive_kernel=kernel,null_multiple=mu))
   assert sum(len(b['incident'])-1 for b in blocks)==12
   # E7 supplies a genuine spherical recession in the ambient link;
   # it does not disappear when the intrinsic12-dimensional link closes.
   recession=K.vscale(omega,K.scalar(-1))
   assert dot(recession,point)==Z and dot(recession,recession)==K.scalar(F(3,2))
   values=[dot(recession,N[i])for i in incident]
   assert all(K.sign(g)<=0 for g in values)and any(K.sign(g)<0 for g in values)
   row.update(intrinsic_affine_blocks=blocks,intrinsic_quotient_rank=12,ambient_spherical_recession=recession,recession_norm=F(3,2),recession_pairings=values)
   ideal.append(row)
 assert len(finite)==source['finite_vertex_count']==44
 assert len(ideal)==source['ideal_vertex_count']==2
 # The inherited global escape survives this changed chamber exactly.
 ep=[dot(escape,n)for n in N];assert all(K.sign(v)<0 for v in ep)
 assert dot(escape,escape)==K.field(d['strict_escape']['norm'])and K.sign(dot(escape,H))<0
 paths=[BASE,FACE,SOURCE,f'scripts/{TAG}.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
 return enc(dict(status='PASS_85_WALL_REPLACEMENT_COMPLETES_H13_FACE_GLOBAL_ESCAPE_REMAINS',sources={p:sha(p)for p in paths},rejected_fixed87_transfers=transfer,removed_old87_walls=removed,retained_old87_walls=retained,old_to_new=old_to_new,form=Q,normals=N,norms=norms,new_pair_counts=dict(counts),new_pairs=newpairs,strict_center=H,facet_points=facets,E7_support=E7,face_facet_walls=facephysical,face_basis=basis,face_form=C,face_facet_coordinates=reduced,face_gram=FG,source_normal_permutation=source_permutation,causal_nonfacet_old87_walls=causal_nonfacets,positive_timelike_dual_weights=dual,projective_census=census,face_finite_vertices=finite,face_ideal_vertices=ideal,strict_escape=dict(point=escape,norm=dot(escape,escape),all85_pairings=ep),scope='All4 A11 reversal/A2 swap transfers and both necessary lifts fail when retaining87. Removing exactly4 specified walls permits both right sourcecaps, giving85actualfacets and a complete finite-volume H13 E7 face, certified by all105 projective supports and bounded intrinsic cusp links. Its ideal vertices retain an ambient spherical E7 link factor with explicit recession. The unchanged strict spacelike escape proves infinite ambient volume. No full85 ambient corner/cusp census and no universal completion obstruction.'))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'],'face finite/ideal',len(d['face_finite_vertices']),len(d['face_ideal_vertices']))
