#!/usr/bin/env python3
"""Independent complete support-cost census and continuous-axis gate replay.

No uniform-gate producer or producer field/interval imports. Physical data use
an earlier independent sparse radical-tower review engine, explicitly bound.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from math import gcd,isqrt
from pathlib import Path
from fresh_h20_bridge_edge_closure_20261014_review import Field,parse
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_uniform_chain_cap_gate_20261016_review'
PRIMARY='results/fresh_h20_uniform_chain_cap_gate_20261016.json'
CHAIN='results/fresh_h20_s61_multicap_20261014_exact.json'
CHAIN_REVIEW='results/fresh_h20_s61_multicap_20261014_exact_review.json'
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json','results/fresh_h20_e7_bridge_face_20261014_exact.json','results/fresh_h20_bridge_edge_closure_20261014.json']
RAD=(1,2,3,5,6,10,15,30)
PRIMARY_RAD=(1,2,3,6,5,10,15,30)
BITS=320;DEN=1<<BITS

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d['sources'].items():assert sha(f)==h,f
 return d
def inv(A):
 # Reverse-coordinate rational Gauss-Jordan reconstruction.
 n=len(A);order=list(reversed(range(n)))
 B=[[F(A[i][j])for j in order]+[F(k==l)for l in range(n)]for k,i in enumerate(order)]
 for k in range(n):
  i=next(i for i in range(k,n)if B[i][k]);B[k],B[i]=B[i],B[k]
  q=B[k][k];B[k]=[x/q for x in B[k]]
  for i in range(n):
   if i!=k:
    q=B[i][k]
    if q:B[i]=[x-q*y for x,y in zip(B[i],B[k])]
 return [[B[n-1-i][n+n-1-j]for j in range(n)]for i in range(n)]
def pd(A):
 B=[list(row)for row in A];p=[]
 for k in reversed(range(len(B))):
  q=B[k][k];assert q>0;p.append(q)
  for i in range(k):
   for j in range(k):B[i][j]-=B[i][k]*B[k][j]/q
 return p
class I:
 __slots__=('lo','hi')
 def __init__(self,a=0,b=None):
  if isinstance(a,I):self.lo,self.hi=a.lo,a.hi
  else:self.lo,self.hi=F(a),F(a if b is None else b)
  assert self.lo<=self.hi
 def __add__(self,x):x=I(x);return I(self.lo+x.lo,self.hi+x.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,x):return self+-I(x)
 def __rsub__(self,x):return I(x)+-self
 def __mul__(self,x):
  x=I(x);v=[a*b for a in(self.lo,self.hi)for b in(x.lo,x.hi)];return I(min(v),max(v))
 __rmul__=__mul__
 def __truediv__(self,x):
  x=I(x);assert x.lo>0 or x.hi<0
  return self*I(1/x.hi,1/x.lo)
 def sqrt(self):
  assert self.lo>0
  a=isqrt(self.lo.numerator*DEN*DEN//self.lo.denominator)
  b=isqrt(self.hi.numerator*DEN*DEN//self.hi.denominator)+1
  assert F(a*a,DEN*DEN)<=self.lo and F(b*b,DEN*DEN)>self.hi
  return I(F(a,DEN),F(b,DEN))
 def row(self):return[str(self.lo),str(self.hi)]
@lru_cache(None)
def sqrt_interval(d):
 a=isqrt(d*DEN*DEN);return F(a,DEN),F(a+int(a*a!=d*DEN*DEN),DEN)
def sign(q):
 if not any(q.values()):return 0
 lo=hi=F(0)
 for d,c in q.items():
  a,b=sqrt_interval(d);lo+=c*(a if c>=0 else b);hi+=c*(b if c>=0 else a)
 assert lo>0 or hi<0,q
 return 1 if lo>0 else-1
def mul(a,b):
 out={}
 for i,x in a.items():
  for j,y in b.items():
   g=gcd(i,j);r=i*j//(g*g);out[r]=out.get(r,F(0))+g*x*y
 return{d:x for d,x in out.items()if x}
def difference(q,c):
 out=dict(q);out[1]=out.get(1,F(0))-c;return out

def census(J,S):
 values={3:{1:F(1)},4:{2:F(1)},5:{1:F(1,2),5:F(1,2)},6:{3:F(1)}}
 products={(a,b):mul(values[a],values[b])for a in values for b in values}
 records=[];counts=Counter();mins=[];minimum=None
 # Different traversal and diagonal-plus-unordered-pair evaluation.
 for k in(3,2,1):
  for ids in reversed(list(combinations(range(18),k))):
   for labels in product((6,5,4,3),repeat=k):
    q={}
    def add_scaled(poly,c):
     for d,x in poly.items():q[d]=q.get(d,F(0))+c*x
    for i,m in zip(ids,labels):add_scaled(products[m,m],J[i][i])
    for a,b in combinations(range(k),2):add_scaled(products[labels[a],labels[b]],2*J[ids[a]][ids[b]])
    side=sign(difference(q,2));counts[{-1:'below_two',0:'equal_two',1:'above_gap'}[side]]+=1
    if side>0:
     assert sign(difference(q,F(137,66)))>=0
     if sign(difference(q,F(137,66)))==0:mins.append(dict(walls=[S[i]for i in ids],labels=list(labels)))
    records.append((ids,labels,q))
 digest=sha256()
 for ids,labels,q in sorted(records,key=lambda r:(len(r[0]),r[0],r[1])):
  digest.update(json.dumps([list(ids),list(labels),[str(q.get(d,F(0)))for d in PRIMARY_RAD]],separators=(',',':')).encode())
 assert len(records)==54744 and counts==dict(below_two=97,equal_two=8,above_gap=54639)
 assert len(mins)==8 and all(x['labels']==[3,3,3]for x in mins)
 expected={tuple(x)for x in product((1,10),(17,31),(61,87))}
 assert {tuple(x['walls'])for x in mins}==expected
 return dict(profiles=len(records),counts=dict(counts),minimum_above_two='137/66',minimizing_profiles=sorted(mins,key=lambda x:x['walls']),ordered_profile_digest=digest.hexdigest())

def run():
 primary,chain,chain_review=map(load,(PRIMARY,CHAIN,CHAIN_REVIEW))
 old,bridge,edge=map(load,FILES)
 assert primary['status']=='PASS_ALL_LABEL_UNIFORM_FINITE_CHAIN_ANGLE_CAP_OBSTRUCTION'
 assert chain_review['sources'][CHAIN]==sha(CHAIN)
 ff=Field([parse(old['field']['discriminant'],0),parse(bridge['field']['new_positive_sqrt_radical'],1),parse(edge['field']['new_positive_sqrt_radical'],2)])
 c,add,mul0,sc,total=ff.const,ff.add,ff.mul,ff.scale,ff.total
 zero=ff.zeros[3]
 Q=[list(map(F,row))for row in old['form']];sparse=[[(j,q)for j,q in enumerate(row)if q]for row in Q];cache={}
 def vp(raw,k):return tuple(ff.lift(parse(x,k),k)for x in raw)
 @lru_cache(None)
 def dot(v,w):
  if w not in cache:cache[w]=tuple(total(sc(w[j],q)for j,q in row)for row in sparse)
  return total(mul0(x,y)for x,y in zip(v,cache[w]))
 def rational(x):
  for k in(3,2,1):assert x[1]==ff.zeros[k-1];x=x[0]
  assert not any(x[1:]);return x[0]
 def interval(x):return I(*ff.interval(x,3,BITS))
 S=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,61,87]
 assert chain['support']==primary['support']==edge['recipe']['support']==S
 N=[vp(v,1)for v in old['normals']]
 G=[[rational(dot(N[i],N[j]))for j in S]for i in S]
 assert all(G[i][i]==2 for i in range(18))
 assert G==[list(map(F,row))for row in chain['support_gram']]
 piv=pd(G);J=inv(G)
 assert J==[list(map(F,row))for row in chain['support_inverse']]
 assert all(x>=0 for row in J for x in row)
 diags=[J[i][i]for i in range(18)]
 assert sorted(diags)[:4]==[F(1,2),F(1,2),F(2,3),F(2,3)]
 assert [S[i]for i,x in enumerate(diags)if x==F(1,2)]==[17,31]
 for wall in(17,31):
  i=S.index(wall);assert J[i]==[F(1,2)if i==j else F(0)for j in range(18)]
 # Literal graph components, independent of the saved type string.
 todo=set(range(18));types=[]
 while todo:
  seen={min(todo)};stack=list(seen)
  while stack:
   i=stack.pop()
   for j in todo-seen:
    if G[i][j]:seen.add(j);stack.append(j)
  todo-=seen;degrees=[sum(G[i][j]!=0 for j in seen if j!=i)for i in seen]
  assert max(degrees)<=2 and all(G[i][j]in(0,-1)for i in seen for j in seen if i!=j)
  types.append(len(seen))
 assert sorted(types)==[1,1,2,4,10]
 cc=census(J,S)
 assert cc['ordered_profile_digest']==primary['finite_support_cost_census']['ordered_profile_digest']
 assert cc['counts']==primary['finite_support_cost_census']['counts']
 print('PASS independent all54744 support profiles',cc['counts'],flush=True)
 # All-label and all-support reduction, including continuous nonfinite cases.
 ceiling=F(141,70);z=F(31,25);poly=z**3+z*z-2*z-1
 bounds={'four_contacts':sum(sorted(diags)[:4]),'large_nonisolated':F(2,3)*F(81,25),'large_isolated_plus_contact':F(81,50)+F(1,2),'nonfinite_nonisolated':F(8,3),'nonfinite_isolated_plus_contact':F(5,2),'finite_gap':F(137,66)}
 assert poly<0 and all(x>ceiling for x in bounds.values())
 # Reconstruct X from the previous edge recipe, not the primary gate.
 X=vp(edge['recipe']['X'],2);A=ff.lift(parse(edge['recipe']['A'],2),2)
 assert X==vp(chain['frame']['X'],2) and A==ff.lift(parse(chain['frame']['A'],2),2)
 assert dot(X,X)==sc(A,-1) and ff.sign(A)>0
 assert all(dot(X,N[i])==zero for i in S)
 alo,ahi=map(F,chain['parameter']['isolating_interval']);assert ahi>alo>0 and alo*alo>280
 tt=I(1/ahi,2/alo);dd=(I(4)+2*tt*tt).sqrt()
 assert 2+4/(alo*alo)<ceiling
 axis=[]
 for wall,block in[(17,83),(31,11)]:
  q=rational(dot(N[block],N[block]));assert q==F(old['norms'][block])>0
  den=I(q/2).sqrt()
  xp=interval(dot(X,N[block]))/(interval(A).sqrt()*den)
  wp=interval(dot(N[wall],N[block]))/(2*den)
  result=tt*xp-dd*wp;assert result.lo>0
  axis.append(dict(wall=wall,blocker=block,old_norm=str(q),time_pairing=xp.row(),dual_pairing=wp.row(),uniform_positive_pairing=result.row()))
 aa=I(alo,ahi);a2=aa*aa
 co=a2/(a2+2);si=2*(a2+1).sqrt()/(a2+2);area=2*si*(1-co)
 assert co.lo>0 and co.hi<1 and si.lo>0 and area.lo>0
 print('PASS separate physical-frame arithmetic and continuous blockers',[(x['wall'],x['blocker'])for x in axis],flush=True)
 paths=[f'scripts/{TAG}.py','scripts/fresh_h20_uniform_chain_cap_gate_20261016.py',PRIMARY,CHAIN,CHAIN_REVIEW]+FILES+['scripts/fresh_h20_bridge_edge_closure_20261014_review.py','scripts/fresh_h20_bridge_faces_20261014_review.py']
 return dict(status='PASS_INDEPENDENT_ALL_LABEL_UNIFORM_CHAIN_CAP_GATE',sources={p:sha(p)for p in paths},
  independent_support=dict(walls=S,gram=[[str(x)for x in row]for row in G],inverse=[[str(x)for x in row]for row in J],reverse_LDL_pivots=list(map(str,piv)),component_ranks=sorted(types)),
  finite_cost_census=cc,all_label_lower_bounds={k:str(v)for k,v in bounds.items()},label7_cubic_at_31_over25=str(poly),
  chain_parameter_interval=[str(alo),str(ahi)],whole_time_interval=tt.row(),whole_nonfinite_contact_interval=dd.row(),noncollinearity_bound=area.row(),axis_blockers=axis,
  scope='Independent squarefree-dictionary Fraction evaluation of every54744 support profile, reversed rational support inversion, exact analytic all-label reduction, and both whole-interval physical blockers. Physical radical arithmetic reuses the earlier independent bridge-edge review engine, not the producer Field or interval classes. The complete133-wall construction and parameter isolation are inherited from their frozen independent review; no full133 replay. Exclusion is only one common finite chain angle in a coherent acute extension, not arbitrary H20 completions.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved independent uniform-chain gate replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
