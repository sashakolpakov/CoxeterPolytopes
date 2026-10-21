#!/usr/bin/env python3
"""Self-contained exact rank11 obstruction to the first star label state.

No probe-domain dependency. Every finite-contact bound is guarded by a
literal actual-vertex Gram failure or an explicit positive Schur budget.
"""
import argparse,json
from itertools import combinations,product
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_cube_metric_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json';STATE='results/fresh_h9_tree_cube_wedges_20261014_case_00.json'
Z=(F(0),F(0));O=(F(1),F(0));H=(F(1,2),F(0));A=(F(0),F(1,2))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return-x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def mul(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def inv(x):
 q=x[0]*x[0]-2*x[1]*x[1];assert q;return x[0]/q,-x[1]/q
def div(x,y):return mul(x,inv(y))
def scale(x,s):return x[0]*s,x[1]*s
def total(xs):
 q=Z
 for x in xs:q=add(q,x)
 return q
def sign(x):
 a,b=x;s=lambda v:int(v>0)-int(v<0)
 if not b:return s(a)
 if not a:return s(b)
 if s(a)==s(b):return s(a)
 return s(a)*s(a*a-2*b*b)
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(tuple,list)):return[enc(v)for v in x]
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 return x
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def project(i,extra):
 x=[Z]*8;x[2*i]=x[2*i+1]=neg(H)
 if extra is not None:x[extra]=neg(H)
 return tuple(x)
def pdot(x,y):
 # Inverse of each normalized B2 block: [[2,sqrt2],[sqrt2,2]].
 return total(add(scale(add(mul(x[k],y[k]),mul(x[k+1],y[k+1])),2),
  mul((F(0),F(1)),add(mul(x[k],y[k+1]),mul(x[k+1],y[k]))))for k in range(0,8,2))
def pivots(M):
 M=[r[:]for r in M];out=[]
 for k in range(len(M)):
  p=M[k][k];out.append(p)
  if sign(p)<=0:return out
  for i in range(k+1,len(M)):
   for j in range(k+1,len(M)):M[i][j]=sub(M[i][j],div(mul(M[i][k],M[k][j]),p))
 return out
def determinant3(M):
 a,b,c=M[0];d,e,f=M[1];g,h,i=M[2]
 return add(sub(mul(a,sub(mul(e,i),mul(f,h))),mul(b,sub(mul(d,i),mul(f,g)))),mul(c,sub(mul(d,h),mul(e,g))))
def run():
 geo=json.loads((ROOT/GEOM).read_text());state=json.loads((ROOT/STATE).read_text());g=geo['representatives'][0]['geometry']
 for d in [geo,state]:
  for p,h in d['sources'].items():assert sha(p)==h,p
 assert state['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
 V=[v['incident']for v in g['vertices']];S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
 labels={(i,j):m for i,j,m in state['solution']['pair_labels']}
 assert len(labels)==20
 assert all(labels[tuple(sorted((v,S[2*k+j])))]==3 for k in range(4)for j in range(2)for v in[P[k],C[k]])
 assert all(labels[S[2*k],S[2*k+1]]==4 for k in range(4))
 COS={2:Z,3:H,4:A}
 def gram(v,lab):return[[O if i==j else neg(COS[lab.get(tuple(sorted((i,j))),2)])for j in v]for i in v]
 def failure(pair,m):
  lab=labels|{pair:m}
  for k,v in enumerate(V):
   if not set(pair)<=set(v):continue
   ps=pivots(gram(v,lab))
   if sign(ps[-1])<=0:return dict(pair=pair,test_label=m,actual_vertex=k,leading_order=len(ps),last_pivot=ps[-1])
  raise AssertionError(('no exact upper-bound witness',pair,m))
 cross=[tuple(sorted((a,b)))for i,j in combinations(range(4),2)for a in S[2*i:2*i+2]for b in S[2*j:2*j+2]]
 forced=[failure(p,3)for p in cross+list(zip(P,C))]
 assert len(forced)==28
 SB=gram(S,labels);assert all(sign(q)>0 for q in pivots(SB))
 SI=[[pdot(tuple(O if k==i else Z for k in range(8)),tuple(O if k==j else Z for k in range(8)))for j in range(8)]for i in range(8)]
 assert all(total(mul(SB[i][k],SI[k][j])for k in range(8))==(O if i==j else Z)for i in range(8)for j in range(8))
 budgets=[]
 for i in range(4):
  for v in[P[i],C[i]]:
   support=sorted((set(S)-{S[2*i]})|{v});assert len(support)==8
   k=next(j for j,w in enumerate(V)if set(support)<=set(w))
   budgets.append(dict(leaf=v,actual_vertex=k,support=support,own_contact_cost=F(1,4),
    foreign_diagonal_inverse=2,two_minimal_foreign_cost=F(5,4),single_label4_cost=F(5,4)))
 assert F(1,4)+2*F(1,4)*2==F(5,4)>1
 assert F(1,4)+2*F(1,2)==F(5,4)>1
 # Any foreign nonright contact is consequently a single label3.
 patterns=[[project(i,e)for e in[None]+[j for j in range(8)if j//2!=i]]for i in range(4)]
 assert all(sub(O,pdot(p,p))in[neg(A),neg(add(H,A))]for arr in patterns for p in arr)
 # C0 meets both final anchors. Label4 already fails, hence all larger
 # labels do too by monotonicity. No finite-label ceiling is assumed.
 upper=[failure(tuple(sorted((C[0],j))),4)for j in[P[3],C[3]]]
 rows=[];anchors=[];signs={-1:0,1:0}
 for ia,ib in product(range(7),repeat=2):
  aa,bb=patterns[3][ia],patterns[3][ib];qa=sub(O,pdot(aa,aa));qb=sub(O,pdot(bb,bb));ab=neg(pdot(aa,bb));ad=sub(mul(qa,qb),mul(ab,ab));assert sign(ad)<0
  anchors.append(dict(patterns=[ia,ib],residual_gram=[[qa,ab],[ab,qb]],determinant=ad))
  for ic,cc in enumerate(patterns[0]):
   qc=sub(O,pdot(cc,cc))
   for ma,mb in product([2,3],repeat=2):
    ac=sub(neg(COS[ma]),pdot(aa,cc));bc=sub(neg(COS[mb]),pdot(bb,cc))
    M=[[qa,ab,ac],[ab,qb,bc],[ac,bc,qc]];det=determinant3(M);assert sign(det)!=0
    signs[sign(det)]+=1
    rows.append(dict(anchor_patterns=[ia,ib],third_pattern=ic,contact_labels=[ma,mb],residual_determinant=det))
 assert len(rows)==1372 and len(anchors)==49
 paths=[GEOM,STATE,f'scripts/{TAG}.py']
 return enc(dict(status='PASS_FIRST_STAR_MANDATORY_ASSIGNMENT_RANK_AT_LEAST_ELEVEN',sources={p:sha(p)for p in paths},
  fixed_positive_support=S,positive_leaves=P,cut_leaves=C,positive_gram=SB,positive_inverse=SI,
  forced_right_witnesses=forced,foreign_contact_budget_witnesses=budgets,third_contact_upper_witnesses=upper,
  projection_patterns=patterns,anchors=anchors,determinant_sign_counts=signs,rank_obstructions=rows,
  scope='First saved star mandatory assignment only, not the star incidence. Every omitted contact bound is proved from an actual-vertex Gram and monotonicity or an explicit Schur budget. The PD8 block leaves at most a two-dimensional residual. All49 anchor patterns have negative determinant, so none is degenerate. In all1372 possible third-pattern/contact cases the residual3x3 determinant is nonzero; the eleven-normal full Gram therefore has rank11. Arbitrary dotted distances and the two central normals cannot repair this principal-rank obstruction.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['rank_obstructions']),d['determinant_sign_counts'],flush=True)
