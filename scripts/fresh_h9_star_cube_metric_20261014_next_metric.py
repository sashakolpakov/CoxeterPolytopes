#!/usr/bin/env python3
"""Exact metric gate for the next orientation-search mandatory assignment."""
import argparse,json
from itertools import product
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from sympy import QQ,sqrt
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_cube_metric_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json';FIRST=f'results/{TAG}.json';STATE=f'results/{TAG}_orientation_search_case_00.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 first,src=[json.loads((ROOT/p).read_text())for p in[FIRST,GEOM]]
 for d in[first,src]:
  for p,h in d['sources'].items():assert sha(p)==h,p
 assert first['status']=='PASS_FIRST_STAR_MANDATORY_ASSIGNMENT_RANK_AT_LEAST_ELEVEN'
 g=src['representatives'][0]['geometry'];V=[z['incident']for z in g['vertices']]
 S=first['fixed_positive_support'];P=first['positive_leaves'];C=first['cut_leaves'];guards=[]
 for i in range(4):
  for j in range(2):
   root=S[2*i+j];support=sorted((set(S)-{S[2*i+1-j]})|{C[i]})
   vi=next(k for k,v in enumerate(V)if set(support)<=set(v))
   tri=sorted([P[i],root,C[i]]);ti=next(k for k,v in enumerate(V)if set(tri)<=set(v))
   guards.append(dict(cut=C[i],increased_own_root=root,eight_support=support,actual_vertex=vi,
    label6_path=tri,path_actual_vertex=ti,label6_determinant='1-1/4-3/4=0'))
 assert len(guards)==8 and F(1,2)+2*F(1,4)==1 and 1-F(1,4)-F(3,4)==0
 K=QQ.algebraic_field(sqrt(2),sqrt(5));z=K.zero;o=K.one;h=K.convert(QQ(1,2));r2=K.from_sympy(sqrt(2))
 co={3:h,4:r2*h,5:K.from_sympy((1+sqrt(5))/4)}
 def enc(x):return[str(q)for q in x.to_list()]
 def pdot(x,y):return sum((2*(x[k]*y[k]+x[k+1]*y[k+1])+r2*(x[k]*y[k+1]+x[k+1]*y[k])for k in range(0,8,2)),z)
 def patterns(group,profile):
  allowed=[None]+[k for k in range(8)if k//2!=group]if profile==(3,3)else[None]
  ans=[]
  for extra in allowed:
   p=[z]*8;p[2*group]=-co[profile[0]];p[2*group+1]=-co[profile[1]]
   if extra is not None:p[extra]=-h
   ans.append((extra,p))
  return ans
 def det3(a,b,c,d,e,f):return a*(d*f-e*e)-b*(b*f-c*e)+c*(b*e-c*d)
 # Recheck all domain witnesses directly: C3's new right contact is not
 # a monotone change of the first certificate.
 state=json.loads((ROOT/STATE).read_text())
 for p,hh in state['sources'].items():assert sha(p)==hh,p
 labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
 expected={(a,b):m for a,b,m in json.loads((ROOT/'results/fresh_h9_tree_cube_wedges_20261014_case_00.json').read_text())['solution']['pair_labels']}
 expected[11,13]=4;expected[12,17]=5;expected[13,17]=2;assert labels==expected
 def sign2(a,b):
  sg=lambda t:int(t>0)-int(t<0)
  if not b:return sg(a)
  if not a:return sg(b)
  if sg(a)==sg(b):return sg(a)
  return sg(a)*sg(a*a-2*b*b)
 def sign(x):
  pp=[F(str(t))for t in x.to_list()];pp=[F(0)]*(4-len(pp))+pp
  a3,a2,a1,a0=pp;a=a0+7*a2;b=a1+17*a3;c=a1+11*a3;d=2*a2
  su=sign2(a,b);sv=sign2(c,d)
  if not sv:return su
  if not su:return sv
  if su==sv:return su
  return su*sign2(a*a+2*b*b-5*c*c-10*d*d,2*a*b-10*c*d)
 assert sign(r2)>0 and K.mod.to_list()==[1,0,-14,0,9]
 co[2]=z
 def gram(v,lab):return [[o if a==b else -co[lab.get(tuple(sorted((a,b))),2)]for b in v]for a in v]
 def pivots(M):
  M=[row[:]for row in M];out=[]
  for k in range(len(M)):
   q=M[k][k];out.append(q)
   if sign(q)<=0:break
   for i in range(k+1,len(M)):
    for j in range(k+1,len(M)):M[i][j]-=M[i][k]*M[k][j]/q
  return out
 def failure(pair,m):
  for vi,v in enumerate(V):
   if not set(pair)<=set(v):continue
   pp=pivots(gram(v,labels|{pair:m}))
   if sign(pp[-1])<=0:return dict(pair=pair,test_label=m,actual_vertex=vi,leading_order=len(pp),last_pivot=enc(pp[-1]))
  raise AssertionError(('no failure',pair,m))
 from itertools import combinations
 cross=[tuple(sorted((a,b)))for i,j in combinations(range(4),2)for a in S[2*i:2*i+2]for b in S[2*j:2*j+2]]
 forced=[failure(p,3)for p in cross+list(zip(P,C))]
 upper=[failure(tuple(sorted((C[0],j))),4)for j in[P[3],C[3]]]
 assert sign(co[5]**2+h-o)>0
 p3=patterns(3,(3,4));rows=[];zeros=[];distinct=set();total=0;anchorzeros=0
 left=(3,3);right=(5,2);cases=[]
 for ia,x in p3:
  for ib,y in patterns(3,right):
   qx=o-pdot(x,x);qy=o-pdot(y,y);xy=-pdot(x,y);ad=qx*qy-xy*xy
   if ad==z:anchorzeros+=1
   for ic,v in patterns(0,left):
    qv=o-pdot(v,v)
    for ma,mb in product([2,3],repeat=2):
     xv=-(z if ma==2 else h)-pdot(x,v);yv=-(z if mb==2 else h)-pdot(y,v)
     det=det3(qx,xy,xv,qy,yv,qv);total+=1;code=enc(det);distinct.add(tuple(code))
     record=dict(foreign_contacts=[ia,ib,ic],anchor_labels=[ma,mb],residual_determinant=code)
     cases.append(record)
     if det==z:zeros.append(record)
 rows.append(dict(C0_profile=left,C3_profile=right,cases=cases))
 assert total==28
 paths=[GEOM,FIRST,STATE,f'scripts/{TAG}_next_metric.py']
 return dict(status='PASS_NEXT_ORIENTATION_ASSIGNMENT_RANK_ELEVEN'if not zeros else'EXACT_NEXT_ORIENTATION_RANK_SURVIVORS',
  sources={p:sha(p)for p in paths},own_contact_guards=guards,forced_right_witnesses=forced,anchor_label_upper_witnesses=upper,
  field='Q(sqrt2,sqrt5), primitive generator sqrt2+sqrt5; descending rational coefficient lists',
  residual_cases=total,anchor_zero_cases=anchorzeros,records=rows,surviving_rank_cases=zeros,
  scope='Exactly the source-bound orientation-search mandatory assignment. All30 right/contact-bound witnesses are recomputed without monotonicity across the newly right contact. C3 has no foreign nonright contact since cos(pi/5)^2+1/2>1 on an actual spherical eight-support. P3 also has no foreign contact because its increased own label4 uses cost1/2; any foreign label3 would add1/2. C0 has seven patterns. Thus the two anchors are fixed and all28 third-leaf/contact cases are covered. Every actual rank-ten realization must have zero residual determinant. All omitted and dotted contacts remain free.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_next_metric.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['residual_cases'],len(d['surviving_rank_cases']),flush=True)
