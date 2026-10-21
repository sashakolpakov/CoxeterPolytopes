#!/usr/bin/env python3
"""Complete cut-contact metric gate with fixed positive leaves and B2 blocks."""
import argparse,json
from itertools import product
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from sympy import QQ,sqrt
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_cube_metric_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json';FIRST=f'results/{TAG}.json';STATE=f'results/{TAG}_family_search_case_00.json'
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
 # Relax every cut contact to right except the two explicitly enumerated
 # profiles. This avoids transferring witnesses through another cut's
 # decreased contact.
 labels={}
 for i in range(4):
  labels[S[2*i],S[2*i+1]]=4
  for j in range(2):labels[tuple(sorted((P[i],S[2*i+j])))]=3
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
 forced=[failure(p,3)for p in cross]
 profiles=[pr for pr in product(range(2,6),repeat=2)if F(1,4)+F(1,pr[0])+F(1,pr[1])<1]
 assert len(profiles)==11
 own=[]
 for i in range(4):
  for profile in profiles:
   j=next(j for j,m in enumerate(profile)if m>=3)
   tri=sorted([P[i],C[i],S[2*i+j]]);vi=next(k for k,v in enumerate(V)if set(tri)<=set(v))
   witness=o-2*co[profile[j]];assert sign(witness)<=0
   own.append(dict(group=i,profile=profile,three_support=tri,actual_vertex=vi,
    tested_PiCi_label=3,ones_quadratic=enc(witness)))
 assert sign(co[5]**2+h-o)>0
 base=labels.copy();p3=patterns(3,(3,3));rows=[];zeros=[];distinct=set();total=0;anchorzeros=0;upper=[]
 for left,right in product(profiles,repeat=2):
  labels=base.copy()
  for i,profile in [(0,left),(3,right)]:
   for j in range(2):labels[tuple(sorted((C[i],S[2*i+j])))]=profile[j]
  upper.append(dict(C0_profile=left,C3_profile=right,witnesses=[failure(tuple(sorted((C[0],j))),4)for j in[P[3],C[3]]]))
  cases=[]
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
      if det==z:zeros.append(dict(C0_profile=left,C3_profile=right,**record))
  rows.append(dict(C0_profile=left,C3_profile=right,cases=cases))
 assert total==8092
 affine=[];survivors=[]
 for row in zeros:
  labels=base.copy()
  for i,profile in [(0,row['C0_profile']),(3,row['C3_profile'])]:
   for j in range(2):labels[tuple(sorted((C[i],S[2*i+j])))]=profile[j]
  for leaf,extra in zip([P[3],C[3],C[0]],row['foreign_contacts']):
   if extra is not None:labels[tuple(sorted((leaf,S[extra])))]=3
  for anchor,m in zip([P[3],C[3]],row['anchor_labels']):labels[tuple(sorted((C[0],anchor)))]=m
  for vi,v in enumerate(V):
   pp=pivots(gram(v,labels))
   if sign(pp[-1])<=0:
    affine.append(dict(rank_survivor=row,actual_vertex=vi,leading_order=len(pp),last_pivot=enc(pp[-1])));break
  else:survivors.append(row)
 paths=[GEOM,FIRST,f'scripts/{TAG}_all_cut_profiles.py']
 return dict(status='PASS_ALL_CUT_CONTACT_PROFILES_EXCLUDED'if not survivors else'EXACT_ALL_CUT_CONTACT_SURVIVORS',
  sources={p:sha(p)for p in paths},own_contact_guards=guards,forced_cross_right_witnesses=forced,
  own_pair_right_witnesses=own,profiles=profiles,anchor_label_upper_witnesses=upper,
  field='Q(sqrt2,sqrt5), primitive generator sqrt2+sqrt5; descending rational coefficient lists',
  residual_cases=total,anchor_zero_cases=anchorzeros,records=rows,zero_rank_cases=zeros,
  actual_vertex_exclusions=affine,surviving_cases=survivors,
  scope='Star case0: B2 internal labels4 and all positive leaves ownlabels3/3 are fixed, but all cut ownlabels are free. PiCi must be right: every Lanner cut triangle has at least one ownlabel>=3, and testing PiCi=3 gives a nonpositive all-ones quadratic on its actual three-face. Then the actual 3-m paths exclude each ownlabel>=6. Hyperbolicity of (4,a,b) leaves the eleven ordered profiles. Both actual eight-support orientations restrict (3,3) cuts to seven projection patterns; every other cut profile has no foreign contact. The C0-anchor label bounds are recomputed after relaxing all unused cut contacts to right. All8092 residual cases are exact. Actual finite-vertex tests additionally reject rank survivors. Cut profiles at C1/C2 cannot alter these necessary principal tests. This is a compact-realization obstruction, not an exclusion of changed B2 or positive-leaf labels.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_all_cut_profiles.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['residual_cases'],len(d['zero_rank_cases']),len(d['surviving_cases']),flush=True)
