#!/usr/bin/env python3
"""Monotone cut-contact rank gate over Q(sqrt2,sqrt5), no numerical fits."""
import argparse,json
from itertools import product
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from sympy import QQ,sqrt
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_cube_metric_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json';FIRST=f'results/{TAG}.json'
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
 p3=patterns(3,(3,3));profiles=list(product([3,4,5],repeat=2));rows=[];zeros=[];distinct=set();total=0;anchorzeros=0
 for left,right in product(profiles,repeat=2):
  cases=[]
  for ia,x in p3:
   for ib,y in patterns(3,right):
    qx=o-pdot(x,x);qy=o-pdot(y,y);xy=-pdot(x,y);ad=qx*qy-xy*xy
    if ad==z:anchorzeros+=1
    for ic,v in patterns(0,left):
     qv=o-pdot(v,v)
     for ma,mb in product([2,3],repeat=2):
      xv=-(z if ma==2 else h)-pdot(x,v);yv=-(z if mb==2 else h)-pdot(y,v)
      d=det3(qx,xy,xv,qy,yv,qv);total+=1;code=enc(d);distinct.add(tuple(code))
      record=dict(foreign_contacts=[ia,ib,ic],anchor_labels=[ma,mb],residual_determinant=code)
      cases.append(record)
      if d==z:zeros.append(dict(C0_profile=left,C3_profile=right,**record))
  rows.append(dict(C0_profile=left,C3_profile=right,cases=cases))
 assert len(rows)==81 and total==6300
 affine=[]
 for row in zeros:
  ia,ib,ic=row['foreign_contacts'];ma,mb=row['anchor_labels']
  x=dict(p3)[ia];v=dict(patterns(0,tuple(row['C0_profile'])))[ic]
  assert ia is not None and ic is not None and ia//2==0 and ic//2==3
  support=[P[3],C[0],S[ia],S[ic]]
  gram=[[o,-(z if ma==2 else h),x[ia],x[ic]],
        [-(z if ma==2 else h),o,v[ia],v[ic]],
        [x[ia],v[ia],o,z],[x[ic],v[ic],z,o]]
  assert all(sum(rr,z)==z for rr in gram)
  vi=next(k for k,vert in enumerate(V)if set(support)<=set(vert))
  affine.append(dict(rank_survivor=row,support=support,actual_vertex=vi,
   gram=[[enc(q)for q in rr]for rr in gram],null_vector=[1]*4))
 assert len(affine)==4
 paths=[GEOM,FIRST,f'scripts/{TAG}_cut_profiles.py']
 return dict(status='PASS_MONOTONE_CUT_FAMILY_EXCLUDED',
  sources={p:sha(p)for p in paths},own_contact_guards=guards,profile_pairs=81,residual_cases=6300,
  field='Q(sqrt2,sqrt5), primitive generator sqrt2+sqrt5; descending rational coefficient lists',
  anchor_zero_cases=anchorzeros,distinct_residual_determinants=len(distinct),records=rows,surviving_rank_cases=zeros,
  rank_excluded_cases=total-len(zeros),affine_vertex_excluded_cases=affine,
  scope='Star case0 with four fixed B2 blocks, all positive leaves owncontacts3/3, and each cut owncontacts>=3. First-state forced-right and anchor-angle bounds persist under all monotone increases. Actual three-node paths exclude cut ownlabels>=6; both orientations of actual eight-root supports forbid every foreign contact when any ownlabel>=4. Thus the81 C0/C3 profiles and6300 projection/contact cases are complete for this monotone family. C1/C2 may vary because they are not needed for the principal-rank obstruction. Ownright contacts, changed positive-leaf labels and changed internal B2 labels are outside the gate.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_cut_profiles.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['residual_cases'],len(d['surviving_rank_cases']),flush=True)
