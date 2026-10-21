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
 budgets=[]
 for i in range(4):
  for leaf in [P[i],C[i]]:
   for j in range(2):
    support=sorted((set(S)-{S[2*i+1-j]})|{leaf})
    vi=next(k for k,v in enumerate(V)if set(support)<=set(v))
    budgets.append(dict(leaf=leaf,retained_own_root=S[2*i+j],
     support=support,actual_vertex=vi))
 assert len(budgets)==16 and 3-6*F(1,2)==0
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
 orientations=[];detcache={}
 for orientation in range(16):
  base={(S[2*i],S[2*i+1]):4 for i in range(4)}
  fixed=[C[i]if orientation&(1<<i)else P[i]for i in range(4)]
  for i,leaf in enumerate(fixed):
   for j in range(2):base[tuple(sorted((leaf,S[2*i+j])))]=3
  rows=[];unguarded=[];zeros=[];cases_count=0
  thirds=[(3,3)]if orientation&1 else profiles
  for third,free in product(thirds,profiles):
   ap,bp=(free,(3,3))if orientation&8 else((3,3),free)
   labels=base.copy()
   for leaf,group,profile in [(C[0],0,third),(P[3],3,ap),(C[3],3,bp)]:
    for j in range(2):labels[tuple(sorted((leaf,S[2*group+j])))]=profile[j]
   upper=[]
   for anchor in [P[3],C[3]]:
    pair=tuple(sorted((C[0],anchor)))
    try:upper.append(failure(pair,4))
    except AssertionError:unguarded.append(dict(C0_profile=third,P3_profile=ap,C3_profile=bp,pair=pair))
   if len(upper)!=2:continue
   cases=[]
   for ia,x in patterns(3,ap):
    for ib,y in patterns(3,bp):
     qx=o-pdot(x,x);qy=o-pdot(y,y);xy=-pdot(x,y)
     for ic,v in patterns(0,third):
      qv=o-pdot(v,v)
      for ma,mb in product([2,3],repeat=2):
       key=(tuple(ap),tuple(bp),tuple(third),ia,ib,ic,ma,mb)
       if key not in detcache:
        xv=-co[ma]-pdot(x,v);yv=-co[mb]-pdot(y,v)
        detcache[key]=det3(qx,xy,xv,qy,yv,qv)
       det=detcache[key];cases_count+=1
       rec=dict(foreign_contacts=[ia,ib,ic],anchor_labels=[ma,mb],residual_determinant=enc(det))
       cases.append(rec)
       if det==z:
        lab=labels.copy()
        for leaf,extra in zip([P[3],C[3],C[0]],[ia,ib,ic]):
         if extra is not None:lab[tuple(sorted((leaf,S[extra])))]=3
        for anchor,m in zip([P[3],C[3]],[ma,mb]):lab[tuple(sorted((C[0],anchor)))]=m
        witness=None
        for vi,vert in enumerate(V):
         pp=pivots(gram(vert,lab))
         if sign(pp[-1])<=0:
          witness=dict(actual_vertex=vi,leading_order=len(pp),last_pivot=enc(pp[-1]));break
        zeros.append(dict(C0_profile=third,P3_profile=ap,C3_profile=bp,**rec,vertex_exclusion=witness))
   rows.append(dict(C0_profile=third,P3_profile=ap,C3_profile=bp,upper_witnesses=upper,cases=cases))
  unresolved=[v for v in zeros if v['vertex_exclusion']is None]
  closed=not unguarded and not unresolved
  orientations.append(dict(orientation=orientation,fixed_33_leaves=fixed,status='EXACT_ORIENTATION_EXCLUDED'if closed else'UNRESOLVED_ORIENTATION',
   residual_cases=cases_count,records=rows,missing_label4_guards=unguarded,zero_rank_cases=zeros))
  print('orientation',orientation,'cases',cases_count,'unguarded',len(unguarded),'zeros',len(zeros),'unresolved',len(unresolved),flush=True)
 paths=[GEOM,FIRST,f'scripts/{TAG}_orientations.py']
 return dict(status='PASS_ALL_SIXTEEN_ORIENTATIONS_EXCLUDED'if all(r['status']=='EXACT_ORIENTATION_EXCLUDED'for r in orientations)else'BOUNDED_ORIENTATION_FRONTIER',
  sources={p:sha(p)for p in paths},base_own_contact_guards=guards,
  all_leaf_both_orientation_budget_guards=budgets,forced_cross_right_witnesses=forced,
  own_pair_right_guard='Each free leaf has a nonright own contact by its Lanner triangle. The partner contact is 3. Testing their mutual contact at 3 gives the all-ones quadratic at most 3-6*(1/2)=0 on the saved actual three-face. Thus the pair is right before the label>=6 path bound is applied.',
  field='Q(sqrt2,sqrt5), primitive generator sqrt2+sqrt5; descending rational coefficient lists',
  orientations=orientations,unique_determinants_computed=len(detcache),
  scope='Sixteen choices of one fixed own-3/3 leaf in each B2 block. Every other own profile must be among eleven as in the prior exact gate, but actual finite-vertex upper-contact witnesses are rebuilt for each orientation with unused free leaves relaxed to right. A missing label4 guard or a surviving rank case remains unresolved. No global incidence automorphism is assumed.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_orientations.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],sum(r['status']=='EXACT_ORIENTATION_EXCLUDED'for r in d['orientations']),flush=True)
