#!/usr/bin/env python3
"""Complete last-own-profile gate for star0, last internal4/5."""
import argparse,json
from fractions import Fraction as F
from itertools import product,combinations
from pathlib import Path
from hashlib import sha256
from sympy import symbols,expand
from fresh_h9_star_zero_20261014_own_three import K,z,o,h,r2,co,sign,enc,det3
from independent_h9_14_record30_a1_labels import finite_code
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_lastblock_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json';GATES='results/fresh_h9_star_zero_20261014_gates.json'
OWN='results/fresh_h9_star_zero_20261014_own_three.json';OPP='results/fresh_h9_star_zero_20261014_opposite_single.json'
S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 paths=[GEOM,GATES,OWN,OPP,'results/fresh_h9_star_zero_20261014_own_three_review.json','results/fresh_h9_star_zero_20261014_opposite_single_review.json','scripts/fresh_h9_star_zero_20261014_own_three.py','scripts/independent_h9_14_record30_a1_labels.py'];seen=set()
 def read(p):
  d=json.loads((ROOT/p).read_text())
  if p not in seen:
   seen.add(p)
   for q,hh in d.get('sources',{}).items():
    assert sha(q)==hh,q
    if q.endswith('.json'):read(q)
  return d
 geo,gates,owncert,oppcert=[read(p)for p in paths[:4]]
 for p in paths[4:6]:read(p)
 assert owncert['status']=='PASS_ALL_OWN_THREE_ARBITRARY_BLOCK_LABELS_EXCLUDED'
 assert oppcert['status']=='PASS_ALL_LABEL_OPPOSITE_SINGLE_FAMILY_EXCLUDED'
 g=geo['representatives'][0]['geometry'];V=[v['incident']for v in g['vertices']];sets=list(map(set,V));missing={frozenset(m)for m in g['minimal_nonfaces']}
 def face(s):
  support=sorted(s);vi=next(i for i,v in enumerate(sets)if set(support)<=v)
  return dict(support=support,actual_vertex=vi)
 own={tuple(sorted((leaf,a))):3 for i in range(3)for leaf in(P[i],C[i])for a in S[2*i:2*i+2]}
 assert all(frozenset(S[2*i:2*i+2]+[leaf])in missing for i in range(4)for leaf in(P[i],C[i]))
 high=[]
 for i in range(3):
  lab=own|{tuple(S[2*j:2*j+2]):(6 if i==j else 4)for j in range(4)}
  hits=[f['id']for f in gates['families']if f['origin']=='two_dense'and all(lab.get((a,b),0)>=m for a,b,m in f['bounds'])]
  assert hits;high.append(dict(high_block=i,reviewed_gate=hits[0]))
 cross=[]
 for i,j in combinations(range(4),2):
  for a in S[2*i:2*i+2]:
   for b in S[2*j:2*j+2]:cross.append(dict(pair=sorted((a,b)),**face(S[2*i:2*i+2]+S[2*j:2*j+2]),relaxed_path=[4,3,4],determinant='0'))
 assert 1-F(1,2)-F(1,4)-F(1,2)+F(1,4)==0
 budgets=[]
 for leaf,group in[(11,3),(17,3),(14,0)]:
  for j in range(2):
   keep=S[2*group+j];drop=S[2*group+1-j]
   budgets.append(dict(leaf=leaf,retained_own_root=keep,**face((set(S)-{drop})|{leaf})))
 high_own=[]
 for anchor in(11,17):
  for root in(12,13):high_own.append(dict(anchor=anchor,high_root=root,**face([14,anchor,root])))
 # All high-own labels >=6, including dense profiles, are excluded by a
 # negative residual plane. This identity is polynomial after D(1-c²)=1.
 D,c,alpha,beta,Q=symbols('D c alpha beta Q');X=D*(alpha**2+beta**2+2*c*alpha*beta)
 lhs=(X-1)*(Q-1+D/4)-(D*c*alpha+D*beta)**2/4
 rhs=(X-1)*(Q-1)+D*(alpha**2-1)/4
 assert expand(lhs-rhs).subs(c**2,1-1/D).expand()==0
 assert F(3,4)+2*F(1,4)>1
 assert F(1,4)+2*2*F(1,4)>1 and F(1,4)+2*F(1,2)>1
 high_lower=r2/4-o/8;assert sign(high_lower)>0
 # The lower bound is linear increasing in D>=2 after X>=3D/4.
 assert sign(3*r2/8-o/16)>0
 high_algebra=dict(identity='(X-1)*(Q-1+D/4)-(D*c*alpha+D*beta)^2/4=(X-1)*(Q-1)+D*(alpha^2-1)/4',relation='D*(1-c^2)=1',bounds=['D>=2','alpha^2>=3/4','beta>=0','Q>=1+sqrt2/2'],strict_determinant_lower=enc(high_lower),scope='Q excludes the possible last-block contribution D/4; no restriction on beta or the high finite label.')
 def exact_failure(lab,required=None):
  # The Dynkin classifier only locates a witness. Every accepted rejection
  # is independently certified here by exact field LDL pivots.
  for vi,v in enumerate(V):
   if required is not None and not set(required)<=sets[vi]:continue
   code=tuple(lab.get(tuple(sorted(p)),2)for p in combinations(v,2))
   if finite_code(9,code):continue
   B=[[o if a==b else-co[lab.get(tuple(sorted((a,b))),2)]for b in v]for a in v]
   piv=[]
   for k in range(9):
    pivot=B[k][k];piv.append(enc(pivot))
    if sign(pivot)<=0:return dict(actual_vertex=vi,prefix_order=k+1,pivots=piv)
    for a in range(k+1,9):
     for b in range(k+1,9):B[a][b]-=B[a][k]*B[k][b]/pivot
   raise AssertionError(('classifier witness not verified',vi,code))
  return None
 profiles={n:[q for q in product(range(2,6),repeat=2)if F(1,n)+sum(F(1,m)for m in q)<1]for n in(4,5)}
 assert [len(profiles[n])for n in(4,5)]==[11,13]
 assert all(pr==(3,3) or max(pr)>=4 for pp in profiles.values()for pr in pp)
 def pattern(pr,group):
  ans=[]
  for extra in([None]+[j for j in range(8)if j//2!=group]if pr==(3,3)else[None]):
   v=[z]*8;v[2*group]=-co[pr[0]];v[2*group+1]=-co[pr[1]]
   if extra is not None:v[extra]=-h
   ans.append((extra,v))
  return ans
 tr=pattern((3,3),0);domain=[];active={4:[],5:[]};upper=[]
 for n in(4,5):
  for a,b in product(profiles[n],repeat=2):
   row=dict(last_internal=n,anchor_profiles=[a,b])
   common=[j for j in range(2)if a[j]>=3 and b[j]>=3]
   if not common:
    assert a.count(2)==b.count(2)==1 and a.index(2)!=b.index(2)
    assert max(a)>=4 and max(b)>=4
    row['status']='INHERITED_OPPOSITE_SINGLE';domain.append(row);continue
   j=common[0];row['mutual_right_guard']=dict(**face([11,17,12+j]),ones_quadratic_at_test3=enc(2*o-2*(co[a[j]]+co[b[j]])))
   assert sign(2*o-2*(co[a[j]]+co[b[j]]))<=0
   bad=next((j for j in range(2)if sign(o-co[a[j]]**2-co[b[j]]**2)<=0),None)
   if bad is not None:
    row['status']='NONPOSITIVE_OWN_TRIANGLE';row['triangle']=face([11,17,12+bad]);row['determinant']=enc(o-co[a[bad]]**2-co[b[bad]]**2);domain.append(row);continue
   if a==b==(3,3):row['status']='INHERITED_ALL_OWN_THREE';domain.append(row);continue
   row['status']='RANK_CENSUS';domain.append(row);active[n].append((a,b))
   lab=own|{tuple(S[2*i:2*i+2]):(n if i==3 else 4)for i in range(4)}|{tuple(sorted((leaf,root))):m for leaf,pr in zip((11,17),(a,b))for root,m in zip((12,13),pr)}
   for pair in((11,14),(14,17)):
    fail=exact_failure(lab|{pair:4},pair);assert fail,(n,a,b,pair)
    upper.append(dict(last_internal=n,anchor_profiles=[a,b],pair=pair,test_label=4,**fail))
 print('domains',len(domain),'active',list(map(len,active.values())),'upper',len(upper),flush=True)
 records=[];zeros=[];unresolved=[];total=0;anchor_signs={};minors={}
 for blocks in product((4,5),repeat=4):
  inv=[o/(o-co[m]**2)for m in blocks]
  def pdot(x,y):return sum((inv[i]*(x[2*i]*y[2*i]+x[2*i+1]*y[2*i+1]+co[blocks[i]]*(x[2*i]*y[2*i+1]+x[2*i+1]*y[2*i]))for i in range(4)),z)
  base=own|{tuple(S[2*i:2*i+2]):blocks[i]for i in range(4)}
  for a,b in active[blocks[3]]:
   ar=pattern(a,3);br=pattern(b,3);values=[]
   aa=[o-pdot(v,v)for extra,v in ar];bb=[o-pdot(v,v)for extra,v in br];tt=[o-pdot(v,v)for extra,v in tr]
   ab=[[-pdot(v,w)for ey,w in br]for ex,v in ar];at=[[-pdot(v,w)for ey,w in tr]for ex,v in ar];bt=[[-pdot(v,w)for ey,w in tr]for ex,v in br]
   for ia,ib,ic in product(range(len(ar)),range(len(br)),range(7)):
    ad=aa[ia]*bb[ib]-ab[ia][ib]**2;sg=sign(ad);anchor_signs[sg]=anchor_signs.get(sg,0)+1
    for ma,mb in product((2,3),repeat=2):
     determinant=det3(aa[ia],ab[ia][ib],at[ia][ic]-co[ma],bb[ib],bt[ib][ic]-co[mb],tt[ic]);total+=1
     value=enc(determinant);values.append(value)
     if determinant==z:
      lab=base|{tuple(sorted((leaf,root))):m for leaf,pr in zip((11,17),(a,b))for root,m in zip((12,13),pr)}|{(11,14):ma,(14,17):mb}
      for leaf,(extra,v)in zip((11,17,14),(ar[ia],br[ib],tr[ic])):
       if extra is not None:lab[tuple(sorted((leaf,S[extra])))]=3
      fail=exact_failure(lab)
      zr=dict(block_labels=blocks,anchor_profiles=[a,b],foreign_contacts=[ar[ia][0],br[ib][0],tr[ic][0]],third_labels=[ma,mb],actual_vertex_failure=fail)
      zeros.append(zr)
      if fail is None:unresolved.append(zr)
   records.append(dict(block_labels=blocks,anchor_profiles=[a,b],pattern_order=[[x for x,v in ar],[x for x,v in br],[x for x,v in tr]],contact_order=[[2,2],[2,3],[3,2],[3,3]],determinants=values))
  print('blocks',blocks,'rank',total,'zeros',len(zeros),'unresolved',len(unresolved),flush=True)
 assert total==sum(8*7*4*sum(len(pattern(a,3))*len(pattern(b,3))for a,b in active[n])for n in(4,5))
 assert len(domain)==290 and len(upper)==168 and total==77952 and len(zeros)==32
 paths += [f'scripts/{TAG}.py']
 return dict(status='PASS_COMPLETE_LAST_BLOCK_4_5_OWN_PROFILE_EXCLUSION'if not unresolved else'COMPLETE_LAST_BLOCK_4_5_RANK_FRONTIER',sources={p:sha(p)for p in paths},case=0,fixed_first_three_own=sorted([list(p)+[m]for p,m in own.items()]),internal_lower_bounds=[S[2*i:2*i+2]+[4]for i in range(3)],last_internal_values=[4,5],first_three_high_reductions=high,cross_right_guards=cross,retained_root_budgets=budgets,high_own_triangle_guards=high_own,high_own_negative_plane=high_algebra,ordered_own_profiles=profiles,own_pair_domains=domain,third_upper_guards=upper,records=records,residual_cases=total,anchor_determinant_sign_counts=anchor_signs,zero_rank_cases=zeros,unresolved=unresolved,field='Q(sqrt2,sqrt5); descending ANP at sqrt2+sqrt5',scope='Star0 first3 blocks have all six leaves own3/3 and arbitrary internal labels; last internal is4 or5, all four last own labels arbitrary finite. Unlisted finite labels and dotted distances free. This is not whole-incidence exclusion and does not address other last internal labels.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['residual_cases'],len(d['zero_rank_cases']),len(d['unresolved']),flush=True)
