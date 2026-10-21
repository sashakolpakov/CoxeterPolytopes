#!/usr/bin/env python3
"""Independent saved arithmetic; not an old-matrix or incidence replay.

Checks the exact face-mass deletion and ratio, and the seven-profile saved
Gram/inverse, cap pair signs, complete saved type blocks and coface loads.
Full old signed feasibility is discharged by the producer's separate
full-matrix replay. Historical native/type completeness is inherited.
"""
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
W='fresh_raw23_surviving_coface_witness_20261018'
T='fresh_core1_mixed_coface_capacity_pilot_20261018'
def read(t):return json.loads((ROOT/f'results/{t}.json').read_text())
def bind(d):
 for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
def field(x):return tuple(map(F,x))
zero=(F(),F());one=(F(1),F())
def plus(a,b):return(a[0]+b[0],a[1]+b[1])
def times(a,b):return(a[0]*b[0]+5*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def total(xs):
 z=zero
 for x in xs:z=plus(z,x)
 return z
def sign(x):
 a,b=x
 if a==0:return(b>0)-(b<0)
 if b==0:return(a>0)-(a<0)
 if a*b>0:return(a>0)-(a<0)
 d=a*a-5*b*b
 return((a>0)-(a<0))*((d>0)-(d<0))
def inner(a,g,b):return total(times(times(x,g[i][j]),y)for i,x in enumerate(a)for j,y in enumerate(b))

w,t=read(W),read(T);bind(w);bind(t)
assert w['status']=='PASS_RAW23_ONLY_FULL_MATRIX_LOWER_WITNESS'
assert t['status']=='PASS_CORE1_SEVEN_PROFILE_MIXED_COFACE_CAPACITY_PILOT'
old=read('fresh_raw357_coface_full_model_20261018');bind(old)
prior=old['candidate'];c=w['candidate']
oa={f:F(a)for f,a in prior['face_multipliers']};a={f:F(x)for f,x in c['face_multipliers']}
assert len(w['deleted_Q_sources'])==2331 and len(w['kept_Q_sources'])==50
gone={r['face']:r for r in w['deleted_Q_sources']};kept={r['face']:r for r in w['kept_Q_sources']}
assert len(gone)==2331 and len(kept)==50 and not(set(gone)&set(a))
assert set(a)==set(oa)-set(gone) and all(a[f]==oa[f]for f in a)
assert all(r['kind']=='Q'and r.get('raw_family')!=23 and F(r['mass'])==oa[f]for f,r in gone.items())
assert all(r['kind']=='Q'and r['raw_family']==23 and r['required']==10 and F(r['mass'])==a[f]for f,r in kept.items())
assert c['vertex_multipliers']==prior['vertex_multipliers']
mu={v:F(x)for v,x in c['vertex_multipliers']};normal=sum(mu.values(),F())
deleted=sum((F(r['mass'])*r['required']for r in gone.values()),F())
num=F(prior['objective_numerator'])-deleted
assert str(deleted)==w['deleted_required_mass']and num==F(c['objective_numerator'])
assert normal==F(c['normalizer'])and num/normal==F(c['bound'])>12
assert F(w['retained_simplex_required_mass'])+F(w['retained_raw23_required_mass'])==num
assert F(w['retained_raw23_required_mass'])==10*sum((a[f]for f in kept),F())
assert w['prepared']==old['prepared']and w['full_matrix_column_checks']==442226
assert not w['deletion_violations']and not w['deleted_rows_with_negative_coefficients']
assert w['certificates'][0]==old['certificates'][0]and w['certificates'][1]==c
assert w['optimizer_invocations']==0 and not w['model_changed']
g=[[field(x)for x in r]for r in t['gram']];j=[[field(x)for x in r]for r in t['inverse']]
assert len(g)==len(j)==5 and all(len(r)==5 for r in g+j)
assert all(g[i][k]==g[k][i]and j[i][k]==j[k][i]for i in range(5)for k in range(5))
assert all(total(times(g[i][k],j[k][h])for k in range(5))==(one if i==h else zero)for i in range(5)for h in range(5))
assert len(t['profiles'])==7 and all(len(p['contacts'])==5 for p in t['profiles'])
profiles={p['raw']:p for p in t['profiles']};assert set(profiles)=={23,24,25,33,36,38,39}
x={r:list(map(field,p['contacts']))for r,p in profiles.items()};tau={r:field(p['tau'])for r,p in profiles.items()}
for r,p in profiles.items():
 norm=inner(x[r],j,x[r]);assert norm==field(p['projection_norm'])
 assert plus(norm,times((F(1,2),F()),tau[r]))==(F(2),F())
expected={(r,s)for r in profiles for s in profiles if r<s}
assert {tuple(p['raw'])for p in t['same_coface_pairs']}==expected and len(t['same_coface_pairs'])==21
for p in t['same_coface_pairs']:
 r,s=p['raw'];dot=inner(x[r],j,x[s]);gap=plus(times(tau[r],tau[s]),times((F(-4),F()),times(dot,dot)))
 assert dot==field(p['projection_product'])and gap==field(p['fourfold_positive_product_gap'])
 assert sign(dot)<0 and sign(gap)>0 and p['same_actual_A1_coface_excluded']
expected={(r,s)for r in profiles for s in profiles if r<=s}
assert {tuple(p['raw'])for p in t['different_coface_pairs']}==expected and len(t['different_coface_pairs'])==28
for p in t['different_coface_pairs']:
 r,s=p['raw'];dot=inner(x[r],j,x[s]);assert dot==field(p['product'])
 gap=sign(dot)<0 and sign(plus(dot,one))>0
 assert gap==p['excluded_by_gap']==(r==s==23)
allcols=t['existing_target_columns']+t['absent_target_rows']
assert len(allcols)==len({z['base_type']for z in allcols})==235
sources=[r for z in allcols for r in z['sources']]
assert len(sources)==len({r['context']for r in sources})==len({r['face']for r in sources})==1645
assert sorted(sources,key=lambda r:r['context'])==sorted(t['typed_sources'],key=lambda r:r['context'])
violations=[]
for z in allcols:
 assert sorted(r['raw']for r in z['sources'])==sorted(profiles)and all(r['multiplicity']==1 for r in z['sources'])
 assert z['capacity']==z['forest'].count('A1')>=1
 coeff=sorted([[r['face'],1]for r in z['sources']]);assert coeff==z['positive_face_coefficients']
 load=sum((oa.get(f,F())for f,k in coeff),F());newload=sum((a.get(f,F())for f,k in coeff),F())
 sf=z['simplex_face']
 if sf is None:
  assert newload==0 and load==F(z['current_source_mass']);continue
 assert z['integer_face_coefficients']==sorted(coeff+[[sf,-z['capacity']]])
 cost=z['capacity']*oa.get(sf,F());newcost=z['capacity']*a.get(sf,F())
 assert load==F(z['current_load'])and cost==F(z['current_capacity'])and load-cost==F(z['current_excess'])
 assert newload==F(z['raw23_witness_load'])<=newcost==F(z['raw23_witness_capacity'])
 if load>cost:violations.append(z)
assert violations==t['violating_existing_columns']and len(violations)==3
assert len(t['raw23_capacity1_base_checks'])==100
assert len({(r['face'],r['base_type'])for r in t['raw23_capacity1_base_checks']})==100
assert {r['face']for r in t['raw23_capacity1_base_checks']}==set(kept)
assert all(sum(r['face']==f for r in t['raw23_capacity1_base_checks'])==2 for f in kept)
for r in t['raw23_capacity1_base_checks']:
 assert a[r['face']]==F(r['mass'])<=a[r['simplex_face']]==F(r['simplex_mass'])
print('PASS_SAVED_RAW23_WITNESS_AND_SEVEN_PROFILE_COFACE_ARITHMETIC',str(num/normal),len(allcols),len(sources),len(violations))
