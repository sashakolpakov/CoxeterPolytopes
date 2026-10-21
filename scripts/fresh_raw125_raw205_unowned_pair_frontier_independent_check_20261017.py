#!/usr/bin/env python3
"""Independent saved-data arithmetic review; no graph or physical census replay."""
import json
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
from hashlib import sha256

ROOT=Path(__file__).resolve().parents[1]
def read(t): return json.loads((ROOT/'results'/f'{t}_20261017.json').read_text())
def digest(x): return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
a=read('fresh_raw125_raw205_unowned_pair_frontier')
b=read('fresh_raw125_raw205_actual_corner_frontier')
o=read('fresh_unowned_pair_optimistic_frontier')
for path,h in a['sources'].items(): assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
assert a['face_multipliers']==b['face_multipliers']==o['face_multipliers']
assert a['baseline_full_repair']==b['repairs'][0]
alpha={i:F(x) for i,x in a['face_multipliers']}
mu={i:F(x) for i,x in a['baseline_full_repair']['vertex_multipliers']}
omu={i:F(x) for i,x in o['baseline_full_repair']['vertex_multipliers']}
prior={r['vertex']:r for r in o['forest_records']}
changed={r['vertex']:r for r in b['forest_records']}
rows={r['vertex']:r for r in a['forest_records']}
assert len(rows)==len(changed)==96 and set(rows)==set(changed)
assert len(a['inherited_forests'])==475
for z in a['inherited_forests']:
 v=z['vertex']; assert v not in rows and digest(prior[v])==z['source_record_sha256']
 assert mu.get(v,F())==omu.get(v,F());rows[v]=prior[v]
assert len(rows)==571 and set(rows)==set(prior)|set(changed)
newpairs=unknown=0
for v,brow in changed.items():
 row=rows[v];pp=brow['newly_owned_pair_domain'];n=len(brow['literal_stars']);k=brow['old_literal_stars']
 assert len(pp)==n*(n-1)//2-k*(k-1)//2
 assert len({(i,j) for i,j,c in pp})==len(pp)
 assert all(0<=i<j<n and j>=k for i,j,c in pp)
 assert digest(pp)==row['new_owned_domain_sha256']
 assert row['corner_orbits']==brow['corner_orbits']
 newpairs+=len(pp);unknown+=sum(c is None for i,j,c in pp)
assert (newpairs,unknown)==(508351,8974)

counts=Counter();scenario_count=cert_count=charge_inequalities=0

def check_scenario(r,q):
 global scenario_count,cert_count,charge_inequalities
 scenario_count+=1
 v=r['vertex'];d=list(map(F,r['demands']));owner={i:j for j,g in enumerate(r['corner_orbits']) for i in g['literal_stars']}
 cover=q['split_cover'];total=F();marg=[F() for _ in d]
 for ids,amount in cover['terms']:
  x=F(amount);assert x>0 and len(ids)==len(set(ids));total+=x
  for i in ids:marg[owner[i]]+=x
 assert total==F(cover['total_mass']) and all(x>=y for x,y in zip(marg,d))
 lo=F(q['optimal_mass_lower_bound']);hi=F(q['optimal_mass_upper_bound'])
 assert 0<=F(q['baseline_optimal_mass'])<=lo<=hi<=total
 assert F(q['baseline_mu'])==mu.get(v,F())
 assert F(q['required_mu_increase_lower_bound'])==max(F(),lo-mu.get(v,F()))
 assert F(q['required_mu_increase_upper_bound'])==max(F(),hi-mu.get(v,F()))
 if q['method']=='exact_charge_and_cover':
  cert_count+=1;c=q['certificate'];vv=q['maximal_count_vectors'];ww=list(map(F,c['weights']));mix=[(j,F(x)) for j,x in c['mixture']]
  assert len(ww)==len(d) and all(x>=0 for x in ww)
  assert len({j for j,x in mix})==len(mix)
  assert all(len(t)==len(d) and all(isinstance(x,int) and x>=0 for x in t) for t in vv)
  assert all(sum((F(t[i])*ww[i] for i in range(len(d))),F())<=1 for t in vv)
  charge_inequalities+=len(vv)
  assert all(0<=j<len(vv) and x>0 for j,x in mix)
  assert all(sum((vv[j][i]*x for j,x in mix),F())>=d[i] for i in range(len(d)))
  val=sum((x for j,x in mix),F());assert val==sum((x*y for x,y in zip(ww,d)),F())==F(c['optimal_mass'])==lo==hi
  assert len(vv)==len(q['literal_representatives'])
  for vec,ids in zip(vv,q['literal_representatives']):
   assert len(ids)==len(set(ids));cc=Counter(owner[i] for i in ids)
   assert vec==[cc.get(i,0) for i in range(len(d))]

alllo=dict(mu);allhi=dict(mu);ranked={}
obj=F(a['baseline_full_repair']['objective_numerator']);den=sum(mu.values(),F())
assert den==F(a['baseline_full_repair']['normalizer'])
for v,r in rows.items():
 demand=[sum((alpha.get(fi,F())*F(m) for fi,m in g['positive_face_coefficients']),F()) for g in r['corner_orbits']]
 assert list(map(str,demand))==r['demands'] and all(x>0 for x in demand)
 q=r['all_unowned_deletion'];check_scenario(r,q)
 assert sum(c['complete_domain_pairs'] for c in r['candidate_domains'])==r['complete_unowned_pairs']==q['deleted_pairs']
 alllo[v]=max(mu.get(v,F()),F(q['optimal_mass_lower_bound']));allhi[v]=max(mu.get(v,F()),F(q['optimal_mass_upper_bound']))
 counts.update(forests=1,unowned_pairs=r['complete_unowned_pairs'],candidate_domains=len(r['candidate_domains']),positive_grand_forests=int(F(q['required_mu_increase_lower_bound'])>0),zero_gain_forests=int(F(q['required_mu_increase_upper_bound'])==0),exact_grand_cones=int(q['method']=='exact_charge_and_cover'))
 for c in r['candidate_domains']:
  s=c.get('scenario',c)
  if 'scenario' in c:
   check_scenario(r,s);assert s['deleted_pairs']==c['complete_domain_pairs']
  else: assert F(q['required_mu_increase_upper_bound'])==0 and F(s['required_mu_increase_lower_bound'])==F(s['required_mu_increase_upper_bound'])==0
  l=F(s['required_mu_increase_lower_bound']);h=F(s['required_mu_increase_upper_bound'])
  if h:ranked[v,tuple(c['corner_orbits'])]=(c['complete_domain_pairs'],l,h)
assert dict(counts)==a['counts'] and len(ranked)==len(a['ranked_candidate_domains'])==167
assert alllo==allhi
opt=a['optimistic_all_unowned'];assert F(opt['required_normalizer_lower'])==sum(alllo.values(),F())==F(opt['required_normalizer_upper'])
assert F(opt['hypothetical_ratio_lower'])==obj/sum(allhi.values(),F())==F(opt['hypothetical_ratio_upper'])
for r in a['ranked_candidate_domains']:
 cost,l,h=ranked[r['vertex'],tuple(r['corner_orbits'])]
 assert (r['pairs'],F(r['increase_lower']),F(r['increase_upper']))==(cost,l,h)
 assert r['certified_exact']==(l==h) and l==h
 assert F(r['hypothetical_ratio_lower'])==obj/(den+h)
 assert F(r['hypothetical_ratio_upper'])==obj/(den+l)
need=obj/12-den;s=a['shortlist'];assert need==F(s['required_increase'])>0
chosen=s['chosen'];assert len({r['vertex'] for r in chosen})==len(chosen)==8
assert sum(r['pairs'] for r in chosen)==s['complete_pair_count']==118
assert sum((F(r['increase_lower']) for r in chosen),F())==F(s['certified_increase_lower'])
assert F(s['hypothetical_ratio_upper'])==obj/(den+F(s['certified_increase_lower']))<12
assert all(l<need for cost,l,h in ranked.values())
# Independent exact-cost recurrence; no dominance pruning and no producer imports.
# The saved attaining witness bounds the relevant integer cost by118.
groups=defaultdict(list)
for (v,pair),(cost,l,h) in ranked.items():
 if l>0:groups[v].append((cost,l))
dp={0:F()}
for v,choices in sorted(groups.items()):
 nd=dict(dp)
 for c,g in dp.items():
  for cost,gain in choices:
   if c+cost<=118:nd[c+cost]=max(nd.get(c+cost,F(-1)),g+gain)
 dp=nd
assert len(groups)==59 and min(c for c,g in dp.items() if g>=need)==118
assert a['actual_new_incompatibilities']==0 and a['actual_exported_columns']==[] and a['actual_witness_unchanged'] is True
print('PASS independent saved-data arithmetic; source bindings',len(a['sources']),'new owned/unknown',newpairs,unknown,'forests',len(rows),'inherited',len(a['inherited_forests']))
print('saved scenarios',scenario_count,'exact charge/cover certificates',cert_count,'charge inequalities',charge_inequalities,'DP forests/cost states',len(groups),len(dp),'minimum pairs',118)
print('all-unowned ratio',opt['hypothetical_ratio_upper'],'shortlist ratio',s['hypothetical_ratio_upper'])
