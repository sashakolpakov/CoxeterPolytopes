#!/usr/bin/env python3
"""No-project-import exact saved-data accounting check; no graph census."""
import json,hashlib
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
root=Path(__file__).resolve().parents[1];p=json.loads((root/'results/fresh_restricted_290_unowned_pair_frontier_20261017.json').read_text())
c=Counter();mu=dict((i,F(x))for i,x in p['baseline_restricted_repair']['vertex_multipliers']);lo=dict(mu);hi=dict(mu);alpha=dict((i,F(x))for i,x in p['face_multipliers'])
def cert(vs,d,cc):
 w=list(map(F,cc['weights']));assert min(w)>=0
 for v in vs:assert sum((x*y for x,y in zip(v,w)),F())<=1;c['charge_inequalities']+=1
 cover=[F()for _ in d];mass=F();seen=set()
 for i,x in cc['mixture']:
  x=F(x);assert x>0 and i not in seen;seen.add(i);mass+=x
  for j,a in enumerate(vs[i]):cover[j]+=x*a
  c['mixture_terms']+=1
 assert all(a>=b for a,b in zip(cover,d));assert mass==sum((a*b for a,b in zip(w,d)),F())==F(cc['optimal_mass']);c['certificates']+=1
for path,h in p['sources'].items():assert hashlib.sha256((root/path).read_bytes()).hexdigest()==h;c['bindings']+=1
for r in p['forest_records']:
 d=list(map(F,r['demands']));groups=r['corner_orbits'];owner={i:j for j,g in enumerate(groups)for i in g['literal_stars']}
 assert d==[sum((alpha.get(i,F())*F(a)for i,a in g['positive_face_coefficients']),F())for g in groups]
 pc=r['projected_baseline_cover'];cert(pc['count_vectors'],d,pc)
 for cl,vec in zip(pc['cliques'],pc['count_vectors']):assert vec==[sum(owner[i]==j for i in cl)for j in range(len(d))]
 q=r['all_unowned_deletion'];scenarios=[q]+[z['scenario']for z in r['candidate_domains']if 'scenario'in z]
 for s in scenarios:
  cover=[F()for _ in d];mass=F()
  for cl,x in s['split_cover']['terms']:
   x=F(x);mass+=x
   for i in cl:cover[owner[i]]+=x
  assert all(a>=b for a,b in zip(cover,d));assert mass==F(s['split_cover']['total_mass'])
  a,b=map(F,(s['optimal_mass_lower_bound'],s['optimal_mass_upper_bound']));assert a<=b<=mass
  if 'certificate'in s:cert(s['maximal_count_vectors'],d,s['certificate']);assert a==b==F(s['certificate']['optimal_mass'])
  assert F(s['required_mu_increase_lower_bound'])==max(F(),a-mu.get(r['vertex'],F()))
  assert F(s['required_mu_increase_upper_bound'])==max(F(),b-mu.get(r['vertex'],F()));c['scenario_bounds']+=1
 lo[r['vertex']]=max(mu.get(r['vertex'],F()),F(q['optimal_mass_lower_bound']));hi[r['vertex']]=max(mu.get(r['vertex'],F()),F(q['optimal_mass_upper_bound']))
 c.update(r['counts']);c['forests']+=1
n=F(p['baseline_restricted_repair']['objective_numerator']);o=p['optimistic_all_unowned'];assert lo==hi
assert sum(lo.values(),F())==F(o['required_normalizer_lower'])==F(o['required_normalizer_upper'])
assert n/sum(lo.values(),F())==F(o['hypothetical_ratio_lower'])==F(o['hypothetical_ratio_upper'])>12
assert c['unowned_pairs']==570398 and c['forests']==571
print('PASS_NO_PRODUCER_SAVED_FRACTION_CHECK',dict(c),flush=True)
