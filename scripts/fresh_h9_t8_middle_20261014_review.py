#!/usr/bin/env python3
"""Bind the middle-label extension to independent quartic algebra/domains."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_middle_20261014_review'
PRIMARY='results/fresh_h9_t8_middle_20261014.json'
BASE='results/fresh_h9_t8_asymmetric_20261014.json'
REVIEW='results/fresh_h9_t8_asymmetric_20261014_review.json'
ALGEBRA='results/fresh_h9_t8_middle_20261014_review_algebra.json'
NAMES=('e','a','b','d','c','delta')
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def value(x):return tuple(map(F,x))+(F(0),)*(4-len(x))
def run():
 docs={p:json.loads((ROOT/p).read_text())for p in[PRIMARY,BASE,REVIEW,ALGEBRA]};bound=set()
 def bind(p):
  if p in bound:return
  bound.add(p)
  for q,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():
   assert sha(q)==h,q
   if q.endswith('.json'):bind(q)
 for p in docs:bind(p)
 data,base,review,algebra=[docs[p]for p in[PRIMARY,BASE,REVIEW,ALGEBRA]]
 assert data['status']=='PASS_T8_TWO_MIDDLE_AND_LAST_BLOCK_LABELS_GE_4_ALL_OWN_PATTERNS_EXCLUDED'
 assert review['status']=='PASS_INDEPENDENT_FULL_T8_OWN_PATTERN_FAMILY_REVIEW'
 assert algebra['status']=='PASS_INDEPENDENT_256_QUARTIC_MIDDLE_BLOCK_CASES'
 free={tuple(x)for x in base['free_mandatory_pairs']}|{(7,8)}
 assert set(map(tuple,data['free_mandatory_pairs']))==free
 old={(a,b):m for a,b,m in base['source_mandatory_labels']};new={(a,b):m for a,b,m in data['source_mandatory_labels']}
 assert set(old)==set(new)and new[7,8]==5 and all(new[k]==v for k,v in old.items()if k not in free)
 # For p>=4, inverse I2(p) entries are1/(1-c²) and c/(1-c²), both
 # increasing for c>=0. Contact magnitude vectors are nonnegative.
 # Thus every inherited lower Schur-cost budget persists. Every saved
 # nonpositive actual-vertex witness also persists: increasing this
 # label lowers an off-diagonal pair of a symmetric Z-matrix.
 assert len(base['anchor_budgets'])==4 and len(base['third_budgets'])==2
 assert len(base['cross_right_guards'])==24 and len(base['own_support_patterns'])==9
 # A nonright foreign contact at p>=6 has cost>=1 before the positive
 # own contact retained by one of the four actual-vertex budgets.
 assert F(1,4)/(1-F(3,4))==1 and F(1,4)+1>1
 full={(b['id'],tuple(r['foreign_contacts']),tuple(r['third_contact_labels']))for b in base['blocks']for r in b['cases']}
 changed={k for k in full if any(v in[7,8]for v in k[1])};inherited=full-changed
 got_inherited={(r['block_id'],tuple(r['foreign_contacts']),tuple(r['third_contact_labels']))for r in data['inherited_cases']}
 assert got_inherited==inherited and len(inherited)==152
 separate={(r['block'],tuple(r['foreign']),tuple(r['third_labels'])):r for r in algebra['cases']}
 assert set(separate)==changed and len(changed)==256
 seen=set();count=0
 for row in data['cases']:
  key=row['block_id'],tuple(row['foreign_contacts']),tuple(row['third_contact_labels']);assert key not in seen;seen.add(key)
  own=separate[key];variables=row['free_variables'];indices=[NAMES.index(x)for x in variables]
  assert len(set(indices))==len(indices)
  for v,interval in zip(variables,row['variable_intervals']):assert list(map(value,interval))==[value(own['lower_bounds'][v]),(F(1),F(0),F(0),F(0))]
  for label,saved in [('positive_residual_numerator',own['certificates'][0]),('negative_anchor_numerator',own['certificates'][1])]:
   assert all(saved['degree'][j]==0 for j in range(6)if j not in indices)
   cert=row[label];assert cert['degrees']==[saved['degree'][j]for j in indices]
   for target,source in [('power_terms','unit_power_terms'),('positive_coefficients','positive_Bernstein_coefficients')]:
    a={tuple(r['index']):value(r['value'])for r in cert[target]};b={tuple(r['index'][j]for j in indices):value(r['value'])for r in saved[source]}
    assert a==b and len(a)==len(cert[target])
    if target=='positive_coefficients':count+=len(a)
 assert seen==changed and count==4320
 assert data['positive_coefficient_count']==algebra['positive_coefficient_count']==count
 minimum=(F(7,50),F(0),F(-31,2000),F(0));assert value(data['minimum_positive_coefficient'])==value(algebra['minimum_positive_coefficient'])==minimum
 assert F(7,50)**2>5*F(31,2000)**2
 paths=[PRIMARY,BASE,REVIEW,ALGEBRA,f'scripts/{TAG}.py']
 return dict(status='PASS_INDEPENDENT_T8_UNRESTRICTED_MIDDLE_LABEL_EXTENSION',sources={p:sha(p)for p in paths},
  recursively_bound_manifests=sorted(bound),changed_cases=256,inherited_cases=152,positive_coefficients=count,
  minimum_positive_coefficient=list(map(str,minimum)),
  domain_argument='Increasing label7-8 from4 preserves every literal non-PD symmetric-Z-matrix witness and every nonnegative inverse-Gram Schur-cost lower bound. For p>=6 an anchor foreign nonright contact to this block costs>=1, before a retained own contact, and is impossible. p=4 uses the408-case base. p=5 has256 independently reconstructed quartic cases and152 unchanged cases. All remaining rank/monotonicity arguments are inherited from the independently checked base theorem.',
  scope='Only the source-bound mandatory family: first internal label>=7, middle7-8 and last12-13 labels>=4, and all four anchor-own finite labels>=2. Other mandatory labels fixed. No whole-incidence exclusion or new finite-volume example.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['changed_cases'],d['positive_coefficients'],flush=True)
