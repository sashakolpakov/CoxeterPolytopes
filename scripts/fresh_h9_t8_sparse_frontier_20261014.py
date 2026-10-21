#!/usr/bin/env python3
"""Narrow literal dependency audit; no new parameter bank or search.

Two actual four-node affine stars replace the old full-vertex upper-label
witnesses. Every other domain guard uses only the positive eight-block,
the selected third wall and the anchors. All earlier files stay frozen.
"""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_sparse_frontier_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
BASE='results/fresh_h9_t8_asymmetric_20261014.json'
MIDDLE='results/fresh_h9_t8_middle_20261014.json'
STATE='results/fresh_h9_t8_middle_20261014_search_case_08.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):return(x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def run():
 docs=[json.loads((ROOT/p).read_text())for p in[GEOM,BASE,MIDDLE,STATE]]
 for obj in docs:
  for p,h in obj['sources'].items():assert sha(p)==h,p
 geo,base,middle,state=docs
 assert base['status']=='PASS_T8_ALL_OWN_CONTACT_PATTERNS_N_GE_4_EXCLUDED'
 assert middle['status']=='PASS_T8_TWO_MIDDLE_AND_LAST_BLOCK_LABELS_GE_4_ALL_OWN_PATTERNS_EXCLUDED'
 vertices=[r['incident']for r in geo['representatives'][8]['geometry']['vertices']]
 z=(F(0),F(0));o=(F(1),F(0));h=(F(1,2),F(0));r=(F(0),F(1,2))
 gram=[[o if i==j else z for j in range(4)]for i in range(4)]
 for j,co in[(0,r),(2,h),(3,h)]:gram[1][j]=gram[j][1]=(-co[0],-co[1])
 w=[r,o,h,h]
 assert all(sum((mul(gram[i][j],w[j])[0]for j in range(4)),F(0))==0 and sum((mul(gram[i][j],w[j])[1]for j in range(4)),F(0))==0 for i in range(4))
 rows=[]
 for anchor in[11,17]:
  support=[anchor,14,2,4];vi=next(i for i,v in enumerate(vertices)if set(support)<=set(v))
  rows.append(dict(anchor=anchor,support=support,actual_vertex=vi,containing_vertex=vertices[vi],relaxed_gram=gram,positive_null_vector=w))
 # All remaining literal supports use only the eleven selected roots.
 S=base['positive_support'];core=set(S)|{11,14,17}
 assert len(S)==8 and len(core)==11
 for guard in base['cross_right_guards']+base['anchor_budgets']+base['third_budgets']:
  assert set(guard['support'])<=core
  assert set(guard['support'])<=set(vertices[guard['actual_vertex']])
 for pattern in base['own_support_patterns']:
  mutual=pattern['anchor_pair_guard'];support=mutual.get('support',[11,17])
  assert set(support)<=core and set(support)<=set(vertices[mutual['actual_vertex']])
 assert [11,12,13]in geo['representatives'][8]['geometry']['minimal_nonfaces']
 assert [12,13,17]in geo['representatives'][8]['geometry']['minimal_nonfaces']
 fixed={(1,14):2,(2,14):3,(4,14):3,(5,14):3,(4,5):4}
 bounds={(1,2):7,(7,8):4,(12,13):4}
 labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
 assert all(labels[p]==m for p,m in fixed.items())
 assert all(labels[p]>=m for p,m in bounds.items())
 assert all(any(labels[tuple(sorted((a,b)))]>2 for b in[12,13])for a in[11,17])
 def enc(x):
  if isinstance(x,F):return str(x)
  if isinstance(x,dict):return{k:enc(v)for k,v in x.items()}
  if isinstance(x,(tuple,list)):return[enc(v)for v in x]
  return x
 paths=[GEOM,BASE,MIDDLE,STATE,f'scripts/{TAG}.py']
 return enc(dict(status='PASS_SPARSE_LITERAL_DEPENDENCY_AND_SAVED_STATE_EXCLUSION',sources={p:sha(p)for p in paths},
  core_roots=sorted(core),fixed_labels=[[a,b,m]for(a,b),m in sorted(fixed.items())],finite_label_lower_bounds=[[a,b,m]for(a,b),m in sorted(bounds.items())],
  own_contacts='All four arbitrary finite>=2; each anchor has at least one nonright by its literal missing triple.',
  replacement_third_contact_guards=rows,other_guard_supports_all_within_core=True,
  unused_mandatory_labels=[[a,b,m]for(a,b),m in labels.items()if(a,b)not in fixed and(a,b)not in bounds and(a,b)not in[(11,12),(11,13),(12,17),(13,17)]],
  saved_state_matches_sparse_hypotheses=True,
  scope='A dependency reduction of the two frozen exact polynomial theorems; no new coefficient bank, metric search or parameter cutoff. The displayed actual affine stars replace their only full-vertex guard dependencies on outside leaves. All other contacts may vary within compact Coxeter requirements. The saved next mandatory state satisfies this sparse theorem and is excluded; the whole T incidence remains open and no later state was searched.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],[(r['anchor'],r['actual_vertex'])for r in d['replacement_third_contact_guards']])
