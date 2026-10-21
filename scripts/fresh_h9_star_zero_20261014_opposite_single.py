#!/usr/bin/env python3
"""Continuous opposite-single own-contact family on the last star block.

All other own contacts are3. High internal labels are removed by existing
literal gates;4/5 blocks leave polynomial rectangles and exact vertex guards.
"""
import argparse,json,math
from fractions import Fraction as F
from hashlib import sha256
from itertools import product,combinations
from pathlib import Path
from sympy import QQ,sqrt
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_zero_20261014_opposite_single'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json';GATES='results/fresh_h9_star_zero_20261014_gates.json'
S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
K=QQ.algebraic_field(sqrt(2),sqrt(5));z=K.zero;o=K.one;h=K.convert(QQ(1,2))
assert K.mod.to_list()==[1,0,-14,0,9]
co={2:z,3:h,4:K.from_sympy(sqrt(2)/2),5:K.from_sympy((1+sqrt(5))/4)}
def enc(x):return[str(t)for t in x.to_list()]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def sign2(a,b):
 sg=lambda x:(x>0)-(x<0)
 if not b:return sg(a)
 if not a:return sg(b)
 if sg(a)==sg(b):return sg(a)
 return sg(a)*sg(a*a-2*b*b)
def sign(x):
 p=[F(str(t))for t in x.to_list()];p=[F(0)]*(4-len(p))+p
 a3,a2,a1,a0=p;a=a0+7*a2;b=a1+17*a3;c=a1+11*a3;d=2*a2
 su=sign2(a,b);sv=sign2(c,d)
 if not sv:return su
 if not su:return sv
 if su==sv:return su
 return su*sign2(a*a+2*b*b-5*c*c-10*d*d,2*a*b-10*c*d)
def add(*ps):
 r={}
 for p in ps:
  for k,v in p.items():r[k]=r.get(k,z)+v
 return{k:v for k,v in r.items()if v!=z}
def mul(p,q):
 r={}
 for(a,b),v in p.items():
  for(c,d),w in q.items():r[a+c,b+d]=r.get((a+c,b+d),z)+v*w
 return{k:v for k,v in r.items()if v!=z}
def neg(p):return{k:-v for k,v in p.items()}
def det(a,b,c,d,e,f):
 return add(mul(a,add(mul(d,f),neg(mul(e,e)))),neg(mul(b,add(mul(b,f),neg(mul(c,e))))),mul(c,add(mul(b,e),neg(mul(c,d)))))
def bernstein(p,lo):
 n=max(k[0]for k in p);m=max(k[1]for k in p);power={}
 for(a,b),v in p.items():
  for i in range(a+1):
   for j in range(b+1):power[i,j]=power.get((i,j),z)+v*math.comb(a,i)*math.comb(b,j)*lo**(a+b-i-j)*(o-lo)**(i+j)
 out={}
 for i in range(n+1):
  for j in range(m+1):out[i,j]=sum((v*K.convert(QQ(math.comb(i,a),math.comb(n,a)))*K.convert(QQ(math.comb(j,b),math.comb(m,b)))for(a,b),v in power.items()if a<=i and b<=j),z)
 return out

def run():
 paths=[GEOM,GATES];seen=set()
 def read(p):
  d=json.loads((ROOT/p).read_text())
  if p not in seen:
   seen.add(p)
   for q,hh in d.get('sources',{}).items():
    assert sha(q)==hh,q
    if q.endswith('.json'):read(q)
  return d
 geo,gates=[read(p)for p in paths]
 g=geo['representatives'][0]['geometry'];V=[r['incident']for r in g['vertices']];sets=list(map(set,V))
 missing={frozenset(t)for t in g['minimal_nonfaces']}
 perm=list(range(18));perm[12],perm[13]=13,12
 assert {tuple(sorted(perm[x]for x in v))for v in V}==set(map(tuple,V))
 assert {frozenset(perm[x]for x in m)for m in missing}==missing
 own={tuple(sorted((leaf,a))):3 for i in range(3)for leaf in(P[i],C[i])for a in S[2*i:2*i+2]}
 own[11,12]=own[13,17]=2
 def witness(support):return next(i for i,v in enumerate(sets)if set(support)<=v)
 high=[]
 for i in range(4):
  lab=own|{(11,13):4,(12,17):4}|{tuple(S[2*j:2*j+2]):(6 if i==j else 4)for j in range(4)}
  hits=[f['id']for f in gates['families']if f['origin']in('two_dense','strong_contact')and all(lab.get((a,b),0)>=m for a,b,m in f['bounds'])]
  assert hits;high.append(dict(high_block=S[2*i:2*i+2],reviewed_gate=hits[0]))
 assert all(frozenset(S[2*i:2*i+2]+[leaf])in missing for i in range(4)for leaf in(P[i],C[i]))
 # Every cross-block contact has the4-3-4 affine path as a relaxed principal
 # minor. No leaf contact is used in these guards.
 cross=[]
 for i,j in combinations(range(4),2):
  for a in S[2*i:2*i+2]:
   for b in S[2*j:2*j+2]:
    support=S[2*i:2*i+2]+S[2*j:2*j+2]
    assert 1-F(1,2)-F(1,4)-F(1,2)+F(1,4)==0
    cross.append(dict(pair=sorted((a,b)),support=support,actual_vertex=witness(support),relaxed_path=[4,3,4],determinant='0'))
 budgets=[]
 for leaf,keep,drop in [(11,13,12),(17,12,13),(14,4,3)]:
  support=sorted((set(S)-{drop})|{leaf});budgets.append(dict(leaf=leaf,own_retained=keep,support=support,actual_vertex=witness(support)))
 assert F(1,2)+2*F(1,4)==1
 assert F(1,4)+2*2*F(1,4)>1 and F(1,4)+2*F(1,2)>1
 mutual_support=[11,13,17];mutual_vi=witness(mutual_support)
 assert 1-F(1,2)-F(1,2)==0
 def vertex_failure(labels,required_pair=None):
  for vi,v in enumerate(V):
   if required_pair is not None and not set(required_pair)<=sets[vi]:continue
   B=[[o if a==b else-co[labels.get(tuple(sorted((a,b))),2)]for b in v]for a in v]
   for k in range(9):
    pivot=B[k][k]
    if sign(pivot)<=0:return dict(actual_vertex=vi,prefix_order=k+1,last_pivot=enc(pivot))
    for a in range(k+1,9):
     for b in range(k+1,9):B[a][b]-=B[a][k]*B[k][b]/pivot
  return None
 upper=[]
 for last in(4,5):
  low=5 if last==4 else 4
  labels=own|{(11,13):low,(12,17):low}|{tuple(S[2*j:2*j+2]):(last if j==3 else 4)for j in range(4)}
  for pair in((11,14),(14,17)):
   failure=vertex_failure(labels|{pair:4},pair);assert failure
   upper.append(dict(last_internal=last,minimum_own_label=low,pair=pair,test_label=4,**failure))
 rows=[];positive=0;rejected=0;coefficients=0;minimum=None
 for blocks in product((4,5),repeat=4):
  inv=[o/(o-co[m]**2)for m in blocks];D=inv[3];C0=co[blocks[3]]*D
  low=5 if blocks[3]==4 else 4;lo=co[low]
  assert all(sign(v-K.convert(2))>=0 for v in inv)
  for extra in(None,2,3,4,5,6,7):
   cost=inv[0]*(o+co[blocks[0]])/2
   if extra is not None:cost+=inv[extra//2]/4
   Q=o-cost;wa=wb=z
   if extra in(6,7):wa=(D if extra==7 else C0)/2;wb=(D if extra==6 else C0)/2
   for mutual,ma,mb in product((2,3),repeat=3):
    # alpha is11-root13, beta is17-root12, both positive cosines.
    aa={(0,0):o,(2,0):-D};bb={(0,0):o,(0,2):-D}
    ab={(0,0):-co[mutual],(1,1):-C0}
    ac={(0,0):-co[ma],(1,0):-wa};bc={(0,0):-co[mb],(0,1):-wb}
    poly=det(aa,ab,ac,bb,bc,{(0,0):Q});bern=bernstein(poly,lo)
    row=dict(block_labels=blocks,foreign_third=extra,mutual_label=mutual,third_labels=[ma,mb],
       lower_own_label=low,variable_order=['alpha','beta'],interval_lower=enc(lo),interval_upper=['1'],
       polynomial=[dict(exponents=k,coefficient=enc(v))for k,v in sorted(poly.items())])
    if all(sign(v)>0 for v in bern.values()):
     row['status']='STRICT_POSITIVE_BERNSTEIN';row['bernstein']=[dict(index=k,coefficient=enc(v))for k,v in sorted(bern.items())]
     positive+=1;coefficients+=len(bern)
     for v in bern.values():
      if minimum is None or sign(v-minimum)<0:minimum=v
    else:
     lab=own|{tuple(S[2*i:2*i+2]):blocks[i]for i in range(4)}|{(11,13):low,(12,17):low,(11,17):mutual,(11,14):ma,(14,17):mb}
     if extra is not None:lab[tuple(sorted((14,S[extra])))]=3
     failure=vertex_failure(lab);assert failure,(blocks,extra,mutual,ma,mb)
     row['status']='NONPOSITIVE_ACTUAL_VERTEX_AT_MINIMUM_OWN_LABELS';row['vertex_failure']=failure;rejected+=1
    rows.append(row)
  print('blocks',blocks,'positive',positive,'vertex_rejected',rejected,flush=True)
 assert len(rows)==896 and positive==832 and rejected==64
 paths += [f'scripts/{TAG}.py']
 return dict(status='PASS_ALL_LABEL_OPPOSITE_SINGLE_FAMILY_EXCLUDED',sources={p:sha(p)for p in paths},case=0,
  fixed_other_own_contacts=sorted([list(p)+[m]for p,m in own.items()]),nonright_own_contacts=[[11,13,4],[12,17,4]],
  internal_label_lower_bounds=[S[2*i:2*i+2]+[4]for i in range(4)],
  literal_twin_swap=perm,arbitrary_high_block_reduction=high,cross_right_guards=cross,projection_budgets=budgets,
  mutual_upper_guard=dict(support=mutual_support,actual_vertex=mutual_vi,test_label=4,determinant_at_minimum='0'),
  third_contact_upper_guards=upper,cases=rows,positive_polynomial_cases=positive,actual_vertex_rejections=rejected,
  positive_bernstein_coefficients=coefficients,minimum_bernstein_coefficient=enc(minimum),
  field='Q(sqrt2,sqrt5); ANP descending rational coefficients at sqrt2+sqrt5',
  scope='Case0 first3 blocks have all six leaves own3/3; last two leaves have opposite single own nonright contacts of arbitrary finite labels>=4. All four internal labels are arbitrary finite integers>=4; last internal3 is outside this certificate. Both last-root orientations follow by a literal incidence swap. Every unlisted finite label and dotted distance remains free. Not a whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['positive_polynomial_cases'],d['actual_vertex_rejections'],flush=True)
