#!/usr/bin/env python3
"""All-label dense-third/single-anchor core obstruction, with six polynomials."""
import argparse,json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_third34_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
BATCH='results/fresh_h9_star_core_transfers_20261014_search_batch_011.json'
CHECKPOINT='results/full_norm_vinberg_and_star_contacts_2026_10_14_audit.json'
S=[3,4,6,7,9,10,12,13];CORE=sorted(S+[11,14,17]);R=s.sqrt(2)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(e):
 e=s.expand(e);b=e.coeff(R);a=s.expand(e-b*R)
 assert a.is_Rational and b.is_Rational,(a,b)
 return [str(a),str(b)]
def sign(e):
 a,b=map(Fraction,enc(e))
 if not a:return (b>0)-(b<0)
 if not b:return (a>0)-(a<0)
 if a>0 and b>0:return 1
 if a<0 and b<0:return -1
 q=a*a-2*b*b;assert q
 return ((q>0)-(q<0))*((a>0)-(a<0))
def run():
 seen={}
 def read(p):
  if p in seen:return seen[p]
  d=json.loads((ROOT/p).read_text());seen[p]=d
  for q,h in d.get('sources',{}).items():
   assert sha(q)==h,q
   if q.endswith('.json'):read(q)
  return d
 geo,batch,checkpoint=[read(p)for p in (GEOM,BATCH,CHECKPOINT)]
 assert checkpoint['status']=='PASS_RECURSIVE_SOURCES_FROZEN_PREDECESSOR_AND_LOCAL_LINKS'
 assert batch['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'and batch['cumulative_additional_nodes']==20837
 assert len(batch['pending_frames'])==9
 g=geo['representatives'][0]['geometry'];V=[sorted(v['incident'])for v in g['vertices']];VS=list(map(set,V))
 M={frozenset(m)for m in g['minimal_nonfaces']}
 assert any(set(S)<=v for v in VS)
 assert all(any(set(p)<=v for v in VS)for p in combinations(CORE,2))
 def face(nodes):
  nodes=sorted(nodes);vi=next(i for i,v in enumerate(VS)if set(nodes)<=v)
  return dict(support=nodes,actual_vertex=vi)
 guards=[];cos={3:s.Rational(1,2),4:R/2}
 def nullguard(nodes,labels,w,kind,test,cycle=False):
  n=len(nodes);G=s.eye(n)
  for i,m in enumerate(labels):G[i,i+1]=G[i+1,i]=-cos[m]
  if cycle:G[0,n-1]=G[n-1,0]=-s.Rational(1,2)
  w=s.Matrix(w);assert all(sign(x)>0 for x in w)
  assert all(s.expand(x)==0 for x in G*w)
  guards.append(dict(kind=kind,ordered_support=nodes,labels=labels,cycle=cycle,
    test_pair=sorted(test),positive_kernel=[enc(x)for x in w],gram=[[enc(x)for x in row]for row in G.tolist()],**face(nodes)))
 for i,j in combinations(range(4),2):
  for a,b in product(S[2*i:2*i+2],S[2*j:2*j+2]):
   aa=next(x for x in S[2*i:2*i+2]if x!=a);bb=next(x for x in S[2*j:2*j+2]if x!=b)
   nullguard([aa,a,b,bb],[4,3,4],[1,R,R,1],'cross_block_right',[a,b])
 for leaf in (11,17):
  nullguard([4,14,leaf],[4,4],[1,R,1],'third_anchor_label_at_most_three',[14,leaf])
  for root in (3,4):
   nullguard([root,14,leaf],[3,3],[1,1,1],'owner_contact_forces_third_contact_right',[14,leaf],True)
 assert len(guards)==30
 budgets=[]
 for leaf,drop,keep,cost in [(11,12,13,s.Rational(1,4)),(17,12,13,s.Rational(1,4)),(17,13,12,s.Rational(1,4)),(14,3,4,s.Rational(1,2))]:
  row=dict(leaf=leaf,removed_own_root=drop,retained_own_root=keep,own_cost_lower=str(cost),**face((set(S)-{drop})|{leaf}))
  assert len(row['support'])==8;budgets.append(row)
 assert s.Rational(1,4)+2*s.Rational(1,2)>1
 assert s.Rational(1,4)+2*s.Rational(1,2)>1
 assert s.Rational(1,2)+2*s.Rational(1,4)==1
 A,B,c,x,y,z=s.symbols('alpha beta c x y z',real=True)
 H=s.Matrix([[1,-c],[-c,1]]);J=H.inv();D=1/(1-c*c)
 assert s.simplify(J-D*s.Matrix([[1,c],[c,1]]))==s.zeros(2)
 p=s.Matrix([A,B]);X=(p.T*J*p)[0]
 assignment_counts=Counter((a.count(0),a.count(1))for a in product((None,0,1),repeat=2))
 assert set(assignment_counts)=={(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)}
 cases=[];count=0
 for n0,n1 in sorted(assignment_counts):
  u=s.Matrix([s.Rational(n0,2),s.Rational(n1,2)])
  F=(u.T*J*u)[0];T=(p.T*J*u)[0];delta=s.Rational(2-n0-n1,2)
  lower=(X-1)*(s.Rational(3,2)+F)-(delta+T)**2
  P=s.cancel((1-c*c)*lower)
  assert s.denom(P)==1 or not(s.denom(P).free_symbols)
  shifted=s.Poly(s.expand(P.subs({A:s.Rational(1,2)+x,B:R/2+y,c:R/2+z})),x,y,z)
  terms=[dict(exponents=list(mon),coefficient=enc(a))for mon,a in shifted.terms()]
  assert len(terms)==11 and all(sign(a)>0 for mon,a in shifted.terms())
  constant=shifted.coeff_monomial(1);assert sign(constant)>0
  # The last projection is Y_last>7/2, for both shared and opposite anchors.
  L,FM,E,MU=s.symbols('Y_last F_middle delta_actual anchor_cosine',real=True)
  residual=s.Matrix([[1-X,-E-T],[-E-T,2-2*MU-L-F-FM]])
  exact=residual.det()
  difference=s.cancel(exact-lower)
  target=(X-1)*(L-s.Rational(7,2)+FM+2*MU)+(delta-E)*(delta+E+2*T)
  assert s.cancel(difference-target)==0
  cases.append(dict(owner_contact_counts=[n0,n1],anchor_assignment_multiplicity=assignment_counts[n0,n1],
    ambient_contact_sum_upper=str(delta),polynomial=str(s.expand(P)),shifted_terms=terms,
    shifted_constant=enc(constant),strict_residual_difference_identity='(X-1)*(Y_last-7/2+F_middle+2*anchor_cosine)+(delta_max-delta_actual)*(delta_max+delta_actual+2*T)'))
  count+=len(terms)
 assert count==66 and sum(r['anchor_assignment_multiplicity']for r in cases)==9
 # Even omitting the positive cross term, X>=2*(1/4+1/2)>1.
 assert 2*(s.Rational(1,4)+s.Rational(1,2))>1
 assert s.Rational(3,4)+s.Rational(1,4)==1
 assert s.simplify(1/(2*(1-s.sqrt(3)/2))-(2+s.sqrt(3)))==0
 assert s.Rational(3)>s.Rational(3,2)**2
 perms=[];families=[]
 for opposite,owner_reverse,last_reverse in product((False,True),repeat=3):
  perm=list(range(18))
  if owner_reverse:perm[3],perm[4]=4,3
  if last_reverse:perm[12],perm[13]=13,12
  assert {tuple(sorted(perm[i]for i in v))for v in V}==set(map(tuple,V))
  assert {frozenset(perm[i]for i in m)for m in M}==M
  fixed=[[11,12,2],[11,13,3],[12,17,3 if opposite else 2],[13,17,2 if opposite else 3]]
  bounds=[[3,4,4],[6,7,4],[9,10,4],[12,13,7],[3,14,3],[4,14,4]]
  image={k:sorted([sorted((perm[a],perm[b]))+[m]for a,b,m in vals])for k,vals in [('fixed',fixed),('bounds',bounds)]}
  ident='dense_third_single_'+('opposite'if opposite else'shared')+('_owner_swap'if owner_reverse else'')+('_last_swap'if last_reverse else'')
  families.append(dict(id=ident,origin='all_label_dense_third_negative_plane',**image));perms.append(perm)
 lab={tuple(r[:2]):r[2]for r in batch['solution']['pair_labels']}
 assert any(all(lab[a,b]==m for a,b,m in f['fixed'])and all(lab[a,b]>=m for a,b,m in f['bounds'])for f in families)
 faces={tuple(sorted(set(v)&set(CORE)))for v in V}
 maximal=sorted(f for f in faces if not any(set(f)<set(g)for g in faces))
 assert len(maximal)==9
 paths=[GEOM,BATCH,CHECKPOINT,f'scripts/{TAG}.py']
 return dict(status='PASS_ALL_LABEL_DENSE_THIRD_SINGLE_ANCHOR_CORE',case=0,sources={p:sha(p)for p in paths},rank_core=CORE,
  complete_maximal_core_faces=[list(f)for f in maximal],positive_kernel_guards=guards,projection_budgets=budgets,
  foreign_domain='Anchors have at most one foreign contact, exactlypi/3 if nonright. Third retains owncost>=1/2, so any foreign nonright contact would give cost>=1 and is impossible.',
  angle_domain='Each third-anchor contact has label2or3. If that anchor has an owner-block nonright contact, its third contact is right by an actual affine-A2 relaxation.',
  variables=['alpha','beta','c'],intervals=['alpha>=1/2','beta>=sqrt2/2','sqrt2/2<=c<1'],
  physical_formulas=dict(owner_inverse='[[1,c],[c,1]]/(1-c^2)',X='(alpha^2+beta^2+2*c*alpha*beta)/(1-c^2)',
    F='(n0^2+n1^2+2*c*n0*n1)/(4*(1-c^2))',T='[n0*(alpha+c*beta)+n1*(beta+c*alpha)]/[2*(1-c^2)]',
    delta_max='(2-n0-n1)/2',residual_norms=['1-X','2-2*anchor_cosine-Y_last-F-F_middle'],residual_pair='-delta_actual-T'),
  continuous_cases=cases,positive_coefficient_count=count,
  tail='Every finite last internal label>=7 has t=cos(pi/m)>sqrt3/2. The anchor-sum last projection is1/(1-t^2)>4 when shared, or1/[2*(1-t)]>2+sqrt3>7/2 when opposite. Thus Y_last>7/2 in both cases; anchor_cosine>=0. No finite cutoff or coefficient-field assumption.',
  negative_plane='X>1, Y_last>7/2, anchor_cosine>=0 and F,F_middle>=0 make both residual norms negative. Each of six cleared lower polynomials has a strictly positive constant and positive coefficients after nonnegative shifts. Their positive determinant forces two negative directions, impossible in Lorentz signature(9,1).',
  literal_orientations=perms,families=families,covered_saved_cursor=20837,
  scope='Star0 first3internal labels>=4,last>=7; each anchor11/17 has exactly one ownpi/3 contact and one right contact, shared or opposite; third14ownlabels>=3and>=4 in eitherorientation. Every noncore contact and all remaining finite contacts/dotted distances free. All finite integer labels are covered. Not a whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['continuous_cases']),d['positive_coefficient_count'])
