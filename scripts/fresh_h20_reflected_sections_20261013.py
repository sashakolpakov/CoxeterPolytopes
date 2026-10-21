#!/usr/bin/env python3
"""Exact obstructions to finite reflection stars of twelve fixed H20 sections.

Floating point chooses actual-ridge proposals. A finite UNSAT core then
retains only exactly certified ridges and rational nonspherical root sets.
A plain Boolean propagation tree makes the final obstruction independently
replayable without a SAT/SMT solver.
"""
from collections import Counter
from itertools import combinations
import json
import numpy as np
import sympy as sp
import z3
import fresh_h20_joint_cusp_irrational_20261013 as prior
O=prior.O;old=O.old;P=old.P;t=old.t;ROOT=O.ROOT
SOURCE=prior.OUT;OUT=ROOT/'data/enumerations/fresh_h20_reflected_sections_20261013.json'

def scouts(Q,W,G,e,C):
 Qf=np.asarray(Q,float);Wf=np.asarray(W,float);ef=np.array([float(p.eval(2**.5)) for p in e]);cf=np.array([float(p.eval(2**.5)) for p in C]);q=ef@Qf@ef;v=Wf@Qf@ef
 H=np.asarray(G,float)-np.outer(v,v)/q;cp=Wf@Qf@cf;bad=[]
 for a,b in combinations(range(210),2):
  aa,bb,ab=H[a,a],H[b,b],H[a,b];D=aa*bb-ab*ab
  if min(aa,bb,D)<1e-8:continue
  ratio=ab*ab/(aa*bb)
  if ab<=0 and(ratio>=.749999 or min(abs(ratio-z) for z in [0,.25,.5,(3+5**.5)/8])<1e-8):continue
  ca=bb*cp[a]-ab*cp[b];cb=aa*cp[b]-ab*cp[a];values=D*cp-ca*H[:,a]-cb*H[:,b];scale=max(1.,max(abs(values)))
  if max(values)>1e-8*scale:continue
  if list(np.flatnonzero(abs(values)<1e-7*scale))!=[a,b]:continue
  bad.append((a,b))
 return bad

def symbolic_literal(x,variables):return variables[abs(x)-1] if x>0 else z3.Not(variables[-x-1])
def sat_core(clauses):
 variables=[z3.Bool('x%d'%i) for i in range(210)];solver=z3.Solver()
 for i,c in enumerate(clauses):solver.assert_and_track(z3.Or([symbolic_literal(x,variables) for x in c]),z3.Bool('c%d'%i))
 if solver.check()!=z3.unsat:return None
 core=[int(str(x)[1:]) for x in solver.unsat_core()]
 # Shrink the solver core to make the archived propositional proof small.
 for i in list(core):
  trial=[j for j in core if j!=i];S=z3.Solver();S.add(*[z3.Or([symbolic_literal(x,variables) for x in clauses[j]]) for j in trial])
  if S.check()==z3.unsat:core=trial
 return sorted(core)

def dpll(clauses,assignment=None):
 assignment={} if assignment is None else dict(assignment);units=[]
 while True:
  pending=[];unit=None
  for index,c in enumerate(clauses):
   if any(abs(l) in assignment and assignment[abs(l)]==(l>0) for l in c):continue
   left=[l for l in c if abs(l) not in assignment]
   if not left:return {'units':units,'conflict_clause':index}
   pending.append(left)
   if len(left)==1:unit=(left[0],index);break
  if unit is None:break
  literal,index=unit;units.append([literal,index]);assignment[abs(literal)]=literal>0
 if not pending:raise AssertionError('SAT model in purported unsatisfiable core')
 occurs=Counter(abs(l) for c in pending for l in c);variable=max(occurs,key=lambda v:(occurs[v],-v))
 return {'units':units,'branch_variable':variable,'false':dpll(clauses,{**assignment,variable:False}),'true':dpll(clauses,{**assignment,variable:True})}

def nonspherical_vector(G,ids):
 M=sp.Matrix(G[np.ix_(ids,ids)])
 for v in M.nullspace():return list(map(str,v)),str((v.T*M*v)[0])
 # Find a rational vector of nonpositive square by exact congruence.
 B=M.copy();T=sp.eye(len(ids))
 for k in range(len(ids)):
  if B[k,k]<=0:return list(map(str,T[:,k])),str(B[k,k])
  for j in range(k+1,len(ids)):
   c=B[k,j]/B[k,k];T[:,j]-=c*T[:,k]
  B=T.T*M*T
 raise AssertionError('positive definite block passed nonspherical gate')

def augment_cover(Q,W,G,h,sd,F,ids,only_pair=None):
 from scipy.optimize import linprog
 e,q,v,gr=sd;Qf=np.asarray(Q,float);Wf=np.asarray(W,float);hf=np.asarray(h,float)
 ef=np.array([float(p.eval(2**.5)) for p in e]);vf=Wf@Qf@ef;H=np.asarray(G,float)-np.outer(vf,vf)/(ef@Qf@ef)
 covered=set(ids)
 for a,b in combinations(range(210),2):
  if a in covered or b in covered or (only_pair is not None and (a,b)!=only_pair):continue
  aa,bb,ab=H[a,a],H[b,b],H[a,b];D=aa*bb-ab*ab
  if min(aa,bb,D)<1e-8:continue
  ratio=ab*ab/(aa*bb)
  if ab<=0 and(ratio>=.749999 or min(abs(ratio-z) for z in [0,.25,.5,(3+5**.5)/8])<1e-8):continue
  others=[i for i in range(210) if i not in (a,b)]
  lp=linprog(np.r_[np.zeros(22),-1.],A_ub=np.c_[Wf[others]@Qf,np.ones(len(others))],b_ub=np.zeros(len(others)),
    A_eq=np.c_[np.array([hf@Qf,ef@Qf,Wf[a]@Qf,Wf[b]@Qf]),np.zeros(4)],b_eq=[-1.,0.,0.,0.],bounds=[(None,None)]*22+[(0,None)],method='highs')
  if not lp.success or lp.x[-1]<1e-8:continue
  for denominator in (10000,1000000,100000000):
   raw=[P(sp.Rational(str(float(x))).limit_denominator(denominator)) for x in lp.x[:-1]]
   xe=old.pdot(Q,raw,e,F);point=[F.rem(q*x-xe*y) for x,y in zip(raw,e)]
   proof=old.witness(Q,W,h,sd,point,F,a,b)
   if proof is not None and proof['zero_walls_one_based']==[a+1,b+1]:return (a,b,proof)
 return None

def build():
 Q,W,G,h=O.base.load();source=json.loads(SOURCE.read_text());results=[]
 for ci,case in enumerate(source['cases']):
  for ri,row in enumerate(case['candidates']):
   F=O.Field(P(t*t-2),*map(sp.Rational,case['field_interval']));e=[sp.Poly.from_list(list(map(sp.Rational,p)),t,domain=sp.QQ) for p in row['normal']]
   C=[sp.Poly.from_list(list(map(sp.Rational,p)),t,domain=sp.QQ) for p in row['strict_center']]
   sd=prior.sd_from_e(Q,W,G,e,F);cp=[old.pdot(Q,C,[P(int(x)) for x in w],F) for w in W]
   bad=scouts(Q,W,G,e,C);clauses=[];metadata=[]
   for a,b in bad:clauses.append([a+1,b+1]);metadata.append({'kind':'actual_bad_ridge','walls':[a+1,b+1]})
   for a,b in combinations(range(210),2):
    if G[a,b]**2>=G[a,a]*G[b,b]:clauses.append([-a-1,-b-1]);metadata.append({'kind':'nonelliptic_pair','walls':[a+1,b+1]})
   extra=[];extra_ridges={};rejected_scouts=[];core=None
   for step in range(100):
    core=sat_core(clauses)
    if core is None:
     variables=[z3.Bool('x%d'%i) for i in range(210)];solver=z3.Optimize()
     solver.add(*[z3.Or([symbolic_literal(x,variables) for x in c]) for c in clauses]);solver.minimize(z3.Sum([z3.If(x,1,0) for x in variables]));assert solver.check()==z3.sat
     model=solver.model();ids=[i for i,x in enumerate(variables) if z3.is_true(model.eval(x))]
     if sp.Matrix(G[np.ix_(ids,ids)]).is_positive_definite:
      new=augment_cover(Q,W,G,h,sd,F,ids)
      if new is None:raise AssertionError(('Spherical star survives exact additional LP scout',ci,ri,ids))
      a,b,proof=new;extra_ridges[(a,b)]=proof;clauses.append([a+1,b+1]);metadata.append({'kind':'actual_bad_ridge','walls':[a+1,b+1]});continue
     smaller=ids[:]
     for i in ids:
      trial=[j for j in smaller if j!=i]
      if trial and not sp.Matrix(G[np.ix_(trial,trial)]).is_positive_definite:smaller=trial
     vector,square=nonspherical_vector(G,smaller)
     meta={'kind':'nonspherical_subset','walls':[i+1 for i in smaller],'nonpositive_vector':vector,'square':square}
     clauses.append([-i-1 for i in smaller]);metadata.append(meta);extra.append(meta);continue
    selected=[];selected_clauses=[];failed=None
    for index in core:
     meta=metadata[index]
     if meta['kind']=='actual_bad_ridge':
      a,b=[x-1 for x in meta['walls']];proof=extra_ridges.get((a,b)) or old.witness(Q,W,h,sd,C,F,a,b,cp)
      if proof is None or proof['zero_walls_one_based']!=sorted(meta['walls']):
       improved=augment_cover(Q,W,G,h,sd,F,[],(a,b));proof=None if improved is None else improved[2]
      if proof is None:failed=index;break
      extra_ridges[(a,b)]=proof;meta={**meta,'certificate':proof}
     selected.append(meta);selected_clauses.append(clauses[index])
    if failed is None:break
    print('discard unverified numerical ridge',ci,ri,metadata[failed]['walls'],flush=True);rejected_scouts.append(metadata[failed]['walls']);del clauses[failed];del metadata[failed]
   else:raise AssertionError('bounded exact-core refinement exhausted')
   assert core is not None
   tree=dpll(selected_clauses)
   result={'source_case':ci,'source_candidate':ri,'normal':row['normal'],'field_interval':[str(F.lo),str(F.hi)],
      'clauses':selected_clauses,'clause_certificates':selected,'boolean_unsat_tree':tree,
      'scout_actual_ridge_proposals':len(bad),'extra_nonspherical_subsets_scouted':extra,'exact_cached_ridges':len(extra_ridges),'discarded_numerical_ridge_proposals':rejected_scouts,
      'status':'NO_FINITE_SPHERICAL_REFLECTION_STAR'}
   results.append(result);print('certified',ci,ri,'core',len(core),'ridges',sum(x['kind']=='actual_bad_ridge' for x in selected),'extras',len(extra),flush=True)
   OUT.write_text(json.dumps({'complete':False,'cases':results},indent=2)+'\n')
 OUT.write_text(json.dumps({'schema':'finite-reflection-star-section-obstruction-v1','complete':True,'source':str(SOURCE.relative_to(ROOT)),
   'cases':results,'scope':'All finite spherical subgroups generated by subsets of the 210 original wall reflections, for the twelve literal fixed section planes. Covers their full stars and any convex chamber subunion within one such star containing the original chamber. Does not cover arbitrary finite chamber unions from infinite reflection subgroups or changing the section plane.'},indent=2)+'\n')
if __name__=='__main__':build()
