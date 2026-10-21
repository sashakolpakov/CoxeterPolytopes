#!/usr/bin/env python3
"""Irrational cusp cuts with an exact global right-angle target.

The quotient direction already solves all local Coxeter-link equations.
Translating by the cusp vector then solves a second, global angle condition.
"""
from collections import Counter
import json
import numpy as np
import sympy as sp
import fresh_h20_joint_cusp_directions_20261013 as D
O=D.J.old;old=O.old;P=old.P;t=old.t;ROOT=D.ROOT
OUT=ROOT/'data/enumerations/fresh_h20_joint_cusp_irrational_20261013.json'

def sd_from_e(Q,W,G,e,F):
 q=old.pdot(Q,e,e,F);v=[old.pdot(Q,e,[P(int(x)) for x in w],F) for w in W]
 gr=lambda a,b:F.rem(q*int(G[a,b])-v[a]*v[b])
 return e,q,v,gr

def component_coefficients(W,e,comp,F):
 ids=[i-1 for i in comp];B=sp.Matrix(W[ids].T);rational=[]
 for power in range(F.p.degree()):
  x=sp.Matrix([p.nth(power) for p in e]);sol=B.gauss_jordan_solve(x)[0];rational.append(sol)
 return [F.rem(sum((t**power*rational[power][i] for power in range(F.p.degree())),0)) for i in range(len(ids))]


def cusp_proof(Q,W,G,u,sd,F,dropped):
 e,q,v,gr=sd;ids=[i for i,p in enumerate(W@Q@u) if not p and i+1!=dropped]
 H=[[gr(i,j) for j in ids] for i in ids]
 assert all(F.sign(H[i][i])>0 for i in range(len(ids)))
 for i in range(len(ids)):
  for j in range(i+1,len(ids)):
   assert F.sign(H[i][j])<=0
   N=F.rem(H[i][j]**2);den=F.rem(H[i][i]*H[j][j])
   assert any(F.rem(N-z*den).is_zero for z in D.J.ALLOWED)
 adj=np.asarray([[int(not p.is_zero) for p in row] for row in H]);result=[]
 for ks in O.components(adj):
  marks=D.kernel_field([[H[k][l] for l in ks] for k in ks],F);assert marks is not None
  walls=[ids[k] for k in ks]
  null=[F.rem(sum((c*(q*int(W[i,d])-v[i]*e[d]) for c,i in zip(marks,walls)),P(0))) for d in range(22)]
  index=next(d for d,z in enumerate(u) if z);mu=F.rem(null[index]/int(u[index]));assert F.sign(mu)>0
  assert all(F.rem(null[d]-mu*int(u[d])).is_zero for d in range(22))
  result.append({'walls':[i+1 for i in walls],'marks':[D.algebra.coeffs(p) for p in marks],
                 'cusp_multiple':D.algebra.coeffs(mu),'rank':len(ks)-1})
 assert sum(c['rank'] for c in result)==19
 return result

def run(limit=6):
 Q,W,G,h=O.base.load();local=json.loads(D.OUT.read_text());records=[];counts=Counter()
 for r in local['survivors']:
  if len(r['root']['polynomial'])!=3:continue
  F=D.algebra.Field(P(t*t-2),sp.Rational(1414213562373095,10**15),sp.Rational(1414213562373096,10**15))
  sign=1 if sp.Rational(r['root']['interval'][0])>0 else -1
  e0=[P(sp.Rational(a)+sign*t*sp.Rational(b)) for a,b in zip(r['e0'],r['e1'])];u=np.array(r['u'],dtype=object)
  cs=component_coefficients(W,e0,r['component'],F);di=r['component'].index(r['dropped']);sg=F.sign(cs[di])
  if not sg or any(F.sign(c)*sg>=0 for k,c in enumerate(cs) if k!=di):continue
  sd=sd_from_e(Q,W,G,e0,F);e,q,v,gr=sd;C=[F.rem(q*int(x)) for x in h]
  up=W@Q@u;targets=[];seen=set()
  for a,b in old.candidates(Q,W,h,sd,C,F):
   if up[a]*up[b] or not (up[a] or up[b]):continue
   den=F.rem(v[a]*int(up[b])+v[b]*int(up[a]))
   if den.is_zero:continue
   tau=F.rem(gr(a,b)*sp.invert(den,F.p));key=tuple(D.algebra.coeffs(tau))
   if key in seen:continue
   seen.add(key);e1=[F.rem(x+tau*int(z)) for x,z in zip(e0,u)];cs1=component_coefficients(W,e1,r['component'],F);sg1=F.sign(cs1[di])
   if not sg1 or any(F.sign(c)*sg1>=0 for k,c in enumerate(cs1) if k!=di):continue
   targets.append((a,b,tau,e1,cs1))
   if len(targets)==limit:break
  case={'cusp_index':r['cusp_index'],'u':r['u'],'component':r['component'],'dropped':r['dropped'],
        'source_path':r['path'],'source_direction_sign':sign,'source_e0':r['e0'],'source_e1':r['e1'],
        'base_normal':[D.algebra.coeffs(p) for p in e0],'base_component_coefficients':[D.algebra.coeffs(p) for p in cs],
        'field_polynomial':['1','0','-2'],'field_interval':[str(F.lo),str(F.hi)],'candidates':[]}
  for a,b,tau,e1,cs1 in targets:
   s1=sd_from_e(Q,W,G,e1,F);ee,qq,vv,gg=s1;he=old.pdot(Q,ee,[P(int(x)) for x in h],F)
   center=[F.rem(qq*int(h[d])-he*ee[d]) for d in range(22)]
   cp=[old.pdot(Q,center,[P(int(x)) for x in w],F) for w in W]
   if any(F.sign(p)>=0 for p in cp):center=O.four.alternate_center(Q,W,h,s1,F)
   row={'target_pair':[a+1,b+1],'translation':D.algebra.coeffs(tau),'normal':[D.algebra.coeffs(p) for p in ee],
        'component_coefficients':[D.algebra.coeffs(p) for p in cs1]}
   if center is None:row['status']='UNRESOLVED_CENTER';case['candidates'].append(row);continue
   cp=[old.pdot(Q,center,[P(int(x)) for x in w],F) for w in W]
   assert F.sign(gg(a,b))==0
   target=O.ridge.actual_ridge(Q,W,h,s1,center,F,a,b,cp)
   bad=None
   for x,y in old.candidates(Q,W,h,s1,center,F):
    bad=old.witness(Q,W,h,s1,center,F,x,y,cp)
    if bad is not None:break
   extra=[]
   # Record a second actual right ridge when available, avoiding the target.
   for x in range(210):
    if extra:break
    for y in range(x+1,210):
     if (x,y)==(a,b) or not gg(x,y).is_zero:continue
     proof=O.ridge.actual_ridge(Q,W,h,s1,center,F,x,y,cp)
     if proof is not None:extra.append(proof);break
   row.update(status='NON_COXETER_ACTUAL_RIDGE' if bad else 'OPEN_COXETER_CANDIDATE',
              strict_center=[D.algebra.coeffs(p) for p in center],normal_square=D.algebra.coeffs(qq),
              repaired_target_ridge=target,additional_right_ridges=extra,actual_forbidden_ridge=bad,
              cusp_link=cusp_proof(Q,W,G,u,s1,F,r['dropped']))
   case['candidates'].append(row);print('candidate',r['path'],r['dropped'],a+1,b+1,'tau',tau.as_expr(),row['status'],flush=True)
  records.append(case)
 for case in records:
  for row in case['candidates']:counts[row['status']]+=1
 OUT.write_text(json.dumps({'schema':'irrational-full-cusp-joint-repair-v1','cusp_source':str(O.CUSP_SOURCE.relative_to(ROOT)),
  'local_source':str(D.OUT.relative_to(ROOT)),'cases':records,'counts':dict(counts),
  'scope':'Two literal Q(sqrt2) cusp quotient directions and up to six deterministic global right-angle targets each. Complete exact roots of these linear translation equations only; no all-translations exclusion.'},indent=2)+'\n')
 print('DONE',dict(counts),flush=True)
if __name__=='__main__':run()
