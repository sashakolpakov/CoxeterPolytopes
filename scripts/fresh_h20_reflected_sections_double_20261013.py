#!/usr/bin/env python3
"""Fresh cusp-preserving sections after doubling the ambient chamber at wall169."""
from itertools import combinations
from collections import Counter
import json
import numpy as np
import sympy as sp
import fresh_h20_joint_cusp_irrational_20261013 as prior
O=prior.O;old=O.old;P=old.P;t=old.t;ROOT=O.ROOT
OUT=ROOT/'data/enumerations/fresh_h20_reflected_sections_double_20261013.json'

def normals(Q,W,G,wall):
 wr=W[wall];q=G[wall,wall];lookup={};rows=[]
 for tile in (0,1):
  for j,w in enumerate(W):
   if j==wall:continue
   v=w if tile==0 else w-sp.Rational(2*G[wall,j],q)*wr
   assert all(x==int(x) for x in v);key=tuple(map(int,v))
   if key not in lookup:lookup[key]=len(rows);rows.append({'normal':list(key),'owners':[]})
   rows[lookup[key]]['owners'].append([tile,j+1])
 return np.array([r['normal'] for r in rows],dtype=object),rows

def candidates(Q,R,G,e,C):
 Qf=np.asarray(Q,float);Rf=np.asarray(R,float);ef=np.array([float(p.eval(2**.5)) for p in e]);cf=np.array([float(p.eval(2**.5)) for p in C]);q=ef@Qf@ef;v=Rf@Qf@ef
 H=np.asarray(G,float)-np.outer(v,v)/q;cp=Rf@Qf@cf
 for a,b in combinations(range(len(R)),2):
  aa,bb,ab=H[a,a],H[b,b],H[a,b];det=aa*bb-ab*ab
  if min(aa,bb,det)<1e-8:continue
  ratio=ab*ab/(aa*bb)
  if ab<=0 and(ratio>=.749999 or min(abs(ratio-z) for z in [0,.25,.5,(3+5**.5)/8])<1e-8):continue
  ca=bb*cp[a]-ab*cp[b];cb=aa*cp[b]-ab*cp[a];vals=det*cp-ca*H[:,a]-cb*H[:,b]
  if max(vals)<=1e-8*max(1.,max(abs(vals))):yield a,b

def ridge(Q,R,h,sd,C,F,a,b,forbidden=False):
 e,q,v,gr=sd;aa,bb,ab=map(F.rem,(gr(a,a),gr(b,b),gr(a,b)));D=F.rem(aa*bb-ab*ab)
 if min(F.sign(aa),F.sign(bb),F.sign(D))<=0:return None
 N=F.rem(ab*ab);den=F.rem(aa*bb)
 reason=None
 if forbidden:
  if F.sign(ab)>0:reason='obtuse'
  elif F.sign(3*den-4*N)>0 and F.sign(N)!=0 and all(F.sign(p)!=0 for p in (4*N-den,2*N-den,16*N*N-12*N*den+den*den)):reason='acute_cosine_squared_below_3/4_but_not_Coxeter'
  else:return None
 cp=[old.pdot(Q,C,[P(int(x)) for x in w],F) for w in R]
 A=F.rem(bb*cp[a]-ab*cp[b]);B=F.rem(aa*cp[b]-ab*cp[a]);vals=[F.rem(D*cp[l]-A*gr(l,a)-B*gr(l,b)) for l in range(len(R))]
 signs=[F.sign(x) for x in vals]
 if max(signs)>0:return None
 zero=[i for i,x in enumerate(signs) if x==0];cones=[]
 for l in zero:
  if l in (a,b):continue
  alpha=F.rem(bb*gr(l,a)-ab*gr(l,b));beta=F.rem(aa*gr(l,b)-ab*gr(l,a))
  if min(F.sign(alpha),F.sign(beta))<0 or not F.rem(D*gr(l,l)-alpha*gr(l,a)-beta*gr(l,b)).is_zero:return None
  cones.append({'wall_one_based':l+1,'alpha_numerator':O.coeffs(alpha),'beta_numerator':O.coeffs(beta),'denominator':O.coeffs(D)})
 point=[F.rem(D*C[d]-A*(q*int(R[a,d])-v[a]*e[d])-B*(q*int(R[b,d])-v[b]*e[d])) for d in range(22)]
 square=old.pdot(Q,point,point,F);hp=old.pdot(Q,point,[P(int(x)) for x in h],F)
 if F.sign(square)>=0 or F.sign(hp)>=0:return None
 out={'ridge_walls_one_based':[a+1,b+1],'projected_gram':[O.coeffs(x) for x in (aa,ab,bb)],'point_coordinate_polynomials':[O.coeffs(x) for x in point],
 'point_square':O.coeffs(square),'point_h_pairing':O.coeffs(hp),'zero_walls_one_based':[i+1 for i in zero],'redundant_zero_wall_cones':cones}
 if reason:out['angle_rejection']=reason
 return out

def run(limit=6):
 Q,W,G,h=O.base.load();source=json.loads(prior.OUT.read_text());case=source['cases'][0];wall=168;R,record=normals(Q,W,G,wall);H=R@Q@R.T
 F=O.Field(P(t*t-2),*map(sp.Rational,case['field_interval']));e0=[sp.Poly.from_list(list(map(sp.Rational,p)),t,domain=sp.QQ) for p in case['base_normal']];u=np.array(case['u'],dtype=object)
 assert W[wall]@Q@u<0
 sd=prior.sd_from_e(Q,R,H,e0,F);e,q,v,gr=sd;C=[F.rem(q*int(x)) for x in h];up=R@Q@u;seen=set();selected=[]
 old_translations={tuple(row['translation']) for row in case['candidates']}
 proposals=[]
 for a in range(len(R)):
  if any(owner[0]==0 for owner in record[a]['owners']) or not up[a]:continue
  tau=F.rem(-v[a]/int(up[a]));key=tuple(O.coeffs(tau))
  if key in seen or key in old_translations:continue
  seen.add(key);ee=[F.rem(p+tau*int(z)) for p,z in zip(e0,u)];co=prior.component_coefficients(W,ee,case['component'],F);di=case['component'].index(case['dropped']);sg=F.sign(co[di])
  if not sg or any(F.sign(c)*sg>=0 for k,c in enumerate(co) if k!=di):continue
  proposals.append((abs(float(tau.eval(2**.5))),a,tau,ee,co))
 for _,a,tau,ee,co in sorted(proposals,key=lambda x:x[:2])[:limit]:selected.append((a,None,tau,ee,co))
 rows=[]
 for a,b,tau,ee,co in selected:
  sd=prior.sd_from_e(Q,R,H,ee,F);_,q,v,gr=sd;he=old.pdot(Q,ee,[P(int(x)) for x in h],F);C=[F.rem(q*int(h[d])-he*ee[d]) for d in range(22)]
  if any(F.sign(old.pdot(Q,C,[P(int(x)) for x in w],F))>=0 for w in W):
   original_sd=prior.sd_from_e(Q,W,G,ee,F);C=O.four.alternate_center(Q,W,h,original_sd,F)
  assert C is not None and F.sign(old.pdot(Q,C,C,F))<0
  assert all(F.sign(old.pdot(Q,C,[P(int(x)) for x in w],F))<0 for w in R)
  target=None
  for possible in range(len(R)):
   if possible==a or H[a,possible]!=0:continue
   target=ridge(Q,R,h,sd,C,F,a,possible)
   if target is not None:b=possible;break
  bad=None
  for x,y in candidates(Q,R,H,ee,C):
   bad=ridge(Q,R,h,sd,C,F,x,y,True)
   if bad is not None:break
  original_sd=prior.sd_from_e(Q,W,G,ee,F);link=prior.cusp_proof(Q,W,G,u,original_sd,F,case['dropped'])
  repaired=None
  aa=next(i for i,r in enumerate(record) if [0,1] in r['owners']);bb=next(i for i,r in enumerate(record) if [1,1] in r['owners'])
  assert gr(aa,bb).is_zero;repaired=ridge(Q,R,h,sd,C,F,aa,bb)
  rows.append({'target_pair':None if b is None else [a+1,b+1],'orthogonal_new_boundary_index':a+1,'translation':O.coeffs(tau),'normal':[O.coeffs(p) for p in ee],
    'component_coefficients':[O.coeffs(p) for p in co],'strict_center':[O.coeffs(p) for p in C],
    'normal_square':O.coeffs(q),'section_not_fixed_by_reflection_pairing':O.coeffs(old.pdot(Q,ee,[P(int(x)) for x in W[wall]],F)),
    'cusp_link_in_original_wall_labels':link,'target_ridge':target,'doubled_old_ridge_right_angle':repaired,
    'actual_forbidden_ridge':bad,'status':'NON_COXETER_ACTUAL_RIDGE' if bad else 'OPEN_COXETER_CANDIDATE'})
  print('double candidate',a+1,None if b is None else b+1,'tau',tau.as_expr(),'target',target is not None,'doubledrepair',repaired is not None,'status',rows[-1]['status'],flush=True)
 out={'schema':'changed-ambient-double-cusp-sections-v1','source':str(prior.OUT.relative_to(ROOT)),'source_case':0,'doubling_wall_one_based':wall+1,
   'ambient_boundary_normals':record,'ambient_boundary_count':len(R),'field_polynomial':['1','0','-2'],'field_interval':[str(F.lo),str(F.hi)],
   'candidates':rows,'scope':'The literal new global right-angle equations after ambient doubling at wall169; one irrational quotient direction and its declared translations. Not all ambient doublings, translations, or finite chamber unions.'}
 OUT.write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(R),len(rows),flush=True)
def save_control():
 Q,W,G,h=O.base.load();data=json.loads(OUT.read_text());source=json.loads(prior.OUT.read_text());case=source['cases'][0];oldrow=case['candidates'][0]
 F=O.Field(P(t*t-2),*map(sp.Rational,data['field_interval']));R=np.asarray([r['normal'] for r in data['ambient_boundary_normals']],dtype=object);H=R@Q@R.T
 e=[sp.Poly.from_list(list(map(sp.Rational,p)),t,domain=sp.QQ) for p in oldrow['normal']];C=[sp.Poly.from_list(list(map(sp.Rational,p)),t,domain=sp.QQ) for p in oldrow['strict_center']];sd=prior.sd_from_e(Q,R,H,e,F)
 a=next(i for i,r in enumerate(data['ambient_boundary_normals']) if [0,1] in r['owners']);b=next(i for i,r in enumerate(data['ambient_boundary_normals']) if [1,1] in r['owners'])
 repaired=ridge(Q,R,h,sd,C,F,a,b);assert repaired is not None and sd[3](a,b).is_zero
 bad=None
 for i,j in candidates(Q,R,H,e,C):
  bad=ridge(Q,R,h,sd,C,F,i,j,True)
  if bad is not None:break
 assert bad is not None
 data['control']={'source_candidate':0,'normal':oldrow['normal'],'strict_center':oldrow['strict_center'],
  'translation':oldrow['translation'],'component_coefficients':oldrow['component_coefficients'],
  'normal_square':O.coeffs(sd[1]),'cusp_link_in_original_wall_labels':oldrow['cusp_link'],
  'doubled_old_ridge_right_angle':repaired,'actual_forbidden_ridge':bad,'status':'NON_COXETER_ACTUAL_RIDGE'}
 OUT.write_text(json.dumps(data,indent=2)+'\n');print('CONTROL original 1/169 defect repaired to actual right ridge',a+1,b+1,flush=True)
if __name__=='__main__':run();save_control()
