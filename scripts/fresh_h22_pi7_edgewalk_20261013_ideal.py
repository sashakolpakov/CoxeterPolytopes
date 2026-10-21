#!/usr/bin/env python3
"""Two actual ideal endpoints and their unbounded Euclidean links."""
import argparse,json,math
from fractions import Fraction as F
from pathlib import Path
from sympy import Matrix
import fresh_h22_incidence_prism_20261013_verify as V
P=V.P;TAG='fresh_h22_pi7_edgewalk_20261013';OUT=P/'results'/f'{TAG}_ideal.json'
def inverse(A):return[[F(x)for x in row]for row in Matrix(A).inv().tolist()]
def run(verify=False):
 M,R,h=V.setup();cert=V.read(V.OUT);front=V.read('results/'+TAG+'.json');recipes=V.read(V.REC)['recipes'];rows=[]
 for obj in(cert,front):
  for p,z in obj['sources'].items():assert V.sha(p)==z
 for row in front['rows']:
  if row['status']!='IDEAL_ENDPOINT_AUDIT_PENDING':continue
  edge=row['edge'];caps=[i for i in edge if i>=812];assert len(caps)==1;cap=caps[0];base=[i for i in edge if i<812];assert len(base)==20
  finite=cert['new_caps'][cap-812]['finite_contacts'];assert all(finite[str(i)]==2 for i in base)
  j=row['near_zero_candidates'][0]['wall'];assert finite[str(j)]==2
  B=[[V.pair(M,R[i],R[k])for k in base]for i in base];assert V.inertia(B)==[20,0,0];co=V.mv(inverse(B),[V.pair(M,R[j],R[i])for i in base]);p=[F(R[j][k])-sum(c*R[i][k]for c,i in zip(co,base))for k in range(23)]
  assert any(p)and V.pair(M,p,p)==0 and all(V.pair(M,p,R[i])==0 for i in base)
  if V.pair(M,p,h)>0:p=[-x for x in p]
  assert V.pair(M,p,h)<0;den=math.lcm(*(x.denominator for x in p));nums=[int(x*den)for x in p];gg=math.gcd(*nums);nums=[x//gg for x in nums];p=list(map(F,nums))
  oldpairs=[V.pair(M,p,r)for r in R];assert max(oldpairs)<=0 and oldpairs[row['released']]<0;active=[i for i,x in enumerate(oldpairs)if x==0];checks=[];activecaps=[]
  for k,rec in enumerate(recipes):
   u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm']);qu=V.pair(M,u,u);qv=V.pair(M,v,v);rad=(q/qv,F(0),-q*q*qu/(4*qv));L=-q*V.pair(M,u,p)/2;B0=rec['radical_branch']*V.pair(M,v,p)
   if L==B0==0:value=(0,0);activecaps.append(812+k)
   else:value=V.add(V.scale(V.T,L),V.scale(V.sqrt(V.ki(rad)),B0));assert value[1]<0
   checks.append({'wall':812+k,'tau_coefficient':str(L),'radical_coefficient':str(B0),'pairing_interval':list(map(str,value))})
  assert activecaps==[cap]and len(active)==22
  assert all(finite[str(i)]==2 for i in active)
  assert Matrix([R[i]for i in active]).rank()==21
  gram=[[V.pair(M,R[i],R[j])for j in active]for i in active];assert V.inertia(gram)==[20,0,2]
  comps=[];unseen=set(active)
  while unseen:
   todo=[unseen.pop()];cc=[]
   while todo:
    a=todo.pop();cc.append(a);near={b for b in unseen if V.pair(M,R[a],R[b])!=0};unseen-=near;todo.extend(near)
   gm=[[V.pair(M,R[i],R[j])for j in cc]for i in cc];inertia=V.inertia(gm);assert inertia==[len(cc)-1,0,1];ker=Matrix(gm).nullspace();assert len(ker)==1;vv=[F(x)for x in ker[0]]
   if vv[0]<0:vv=[-x for x in vv]
   assert min(vv)>0;comps.append({'walls':sorted(cc),'rank':len(cc)-1,'positive_null_coefficients_order':cc,'positive_null_coefficients':list(map(str,vv))})
  q=F(recipes[cap-812]['norm']);assert q>0
  rows.append({'vertex':row['vertex'],'released':row['released'],'edge':edge,'status':'CERTIFIED_IDEAL_ENDPOINT_WITH_UNBOUNDED_LINK','primitive_integer_null_ray':nums,'center_pairing':str(V.pair(M,p,h)),'released_wall_pairing':str(oldpairs[row['released']]),'active_base_walls':active,'active_new_caps':activecaps,'old_active_span_rank':21,'old_active_Gram_inertia':[20,0,2],'affine_components':comps,'all_cap_pairings':checks,'link_recession_direction':{'negative_of_wall':cap,'norm':str(q),'pairing_with_cap':str(-q),'pairings_with_all_other_active_walls':0},'geometry':'The23 active walls span22 ambient dimensions, so their common future null ray is an actual ideal endpoint. The cap is orthogonal to all22 active parent walls. Its negative is a nonzero Euclidean-link recession direction, so this ideal endpoint has an unbounded link and is not a finite-volume cusp.'})
 assert len(rows)==2
 sources=[V.OUT,V.REC,'results/'+TAG+'.json','scripts/'+Path(__file__).name,'scripts/'+V.PRE+'_verify.py']
 out={'status':'PASS_TWO_EXACT_IDEAL_ENDPOINTS_UNBOUNDED_LINKS','sources':{p:V.sha(p)for p in sources},'rows':rows,'scope':'These two edge directions have genuine feasible ideal endpoints; they are not unclosed hyperbolic rays. Their complete active-wall sets have an explicit Euclidean-link recession vector, so neither has a bounded cusp link. All835 inequalities and every additional-cap equality/strict sign are checked.'}
 if verify:assert out==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],[(r['vertex'],r['released'])for r in rows],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
