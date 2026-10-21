#!/usr/bin/env python3
"""Exact conditional affine mechanism; numerical realization and boundary cost."""
import argparse,collections,hashlib,json,math
import numpy as np
import sympy as s
import fresh_h20_replacement150_boundary_20261017 as A
O,L,ROOT=A.O,A.L,A.ROOT
TAG='fresh_h20_replacement150_parabolic_cusp_20261017'
TRI='results/fresh_h20_replacement150_triangles_20261017.json'
CONT='results/fresh_h20_replacement150_continuation_20261017.json'
COMPONENTS=[list(range(1,11)),[12,13,14,15],[61,87],[151],[174,176]]
def metrics(b,ray=None):
 counts=collections.Counter();proper=[]
 if ray is not None:
  proper=[i for i,z in enumerate(b['ideal_edge_endpoints'])if np.linalg.norm(np.array(z['point'])-ray[1:])<1e-6];assert len(proper)==1
 for e in b['poly_edges']:
  if not e['actual_hyperbolic_edge']:continue
  n=sum(x['kind']=='ideal'and x['ideal_index']not in proper for x in e['hyperbolic_endpoints']);counts[int(n)]+=1
 return dict(finite_vertices=b['poly_vertex_counts'].get('finite',0),proper_ideal_vertices=len(proper),closed_edges=counts[0],one_open_edges=counts[1],two_open_edges=counts[2],open_ends=counts[1]+2*counts[2],exterior_vertices=b['poly_vertex_counts'].get('exterior',0),exterior_skeleton_components=len(b['exterior_skeleton_components']),maximum_squared_radius=b['maximum_squared_radius'])
def run():
 D,F,Y,CT,J,q,m,T,ret=A.geometry();ret=[j for j in ret if j!=150];R=O.P.P.load(TRI)['candidates'][20];C=O.P.P.load(f'results/{A.TAG}.json')['candidates'][14]['recipe'];E=O.P.P.load(CONT)['candidates'][112]
 assert R['profile']['contacts']==[]and C['contacts']==[]and E['contacts']==[[152,3],[172,4]]and E['label']==2
 for z in[R,C,E]:Y=np.r_[Y,[z['residual4']]];CT=np.r_[CT,[z['T_pairings']]];q=np.r_[q,2]
 ret += [175,176,177];assert len(ret)==175;old=E['boundary'];e=next(z for z in old['open_edges']if z['pair']==[174,176]);nu=np.array([1]+next(z['point']for z in e['hyperbolic_endpoints']if z['kind']=='ideal'));W=T+[174,176];whole=Y@m@Y.T+CT@J@CT.T;cen=O.census(W,whole[np.ix_(W,W)]);GW=2*s.eye(19)
 for i,j,label in cen['old_edges']:
  assert label==3;GW[W.index(i),W.index(j)]=GW[W.index(j),W.index(i)]=-1
 assert sorted(map(sorted,cen['components']))==sorted(map(sorted,COMPONENTS));assert abs(nu@m@nu)<1e-10 and nu[0]>0
 # Exact abstract normalized-root calculation, conditional on the base and nu.
 coeff=s.eye(20)[:19,:];rows=[];heights=[];records=[];pairs=Y@m@nu;outside=[j for j in ret if j not in W]
 assert max(abs(pairs[j])for j in W)<1e-8 and max(pairs[j]for j in outside)<-0.3
 V=np.array([[1]+x['point']for x in old['poly_vertices']if x['kind']=='finite']);V/=np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None]
 H=V.sum(axis=0);H/=math.sqrt(-H@m@H);assert H[0]>0 and max(Y[j]@m@H for j in ret if j not in T)<-1e-8
 for comp in COMPONENTS:
  wi=[W.index(j)for j in comp];v=s.zeros(19,1)
  for j in wi:v[j]=-1
  assert (v.T*GW*v)[0]==2
  pp=-np.sum(Y[comp],axis=0);dt=-np.sum(CT[comp],axis=0);base=pp@m@Y.T+CT@J@dt;bounds=[((base[j]+math.sqrt(2*q[j]))/(-pairs[j]),j)for j in outside];bound,block=max(bounds);h=max(1,math.ceil(bound)+1);heights.append(h);r=pp+h*nu
  assert max(r@m@V.T)<-1e-6
  margins=[float(-math.sqrt(2*q[j])-base[j]-h*pairs[j])for j in outside];assert min(margins)>0
  row=list(v)+[s.Integer(h)];rows.append(row);records.append(dict(component=comp,height=h,all_dotted_height_bound=float(bound),blocking_wall=block,minimum_outside_dotted_margin=min(margins),base_native_pairs=base.tolist(),nu_native_pairs=pairs.tolist(),residual4=r.tolist(),T_pairings=dt.tolist(),T_projected_square_exact='0'if set(comp)<=set(T)else'2',T_projected_square_numerical=float(2-dt@J@dt),maximum_old_finite_vertex_product=float(max(r@m@V.T))))
 assert heights==[8,7,7,7,24]
 coeff=coeff.col_join(s.Matrix(rows));Q=GW.row_join(s.zeros(19,1)).col_join(s.zeros(1,20));AG=coeff*Q*coeff.T
 assert coeff.rank()==20 and AG.rank()==19 and all(AG[i,i]==2 for i in range(24))
 relations=[]
 for k,comp in enumerate(COMPONENTS):
  z=s.zeros(24,1)
  for j in comp:z[W.index(j)]=1
  z[19+k]=1;assert AG*z==s.zeros(24,1) and (z.T*coeff)==s.Matrix([[0]*19+[heights[k]]]);relations.append(list(map(str,z)))
 RR=np.array([z['residual4']for z in records]);TT=np.array([z['T_pairings']for z in records]);YY=np.r_[Y,RR]
 # Four roots project to positive multiples of the future null nu. Their
 # inequalities are strict inside H3 and can be omitted from its side ledger.
 for z in records[:4]:assert np.linalg.norm(np.array(z['residual4'])-z['height']*nu)<1e-6
 full=L.boundary(YY,m,ret+[182],T);fullmetric=metrics(full,nu)
 assert fullmetric['proper_ideal_vertices']==1
 # One affine-A2 height family. At c=0 the160 contact is right; other
 # unprescribed native products are already strictly dotted and decrease.
 last=records[-1];base=np.array(last['base_native_pairs']);slope=np.array(last['nu_native_pairs']);threshold=np.sqrt(2*q);lo=max(0,max(base[j]/-slope[j]for j in outside));hi=max((base[j]+threshold[j])/-slope[j]for j in outside);h0=base[160]/-slope[160]
 other_margin=min(-threshold[j]-base[j]-h0*slope[j]for j in outside if j!=160);assert other_margin>0.7
 events=[];trials=[]
 for anchor in outside:
  for lab in list(range(2,13))+['parallel']:
   cv=1 if lab=='parallel' else(0 if lab==2 else math.cos(math.pi/lab));h=(base[anchor]+threshold[anchor]*cv)/-slope[anchor]
   if h<lo-1e-8 or h>=hi-1e-8:continue
   co=-(base+h*slope)/threshold;ls=[O.R.classify(float(x))for x in co];bad=[j for j in ret if ls[j]in('obtuse','angle_gap')];trials.append(dict(anchor=anchor,label=lab,height=float(h),bad=bad))
   if bad:continue
   r=-Y[174]-Y[176]+h*nu;assert max(r@m@V.T)<-1e-6
   if any(abs(h-z['height'])<1e-7 for z in events):continue
   b=L.boundary(np.r_[Y,[r]],m,ret+[178],T);events.append(dict(anchor=anchor,label=lab,height=float(h),residual4=r.tolist(),retained_native_labels=[ls[j]for j in ret],boundary=b,metrics=metrics(b,nu)))
 assert len(events)==11 and all(z['anchor']==160 for z in events)
 paths=[f'scripts/{TAG}.py',TRI,CONT,f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',f'scripts/{O.TAG}.py',f'scripts/{L.TAG}.py']
 return dict(status='EXACT_CONDITIONAL_AFFINE_COMPLETION_WITH_NUMERICAL_BOUNDARY_COST',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},base_retained_ids=ret,base_continuation_index=112,edge=[174,176],support=W,endpoint_census=cen,future_null_ray=nu.tolist(),maximum_nonincident_ray_product=float(max(pairs[j]for j in outside)),components=COMPONENTS,completion_records=records,exact_conditional=dict(normalized_old_gram=[[str(x)for x in row]for row in GW.tolist()],coefficient_matrix=[[str(x)for x in row]for row in coeff.tolist()],completed_gram=[[str(x)for x in row]for row in AG.tolist()],positive_relations=relations,euclidean_gram_rank=19,normal_span_rank=20,relation_multipliers=heights),common_face_center=dict(point=H.tolist(),norm=float(H@m@H),maximum_old_pair=float(max(Y[j]@m@H for j in ret if j not in T)),maximum_new_pair=float(max(RR@m@H))),base_metrics=metrics(old),full_completion_boundary=full,full_completion_metrics=fullmetric,height_family=dict(parameter='c=cos(pi/m); h(c)=(p.N160+sqrt(2*q160)*c)/(-nu.N160)',right_height=float(h0),parallel_height=float(hi),height_slope=float(threshold[160]/-slope[160]),minimum_other_dotted_margin_at_right=float(other_margin),numerical_acute_lower=float(lo),bounded_trials=trials,candidates=events),scope='Exact normalized affine Gram/rank/positive-relation identities are conditional on an exact coherent realization of the numerical175-wall base and its selected future null ray. Physical old-wall margins, chosen heights, and H3 boundary ledgers are numerical. Four new roots have future-null T projections and impose no interior H3 side. The fifth caps the174/176 endpoint. No exact promotion, full ambient vertex census, or finite-volume conclusion.')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');z=a.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if z.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(p.read_text()))
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['base_metrics'],d['full_completion_metrics'],flush=True)
