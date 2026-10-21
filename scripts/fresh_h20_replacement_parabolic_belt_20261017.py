#!/usr/bin/env python3
"""Small explicit affine-triangle completion test at the forward173/174 end."""
import argparse,hashlib,itertools,json,math
from collections import Counter
import numpy as np
import fresh_h20_replacement_outer_patch_20261017 as S
O,L,ROOT=S.O,S.L,S.ROOT;TAG='fresh_h20_replacement_parabolic_belt_20261017';PAIR='results/fresh_h20_replacement_outer_patch_pairs_20261017.json'
def run():
 D=O.P.P.load(S.INPUT);E=O.P.P.load(f'results/{S.TAG}.json');F=O.P.P.load(PAIR);GG=D['geometry'];Y0=np.array(GG['residual_normals']);CT0=np.array(GG['T_products']);J=np.array(GG['T_inverse']);q0=np.array(GG['native_norms']);m=np.array(GG['residual_form']);T=D['support'];ret=D['stages'][-1]['retained_ids']+[173,174];first=E['candidates'][0];choices=[(i,z)for i,z in enumerate(F['candidates'])if z['first_candidate']==0 and not any(j in T for j,lab in z['contacts'])];patterns=sorted(set(itertools.permutations((3,3,3)))|set(itertools.permutations((2,4,4)))|set(itertools.permutations((2,3,6))));out=[];records=[];stats=Counter()
 for idx,second in choices:
  Y=np.r_[Y0,[first['residual4'],second['residual4']]];q=np.r_[q0,2,2];bond=O.R.classify(float(-Y[173]@m@Y[174]/2));edge=next(z for z in second['boundary']['open_edges']if z['pair']==[173,174]);ends=[z for z in edge['hyperbolic_endpoints']if z['kind']=='ideal'];assert len(ends)==1;nu=np.r_[1.,ends[0]['point']];assert abs(nu@m@nu)<1e-7 and max(abs(Y[[173,174]]@m@nu))<1e-7;U=Y[[173,174]];inv=np.linalg.inv(U@m@U.T);V=np.array([[1.]+v['point']for v in second['boundary']['poly_vertices']if v['kind']=='finite'])
  for old,a,b in patterns:
   if old!=bond:continue
   f=np.array([-2*math.cos(math.pi/a),-2*math.cos(math.pi/b)]);p=(inv@f)@U;assert abs(p@m@p-2)<1e-7
   for lab in list(range(2,13))+['parallel']:
    target=-2 if lab=='parallel'else -2*(0 if lab==2 else math.cos(math.pi/lab));slope=float(nu@m@Y[171]);assert slope<0;h=float((target-p@m@Y[171])/slope);r=p+h*nu;co=-(r@m@Y.T)/np.sqrt(2*q);labels=[O.R.classify(float(x))for x in co];bad=[j for j in ret if labels[j]in('obtuse','angle_gap')];oldcorn=float(max(r@m@V.T));gram=np.r_[Y[[173,171]],[r]]@m@np.r_[Y[[173,171]],[r]].T;ev=np.linalg.eigvalsh(gram);stats['recipes']+=1
    rec=dict(second_candidate=idx,old_pair_label=bond,new_pair_labels=[a,b],label_to171=lab,height=h,particular=p.tolist(),forward_null=nu.tolist(),residual4=r.tolist(),old_vertex_maximum=oldcorn,second_endpoint_eigenvalues=ev.tolist(),failed_gates=[dict(wall=j,kind=labels[j],cosine=float(co[j]))for j in bad])
    if ev[0]<-1e-7:rec['rejection']='second endpoint indefinite';stats['second_endpoint_indefinite']+=1
    elif oldcorn>=-1e-8:rec['rejection']='old finite vertex not strict';stats['old_vertex_violation']+=1
    elif bad:rec['rejection']='retained native gate';stats['native_gate_failure']+=1
    else:
     rec['all_retained_labels']=[labels[j]for j in ret];rec['boundary']=L.boundary(np.r_[Y,[r]],m,ret+[175],T);out.append(rec);stats['accepted']+=1
    records.append(rec)
 paths=[f'scripts/{TAG}.py',S.INPUT,f'results/{S.TAG}.json',PAIR,f'scripts/{L.TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_AFFINE_TRIANGLE_BELT_ENDPOINT_TEST',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},second_candidate_ids=[i for i,z in choices],affine_triangle_patterns=[list(z)for z in patterns],statistics=dict(stats),records=records,candidates=out,scope='Explicit right-T affine triangle profiles at the unique forward ideal173/174 endpoint of22 right-T pair variants. Height fixed by171 label2..12 or parallel. All retained173 native gates, old finite vertices and second endpointPSD tested numerically. No other affine-component profile or arbitrary height exhaustion; no proper cusp or exact construction assertion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],flush=True)
