#!/usr/bin/env python3
"""Coherent repair of17 inherited H20 cusp links and an all-height limitation."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from fresh_h20_changed_assembly_20261013_exact import (
 setup,read,sha,ROOT,OUT as SOURCE,B,E,Z,O,vec,dot,mv,vp,pairing,
 norm,e,label,encoded,bind_dependencies,strict_point)
from fresh_h20_changed_assembly_20261013_cusps import OUT as CUSPS
from fresh_dimension_bridge_20261013_verify import rank
OUT='data/enumerations/fresh_h20_changed_cusp_repair_20261013.json'

def state():
 Q,_,_,old,_=setup();d=read(SOURCE);h=vec(d['strict_center']);hc=mv(Q,h)
 roots=[(B,vec(r))for r in d['new_core_normals']]+[old[i]for i in d['retained_old_indices']]+[(B,vec(c['normal']))for c in d['new_polars']]
 assert len(roots)==179
 return Q,h,hc,roots

def run():
 dependencies=bind_dependencies(SOURCE);cd=read(CUSPS);assert cd['source']==SOURCE and cd['source_sha256']==sha(SOURCE)
 Q,h,hc,roots=state();old=roots[:];target=vec(read(SOURCE)['strict_escape']['vector']);tc=mv(Q,target)
 strict_point(target,Q,hc,roots);assert target[-1].sign()<0
 additions=[];family=[];links=[]
 for ci in cd['deficient_cusp_indices']:
  row=cd['rays'][ci];delta=tuple(E(x)for x in row['ray']);dc=mv(Q,delta)
  assert dot(delta,dc)==Z and dot(delta,hc).sign()<0 and delta[-1]==Z
  xd=dot(delta,tc);assert xd.sign()<0
  for sign in(1,-1):
   # Any spacelike normal in this plane has positive-rescaling form tδ±z.
   # An original core inequality gives the exact necessary lower bound t>=A.
   choices=[]
   for j,(f,r)in enumerate(old[:71]):
    assert f is B
    a=-dot(r,dc)
    if a.sign()>0:choices.append((E(3*sign)*r[-1]/a,j,a))
   A,j,a=max(choices,key=lambda x:x[0].a)
   assert A.b==0 and A.sign()>0
   cut=E(3*sign)*target[-1]/(-xd);assert cut.b==0
   proof={'cusp':ci,'sign':sign,'acute_wall':j,'acute_slope':a.enc(),'acute_lower_bound':A.enc(),'escape_nonnegative_upper_bound':cut.enc()}
   if (A-cut).sign()>0:proof['exclusion']='acuteness forces t strictly above the escape-cut threshold'
   else:
    k=next(k for k,(f,r)in enumerate(old[:27])if (-dot(r,dc)).sign()>0 and (E(3)*norm(f,r,Q)/E(4)-dot(r,dc)*dot(r,dc)*cut*cut).sign()>0)
    f,r=old[k];aa=-dot(r,dc);qq=norm(f,r,Q)
    assert r[-1]==Z and aa.sign()>0 and cut.sign()>0
    assert (E(3)*qq/E(4)-aa*aa*cut*cut).sign()>0
    proof.update(exclusion='for every A<=t<=cut, the indicated short-wall cosine lies strictly between0and1/2',angle_wall=k,angle_slope=aa.enc(),angle_wall_norm=qq.enc(),upper_cosine_square=(aa*aa*cut*cut/(E(3)*qq)).enc())
   family.append(proof)
   r=tuple(E(5)*x for x in delta[:-1])+(E(sign),);rc=mv(Q,r);q=dot(r,rc)
   assert q==E(3)and dot(r,hc).sign()<0 and dot(r,tc).sign()<0
   pairs=[]
   for k,(f,v)in enumerate(roots):
    H,p=pairing(B,rc,f,v);kind=label(H,p,e(H,q),e(H,norm(f,v,Q))if f is B else norm(f,v,Q))
    pairs.append({'old':k,'pairing':H.encode(p),'kind':kind})
   cusp_signs=[]
   for c in cd['rays']:
    v=tuple(E(x)for x in c['ray']);p=dot(r,mv(Q,v));assert p.sign()==(0 if c is row else-1);cusp_signs.append(p.enc())
   additions.append({'wall':len(roots),'cusp':ci,'sign':sign,'normal':encoded(r),'norm':q.enc(),'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(p['kind']for p in pairs).items())),'source_cusp_pairings':cusp_signs})
   roots.append((B,r))
  plus=roots[-2][1];minus=roots[-1][1]
  assert vp(plus,minus)==tuple(E(10)*x for x in delta)and dot(plus,mv(Q,minus))==E(-3)
  normals=[old[i][1]for i in row['incident']]+[plus,minus]
  assert all(dot(r,mv(Q,plus))==dot(r,mv(Q,minus))==Z for r in normals[:-2])
  assert all(x.b==0 for r in normals for x in r)
  assert rank([[x.a for x in r]for r in normals])==20
  # Previous positive affine components + the new parallel pair positively
  # span the whole19-dimensional Euclidean quotient, so recession is zero.
  weights=[]
  for component in row['deficient_link']['affine_components']:
   weights.extend((i,F(w))for i,w in zip(component['incident_walls'],component['positive_kernel'],strict=True))
  assert len(weights)==len(row['incident'])and {i for i,w in weights}==set(row['incident'])
  assert all(w>0 for i,w in weights)
  total=tuple(sum((E(w)*old[i][1][k]for i,w in weights),Z)+plus[k]+minus[k]for k in range(21))
  nz=next(k for k,x in enumerate(delta)if x!=Z);mu=total[nz]/delta[nz]
  assert mu.sign()>0 and total==tuple(mu*x for x in delta)
  links.append({'cusp':ci,'old_incident':row['incident'],'new_pair':[len(roots)-2,len(roots)-1],'old_positive_weights':[[i,str(w)]for i,w in weights],'positive_delta_multiple':mu.enc(),'normal_span_rank':20,'quotient_rank':19,'status':'bounded'})
 assert len(roots)==213 and len(additions)==34 and len(links)==17
 escape=strict_point(target,Q,hc,roots)
 return {'status':'PASS_EXACT_213_WALL_REPAIRED_LINKS_WITH_UNIVERSAL_CUSP_FAMILY_ESCAPE','sources':{SOURCE:sha(SOURCE),CUSPS:sha(CUSPS)},'bound_prior_dependencies':dependencies,'wall_count':213,'formula':'C_delta,sign=5*delta+sign*z, z²=3','new_caps':additions,'repaired_links':links,'all_real_family_proof':family,'strict_escape':escape,
         'scope':'The34 coherent caps close all17 specified deficient inherited links and retain179 actual facets. The saved179 strict escape survives213. More strongly, every acute-compatible Coxeter normal in any of the17 planes span(delta,z), with arbitrary real coefficient and either orientation, is strictly negative at that same escape. Thus any finite additions from those planes cannot complete retained179 to finite volume. Other cusp directions, wall replacements and general normals remain unrestricted.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print('walls',d['wall_count'],'links',len(d['repaired_links']))
