#!/usr/bin/env python3
"""Two explicit unbounded cusp links after the changed-profile225/226 caps."""
import argparse,json
from fractions import Fraction as F
from fresh_h20_four_cusp_completion_20261013_next_exact import *
from fresh_h20_changed_assembly_20261013_cusps import kernel1
SOURCE=OUT
OUT='data/enumerations/fresh_h20_four_cusp_completion_20261013_new_cusps.json'
def run():
 dependencies=bind_dependencies(SOURCE);Q,h,hc,roots,f,u,r=construct();roots.append((f,r));Qf=tuple(tuple(e(f,x)for x in row)for row in Q)
 d0=tuple(e(f,a+b/E(3))for a,b in zip(roots[224][1],roots[48][1],strict=True))
 d1=tuple(f.add(x,f.mul(f.rational(F(1,2)),e(f,roots[17][1][i]+roots[25][1][i])))for i,x in enumerate(r))
 data=[]
 for ci,delta in enumerate((d0,d1)):
  dc=fmv(f,Qf,delta);assert f.sign(fdot(f,delta,dc))==0 and f.sign(fdot(f,delta,tuple(e(f,x)for x in hc)))<0
  inc=[];signs=[]
  for j,(g,v)in enumerate(roots):
   H,p=pairing(f,dc,g,v);sg=H.sign(p);assert sg<=0;signs.append(sg)
   if not sg:inc.append(j)
  assert all(roots[i][0]is B or i==225 for i in inc)
  normals=[tuple(e(f,x)for x in roots[i][1])if roots[i][0]is B else roots[i][1]for i in inc]
  GF=[[fdot(f,a,fmv(f,Qf,b))for b in normals]for a in normals];assert all(x[1]==Z and x[0].b==0 for row in GF for x in row)
  G=[[x[0].a for x in row]for row in GF];sig=inertia(G);assert sig==(len(inc)-1,0,1)
  todo=set(range(len(inc)));components=[];affine_rank=spherical_rank=0;recession=None
  while todo:
   seed=min(todo);todo.remove(seed);cc=[seed]
   for i in cc:
    for j in sorted(todo.copy()):
     if G[i][j]<0:cc.append(j);todo.remove(j)
   M=[[G[i][j]for j in cc]for i in cc];inert=inertia(M);row={'walls':[inc[i]for i in cc],'gram':[[str(x)for x in a]for a in M],'inertia':list(inert)}
   if inert==(len(cc)-1,0,1):
    weights=kernel1(M);assert all(w>0 for w in weights)
    total=tuple(sumf(f,[f.mul(f.rational(w),normals[i][k])for w,i in zip(weights,cc,strict=True)])for k in range(21))
    nz=next(k for k,x in enumerate(delta)if x[1]==Z and x[0]!=Z);mu=tuple(y/delta[nz][0]for y in total[nz]);assert f.sign(mu)>0 and total==tuple(f.mul(mu,x)for x in delta)
    row.update(kind='affine',positive_kernel=list(map(str,weights)),positive_delta_multiple=f.encode(mu));affine_rank+=len(cc)-1
   else:
    assert inert==(len(cc),0,0);row['kind']='spherical';spherical_rank+=len(cc)
    if recession is None:
     I=inverse(M);recession=tuple(f.neg(sumf(f,[f.mul(f.rational(I[k][0]),normals[i][j])for k,i in enumerate(cc)]))for j in range(21))
   components.append(row)
  assert recession is not None
  wc=fmv(f,Qf,recession);wq=fdot(f,recession,wc);wp=[fdot(f,v,wc)for v in normals]
  assert f.sign(wq)>0 and f.sign(fdot(f,delta,wc))==0 and all(f.sign(x)<=0 for x in wp)and any(f.sign(x)<0 for x in wp)
  assert (affine_rank,spherical_rank)==((1,18)if ci==0 else(2,15))
  data.append({'cusp':ci,'defining_walls':[48,224]if ci==0 else[17,25,225],'null_ray':list(map(f.encode,delta)),'null_ray_center_pairing':f.encode(fdot(f,delta,tuple(e(f,x)for x in hc))),'incident_walls':inc,'all226_pairing_signs':signs,'components':components,'positive_quotient_rank':sig[0],'normal_span_rank':sig[0]+1,'positive_affine_rank':affine_rank,'spherical_rank':spherical_rank,'recession_vector':list(map(f.encode,recession)),'recession_square':f.encode(wq),'incident_recession_pairings':list(map(f.encode,wp)),'status':'unbounded'})
 return {'status':'PASS_TWO_NEW_EXACT_UNBOUNDED_CUSP_LINKS_AT226','source':SOURCE,'source_sha256':sha(SOURCE),'bound_prior_dependencies':dependencies,'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':['7','0']},'cusps':data,'scope':'These two actual feasible null rays have explicit positive-norm horospherical recession directions and hence unbounded links. Their affine/spherical ranks are1+18 and2+15; quotientrank19 alone is insufficient for the first. This is a scoped diagnosis, not a complete226cusp census.'}
def sumf(f,xs):
 out=f.zero
 for x in xs:out=f.add(out,x)
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print([(c['positive_affine_rank'],c['spherical_rank'])for c in d['cusps']])
