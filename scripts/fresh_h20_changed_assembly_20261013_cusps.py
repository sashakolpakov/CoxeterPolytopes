#!/usr/bin/env python3
"""Exact doubled-source cusp incidence and seventeen deficient links.

This is a separate diagnostic of the frozen179-wall chamber, not an exhaustive
enumeration of all its ideal vertices. All arithmetic is rational/radical.
"""
import argparse,json,math
from fractions import Fraction as F
from fresh_h20_changed_assembly_20261013_exact import (
 setup,read,sha,ROOT,BASE,OUT as SOURCE,B,E,Z,O,vec,dot,mv,vp,
 pairing,encoded,bind_dependencies)
from fresh_dimension_bridge_20261013_verify import rank

OUT='data/enumerations/fresh_h20_changed_assembly_20261013_cusps.json'

def kernel1(M):
 a=[list(map(F,row))for row in M];n=len(a);piv=[];rr=0
 for col in range(n):
  k=next((i for i in range(rr,n)if a[i][col]),None)
  if k is None:continue
  a[rr],a[k]=a[k],a[rr];v=a[rr][col];a[rr]=[x/v for x in a[rr]]
  for i in range(n):
   if i!=rr:
    v=a[i][col];a[i]=[x-v*y for x,y in zip(a[i],a[rr],strict=True)]
  piv.append(col);rr+=1
 free=[i for i in range(n)if i not in piv];assert len(free)==1
 v=[F(0)]*n;v[free[0]]=F(1)
 for i,j in enumerate(piv):v[j]=-a[i][free[0]]
 if all(x<0 for x in v):v=[-x for x in v]
 assert all(x>0 for x in v)
 assert all(sum(x*y for x,y in zip(row,v,strict=True))==0 for row in M)
 return v

def run():
 dependencies=bind_dependencies(SOURCE);d=read(SOURCE)
 Q,_,_,old,_=setup();h=vec(d['strict_center']);hc=mv(Q,h)
 roots=[(B,vec(r))for r in d['new_core_normals']]+[old[i]for i in d['retained_old_indices']]+[(B,vec(c['normal']))for c in d['new_polars']]
 assert len(roots)==179
 a=vec(read(BASE)['old_normals'][24]);rays={}
 for ci,row in enumerate(read(BASE)['cusp_rows']):
  v=vec(row['ray'])
  for branch,w in enumerate((v,vp(v,a,-dot(v,mv(Q,a))/E(3)))):
   assert all(x.b==0 for x in w)
   den=math.lcm(*(x.a.denominator for x in w));z=[int(x.a*den)for x in w];gg=math.gcd(*z);z=tuple(x//gg for x in z)
   rays.setdefault(z,[]).append([ci,branch])
 assert len(rays)==204
 rows=[];defects=[];u=(Z,)*20+(O,);assert dot(u,mv(Q,u))==E(3)
 for z,provenance in rays.items():
  v=tuple(E(x)for x in z);vc=mv(Q,v);assert dot(v,vc)==Z and dot(v,hc).sign()<0
  inc=[];signs=[]
  for i,(f,r)in enumerate(roots):
   H,p=pairing(B,vc,f,r);sg=H.sign(p);assert sg<=0;signs.append(sg)
   if not sg:inc.append(i)
  row={'ray':list(z),'provenance':provenance,'incident':inc,'all_179_pairing_signs':signs}
  if all(i<27 for i in inc):
   rs=[roots[i][1]for i in inc];assert rank([[x.a for x in r]for r in rs])==19
   assert all(x.b==0 for r in rs for x in r)
   assert all(dot(r,mv(Q,u))==Z for r in rs)and dot(v,mv(Q,u))==Z
   gram=[[dot(r,mv(Q,s)).a for s in rs]for r in rs];todo=set(range(len(inc)));components=[];totalrank=0
   while todo:
    seed=min(todo);todo.remove(seed);comp=[seed]
    for i in comp:
     for j in sorted(todo.copy()):
      if gram[i][j]<0:comp.append(j);todo.remove(j)
    G=[[gram[i][j]for j in comp]for i in comp];weights=kernel1(G);rr=rank(G);assert rr==len(comp)-1;totalrank+=rr
    null=tuple(sum((E(w)*rs[i][k]for i,w in zip(comp,weights,strict=True)),Z)for k in range(21));nz=next(k for k,x in enumerate(v)if x!=Z);mu=null[nz]/v[nz]
    assert mu.sign()>0 and null==tuple(mu*x for x in v)
    components.append({'incident_walls':[inc[i]for i in comp],'rank':rr,'positive_kernel':list(map(str,weights)),'positive_delta_multiple':mu.enc()})
   assert totalrank==18
   row['deficient_link']={'normal_span_rank':19,'affine_quotient_rank':18,'required_quotient_rank':19,'affine_components':components,'recession_vector':encoded(u),'recession_norm':['3','0']}
   defects.append(len(rows))
  rows.append(row)
 assert len(defects)==17
 return {'status':'PASS_204_INHERITED_CUSP_SIGNS_AND_17_DEFICIENT_LINKS','source':SOURCE,'source_sha256':sha(SOURCE),'bound_dependencies':dependencies,
         'wall_count':179,'inherited_cusp_count':204,'deficient_cusp_indices':defects,'rays':rows,
         'scope':'All204 original/reflected source cusp rays are feasible. Seventeen have only short height-zero incident normals, positive affine quotient rank18 and the explicit spacelike vertical recession direction. Their H20 horospherical links are unbounded. This is not a complete cusp census of179 walls; no boundedness claim is made for all remaining source links here.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print('deficient',d['deficient_cusp_indices'])
