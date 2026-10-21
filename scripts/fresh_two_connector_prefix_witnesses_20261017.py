#!/usr/bin/env python3
"""Numerical coherent A-X-Y-R prefixes; no exact geometry promotion."""
import argparse,json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
import fresh_two_connector_native53_deformation_20261017 as D
T,N,S,H,B,G,M,ROOT=D.T,D.N,D.S,D.H,D.B,D.G,D.M,D.ROOT
TAG='fresh_two_connector_prefix_witnesses_20261017'

def run():
 z,labels,ns,h,pro,pr=H.setup();old=[(i,ns[i])for i in labels if i not in S.DELETE+D.D.WHITE];actual=[i for i in z['actual_sides']if i not in S.DELETE]
 k0,k1=T.a_curve(ns,h);aa=2**-.5*(k0+k1);bb=2**-.5*(-k0+k1)
 t=max(np.roots([aa@M@ns[53],1,bb@M@ns[53]]));a=aa*t+bb/t;out=[]
 sq={v['wall']:B.mid(v['square'])for v in z['projections']}
 for mr,ma in [(6,5),(7,3),(8,3),(10,2),(15,2)]:
  rr,info=N.norm_roots([ns[8],ns[44],ns[63]],[0,-B.theta(mr),0],2)
  r=[r for sg,r in rr if r@M@h<0][0]
  xs=D.conics(a,B.theta(ma),ns,h,old,r);ys=D.conics(r,B.theta(ma),ns,h,old,a)
  for x in xs:
   for y in ys:
    for ix in x['intervals']:
     for iy in y['intervals']:
      for factor in [1.001,1.01,1.1]:
       s=ix[0]*factor;xx=x['part']+x['aa']*s+x['bb']/s
       co=[xx@M@y['aa'],xx@M@y['part']+1,xx@M@y['bb']]
       for zz in np.roots(co):
        if abs(zz.imag)>1e-8 or not iy[0]+1e-8<zz.real<iy[1]:continue
        zz=float(zz.real);yy=y['part']+y['aa']*zz+y['bb']/zz;new=[a,r,xx,yy]
        assert max(abs(v@M@v-2)for v in new)<2e-6
        assert max(v@M@h for v in new)<-1e-8
        assert abs(a@M@ns[53]+1)<1e-6
        assert max(abs(float(a@M@ns[j]))for j in [7,8])<1e-6
        assert max(abs(float(r@M@ns[j]))for j in [8,63])<1e-6
        assert abs(r@M@ns[44]+B.theta(mr))<1e-6
        assert max(abs(v@M@ns[8])for v in [xx,yy])<1e-6
        assert max(abs(a@M@xx+B.theta(ma)),abs(r@M@yy+B.theta(ma)),abs(xx@M@yy+1))<1e-6
        assert not any(T.gates(T.native_values(v,0,labels,ns,pro))for v in new)
        pairs=[dict(pair=[i,j],product=float(new[i]@M@new[j]),relation=B.relation(float(new[i]@M@new[j])))for i,j in combinations(range(4),2)]
        assert all(p['relation']is not None and p['relation']!='parallel_or_unresolved'for p in pairs)
        assert T.corner(0,2,ma,0,2)is not None and T.corner(0,2,3,0,2)is not None
        nd=dict(ns);nd.update(dict(zip(range(67,71),new)));selected=actual+list(range(67,71))
        bd=S.C.detailed_cone(nd,selected)
        ic=H.incidence_score([nd[i]for i in selected],selected,[sq[i]for i in actual]+[2]*4)
        out.append(dict(Rlabel=mr,AXYR=ma,factor=factor,A_t=float(t),X_s=float(s),Y_z=float(zz),normals=[v.tolist()for v in new],pair_gates=pairs,
          center_products=[float(v@M@h)for v in new],boundary=bd,PD_incidence=ic))
        print('NUMERICAL_PREFIX',mr,ma,factor,bd['counts'],'ends',bd['open_ends'],'warnings',bd['warnings'],flush=True)
 assert len(out)==9 and all(not r['boundary']['warnings']for r in out)
 parent=json.loads((ROOT/'results'/f'{S.TAG}.json').read_text())
 paths=[Path(__file__),Path(D.__file__),Path(T.__file__),Path(N.__file__),Path(H.__file__),ROOT/'results'/f'{S.TAG}.json',H.PARENT,H.BOUNDARY]
 return dict(status='NINE_COHERENT_NUMERICAL_TWO_CONNECTOR_PREFIXES',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},count=len(out),records=out,
  deleted_walls=S.DELETE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],
  scope='Numerical, not exact. Full native/contact/common-center gates and complete numerical S6 cone checked; independent numerical PD incidence also recorded. Each prefix has21finite/9exterior intrinsic rays and13openends, worse than the frozen10-end control. B,C remain unsolved. No ambient finite census, geometry promotion, or finite-volume claim.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],out['count'],flush=True)
