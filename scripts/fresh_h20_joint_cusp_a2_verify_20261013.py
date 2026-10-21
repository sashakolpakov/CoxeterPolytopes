#!/usr/bin/env python3
"""Independent stdlib/Fraction proof of the full affine-A2 section obstruction."""
from copy import deepcopy
from fractions import Fraction as F
import json
import fresh_h20_bivariate_verify_20261013 as B
import fresh_h20_new_route_20261013_verify as K
A=B.A;ROOT=B.ROOT
SOURCE=ROOT/'data/enumerations/fresh_h20_joint_cusp_a2_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_joint_cusp_a2_verify_20261013.json'

def certify(row,Q,W,G,h):
 i,j=[x-1 for x in row['positive_walls_one_based']];k=row['negative_wall_one_based']-1;a,b=[x-1 for x in row['ridge_walls_one_based']]
 e,q,v,gr=B.section(Q,W,G,i,j,k)
 wantedA=B.summ([B.power(B.S,2),B.mul(B.S,B.T),B.scale(B.S,-1),B.power(B.T,2),B.T,B.ONE])
 assert q==B.scale(wantedA,4)==B.read(row['normal_square'])
 complete=B.summ([B.power(B.sub(B.S,B.const(F(1,2))),2),B.const(F(3,4)),B.mul(B.S,B.T),B.power(B.T,2),B.T])
 assert wantedA==complete # Strictly positive for s,t>0.
 zplus,zminus=[list(map(F,row[key]['point'])) for key in ('positive_face','negative_face')]
 for key,z,ids in [('positive_face',zplus,[a,b,i,j]),('negative_face',zminus,[a,b,k])]:
  assert sorted(row[key]['walls_one_based'])==sorted(x+1 for x in ids);B.check_face(z,ids,Q,W,h)
 up=B.pdot(e,B.cv(zplus),Q);um=B.pdot(e,B.cv(zminus),Q)
 assert [up,um]==[B.read(x) for x in row['face_normal_pairings']]
 B.positive(up);B.positive(B.scale(um,-1))
 p=[B.sub(B.scale(up,zminus[d]),B.scale(um,zplus[d])) for d in range(22)]
 assert p==[B.read(x) for x in row['ridge_point']]
 assert B.pdot(p,e,Q)==B.ZERO
 square=B.pdot(p,p,Q);assert square==B.read(row['ridge_square']);B.positive(B.scale(square,-1))
 B.positive(B.scale(B.pdot(p,B.cv(h),Q),-1))
 for l,w in enumerate(W):
  value=B.pdot(p,B.cv(w),Q)
  if l in (a,b):assert not value
  else:B.positive(B.scale(value,-1))
 # The normal e and the two wall normals are independent: e has strictly
 # positive pairing with zplus, whereas both wall normals vanish there.
 # All are in timelike p-perp, hence their Gram is positive definite.
 assert G[a][a]*G[b][b]>G[a][b]**2
 aa,ab,bb=gr(a,a),gr(a,b),gr(b,b)
 assert [aa,ab,bb]==[B.read(x) for x in row['projected_gram']]
 assert ab==B.scale(q,G[a][b]) and G[a][b]<0
 N=B.mul(ab,ab);den=B.mul(aa,bb)
 assert B.sub(B.scale(N,4),den)==B.scale(B.mul(q,B.power(B.T,2)),8)
 gap=B.summ([B.power(B.sub(B.S,B.const(F(1,2))),2),B.const(F(3,4)),B.mul(B.S,B.T),B.T])
 assert B.sub(den,B.scale(N,2))==B.scale(B.mul(q,gap),8)
 # Construct a strict interior independently of the ridge.
 cp,cm=[list(map(F,r['point'])) for r in row['strict_center_faces']]
 B.check_face(cp,[i,j],Q,W,h);B.check_face(cm,[k],Q,W,h)
 xp=B.pdot(e,B.cv(cp),Q);xm=B.pdot(e,B.cv(cm),Q);B.positive(xp);B.positive(B.scale(xm,-1))
 center=[B.sub(B.scale(xp,cm[d]),B.scale(xm,cp[d])) for d in range(22)]
 assert center==[B.read(x) for x in row['strict_center']];assert B.pdot(center,e,Q)==B.ZERO
 B.positive(B.scale(B.pdot(center,center,Q),-1));B.positive(B.scale(B.pdot(center,B.cv(h),Q),-1))
 for w in W:B.positive(B.scale(B.pdot(center,B.cv(w),Q),-1))
 # Omitted wall k is the positive combination (ri+s*rj)/t in the section.
 assert all(B.sub(B.add(B.const(W[i][d]),B.scale(B.S,W[j][d])),B.scale(B.T,W[k][d]))==e[d] for d in range(22))
 # Active A1 kernel: (1+t) Pi(ri)+(s+t) Pi(rj)=2t*u.
 ri=[B.sub(B.scale(q,W[i][d]),B.mul(v[i],e[d])) for d in range(22)]
 rj=[B.sub(B.scale(q,W[j][d]),B.mul(v[j],e[d])) for d in range(22)]
 U=[W[i][d]+W[j][d]+W[k][d] for d in range(22)]
 for d in range(22):
  left=B.add(B.mul(B.add(B.ONE,B.T),ri[d]),B.mul(B.add(B.S,B.T),rj[d]))
  assert left==B.scale(B.mul(B.T,q),U[d])
 return True

def run():
 raw=json.loads(SOURCE.read_text());Q,W,G,h=A.source();u=list(map(F,raw['cusp_vector']));support=[x-1 for x in raw['support_one_based']]
 assert support==[59,98,124] and A.dot(u,u,Q)==0 and A.dot(u,h,Q)<0
 assert all(A.dot(w,u,Q)<=0 for w in W)
 assert all(sum(W[i][d] for i in support)==2*u[d] for d in range(22))
 assert [[G[i][j] for j in support] for i in support]==[[4,-2,-2],[-2,4,-2],[-2,-2,4]]
 incident=[i for i,w in enumerate(W) if A.dot(w,u,Q)==0];components=K.groups(G,incident)
 assert [sorted(x-1 for x in c['walls_one_based']) for c in raw['cusp_components']]==components
 for row,ids in zip(raw['cusp_components'],components):
  M=[[G[i][j] for j in ids] for i in ids];marks=K.kernel(M);assert marks==row['marks'] and len(ids)-1==row['rank']
  assert all(G[i][j]<=0 and F(G[i][j]**2,G[i][i]*G[j][j]) in (0,F(1,4),F(1,2),F(3,4),1) for ni,i in enumerate(ids) for j in ids[ni+1:])
  null=[sum(m*W[i][d] for m,i in zip(marks,ids)) for d in range(22)];index=next(d for d,x in enumerate(u) if x);mu=null[index]/u[index]
  assert mu>0 and all(null[d]==mu*u[d] for d in range(22))
 assert sorted(c['rank'] for c in raw['cusp_components'])==[2,7,11]
 assert {r['negative_wall_one_based'] for r in raw['families']}=={60,99,125}
 for row in raw['families']:certify(row,Q,W,G,h)
 forged=deepcopy(raw['families'][0]);forged['positive_face']['point'][0]=str(F(forged['positive_face']['point'][0])+1)
 try:certify(forged,Q,W,G,h)
 except AssertionError:pass
 else:raise AssertionError('forged ambient face accepted')
 result={'status':'PASS','standard_library_only':True,'all_positive_parameter_families':3,'actual_forbidden_ridge_interval':'1/4 < cosine_squared < 1/2','whole_cusp_rank':19,'forged_face_rejected':True,
         'scope':'All mixed nonzero coefficients in the one displayed affine-A2 span, not all H20 sections. Zero-coefficient sections reduce to prior universal <=2-wall theorem.'}
 OUTPUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':run()
