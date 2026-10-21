#!/usr/bin/env python3
"""Certify all positive two-parameter sections through the sampled A2 cusp span."""
import json
import numpy as np
import sympy as sp
import fresh_h20_new_route_20261013 as old
ROOT=old.ROOT;OUT=ROOT/'data/enumerations/fresh_h20_joint_cusp_a2_20261013.json'
s,t=sp.symbols('s t')
def enc(p):return [[int(a),int(b),str(c)] for (a,b),c in sp.Poly(p,s,t,domain=sp.QQ).terms()]
def face(Q,W,G,h,ids):
 M=sp.Matrix(G[np.ix_(ids,ids)]);assert M.is_positive_definite
 z=sp.Matrix(h)-sp.Matrix(W[ids].T)*M.inv()*sp.Matrix(W[ids]@Q@h)
 vals=sp.Matrix(W)*sp.Matrix(Q)*z
 assert [i for i,x in enumerate(vals) if not x]==sorted(ids) and max(vals)==0
 assert (z.T*sp.Matrix(Q)*z)[0]<0
 return z

def build():
 Q,W,G,h=old.base.load();support=[59,98,124];rows=[]
 for k,a,b in [(59,2,26),(98,8,7),(124,0,21)]:
  i,j=[x for x in support if x!=k];e=sp.Matrix(W[i])+s*sp.Matrix(W[j])-t*sp.Matrix(W[k]);q=sp.expand((e.T*sp.Matrix(Q)*e)[0]);v=sp.Matrix(W)*sp.Matrix(Q)*e
  gram=lambda l,m:sp.expand(q*G[l,m]-v[l]*v[m])
  zplus=face(Q,W,G,h,[a,b,i,j]);zminus=face(Q,W,G,h,[a,b,k]);up=sp.expand((zplus.T*sp.Matrix(Q)*e)[0]);um=sp.expand((zminus.T*sp.Matrix(Q)*e)[0]);point=up*zminus-um*zplus
  assert sp.expand((point.T*sp.Matrix(Q)*e)[0])==0
  vals=sp.Matrix(W)*sp.Matrix(Q)*point
  assert [l for l,x in enumerate(vals) if sp.expand(x)==0]==sorted([a,b])
  for l in range(210):
   if l not in (a,b):assert all(c<=0 for c in sp.Poly(vals[l],s,t).coeffs())
  aa,ab,bb=gram(a,a),gram(a,b),gram(b,b);N=sp.expand(ab**2);den=sp.expand(aa*bb)
  A=s*s+s*t-s+t*t+t+1;B=2*s*s+2*s*t-2*s+t*t+2*t+2
  assert sp.expand(2*B*N-A*den)==0
  assert sp.expand(4*N-den-8*q*t*t)==0
  assert sp.expand(den-2*N-8*q*((s-sp.Rational(1,2))**2+sp.Rational(3,4)+s*t+t))==0
  cplus=face(Q,W,G,h,[i,j]);cminus=face(Q,W,G,h,[k]);cpair=sp.expand((cplus.T*sp.Matrix(Q)*e)[0]);cm=sp.expand((cminus.T*sp.Matrix(Q)*e)[0]);center=cpair*cminus-cm*cplus
  rows.append({'strict_center_faces':[{'walls_one_based':[i+1,j+1],'point':list(map(str,cplus))},{'walls_one_based':[k+1],'point':list(map(str,cminus))}],
    'strict_center':[enc(x) for x in center],'positive_walls_one_based':[i+1,j+1],'negative_wall_one_based':k+1,'ridge_walls_one_based':[a+1,b+1],
    'positive_face':{'walls_one_based':[a+1,b+1,i+1,j+1],'point':list(map(str,zplus))},
    'negative_face':{'walls_one_based':[a+1,b+1,k+1],'point':list(map(str,zminus))},
    'face_normal_pairings':[enc(up),enc(um)],'ridge_point':[enc(x) for x in point],
    'ridge_square':enc((point.T*sp.Matrix(Q)*point)[0]),'normal_square':enc(q),
    'projected_gram':[enc(x) for x in (aa,ab,bb)],'status':'ALL_POSITIVE_S_T_NON_COXETER_ACTUAL_RIDGE'})
 cusp=json.loads(old.CUSP_SOURCE.read_text())['cusps'][1]['vector'];u=np.asarray(cusp,dtype=object);ids=[i for i,p in enumerate(W@Q@u) if not p]
 comps=[]
 for cs in old.components(G,ids):
  marks=old.positive_kernel(G[np.ix_(cs,cs)]);comps.append({'walls_one_based':[i+1 for i in cs],'marks':marks,'rank':len(cs)-1})
 assert np.array_equal(sum(W[support]),2*u)
 result={'schema':'universal-affine-A2-cusp-section-v1','source':str(old.base.SOURCE.relative_to(ROOT)),
   'cusp_vector':cusp,'cusp_components':comps,'support_one_based':[i+1 for i in support],
   'families':rows,'scope':'Every interior section with normal in the rational or real span of the three displayed walls. Nonzero mixed coefficients are covered up to overall sign and scaling; zero coefficients use the previous universal two-wall obstruction. No other cusp normal span is excluded.'}
 OUT.write_text(json.dumps(result,indent=2)+'\n');print('Archived three all-positive-parameter ridge certificates')
if __name__=='__main__':build()
