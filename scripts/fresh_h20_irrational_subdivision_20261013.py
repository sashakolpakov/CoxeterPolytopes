#!/usr/bin/env python3
"""Exact local obstructions to finite face-to-face Coxeter subdivisions.

All 59 frozen twin-leaf section planes are covered, with no parameter or
angle-label cutoff.  The geometric implications are in the matching report.
This uses the prior elementary Fraction matrix helpers, not a numerical scout.
"""
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
import argparse, hashlib, json
import fresh_h20_cusp_sections_20261013 as P

ROOT=P.ROOT
TAG='fresh_h20_irrational_subdivision_20261013'
OLD=f'results/{P.TAG}.json'
OUT=ROOT/f'results/{TAG}.json'
def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F): return str(x)
 if isinstance(x,(list,tuple)): return [enc(y) for y in x]
 if isinstance(x,dict): return {str(k):enc(v) for k,v in x.items()}
 return x
def qadd(x,y): return (x[0]+y[0],x[1]+y[1])
def qmul(x,y): return (x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def qdiv(x,y):
 d=y[0]**2-5*y[1]**2
 return qmul(x,(y[0]/d,-y[1]/d))
def tile_data():
 one=(F(1),F(0)); B=(F(1,4),F(0)); out={}
 for typ,C in [('A3',B),('B3',(F(1,2),F(0))),('H3',(F(3,8),F(1,8)))]:
  sb=qadd(one,(-B[0],-B[1]));sc=qadd(one,(-C[0],-C[1]))
  out[typ]=[qdiv(qmul(B,C),qmul(sb,sc)),qdiv(B,sc),qdiv(C,sb)]
 assert set(out['A3'])=={(F(1,9),F(0)),(F(1,3),F(0))}
 assert set(out['B3'])=={(F(1,3),F(0)),(F(1,2),F(0)),(F(2,3),F(0))}
 assert set(out['H3'])=={(F(1,3),F(2,15)),(F(1,2),F(1,10)),(F(1,2),F(1,6))}
 assert all(x[1] for x in out['H3'])
 available=[F(1,2),F(1,3),F(1,4),F(1,5)]
 sums=sorted({a+b for a in available for b in available if F(1,3)<a+b<F(1,2)})
 assert sums==[F(2,5),F(9,20)]
 assert 3*F(1,5)>F(1,2)
 assert (F(2,5)+F(1,4)+F(1,2)-1)/(F(1,2)+F(1,3)+F(1,5)-1)==F(9,2)
 return {'side_cosine_squares_in_Qsqrt5':out,'possible_exceptional_corner_angles_over_pi':sums,
         'all_H3_area_ratio_at_2pi_over5':F(9,2)}

def run():
 raw=json.loads((ROOT/P.RAW).read_text())['coxeter_input']
 Q=[[-int(x) for x in r] for r in raw['basis_gram']]
 W=[list(map(int,r)) for r in raw['wall_coordinates_as_columns_transposed']]
 C=[P.mv(Q,w) for w in W];G=[[P.dot(w,c) for c in C] for w in W]
 old=json.loads((ROOT/OLD).read_text());h=list(map(F,old['strict_center']))
 assert old['status']=='PASS_ALL_REAL_FORK_LEAF_CUSP_SECTION_PENCILS_EXCLUDED'
 assert old['sources'][P.RAW]==sha(P.RAW)
 assert old['sources'][f'scripts/{P.TAG}.py']==sha(f'scripts/{P.TAG}.py')
 assert len(W)==210 and P.pair(h,h,Q)==-14
 assert all(P.dot(h,c)<0 for c in C)
 faces={}
 def face(ids):
  key=tuple(sorted(ids))
  if key not in faces: faces[key]=P.project(h,list(key),W,G,Q,C)
  return faces[key]
 def actual(row,js):
  a,b=row['retained_leaf'],row['omitted_leaf'];delta=list(map(F,row['cusp']))
  e=list(map(F,row['base_normal']));cs=[i for i,c in row['endpoint_positive_combination']]+[b]
  left=sorted([a,b]+js);right=sorted(set(cs)-{b}|set(js))
  zL,cL=face(left);zR,cR=face(right)
  A=-P.pair(delta,zL,Q);B=-P.pair(delta,zR,Q)
  assert A>0 and B>0 and P.pair(e,zL,Q)==0
  assert P.pair(P.add(e,P.scale(F(1,2),delta)),zR,Q)==0
  assert set(left)&set(right)==set([a]+js)
  # x(t)=(1/2-t) B zL+t A zR is strictly future timelike for 0<t<1/2.
  # Its section pairing cancels identically, and its exact zero set is the
  # intersection of the two explicitly checked ambient face zero sets.
  for i in range(210):
   pl=P.dot(zL,C[i]);pr=P.dot(zR,C[i])
   assert pl<=0 and pr<=0
   assert (pl==pr==0)==(i in set([a]+js))
  return {'left_face':cL,'right_face':cR,'A':A,'B':B,
          'actual_section_zero_walls_for_0_lt_t_lt_half':sorted([a]+js),
          'interpolation':'(1/2-t)*B*zL+t*A*zR',
          't_zero_extra_wall':b,'t_zero_duplicate_projected_normal':a}
 rows=[];planes={};counts=Counter()
 for ri,row in enumerate(old['pencils']):
  a,b,j=row['retained_leaf'],row['omitted_leaf'],row['auxiliary_wall']
  delta=list(map(F,row['cusp']));e=list(map(F,row['base_normal']))
  cs=[i for i,c in row['endpoint_positive_combination']]+[b]
  assert G[a][a]==G[b][b]==2 and G[a][b]==0
  assert P.pair(e,e,Q)==1 and P.pair(e,delta,Q)==P.pair(delta,delta,Q)==0
  assert e==P.scale(F(1,2),P.add(W[a],P.scale(-1,W[b])))
  end=row['endpoint_positive_combination']
  assert min(F(v) for i,v in end)>0
  assert P.add(e,P.scale(F(1,2),delta))==[sum(F(v)*W[i][k] for i,v in end) for k in range(22)]
  assert all(P.pair(delta,W[i],Q)==0 for i in cs)
  assert G[j][j]==2 and P.dot(delta,C[j])==-1
  nu=P.dot(e,C[j]);assert nu in [F(-1,2),F(1,2)]
  assert G[a][j]-nu==F(-1,2)
  common=[c for c in cs if c not in [a,b] and G[a][c]!=0]
  assert len(common)==1;c=common[0]
  assert G[a][c]==G[b][c] and G[c][j]==0 and P.dot(e,C[c])==0
  rec={'source_row':ri,'cusp_orbit':row['cusp_orbit'],'retained_leaf':a,'omitted_leaf':b,
       'cusp':delta,'base_normal':e,'common_neighbor':c}
  if G[c][c]==2 and G[a][c]==-1:
   counts['spherical_triangle_orientations']+=1
   rec.update({'method':'SPHERICAL_TRIANGLE_NO_FINITE_COXETER_TILING',
       'auxiliary_wall':j,'actual_codimension3_face':actual(row,[c,j]),
       'projected_Gram_order':[a,c,j],
       'projected_Gram_polynomials':[[[1],[-1],[F(-1,2),1]],
          [[-1],[2],[0]],[[F(-1,2),1],[0],[F(7,4),2*nu,-1]]],
       'link_angles':['pi/4','pi/2','theta'],
       'theta_cosine_squared':{'numerator':[F(1,4),-1,1],
                                'denominator':[F(7,4),2*nu,-1]}})
   # Denominator-4*numerator is concave, strictly positive on [0,1/2].
   assert F(3,4)>0 and F(7,4)+nu-F(1,4)>0
   # Three-wall Gram determinant D-2(1/2-t)^2: concave, positive endpoints.
   assert F(5,4)>0 and F(3,2)+nu>0
  else:
   assert G[c][c]==4 and G[a][c]==-2 and len(cs)==3
   counts['cyclotomic_two_ridge_orientations']+=1
   def pick(q,nu,d,g):
    return next(k for k in range(210) if k not in cs and
      (G[k][k],P.dot(e,C[k]),P.dot(delta,C[k]),G[a][k])==(q,nu,d,g))
   jA=pick(2,F(1,2),-1,0);jD=pick(4,1,-4,0)
   rec.update({'method':'TWO_ACTUAL_RIDGES_CYCLOTOMIC_NORM_CONTRADICTION',
      'A_ridge':actual(row,[jA]),'D_ridge':actual(row,[jD]),
      'A_auxiliary_wall':jA,'D_auxiliary_wall':jD,
      'A_4cosine_squared':{'numerator':[1,-4,4],'denominator':[F(7,4),1,-1]},
      'D_4cosine_squared':{'numerator':[4,-32,64],'denominator':[3,8,-16]}})
  key=P.plane_key(e,delta);planes.setdefault(key,[]).append(ri);rows.append(rec)
 assert len(rows)==118 and len(planes)==59 and all(len(v)==2 for v in planes.values())
 assert counts=={'spherical_triangle_orientations':114,'cyclotomic_two_ridge_orientations':4}
 for ids in planes.values():
  u,v=(rows[i] for i in ids)
  assert u['method']==v['method'] and u['cusp']==v['cusp']
  assert u['base_normal']==P.scale(-1,v['base_normal'])
 lo,hi=F(-1,8),F(5,8)
 # sqrt(2)<3/2 gives strict containment in this rational interval.
 assert 2<F(9,4)
 dp=lambda t:F(3,4)+5*t-5*t*t
 den=lambda t:F(7,4)+t-t*t
 assert min(dp(lo),dp(hi))==F(3,64)>0 and min(den(lo),den(hi))>0
 sources=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',OLD,P.RAW]
 return enc({'status':'PASS_ALL_59_SECTION_PLANES_NO_FINITE_FACE_TO_FACE_COXETER_SUBDIVISION',
    'sources':{p:sha(p) for p in sources},'counts':dict(counts),
    'unoriented_plane_count':59,'parameter_scope':'all real t giving a full-dimensional section; t=0 included',
    'rows':rows,'spherical_tile_arithmetic':tile_data(),
    'exceptional_two_plane_interval':{'rational_outer_interval':[lo,hi],
       'denominator_minus_A_numerator_endpoint_values':[dp(lo),dp(hi)],
       'denominator_endpoint_values':[den(lo),den(hi)]},
    'scope':'Finite face-to-face subdivision of the entire fixed section into full-dimensional Coxeter polytopes; no arbitrary changed section or non-face-to-face claim.'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 result=run()
 if args.verify: assert result==json.loads(OUT.read_text())
 else: OUT.write_text(json.dumps(result,indent=2)+'\n')
 print(result['status']);print(result['counts'])
