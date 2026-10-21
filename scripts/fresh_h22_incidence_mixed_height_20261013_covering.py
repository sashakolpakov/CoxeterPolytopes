#!/usr/bin/env python3
"""Stdlib-only check of an untouched old cusp and a simpler strict escape.

This is independent of the construction arithmetic. It verifies the supplied
literal null ray and all projection pairings directly with integer arithmetic.
"""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_incidence_mixed_height_20261013'
BASE='data/enumerations/'+PRE+'.json';GAP='data/enumerations/'+PRE+'_gap.json';RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';OUT='data/enumerations/'+PRE+'_covering.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mv(A,x):return[dot(r,x)for r in A]
def rank(A):
 A=[[F(x)for x in row]for row in A];r=0
 for j in range(len(A[0])):
  k=next((i for i in range(r,len(A))if A[i][j]),None)
  if k is None:continue
  A[r],A[k]=A[k],A[r];p=A[r][j];A[r]=[x/p for x in A[r]]
  for i in range(r+1,len(A)):
   p=A[i][j];A[i]=[x-p*y for x,y in zip(A[i],A[r])]
  r+=1
  if r==len(A):break
 return r
def run():
 b=read(BASE);g=read(GAP);raw=read(RAW)['coxeter_input']
 for d in(b,g):
  for p,s in d['sources'].items():assert sha(p)==s
 Q=[[-int(x)for x in r]for r in raw['basis_gram']];W=[[int(x)for x in r]for r in raw['wall_coordinates_as_columns_transposed']]
 D=b['null_projections'];X2=b['sqrt3_twice_projections'];u=g['new_old_projection'];h=b['center']
 delta=[20,-47,-32,-63,-93,-76,-58,-40,-20,80,52,106,160,135,111,87,64,-20,-40,-29,-16,42]
 cd=mv(Q,delta);ch=mv(Q,h);assert dot(delta,cd)==0 and dot(h,ch)==-14 and dot(h,cd)==-8
 pw=[dot(w,cd)for w in W];assert max(pw)==0;I=[i for i,p in enumerate(pw)if p==0]
 assert len(I)==24
 assert max(dot(w,cd)for w in D)==-2 and max(dot(w,cd)for w in X2)==-4 and dot(u,cd)==-111
 # Every new wall is strict at(delta,0). The added-coordinate vector has
 # norm1 and is perpendicular to delta and to every incident old wall.
 G=[[dot(W[i],mv(Q,W[j]))for j in I]for i in I]
 nr=rank([W[i]for i in I]);gr=rank(G);assert(nr,gr)==(21,20)
 v=[100*a+c for a,c in zip(delta,h)];z=100;cv=mv(Q,v)
 assert dot(v,cv)+z*z==8386 and dot(h,cv)==-814
 op=max(dot(w,cv)for w in W);dp=max(dot(w,cv)for w in D);xp=max(dot(w,cv)for w in X2);up=dot(u,cv)
 assert(op,dp,xp,up)==(-1,-206,-414,-11260)
 # a^2=(9+6sqrt2)/4>4: 6sqrt2>7 follows from72>49.
 assert 72>49
 assert(2*dp)**2>2*z*z and xp*xp>12*z*z and up*up>22*z*z
 return {'status':'PASS_INDEPENDENT_UNTOUCHED_CUSP_AND_STRICT_ESCAPE_ALL2296','sources':{p:sha(p)for p in(BASE,GAP,RAW)},'old_null_ray':delta,'old_null_square':0,'old_center_pair':-8,'incident_old_walls':I,'incident_normal_rank':nr,'incident_gram_rank':gr,'new_projection_maxima':{'delta_bank':-2,'twice_x_bank':-4,'sqrt22_u':-111},'new_caps_incident':0,'link_recession':'The unit added-coordinate vector is orthogonal to every incident normal and to(delta,0), and has square1. Thus the horospherical quotient has a free vertical direction.','independent_strict_escape':{'old_coordinates':v,'last_coordinate':100,'square':8386,'center_pair':-814,'max_old_pair':op,'max_delta_pair':dp,'max_twice_x_pair':xp,'sqrt22_projection_pair':up},'all_finite_timelike_additions_obstruction':'At any old future null ray delta, every future timelike old projection has strictly negative pairing. Any finite collection of such additional caps remains strict at(delta,0); small vertical displacement plus a small move toward a strict center gives a strict spacelike escape. New null projections not proportional to delta are also strict there.','scope':'Explicit missing old cusp in allthree1336 chambers; independent integer verification covers the full2296 relaxed union. This does not exclude a new null cap proportional to delta, pure vertical walls, or replacement of old facets.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
