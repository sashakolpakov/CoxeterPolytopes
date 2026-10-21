#!/usr/bin/env python3
"""Independent source-ray and192-bit interval review of the radial gate."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
from fresh_h9_h8_facet_cover_20261014_line_review import I,dec,dv,pair,s5,BITS

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_global_cap_gap_20261014_radial_review'
PRIMARY='results/fresh_h9_h8_global_cap_gap_20261014_radial.json'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
SOURCE='data/known_examples/bugaenko_compact_h8.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/PRIMARY).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 base=json.loads((ROOT/BASE).read_text());src=json.loads((ROOT/SOURCE).read_text())
 S5=s.sqrt(5);R=s.sqrt((16+4*S5)/11)
 K=s.QQ.algebraic_field(S5,R)
 def ex(v):return s.Rational(v[0][0])+s.Rational(v[0][1])*S5+(s.Rational(v[1][0])+s.Rational(v[1][1])*S5)*R
 def kv(v):return K.from_sympy(ex(v))
 Qk=[[kv(x)for x in row]for row in base['form']]
 Nk=[[kv(x)for x in row]for row in base['final_normals']]
 def kp(x,y):return sum((x[i]*sum((Qk[i][j]*y[j]for j in range(10)),K.zero)for i in range(10)),K.zero)
 nQ=[[dec(x)for x in row]for row in base['form']];nN=[dv(v)for v in base['final_normals']]
 norms=[pair(v,v,nQ)for v in nN];nk=[kp(v,v)for v in Nk]
 assert all(x.l>0 for x in norms)
 graph=[];products=[];vals={}
 for i,j in combinations(range(12),2):
  exact=kp(Nk[i],Nk[j]);test=4*exact*exact-nk[i]*nk[j]
  g=pair(nN[i],nN[j],nQ)
  if exact==K.zero:g=I(0)
  else:assert g.h<0
  C=-2*g/(norms[i]*norms[j]).sqrt();vals[i,j]=C
  if test==K.zero:assert C.l<=1<=C.h;sign=0
  elif C.l>1:sign=1;graph.append([i,j])
  else:assert C.h<1;sign=-1
  products.append(dict(pair=[i,j],comparison_to_one=sign,normalized_negative_product=C.enc()))
 assert graph==[[0,1],[4,11],[5,6],[8,11]]==d['greater_than_one_graph']
 assert [r['comparison_to_one']for r in products]==[r['comparison_to_one']for r in d['old_pair_products']]
 # Rebuild each source ray as the literal8-row physical nullspace, then
 # normalize its minimum positive source contact to1.
 sf=lambda x:s.Rational(x['a'])+s.Rational(x['b'])*(1+S5)/2
 SQ=s.Matrix([[sf(x)for x in row]for row in src['lattice_gram']])
 SN=s.Matrix([[sf(x)for x in row]for row in src['simple_roots']])
 SH=s.Matrix([sf(x)for x in src['geometry']['interior_point']]);rays={};records=[]
 simp=lambda x:s.simplify(s.expand(x))
 for vi in[0,32,33]:
  ids=src['geometry']['vertices'][vi];M=SN[ids,:]*SQ;ker=M.nullspace();assert len(ker)==1
  v=ker[0];orient=simp((v.T*SQ*SH)[0]);assert orient!=0
  if orient>0:v=-v
  contacts=[simp(x)for x in -SN*SQ*v];assert all(x>=0 for x in contacts)
  positive=[x for x in contacts if x>0];m=min(positive);assert m>0
  v=v.applyfunc(lambda x:simp(x/m));contacts=[simp(x/m)for x in contacts]
  assert simp((v.T*SQ*v)[0]+(S5-1)/4)==0
  old=d['source_exceptional_points'][str(vi)]
  assert all(simp(v[i]-ex(old['point'][i]))==0 for i in range(9))
  assert all(simp(contacts[i]-ex(old['contacts'][i]))==0 for i in range(11))
  assert [i for i,x in enumerate(contacts)if x==0]==ids
  rays[vi]=v;records.append(dict(source_vertex=vi,point=[str(x)for x in v],contacts=[str(x)for x in contacts]))
 sourcepairs=[]
 for i,j in combinations([0,32,33],2):
  C=simp(-(rays[i].T*SQ*rays[j])[0]/((S5-1)/4))
  assert simp(C-(9+4*S5 if i==0 else 17+8*S5))==0
  sourcepairs.append(dict(pair=[i,j],cosh=str(C)))
 a=(s5-1)/4;B=4*a;DA=9+4*s5;DB=17+8*s5;phi=(1+s5)/2
 def reach(x,y,C):
  dc=(((x+2)*(y+2)).sqrt()+C)/(x*y).sqrt()
  return dc*((y+2)/y).sqrt()+(dc*dc-1).sqrt()*(2/y).sqrt()
 def union(x,C):
  dc=(x+2+C)/x;cr=((x+2)/x).sqrt();sr=(2/x).sqrt()
  return dc*(cr*cr+sr*sr)+2*(dc*dc-1).sqrt()*cr*sr
 bounds=dict(ordinary_neighbor_max_reach_cosh=reach(a,B,I(1)),
             own_ball_max_radius_cosh=((a+2)/a).sqrt(),high_ball_max_diameter_cosh=1+4/B,
             high_balls_4_8_union_diameter_cosh=union(B,I(0)),
             scale_at_least_phi_max_reach_cosh=reach(a*phi*phi,B,vals[4,11]),
             sqrt2_scale_ball8_max_reach_cosh=reach(2*a,B,vals[8,11]))
 margins={}
 scale=int(d['interval_scale'])
 for k,value in bounds.items():
  target=DA if k in['ordinary_neighbor_max_reach_cosh','own_ball_max_radius_cosh','high_ball_max_diameter_cosh']else DB
  margin=target-value;assert margin.l>0;margins[k]=margin.enc()
  old=I(F(int(d[k][0]),scale),F(int(d[k][1]),scale))
  assert max(old.l,value.l)<=min(old.h,value.h)
 assert vals[4,11].l>vals[8,11].h
 assert kp(Nk[4],Nk[8])==K.zero
 # Independent analytic checks behind uniformity over all cap directions.
 x,y,c=s.symbols('x y c',positive=True)
 f=(s.sqrt((x+2)*(y+2))+c)/s.sqrt(x*y)
 dx=-(2*s.sqrt(y+2)/s.sqrt(x+2)+c)/(2*x*s.sqrt(x*y))
 assert s.simplify(s.diff(f,x)-dx)==0
 assert s.simplify((1+S5)**2/4-(3+S5)/2)==0
 assert d['eliminated_low_old_walls']==list(range(11))and d['eliminated_low11_source_corner']==0
 assert [(r['source_corner'],r['scale'])for r in d['surviving_low_branches']]==[(32,'1'),(32,'sqrt2'),(33,'1'),(33,'sqrt2')]
 paths=[PRIMARY,'scripts/fresh_h9_h8_global_cap_gap_20261014_radial.py',BASE,SOURCE,
  'results/fresh_h9_h8_facet_cover_20261014_gap_gate.json',
  'results/fresh_h9_h8_corner_anchors_20261014_pair_gate.json',
  'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py',f'scripts/{TAG}.py']
 return dict(status='PASS_INDEPENDENT_UNIFORM_RADIAL_REVIEW',sources={p:sha(p)for p in paths},
  interval_bits=BITS,physical_pair_count=66,greater_than_one_graph=graph,
  independent_old_products=products,independent_source_rays=records,independent_source_pair_cosh=sourcepairs,
  independent_uniform_bounds={k:v.enc()for k,v in bounds.items()},strict_margins=margins,
  scope='No producer arithmetic imports. Separate SymPy ANP physical Gram, literal SymPy source nullspaces and192-bit rational intervals recheck66oldpairs, three source rays and allsix radial bounds. Source/proof review covers the two-remote-vertex argument and all-label scales; the exhaustive source-gap and at-most-one-low theorems remain inherited. Raw twelve-ball attachment scope, not arbitrary additional nonperpendicular walls.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run()
 p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
