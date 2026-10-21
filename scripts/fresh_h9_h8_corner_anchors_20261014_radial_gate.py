#!/usr/bin/env python3
"""Exact anchored radial-cover obstruction for one isolated cap direction.

Reconstructs the physical ball geometry with 192-bit rational intervals.
The source contact-gap classification is a named frozen dependency.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as s
from fresh_h9_h8_facet_cover_20261014_line_review import I,dec,dv,pair,s5,BITS

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_h8_corner_anchors_20261014_radial_gate'
CAP='results/fresh_h9_h8_global_cap_gap_20261014.json'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
GAP='results/fresh_h9_h8_facet_cover_20261014_gap_gate.json'
SOURCE='data/known_examples/bugaenko_compact_h8.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for path,h in d.get('sources',{}).items():assert sha(path)==h,path
 return d
def run():
 cap,b,gap,src=[load(p)for p in[CAP,BASE,GAP,SOURCE]]
 assert cap['unique_low_wall']==8
 assert sorted(r['source_vertex']for r in gap['exceptional_rays'])==[0,32,33]
 assert gap['exceptional_minimum_norm']==['-1/4','1/4']
 assert gap['ordinary_minimum_norm']==['-1','1']
 Q=[[dec(x)for x in row]for row in b['form']]
 N=[dv(v)for v in b['final_normals']];H=dv(b['center'])
 recipe=cap['cap_recipe'];lo,hi=map(F,recipe['parameter_interval']);assert lo<hi
 t=I(lo,hi);V=[I(x)for x in recipe['rational_endpoint']]
 delta=[v-h for v,h in zip(V,H)];z=[h+t*d for h,d in zip(H,delta)]
 a=(s5-1)/4;q=pair(z,z,Q);assert q.l>0
 g=[pair(z,n,Q)for n in N];nn=[pair(n,n,Q)for n in N]
 assert all(x.h<0 for x in g)and all(x.l>0 for x in nn)
 assert pair(z,H,Q).h<0
 # Reconstruct the defining polynomial, and bind its unique root independently
 # of the 128-bit coordinate boxes in the cap certificate.
 q0,q1,q2=pair(H,H,Q),2*pair(H,delta,Q),pair(delta,delta,Q)
 g0,g1=pair(H,N[8],Q),pair(delta,N[8],Q)
 polys=[g0*g0-(2+a)*q0,2*g0*g1-(2+a)*q1,g1*g1-(2+a)*q2]
 for p,old in zip(polys,recipe['polynomial_ascending']):
  v=dec(old);assert max(p.l,v.l)<=min(p.h,v.h)
 def pe(x):return polys[0]+x*polys[1]+x*x*polys[2]
 end=[pe(lo),pe(hi)]
 assert end[0].h<0<end[1].l or end[1].h<0<end[0].l
 deriv=[polys[1]+2*x*polys[2]for x in[lo,hi]]
 assert all(x.l>0 for x in deriv)or all(x.h<0 for x in deriv)
 k=[g[i]*g[i]/q-nn[i]for i in range(12)]
 assert all(x.l>0 for x in k)
 assert max(k[8].l,a.l)<=min(k[8].h,a.h)
 # n8 has square2 and the defining polynomial makes k8 exactly a.
 assert nn[8].l<=2<=nn[8].h
 assert all((2*k[i]/nn[i]-4*a).l>0 for i in range(12)if i!=8)
 k[8]=a
 balls=[]
 for i in range(12):
  cr=((k[i]+nn[i])/k[i]).sqrt();sr=(nn[i]/k[i]).sqrt()
  if i==8:dc,ds=I(1),I(0)
  else:
   dc=(g[8]*g[i]/q-pair(N[8],N[i],Q))/(a*k[i]).sqrt()
   assert dc.l>1;ds=(dc*dc-1).sqrt()
  reach=dc*cr+ds*sr
  balls.append(dict(old_wall=i,center_distance_cosh=dc.enc(),radius_cosh=cr.enc(),
                    anchored_reach_cosh=reach.enc()))
 # Exact source ray Gram from the literal source Gram and contact vectors.
 S5=s.sqrt(5);sf=lambda x:s.Rational(x['a'])+s.Rational(x['b'])*(1+S5)/2
 SG=s.Matrix([[sf(x)for x in row[:9]]for row in src['wall_gram'][:9]])
 ray={r['source_vertex']:r for r in gap['exceptional_rays']};Ginv=SG.inv()
 hs={v:s.Matrix([s.Rational(x[0])+s.Rational(x[1])*S5 for x in r['contacts'][:9]])for v,r in ray.items()}
 for v in[0,32,33]:
  assert ray[v]['zero_set']==src['geometry']['vertices'][v]
  assert s.simplify((hs[v].T*Ginv*hs[v])[0]+(S5-1)/4)==0
 target_pairs=[(0,32),(32,33),(33,32)];rows=[]
 for v,w in target_pairs:
  product=s.simplify((hs[v].T*Ginv*hs[w])[0])
  bval=(11+5*s5)/4 if v==0 else(23+9*s5)/4
  expected=-(11+5*S5)/4 if v==0 else-(23+9*S5)/4
  assert s.simplify(product-expected)==0
  dc=bval/a;margins=[]
  for ball in balls:
   reach=I(*map(F,ball['anchored_reach_cosh']));margin=dc-reach
   assert margin.l>0;margins.append(margin.enc())
  rows.append(dict(anchor_source_vertex=v,uncovered_source_vertex=w,
                   source_distance_cosh=dc.enc(),strict_cosh_margins=margins))
 # The maximum is identified only for the explanatory diagnostic, not needed
 # for the 36 individual strict comparisons.
 assert all(F(balls[7]['anchored_reach_cosh'][0])>F(row['anchored_reach_cosh'][1])for i,row in enumerate(balls)if i!=7)
 paths=[CAP,BASE,GAP,SOURCE,'scripts/fresh_h9_h8_global_cap_gap_20261014.py',
        'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py',f'scripts/{TAG}.py']
 return dict(status='PASS_SINGLE_ANCHOR_CAP_ALL_THREE_SOURCE_PLACEMENTS_EXCLUDED',
  sources={p:sha(p)for p in paths},interval_bits=BITS,
  parameter_interval=[str(lo),str(hi)],root_endpoint_bounds=[x.enc()for x in end],
  derivative_endpoint_bounds=[x.enc()for x in deriv],raw_cap_norm=q.enc(),
  anchor_old_wall=8,anchor_exact_normalized_projection_norm=['-1/4','1/4'],
  balls=balls,maximum_reach_old_ball=7,anchor_obstructions=rows,strict_comparisons=36,
  conclusion='For this exact cap direction, the raw attachment formed by the old12 walls, this cap and any coherent perpendicular placement of the fixed compact11-facet H8 source has infinite volume: every possible source-corner anchor leaves a source vertex outside all twelve closed old balls.',
  scope='One isolated algebraic cap direction, fixed old12 chamber and fixed H8 source. Inherits the exact source contact-gap classification and warped-end necessary coverage lemma. Additional nonperpendicular walls can cut the uncovered end and are outside this raw-attachment theorem. It also does not exclude other cap directions, changed old walls or different H8 sources. No numerical optimizer is a proof dependency.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run()
 p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['strict_comparisons'],'strict comparisons',flush=True)
