#!/usr/bin/env python3
"""All-label two-low-projection anchor obstruction for the frozen mixed12."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from fresh_h9_h8_facet_cover_20261014_line_review import I,dec,dv,pair,s5
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_corner_anchors_20261014_pair_gate'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 base='results/fresh_h9_mixed_d4_coupling_20261014.json';gap='results/fresh_h9_h8_facet_cover_20261014_gap_gate.json'
 d,c=[json.loads((ROOT/p).read_text())for p in[base,gap]]
 for doc in[d,c]:
  for p,h in doc['sources'].items():assert sha(p)==h
 assert sorted(x['source_vertex']for x in c['exceptional_rays'])==[0,32,33]
 Q=[[dec(x)for x in row]for row in d['form']];N=[dv(v)for v in d['final_normals']];norm=[pair(v,v,Q)for v in N];assert all(x.l>0 for x in norm)
 a=(s5-1)/4;phi=(1+s5)/2;root2=I(2).sqrt();ba=(11+5*s5)/4;bb=(23+9*s5)/4
 def f(b,s,t):return b*s*t-((2+a*s*s)*(2+a*t*t)).sqrt()
 fa11=f(ba,I(1),I(1));fb11=f(bb,I(1),I(1))
 assert (ba*ba-a*a-2*a).l>0 and(bb*bb-a*a-2*a).l>0
 ordinary=[];special={}
 for i in range(12):
  for j in range(i+1,12):
   C=-2*pair(N[i],N[j],Q)/(norm[i]*norm[j]).sqrt()
   if(i,j)in[(4,11),(8,11)]:special[i,j]=C
   else:assert C.h<fa11.l;ordinary.append(dict(pair=[i,j],normalized_negative_product=C.enc(),strict_minimum_anchor_margin=(fa11-C).enc()))
 C4=special[4,11];C8=special[8,11]
 cases={
  'pair8_11_A_lower':C8-fa11,
  'pair8_11_A_next':f(ba,I(1),root2)-C8,
  'pair4_11_A_no_scale1':C4-f(ba,I(1),I(2)),
  'pair4_11_A_44_lower':C4-f(ba,root2,root2),
  'pair4_11_A_45_upper':f(ba,root2,phi)-C4,
  'pair4_11_B_lower':C4-fb11,
  'pair4_11_B_next':f(bb,I(1),root2)-C4,
  'pair8_11_B_below_minimum':fb11-C8,
 }
 assert all(x.l>0 for x in cases.values())
 # Reconstruct the3 exception Gram entries from contact vectors by a
 # separate exact symbolic linear solve in the source Gram; no saved b used.
 import sympy as s
 srcpath='data/known_examples/bugaenko_compact_h8.json';src=json.loads((ROOT/srcpath).read_text());S5=s.sqrt(5)
 def sx(x):return s.Rational(x['a'])+s.Rational(x['b'])*(1+S5)/2
 G=s.Matrix([[sx(x)for x in row[:9]]for row in src['wall_gram'][:9]])
 hs=[s.Matrix([s.Rational(x[0])+s.Rational(x[1])*S5 for x in row['contacts'][:9]])for row in c['exceptional_rays']]
 gi=G.inv();vals=[[s.simplify((x.T*gi*y)[0])for y in hs]for x in hs]
 expected=[[(1-S5)/4,-(11+5*S5)/4,-(11+5*S5)/4],[-(11+5*S5)/4,(1-S5)/4,-(23+9*S5)/4],[-(11+5*S5)/4,-(23+9*S5)/4,(1-S5)/4]]
 assert all(s.simplify(vals[i][j]-expected[i][j])==0 for i in range(3)for j in range(3))
 paths=[base,gap,srcpath,'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py',f'scripts/{TAG}.py']
 return dict(status='PASS_ALL_LABEL_TWO_EXCEPTIONAL_CENTER_ANCHORS_EXCLUDED',sources={p:sha(p)for p in paths},
  old_walls=12,source_exceptional_vertices=[0,32,33],source_exceptional_gram=[[str(x)for x in row]for row in vals],
  ordinary_pairs=ordinary,special_normalized_negative_products={f'{i},{j}':x.enc()for(i,j),x in special.items()},
  strict_interval_margins={k:v.enc()for k,v in cases.items()},
  conclusion='For any all-dotted cap direction retaining this mixed12, any coherent perpendicular placement of the fixedH8source has at mostONE projected old center with normalized kappa²<sqrt5−1. Allfinite labels included.',
  scope='Necessary anchor obstruction; it does not exclude single exceptional-center anchors, all-high-norm projections, different source polytopes, nonperpendicular attachments or changed old walls. No cap vector or finite-volume construction is asserted.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['ordinary_pairs']),'ordinarypairs',flush=True)
