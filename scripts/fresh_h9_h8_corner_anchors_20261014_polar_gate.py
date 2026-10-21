#!/usr/bin/env python3
"""Uniform low-polar anchor obstruction for raw perpendicular H8 attachment."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as s
import fresh_h9_mixed_d4_coupling_20261014_review as A
from fresh_h9_h8_facet_cover_20261014_line_review import I,dec,dv,pair,s5,BITS

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_h8_corner_anchors_20261014_polar_gate'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
GAP='results/fresh_h9_h8_facet_cover_20261014_gap_gate.json'
PAIR='results/fresh_h9_h8_corner_anchors_20261014_pair_gate.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for path,h in d.get('sources',{}).items():assert sha(path)==h,path
 return d
def run():
 b,gap,pg=[load(p)for p in[BASE,GAP,PAIR]]
 assert pg['source_exceptional_vertices']==[0,32,33]
 assert len(pg['ordinary_pairs'])==64 and len(pg['strict_interval_margins'])==8
 assert gap['exceptional_minimum_norm']==['-1/4','1/4']
 assert gap['ordinary_minimum_norm']==['-1','1']
 assert sorted(r['source_vertex']for r in gap['exceptional_rays'])==[0,32,33]
 # Literal physical Gram guards; zero identities are exact, not intervals.
 Q,N=A.dm(b['form']),A.dm(b['final_normals'])
 nn=[A.pair(v,v,Q)for v in N]
 assert all(nn[i]==2 for i in range(11))and nn[11].sign()>0
 zero_pairs=[[0,2],[0,4],[2,4],[0,11],[2,11]]
 assert all(not A.pair(N[i],N[j],Q)for i,j in zero_pairs)
 polar_right=[i for i in range(11)if not A.pair(N[11],N[i],Q)]
 assert polar_right==[0,1,2,3,5,6,7,9,10]
 g4=A.pair(N[11],N[4],Q);g8=A.pair(N[11],N[8],Q)
 assert g4.sign()<0 and g8.sign()<0
 C4=-dec(A.enc(g4))*(2/dec(A.enc(nn[11]))).sqrt()
 assert C4.l>F(44,5)and C4.h<9
 a=(s5-1)/4;B=4*a;L2=2*(B+2)
 DA=9+4*s5;DB=17+8*s5
 assert DB.l>DA.h and L2.l>6 and L2.h<F(169,25)
 # The three exceptional source rays all have norm-a. Their distinct pair
 # distances have cosh DA or DB, as rebuilt in the bound pair certificate.
 S5=s.sqrt(5)
 expected=[[(1-S5)/4,-(11+5*S5)/4,-(11+5*S5)/4],
           [-(11+5*S5)/4,(1-S5)/4,-(23+9*S5)/4],
           [-(11+5*S5)/4,-(23+9*S5)/4,(1-S5)/4]]
 assert all(s.simplify(s.sympify(pg['source_exceptional_gram'][i][j])-expected[i][j])==0 for i in range(3)for j in range(3))
 def reach(x,y,C):
  dc=(((x+2)*(y+2)).sqrt()+C)/(x*y).sqrt();assert dc.l>1
  return dc*((y+2)/y).sqrt()+(dc*dc-1).sqrt()*(2/y).sqrt()
 # A low anchor has kappa²>=a. Every other center has kappa²>=B by the
 # all-label at-most-one theorem. Right old pairs cannot cover either other
 # exceptional vertex; neither can the anchor's own ball.
 right_reach=reach(a,B,I(0))
 own_reach=((a+2)/a).sqrt()
 high_double_radius=1+4/B
 assert right_reach.h<DA.l and own_reach.h<DA.l and high_double_radius.h<DA.l
 # Hence both old balls4 and8 are required. In particular ball4 must reach
 # distance>=DA. If its ambient contact magnitude were>=4, kappa4²>=6.
 # The following upper bound uses C4<9, and is already too small.
 four_reach=reach(a,I(6),I(9));assert four_reach.h<DA.l
 # With cap square2, low11 gives contact A<13/5, while coverage gives
 # contact D<4. The Lorentz block {n11,n4} is orthogonal to n0,n2.
 # Both latter contacts have square>=L2>6, so their combined positive
 # projection cost is>6. Bound the negative Lorentz projection uniformly.
 AA,DD,CC=F(13,5),F(4),F(44,5)
 penalty=2*(AA*AA+DD*DD+CC*AA*DD)/(CC*CC-4)
 lower=6-penalty;assert penalty==F(2857,918)and lower==F(2651,918)>2
 # Separate symbolic checks of both monotonicity formulas and decomposition.
 x,y,c=s.symbols('x y c',positive=True)
 dc=(s.sqrt((x+2)*(y+2))+c)/s.sqrt(x*y)
 dx=-(2*s.sqrt(y+2)/s.sqrt(x+2)+c)/(2*x*s.sqrt(x*y))
 assert s.simplify(s.diff(dc,x)-dx)==0
 aa,dd,cc=s.symbols('A D C',positive=True)
 penalty_expr=2*(aa*aa+dd*dd+cc*aa*dd)/(cc*cc-4)
 deriv=-2*(aa*dd*(cc*cc+4)+2*cc*(aa*aa+dd*dd))/(cc*cc-4)**2
 assert s.simplify(s.diff(penalty_expr,cc)-deriv)==0
 gram=s.Matrix([[2,-cc],[-cc,2]]);v=s.Matrix([-aa,-dd])
 assert s.simplify((v.T*gram.inv()*v)[0]+penalty_expr)==0
 paths=[BASE,GAP,PAIR,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
        'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py',f'scripts/{TAG}.py']
 return dict(status='PASS_ALL_LOW_POLAR11_SOURCE_ANCHORS_EXCLUDED',
  sources={p:sha(p)for p in paths},interval_bits=BITS,
  normalized_old_polar4_negative_product=C4.enc(),exact_zero_pairs=zero_pairs,
  old_polar_right_neighbors=polar_right,source_anchor_vertices=[0,32,33],
  universal_reach_bounds=dict(right_pair=right_reach.enc(),own_ball=own_reach.enc(),
    high_ball_diameter=high_double_radius.enc(),wall4_if_contact_at_least4=four_reach.enc()),
  required_source_distance_cosh=DA.enc(),high_contact_square_lower=L2.enc(),
  strict_margins=dict(right_pair= (DA-right_reach).enc(),own_ball=(DA-own_reach).enc(),
    high_ball_diameter=(DA-high_double_radius).enc(),wall4_contact_bound=(DA-four_reach).enc(),
    old_polar4_lower=(C4-F(44,5)).enc(),old_polar4_upper=(9-C4).enc(),
    contact_square_above6=(L2-6).enc(),polar_contact_below13_5=(F(169,25)-L2).enc()),
  rational_lorentz_penalty_upper=str(penalty),rational_total_norm_lower=str(lower),
  norm2_contradiction_margin=str(lower-2),symbolic_identity_count=3,
  conclusion='For every all-dotted cap direction of the frozen mixed12, no raw coherent perpendicular attachment of this fixed H8 source can have low normalized projected center11 (kappa11²<sqrt5−1). Allthree exceptional sourcecorners and allcontact scales are covered.',
  scope='Inherits the exact source-gap and at-most-one-low-center theorems. Coverage concerns only the twelve old-wall balls, so added nonperpendicular caps, changed old walls and other source polytopes are outside scope. No-low directions remain open here. No numerical optimizer or sampled cap direction is a proof dependency.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run()
 p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['rational_total_norm_lower'],'>2',flush=True)
