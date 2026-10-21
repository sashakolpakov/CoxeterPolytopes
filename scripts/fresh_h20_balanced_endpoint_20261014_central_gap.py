#!/usr/bin/env python3
"""Uniform exact gap for central plus every balanced endpoint cap."""
import argparse,json,hashlib
from fractions import Fraction as F
from pathlib import Path
import fresh_h20_balanced_endpoint_20261014 as E
ROOT=E.ROOT;TAG='fresh_h20_balanced_endpoint_20261014_central_gap'
SRC='results/fresh_h20_balanced_endpoint_20261014_bounds.json'
sub=lambda a,b:E.add(a,E.scale(b,-1))
def dec(x):return tuple(map(F,x))
def run():
 d=json.loads((ROOT/SRC).read_text())
 for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 assert d['uniform_radius_upper_bound']==14 and d['uniform_incenter_ball_power_lower_bound']==68
 assert F(264)<F(65,4)**2 and F(271)<F(33,2)**2
 records=[]
 for row in d['scales']:
  unit=row['name']=='unit';pmin=68 if unit else 75
  for p in row['profiles']:
   assert E.sign(sub(dec(p['all_radius_upper_bound']),E.c(14)))<=0
   assert E.sign(sub(dec(p['incenter_ball_power_lower_bound']),E.c(pmin)))>=0
  assert F(pmin-4,56)>1
  lam=F(3,20)if unit else F(1,8);N=dec(row['all_vertex_square']);q=E.scale(N,lam*lam)
  if unit:
   assert q==E.c(F(441,100))
   endpoint_power=F(68)+F(441,100)-F(21,5)*F(65,4)
   assert endpoint_power==F(104,25)
  else:
   assert E.sign(sub(q,E.c(5)))>0 and E.sign(sub(E.c(F(529,100)),q))>0
   endpoint_power=F(75)+5-2*F(23,10)*F(33,2)
   assert endpoint_power==F(41,10)
  assert endpoint_power>F(1,4)and E.sign(sub(q,E.c(F(17,4))))>0
  # Monotonicity in the unknown ball power P follows from P>=68>|x|^2.
  assert E.sign(sub(E.c(pmin),q))>0
  D=E.scale(sub(q,E.c(F(1,4))),F(1,2))
  assert sub(q,E.scale(D,2))==E.c(F(1,4))
  central_pair=sub(E.c(2),D);future_pair=E.add(E.c(-128),E.scale(D,-1))
  assert E.sign(central_pair)<0 and E.sign(future_pair)<0
  offsets=[dec(x)for x in row['factor_scales']];gaps=[]
  for n,scale in zip([6,6,7],offsets):
   # z0 uses the affine-zero omission of each factor.
   fg=[E.scale(scale,1-lam+lam*(n+1))]+[E.scale(scale,1-lam)]*n
   assert all(E.sign(x)>0 for x in fg);gaps+=fg
  assert len(gaps)==22
  records.append(dict(scale=row['name'],spatial_point='lambda times vertex omitting0,7,14',
     vertex_multiplier=lam,spatial_square=q,D=D,U=1,Lorentz_square='1/4',
     all22_seed_gaps=gaps,central_pairing=central_pair,strict_center_pairing=future_pair,
     uniform_endpoint_ball_power_strict_lower_bound=endpoint_power,
     normalized_negative_central_endpoint_pairing_lower_bound=F(pmin-4,56)))
 paths=[SRC,'scripts/fresh_h20_balanced_endpoint_20261014_bounds.py',
        'scripts/fresh_h20_balanced_endpoint_20261014.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return E.enc(dict(status='PASS_UNIFORM_CENTRAL_PLUS_ALL_ENDPOINT_CAP_STRICT_ESCAPE',
  sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},records=records,
  proof='For every endpoint ball, P=|w|^2-R^2>=Pmin and R<=14. '
        'Cauchy gives power(x)>=P+|x|^2-2|x|sqrt(P+R^2). '
        'This increases withP sinceP>=Pmin>|x|^2, and decreases withR. '
        'UsePmin=68 or75 andR=14, then rational square-root upper bounds.',
  scope='No finite-volume completion retaining the22 seed facets and radius-two central facet using only caps ending an original cusp edge, at unit or(1,1,phi)scales. The bound ignores mutual cap compatibility and is therefore stronger than any such subset selection. Interior-face caps and changed seed walls remain outside the theorem.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
