#!/usr/bin/env python3
"""All-label rigidity of side lifts of the literal compact H8 source."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_side_tilts_20261014'
SOURCE='data/known_examples/bugaenko_compact_h8.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 src=json.loads((ROOT/SOURCE).read_text());S5=s.sqrt(5);phi=(1+S5)/2
 f=lambda x:s.Rational(x['a'])+s.Rational(x['b'])*phi
 G=s.Matrix([[f(x)for x in row]for row in src['wall_gram']])
 N=s.Matrix([[f(x)for x in row]for row in src['simple_roots']])
 Q=s.Matrix([[f(x)for x in row]for row in src['lattice_gram']])
 assert G.shape==(11,11)and all(s.simplify(x)==0 for x in N*Q*N.T-G)
 assert all(G[i,i]==2 for i in range(11))
 edges=[[i,j]for i in range(11)for j in range(i+1,11)if G[i,j]==-1]
 unseen=set(range(11));components=[]
 while unseen:
  start=min(unseen);unseen.remove(start);c={start};todo=[start]
  while todo:
   x=todo.pop()
   for i,j in edges:
    y=j if i==x else i if j==x else None
    if y in unseen:unseen.remove(y);c.add(y);todo.append(y)
  components.append(sorted(c))
 assert components==[[0],[1,2,3,4,5,6,7,9,10],[8]]
 assert G[1,3]==0 and G[0,1]==-phi and G[8,7]==-phi
 # A label3 source edge with one right lift has ambient cosine sin(theta)/2,
 # strictly between0 and1/2 for every finite nonright theta=pi/m.
 # Thus right/nonright status is constant on each component above. The
 # central component cannot be nonright because source-right pair1,3 would
 # acquire positive ambient pairing2cos(theta1)cos(theta3).
 c5sq=phi*phi/4
 gaps={
  'm3_above_pi3':s.simplify(3*c5sq/4-s.Rational(1,4)),
  'm3_below_pi4':s.simplify(s.Rational(1,2)-3*c5sq/4),
  'm4_above_pi3':s.simplify(c5sq/2-s.Rational(1,4)),
  'm4_below_pi4':s.simplify(s.Rational(1,2)-c5sq/2),
  'm_at_least5_below_pi3':s.simplify(s.Rational(1,4)-c5sq*(1-c5sq)),
 }
 assert all(v>0 for v in gaps.values())
 # Formal normalized lift identity, with angle cosine c and sine sqrt(1-c²).
 ci,cj,u=s.symbols('ci cj u',real=True)
 ambient=s.sqrt(1-ci*ci)*s.sqrt(1-cj*cj)*u+2*ci*cj
 assert s.simplify(ambient.subs({cj:0,u:-1})+s.sqrt(1-ci*ci))==0
 assert s.simplify(ambient.subs(u,0)-2*ci*cj)==0
 paths=[SOURCE,f'scripts/{TAG}.py']
 return dict(status='PASS_ALL_LABEL_H8_SOURCE_SIDE_LIFTS_FORCED_RIGHT',sources={p:sha(p)for p in paths},
  source_dimension=8,source_facets=11,pi3_edges=edges,pi3_components=components,
  central_component_right_pair=[1,3],remaining_pi5_edges=[[0,1],[8,7]],
  exact_squared_cosine_gap_margins={k:str(v)for k,v in gaps.items()},
  conclusion='Every actual ambient Coxeter sidewall inducing one of the11facets of this literal H8 source must meet its containing facet at pi/2. All finite labels are covered.',
  scope='Rigid side-lift theorem for this exact compact H8 facet geometry. It does not prohibit extra ambient walls disjoint from the facet, other H8sources, or changing the face geometry; it does not by itself prove nonexistence of a higher-dimensional completion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
