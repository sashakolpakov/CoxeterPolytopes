#!/usr/bin/env python3
"""Independent small arithmetic review; not a second 48,216-case power census."""
import argparse,json,hashlib
from pathlib import Path
from fractions import Fraction as F
import sympy as s
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_balanced_half_20261014_review';SRC='results/fresh_h20_balanced_half_20261014.json';PY='scripts/fresh_h20_balanced_half_20261014.py'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def val(row):return sum(s.Rational(c)*s.sqrt(2 if i&1 else 1)*s.sqrt(3 if i&2 else 1)*s.sqrt(5 if i&4 else 1)for i,c in enumerate(row))
def run():
 d=json.loads((ROOT/SRC).read_text());assert all(sha(p)==h for p,h in d['sources'].items());T=s.Symbol('T');checks=[]
 rows=d['all_right_forbidden_finite_cosines']+d['cross_A6_A7_forbidden_finite_cosines']
 rows += [r['witness']for p in d['radius_bounds']for r in p['forbidden_exceptional_radii']]
 for row in rows:
  cc=val(row.get('other_cosine',row.get('cosine')));poly=s.Poly(s.minpoly(2*cc,T),T).monic();expected=s.Poly(sum(s.Rational(a)*T**i for i,a in enumerate(row['minpoly_2cos'])),T)
  assert poly==expected and any(v.q!=1 for v in poly.all_coeffs());checks.append(str(poly.as_expr()))
 assert s.Poly(T-s.Rational(7,4),T)==s.Poly(sum(s.Rational(a)*T**i for i,a in enumerate(d['all_right_boundary_minpoly'])),T)
 assert val(d['radial_probe_evaluations'][0]['minimum_power_lower_bound'])==s.Rational(164,25)
 # Independently reproduce the all-right own-vertex minimum at relaxed R=8.
 assert F(196)*F(3,5)**2-64==F(164,25)
 r=16*(s.sqrt(2)-1);assert s.simplify(r*r-42).is_positive and s.simplify(44-r*r).is_positive
 for k in[12,15,16]:
  def power(R):return 49-s.Rational(k,2)*R+(s.Rational(k,16)-1)*R*R
  assert power(s.Rational(49,8))>=0
  assert s.simplify(power(r)).is_negative
 assert s.simplify((1+s.sqrt(5))/4-s.Rational(79,98)).is_positive
 assert F(42,44)-1+F(7,32)>0 and F(35,84)<F(1,2)
 x2=F(784,25);D=(x2-F(1,4))/2;assert D==F(3111,200)and x2-2*D==F(1,4)
 assert F(9,25)>F(1,4)and F(164,25)>F(1,4)
 return dict(status='PASS_SOURCE_PROOF_AND_SMALL_INDEPENDENT_ARITHMETIC_REVIEW',sources={p:sha(p)for p in[SRC,PY,f'scripts/{TAG}.py']},independent_nonintegral_minpolynomials=len(checks)+1,all_right_own_vertex_relaxed_power='164/25',union_escape_norm='1/4',scope='Independent SymPy minimal-polynomial checks and Fraction radial constants plus source/proof review. Not a second endpoint-profile enumeration or 48,216-value power replay.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/'results'/f'{TAG}.json';d=run()
 if a.verify:assert json.loads(p.read_text())==d
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
