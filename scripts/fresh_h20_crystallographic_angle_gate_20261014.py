#!/usr/bin/env python3
"""Exact tables excluding all section ridge labels except2,3,4,6.

The compact-face argument is in the accompanying note. This script checks
the literal spherical-triple catalogs and new uniform pair inequalities;
the previous H3 and short-A2 theorems are bound by their source hashes.
"""
from pathlib import Path
from itertools import combinations_with_replacement
from collections import defaultdict
import argparse, hashlib, json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_crystallographic_angle_gate_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUT=ROOT/f'results/{TAG}.json'
q,C=s.symbols('q C');w=s.sqrt(5);C5=(3+w)/8
POOLS={(2,2):[-1,0],(2,4):[-2,0],(4,4):[-8,-6,-4,-2]}

def positive_open(poly,lo,hi):
 p=s.Poly(s.expand(poly),q,domain=s.QQ)
 assert p.degree()<=2 and not p.is_zero
 assert p.eval(lo)>=0 and p.eval(hi)>=0
 if p.degree()==2:
  a,b,_=p.all_coeffs()
  if a>0:
   t=-b/(2*a)
   if lo<t<hi:assert p.eval(t)>0

def table(name,types,lo,hi,nonright):
 rows=[]
 for a,b in combinations_with_replacement(types,2):
  da,ka=a;db,kb=b
  if ka==kb==0:continue
  na=da-q*ka*ka;nb=db-q*kb*kb
  positive_open(na,lo,hi);positive_open(nb,lo,hi)
  for g in POOLS[tuple(sorted([da,db]))]:
   p=g-q*ka*kb;det=s.expand(p*p-na*nb)
   if p==0:verdict='RIGHT_TO_UNCHANGED_SOURCE_ROOT';assert (ka==0)!=(kb==0)
   elif nonright(a,b,g):
    ratio=s.cancel(p*p/(na*nb));verdict='SPECIAL_FINITE_RATIO'
   else:
    positive_open(-p,lo,hi);positive_open(det,lo,hi)
    verdict='STRICTLY_DIVERGENT';ratio=None
   rows.append({'types':[list(a),list(b)],'old_pairing':g,'new_pairing':str(p),
    'norms':[str(na),str(nb)],'verdict':verdict,
    'squared_cosine':str(s.cancel(p*p/(na*nb)))if verdict=='SPECIAL_FINITE_RATIO'else None})
 return {'name':name,'q_interval':list(map(str,[lo,hi])),'rows':rows}

def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
 Q=-s.Matrix(raw['basis_gram']).applyfunc(s.Rational)
 W=s.Matrix(raw['wall_coordinates_as_columns_transposed']).applyfunc(s.Rational)
 CW=W*Q;GG=CW*W.T;G=[[int(v)for v in GG.row(i)]for i in range(210)]
 D=[G[i][i]for i in range(210)]
 assert all((2*CW[i,j]/D[i]).q==1 for i in range(210)for j in range(22))
 assert {d:D.count(d)for d in set(D)}=={2:42,4:168}
 for ds,pool in POOLS.items():
  assert sorted({G[i][j]for i in range(210)for j in range(i)if tuple(sorted([D[i],D[j]]))==ds})==pool
 one=defaultdict(set);two=defaultdict(set);triangles=0
 for a in range(210):
  for b in range(210):
   if a==b:continue
   da,db,g=D[a],D[b],G[a][b]
   if da*db<=g*g:continue
   for l in range(210):
    if l==a or l==b:continue
    dl,u,v=D[l],G[a][l],G[b][l]
    determinant=da*db*dl+2*g*u*v-da*v*v-db*u*u-dl*g*g
    if determinant<=0:continue
    if g and v==0:one[(da,db,g)].add((dl,u))
    if u and v:
     if g:triangles+=1
     elif da<=db:two[(da,db)].add((dl,u,v))
 assert triangles==0
 assert dict(one)=={(2,2,-1):{(2,0),(4,0),(2,-1),(4,-2)},
  (2,4,-2):{(2,0),(2,-1)},(4,2,-2):{(2,0),(4,-2)},(4,4,-2):{(2,0),(2,-2)}}
 assert dict(two)=={(2,2):{(2,-1,-1)},(2,4):{(4,-2,-2),(2,-1,-2)}}
 mixed=lambda a,b,g:(a[1]==0)!=(b[1]==0)and g==-1
 none=lambda a,b,g:False
 tables=[table('one_contact_long_A2_m_ge7',[(2,0),(2,-2)],s.Rational(1,3),s.Rational(3,8),mixed),
  table('one_contact_B2_short_anchor_m_ge7',[(2,0),(2,-1)],s.Integer(1),s.Rational(4,3),mixed),
  table('one_contact_B2_long_anchor_m_ge7',[(2,0),(4,-2)],s.Rational(1,2),s.Rational(2,3),none),
  table('two_contact_short_orthogonal_pair',[(2,0),(4,0),(2,-1)],s.Integer(1),s.Rational(6,5),mixed),
  table('two_contact_mixed_orthogonal_pair_ratio1',[(2,0),(4,-2)],s.Rational(3,4),s.Integer(1),none),
  table('two_contact_mixed_orthogonal_pair_ratio2',[(2,0),(2,-1)],s.Rational(3,2),s.Integer(2),none)]
 # Long A2 exceptional pair gives exactly the previously certified m=n=8 equation.
 qlong=C/(4*C-1);N=(4*C-1)/(8*C-4)
 assert s.cancel(1/(4-8*qlong)-N)==0
 assert s.cancel((2*C-1)*(2*N-1)-s.Rational(1,2))==0
 # Short-anchor B2 ratios cannot be Coxeter for any integer m>=7.
 f=(2*C-1)/(6*C-4)
 assert s.factor(s.diff(f,C))==-2/(3*C-2)**2/4
 assert f.subs(C,s.Rational(13,16))==s.Rational(5,7)
 assert f.subs(C,s.Rational(4,5))==s.Rational(3,4)
 assert f.subs(C,s.Rational(16,19))==s.Rational(13,20)
 assert s.sign(s.Rational(5,7)-C5)==1 and s.sign(C5-s.Rational(13,20))==1
 C8=(2+s.sqrt(2))/4;assert s.sign(C8-s.Rational(16,19))==1
 # Orthogonal short-short: q=(1+sec(pi/m))/2 <=sqrt5/2<6/5.
 assert s.sign(s.Rational(6,5)-w/2)==1 and s.sign(C5-s.Rational(5,8))==1
 # At m=5 B2 modified boundary norms are negative; long A2 finite
 # nonright pair candidates are divergent, exactly in Q(sqrt5).
 qBshort=s.simplify(C5/(2*C5-1));qBlong=s.simplify(C5/(4*C5-2));qAlong=s.simplify(C5/(4*C5-1))
 assert s.sign(2-qBshort)==-1 and s.sign(4-4*qBlong)==-1
 m5=[]
 for a,b in combinations_with_replacement([(2,0),(2,-2)],2):
  da,ka=a;db,kb=b
  if ka==kb==0:continue
  for g in[-1,0]:
   na=da-qAlong*ka*ka;nb=db-qAlong*kb*kb;p=g-qAlong*ka*kb
   assert s.sign(na)==s.sign(nb)==1
   if p==0:verdict='RIGHT_TO_UNCHANGED_SOURCE_ROOT'
   else:assert s.sign(s.simplify(p*p-na*nb))==1;verdict='STRICTLY_DIVERGENT'
   m5.append({'types':[list(a),list(b)],'old_pairing':g,'verdict':verdict})
 # Bind exact premises, including the complete integer (m,n)=(8,8) gate.
 premises=['results/fresh_h20_coupled_h3_obstruction_20261014.json',
  'results/fresh_h20_dihedral_edge_sections_20261014.json',
  'results/fresh_h20_dihedral_edge_sections_20261014_m5.json',
  'results/fresh_h20_h4_face_closure_20261014.json']
 for name in premises:
  z=json.loads((ROOT/name).read_text())
  for k in['sources','source_sha256']:
   if isinstance(z.get(k),dict):
    for p,h in z[k].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
 files=[RAW,f'scripts/{TAG}.py']+premises
 return {'status':'PASS_COMPLETE_SECTION_FINITE_ANGLE_GATE_ALGEBRA',
  'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in files},
  'one_contact_boundary_catalog':{str(k):[list(v)for v in sorted(vals)]for k,vals in sorted(one.items())},
  'two_contact_orthogonal_modified_boundary_catalog':{str(k):[list(v)for v in sorted(vals)]for k,vals in sorted(two.items())},
  'positive_definite_source_triangles':triangles,'uniform_pair_tables':tables,
  'long_A2_m5_table':m5,'B2_short_anchor_ratio':str(f),
  'allowed_actual_finite_ridge_labels':[2,3,4,6],
  'scope':'Any finite-volume Coxeter H20 hyperplane section of the literal210-wall source. Angle labels only; neither arithmeticity nor existence is asserted.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();z=run()
 if args.verify:assert z==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status']);print('tables',len(z['uniform_pair_tables']),'allowed labels',z['allowed_actual_finite_ridge_labels'])
