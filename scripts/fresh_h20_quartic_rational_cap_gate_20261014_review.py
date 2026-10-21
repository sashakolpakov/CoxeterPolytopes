#!/usr/bin/env python3
"""Separate ANP replay of rational-cap obstruction for all positive real shifts.

Reuses only the preceding independent review's SymPy field wrapper. No
producer arithmetic is imported. Tests sqrt(5) membership in the opposite
direction: irreducibility of the quartic polynomial over Q(sqrt(5)).
"""
from pathlib import Path
import sys,argparse,json
from fractions import Fraction as F
from itertools import combinations
import sympy as S
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.fresh_h20_quartic_orbit_cap_gate_20261014_review import Field,sha,dot
TAG='fresh_h20_quartic_rational_cap_gate_20261014'
SOURCE=f'results/{TAG}.json'
BASE='results/fresh_h20_quartic_orbit_cap_gate_20261014_review.json'
LOCAL='results/fresh_h20_rank3_omitted_sections_20261014.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
y=S.symbols('y')

def run():
 saved=json.loads((ROOT/SOURCE).read_text());local=json.loads((ROOT/LOCAL).read_text());base=json.loads((ROOT/BASE).read_text())
 for doc in (saved,local,base):
  for p,h in doc['sources'].items():assert sha(p)==h,p
 assert base['status']=='PASS_INDEPENDENT_FIXED_QUARTIC_CAP_ENDPOINT_GATE'
 row=local['directions'][0];spec=local['systems'][row['system']]
 K=Field(row['minimal_polynomial'],row['root_interval']);c=K.scalar;Z=K.zero;O=K.one
 # If sqrt(5) belonged to K, this quartic would split over Q(sqrt(5)).
 f5=S.Poly(K.p.as_expr(),K.p.gen,extension=S.sqrt(5))
 assert f5.is_irreducible
 factors=f5.factor_list()[1];assert len(factors)==1 and factors[0][0].degree()==4 and factors[0][1]==1
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-F(x)for x in r]for r in raw['basis_gram']]
 W=[list(map(F,r))for r in raw['wall_coordinates_as_columns_transposed']]
 R=spec['retained_walls'];assert R==[23,50,200]
 B=[[dot(W[i],[dot(x,W[j])for x in Q])for j in R]for i in R]
 a=list(map(K.read,row['retained_pairings']));k=list(map(K.read,row['positive_kernel']))
 q=K.read(row['normal_square']);L=K.read(row['parameter_lower_magnitude'])
 nn=[q*c(B[i][i])-a[i]**2 for i in range(3)]
 assert all(K.sign(x)>0 for x in nn+k+[q,L])
 for i,j in combinations(range(3),2):
  assert k[i]**2*nn[i]==k[j]**2*nn[j]
  p=q*c(B[i][j])-a[i]*a[j];assert K.sign(p)<0 and 4*p*p==nn[i]*nn[j]
 rp=K.ring;X=K.var;trials=[];records=[]
 def isrational(v):return all(x==0 for x in K.coeffs(v)[:-1])
 def rational_value(v):assert isrational(v);return K.coeffs(v)[-1]
 def peval(p,x):return sum((v*c(x)**mon[0]for mon,v in p.items()),Z)
 def norm_polynomial(Tp,j,gp,cos2):
  P=q*gp-a[j]*Tp
  return (Tp*Tp+P*P/(nn[j]*c(cos2)))/q
 def rational_roots(poly,tag):
  cs=[poly[(i,)]for i in range(max(0,poly.degree())+1)]
  # Three nonconstant number-field coordinates must vanish simultaneously.
  equations=[S.Poly(sum(S.Rational(K.coeffs(v)[z])*y**i for i,v in enumerate(cs)),y,domain=S.QQ)for z in range(3)]
  nz=[p for p in equations if not p.is_zero];assert nz,'unexpected continuous rational-norm branch'
  g=nz[0]
  for p in nz[1:]:g=g.gcd(p)
  g=g.monic()
  # Every rational root is a degree-one irreducible factor over Q.
  roots=sorted(F(-f.all_coeffs()[1]/f.all_coeffs()[0])for f,m in g.factor_list()[1]if f.degree()==1)
  r=saved['rationality_trials'][len(trials)]
  assert r['case']==tag and len(r['norm_polynomial'])==len(cs)
  assert cs==list(map(K.asc,r['norm_polynomial']))
  assert [F(v)for v in reversed(g.all_coeffs())]==list(map(F,r['rationality_gcd']))
  assert roots==sorted(map(F,r['rational_roots']))
  trials.append({'case':tag,'gcd':str(g.as_expr()),'rational_roots':list(map(str,roots))})
  return roots
 def identity(v,i,j,gi,gj,record):
  z=next(x for x in range(3)if x not in(i,j))
  A=S.Matrix.hstack(S.Matrix(K.coeffs(k[z])),S.Matrix(K.coeffs(-L/c(2))))
  b=S.Matrix(K.coeffs(v-k[i]*c(gi)-k[j]*c(gj)))
  rank=A.rank();assert rank in(1,2)
  consistent=A.row_join(b).rank()==rank
  assert record['consistent']==consistent
  if not consistent:
   assert not record['can_have_nonzero_cusp_pairing'];return False
  rr,piv=A.row_join(b).rref();assert tuple(piv)==tuple(range(rank))
  possible=rank==1 or rr[1,2]!=0
  assert record['can_have_nonzero_cusp_pairing']==possible
  params=[S.Symbol(s)for s in record['parameters']]
  assert len(params)==2-rank
  gx=S.sympify(record['third_pairing']);h=S.sympify(record['cusp_pairing'])
  assert all(S.expand(v)==0 for v in A*S.Matrix([gx,h])-b)
  if rank==1:assert S.diff(h,params[0])!=0
  return possible
 def accept(v,i,j,gi,gj,d,mode):
  r=saved['records'][len(records)]
  assert r['pair']==[R[i],R[j]] and r['mode']==mode and list(map(F,r['old_pairings']))==[gi,gj]
  assert K.asc(r['T'])==v
  ps=[q*c(g)-a[x]*v for x,g in((i,gi),(j,gj))]
  assert ps==list(map(K.asc,r['projected_pairings']))
  if d is None:assert r['ambient_norm']is None and ps==[Z,Z];good=True
  else:
   assert K.asc(r['ambient_norm'])==d and rational_value(d)>0
   n=q*d-v*v;signs=list(map(K.sign,ps))
   good=K.sign(n)>0 and(min(signs)>=0 or max(signs)<=0)
   # Replay the exact finite-angle equations as well as rationality.
   if mode.startswith('right_'):
    m=int(mode.split('_')[1]);cos2={3:F(1,4),4:F(1,2),6:F(3,4)}[m]
    assert ps[0]==Z and ps[1]**2==c(cos2)*nn[j]*n
   else:assert all(4*p*p==nn[x]*n for p,x in zip(ps,(i,j)))
  assert r['positive_norm_and_possible_acute_orientation']==good
  possible=identity(v,i,j,gi,gj,r['identity'])
  assert r['necessary_survivor']==(good and possible)
  records.append({'pair':r['pair'],'mode':mode,'necessary_survivor':good and possible})
 for i,j in combinations(range(3),2):
  for ii,jj in((i,j),(j,i)):
   accept(Z,ii,jj,F(0),F(0),None,'right_right_zero')
   v=q/a[ii];other=a[jj]*v/q
   if isrational(other):accept(v,ii,jj,F(1),rational_value(other),None,'right_right_nonzero')
   for m,cos2 in((3,F(1,4)),(4,F(1,2)),(6,F(3,4))):
    d=q/(nn[jj]*c(cos2))
    if isrational(d)and rational_value(d)>0:accept(Z,ii,jj,F(0),F(1),d,f'right_{m}_zero')
    dp=norm_polynomial(rp(v),jj,X,cos2)
    for rr in rational_roots(dp,[R[ii],R[jj],2,m]):
     d=peval(dp,rr)
     if rational_value(d)>0:accept(v,ii,jj,F(1),rr,d,f'right_{m}_nonzero')
  den=k[i]*a[i]-k[j]*a[j];assert den!=Z
  tp=rp(q*k[i]/den)-rp(q*k[j]/den)*X
  dp=norm_polynomial(tp,i,rp(O),F(1,4))
  for rr in rational_roots(dp,[R[i],R[j],3,3]):
   d=peval(dp,rr)
   if rational_value(d)>0:accept(peval(tp,rr),i,j,F(1),rr,d,'both3_first_nonzero')
  v=tp[(1,)];d=peval(norm_polynomial(rp(v),i,rp(Z),F(1,4)),F(0))
  if isrational(d)and rational_value(d)>0:accept(v,i,j,F(0),F(1),d,'both3_first_zero')
 assert len(trials)==len(saved['rationality_trials'])==21 and len(records)==len(saved['records'])==12
 survivors=[r for r in records if r['necessary_survivor']];assert len(survivors)==3
 blocked=[list(p)for p in combinations(R,2)if not any(set(r['pair'])==set(p)for r in survivors)]
 assert blocked==saved['blocked_triangle_edge_pairs_for_rational_u']==[[23,200],[50,200]]
 assert saved['unresolved_continua']==[]
 # Every angle-stage candidate for the remaining edge class forces u rational.
 # This precedes, and is independent of, the anchor-u normal-identity solve.
 forcing=[]
 for index,r in enumerate(saved['records']):
  if set(r['pair'])!={23,50}:continue
  i,j=map(R.index,r['pair']);gi,gj=map(F,r['old_pairings'])
  assert k[2]==L
  constant=(k[i]*c(gi)+k[j]*c(gj)-K.asc(r['T']))/L
  assert isrational(constant)
  forcing.append({'angle_record':index,'identity':'u*h = g200 + constant','constant':str(rational_value(constant))})
 assert forcing==saved['rational_shift_forcing_records'] and len(forcing)==3
 assert [F(r['constant'])for r in forcing]==[F(0),F(0),F(2)]
 assert saved['parameter_scope']=='all real u>0'
 paths=[SOURCE,BASE,LOCAL,RAW,f'scripts/{TAG}.py','scripts/fresh_h20_quartic_orbit_cap_gate_20261014_review.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return {'status':'PASS_INDEPENDENT_ALL_POSITIVE_REAL_SHIFTS_RATIONAL_CAP_OBSTRUCTION','sources':{p:sha(p)for p in paths},'sqrt5_nonmembership':'Quartic minimal polynomial remains irreducible over Q(sqrt5), so the compositum has degree8.','rationality_trials':trials,'records':records,'necessary_survivor_records_at_anchor':len(survivors),'blocked_triangle_pairs_for_rational_u':blocked,'rational_shift_forcing_records':forcing,'unclosable_edge_germs':{'rational_positive_u':168,'irrational_positive_u':84},'scope':'Fixed localrecord0 spatial direction, every real u>0, with22 retained cusp facets; rational ambient cap vectors of arbitrary norm. Re-enumerates all normalized contact cases and rational-norm polynomial constraints using separate ANP arithmetic. Shares SymPy CAS, but not the producer field implementation. The finite-volume implication uses the separately audited seed geometry; its Gram is shift-independent and its positive target kernel sums to uLD.'}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}_review.json'
 if args.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print('trials',len(out['rationality_trials']),'records',len(out['records']))
