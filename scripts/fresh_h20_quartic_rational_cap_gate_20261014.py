#!/usr/bin/env python3
"""Necessary endpoint gates for rational ambient caps, one full quartic pencil.

Uses exact quartic arithmetic and rational polynomial gcds. This is not a
construction or an enumeration of the full reflection orbit.
"""
from pathlib import Path
import sys,json,argparse
from itertools import combinations
import sympy as sp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.fresh_h20_quartic_orbit_cap_gate_20261014 import KField,F,dot,enc,sha
TAG='fresh_h20_quartic_rational_cap_gate_20261014'
LOCAL='results/fresh_h20_rank3_omitted_sections_20261014.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
SEED='results/fresh_h20_quartic_orbit_cap_gate_20261014.json'


def run():
 packet=json.loads((ROOT/LOCAL).read_text());row=packet['directions'][0];system=packet['systems'][row['system']]
 for p,h in packet['sources'].items():assert sha(p)==h
 seed=json.loads((ROOT/SEED).read_text())
 for p,h in seed['sources'].items():assert sha(p)==h
 assert seed['local_record']==0 and seed['u']=='1/2'
 assert len(seed['seed_cusp']['seed_walls'])==22 and seed['seed_cusp']['projected_span_rank']==20
 K=KField(row['minimal_polynomial'],row['root_interval']);Z,O=K.zero,K.one
 t,r,x=sp.symbols('t r x')
 f=sp.Poly(sum(sp.Rational(c)*t**i for i,c in enumerate(reversed(row['minimal_polynomial']))),t)
 SK=sp.QQ.alg_field_from_poly(f)
 sf=sp.Poly(x*x-5,x,domain=SK).factor_list()[1]
 assert len(sf)==1 and sf[0][0].degree()==2 and sf[0][1]==1
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in rr]for rr in raw['basis_gram']];W=[list(map(int,rr))for rr in raw['wall_coordinates_as_columns_transposed']]
 R=system['retained_walls'];assert R==[23,50,200] and system['target']=='A2'
 B=[[dot(W[i],[dot(xx,W[j])for xx in Q])for j in R]for i in R]
 a=[K.read(p)for p in row['retained_pairings']];k=[K.read(p)for p in row['positive_kernel']];q=K.read(row['normal_square']);L=K.read(row['parameter_lower_magnitude'])
 nn=[K.sub(K.scale(q,B[i][i]),K.mul(a[i],a[i]))for i in range(3)]
 assert all(K.sign(v)>0 for v in nn+k+[q,L])
 for i,j in combinations(range(3),2):
  assert K.mul(K.mul(k[i],k[i]),nn[i])==K.mul(K.mul(k[j],k[j]),nn[j])
  pp=K.sub(K.scale(q,B[i][j]),K.mul(a[i],a[j]))
  assert K.sign(pp)<0 and K.scale(K.mul(pp,pp),4)==K.mul(nn[i],nn[j])
 records=[];continuum=[];trials=[]
 def polyexpr(p,coord):return sp.Poly(sum(sp.Rational(c[coord])*r**z for z,c in enumerate(p)),r,domain=sp.QQ)
 def rational_norm_roots(p,tag):
  pp=[polyexpr(p,z)for z in range(1,K.d)];nonzero=[v for v in pp if not v.is_zero]
  if not nonzero:
   continuum.append({'case':tag,'norm_polynomial':p});return []
  g=nonzero[0]
  for v in nonzero[1:]:g=sp.gcd(g,v)
  rr=[F(str(v))for v in sp.polys.polytools.ground_roots(g)]
  trials.append({'case':tag,'norm_polynomial':p,'rationality_gcd':list(reversed([str(v)for v in g.monic().all_coeffs()])),'rational_roots':rr})
  return rr
 def identity(T,i,j,gi,gj):
  z=next(v for v in range(3)if v not in[i,j]);v=k[z];w=K.scale(L,F(-1,2));b=K.sub(T,K.add(K.scale(k[i],gi),K.scale(k[j],gj)))
  A=sp.Matrix([[sp.Rational(v[z]),sp.Rational(w[z])]for z in range(K.d)]);bb=sp.Matrix(list(map(sp.Rational,b)))
  try:sol,params=A.gauss_jordan_solve(bb)
  except ValueError:return {'consistent':False,'can_have_nonzero_cusp_pairing':False}
  h=sp.expand(sol[1]);ok=h!=0
  return {'consistent':True,'third_pairing':str(sol[0]),'cusp_pairing':str(h),'parameters':[str(v)for v in params],'can_have_nonzero_cusp_pairing':ok}
 def accept(T,i,j,gi,gj,d,mode):
  ident=identity(T,i,j,gi,gj)
  if d is None:
   # Both pair contacts are right; a positive projected norm is left free.
   assert K.sub(K.scale(q,gi),K.mul(a[i],T))==Z and K.sub(K.scale(q,gj),K.mul(a[j],T))==Z
   good=True;ps=[Z,Z]
  else:
   assert all(v==0 for v in d[1:]);assert d[0]>0
   n=K.sub(K.mul(q,d),K.mul(T,T));good=K.sign(n)>0
   ps=[K.sub(K.scale(q,g),K.mul(a[z],T))for z,g in[(i,gi),(j,gj)]]
   signs=[K.sign(v)for v in ps];good=good and (min(signs)>=0 or max(signs)<=0)
  records.append({'pair':[R[i],R[j]],'mode':mode,'old_pairings':[gi,gj],'T':T,'ambient_norm':d,'projected_pairings':ps,'positive_norm_and_possible_acute_orientation':good,'identity':ident,'necessary_survivor':good and ident['can_have_nonzero_cusp_pairing']})
 def dpoly(Tp,i,gip,c):
  ap=K.psub(K.pfscale(gip,q),K.pfscale(Tp,a[i]))
  return K.padd(K.pfscale(K.pmul(Tp,Tp),K.inverse(q)),K.pfscale(K.pmul(ap,ap),K.inverse(K.scale(K.mul(nn[i],q),c))))
 for i,j in combinations(range(3),2):
  for ii,jj in[(i,j),(j,i)]:
   # Both right: g_ii=0 -> T=0 -> g_jj=0; otherwise normalize g_ii=1.
   accept(Z,ii,jj,F(0),F(0),None,'right_right_zero')
   T=K.div(q,a[ii]);other=K.div(K.mul(a[jj],T),q)
   if all(v==0 for v in other[1:]):accept(T,ii,jj,F(1),other[0],None,'right_right_nonzero')
   # One right, one m=3,4,6. Label5 is unavailable since sqrt5 is not in K.
   for m,c in[(3,F(1,4)),(4,F(1,2)),(6,F(3,4))]:
    # g_ii=0, g_jj !=0 can be scaled to 1.
    dd=K.div(q,K.scale(nn[jj],c))
    if all(v==0 for v in dd[1:]) and dd[0]>0:accept(Z,ii,jj,F(0),F(1),dd,f'right_{m}_zero')
    dp=dpoly([T],jj,[Z,O],c)
    for rr in rational_norm_roots(dp,[R[ii],R[jj],2,m]):
     d=K.peval(dp,rr)
     if d[0]>0:accept(T,ii,jj,F(1),rr,d,f'right_{m}_nonzero')
  # Both labels3. The affine determinant vanishes; length ratios give T linearly.
  den=K.sub(K.mul(k[i],a[i]),K.mul(k[j],a[j]));assert den!=Z
  tp=[K.div(K.mul(q,k[i]),den),K.scale(K.div(K.mul(q,k[j]),den),-1)]
  dp=dpoly(tp,i,[O],F(1,4))
  for rr in rational_norm_roots(dp,[R[i],R[j],3,3]):
   d=K.peval(dp,rr)
   if d[0]>0:accept(K.peval(tp,rr),i,j,F(1),rr,d,'both3_first_nonzero')
  T=tp[1];d=K.peval(dpoly([T],i,[Z],F(1,4)),F(0))
  if all(v==0 for v in d[1:]) and d[0]>0:accept(T,i,j,F(0),F(1),d,'both3_first_zero')
 # All-zero pairings in the both3 branch imply T=0 and projected cap norm0.
 assert not continuum,continuum
 survivors=[v for v in records if v['necessary_survivor']]
 blocked=[list(p)for p in combinations(R,2)if not any(set(v['pair'])==set(p)for v in survivors)]
 # The entire angle stage above is independent of the real shift u. At any
 # PSD3 endpoint, a right contact or the both3 ratio forces T into K first.
 # Every angle candidate on pair23,50 forces u rational by the normal identity.
 forcing=[]
 for index,rec in enumerate(records):
  if set(rec['pair'])!={23,50}:continue
  i,j=[R.index(v)for v in rec['pair']];gi,gj=rec['old_pairings']
  z=next(v for v in range(3)if v not in[i,j]);assert R[z]==200 and k[z]==L
  residual=K.sub(K.add(K.scale(k[i],gi),K.scale(k[j],gj)),rec['T'])
  constant=K.div(residual,L);assert all(v==0 for v in constant[1:])
  forcing.append({'angle_record':index,'identity':'u*h = g200 + constant','constant':constant[0]})
 assert len(forcing)==3 and blocked==[[23,200],[50,200]]
 # For nonzero rational u, h*=2*u*h is an arbitrary nonzero rational number.
 # Substituting it into the anchor identity reuses the two blocked pair gates.
 src=[LOCAL,RAW,SEED,'scripts/fresh_h20_quartic_orbit_cap_gate_20261014.py','scripts/fresh_h20_rank3_global_sections_20261014.py',f'scripts/{TAG}.py']
 return enc({'status':'PASS_ALL_POSITIVE_REAL_SHIFTS_NO_RATIONAL_CAP_FINITE_VOLUME_COMPLETION','scope':'Localrecord0 spatial direction, every real u>0, retaining the22 cusp walls and adding projections of rational ambient normals of arbitrary positive norm. Does not cover other spatial directions, nonrational ambient caps, or altered seed walls.','sources':{p:sha(p)for p in src},'sympy_version':sp.__version__,'local_record':0,'anchor_u':'1/2','parameter_scope':'all real u>0','sqrt5_not_in_K_factor_degrees':[v.degree()for v,m in sf],'rationality_trials':trials,'records':records,'survivors_at_anchor':survivors,'blocked_triangle_edge_pairs_for_rational_u':blocked,'rational_shift_forcing_records':forcing,'unresolved_continua':continuum})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'records',len(out['records']),'anchor_survivors',len(out['survivors_at_anchor']),'blocked_for_rational_u',out['blocked_triangle_edge_pairs_for_rational_u'])
