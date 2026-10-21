#!/usr/bin/env python3
"""Exact all-parameter global gates for the rank>=4 omitted-factor census."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import scripts.fresh_h20_omitted_factor_classification_20261014 as E
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache
import argparse,json,math
ROOT=E.ROOT;TAG=E.TAG;Z=E.ZERO;O=E.ONE
add,sub,mul,scale,ff,sign=E.add,E.sub,E.mul,E.scale,E.ff,E.sign
RAW=E.RAW;LOCAL=f'results/{TAG}.json';REC=f'results/{TAG}_global_recipes.json'
DIFF='results/fresh_h20_cusp_sections_20261013.json'
ANGLE='results/fresh_h20_crystallographic_angle_gate_20261014.json'
def dec(a):return tuple(map(F,a))
def fv(v):return [dec(a)for a in v]
def invf(a):
 assert a!=Z;v=O
 for mask in range(1,8):
  v=mul(v,tuple((-x if (i&mask).bit_count()%2 else x)for i,x in enumerate(a)))
 norm=mul(a,v);assert norm[0]and not any(norm[1:]);return scale(v,1/norm[0])
def trim(p):
 p=list(p)
 while len(p)>1 and p[-1]==Z:p.pop()
 return p
def padd(a,b):return trim([add(a[i]if i<len(a)else Z,b[i]if i<len(b)else Z)for i in range(max(len(a),len(b)))])
def pscale(a,c):return trim([scale(x,c)for x in a])
def psub(a,b):return padd(a,pscale(b,-1))
def pmul(a,b):
 z=[Z]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i+j]=add(z[i+j],mul(x,y))
 return trim(z)
def peval(p,t):
 z=Z
 for a in reversed(p):z=add(mul(z,t),a)
 return z
def affine(p,l,r):
 q=[Z]
 for i,c in enumerate(p):
  for j in range(i+1):
   while len(q)<=j:q.append(Z)
   q[j]=add(q[j],scale(c,F(math.comb(i,j))*l**(i-j)*(r-l)**j))
 return trim(q)
def positive(p):
 """Strict on (0,1), with all endpoint zeros recorded explicitly."""
 p=trim(p);assert p!=[Z];left=right=0
 while p[0]==Z:p=p[1:];left+=1
 while peval(p,O)==Z:
  q=[p[0]]
  for x in p[1:-1]:q.append(add(x,q[-1]))
  assert p[-1]==E.neg(q[-1]);p=trim(q);right+=1
 def bern(q,lo=F(0),hi=F(1),depth=0):
  n=len(q)-1;b=[E.flinear(q[:k+1],[F(math.comb(k,i),math.comb(n,i))for i in range(k+1)])for k in range(n+1)]
  if all(sign(x)>0 for x in b):return [{'interval':[lo,hi],'Bernstein':b}]
  assert depth<12,(q,lo,hi,b)
  mid=(lo+hi)/2
  return bern(affine(q,F(0),F(1,2)),lo,mid,depth+1)+bern(affine(q,F(1,2),F(1)),mid,hi,depth+1)
 proof=bern(p)
 return {'zero_order_at_0':left,'zero_order_at_1':right,'reduced_polynomial':p,'positive_Bernstein_intervals':proof}

def run():
 packet=json.loads((ROOT/LOCAL).read_text());recipes=json.loads((ROOT/REC).read_text());prior=json.loads((ROOT/DIFF).read_text());angle=json.loads((ROOT/ANGLE).read_text())
 assert angle['status']=='PASS_COMPLETE_SECTION_FINITE_ANGLE_GATE_ALGEBRA'
 assert angle['allowed_actual_finite_ridge_labels']==[2,3,4,6]
 for doc in[packet,prior,angle]:
  for p,h in doc['sources'].items():assert E.sha(p)==h
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']];W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 C=[[E.dot(q,w)for q in Q]for w in W];G=[[E.dot(w,c)for c in C]for w in W]
 h=[sum(w[d]for w in W)for d in range(22)];gg=math.gcd(*h);h=[F(x,gg)for x in h]
 def pair(a,b):return E.dot(a,[E.dot(q,b)for q in Q])
 def fp(a,b):return E.fdot(a,[E.flinear(b,q)for q in Q])
 def mixed(a,b):return E.flinear(a,[E.dot(q,b)for q in Q])
 @lru_cache(None)
 def face(ids):
  A=[[F(G[i][j])for j in ids]for i in ids];piv=E.pd(A);I=E.inv(A);c=[E.dot(row,[E.dot(h,C[i])for i in ids])for row in I]
  x=[h[d]-sum(z*W[i][d]for z,i in zip(c,ids,strict=True))for d in range(22)];ps=[E.dot(x,c)for c in C]
  assert pair(x,x)<0 and pair(x,h)<0 and max(ps)==0 and [i for i,v in enumerate(ps)if not v]==list(ids)
  return {'walls':ids,'point':x,'norm':pair(x,x),'future':pair(x,h),'LDL':piv,'max_other_pairing':max(v for v in ps if v)}
 rows=packet['local_directions'];oldmatches=[];new=[]
 for ix,row in enumerate(rows):
  e=fv(row['unit_normal_constant']);D=list(map(F,row['null_shift']));ds=next(x for x in D if x);canonD=[x/ds for x in D];ep=[mixed(e,W[i])for i in row['factor']];matched=None
  for j,r in enumerate(prior['pencils']):
   d=list(map(F,r['cusp']));d0=next(x for x in d if x)
   if [x/d0 for x in d]!=canonD:continue
   ee=list(map(F,r['base_normal']));p=[pair(ee,W[i])for i in row['factor']]
   if not any(p):continue
   k=next(k for k,v in enumerate(p)if v)
   if all(scale(v,p[k])==scale(ep[k],q)for v,q in zip(ep,p,strict=True)):
    lam=scale(ep[k],1/p[k]);difference=[sub(x,scale(lam,y))for x,y in zip(e,ee,strict=True)]
    dindex=next(k for k,x in enumerate(D)if x)
    assert all(scale(v,D[dindex])==scale(difference[dindex],q)for v,q in zip(difference,D,strict=True))
    matched=j;break
  if matched is None:new.append(ix)
  else:oldmatches.append({'row':ix,'old_difference_plane':matched})
 assert len(oldmatches)==112 and len(new)==63 and sorted(new)==sorted(r['row']for r in recipes['families'])
 assert len({r['old_difference_plane']for r in oldmatches})==56
 global_proofs=[];endpoint_exceptions=[]
 for rec in recipes['families']:
  ix=rec['row'];row=rows[ix];e=fv(row['unit_normal_constant']);D=list(map(F,row['null_shift']));L=dec(row['parameter_lower_magnitude']);i,j=rec['ridge'];R=row['retained_walls'];omit=row['omitted_wall']
  assert omit not in[i,j]and sign(L)>0 and fp(e,e)==O and mixed(e,D)==Z
  A0=tuple(sorted(set(R+[i,j])));F0=face(A0)
  if rec['left_face_uses_only_omitted']:left=[omit]
  else:left=[k for k in row['factor']if k!=row['lower_bound_witness_wall']]
  A1=tuple(sorted(set(left+[i,j])));F1=face(A1);K=sorted(set(A0)&set(A1));assert set(A1)-set(K)=={omit}
  p0=F0['point'];p1=F1['point'];d0=-pair(p0,D);d1=-pair(p1,D);assert d0>0 and d1>0
  p0=[x/d0 for x in p0];p1=[x/d1 for x in p1]
  assert pair(p0,D)==pair(p1,D)==-1 and mixed(e,p0)==Z and mixed(e,p1)==E.neg(L)
  assert pair(p0,p1)<0
  # At an interior parameter, interpolation is timelike and has precisely
  # the common zero set K. e(u).p0=uL!=0 proves projected independence.
  E.pd([[F(G[a][b])for b in K]for a in K]);assert {i,j}<=set(K)and len(K)<=20
  aa=[mixed(e,W[k])for k in[i,j]];bb=[scale(L,-pair(D,W[k]))for k in[i,j]]
  pp=psub([ff(G[i][j])],pmul([aa[0],bb[0]],[aa[1],bb[1]]))
  ni=psub([ff(G[i][i])],pmul([aa[0],bb[0]],[aa[0],bb[0]]));nj=psub([ff(G[j][j])],pmul([aa[1],bb[1]],[aa[1],bb[1]]))
  NN=pmul(pp,pp);DD=pmul(ni,nj);tests={'first_norm':positive(ni),'second_norm':positive(nj)}
  if rec['mode']=='obtuse':tests['positive_pairing']=positive(pp);endpoint_good=sign(peval(pp,O))>0
  elif rec['mode']=='gap':
   lo,hi=map(F,rec['gap']);tests['negative_pairing']=positive(pscale(pp,-1));tests['above_lower']=positive(psub(NN,pscale(DD,lo)));tests['below_upper']=positive(psub(pscale(DD,hi),NN))
   endpoint_good=sign(peval(psub(NN,pscale(DD,lo)),O))>0 and sign(peval(psub(pscale(DD,hi),NN),O))>0
   if lo==F(1,2):assert ix==102 # Uses the separate all-angle source-section theorem.
  else:
   assert ix==110 and rec['mode']=='one_exception'
   assert NN==[Z,Z,O]and DD==[ff(4),Z,E.neg(O)]
   tests['only_possible_parameter']='u=2/sqrt(5), where cosine squared is1/4';endpoint_good=False
  # The sole extra endpoint-zero wall is redundant: on the section it is
  # a nonnegative combination of retained walls whose coefficients vanish
  # off K. Recover the relation directly in ambient coordinates.
  walk=json.loads((ROOT/E.WALK).read_text());fac=next(f for c in walk['corner_orbits']if c['id']==row['cusp']for f in c['affine_components']if f['walls']==row['factor'])
  marks=dict(zip(fac['walls'],map(F,fac['positive_kernel']),strict=True));kk=fv(row['positive_kernel']);coeff=[sub(k,scale(L,marks[r]))for r,k in zip(R,kk,strict=True)]
  assert all(sign(c)>=0 for c in coeff)and all(c==Z or r in K for r,c in zip(R,coeff,strict=True))
  e1=[sub(x,scale(L,d))for x,d in zip(e,D,strict=True)]
  assert e1==[sub(E.flinear(coeff,[W[r][d]for r in R]),scale(L,marks[omit]*W[omit][d]))for d in range(22)]
  boundary_only=all(c==Z for c in coeff)
  if not boundary_only and not endpoint_good:endpoint_exceptions.append(ix)
  global_proofs.append({'row':ix,'ridge':[i,j],'mode':rec['mode'],'angle_gap':rec['gap'],'parameter':'e=e0-u*L*D, 0<u<1',
   'face_at_0':F0,'face_at_1':F1,'normalized_point_0':p0,'normalized_point_1':p1,'common_zero_set':K,
   'point_cross_pair':pair(p0,p1),'pairing_polynomial':pp,'first_norm_polynomial':ni,'second_norm_polynomial':nj,'sign_certificates':tests,
   'u1_omitted_wall_redundancy_coefficients':coeff,'u1_is_old_facet_hyperplane':boundary_only,'u1_bad_angle_verified':endpoint_good})
 assert endpoint_exceptions==[83]
 exceptions=[]
 for rec in recipes['exceptions']:
  ix=rec['row'];row=rows[ix];u=ff(1)if rec['u']=='1'else E.scale(E.sqrtq(5),F(2,5));assert ix in[83,110]
  e0=fv(row['unit_normal_constant']);D=list(map(F,row['null_shift']));L=dec(row['parameter_lower_magnitude']);e=[sub(v,scale(mul(u,L),d))for v,d in zip(e0,D,strict=True)]
  ids=rec['old_edge_walls'];assert len(ids)==20;A=[[F(G[i][j])for j in ids]for i in ids];piv=E.pd(A);I=E.inv(A)
  ch=[E.dot(r,[E.dot(h,C[i])for i in ids])for r in I];p=[h[d]-sum(c*W[i][d]for c,i in zip(ch,ids,strict=True))for d in range(22)]
  ce=[E.flinear([mixed(e,W[i])for i in ids],r)for r in I];z=[sub(e[d],E.flinear(ce,[W[i][d]for i in ids]))for d in range(22)];qz=fp(z,z);assert sign(qz)>0
  alpha=mul(mixed(z,p),invf(qz));x=[sub(ff(p[d]),mul(alpha,z[d]))for d in range(22)]
  vals=[mixed(x,w)for w in W];assert fp(x,e)==Z and all(sign(v)<=0 for v in vals)and sorted(i for i,v in enumerate(vals)if v==Z)==sorted(ids)
  qx=fp(x,x);future=mixed(x,h);assert sign(qx)<0 and sign(future)<0
  i,j=rec['bad_pair'];assert i in ids and j in ids;ai=mixed(e,W[i]);aj=mixed(e,W[j]);pij=sub(ff(G[i][j]),mul(ai,aj));ni=sub(ff(G[i][i]),mul(ai,ai));nj=sub(ff(G[j][j]),mul(aj,aj));n=mul(pij,pij);d=mul(ni,nj)
  assert sign(pij)<0 and sign(n)>0 and sign(sub(scale(d,F(1,4)),n))>0
  exceptions.append({'row':ix,'u':u,'normal':e,'old_edge_walls':ids,'old_Gram_LDL':piv,'orthogonal_normal_square':qz,'actual_vertex':x,'norm':qx,'future_pair':future,'all210_pairings':vals,
   'bad_pair':[i,j],'pairing':pij,'squared_cosine_numerator':n,'squared_cosine_denominator':d,'gap':'0<cosine squared<1/4'})
 sources=[RAW,E.WALK,LOCAL,REC,DIFF,ANGLE,f'scripts/{TAG}.py',f'scripts/{TAG}_global.py']
 return E.enc({'status':'PASS_ALL_RANK_AT_LEAST_FOUR_ONE_OMITTED_GLOBAL_EXCLUSIONS','sources':{p:E.sha(p)for p in sources},
  'old_difference_plane_matches':oldmatches,'new_plane_global_certificates':global_proofs,'isolated_parameter_actual_vertices':exceptions,
  'counts':{'old_matched_omission_directions':112,'old_distinct_planes':56,'new_planes':63,'new_sqrt6_planes':48,'exceptional_actual_vertices':2},
  'scope':'All real one-omitted single-factor directions of rank at least four in the local census are excluded as finite-volume Coxeter H20 interior sections. Rank-three and other local facet counts are separate. One angle gap invokes the separately archived crystallographic-angle theorem.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}_global.json'
 if args.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['counts'])
