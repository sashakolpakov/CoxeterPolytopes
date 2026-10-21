#!/usr/bin/env python3
"""Exact one-factor cusp link families and actual-ridge obstructions.

All predicates use standard-library Fraction arithmetic. The report supplies
continuous parameter arguments; no coefficient or angle-label grid is used.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from collections import Counter
import argparse,hashlib,json,math
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_single_factor_links_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';WALK='results/fresh_h21_source_corner_walk_20261013.json';OLD='results/fresh_h20_cusp_sections_20261013.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mv(A,x):return [dot(r,x)for r in A]
def pair(a,b,Q):return dot(a,mv(Q,b))
def inv(A):
 n=len(A);M=[[F(x)for x in r]+[F(i==j)for j in range(n)]for i,r in enumerate(A)]
 for j in range(n):
  i=next(i for i in range(j,n)if M[i][j]);M[i],M[j]=M[j],M[i];q=M[j][j];M[j]=[x/q for x in M[j]]
  for i in range(n):
   if i!=j:q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[j],strict=True)]
 return [r[n:]for r in M]
def pd(A):
 A=[[F(x)for x in r]for r in A];out=[]
 while A:
  q=A[0][0];assert q>0;out.append(q)
  A=[[A[i][j]-A[i][0]*A[0][j]/q for j in range(1,len(A))]for i in range(1,len(A))]
 return out
def groups(A):
 left=set(range(len(A)));out=[]
 while left:
  todo=[min(left)];left.remove(todo[0])
  for i in todo:
   nn=[j for j in sorted(left)if A[i][j]];todo+=nn;left.difference_update(nn)
  out.append(sorted(todo))
 return out
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def proportional(a,b):
 i=next(i for i,x in enumerate(a)if x)
 return bool(b[i])and all(x*b[i]==y*a[i]for x,y in zip(a,b,strict=True))
def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']];W=[[int(x)for x in r]for r in raw['wall_coordinates_as_columns_transposed']]
 C=[mv(Q,w)for w in W];G=[[dot(w,c)for c in C]for w in W];assert len(W)==210
 h=[sum(w[i]for w in W)for i in range(22)];g=math.gcd(*h);h=[x//g for x in h]
 assert pair(h,h,Q)==-14 and all(dot(h,c)<0 for c in C)
 assert all(G[i][j]<=0 for i,j in combinations(range(210),2))
 short=[i for i in range(210)if G[i][i]==2];assert len(short)==42
 walk=json.loads((ROOT/WALK).read_text());assert walk['status']=='PASS_COMPLETE_CORNER_ORBIT_CENSUS'and not walk['pending_orbit_ids']and walk['sources'][RAW]==sha(RAW)
 cusps=[r for r in walk['corner_orbits']if r['type']=='ideal'];assert len(cusps)==18
 old=json.loads((ROOT/OLD).read_text());assert old['distinct_normal_planes']==59
 faces={}
 def face(ids):
  ids=tuple(sorted(set(ids)))
  if ids in faces:return faces[ids]
  A=[[G[i][j]for j in ids]for i in ids];piv=pd(A);B=inv(A);assert all(x>=0 for r in B for x in r)
  cc=mv(B,[dot(h,C[i])for i in ids]);p=[F(h[d])-sum(c*W[i][d]for c,i in zip(cc,ids,strict=True))for d in range(22)]
  ps=[dot(p,c)for c in C];q=pair(p,p,Q);future=pair(p,h,Q)
  assert q<0 and future<0 and [i for i,x in enumerate(ps)if not x]==list(ids)and max(ps)==0
  out={'walls':list(ids),'point':p,'norm':q,'future_pair':future,'maximum_other_pairing':max(x for x in ps if x),'LDL':piv};faces[ids]=out;return out
 census=[];sumplanes=[];rank2=[];distinctness=0
 for cusp in cusps:
  delta=list(map(F,cusp['primitive_ray']));active=cusp['active_walls'];assert pair(delta,delta,Q)==0 and pair(delta,h,Q)<0
  values=[dot(delta,c)for c in C];assert max(values)==0 and [i for i,x in enumerate(values)if not x]==active
  for comp in cusp['affine_components']:
   ids=comp['walls'];marks=list(map(F,comp['positive_kernel']));mu=F(comp['ray_multiple']);A=[[G[i][j]for j in ids]for i in ids];n=comp['rank']
   assert n==len(ids)-1 and min(marks)>0 and mu>0 and all(dot(r,marks)==0 for r in A)
   assert [sum(v*W[i][d]for v,i in zip(marks,ids,strict=True))for d in range(22)]==[mu*x for x in delta]
   for j in range(len(ids)):pd([[A[a][b]for b in range(len(ids))if b!=j]for a in range(len(ids))if a!=j])
   census.append({'cusp':cusp['id'],'walls':ids,'rank':n,'marks':marks,'null_multiple':mu})
   if n==2:rank2.append((cusp,ids,delta))
   for k in[2,3]:
    for selected in combinations(ids,k):
     if not all(G[i][i]==2 and sum(G[i][j]!=0 for j in ids if j!=i)==1 for i in selected):continue
     centers=[next(j for j in ids if j!=i and G[i][j])for i in selected]
     if len(set(centers))!=1:continue
     center=centers[0]
     if not all(G[i][j]==G[selected[0]][j]for i in selected for j in ids if j not in selected):continue
     if not all(G[i][j]==0 for i,j in combinations(selected,2)):continue
     mark=marks[ids.index(selected[0])];assert all(marks[ids.index(i)]==mark for i in selected)
     Delta=[mu/mark*x for x in delta];e0=[sum(W[i][d]for i in selected)for d in range(22)];q=pair(e0,e0,Q);assert q==2*k
     ep=[dot(e0,c)for c in C];assert all(ep[i]==0 for i in active if i not in ids)
     weights={i:marks[ids.index(i)]/mark for i in ids};rest=[i for i in ids if i not in selected]
     assert Delta==[e0[d]+sum(weights[i]*W[i][d]for i in rest)for d in range(22)]
     H=[[F(G[i][j])-F(ep[i]*ep[j],q)for j in ids]for i in ids]
     zero=[i for i in ids if H[ids.index(i)][ids.index(i)]==0];keep=[i for i in ids if i not in zero]
     assert all(all(v==0 for v in H[ids.index(i)])for i in zero)
     B=[[H[ids.index(i)][ids.index(j)]for j in keep]for i in keep]
     assert all(B[i][j]<=0 and F(B[i][j]**2,B[i][i]*B[j][j])in[0,F(1,4),F(1,2),F(3,4),1]for i,j in combinations(range(len(keep)),2))
     blocks=[];rank=0
     for positions in groups(B):
      ii=[keep[x]for x in positions];ww=[F(1)if i in selected else weights[i]for i in ii];M=[[B[a][b]for b in positions]for a in positions]
      assert all(dot(r,ww)==0 for r in M)and min(ww)>0;pd([r[:-1]for r in M[:-1]])
      const=[sum(v*(F(W[i][d])-F(ep[i],q)*e0[d])for v,i in zip(ww,ii,strict=True))for d in range(22)]
      slope=-sum(v*F(ep[i],q)for v,i in zip(ww,ii,strict=True))
      expected=(F(0),F(-1))if set(ii)==set(selected)else(F(1),F(1))
      assert slope==expected[1]and const==[expected[0]*x for x in Delta]
      blocks.append({'walls':ii,'positive_weights':ww,'rank':len(ii)-1,'null_multiple_constant_slope':list(expected)});rank+=len(ii)-1
     assert rank==n-1 and len(blocks)==(1 if zero else 2)
     for i in zero:
      assert [F(W[i][d])-F(ep[i],q)*e0[d]for d in range(22)]==Delta and -F(ep[i],q)==1
     # Auxiliary short wall attaches to exactly one selected leaf and no rest.
     auxiliary=next(j for j in short if sum(G[j][i]==-1 for i in selected)==1 and all(G[j][i]in[0,-1]for i in selected)and all(G[j][i]==0 for i in rest))
     attached=next(i for i in selected if G[auxiliary][i]==-1);retained=next(i for i in selected if i!=attached)
     assert G[auxiliary][retained]==0 and ep[auxiliary]==-1 and dot(Delta,C[auxiliary])==-1
     left=face(rest+[retained,auxiliary]);right=face(list(selected)+[auxiliary]);L=left['point'];R=right['point']
     a=-dot(e0,mv(Q,L));b=-dot(Delta,mv(Q,R));assert a>0 and b>0
     assert pair(Delta,L,Q)==-a and pair(e0,R,Q)==0
     X=[b*x for x in L];Y=[a*x for x in R]
     assert pair(e0,X,Q)==pair(Delta,Y,Q)==-a*b and pair(e0,Y,Q)==0 and pair(Delta,X,Q)==-a*b
     assert sorted(set(left['walls'])&set(right['walls']))==sorted([retained,auxiliary])and pair(X,Y,Q)<0
     assert F(4)-F(5,k)>0 # projected ridge determinant for 0<1+t<1.
     # Each plane is new relative to all59 old twin-difference planes.
     new_support=[i for i in active if ep[i]];assert len(new_support)==k+1
     for oldrow in old['pencils']:
      od=list(map(F,oldrow['cusp']))
      if not proportional(Delta,od):continue
      oe=list(map(F,oldrow['base_normal']));ops=[i for i in active if dot(oe,C[i])]
      assert pair(od,od,Q)==0 and pair(oe,od,Q)==0 and pair(oe,oe,Q)>0
      assert len(ops)==2 and len(new_support)>2;distinctness+=1
     sumplanes.append({'cusp':cusp['id'],'factor':ids,'selected_leaves':list(selected),'center':center,'normal_constant':e0,'normal_shift':Delta,'normal_norm':q,'domain':'-1<t<0','local_blocks':blocks,'zero_quotient_walls':zero,'retained_local_facets':len(keep),'auxiliary':auxiliary,'attached_leaf':attached,'ridge_walls':[retained,auxiliary],'left_face':left,'right_face':right,'ridge_point_at_minus_one':X,'ridge_point_at_zero':Y,'ridge_projected_pairing_constant_slope':[F(1,k),F(1,k)],'ridge_gram_determinant_in_u_equals_1_plus_t':[F(4)-F(4,k),F(0),-F(1,k)]})
 c2=[];a2=[]
 for cusp,ids,delta in rank2:
  norms=[G[i][i]for i in ids]
  if sorted(norms)==[2,2,4]:
   a,b=[i for i in ids if G[i][i]==2];c=next(i for i in ids if G[i][i]==4)
   assert [[G[i][j]for j in[a,b,c]]for i in[a,b,c]]==[[2,0,-2],[0,2,-2],[-2,-2,4]]
   D=[W[a][d]+W[b][d]+W[c][d]for d in range(22)];assert proportional(D,delta)
   pairs=[(j,k)for j in short if G[j][a]==-1 and G[j][b]==G[j][c]==0 for k in short if G[k][b]==-1 and G[k][a]==G[k][c]==G[k][j]==0]
   assert pairs;j,k=pairs[0]
   triangle=[]
   for omitted in[a,b,c]:
    z=face([i for i in[a,b,c]if i!=omitted]+[j,k]);p=z['point'];scale=-pair(p,D,Q);assert scale>0
    p=[x/scale for x in p];assert pair(p,D,Q)==-1 and pair(p,W[omitted],Q)==-1
    triangle.append({'only_nonzero_factor_wall':omitted,'face':z,'normalized_point':p})
   assert set.intersection(*(set(r['face']['walls'])for r in triangle))=={j,k}
   same_sign=[]
   for large,small,aux in[(a,b,j),(b,a,k)]:
    l=face([small,c,aux]);r=face([a,b,aux]);lp=l['point'];rp=r['point']
    lp=[x/-pair(lp,D,Q)for x in lp];rp=[x/-pair(rp,D,Q)for x in rp]
    assert pair(lp,W[large],Q)==-1 and pair(rp,W[c],Q)==-1
    assert set(l['walls'])&set(r['walls'])=={small,aux}
    same_sign.append({'large_coefficient_wall':large,'positive_smaller_wall':small,'auxiliary':aux,'left_face':l,'right_face':r,'normalized_left':lp,'normalized_right':rp})
   c2.append({'cusp':cusp['id'],'walls':[a,b,c],'null_sum':D,'auxiliary_pair':[j,k],'triangle_face_points':triangle,'same_sign_faces':same_sign,'normal':'e=x*a+y*b+t*(a+b+c), (x,y)!=0','interior_domain':'-max(x,y,0)<t<-min(x,y,0)','normal_norm':'2*(x*x+y*y)','auxiliary_ridge_pairing':'-(x+t)*(y+t)/(2*(x*x+y*y))','same_sign_ridge_pairing':'2*y*(x+t)/(2*(x*x+y*y)), after x>=y>0'})
  else:
   assert norms==[4,4,4]and all(G[i][j]==-2 for i,j in combinations(ids,2))
   D=[sum(W[i][d]for i in ids)for d in range(22)];assert proportional(D,delta)
   for negative in ids:
    positives=[i for i in ids if i!=negative]
    pairs=[(a,b)for b in short if G[b][negative]==-2 and all(G[b][i]==0 for i in positives)for a in short if G[a][b]==-1 and all(G[a][i]==0 for i in ids)]
    assert pairs;a,b=pairs[0];pl=face([a,b]+positives);mi=face([a,b,negative])
    assert set(pl['walls'])&set(mi['walls'])=={a,b}
    assert pair(pl['point'],W[negative],Q)<0 and all(pair(mi['point'],W[i],Q)<0 for i in positives)
    a2.append({'cusp':cusp['id'],'factor':ids,'positive_walls':positives,'negative_wall':negative,'ridge_walls':[a,b],'positive_face':pl,'negative_face':mi,'parameter_domain':'x,y>=0, x+y>0, t>0','normal':'x*ri+y*rj-t*rk','normal_norm':'q=4*(x*x+y*y+t*t-x*y+x*t+y*t)','projected_ridge_squared_cosine':'q/(4*q-8*t*t)','strict_gap':'q-4*t*t=4*((x-y)^2+x*y+t*(x+y))>0'})
 assert len(census)==53 and len(sumplanes)==79 and len(rank2)==4 and len(c2)==2 and len(a2)==6
 assert Counter(len(r['selected_leaves'])for r in sumplanes)=={2:59,3:20}
 assert sum(bool(r['zero_quotient_walls'])for r in sumplanes)==2
 mutual=0
 for a,b in combinations(sumplanes,2):
  da,db=a['normal_shift'],b['normal_shift']
  if not proportional(da,db):continue
  incident=[i for i in range(210)if dot(da,C[i])==0]
  sa=[i for i in incident if dot(a['normal_constant'],C[i])]
  sb=[i for i in incident if dot(b['normal_constant'],C[i])]
  assert sa!=sb;mutual+=1
 return enc({'status':'PASS_SINGLE_FACTOR_SUM_PLANES_AND_ALL_RANK_TWO_SPANS','sources':{p:sha(p)for p in[RAW,WALK,OLD,f'scripts/{TAG}.py']},'factor_census':census,'counts':{'cusps':18,'factors':53,'sum_planes':79,'two_leaf_sum_planes':59,'three_leaf_sum_planes':20,'locally_two_factor_links':77,'zero_quotient_bound_cases':2,'old_plane_distinctness_checks':distinctness,'mutual_same_null_ray_distinctness_checks':mutual,'C2_full_spans':2,'A2_full_spans':2,'A2_negative_coefficient_cases':6,'exact_ambient_faces':len(faces)},'sum_planes':sumplanes,'full_C2_spans':c2,'full_A2_spans':a2,'scope':'All real shifts in79 explicit positive leaf-sum planes are rejected by actual obtuse ridges. All real interior section normals in every rank-two affine factor span in the18-cusp census are rejected. No classification or exclusion of remaining higher-rank single-factor directions is claimed.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');aa=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if aa.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['counts'])
