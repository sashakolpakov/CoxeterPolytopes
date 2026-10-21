#!/usr/bin/env python3
"""Exact finite B2-closure roots and actual-ridge replay on stated families.

Elements are a+b*sqrt(D) over Q(sqrt(d)), with recursive exact signs.
No independence of the two radicals is assumed; no outer-field division is
used. A numerical scout supplies only candidate defect indices.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from functools import lru_cache
import argparse,gzip,hashlib,json,math
import fresh_h20_cusp_sections_20261013 as P
ROOT=P.ROOT;TAG='fresh_h20_b2_edge_sections_20261014';ZERO=(F(0),F(0));ONE=(F(1),F(0))
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def scale(a,c):return(a[0]*c,a[1]*c)
def mul(a,b,d):return(a[0]*b[0]+d*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def sgn(x):return int(x>0)-int(x<0)
def ksign(a,d):
 x,y=a
 if not y or not d:return sgn(x)
 if not x:return sgn(y)
 if sgn(x)==sgn(y):return sgn(x)
 return sgn(x)*sgn(x*x-d*y*y)
def divide(a,b,d):
 den=b[0]**2-d*b[1]**2;assert den
 return scale(mul(a,(b[0],-b[1]),d),1/den)
class Field:
 def __init__(self,d,D):self.d=d;self.D=D;assert ksign(D,d)>=0
 def k(self,a):return(a,ZERO)
 def c(self,a):return((F(a),F(0)),ZERO)
 def plus(self,a,b):return(add(a[0],b[0]),add(a[1],b[1]))
 def times(self,a,b):return(add(mul(a[0],b[0],self.d),mul(mul(a[1],b[1],self.d),self.D,self.d)),add(mul(a[0],b[1],self.d),mul(a[1],b[0],self.d)))
 def scaled(self,a,c):return(scale(a[0],c),scale(a[1],c))
 def sign(self,a):
  x,y=ksign(a[0],self.d),ksign(a[1],self.d)
  if not ksign(self.D,self.d)or not y:return x
  if not x:return y
  if x==y:return x
  det=add(mul(a[0],a[0],self.d),scale(mul(mul(a[1],a[1],self.d),self.D,self.d),-1))
  return x*ksign(det,self.d)
 def equal(self,a,b):return self.sign(self.plus(a,self.scaled(b,-1)))==0
 def dot(self,a,b):
  z=self.c(0)
  for x,y in zip(a,b,strict=True):z=self.plus(z,self.times(x,y))
  return z
 def rational_dot(self,a,b):
  z=self.c(0)
  for x,y in zip(a,b,strict=True):z=self.plus(z,self.scaled(y,x))
  return z
 def numeric(self,a):
  v=lambda x:float(x[0])+float(x[1])*math.sqrt(self.d)
  return v(a[0])+v(a[1])*math.sqrt(max(0,v(self.D)))
 def enc(self,a):return[list(map(str,a[0])),list(map(str,a[1]))]

def contacts(ga,gb,d):
 out=[]
 if gb==0:
  if d-ga*ga>0:out.append(('ORTHOGONAL',0,(F(ga),F(0))))
  dd,x={(0,2):(10,(F(0),F(1,5))),(0,4):(5,(F(0),F(2,5))),
        (-1,2):(6,(F(-4,5),F(1,5))),(-2,4):(0,(F(-6,5),F(0)))}[ga,d]
  out.append(('B3_AT_MODIFIED_END_1',dd,x))
 if gb<0 and d-ga*ga==2*gb*gb:out.append(('B3_AT_UNCHANGED_END',0,(F(ga),F(0))))
 return out

@lru_cache(None)
def cosine_bounds(m):
 """Rational enclosure for cos²(pi/m), using Machin and alternating series."""
 assert m>=2
 if m==2:return F(0),F(0)
 def atan_bounds(x):
  terms=[(-1)**k*x**(2*k+1)/F(2*k+1)for k in range(18)]
  low=sum(terms);high=low+x**37/F(37);return low,high
 a,b=atan_bounds(F(1,5));c,d=atan_bounds(F(1,239));pl,pu=16*a-4*d,16*b-4*c
 assert F(314159,100000)<pl<pu<F(314160,100000)
 def cos_bounds(t):
  terms=[(-1)**k*t**(2*k)/F(math.factorial(2*k))for k in range(18)]
  low=sum(terms);high=low+t**36/F(math.factorial(36));return low,high
 low=cos_bounds(2*pu/m)[0];high=cos_bounds(2*pl/m)[1]
 return(1+low)/2,(1+high)/2

def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return[enc(y)for y in x]
 if isinstance(x,dict):return{str(k):enc(y)for k,y in x.items()}
 return x

def run(selected):
 recipes=json.loads((ROOT/f'results/{TAG}_exact_recipes.json').read_text())
 for p,h in recipes['source_sha256'].items():assert sha(p)==h,p
 proposed={i:[r for r in recipes['candidate_witnesses']if r['family']==i]for i in selected}
 raw=json.loads((ROOT/P.RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']]
 W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']];C=[P.mv(Q,w)for w in W];G=[[P.dot(w,c)for c in C]for w in W]
 walk=json.loads((ROOT/P.WALK).read_text())['corner_orbits'];h=[sum(w[j]for w in W)for j in range(22)];gg=math.gcd(*h);h=[F(x,gg)for x in h];hc=P.mv(Q,h)
 gamma=F(1,3);families=[];rows=[];counts=Counter()
 for fid in selected:
  f=recipes['families'][str(fid)];corner=walk[f['corner']];edge=corner['edge_orbits'][f['edge_index']];E=f['edge_walls'];a=f['anchor'];b=next(i for i in f['component']if i!=a);ai=E.index(a)
  assert E==edge['edge_walls']and G[a][a]==G[b][b]==2 and G[a][b]==-1
  assert all(G[i][j]==0 for i in[a,b]for j in E if j not in[a,b])
  A=[[F(G[i][j])for j in E]for i in E];piv=P.pd(A);I=P.inv(A);assert I[ai][ai]==F(2,3)
  dual=[[sum(I[k][j]*G[l][i]for j,l in enumerate(E))for i in range(210)]for k in range(20)]
  dh=[sum(I[k][j]*P.dot(W[l],hc)for j,l in enumerate(E))for k in range(20)]
  left=list(map(F,corner['primitive_ray']));right=list(map(F,edge['endpoint_primitive_ray']));left=P.scale(-1/P.dot(left,hc),left);right=P.scale(-1/P.dot(right,hc),right)
  B=[[P.pair(x,y,Q)for y in[left,right]]for x in[left,right]];pv=[P.dot(left,c)for c in C];pb=[P.dot(right,c)for c in C]
  assert max(pv)==max(pb)==0 and all(pv[i]==pb[i]==0 for i in E)
  census=[];family_counts=Counter();matched=set()
  for l in range(210):
   if l in E:continue
   ga,gb,d=G[l][a],G[l][b],G[l][l]
   for kind,base,x in contacts(ga,gb,d):
    p=[pv[l],pb[l]];pp=P.dot(p,p);assert pp>0
    target=add(x,(-F(2*ga+gb,3),F(0)));v0=[F(y,pp)for y in p];y=[p[1],-p[0]]
    aa=P.dot(y,P.mv(B,y));bb=scale(target,2*P.dot(v0,P.mv(B,y)));dd=add(scale(mul(target,target,base),P.dot(v0,P.mv(B,v0))),(-gamma,F(0)))
    D=add(mul(bb,bb,base),scale(dd,-4*aa));record={'wall':l,'kind':kind,'base_radicand':base,'contact':x,'quadratic':[aa,bb,dd],'discriminant':D};census.append(record)
    if aa:
     if ksign(D,base)<0:record['status']='NO_REAL_ROOT';continue
     field=Field(base,D);roots=[(s,(scale(bb,-1/(2*aa)),(F(s,2)/aa,F(0))))for s in([-1,1]if ksign(D,base)>0 else[0])]
    elif ksign(bb,base):
     field=Field(base,ZERO);roots=[(0,field.k(scale(divide(dd,bb,base),-1)))]
    else:assert ksign(dd,base);record['status']='NO_ROOT';continue
    record['branches']=[]
    for branch,root in roots:
     zc=[field.plus(field.k(scale(target,v0[i])),field.scaled(root,y[i]))for i in range(2)]
     vals=[field.rational_dot(B[i],zc)for i in range(2)]
     assert field.equal(field.dot(zc,vals),field.c(gamma))and field.equal(field.rational_dot(p,zc),field.k(target))
     signs=[field.sign(v)for v in vals]
     if signs[0]*signs[1]>=0:
      status='BOUNDARY_ROOT'if 0 in signs else'OUTSIDE_TRANSVERSE_INTERIOR';record['branches'].append(status);family_counts[status]+=1;continue
     record['branches'].append('TRANSVERSE_INTERIOR');family_counts['transverse_recipes_with_multiplicity']+=1
     distances=sorted((max(abs(field.numeric(zc[i])-r['edge_plane_coefficients'][i])for i in range(2)),j)for j,r in enumerate(proposed[fid]))
     assert distances and distances[0][0]<1e-5
     r=proposed[fid][distances[0][1]];matched.add(r['index']);assert r['status']=='NUMERICAL_ACTUAL_RIDGE_DEFECT',r
     coeff=[vals[1],field.scaled(vals[0],-1)];future=field.scaled(field.plus(coeff[0],coeff[1]),-1)
     if field.sign(future)>0:coeff=[field.scaled(v,-1)for v in coeff];future=field.scaled(future,-1)
     ps=[field.plus(field.scaled(coeff[0],pv[i]),field.scaled(coeff[1],pb[i]))for i in range(210)]
     qp=field.dot(coeff,[field.rational_dot(B[i],coeff)for i in range(2)])
     assert field.sign(qp)<0 and field.sign(future)<0 and all(field.sign(ps[i])==(0 if i in E else-1)for i in range(210))
     zp=[field.plus(field.scaled(zc[0],pv[i]),field.scaled(zc[1],pb[i]))for i in range(210)]
     ep=[field.plus(field.c(dual[ai][i]),zp[i])for i in range(210)];assert all(field.equal(ep[i],field.c(int(i==a)))for i in E)
     witness=r['witness'];release=witness['release'];new=witness['new_wall'];ri=E.index(release)
     dp=[field.plus(field.c(dual[ri][i]),field.scaled(zp[i],-I[ri][ai]/gamma))for i in range(210)]
     assert field.sign(dp[new])<0
     xp=[field.plus(field.scaled(field.times(dp[new],ps[i]),-1),field.times(ps[new],dp[i]))for i in range(210)]
     qd=I[ri][ri]+I[ri][ai]**2/gamma
     qX=field.plus(field.times(field.times(dp[new],dp[new]),qp),field.scaled(field.times(ps[new],ps[new]),qd))
     dcenter=field.plus(field.c(dh[ri]),field.scaled(field.plus(zc[0],zc[1]),I[ri][ai]/gamma))
     xcenter=field.plus(field.scaled(field.times(dp[new],future),-1),field.times(ps[new],dcenter))
     ss=[field.sign(v)for v in xp];active=[i for i,v in enumerate(ss)if v==0]
     assert field.sign(qX)<0 and field.sign(xcenter)<0 and max(ss)==0
     assert active==sorted(set(E)-{release}|{new})==witness['active']
     i,j=witness['pair'];pij=field.plus(field.c(G[i][j]),field.scaled(field.times(ep[i],ep[j]),-1))
     ni=field.plus(field.c(G[i][i]),field.scaled(field.times(ep[i],ep[i]),-1));nj=field.plus(field.c(G[j][j]),field.scaled(field.times(ep[j],ep[j]),-1));prod=field.times(ni,nj);sq=field.times(pij,pij)
     det=field.plus(prod,field.scaled(sq,-1));assert field.sign(ni)>0 and field.sign(nj)>0 and field.sign(det)>0
     proof={}
     if witness['kind']=='OBTUSE':assert field.sign(pij)>0
     else:
      assert field.sign(pij)<0;m,n=witness['adjacent_labels'];assert n==m+1
      lower=cosine_bounds(m)[1];upper=cosine_bounds(n)[0]
      assert field.sign(field.plus(sq,field.scaled(prod,-lower)))>0
      assert field.sign(field.plus(field.scaled(prod,upper),field.scaled(sq,-1)))>0
      proof={'adjacent_integer_labels':[m,n],'lower_cosine_upper_bound':lower,'upper_cosine_lower_bound':upper}
     counts[witness['kind']]+=1
     rows.append({'family':fid,'closure_wall':l,'contact_type':kind,'branch':branch,
       'field':{'base_radicand':base,'outer_radicand':D},'edge_plane_coefficients':zc,
       'initial_point_edge_coefficients':coeff,'initial_point_norm':qp,'release':release,'new_wall':new,
       'contact_norm':qX,'contact_future_pairing':xcenter,'all210_wall_signs':ss,'active_walls':active,
       'actual_bad_ridge':[i,j],'projected_pairing':pij,'projected_norms':[ni,nj],
       'positive_two_by_two_determinant':det,'defect':witness['kind'],**proof})
  actual=[r for r in proposed[fid]if r['status']=='NUMERICAL_ACTUAL_RIDGE_DEFECT']
  assert matched=={r['index']for r in actual}
  families.append({'family':fid,'source_family':f,'old_edge_LDL':piv,'complete_contact_quadratics':census,
    'exact_recipe_counts':dict(family_counts),'distinct_scout_normals_exactly_rejected':len(matched)})
  print('exact family',fid,'recipes',family_counts['transverse_recipes_with_multiplicity'],'distinct',len(matched),flush=True)
 sources=[P.RAW,P.WALK,f'results/{TAG}_exact_recipes.json',f'scripts/{TAG}_exact.py',
  'scripts/fresh_h20_cusp_sections_20261013.py']
 return enc({'status':'PASS_COMPLETE_SELECTED_CONTINUOUS_B2_FAMILIES_ACTUAL_RIDGE_DEFECTS',
  'sources':{p:sha(p)for p in sources},'selected_family_ids':selected,'families':families,'counts':dict(counts),'cases':rows,
  'scope':'Exact exhaustive closure roots on the stated families, including both B3 extensions and nonsquare two-stage radicals. Actual first contacts and strict Coxeter-angle gaps are checked against all210 source inequalities. No exclusion of unstated families.'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--families',type=int,nargs='+',default=[0,34,660]);args=ap.parse_args();z=run(args.families)
 out=ROOT/f'results/{TAG}_exact.json'
 if args.verify:assert z==json.loads(out.read_text())
 else:out.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status']);print(z['counts'])
