#!/usr/bin/env python3
"""Exact complete incidence and unbounded links at two new H22 null rays.

Reads the frozen lower placement, enumerates every signed wall family with
checked native integer bounds, then uses independent Fraction/SymPy arithmetic
for the local Gram components and a four-dimensional positive recession space.
The added coordinate in local geometry denotes a multiple of sqrt(11).
"""
import argparse,hashlib,json
from fractions import Fraction as F
from functools import reduce
from math import gcd,lcm
from pathlib import Path
import numpy as np
import sympy as sp
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_new_cusp_completion_20261013'
LOWER='results/'+PRE+'_lower.json';COMPLETE='results/'+PRE+'_verify.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
FRAMES='data/enumerations/'+PRE+'_frames.npz';HIGHEST='data/enumerations/'+PRE+'_highest_roots.npz'
UFILE='data/enumerations/fresh_h22_symmetry_orbit_20261013.npz'
KFILE='data/enumerations/fresh_h22_cusp_orbits_20261013.npz'
RFILE='data/enumerations/fresh_h22_gap_face_orbit_20261013.npz'
OUT='results/fresh_h22_lower_new_cusp_20261013.json'
DEN=26040

def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),F(0))
def mv(M,x):return[dot(row,x)for row in M]
def frac(x):return F(int(x.p),int(x.q))
def primitive(v):
 a=list(map(frac,v));d=lcm(*(x.denominator for x in a));n=[int(d*x)for x in a];g=reduce(gcd,n);n=[x//g for x in n]
 if next(x for x in n if x)<0:n=[-x for x in n]
 return n

def run():
 lower=read(LOWER);complete=read(COMPLETE);assert lower['source']==COMPLETE and lower['source_sha256']==sha(COMPLETE)
 assert all(sha(p)==h for p,h in complete['sources'].items())
 raw=read(RAW)['coxeter_input'];Q=np.array([[-int(x)for x in row]for row in raw['basis_gram']],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
 h=[int(x)for x in W.sum(0)];gg=reduce(gcd,h);h=np.array([x//gg for x in h],dtype=np.int64)
 assert int(h@Q@h)==-14
 fr=np.load(P/FRAMES);D=fr['six_cusp_old'];V=fr['first_tangent_numerator'];B=fr['second_tangent_numerator']
 hi=np.load(P/HIGHEST);A=hi['old_bases'];cusp_ids=hi['cusp_indices'];sizes=hi['component_sizes']
 U=np.load(P/UFILE)['twice_old_projections'];ks=np.load(P/KFILE);K=ks['old_null_projections'];four=ks['four_heights'];R=np.load(P/RFILE)['six_old_projections']
 assert D.shape==V.shape==B.shape==(120960,22) and A.shape==(846720,22) and R.shape==(241920,22)
 assert 212+2*(len(U)+len(K)+len(R))+22*len(D)==lower['total_facets']==5397548
 for a in(Q,W,D,V,B,A,U,K,R):assert a.dtype==np.int64
 # A uniform prebound covers all unscaled geometric contractions below.
 maxraw=max(int(abs(a).max())for a in(Q,W,D,V,B,A,U,K,R));assert 22**2*maxraw**3<2**63
 assert (DEN+200*(DEN//120))*maxraw<2**63
 C=np.array(lower['new_cap_old_numerators'],dtype=np.int64);shifts=list(map(F,lower['shifts']));assert shifts==list(map(F,['3','3','4','4','1','1','2','10','5','1/20','1/20']))
 assert lower['common_denominator']==DEN
 Qf=[[F(int(x))for x in row]+[F(0)]for row in Q]+[[F(0)]*22+[F(11)]];hf=[F(int(x))for x in h]+[F(0)];hcov=mv(Qf,hf)
 outputs=[]
 for typ,original in [(4,32),(5,63)]:
  # The representative is -old_original+delta_old; its extra parallel wall is26.
  c=[F(int(x),DEN)for x in C[typ]]+[2*shifts[typ]];old26=[F(int(x))for x in W[26]]+[F(0)];oldroot=[F(int(x))for x in W[original]]+[F(0)]
  assert c[:-1]==[F(int(x),6)-r for x,r in zip(D[0],oldroot[:-1],strict=True)] and c[-1]==2
  cc=mv(Qf,c);q=dot(c,cc);assert q==2 and dot(old26,cc)==-2 and dot(oldroot,cc)==-2
  parallels=[i for i,w in enumerate(W)if dot([F(int(x))for x in w]+[F(0)],cc)==-2 and dot([F(int(x))for x in w]+[F(0)],mv(Qf,[F(int(x))for x in w]+[F(0)]))==2]
  assert parallels==[26,original]
  delta=[a+b for a,b in zip(c,old26,strict=True)];dc=mv(Qf,delta);assert dot(delta,dc)==0 and dot(delta,hcov)<0 and delta[-1]==2
  six=[6*x for x in delta[:-1]];assert all(x.denominator==1 for x in six);v=np.array([int(x)for x in six],dtype=np.int64)
  cov=Q.astype(object)@v.astype(object);assert all(abs(int(x))<2**63 for x in cov);cov=np.array(cov,dtype=np.int64)
  assert 22*maxraw*int(abs(cov).max())+792<2**63
  oldpair=W@cov;assert max(oldpair)==0;oldids=np.flatnonzero(oldpair==0).tolist()
  # Timelike old family: U/2 +/-sqrt(5/2), pairing(UQv)/12 +/-sqrt110.
  up=U@cov;um=int(max(up));assert um<0 and um*um>15840
  # Old cusp caps: (K,+/-four/4), pairing(KQv)/6 +/-four*sqrt11/2.
  kp=K@cov;nullmax={}
  for f in(1,4):
   m=int(max(kp[four==f]));assert m<0 and m*m>99*f*f;nullmax[str(f)]=m
  # k6 family: R/6 +/-sqrt11; common denominator36.
  rp=R@cov;rplus=rp+792;rminus=rp-792;assert max(rplus)<0 and max(rminus)<0
  # Both center caps, all m>=3: |height|<=4.
  hp=int(h@cov);assert hp<0 and hp*hp>25344
  newids=[];new_bounds={}
  families=[('highest_roots',A,cusp_ids,None),('first+',V,np.arange(len(D)),200),('first-',-V,np.arange(len(D)),100),('second+',B,np.arange(len(D)),1),('second-',-B,np.arange(len(D)),1)]
  for name,base,ix,const in families:
   maxima=[None,None];zeros=[[],[]]
   for start in range(0,len(base),8192):
    stop=min(start+8192,len(base));aa=base[start:stop];ii=ix[start:stop]
    if name=='highest_roots':
     norms=np.einsum('ij,ij->i',aa,aa@Q);ss=sizes[start:stop]
     assert set(np.unique(norms))<={2,4}
     t=np.where(ss>=4,60,np.where(ss==3,80,np.where(norms==2,20,40)))
     bb=aa*DEN
    elif name.startswith('first'):t=np.full(len(aa),const,dtype=np.int64);bb=aa*(DEN//28)
    else:t=np.full(len(aa),const,dtype=np.int64);bb=aa*(DEN//5208)
    roots=bb+t[:,None]*(DEN//120)*D[ii]
    bound=22*int(abs(roots).max())*int(abs(cov).max())+66*int(max(t))*(DEN//5);assert bound<2**63
    old=roots@cov;extra=66*t*(DEN//5)
    for sign,pr in enumerate((old+extra,old-extra)):
     mx=int(max(pr));assert mx<=0;maxima[sign]=mx if maxima[sign]is None else max(mx,maxima[sign])
     for off in np.flatnonzero(pr==0):
      entry=start+int(off);zeros[sign].append(entry);newids.append({'family':name,'entry':entry,'cusp':int(ix[entry]),'height_sign':1 if sign==0 else-1})
   new_bounds[name]={'max_pair_numerators':maxima,'pair_denominator':6*DEN,'zero_entries':zeros,'signed_walls_tested':2*len(base)}
  assert newids==[{'family':'highest_roots','entry':typ,'cusp':0,'height_sign':1}]
  # All incident inequalities have now been enumerated across the entire chamber.
  incident=[[F(int(x))for x in W[i]]+[F(0)]for i in oldids]+[c];labels=['old'+str(i)for i in oldids]+['new_type'+str(typ)]
  covs=[mv(Qf,r)for r in incident];G=[[dot(r,s)for s in covs]for r in incident]
  assert len(incident)==18 and all(G[i][j]<=0 for i in range(18)for j in range(i))
  assert sp.Matrix(G).rank()==17 and sp.Matrix(incident).rank()==18
  unseen=set(range(18));components=[];finite_types=[];affine_count=0
  while unseen:
   cset={min(unseen)};queue=list(cset);unseen-=cset
   for i in queue:
    extra={j for j in unseen if G[i][j]!=0};unseen-=extra;cset|=extra;queue.extend(extra)
   js=sorted(cset);M=sp.Matrix([[G[i][j]for j in js]for i in js]);n=len(js)
   if M.det()==0:
    assert M==sp.Matrix([[2,-2],[-2,2]]) and len(M.nullspace())==1 and M.nullspace()[0]==sp.ones(2,1)
    assert [incident[js[0]][k]+incident[js[1]][k]for k in range(23)]==delta
    kind='affine A1';affine_count+=1;kernel=['1','1']
   else:
    assert all(M[i,i]==2 for i in range(n)) and all(M[i,j]in(-1,0)for i in range(n)for j in range(i))
    degrees=[sum(M[i,j]==-1 for j in range(n))for i in range(n)]
    assert (degrees==[0]if n==1 else sorted(degrees)==[1,1]+[2]*(n-2))
    assert all(M[:j,:j].det()>0 for j in range(1,n+1));kind='finite A'+str(n);finite_types.append(n);kernel=[]
   components.append({'walls':[labels[j]for j in js],'type':kind,'rank':int(M.rank()),'gram':[[str(x)for x in row]for row in M.tolist()],'positive_kernel':kernel})
  assert sorted(finite_types)==[1,3,3,4,5] and affine_count==1
  # These are genuine lineality directions on a horosphere, not only a rank count.
  ns=sp.Matrix(covs+[hcov]).nullspace();assert len(ns)==4;basis=[primitive(x)for x in ns]
  for w in basis:
   assert all(dot(w,z)==0 for z in covs+[hcov,dc])and dot(w,mv(Qf,w))>0
  bg=[[dot(a,mv(Qf,b))for b in basis]for a in basis];assert all(sp.Matrix(bg)[:j,:j].det()>0 for j in range(1,5))
  centerpairs=[dot(hf,row)for row in covs];assert all(x<0 for x in centerpairs)
  outputs.append({'representative_type':typ,'original_parallel_old_wall':original,'extra_parallel_old_wall':26,'null_ray':list(map(str,delta)),'null_ray_center_pairing':str(dot(delta,hcov)),'incident_old_wall_indices':oldids,'incident_new_orbit_walls':newids,'inherited_strict_bounds':{'twice_timelike_max_numerator':um,'null_caps_max_numerators_by_four_height':nullmax,'k6_signed_max_numerators':[int(max(rplus)),int(max(rminus))],'center_old_numerator':hp},'new_signed_family_incidence':new_bounds,'incident_facets':18,'normal_span_rank':18,'quotient_gram_rank':17,'horosphere_dimension':21,'components':components,'lineality_dimension':4,'positive_lineality_basis':basis,'lineality_gram':[[str(x)for x in row]for row in bg],'strict_center_incident_pairings':list(map(str,centerpairs))})
 return {'status':'PASS_EXACT_TWO_NEW_UNBOUNDED_H22_CUSP_LINKS','sources':{n:sha(n)for n in(LOWER,COMPLETE,RAW,FRAMES,HIGHEST,UFILE,KFILE,RFILE)},'chamber_facets':lower['total_facets'],'local_coordinate_field':'last coordinate is the coefficient of sqrt11; form diag(Q,11)','representatives':outputs,'scope':'Complete signed-wall incidence for these two explicit null rays in the lower placement, uniformly for all old center labels m>=3. Each link has finite A5+A4+2A3+A1 components plus affine A1 and a four-dimensional lineality space, so is unbounded. This is not a classification of all new cusps or a new completion search.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(P/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print([(x['representative_type'],x['incident_facets'],x['quotient_gram_rank'],x['lineality_dimension'])for x in d['representatives']])
