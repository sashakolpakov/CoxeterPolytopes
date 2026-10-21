#!/usr/bin/env python3
"""Exact necessary cap gate at one quartic source section."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.fresh_h20_rank3_global_sections_20261014 import Field,F,rpdiv,rpmul,rpadd,rpscale,dot,enc,sha
from itertools import combinations,product
import json,math,argparse
from collections import Counter
TAG='fresh_h20_quartic_orbit_cap_gate_20261014'
LOCAL='results/fresh_h20_rank3_omitted_sections_20261014.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
ENGINE='scripts/fresh_h20_rank3_global_sections_20261014.py'
WALK='results/fresh_h21_source_corner_walk_20261013.json'

class KField(Field):
 def inverse(self,a):
  assert a!=self.zero
  oldr,r=self.p,list(a);oldt,t=[F(0)],[F(1)]
  while any(r):
   quot,rem=rpdiv(oldr,r);oldr,r=r,rem;oldt,t=t,rpadd(oldt,rpscale(rpmul(quot,t),-1))
  assert len(oldr)==1 and oldr[0]
  out=self.reduce(rpscale(oldt,1/oldr[0]));assert self.mul(a,out)==self.one
  return out
 def div(self,a,b):return self.mul(a,self.inverse(b))
 def pfscale(self,p,c):return self.ptrim([self.mul(x,c)for x in p])
 def pdiv(self,a,b):
  a=self.ptrim(a);b=self.ptrim(b);assert b!=[self.zero];q=[self.zero]*max(1,len(a)-len(b)+1);ib=self.inverse(b[-1])
  while a!=[self.zero]and len(a)>=len(b):
   i=len(a)-len(b);c=self.mul(a[-1],ib);q[i]=self.add(q[i],c)
   for j,x in enumerate(b):a[i+j]=self.sub(a[i+j],self.mul(c,x))
   a=self.ptrim(a)
  return self.ptrim(q),a
 def pgcd(self,a,b):
  while b!=[self.zero]:a,b=b,self.pdiv(a,b)[1]
  return self.pfscale(a,self.inverse(a[-1]))

def rank(A):
 M=[list(map(F,r))for r in A];n=len(M);m=len(M[0]);k=0
 for j in range(m):
  p=next((i for i in range(k,n)if M[i][j]),None)
  if p is None:continue
  M[k],M[p]=M[p],M[k];q=M[k][j];M[k]=[x/q for x in M[k]]
  for i in range(k+1,n):
   q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[k],strict=True)]
  k+=1
  if k==n:break
 return k

def rational_pd(A):
 M=[list(map(F,r))for r in A];p=[]
 for i in range(len(M)):
  v=M[i][i];assert v>0;p.append(v)
  for j in range(i+1,len(M)):
   for k in range(j,len(M)):M[j][k]=M[k][j]=M[j][k]-M[j][i]*M[k][i]/v
 return p

def run():
 packet=json.loads((ROOT/LOCAL).read_text());row=packet['directions'][0];system=packet['systems'][row['system']]
 for p,h in packet['sources'].items():assert sha(p)==h
 K=KField(row['minimal_polynomial'],row['root_interval']);Z=K.zero;O=K.one
 assert K.d==4 and system['target']=='A2'and system['retained_walls']==[23,50,200]
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']];W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 R=system['retained_walls'];B=[[dot(W[i],[dot(x,W[j])for x in Q])for j in R]for i in R]
 a=[K.read(p)for p in row['retained_pairings']];k=[K.read(p)for p in row['positive_kernel']];q=K.read(row['normal_square']);L=K.read(row['parameter_lower_magnitude']);u=F(1,2)
 assert all(v!=Z for v in a)and K.sign(q)>0 and K.sign(L)>0
 nn=[K.sub(K.scale(q,B[i][i]),K.mul(a[i],a[i]))for i in range(3)]
 for i,j in combinations(range(3),2):
  pp=K.sub(K.scale(q,B[i][j]),K.mul(a[i],a[j]));assert K.sign(pp)<0 and K.scale(K.mul(pp,pp),4)==K.mul(nn[i],nn[j])
 # The full22-wall seed cusp is a bounded Euclidean product. These checks
 # prove its spatial rank19 and actual span20 after projection.
 walk=json.loads((ROOT/WALK).read_text())
 for p,h in walk['sources'].items():assert sha(p)==h
 cusp=next(c for c in walk['corner_orbits']if c['id']==4);D=list(map(F,row['null_shift']));C=[[dot(v,w)for v in Q]for w in W];G=[[dot(w,c)for c in C]for w in W]
 pair=lambda a,b:dot(a,[dot(v,b)for v in Q])
 center=[sum(w[d]for w in W)for d in range(22)]
 assert pair(center,center)<0 and pair(D,center)<0 and all(pair(center,w)<0 for w in W)
 assert pair(D,D)==0 and all(pair(D,w)<=0 for w in W)
 assert [i for i,w in enumerate(W)if not pair(D,w)]==sorted(cusp['active_walls'])
 e0=[K.read(v)for v in row['normal_constant']];assert all(e0[d]==K.linear(k,[W[r][d]for r in R])for d in range(22))
 assert all(a[i]==K.linear(k,B[i])for i in range(3));checkq=Z
 for aa,kk in zip(a,k,strict=True):checkq=K.add(checkq,K.mul(aa,kk))
 assert checkq==q and all(K.sign(v)>0 for v in k)
 seed=sorted(set(cusp['active_walls'])-{system['omitted_wall']});assert len(seed)==22 and rank([W[i]for i in seed])==21
 e=[K.sub(v,K.scale(L,u*d))for v,d in zip(e0,D,strict=True)]
 untouched=[]
 for fac in cusp['affine_components']:
  ids=fac['walls'];marks=list(map(F,fac['positive_kernel']))
  if ids==system['factor']:
   assert all(K.sign(K.sub(k[j],K.scale(L,u*marks[ids.index(r)])))>0 for j,r in enumerate(R))
   continue
  mat=[[G[i][j]for j in ids]for i in ids]
  assert all(dot(rr,marks)==0 for rr in mat)and all(x>0 for x in marks)
  multiple=F(fac['ray_multiple']);assert multiple>0
  assert [sum(mm*W[r][d]for r,mm in zip(ids,marks,strict=True))for d in range(22)]==[multiple*x for x in D]
  assert all(K.linear(e0,C[r])==Z and pair(D,W[r])==0 for r in ids)
  minors=[rational_pd([[G[i][j]for j in ids if j!=omit]for i in ids if i!=omit])for omit in ids]
  untouched.append({'walls':ids,'positive_kernel':marks,'ray_multiple':multiple,'all_one_omitted_LDL':minors})
 assert sorted(len(x['walls'])for x in untouched)==[7,12]
 HH=[[K.sub(K.scale(q,B[i][j]),K.mul(a[i],a[j]))for j in range(3)]for i in range(3)]
 for i in range(3):
  s=Z
  for j in range(3):s=K.add(s,K.mul(HH[i][j],k[j]))
  assert s==Z and K.sign(HH[i][i])>0
 # Sum k_i projected roots equals u*L*D, a positive null multiple.
 assert K.sign(K.scale(L,u))>0
 records=[];quadratics=[];old_cases=[];counts={'old_PSD3_integer_cases':0,'right_contact_trials':0,'right_contact_angle_pass':0,'both_label3_common_linear':0,'both_label3_common_quadratic':0,'T_identity_consistent':0,'integral_negative_h':0}
 def linidentity(T,i,j,gi,gj):
  z=next(x for x in range(3)if x not in[i,j]);v=k[z];w=K.scale(L,-u);b=K.sub(T,K.add(K.scale(k[i],gi),K.scale(k[j],gj)))
  piv=next(((p,r,v[p]*w[r]-v[r]*w[p])for p,r in combinations(range(K.d),2)if v[p]*w[r]-v[r]*w[p]),None)
  if piv is None:
   p=next(p for p,x in enumerate(v)if x);rho=w[p]/v[p];tau=b[p]/v[p]
   if K.scale(v,tau)!=b:return None
   assert K.scale(v,rho)==w and rho
   counts['T_identity_consistent']+=1
   den=math.lcm(rho.denominator,tau.denominator);A=den;B=int(rho*den);C=int(tau*den);ok=C%math.gcd(A,B)==0
   if ok:counts['integral_negative_h']+=1
   return {'third_plus_rho_h_equals_tau':[rho,tau],'admissible':ok,'h':'negative integer in the indicated congruence class'}
  p,r,d=piv;gz=(b[p]*w[r]-b[r]*w[p])/d;h=(v[p]*b[r]-v[r]*b[p])/d
  if K.add(K.scale(v,gz),K.scale(w,h))!=b:return None
  counts['T_identity_consistent']+=1
  if gz.denominator!=1 or h.denominator!=1 or h>=0:return {'third':gz,'h':h,'admissible':False}
  counts['integral_negative_h']+=1
  return {'third':gz,'h':h,'admissible':True}
 def angle(p,D):
  if p==Z:return 2
  if K.sign(p)>0:return None
  n=K.mul(p,p)
  for m,c in[(3,F(1,4)),(4,F(1,2)),(6,F(3,4))]:
   if n==K.scale(D,c):return m
  f=K.add(K.sub(K.scale(K.mul(n,n),16),K.scale(K.mul(n,D),12)),K.mul(D,D))
  if f==Z and K.sign(K.sub(n,K.scale(D,F(1,2))))>0:return 5
  return None
 def evaluate(T,i,j,gi,gj,norm,mode):
  nl=K.sub(K.scale(q,norm),K.mul(T,T))
  if K.sign(nl)<=0:return
  ps=[K.sub(K.scale(q,g),K.mul(a[x],T))for x,g in[(i,gi),(j,gj)]]
  labels=[angle(p,K.mul(nn[x],nl))for p,x in zip(ps,[i,j])]
  if None in labels:return
  if mode=='right':
   assert 2 in labels;counts['right_contact_angle_pass']+=1
  else:
   if labels!=[3,3]:return
  # Exact PSD3 determinant and principal-minor gates, including ideal cases.
  sij=K.sub(K.scale(q,B[i][j]),K.mul(a[i],a[j]));D=K.sub(K.mul(K.mul(nn[i],nn[j]),nl),K.add(K.add(K.mul(nn[i],K.mul(ps[1],ps[1])),K.mul(nn[j],K.mul(ps[0],ps[0]))),K.mul(nl,K.mul(sij,sij))))
  D=K.add(D,K.scale(K.mul(K.mul(sij,ps[0]),ps[1]),2))
  if K.sign(D)<0:return
  sol=linidentity(T,i,j,gi,gj)
  records.append({'pair':[i,j],'old_g':[gi,gj],'ambient_norm':norm,'mode':mode,'T':T,'labels':labels,'projected_determinant':D,'identity_solution':sol})
 for norm in[2,4]:
  for i,j in combinations(range(3),2):
   di,dj,g=B[i][i],B[j][j],B[i][j];assert di.denominator==dj.denominator==1
   mi=math.isqrt(int(di*norm));mj=math.isqrt(int(dj*norm))
   for gi,gj in product(range(-mi,mi+1),range(-mj,mj+1)):
    det=di*dj*norm+2*g*gi*gj-di*gj*gj-dj*gi*gi-norm*g*g
    if di*norm<gi*gi or dj*norm<gj*gj or det<0:continue
    counts['old_PSD3_integer_cases']+=1
    for ii,gg in[(i,gi),(j,gj)]:
     counts['right_contact_trials']+=1;T=K.div(K.scale(q,gg),a[ii]);evaluate(T,i,j,gi,gj,norm,'right')
    def equation(x,gx):
     return [K.sub(K.scale(K.mul(q,q),4*gx*gx),K.scale(K.mul(nn[x],q),norm)),K.scale(K.mul(q,a[x]),-8*gx),K.add(K.scale(K.mul(a[x],a[x]),4),nn[x])]
    pa,pb=equation(i,gi),equation(j,gj);gg=K.pgcd(pa,pb)
    old_cases.append({'pair':[i,j],'ambient_norm':norm,'old_g':[gi,gj],'old_Gram_determinant':det,'both3_quadratic_gcd':gg})
    if len(gg)==2:
     counts['both_label3_common_linear']+=1;T=K.scale(gg[0],-1);evaluate(T,i,j,gi,gj,norm,'both3')
    elif len(gg)==3:
     counts['both_label3_common_quadratic']+=1;quadratics.append({'pair':[i,j],'old_g':[gi,gj],'ambient_norm':norm,'gcd':gg})
    else:assert len(gg)==1
 assert not quadratics and counts['old_PSD3_integer_cases']==120 and len(records)==15
 survivors=[];keys=set()
 for r in records:
  sol=r['identity_solution']
  if not sol or not sol['admissible']:continue
  key=(tuple(r['pair']),r['ambient_norm'],tuple(r['old_g']),tuple(r['T']))
  if key in keys:continue
  keys.add(key);survivors.append(r)
 assert len(survivors)==3 and all(r['pair']==[0,1]for r in survivors)
 for r in survivors:
  assert r['identity_solution']['third_plus_rho_h_equals_tau'][0]==-F(1,2)
  if r['mode']=='right':
   assert r['ambient_norm']in[2,4]and r['old_g']==[0,0]and r['T']==Z and r['identity_solution']['third_plus_rho_h_equals_tau'][1]==0
  else:
   assert r['ambient_norm']==4 and r['old_g']==[0,-2]and r['T']==a[2]and r['identity_solution']['third_plus_rho_h_equals_tau'][1]==4
 counts['distinct_necessary_survivor_families']=len(survivors)
 counts['seed_cusp_edge_germs']=12*3*7;counts['unclosable_edge_germs']=12*2*7
 print(counts,flush=True)
 return enc({'status':'PASS_FIXED_QUARTIC_CUSP_NO_FINITE_VOLUME_INTEGRAL_NORM2_4_COMPLETION',
  'sources':{p:sha(p)for p in[LOCAL,RAW,WALK,ENGINE,f'scripts/{TAG}.py']},'local_record':0,'u':'1/2','old_seed_Gram':B,'counts':counts,
  'seed_cusp':{'source_cusp':4,'seed_walls':seed,'null_ray':D,'unprojected_span_rank':21,'projected_span_rank':20,'projected_Gram_rank':19,'unchanged_affine_factors':untouched,'scaled_projected_A2_Gram':HH,'A2_positive_kernel':k,'A2_null_sum_multiple':K.scale(L,u),'normal_square':q},
  'old_PSD3_cases':old_cases,'angle_pass_records':records,'necessary_survivor_families':survivors,
  'unclosable_triangle_pairs':[[23,200],[50,200]],
  'scope':'Fixed local record0 at u=1/2. Retain all22 actual seed cusp walls and add only projected spacelike normals of integral ambient norm2 or4 with nonpositive pairings to all seed facets. No finite-volume polytope completion exists:168 seed cusp edge germs cannot have finite or ideal opposite endpoints. Three necessary cap families for the other84 edge germs remain; this is not a proof that every individual cap is impossible. Other section normals, arbitrary real/nonintegral mirrors, different ambient norms, or deletion of seed facets are outside scope.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
