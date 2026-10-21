#!/usr/bin/env python3
"""Separate Fraction replay of the 118 archived cusp-pencil ridge proofs."""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_cusp_sections_20261013_review'
SOURCE=ROOT/'results/fresh_h20_cusp_sections_20261013.json'
SCRIPT=ROOT/'scripts/fresh_h20_cusp_sections_20261013.py'
RAW=ROOT/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUT=ROOT/'results'/f'{TAG}.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mv(A,x):return [dot(r,x)for r in A]
def add(a,b):return [x+y for x,y in zip(a,b,strict=True)]
def scale(a,x):return [a*y for y in x]
def run():
 d=json.loads(SOURCE.read_text());raw=json.loads(RAW.read_text())['coxeter_input']
 for p,h in d['sources'].items():assert sha(ROOT/p)==h
 Q=[[-F(x)for x in row]for row in raw['basis_gram']]
 W=[list(map(F,row))for row in raw['wall_coordinates_as_columns_transposed']]
 C=[mv(Q,w)for w in W];h=list(map(F,d['strict_center']))
 def pair(x,y):return dot(x,mv(Q,y))
 assert pair(h,h)<0 and all(dot(h,c)<0 for c in C)
 count=0;link_count=0;orientations={}
 for row in d['pencils']:
  a,b,j=row['retained_leaf'],row['omitted_leaf'],row['auxiliary_wall']
  e=list(map(F,row['base_normal']));u=list(map(F,row['cusp']))
  assert e==scale(F(1,2),add(W[a],scale(-1,W[b])))
  assert pair(e,e)==1 and pair(u,u)==pair(e,u)==pair(e,h)==0 and pair(u,h)<0
  ew={int(i):F(c)for i,c in row['endpoint_positive_combination']};assert min(ew.values())>0
  assert add(e,scale(F(1,2),u))==[sum(c*W[i][k]for i,c in ew.items())for k in range(22)]
  z0,z1=[list(map(F,row[key]['point']))for key in ['left_face','right_face']]
  zero_sets=[]
  for z,key in [(z0,'left_face'),(z1,'right_face')]:
   assert pair(z,z)<0 and pair(z,h)<0
   ps=[dot(z,c)for c in C];assert max(ps)<=0
   zeros={i for i,p in enumerate(ps)if not p};assert zeros==set(row[key]['walls']);zero_sets.append(zeros)
  assert zero_sets[0]&zero_sets[1]=={a,j} and pair(z0,z1)<0
  A,B=-pair(z0,u),-pair(z1,u);assert A>0 and B>0
  assert A==F(row['interpolation_weights']['A']) and B==F(row['interpolation_weights']['B'])
  assert pair(z0,e)==0 and pair(z1,add(e,scale(F(1,2),u)))==0
  p0=scale(B/2,z0);p1=add(scale(A,z1),scale(-B,z0))
  assert pair(p0,e)==pair(p0,u)+pair(p1,e)==pair(p1,u)==0
  # Positive interpolation is future timelike for 0<t<1/2 and has
  # precisely the two asserted incident inequalities.
  assert pair(W[a],e)==1 and pair(W[a],u)==0 and pair(W[j],u)==-1
  s=1 if row['auxiliary_attached_leaf']==b else -1
  assert row['auxiliary_attached_leaf']in[a,b] and pair(W[j],e)==F(s,2)
  assert pair(W[a],W[a])==pair(W[j],W[j])==2
  assert pair(W[a],W[j])-F(s,2)==F(-1,2)
  assert list(map(F,row['cosine_squared_numerator_coefficients']))==[F(1,4),F(-1),F(1)]
  assert list(map(F,row['cosine_squared_denominator_coefficients']))==[F(7,4),F(s),F(-1)]
  gap=[F(3,4),F(4+s),F(-5)]
  assert list(map(F,row['denominator_minus_four_numerator']))==gap
  assert gap[0]>0 and sum(x*F(1,2)**i for i,x in enumerate(gap))>0
  pa=add(W[a],scale(-1,e));pb=add(W[b],e);pj=add(W[j],scale(-F(s,2),e))
  assert pa==pb and pair(pa,pj)**2/(pair(pa,pa)*pair(pj,pj))==F(1,7)
  totalrank=0
  for link in row['local_link']:
   ids=link['walls'];ks=list(map(F,link['positive_kernel']));lam,slope=map(F,link['null_multiple_coefficients'])
   assert min(ks)>0 and lam>0 and lam+slope/2>=0
   ev=[pair(W[i],e)for i in ids]
   H=[[pair(W[i],W[j])-ev[ii]*ev[jj]for jj,j in enumerate(ids)]for ii,i in enumerate(ids)]
   assert all(H[i][i]>0 for i in range(len(ids)))
   assert all(H[i][j]<=0 for i in range(len(ids))for j in range(i))
   assert all(dot(r,ks)==0 for r in H)
   reached={0}
   while True:
    new=reached|{j for j in range(len(ids))if any(H[i][j]for i in reached)}
    if new==reached:break
    reached=new
   assert len(reached)==len(ids)
   assert [sum(k*(W[i][z]-t*e[z])for k,i,t in zip(ks,ids,ev,strict=True))for z in range(22)]==scale(lam,u)
   assert scale(-sum(k*t for k,t in zip(ks,ev,strict=True)),u)==scale(slope,u)
   for i in range(len(ids)):
    for j in range(i):assert H[i][j]**2/(H[i][i]*H[j][j])in[0,F(1,4),F(1,2),F(3,4),1]
   assert link['rank']==len(ids)-1;totalrank+=len(ids)-1;link_count+=1
  assert totalrank==19
  orientations.setdefault((row['cusp_orbit'],tuple(sorted((a,b)))),set()).add((a,b))
  count+=1
 assert count==118 and len(orientations)==59 and all(len(x)==2 for x in orientations.values())
 return {'status':'PASS_INDEPENDENT_RATIONAL_RIDGE_AND_CUSP_LINK_REPLAY','sources':{str(p.relative_to(ROOT)):sha(p)for p in [SOURCE,SCRIPT,RAW,Path(__file__)]},'oriented_pencils':count,'unoriented_leaf_planes':len(orientations),'cusp_representatives':len({r['cusp_orbit']for r in d['pencils']}),'affine_link_components_checked':link_count,'arithmetic':'Separate standard-library Fraction implementation; no producer imports.','scope':'Replays all 118 archived actual-ridge interpolations, cosine exclusions, supporting endpoint identities, t=0 duplicate-leaf geometry, and full rank19 Coxeter links. Source review checks exhaustiveness of the 59 declared literal normal planes; arbitrary cusp directions or other representatives are not excluded.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run()
 if a.verify:assert json.loads(OUT.read_text())==out
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['oriented_pencils'],flush=True)
