#!/usr/bin/env python3
"""Exact low-type classification and three mixed-height H22 chambers.

No floating point is used. Checked int64 products produce integer Gram data;
Fraction arithmetic checks the witness and the inherited Lorentz signature.
The new null-cap coefficient is a=(sqrt(6)+sqrt(3))/2. The global witness
is strict even for the incompatible union of all three 240-cap choices.
"""
import argparse, hashlib, json
from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
PRE='fresh_h22_incidence_mixed_height_20261013'
SOURCE='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
BANK='data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json'
SYM='data/enumerations/fresh_h20_global_interior_sections_20261013.json'
OUT='data/enumerations/'+PRE+'.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def mm(A,B):
 assert A.dtype==B.dtype==np.dtype('int64')
 assert A.shape[1]*int(abs(A).max())*int(abs(B).max())<2**63
 return A@B
def inertia(M):
 A=[[F(int(x)) for x in r]for r in M];out=[0,0,0]
 while A:
  k=next((i for i in range(len(A)) if A[i][i]),None)
  if k is None:
   ij=next(((i,j)for i in range(len(A))for j in range(i+1,len(A))if A[i][j]),None)
   if ij is None:out[2]+=len(A);break
   i,j=ij
   for t in range(len(A)):A[i][t]+=A[j][t]
   for t in range(len(A)):A[t][i]+=A[t][j]
   k=i
  A[0],A[k]=A[k],A[0]
  for r in A:r[0],r[k]=r[k],r[0]
  p=A[0][0];out[p<0]+=1
  A=[[A[i][j]-A[i][0]*A[0][j]/p for j in range(1,len(A))]for i in range(1,len(A))]
 return out
def gf2_solutions(inc,parity):
 rows=[sum(int(inc[i,j])<<j for j in range(21))|(parity<<21)for i in range(21)]
 piv=[];k=0
 for j in range(21):
  i=next((i for i in range(k,21)if (rows[i]>>j)&1),None)
  if i is None:continue
  rows[k],rows[i]=rows[i],rows[k]
  for i in range(21):
   if i!=k and(rows[i]>>j)&1:rows[i]^=rows[k]
  piv.append(j);k+=1
 assert all(r==0 for r in rows[k:]);free=[j for j in range(21)if j not in piv]
 out=[]
 for bits in range(1<<len(free)):
  m=sum(((bits>>i)&1)<<j for i,j in enumerate(free))
  for i,j in enumerate(piv):m|=(((rows[i]>>21)&1)^((rows[i]&m).bit_count()&1))<<j
  out.append(m)
 assert len(out)==2048
 return out
def run():
 raw=read(SOURCE)['coxeter_input'];Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
 assert inertia(Q)==[21,1,0];G=mm(mm(W,Q),W.T);wn=G.diagonal();short=np.flatnonzero(wn==2).tolist()
 assert Counter(map(int,wn))=={2:42,4:168}
 col={short[0]:0};todo=deque(col)
 while todo:
  i=todo.popleft()
  for j in short:
   if i!=j and G[i,j]==-1:
    if j not in col:col[j]=1-col[i];todo.append(j)
    assert col[j]!=col[i]
 parts=[[i for i in short if col[i]==c]for c in (0,1)];assert list(map(len,parts))==[21,21]
 P,L=W[parts[0]],W[parts[1]];inc=-mm(mm(L,Q),P.T)
 assert np.array_equal(mm(inc,inc.T),4*np.eye(21,dtype=np.int64)+np.ones((21,21),dtype=np.int64))
 SP=P.sum(axis=0);SL=L.sum(axis=0);assert np.all((5*SP+2*SL)%21==0)
 vp=(5*SP+2*SL)//21;assert np.all((SP+SL)%3==0);h=(SP+SL)//3
 assert int(h@Q@h)==-14 and int(vp@Q@vp)==-2 and max(W@Q@h)<0
 records=[];type_counts={}
 for A,weights in ((3,{7,9}),(4,{14,16})):
  for m in gf2_solutions(inc,A%2):
   u=np.array([(m>>i)&1 for i in range(21)],dtype=np.int64);w=int(u.sum())
   if w not in weights:continue
   row=A*vp-u@P;pw=mm(row[None,:],mm(Q,W.T))[0]
   good=True
   for p,n in zip(pw,wn):
    if p>0:good=False;break
    s=F(int(p*p),8*int(n)) # (row/2,w)^2 / (2*w^2)
    if s not in (0,F(1,4),F(1,2),F(3,4)) and s<1:good=False;break
   if good:records.append({'A':A,'mask':m,'weight':w,'twice_old_projection':row.tolist(),'twice_height_square':A*A-w+4})
 for typ in ((3,7),(3,9),(4,14),(4,16)):
  type_counts[str(typ)]=sum((r['A'],r['weight'])==typ for r in records)
 assert type_counts=={'(3, 7)':360,'(3, 9)':280,'(4, 14)':360,'(4, 16)':21}
 source=[r for r in read(BANK)['candidates']if (r['A'],r['weight'])in((3,7),(3,9),(4,14),(4,16))]
 assert sorted(records,key=lambda r:(r['A'],r['mask']))==sorted(source,key=lambda r:(r['A'],r['mask']))
 gens=[np.array(g['matrix'],dtype=np.int64)for g in read(SYM)['ambient_symmetry_generators']]
 for g in gens:assert np.array_equal(mm(mm(g.T,Q),g),Q)
 def orbit(seed):
  found={tuple(seed)};todo=deque(found)
  while todo:
   a=np.array(todo.popleft(),dtype=np.int64)
   for g in gens:
    b=tuple(mm(g,a[:,None])[:,0])
    if b not in found:found.add(b);todo.append(b)
  return found
 orbits={}
 for typ,n in (((3,9),280),((4,16),42),((3,7),720)):
  o=orbit(next(r['twice_old_projection']for r in records if(r['A'],r['weight'])==typ));assert len(o)==n;orbits[str(typ)]=o
 D2=np.array(sorted(orbits['(3, 9)']|orbits['(4, 16)']),dtype=np.int64);assert np.all(D2%2==0);D=D2//2
 X2=np.array(sorted(orbits['(3, 7)']),dtype=np.int64)
 assert set(map(tuple,X2))==set(tuple(r['twice_old_projection'])for r in records if r['twice_height_square']==6)
 DD=mm(mm(D,Q),D.T);XX4=mm(mm(X2,Q),X2.T);DX2=mm(mm(D,Q),X2.T)
 assert set(DD.flat)=={0,-2,-3,-4} and np.all(DD.diagonal()==0)
 assert np.all(XX4.diagonal()==-4) and set(DX2.flat)=={-8,-6,-4}
 assert set(D@Q@h)=={-6} and set(X2@Q@h)=={-14}
 dp=mm(mm(D,Q),W.T);xp=mm(mm(X2,Q),W.T)
 assert set(dp[:,wn==2].flat)=={0,-1,-2} and set(dp[:,wn==4].flat)=={-2,-4,-6}
 # a>2 follows from a^2=(9+6sqrt2)/4>4 (since 72>49).
 assert 72>49
 # Null-null distinct: -2*a^2+2<-6<-2. Opposite signs only decrease it.
 # Null-X same sign, b=2 -> -2a+sqrt6=-sqrt3; b>=3 -> <-2.
 # b=3: (-sqrt6-3sqrt3)/2 < -2 follows already from 3sqrt3>4.
 assert 27>16
 pp=XX4+12;ok=(pp==0)|(pp==-4)|(pp<=-8);np.fill_diagonal(ok,False)
 assert set(ok.sum(axis=1))=={239};components=[];seen=set()
 for i in range(720):
  if i in seen:continue
  c=[i]+np.flatnonzero(ok[i]).tolist();assert len(c)==240
  assert np.all(ok[np.ix_(c,c)]|np.eye(240,dtype=bool));assert not(seen&set(c));seen.update(c);components.append(c)
 assert len(components)==3 and len(seen)==720
 assert set(pp[~ok & ~np.eye(720,dtype=bool)])=={-6,-2,2}
 # Forbidden negative values have rational noninteger 2cos(theta)=3/2,1/2.
 # Positive values cannot be acute-compatible. All opposite-height pairs dotted.
 assert int(XX4.max())-12==-16
 # The full sqrt2 bank is integral for diag(Q,2), including norm-four roots.
 assert np.all((Q@W[wn==4].T)%2==0)
 # Conversely every X cap has a nonintegral length-two cyclic product.
 bad=[]
 for i,row in enumerate(xp):
  j=next(j for j,p in enumerate(row)if F(int(p*p),2*int(wn[j])).denominator!=1)
  bad.append({'cap':i,'old_wall':j,'pair':str(F(int(row[j]),2)),'cyclic_product':str(F(int(row[j]*row[j]),2*int(wn[j])))})
 assert len(bad)==720
 u=np.array([4]+[0]*8+[15,10,20,30,25,20,16,12]+[0]*4+[8],dtype=np.int64)
 v=100*u+h;z=2400
 q=int(v@Q@v)+z*z;assert q==5527986 and int(v@Q@h)==-16014
 oldmax=int(max(W@Q@v));dm=int(max(D@Q@v));xm=int(max(X2@Q@v));assert (oldmax,dm,xm)==(-1,-3206,-8414)
 # D cap: a*dm+2400sqrt2 < 2dm+2400sqrt2 <0.
 assert 2*dm<0 and (2*dm)**2>2*z*z
 # X cap: xm/2+2400sqrt3 <0, both height signs.
 assert xm<0 and xm*xm>12*z*z
 # Actual finite vertices at (x,+/-1/sqrt3): old X incidence is PD21;
 # null-cap value at that point is -ab+sqrt(2/3)<0 for b>=2.
 vertices=[]
 for i in [c[0]for c in components]:
  incids=np.flatnonzero(xp[i]==0).tolist();assert len(incids)==21
  assert inertia(G[np.ix_(incids,incids)])==[21,0,0]
  assert all(p<0 for j,p in enumerate(xp[i])if j not in incids)
  assert all(XX4[i,j]<-4 for j in range(720)if j!=i)
  vertices.append({'sqrt3_cap':i,'old_incident':incids,'point':'(x,+/-1/sqrt3)','square':'-2/3'})
 def counts(M):return dict(sorted(Counter(map(str,M.flat)).items()))
 out={'status':'PASS_EXACT_MIXED_HEIGHT_THREE_1334_FACET_CHAMBERS_WITH_UNIVERSAL_STRICT_ESCAPE',
  'sources':{p:sha(p)for p in(SOURCE,BANK,SYM)},'ambient_signature':[22,1,0],'center':h.tolist(),
  'complete_binary_type_counts':type_counts,'generator_orbit_sizes':{k:len(v)for k,v in orbits.items()},
  'null_projection_count':322,'null_projections':D.tolist(),'sqrt3_twice_projections':X2.tolist(),
  'null_cap_formula':'(a*delta,+/-sqrt(2)); a=(sqrt(6)+sqrt(3))/2; norm2',
  'sqrt3_cap_formula':'(x,+/-sqrt(3)); x=twice_old_projection/2; norm2',
  'null_null_old_gram_counts':counts(DD),'null_sqrt3_twice_old_gram_counts':counts(DX2),
  'sqrt3_four_times_old_gram_counts':counts(XX4),'sqrt3_compatibility_components':components,
  'each_chamber_facets':1334,'each_chamber_composition':{'old':210,'scaled_null_both_signs':644,'sqrt3_both_signs':480},
  'new_cross_angles':'b=2 and equal height signs gives pi/6; all other null-X pairs dotted.',
  'actual_facets':'Every selected normal is spacelike, all distinct selected pairings are nonpositive, and (h,0) is strict timelike. Projecting that point to a wall gives strictness on every other wall.',
  'representative_actual_vertices':vertices,'nonintegral_cycles':bad,
  'pure_unscaled_sqrt2_integrality':'All null coordinates integral; old norm4 Q*w even. Reflections preserve Z^23 for diag(Q,2). This statement is not transferred to arbitrary rational additions.',
  'universal_escape':{'old_coordinates':v.tolist(),'last_coordinate':z,'square':q,'center_pair':-16014,'max_old_pair':oldmax,'max_delta_pair':dm,'max_twice_x_pair':xm,'scope':'Strict on all210 old walls, all644 scaled-null caps and all1440 sqrt3 caps, even across incompatible choices.'},
  'scope':'Three exact infinite-volume Coxeter chambers in H22; no finite-volume, global novelty, or unrestricted impossibility claim. The witness excludes every subset of this fixed combined dictionary.'}
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
