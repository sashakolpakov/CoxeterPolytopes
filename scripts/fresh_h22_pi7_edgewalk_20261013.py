#!/usr/bin/env python3
"""Exact interval geometry of506 outgoing directions of the835-wall chamber."""
import argparse,json,hashlib
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import numpy as np
import fresh_h22_incidence_prism_20261013_verify as V
P=V.P;TAG='fresh_h22_pi7_edgewalk_20261013';OUT=P/'results'/f'{TAG}.json'
def div(a,b):
 assert not b[0]<=0<=b[1]
 if b[1]<0:return div(V.neg(a),V.neg(b))
 return V.mul(a,(V.S*V.S//b[1],-((-V.S*V.S)//b[0])))
def inv(A):
 n=len(A);B=[list(r)+[V.rat(int(i==j))for j in range(n)]for i,r in enumerate(A)];piv=[]
 for i in range(n):
  z=B[i][i];assert z[0]>0;piv.append(z);B[i]=[div(x,z)for x in B[i]]
  for j in range(n):
   if j==i:continue
   z=B[j][i]
   if z==(0,0):continue
   B[j]=[V.add(x,V.neg(V.mul(z,y)))for x,y in zip(B[j],B[i])]
 return[r[n:]for r in B],piv
def ip(a,b):return V.isum(V.mul(x,y)for x,y in zip(a,b))
def enc(v):return list(map(str,v))
def setup():
 M,R,h=V.setup();cert=V.read(V.OUT);geo=V.read('results/'+V.PRE+'_geometry.json');scout=V.read('results/'+TAG+'_scout.json')
 for d in(cert,geo,scout):
  for p,x in d['sources'].items():assert V.sha(p)==x
 mm=np.array(M,dtype=np.int64);rr=np.array(R,dtype=np.int64);assert 23**2*int(abs(mm).max())*int(abs(rr).max())**2<2**63;gg=rr@mm@rr.T;oldh=(rr@(mm@np.array(h,dtype=np.int64))).tolist();hh=[V.rat(x)for x in oldh]+[tuple(map(int,r['center_pairing']))for r in cert['new_caps']]
 def pair(i,j):
  if i<812 and j<812:return V.rat(int(gg[i,j]))
  if i<j:i,j=j,i
  if i==j:return V.rat(F(cert['new_caps'][i-812]['norm']))
  row=cert['new_caps'][i-812]
  if j<812:
   if row['finite_contacts'].get(str(j))==2:return(0,0)
   return tuple(map(int,row['old_pair_intervals'][j]))
  return tuple(map(int,row['previous_cap_pair_intervals'][j-812]))
 return M,R,h,cert,geo,scout,hh,pair

def run(verify=False):
 M,R,h,cert,geo,scout,hh,pair=setup();rows=[];counts=Counter();sources=[V.OUT,V.REC,'results/'+V.PRE+'_geometry.json','results/'+TAG+'_scout.json','scripts/'+Path(__file__).name,'scripts/'+V.PRE+'_verify.py'];vertices=[]
 for vi,v in enumerate(geo['new_finite_vertices']):
  inc=v['base_walls']+[812+v['new_cap_index']];A=[[pair(i,j)for j in inc]for i in inc];AI,piv=inv(A);b=[hh[i]for i in inc];coef=[ip(row,b)for row in AI];pn=V.add(V.rat(-14),V.neg(ip(b,coef)));assert pn[1]<0
  cross=[[pair(i,j)for i in inc]for j in range(835)];pp=[V.add(hh[j],V.neg(ip(coef,cross[j])))if j not in inc else(0,0)for j in range(835)];assert all(x[1]<0 for j,x in enumerate(pp)if j not in inc)
  vertices.append({'vertex':vi,'incident':inc,'positive_Gram_pivots':[enc(x)for x in piv],'point_norm':enc(pn),'maximum_nonincident_pairing':enc(max((x for j,x in enumerate(pp)if j not in inc),key=lambda x:x[1]))})
  for local,released in enumerate(inc):
   row=scout['rows'][vi*22+local];assert(row['vertex'],row['released'])==(vi,released);edge=[i for i in inc if i!=released];col=[AI[i][local]for i in range(22)];un=AI[local][local];assert un[0]>0
   up=[ip(col,cross[j])if j not in inc else V.rat(int(j==released))for j in range(835)]
   rec={'vertex':vi,'released':released,'edge':edge,'dual_square':enc(un)}
   if row['status']=='LIGHT_CONE_ESCAPE_CANDIDATE':
    alpha0=V.sqrt(div(V.neg(pn),un));found=False
    for power in range(4,61,4):
     eps=F(1,2**power);alpha=V.scale(alpha0,1+eps);pairs=[V.add(pp[j],V.neg(V.mul(alpha,up[j])))for j in range(835)if j not in edge]
     if max(x[1]for x in pairs)<0:found=True;break
    assert found,(vi,released,'no exact outward continuation')
    # All retained21 pairings vanish as algebraic identities. p.u=0 and
    # u²=AI_jj, so q(p-alpha*u)=(-p²)(2eps+eps²)>0 exactly.
    q=V.scale(V.neg(pn),2*eps+eps*eps);future=V.add(pn,V.neg(V.mul(alpha,coef[local])));assert q[0]>0 and future[1]<0
    rec.update(status='CERTIFIED_UNCLOSED_EDGE_RAY',beyond_null_factor=str(1+eps),alpha_interval=enc(alpha),spacelike_norm=enc(q),center_pairing=enc(future),maximum_nonedge_pairing=enc(max(pairs,key=lambda x:x[1])),retained_zero_pairings=21)
   elif row['status']=='EXISTING_FINITE_ENDPOINT':
    other=[x['wall']for x in row['finite_endpoints']if x['wall']!=released];assert len(other)==1;j=other[0];assert pp[j][1]<0 and up[j][1]<0;alpha=div(pp[j],up[j]);assert alpha[0]>0
    q=V.add(pn,V.mul(V.mul(alpha,alpha),un));assert q[1]<0;future=V.add(pn,V.neg(V.mul(alpha,coef[local])));assert future[1]<0
    pairs=[V.add(pp[k],V.neg(V.mul(alpha,up[k])))for k in range(835)if k not in edge and k!=j];assert max(x[1]for x in pairs)<0
    # PD22 verifies an actual finite corner independently of endpoint ordering.
    E=edge+[j];EI,epiv=inv([[pair(i,k)for k in E]for i in E]);rec.update(status='CERTIFIED_EXISTING_FINITE_ENDPOINT',next_wall=j,endpoint_incident=E,alpha_interval=enc(alpha),endpoint_norm=enc(q),center_pairing=enc(future),maximum_nonincident_pairing=enc(max(pairs,key=lambda x:x[1])),positive_Gram_pivots=[enc(x)for x in epiv])
   else:
    assert row['status']=='POSSIBLE_IDEAL_ENDPOINT_NEEDS_LINK';rec.update(status='IDEAL_ENDPOINT_AUDIT_PENDING',near_zero_candidates=row['near_zero_schur_candidates'])
   counts[rec['status']]+=1;rows.append(rec)
  print('vertex',vi,'directions',22,flush=True)
 out={'status':'PASS_EXACT_835_WALL_FINITE_VERTEX_FRONTIER','sources':{s:V.sha(s)for s in sources},'interval_denominator':str(V.S),'known_actual_vertices':23,'directed_edges':506,'counts':dict(counts),'vertices':vertices,'rows':rows,'scope':'Exact Gram-inverse identities and rigorous intervals certify all497 proposed unclosed edge rays against all835 existing walls and seven finite endpoint directions. Two possible ideal endpoints are explicitly left to a separate null-ray/link certificate; they are not called open. The starting23 actual finite vertices are not claimed exhaustive. No mirror-existence or full Allcock-algorithm assertion is made.'}
 if verify:assert out==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],dict(counts),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
