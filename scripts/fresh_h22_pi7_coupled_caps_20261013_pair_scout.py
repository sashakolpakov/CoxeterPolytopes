#!/usr/bin/env python3
"""Try finite coupling between two new all-right caps at one old vertex."""
import json,math
from fractions import Fraction
import numpy as np
from scipy.optimize import brentq
import fresh_h22_pi7_edgewalk_20261013_scout as B
P=B.P;TAG='fresh_h22_pi7_coupled_caps_20261013'
def run():
 M,W,h,G,n,geo,cert=B.load();n=n.copy();rows=json.loads((P/'results/fresh_h22_pi7_edgewalk_20261013_scout.json').read_text())['rows']
 for i,rec in enumerate(cert['new_caps']):
  k=812+i
  for j,ab in enumerate(rec['old_pair_intervals']):G[k,j]=G[j,k]=sum(map(int,ab))/(2*B.V.S)
  for j,ab in enumerate(rec['previous_cap_pair_intervals']):G[k,812+j]=G[812+j,k]=sum(map(int,ab))/(2*B.V.S)
  G[k,k]=n[k]=float(rec['norm'])
 hp=W@M@h
 for i,rec in enumerate(cert['new_caps']):hp[812+i]=sum(map(int,rec['center_pairing']))/(2*B.V.S)
 target=np.array(geo['strict_escape']['numerator'],float)/geo['strict_escape']['denominator'];tp=W@M@target;found=[];pairmax=[]
 for vi,vr in enumerate(geo['new_finite_vertices']):
  ids=vr['base_walls']+[812+vr['new_cap_index']];AI=np.linalg.inv(G[np.ix_(ids,ids)]);lam=AI@hp[ids];qp=-14-hp[ids]@lam;pp=(hp-lam@G[ids])/np.sqrt(-qp);pt=(h@M@target-lam@tp[ids])/np.sqrt(-qp);edges=[]
  for ei,row in enumerate(rows):
   if row['vertex']!=vi or row['status']!='LIGHT_CONE_ESCAPE_CANDIDATE' or max(row['edge'])<812:continue
   j=ids.index(row['released']);qu=AI[j,j];up=(AI[j]@G[ids])/np.sqrt(qu);ut=(AI[j]@tp[ids])/np.sqrt(qu);non=np.ones(len(W),bool);non[row['edge']]=False
   def margin(t):return max(np.sqrt(2)*(np.sinh(t)*pp[non]-np.cosh(t)*up[non])+np.sqrt(2*n[non]))
   lo=0. if margin(0)<0 else brentq(margin,0.,16.);boost=Fraction(math.ceil(math.exp(lo)*64)+1,64);ch=float((boost+1/boost)/2);sh=float((boost-1/boost)/2)
   edges.append({'row':ei,'j':j,'lower':lo,'boost':boost,'ch':ch,'sh':sh,'up':up,'ut':ut,'non':non})
  for i,a in enumerate(edges):
   for b in edges[i+1:]:
    corr=AI[a['j'],b['j']]/np.sqrt(AI[a['j'],a['j']]*AI[b['j'],b['j']]);maxpair=2*(-np.sinh(a['lower'])*np.sinh(b['lower'])+corr*np.cosh(a['lower'])*np.cosh(b['lower']));pairmax.append(float(maxpair))
    for aa,bb in((a,b),(b,a)):
     A=aa['sh'];C0=corr*aa['ch'];delta=A*A-C0*C0
     if abs(delta)<1e-10:continue
     for m in(2,3,4,5,6,7,8,10,12,16):
      c=math.cos(math.pi/m);rad=c*c+delta
      if rad<=0:continue
      for branch in(-1,1):
       sh=(A*c+branch*C0*np.sqrt(rad))/delta;ch=(C0*c+branch*A*np.sqrt(rad))/delta
       if sh<=0 or ch<1 or abs(ch*ch-sh*sh-1)>1e-6 or math.asinh(sh)<=bb['lower']+1e-7:continue
       pb=np.sqrt(2)*(sh*pp-ch*bb['up']);margin=max(pb[bb['non']]+np.sqrt(2*n[bb['non']]))
       if margin>=-1e-7:continue
       found.append({'vertex':vi,'first_edge':aa['row'],'second_edge':bb['row'],'first_boost':str(aa['boost']),'label':m,'branch':branch,'second_sinh':float(sh),'second_cosh':float(ch),'first_lower_parameter':aa['lower'],'second_lower_parameter':bb['lower'],'second_dotted_margin':float(margin),'first_target_pairing':float(np.sqrt(2)*(aa['sh']*pt-aa['ch']*aa['ut'])),'second_target_pairing':float(np.sqrt(2)*(sh*pt-ch*bb['ut'])),'dual_correlation':float(corr),'delta':float(delta)})
  print('vertex',vi,'pairs',len(pairmax),'found',len(found),flush=True)
 found.sort(key=lambda r:(r['label'],max(abs(r['first_target_pairing']),abs(r['second_target_pairing']))))
 out={'status':'NUMERICAL_FINITE_COUPLED_ENDPOINT_PAIR_SCOUT','candidates':found,'pair_count':len(pairmax),'maximum_at_dotted_threshold':max(pairmax),'scope':'Bounded numerical angle/boost proposals for two all-right caps; no exact validity or complete-family claim.'}
 (P/'data/enumerations'/f'{TAG}_pair_scout.json').write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(found),'maxpair',max(pairmax),'first',found[:2],flush=True)
if __name__=='__main__':run()
