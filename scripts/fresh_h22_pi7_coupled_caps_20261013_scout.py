#!/usr/bin/env python3
"""Bounded numerical all-right endpoint insertion at the 835-wall frontier."""
import hashlib,json,math
from pathlib import Path
import numpy as np
import fresh_h22_pi7_edgewalk_20261013_scout as B
P=B.P;TAG='fresh_h22_pi7_coupled_caps_20261013'
SRC='results/fresh_h22_pi7_edgewalk_20261013_scout.json'
def run():
 M,W,h,G,n,geo,cert=B.load();n=n.copy();rows=json.loads((P/SRC).read_text())['rows']
 # Use the rigorous pair archives, not cancellation-prone products of long coordinates.
 for i,rec in enumerate(cert['new_caps']):
  k=812+i
  for j,ab in enumerate(rec['old_pair_intervals']):G[k,j]=G[j,k]=(int(ab[0])+int(ab[1]))/(2*B.V.S)
  for j,ab in enumerate(rec['previous_cap_pair_intervals']):G[k,812+j]=G[812+j,k]=(int(ab[0])+int(ab[1]))/(2*B.V.S)
  G[k,k]=float(rec['norm']);n[k]=G[k,k]
 hp=W@M@h
 for i,rec in enumerate(cert['new_caps']):hp[812+i]=sum(map(int,rec['center_pairing']))/(2*B.V.S)
 target=np.array(geo['strict_escape']['numerator'],dtype=float)/geo['strict_escape']['denominator'];tp=W@M@target
 nulls=np.array([r['null_point']for r in rows]);nullpairs=W@M@nulls.T;hnull=h@M@nulls.T
 found=[];miss=[]
 for ei,row in enumerate(rows):
  if row['status']!='LIGHT_CONE_ESCAPE_CANDIDATE' or max(row['edge'])<812:continue
  vrec=geo['new_finite_vertices'][row['vertex']];ids=vrec['base_walls']+[812+vrec['new_cap_index']];j=ids.index(row['released']);A=G[np.ix_(ids,ids)];AI=np.linalg.inv(A);b=hp[ids];lam=AI@b
  qp=-14-b@lam;qu=AI[j,j];assert qp<0 and qu>0
  pp=(hp-lam@G[ids])/np.sqrt(-qp);up=(AI[j]@G[ids])/np.sqrt(qu);non=np.ones(len(W),bool);non[row['edge']]=False
  p=h-lam@W[ids];u=AI[j]@W[ids];tt=(h@M@target-lam@tp[ids])/np.sqrt(-qp);tu=(AI[j]@tp[ids])/np.sqrt(qu)
  for k in range(1,21):
   ch=(2.**k+2.**-k)/2;sh=(2.**k-2.**-k)/2
   ps=np.sqrt(2)*(sh*pp-ch*up)
   if np.max(ps[non]+np.sqrt(2*n[non])) < -1e-7:
    normal=np.sqrt(2)*(sh*p/np.sqrt(-qp)-ch*u/np.sqrt(qu));cut=np.sqrt(2)*(sh*tt-ch*tu)
    allnull=np.sqrt(2)*(sh*(hnull-lam@nullpairs[ids])/np.sqrt(-qp)-ch*(AI[j]@nullpairs[ids])/np.sqrt(qu))
    found.append({'edge_row':ei,'vertex':row['vertex'],'incident':ids,'released':row['released'],'edge':row['edge'],'boost_exponent':k,'center_projection_square':float(qp),'dual_square':float(qu),'maximum_dotted_margin':float(np.max(ps[non]+np.sqrt(2*n[non]))),'target_pairing':float(cut),'normal':normal.tolist(),'cut_edge_limits':np.flatnonzero(allnull>1e-6).tolist()});break
  else:miss.append(ei)
 found.sort(key=lambda r:(r['target_pairing']<=0,r['boost_exponent'],-len(r['cut_edge_limits']),-r['target_pairing']))
 out={'status':'BOUNDED_NUMERICAL_ALL_RIGHT_ENDPOINT_CAP_SCOUT','source_sha256':{s:hashlib.sha256((P/s).read_bytes()).hexdigest()for s in(SRC,B.V.OUT,'results/'+B.V.PRE+'_geometry.json')},'boost_exponents':[1,20],'candidates':found,'no_candidate_in_boost_range':miss,'scope':'Numerical proposals only; every proposed cap retains all835 and is right to21 incident walls including a pi7 cap. No complete cap-family or volume claim.'}
 (P/'data/enumerations'/f'{TAG}_scout.json').write_text(json.dumps(out,indent=2)+'\n')
 print('PASS_SCOUT',len(found),'failed',len(miss),'targetcutters',sum(r['target_pairing']>1e-7 for r in found),'best',[(r['edge_row'],r['boost_exponent'],len(r['cut_edge_limits']),r['target_pairing'])for r in found[:10]],flush=True)
if __name__=='__main__':run()
