#!/usr/bin/env python3
"""Solve B first, then A,X,Y,C for the two-zero-arc six-cap disk.

Bounded numerical discovery with all retained native gates. Unknown labels
and high-label tails are not rejected by this scout.
"""
import argparse,json
from collections import Counter
from itertools import product
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_two_connector_native53_deformation_20261017 as D
T,N,S,H,B,G,M,ROOT=D.T,D.N,D.S,D.H,D.B,D.G,D.M,D.ROOT
TAG='fresh_two_connector_coupled_disk_20261017'

def label(v):
 r=B.relation(float(v))
 return r if isinstance(r,int)and r>=2 else None

def run():
 z,labels,ns,h,profiles,primary=H.setup();oldactual=[i for i in z['actual_sides']if i not in S.DELETE]
 parent=json.loads((ROOT/'results'/f'{S.TAG}.json').read_text());k0,k1=T.a_curve(ns,h)
 counts=Counter();blocks=Counter();bank=[];coupled=[];full=[];unresolved=[]
 # This is a changed high-R geometric branch; finite labels beyond this
 # explicitly bounded set remain open, including isolated high B contacts.
 rlabels=[6,7,8,9,10,12,15]
 for mr in rlabels:
  roots,info=N.norm_roots([ns[8],ns[44],ns[63]],[0.,-B.theta(mr),0.],2.)
  for rsg,r in roots:
   if r@M@h>=0 or T.gates(T.native_values(r,0,labels,ns,profiles),[8,44,63]):continue
   counts['native_R']+=1
   for bp,m30,m63,mrB in product(range(7),range(2,7),range(2,6),range(2,6)):
    if T.corner(bp,m63,m30,0,2)is None or T.corner(bp,m63,mrB,0,2)is None:continue
    bs,bi=N.norm_roots([ns[30],ns[63],r],[-B.theta(m30),-B.theta(m63),-B.theta(mrB)],profiles[bp]['square']);counts['B_equations']+=1
    if bi['status'].startswith('unresolved'):unresolved.append(['B',mr,bp,m30,m63,mrB,bi]);continue
    for bsg,b in bs:
     counts['B_branches']+=1
     if b@M@h+profiles[bp]['center_offset']>=-1e-8:counts['B_wrong_center']+=1;continue
     bad=T.gates(T.native_values(b,bp,labels,ns,profiles),[30,63])
     if bad:blocks[('B',str(bad[0][0]))]+=1;continue
     if b@M@ns[8]>=-2-1e-7:counts['B_incident8']+=1;continue
     br=dict(R44=mr,B_profile=bp,B_labels=[m30,m63,mrB],B_branch=bsg,B_projection=b.tolist());bank.append(br);counts['native_B']+=1
     for asign,mab in product([-1,1],range(2,6)):
      aa=2**-.5*(k0+asign*k1);bb=2**-.5*(-k0+asign*k1)
      ts=np.roots([aa@M@b,B.theta(mab),bb@M@b])
      for t in ts:
       if abs(t.imag)>1e-8 or t.real<=0:continue
       t=float(t.real);a=aa*t+bb/t;counts['A_events']+=1
       if a@M@h>=-1e-8:continue
       bad=T.gates(T.native_values(a,0,labels,ns,profiles),[7,8])
       if bad:blocks[('A',str(bad[0][0]))]+=1;continue
       if a@M@r>=-2-1e-7:counts['A_not_dotted_R']+=1;continue
       counts['native_AB']+=1;xb=[];yb=[]
       for side,endpoint,otherbond,target in [('X',a,mab,xb),('Y',r,mrB,yb)]:
        for mend,mb in product(range(2,6),repeat=2):
         if T.corner(bp,otherbond,mb,0,mend)is None:continue
         vs,vi=N.norm_roots([ns[8],endpoint,b],[0.,-B.theta(mend),-B.theta(mb)],2.)
         counts[side+'_equations']+=1
         if vi['status'].startswith('unresolved'):unresolved.append([side,mr,bp,mab,mend,mb,vi]);continue
         for vsg,v in vs:
          if v@M@h>=-1e-8:continue
          if T.gates(T.native_values(v,0,labels,ns,profiles),[8]):continue
          opposite=r if side=='X'else a
          if v@M@opposite>=-2-1e-7:continue
          counts['native_'+side]+=1;target.append((mend,mb,vsg,v))
       for xr,yr in product(xb,yb):
        mx,mbx,xsg,x=xr;my,mby,ysg,y=yr;xy=float(x@M@y);mxy=label(xy)
        if mxy is None or T.corner(bp,mbx,mby,0,mxy)is None:continue
        counts['coherent_five_cap_prefixes']+=1
        rec=dict(**br,A_AB=mab,A_sign=asign,A_t=t,X_labels=[mx,mbx],Y_labels=[my,mby],XY_label=mxy,XY_product=xy,
          A_projection=a.tolist(),R_projection=r.tolist(),X_projection=x.tolist(),Y_projection=y.tolist())
        nd=dict(ns);nd.update({67:a,68:r,69:x,70:y,71:b})
        rec['boundary']=S.C.detailed_cone(nd,oldactual+[67,68,69,70,71]);coupled.append(rec)
        for cp,m7,mc30,mca in product(range(7),range(2,6),range(2,6),range(2,6)):
         if T.corner(cp,m7,mc30,0,5)is None or T.corner(cp,m7,mca,0,2)is None:continue
         cs,ci=N.norm_roots([ns[7],ns[30],a],[-B.theta(m7),-B.theta(mc30),-B.theta(mca)],profiles[cp]['square']);counts['C_equations']+=1
         if ci['status'].startswith('unresolved'):unresolved.append(['C',mr,cp,m7,mc30,mca,ci]);continue
         for csg,c in cs:
          if c@M@h+profiles[cp]['center_offset']>=-1e-8:continue
          if T.gates(T.native_values(c,cp,labels,ns,profiles),[7,30]):continue
          cb=H.A.native(c,cp,b,bp,profiles);mcb=label(cb)
          if mcb is None:continue
          if T.corner(cp,mca,mcb,bp,mab)is None or T.corner(cp,mc30,mcb,bp,m30)is None:continue
          if any(c@M@v>=-2-1e-7 for v in [r,x,y]):continue
          counts['coherent_complete_disks']+=1;dd=dict(nd);dd[72]=c
          bd=S.C.detailed_cone(dd,oldactual+[67,68,69,70,71,72])
          full.append(dict(prefix=len(coupled)-1,C_profile=cp,C_labels=[m7,mc30,mca,mcb],C_projection=c.tolist(),boundary=bd))
  print('COUPLED_DISK',mr,dict(counts),'full',len(full),flush=True)
 paths=[Path(__file__),Path(D.__file__),Path(T.__file__),Path(N.__file__),Path(H.__file__),ROOT/'results'/f'{S.TAG}.json',H.PARENT,H.BOUNDARY]
 return dict(status='BOUNDED_NUMERICAL_COUPLED_TWO_CONNECTOR_DISK',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},counts=dict(counts),R44_labels=rlabels,B_bank=bank,prefixes=coupled,complete_disks=full,unresolved=unresolved,
  native_blockers=[[list(k),v]for k,v in sorted(blocks.items())],deleted_walls=S.DELETE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],
  scope='B is solved from30,63,R; A from7,8 and its finite AB event; X/Y from8 and endpoint/B; C from7,30,A. All retained native gates/common-center checked before scoring. B uses7 A2 profiles, B30 labels2..6 and other finite labels2..5; A/X/Y profile0. High-label tails and other R seeds remain open. Numerical finite labels and cone decisions are not exact certificates; no finite-volume or all-parameter obstruction claim.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],out['counts'],'full',len(out['complete_disks']),flush=True)
