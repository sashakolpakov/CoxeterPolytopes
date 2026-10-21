#!/usr/bin/env python3
"""Flip interior AB to CX in the six-cap pentagonal disk; bounded scout."""
import argparse,json
from collections import Counter
from itertools import product
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_two_connector_coupled_disk_20261017 as F
D,T,N,S,H,B,G,M,ROOT=F.D,F.T,F.N,F.S,F.H,F.B,F.G,F.M,F.ROOT
TAG='fresh_two_connector_flipped_disk_20261017'

def run():
 z,labels,ns,h,profiles,primary=H.setup();actual=[i for i in z['actual_sides']if i not in S.DELETE]
 old=json.loads((ROOT/'results'/f'{F.TAG}.json').read_text());counts=Counter();blockers=Counter();full=[];near=[];unresolved=[]
 rs={}
 for mr in old['R44_labels']:
  rr,ii=N.norm_roots([ns[8],ns[44],ns[63]],[0.,-B.theta(mr),0.],2.)
  rs[mr]=[r for sg,r in rr if r@M@h<0 and not T.gates(T.native_values(r,0,labels,ns,profiles),[8,44,63])][0]
 def native(v,pi,assigned=()):
  if v@M@h+profiles[pi]['center_offset']>=-1e-8:return [('center',0)]
  return T.gates(T.native_values(v,pi,labels,ns,profiles),assigned)
 for ib,br in enumerate(old['B_bank']):
  bp=br['B_profile'];m30,m63,mrB=br['B_labels'];b=np.array(br['B_projection']);r=rs[br['R44']]
  assert not native(b,bp,[30,63])
  for my,mby in product(range(2,6),repeat=2):
   if T.corner(bp,mrB,mby,0,my)is None:continue
   ys,yi=N.norm_roots([ns[8],r,b],[0.,-B.theta(my),-B.theta(mby)],2.);counts['Y_equations']+=1
   if yi['status'].startswith('unresolved'):unresolved.append(['Y',ib,my,mby,yi]);continue
   for ysg,y in ys:
    bad=native(y,0,[8])
    if bad:blockers[('Y',str(bad[0][0]))]+=1;continue
    counts['native_Y']+=1
    for mxy,mbx in product(range(2,6),repeat=2):
     if T.corner(bp,mbx,mby,0,mxy)is None:continue
     xs,xi=N.norm_roots([ns[8],b,y],[0.,-B.theta(mbx),-B.theta(mxy)],2.);counts['X_equations']+=1
     if xi['status'].startswith('unresolved'):unresolved.append(['X',ib,my,mby,mxy,mbx,xi]);continue
     for xsg,x in xs:
      bad=native(x,0,[8])
      if bad:blockers[('X',str(bad[0][0]))]+=1;continue
      if x@M@r>=-2-1e-7:counts['X_R_not_dotted']+=1;continue
      counts['native_XYB']+=1
      for maxx in range(2,6):
       ars,ai=N.norm_roots([ns[7],ns[8],x],[0.,0.,-B.theta(maxx)],2.);counts['A_equations']+=1
       if ai['status'].startswith('unresolved'):unresolved.append(['A',ib,my,mby,mxy,mbx,maxx,ai]);continue
       for asg,a in ars:
        bad=native(a,0,[7,8])
        if bad:blockers[('A',str(bad[0][0]))]+=1;continue
        if any(a@M@v>=-2-1e-7 for v in [b,r,y]):counts['A_other_not_dotted']+=1;continue
        counts['coherent_five_caps']+=1
        rec=dict(B_index=ib,Y_labels=[my,mby],X_labels=[mxy,mbx],AX=maxx,A_projection=a.tolist(),R_projection=r.tolist(),X_projection=x.tolist(),Y_projection=y.tolist())
        near.append(rec)
        for cp,m7,mc30,mca in product(range(7),range(2,6),range(2,6),range(2,6)):
         if T.corner(cp,m7,mc30,0,5)is None or T.corner(cp,m7,mca,0,2)is None:continue
         cs,ci=N.norm_roots([ns[7],ns[30],a],[-B.theta(m7),-B.theta(mc30),-B.theta(mca)],profiles[cp]['square']);counts['C_equations']+=1
         if ci['status'].startswith('unresolved'):unresolved.append(['C',len(near)-1,cp,m7,mc30,mca,ci]);continue
         for csg,c in cs:
          bad=native(c,cp,[7,30])
          if bad:blockers[('C',str(bad[0][0]))]+=1;continue
          cb=H.A.native(c,cp,b,bp,profiles);cx=float(c@M@x);mcb=F.label(cb);mcx=F.label(cx)
          if mcb is None or mcx is None:counts['C_BX_gap']+=1;continue
          if T.corner(cp,mca,mcx,0,maxx)is None or T.corner(cp,mcx,mcb,bp,mbx)is None or T.corner(cp,mc30,mcb,bp,m30)is None:counts['C_nonfinite_corner']+=1;continue
          if any(c@M@v>=-2-1e-7 for v in [r,y]):counts['C_other_not_dotted']+=1;continue
          counts['coherent_complete_disks']+=1;nd=dict(ns);nd.update({67:a,68:r,69:x,70:y,71:b,72:c})
          bd=S.C.detailed_cone(nd,actual+[67,68,69,70,71,72])
          full.append(dict(prefix=len(near)-1,C_profile=cp,C_labels=[m7,mc30,mca,mcx,mcb],C_projection=c.tolist(),boundary=bd))
  if (ib+1)%50==0:print('FLIPPED_DISK',ib+1,'of',len(old['B_bank']),dict(counts),'full',len(full),flush=True)
 paths=[Path(__file__),Path(F.__file__),Path(D.__file__),Path(T.__file__),Path(N.__file__),Path(H.__file__),ROOT/'results'/f'{F.TAG}.json',H.PARENT,H.BOUNDARY]
 return dict(status='BOUNDED_NUMERICAL_FLIPPED_TWO_CONNECTOR_DISK',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},counts=dict(counts),prefixes=near,complete_disks=full,unresolved=unresolved,native_blockers=[[list(k),v]for k,v in sorted(blockers.items())],
  deleted_walls=S.DELETE,exact_old_deletion_cost=old['exact_old_deletion_cost'],topology=dict(removed_edge=['A','B'],added_edge=['C','X'],removed_triangles=[['A','X','B'],['A','B','C']],added_triangles=[['A','X','C'],['X','B','C']],vertices=11,edges=25,triangles=15),
  scope='Changed triangulation avoids mandatory finite AB. Inherits bounded native B bank; solves Y,X,A,C sequentially, with7 A2 profiles for B,C and zero for A,X,Y. New assigned finite labels2..5; isolated high tails remain open. All retained native gates and all stipulated new contacts checked; full numerical boundary only for complete disks. No exact geometry or finite-volume claim.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],out['counts'],'full',len(out['complete_disks']),flush=True)
