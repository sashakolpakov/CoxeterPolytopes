#!/usr/bin/env python3
"""Bounded two-zero-arc connector discovery, after a second edge subdivision.

A-X-Y-R replaces A-X-R on old8. Parameter seeds come from complete native-
dotted conic intervals rather than an enlarged old A grid. Numerical only.
"""
import argparse,json,math
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_broken_cap_adjacency_20261017 as T
N,S,H,B,G,M,ROOT=T.N,T.S,T.H,T.B,T.G,T.M,T.ROOT
TAG='fresh_zero_arc_cap_two_connectors_20261017'
WHITE=[0,1,2,3,5,6]

def qset(a,b,c):
 roots=np.roots([a,b,c]if abs(a)>1e-12 else[b,c])if max(abs(a),abs(b))>1e-12 else[]
 xs=[0.]+sorted(set(float(x.real)for x in roots if abs(x.imag)<1e-8 and x.real>0))+[float('inf')];out=[]
 for lo,hi in zip(xs,xs[1:]):
  t=(lo+hi)/2 if np.isfinite(hi)else max(1.,2*lo)
  if(a*t+b)*t+c<=1e-10:out.append((lo,hi))
 return out

def meet(aa,bb):return[(max(a,c),min(b,d))for a,b in aa for c,d in bb if max(a,c)<min(b,d)-1e-9]

def run():
 z,labels,ns,h,profiles,primary=H.setup();actual=[i for i in z['actual_sides']if i not in S.DELETE]
 prior=json.loads((ROOT/'results'/f'{S.TAG}.json').read_text());rb=prior['accepted_low_caps']
 k0,k1=T.a_curve(ns,h);at=1.5;ay=np.sqrt(.5)*(at*(k0+k1)+(-k0+k1)/at)
 assert not T.gates(T.native_values(ay,0,labels,ns,profiles),[7,8])and ay@M@h<0
 af=np.array([ns[8],ay]);iv=np.linalg.inv(af@M@af.T)
 def proj(v):return v-(v@M@af.T)@iv@af
 k0=proj(h);k0/=np.sqrt(-k0@M@k0);k1=proj(ns[30]);k1+=float(k1@M@k0)*k0;k1/=np.sqrt(k1@M@k1)
 xbank=[];counts=Counter();blockers=Counter();accepted=[];scores=Counter();unresolved=[]
 for ma in range(2,6):
  assert T.corner(0,2,ma,0,2)is not None
  part=np.array([0.,-B.theta(ma)])@iv@af;d=2-float(part@M@part)
  for sg in[-1,1]:
   aa=np.sqrt(d)/2*(k0+sg*k1);bb=np.sqrt(d)/2*(-k0+sg*k1);ints=[(0.,float('inf'))]
   for j in labels:
    if j in S.DELETE+WHITE+[8]:continue
    a=float(aa@M@ns[j]);b=float(bb@M@ns[j]);c=float(part@M@ns[j])
    if max(abs(a),abs(b))<1e-9:
     if B.relation(c)is None or B.relation(c)=='parallel_or_unresolved':ints=[];break
     continue
    ints=meet(ints,qset(a,c+2,b))
    if not ints:break
   for lo,hi in ints:
    # One explicit interior rational seed per complete open interval.
    st=math.ceil(max(1.,1.1*lo)*1000)/1000
    if st>=hi:continue
    x=part+aa*st+bb/st;bad=T.gates(T.native_values(x,0,labels,ns,profiles),[8])
    assert not bad and x@M@h<0
    xbank.append(dict(index=len(xbank),AX=ma,branch=sg,parameter=st,continuous_interval=[lo,None if not np.isfinite(hi)else hi],projection=x.tolist()))
 for rr in rb:
  r=np.array(rr['projection']);rpi=rr['recipe']['profile'];r8=rr['recipe']['labels'][0]
  for xx in xbank:
   x=np.array(xx['projection']);xr=float(x@M@r);ar=float(ay@M@r)
   if xr>=-2-1e-7 or ar>=-2-1e-7:counts['wrong_dotted_separation']+=1;continue
   for ypi,y8,yr,yx in product(range(7),range(2,6),range(2,6),range(2,6)):
    if T.corner(ypi,y8,yr,rpi,r8)is None or T.corner(ypi,y8,yx,0,2)is None:continue
    off=H.A.native(np.zeros(4),ypi,np.zeros(4),rpi,profiles)
    roots,info=N.norm_roots([ns[8],r,x],[-B.theta(y8),-B.theta(yr)-off,-B.theta(yx)],profiles[ypi]['square']);counts['Y_equations']+=1
    if info['status'].startswith('unresolved'):
     unresolved.append(dict(R=rr['index'],X=xx['index'],recipe=[ypi,y8,yr,yx],norm=info));continue
    for sg,y in roots:
     counts['Y_branches']+=1
     assert abs(y@M@y-profiles[ypi]['square'])<2e-6
     assert max(abs(float(y@M@v)-c)for v,c in zip([ns[8],r,x],[-B.theta(y8),-B.theta(yr)-off,-B.theta(yx)]))<2e-6
     if y@M@h+profiles[ypi]['center_offset']>=-1e-8:counts['Y_wrong_center']+=1;continue
     vals=T.native_values(y,ypi,labels,ns,profiles);bad=T.gates(vals,[8]);ya=float(y@M@ay)
     if B.relation(ya)is None or B.relation(ya)=='parallel_or_unresolved':bad.append(('A',ya))
     if bad:blockers[tuple(str(i)for i,v in bad)]+=1;continue
     counts['coherent_four_cap_prefixes']+=1
     nd=dict(ns);nd.update({67:ay,68:r,69:x,70:y});bd=S.C.detailed_cone(nd,actual+[67,68,69,70]);scores[str((bd['counts'],bd['open_ends'],bool(bd['warnings'])))]+=1
     sq={v['wall']:B.mid(v['square'])for v in z['projections']}
     ic=H.incidence_score([nd[i]for i in actual+[67,68,69,70]],actual+[67,68,69,70],[sq[i]for i in actual]+[2,profiles[rpi]['square'],2,profiles[ypi]['square']])
     accepted.append(dict(R=rr['index'],X=xx['index'],Y_profile=ypi,Y_labels=[y8,yr,yx],Y_branch=sg,Y_projection=y.tolist(),AY_product=ya,XR_product=xr,AR_product=ar,boundary=bd,PD_incidence=ic))
  if(rr['index']+1)%20==0:print('TWO_CONNECTORS',rr['index']+1,'of94',dict(counts),'accepted',len(accepted),flush=True)
 paths=[Path(__file__),Path(T.__file__),Path(N.__file__),Path(S.__file__),Path(H.__file__),ROOT/'results'/f'{S.TAG}.json',H.PARENT,H.BOUNDARY]
 return dict(status='BOUNDED_NUMERICAL_TWO_ZERO_ARC_CONNECTOR_DOMAIN',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},deleted_walls=S.DELETE,boundary_cycle=S.CYCLE,exact_old_deletion_cost=prior['exact_old_deletion_cost'],
  topology=dict(second_subdivided_edge=['X','R'],old_triangles=[[8,'X','R'],['X','R','B']],new_triangles=[[8,'X','Y'],[8,'Y','R'],['X','Y','B'],['Y','R','B']],vertices=11,edges=25,triangles=15,x=5,y=6,z=4,
   scope='A-X-Y-R chain alongold8; both X,Y have zero original boundary arcs. B,C not yet solved.'),
  A_parameter='3/2',A_projection=ay.tolist(),X_conic_seeds=xbank,counts=dict(counts,R_low_seeds=len(rb),X_seeds=len(xbank)),native_blockers=[[list(k),v]for k,v in sorted(blockers.items())],boundary_histogram=dict(scores),accepted_prefixes=accepted,unresolved=unresolved,
  scope='New six-cap incidence after a second stellar subdivision, not the excluded single-X topology. A parameter3/2 and one rational interior X seed per native-dotted interval are bounded discovery choices. R and Y finite labels2..5 only; high labels, alternative parameters and B,C completion remain open. All retained native contacts/common center and complete numerical boundary are checked for any accepted prefix. No exact chamber or finite-volume claim.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],out['counts'],flush=True)
