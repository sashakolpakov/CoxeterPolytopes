#!/usr/bin/env python3
"""Exact finite-profile/all-real-height obstruction on the sixteen new edges.

Every necessary q(t)<=0 is relaxed to q_lower(t)<=0, using rational lower
coefficient bounds and t>0. Rational outward root enclosures suffice; there
is no numerical optimizer or floating point in this checker. Isolated
I2(m>=6) profiles are explicitly outside its finite-angle profile scope.
"""
import argparse,json,sys
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_component_bridge_20261013_verify as N
V=N.V;K=N.K;BASE=N.PREFIX;OUT=BASE+'_continuous_verify.json'
Z4=(Q(0),)*4
def a4(x,y):return tuple(a+b for a,b in zip(x,y,strict=True))
def s4(x,t):return tuple(a*t for a in x)
def sg4(x):
 a,b=x[:2],x[2:];sa=K.sign(a);sb=K.sign(b)
 if not sa:return sb
 if not sb or sa==sb:return sa
 return sa*K.sign(K.sub(K.mul(a,a),K.scale(K.mul(b,b),2)))
def iv4(x):return V.add(V.golden(x[:2]),V.mul(V.golden(x[2:]),V.sqrt(V.rat(2))))
ALPH=[Z4,(Q(1),Q(0),Q(0),Q(0)),(Q(0),Q(0),Q(1),Q(0)),(Q(0),Q(1),Q(0),Q(0))]
PROD=[[Z4 for _ in range(4)]for _ in range(4)]
for i in range(1,4):PROD[1][i]=PROD[i][1]=ALPH[i]
PROD[2][2]=(Q(2),Q(0),Q(0),Q(0));PROD[3][3]=(Q(1),Q(1),Q(0),Q(0));PROD[2][3]=PROD[3][2]=(Q(0),Q(0),Q(0),Q(1))

def profiles(Mi):
 assert all(x[1]==0 and x[0]>=0 for row in Mi for x in row)
 A=[[x[0]for x in row]for row in Mi];rows=[];pruned=0
 def walk(pre,en):
  nonlocal pruned
  i=len(pre)
  if i==8:rows.append((pre,a4((Q(2),Q(0),Q(0),Q(0)),s4(en,-1))));return
  for j in range(4):
   ne=a4(en,s4(PROD[j][j],A[i][i]))
   for k,u in enumerate(pre):ne=a4(ne,s4(PROD[j][u],2*A[i][k]))
   if sg4(a4((Q(2),Q(0),Q(0),Q(0)),s4(ne,-1)))<=0:pruned+=1;continue
   walk(pre+[j],ne)
 walk([],Z4)
 return rows,pruned

def intersect(xs,ys):
 out=[]
 for a,b in xs:
  for c,d in ys:
   lo=max(a,c);hi=d if b is None else b if d is None else min(b,d)
   if hi is None or hi>=lo:out.append((lo,hi))
 return out

def necessary(a,b,c):
 # Coefficients are exact integers: a common positive interval scale is
 # immaterial. The polynomial is a LOWER bound to the actual q(t), t>0.
 if not a:
  if not b:return[(Q(0),None)]if c<=0 else[]
  z=Q(-c,b)
  return[(max(Q(0),z),None)]if b<0 else([(Q(0),z)]if z>=0 else[])
 D=b*b-4*a*c
 if D<0:return[(Q(0),None)]if a<0 else[]
 sd=isqrt(D);su=sd if sd*sd==D else sd+1
 rr=[sorted((Q(-b+sign*sd,2*a),Q(-b+sign*su,2*a)))for sign in(-1,1)];rr.sort(key=lambda r:r[0])
 # For a>0 enlarge the central feasible interval. For a<0 enlarge both
 # external feasible intervals, including uncertain root slivers.
 pieces=[(rr[0][0],rr[1][1])]if a>0 else[(Q(0),rr[0][1]),(rr[1][0],None)]
 return[(max(Q(0),lo),hi)for lo,hi in pieces if hi is None or hi>=max(Q(0),lo)]

def run(verify=False):
 d=V.read(BASE+'_caps.json');certificate=N.check(d);assert certificate==V.read(BASE+'_verify.json')
 G,I,h,g,old,names,center,pairs,_=N.old_state()
 caps=[N.make(s,b,g,old,names,center,pairs)[0]for s,b in zip(N.SPECS,d['branches'],strict=True)]
 for spec,new in zip(N.SPECS,N.NEW,strict=True):
  ff=dict(zip(spec['edge'],spec['labels'],strict=True));ff[spec['anchor']]=spec['anchor_label']
  for n in names:pairs['|'.join(sorted((new,n))) ]={'kind':'finite','label':ff[n]}if n in ff else{'kind':'dotted'}
 pairs['|'.join(sorted(N.NEW))]={'kind':'dotted'};old+=caps;names+=N.NEW
 target=tuple(map(K.parse,d['escape']));yi=list(map(V.golden,target));hi=list(map(V.golden,h));w=list(map(V.golden,[K.neg(x)for x in K.mv(I,target)]))
 edges=[r['edge']for r in certificate['new_vertex_frontier']if len(r['end_walls'])==1];assert len(edges)==16
 counts={'edges':16,'profiles':0,'branches':0,'quadratic_empty':0,'forbidden_seed6':0};stats=[];tails=[];survivors=[]
 for ei,ss in enumerate(edges):
  M=N.C.gram(ss,pairs);Mi=K.inverse(M);assert N.T.P.T.psd(M)=='positive_definite'
  A=[old[names.index(n)]for n in ss];mi=[[V.golden(x)for x in row]for row in Mi]
  def proj(v):
   z=N.combine(V.mv(mi,[V.pair(v,g,r)for r in A]),A)
   return[V.sub(a,b)for a,b in zip(v,z,strict=True)]
  tv=proj(center);qt=V.pair(tv,g,tv);assert qt[1]<0;tv=[V.div(x,V.sqrt(V.neg(qt)))for x in tv]
  ev=proj(w);lam=V.pair(ev,g,tv);ev=[V.add(a,V.mul(lam,b))for a,b in zip(ev,tv,strict=True)];qe=V.pair(ev,g,ev);assert qe[0]>0;ev=[V.div(x,V.sqrt(qe))for x in ev]
  exterior=[j for j,n in enumerate(names)if n not in ss]
  forms=[yi,list(map(V.neg,hi))]+[V.mv(g,old[j])for j in exterior]
  form_names=['current_escape','negative_center']+[names[j]for j in exterior]
  la=[V.dot(tv,f)for f in forms];le=[V.dot(ev,f)for f in forms]
  baseforms=[V.mv(mi,[V.dot(r,f)for r in A])for f in forms]
  rows,pruned=profiles(Mi);counts['profiles']+=len(rows);before=counts.copy()
  for j in range(8):
   if all(M[j][k]==K.Z for k in range(8)if k!=j):tails.append({'edge':ei,'wall':ss[j],'excluded_from_proof':'I2(m), every m>=6'})
  for codes,deficit in rows:
   mags=[iv4(ALPH[j])for j in codes];rho=V.sqrt(iv4(deficit));bvals=[V.neg(V.dot(mags,row))for row in baseforms]
   for branch in(-1,1):
    counts['branches']+=1
    av=[V.mul(V.div(rho,V.rat(2)),V.add(a,V.mul(V.rat(branch),b)))for a,b in zip(la,le,strict=True)]
    cv=[V.mul(V.div(rho,V.rat(2)),V.add(V.neg(a),V.mul(V.rat(branch),b)))for a,b in zip(la,le,strict=True)]
    ranges=[(Q(0),None)];used=[]
    for n,a,b,c in zip(form_names,av,bvals,cv,strict=True):
     new=intersect(ranges,necessary(a[0],b[0],c[0]))
     if new!=ranges:used.append(n)
     ranges=new
     if not ranges:break
    if not ranges:counts['quadratic_empty']+=1;continue
    # All remaining branches are excluded without an angle cutoff: norm2
    # pairings strictly between -1 and 0 are neither right nor pi/m, m>=3.
    assert 'seed6'in form_names;j=form_names.index('seed6');checks=[]
    for lo,up in ranges:
     assert lo>0 and up is not None
     # Subdivision only sharpens an interval enclosure; all endpoints are
     # exact rationals and the union covers the entire necessary interval.
     for k in range(32):
      l=lo+(up-lo)*Q(k,32);u=lo+(up-lo)*Q(k+1,32);t=(V.rat(l)[0],V.rat(u)[1])
      p=V.add(V.add(V.mul(av[j],t),bvals[j]),V.div(cv[j],t))
      assert p[0]>-V.SCALE and p[1]<0,(ei,codes,branch,lo,up,p)
      checks.append(V.enc(p))
    counts['forbidden_seed6']+=1
    survivors.append({'edge':ei,'labels':[[2,3,4,5][j]for j in codes],'branch':branch,'necessary_ranges':[[str(a),str(b)]for a,b in ranges],'restricting_inequalities':used,'seed6_subinterval_pairings':checks})
  stats.append({'edge':ss,'profiles':len(rows),'partial_profiles_pruned':pruned,'branches_empty':counts['quadratic_empty']-before['quadratic_empty'],'branches_forbidden':counts['forbidden_seed6']-before['forbidden_seed6']})
 assert counts=={'edges':16,'profiles':940,'branches':1880,'quadratic_empty':1876,'forbidden_seed6':4}
 out={'status':'PASS_EXACT_FINITE_PROFILE_ALL_HEIGHT_FRONTIER_OBSTRUCTION','sources':{BASE+'_caps.json':V.sha(BASE+'_caps.json'),BASE+'_verify.json':V.sha(BASE+'_verify.json')},'interval_scale':str(V.SCALE),'counts':counts,'edges':stats,'surviving_acute_branches':survivors,'uncovered_I2_tails':tails,
      'proof':'Exact Q(phi,sqrt2) recursion enumerates all PD9 profiles in the declared alphabet; the PD inverse is entrywise nonnegative, so partial energies only increase. Both norm hyperbola branches cover all real residual heights. For t>0 each true necessary quadratic dominates its rational lower-coefficient polynomial; their conservative root intervals have empty intersection on1876 branches. On the four remaining complete necessary intervals, seed6 pairing lies strictly(-1,0), excluding every Coxeter angle and dotted/parallel relations.',
      'scope':'No norm2 acute-compatible cap retaining all104 walls, cutting this exact escape, and making a finite endpoint on one of these16 open PD8 edges can have all8 edge labels in{2,3,4,5}. This does not cover the51 listed isolatedI2(m>=6) profiles, nonspherical endpoints, other edges, wall replacement, or arbitrary new caps. It is not a general compactH9 obstruction.'}
 if verify:assert out==V.read(OUT)
 else:(V.ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],counts,'uncovered_tails',len(tails),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();run(args.verify)
