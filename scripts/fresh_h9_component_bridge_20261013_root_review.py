#!/usr/bin/env python3
"""Separate 110-digit linear-system reconstruction of the two new normals.

The preceding normals come from their rigorous enclosures. This review uses
a coordinate nullspace and a scalar quadratic, not the producer's sequential
orthogonal-projection recipe. It supplements, not replaces, exact replay.
"""
import argparse,json,sys
from pathlib import Path
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_component_bridge_20261013_verify as N
V=N.V;P=Path(__file__).resolve().parents[1]
OUT='results/fresh_h9_component_bridge_20261013_root_review.json'

def run():
 mp.mp.dps=110
 num=lambda iv:mp.mpf(iv[0]+iv[1])/(2*int(V.SCALE))
 d=V.read(N.PREFIX+'_caps.json');saved=V.read(N.PREFIX+'_verify.json')
 _,_,_,gg,rr,names,_,_,_=N.old_state()
 g=mp.matrix([[num(x)for x in row]for row in gg]);old=[mp.matrix(list(map(num,r)))for r in rr]
 check=lambda x,y:(x.T*g*y)[0]
 roots=[];rows=[]
 for spec,branch,enc in zip(N.SPECS,d['branches'],saved['normal_enclosures'],strict=True):
  finite=spec['edge']+[spec['anchor']];labels=spec['labels']+[spec['anchor_label']]
  M=mp.matrix([list((g*old[names.index(n)]).T)for n in finite])
  rhs=mp.matrix([-2*mp.cos(mp.pi/m)for m in labels])
  for free in range(10):
   js=[j for j in range(10)if j!=free];A=mp.matrix([[M[i,j]for j in js]for i in range(9)])
   if abs(mp.det(A))>mp.mpf('1e-30'):break
  else:raise AssertionError('No independent coordinate minor')
  y=mp.lu_solve(A,rhs);z=mp.lu_solve(A,mp.matrix([-M[i,free]for i in range(9)]))
  v=mp.matrix(10,1);n=mp.matrix(10,1);n[free]=1
  for j,a,b in zip(js,y,z,strict=True):v[j]=a;n[j]=b
  a=check(n,n);b=2*check(v,n);c=check(v,v)-2;disc=b*b-4*a*c;assert disc>0
  candidates=[v+((-b+s*mp.sqrt(disc))/(2*a))*n for s in(-1,1)]
  expected=mp.matrix([mp.mpf(int(x[0])+int(x[1]))/(2*int(V.SCALE))for x in enc])
  r=min(candidates,key=lambda x:mp.norm(x-expected));error=mp.norm(r-expected);assert error<mp.mpf('1e-45')
  assert abs(check(r,r)-2)<mp.mpf('1e-70')
  finite_error=max(abs(check(r,old[names.index(k)])+2*mp.cos(mp.pi/m))for k,m in zip(finite,labels,strict=True))
  assert finite_error<mp.mpf('1e-70')
  dotted=max(check(r,w)for name,w in zip(names,old,strict=True)if name not in finite);assert dotted<-2
  roots.append(r);rows.append({'free_coordinate':free,'normal_difference':mp.nstr(error,12),
   'finite_equality_error':mp.nstr(finite_error,12),'largest_dotted_pairing':mp.nstr(dotted,30)})
 cross=check(*roots);assert cross<-2
 parse=lambda pair:mp.mpf(str(pair[0]))+((1+mp.sqrt(5))/2)*mp.mpf(str(pair[1]))
 # K encodings are strings "a,b" or a two-item array depending on ancestry.
 def point(raw):
  out=[]
  for item in raw:
   aa,bb=N.K.parse(item);out.append(mp.mpf(aa.numerator)/aa.denominator+(1+mp.sqrt(5))/2*(mp.mpf(bb.numerator)/bb.denominator))
  return mp.matrix(out)
 escape=point(d['escape']);minimum=min((r.T*escape)[0]for r in old+roots);assert minimum>0
 return {'status':'PASS_SEPARATE_110_DIGIT_LINEAR_RECONSTRUCTION',
  'sources':{p:V.sha(p)for p in(N.PREFIX+'_caps.json',N.PREFIX+'_verify.json')},
  'precision_decimal_digits':110,'new_normals':rows,'cross_pairing':mp.nstr(cross,30),
  'remaining_escape_minimum_pairing':mp.nstr(minimum,30),
  'scope':'Independent coordinate linear-system/quadratic reconstruction and numerical high-precision pair checks, supplementing the separately replayed exact interval certificate. Prior102 normals are shared through their rigorous reconstructions.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=run()
 if args.verify:assert out==V.read(OUT)
 else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
