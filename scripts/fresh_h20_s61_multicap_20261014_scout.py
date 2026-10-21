#!/usr/bin/env python3
"""Bounded symmetric two-cap S61 pentagon scout, with high-precision inputs."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import mpmath as mp
import fresh_h20_coupled_roles_20261014 as C
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_s61_multicap_20261014_scout'
FILES=C.INPUTS+['results/fresh_h20_coupled_roles_20261014_exact.json','scripts/fresh_h20_coupled_roles_20261014.py']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load():
    N,Q,G,H,E,dot=C.load();a,b,edge,polar=[json.loads((ROOT/p).read_text())for p in C.INPUTS];cur=json.loads((ROOT/FILES[-2]).read_text())
    for p,h in cur['sources'].items():assert sha(p)==h,p
    def f(q):q=F(q);return mp.mpf(q.numerator)/q.denominator
    rad=[mp.sqrt(d)for d in a['field']['base_square_root_basis']]
    def v0(x):return sum(f(q)*r for q,r in zip(x,rad))
    D=mp.sqrt(v0(a['field']['discriminant']))
    def v1(x):return v0(x[0])+D*v0(x[1])
    K=mp.sqrt(v1(b['field']['new_positive_sqrt_radical']))
    def v2(x):return v1(x[0])+K*v1(x[1])
    R=mp.matrix([v1(x)for x in cur['cap_R']['normal']]);R*=mp.sqrt(2/dot(R,R));N.append(R)
    bs=[mp.matrix([v1(x)for x in cur['hex_basis'][t]])for t in['C','U','V']];qs=[v1(x)for x in cur['hex_basis']['squared_lengths']]
    for co in cur['hexagon']['frame_coefficients']:N.append(sum((v0(u)*v/mp.sqrt(q)for u,v,q in zip(co,bs,qs)),mp.zeros(21,1)))
    rec=edge['recipe'];A=v2(rec['A']);w=v2(rec['W_norm'])
    X=mp.matrix([v2(x)for x in rec['X']]);W=mp.matrix([v2(x)for x in rec['W']]);frame=[X/mp.sqrt(A),N[95]/mp.sqrt(2),W/mp.sqrt(w)]
    assert max(abs(dot(u,v)-(-1 if i==j==0 else 1 if i==j else 0))for i,u in enumerate(frame)for j,v in enumerate(frame))<mp.mpf('1e-65')
    P=mp.matrix([[dot(u,n)for u in frame]for n in N]);hh=mp.matrix([dot(u,H)for u in frame])
    return N,Q,H,dot,w,frame,P,hh,cur

def relation(g):
    tol=mp.mpf('1e-55');co=-g/2
    if abs(co)<tol:return 'right'
    if co<0:return 'obtuse'
    if co>1+tol:return 'dotted'
    if abs(co-1)<tol:return 'parallel'
    m=mp.pi/mp.acos(co);j=int(mp.nint(m))
    if j>=3 and abs(m-j)<tol:return f'pi/{j}'
    return 'forbidden_angle'
def run():
    N,Q,H,dot,w,frame,P,hh,cur=load();k=mp.sqrt(w-mp.mpf(5)/6);B=2*w-mp.mpf(8)/3
    labels=list(range(2,31))+[40,60,120,1000]
    rows=[];accepted=[];count={}
    for m in[2,3,4,5]:
      cm=mp.cos(mp.pi/m)
      for ell in labels:
       c=mp.sqrt(1+mp.cos(mp.pi/ell))
       for arrangement in[-1,1]:
        t=arrangement*c*mp.sqrt(w)+2*cm
        for branch in[-1,1]:
         aa=(2*k*t+branch*mp.sqrt(2*t*t+B*(2-c*c)))/B;bb=mp.sqrt(2)*(t-aa*k)
         coeff=[mp.matrix([aa,bb,c]),mp.matrix([aa,bb,-c])]
         pairs=[P*u for u in coeff]
         assert abs(-aa*aa+bb*bb+c*c-2)<mp.mpf('1e-65')
         assert abs(-aa*aa+bb*bb-c*c+2*mp.cos(mp.pi/ell))<mp.mpf('1e-65')
         contacts=[96 if arrangement==1 else 97,97 if arrangement==1 else 96]
         assert all(abs(p[j]+2*cm)<mp.mpf('1e-65')for p,j in zip(pairs,contacts))
         bad=[[[i,relation(g),mp.nstr(g,30)]for i,g in enumerate(p)if relation(g)in['obtuse','forbidden_angle']]for p in pairs]
         good=not any(bad);key='accepted'if good else'bad_angles';count[key]=count.get(key,0)+1
         row=dict(side_label=m,mutual_label=ell,arrangement=arrangement,branch=branch,frame_coefficients=[mp.nstr(x,50)for x in[aa,bb,c]],bad_pairs=bad,common_center_signs=[mp.nstr((hh.T*u)[0],30)for u in coeff],accepted=good)
         rows.append(row)
         if good:accepted.append(len(rows)-1)
    paths=FILES+[f'scripts/{TAG}.py']
    return dict(status='BOUNDED_NUMERICAL_SYMMETRIC_PENTAGON_SCOUT',sources={p:sha(p)for p in paths},precision_decimal_digits=mp.mp.dps,side_labels=[2,3,4,5],mutual_labels=labels,trials=len(rows),counts=count,accepted_rows=accepted,rows=rows,scope='Specified symmetric ansatz only; all107 pair tests at100 decimal digits. No exact exclusion, all-real conclusion, or global finite-volume claim. Actual endpoint/compactness checks are required before accepting a construction.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['trials'],d['counts'],d['accepted_rows'],flush=True)
