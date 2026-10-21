#!/usr/bin/env python3
"""Two-vertex finite-edge insertion rather than an isolated new corner."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts.fresh_h9_compact_next_20261013_scout import P,V,K,PREFIX,BLOCKS

def run():
    d=V.read(P.PREFIX+'_caps.json');G,I,h,g,old,names,caps,center,_=P.reconstruct(d);ff=lambda x:sum(x)/(2*V.SCALE)
    C=np.array([[ff(x)for x in r]for r in old+caps]);names+=P.NEW;Q=np.array([[ff(V.golden(x))for x in row]for row in G]);hh=np.array([ff(V.golden(x))for x in h]);cc=np.array([ff(x)for x in center]);yy=np.array([ff(V.golden(K.parse(x)))for x in d['escape']]);w=-np.linalg.solve(Q,yy)
    candidates=[];stats=[];ks=np.geomspace(.05,1e4,600);angles=np.arange(180)*2*np.pi/540
    for ss in BLOCKS:
        S=C[[names.index(n)for n in ss]];Mi=np.linalg.inv(S@Q@S.T)
        def ps(v):return v-(v@Q@S.T)@Mi@S
        for anchor in names:
            if anchor in ss:continue
            a=ps(C[names.index(anchor)]);q=a@Q@a
            if q<=1e-7:continue
            def proj(v):
                z=ps(v);return z-(z@Q@a)/q*a
            t=proj(cc);t/=np.sqrt(-t@Q@t)
            e=proj(w);e+=(e@Q@t)*t;en=e@Q@e
            if en<1e-10:continue
            e/=np.sqrt(en);F=[t,e];sources=[]
            for j in range(10):
                z=proj(np.eye(10)[j])
                for b in F:z-=(z@Q@b)/(b@Q@b)*b
                qn=z@Q@z
                if qn>1e-8:F.append(z/np.sqrt(qn));sources.append(j)
                if len(F)==4:break
            assert len(F)==4
            for m in(3,4,5,6,8,10,12,16):
                if anchor=='seed4'and m==6:continue # Forced positive affine 3-6 path.
                b=2*np.cos(np.pi/m);nu=b*b/q-2
                if nu<=0:continue
                row={'block':ss,'anchor':anchor,'projected_anchor_norm':float(q),'label':m,'boosts':len(ks),'base_compatible':0,'target_cut':0,'full_compatible':0};extra=[i for i,n in enumerate(names)if n not in ss+[anchor]]
                for k in ks:
                    ct=(k+1/k)/2;st=(k-1/k)/2;p=ct*t+st*e;r=-b/q*a+np.sqrt(nu)*p
                    pp=r@Q@C.T
                    if (r@hh)<=0 or np.max(pp[extra])>=-2-1e-8:continue
                    row['base_compatible']+=1
                    if -r@yy<=1e-7:continue
                    row['target_cut']+=1;B=st*t+ct*e
                    if B@hh<0:B=-B
                    for theta in angles:
                        aa=theta+np.arange(3)*2*np.pi/3;R=np.sqrt(2/3)*B+np.sqrt(4/3)*(np.cos(aa)[:,None]*F[2]+np.sin(aa)[:,None]*F[3])
                        if np.min(R@hh)<=0 or np.max((R@Q@C.T)[:,extra])>=-2-1e-7:continue
                        row['full_compatible']+=1
                        candidates.append({'orthogonal_old_walls':ss,'anchor_old':anchor,'anchor_label':m,'frame_sources':sources,'boost':float(k),'triangle_angle':float(theta),'normal_r':r.tolist(),'other_normals':R.tolist(),'target_cut':float(-r@yy)})
                        break
                stats.append(row);print(row,flush=True)
    out={'status':'NUMERICAL_TWO_VERTEX_FINITE_EDGE_SCOUT','source':P.PREFIX+'_caps.json','source_sha256':V.sha(P.PREFIX+'_caps.json'),'stats':stats,'candidates':candidates,
         'scope':'Two spherical rank5 blocks; all old spherical one-wall extensions; labels3,4,5,6,8,10,12,16 with positive timelike remainder;600 boosts and180 triangle phases. The forced affine seed4/label6 profile is omitted. Only r needs cut the saved gap. No exhaustive family claim.'}
    (V.ROOT/(PREFIX+'_edge_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(candidates),flush=True)

if __name__=='__main__':run()
