#!/usr/bin/env python3
"""Couple up to six finite relations across the complete remote nine-wall span."""
from collections import Counter
from itertools import combinations,product
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts import fresh_h9_remote_next_20261013_bridge_verify as B
from scripts.fresh_h9_crossplane_repair_20261013_scout import INCIDENT,PREFIX
V=B.V

def run():
    source=B.PREFIX+'_caps.json';data=V.read(source)
    G,I,h,g,old,names,caps,center,_=B.reconstruct(data);old+=caps;names+=['bridge0','bridge1']
    f=lambda p:sum(p)/(2*V.SCALE)
    C=np.array([[f(p)for p in row]for row in old]);GG=np.array([[f(V.golden(p))for p in row]for row in G])
    hh=np.array([f(V.golden(p))for p in h]);y=np.array([f(V.golden(V.K.parse(p)))for p in data['escape']])
    cov=C@GG;ids=[names.index(n)for n in INCIDENT];A=C[ids];H=A@GG@A.T;Hi=np.linalg.inv(H)
    ev=np.linalg.eigvalsh(H);assert(sum(ev>1e-7),sum(ev< -1e-7))==(8,1)
    cc=np.array([f(x)for x in center]);v=cc-(cc@GG@A.T)@Hi@A;v/=np.sqrt(v@GG@v)
    labels=(3,4,5,6,10);values={m:2*np.cos(np.pi/m)for m in labels}
    gate=np.array([2*np.cos(np.pi/m)for m in(2,3,4,5,6,7,8,9,10,12,15,16,20,24,30)])
    counts=Counter();records=[];seen=set();best=None;subcounts=[]
    for k in range(3,7):
        profiles=list(product(labels,repeat=k));values0=np.array([[-values[m]for m in p]for p in profiles]);row={'nonright_support':k,'profiles':0,'new_caps':0}
        for support in combinations(range(9),k):
            for start in range(0,len(profiles),256):
                chunk=profiles[start:start+256];P=np.zeros((len(chunk),9));P[:,support]=values0[start:start+256];row['profiles']+=len(chunk)
                U=P@Hi@A;q=np.einsum('bi,ij,bj->b',P,Hi,P);rem=2-q;js=np.flatnonzero(rem>1e-8)
                for branch in(-1,1):
                    Z=U[js]+branch*np.sqrt(rem[js,None])*v;counts['norm_branches']+=len(Z)
                    good=(Z@hh>1e-7)&(-Z@y>1e-7);ks=js[good];Z=Z[good];counts['future_cutters']+=len(Z)
                    if not len(Z):continue
                    ps=Z@cov.T;good=ps.max(axis=1)<2e-6;ks=ks[good];Z=Z[good];ps=ps[good];counts['acute_cutters']+=len(Z)
                    if not len(Z):continue
                    err=np.abs(ps[:,:,None]+gate[None,None,:]).min(axis=2);err=np.where(ps<=-2+2e-6,0,err);maximum=err.max(axis=1)
                    j=int(np.argmin(maximum))
                    if best is None or maximum[j]<best['maximum_angle_error']:
                        best={'nonright_support':[INCIDENT[i]for i in support],'labels':list(chunk[ks[j]]),'branch':branch,'maximum_angle_error':float(maximum[j]),
                              'bad_walls':[names[i]for i in np.flatnonzero(err[j]>2e-6)],'coefficients':list(map(float,Z[j])),'pairings':list(map(float,ps[j]))}
                    for j in np.flatnonzero(maximum<2e-6):
                        key=tuple(np.round(Z[j],7))
                        if key in seen:continue
                        seen.add(key);finite=[names[i]for i in range(88)if i not in ids and ps[j,i]>-2+2e-6]
                        records.append({'nonright_support':[INCIDENT[i]for i in support],'labels':list(chunk[ks[j]]),'branch':branch,
                            'coefficients':list(map(float,Z[j])),'pairings':list(map(float,ps[j])),'target_cut':float(-Z[j]@y),'extra_finite_walls':finite});row['new_caps']+=1
        subcounts.append(row);print({'done_support_size':k,'candidates':len(records),'counts':dict(counts),'best_error':None if best is None else best['maximum_angle_error']},flush=True)
    out={'status':'NUMERICAL_FULL_REMOTE_FACE_METRIC_SCOUT','source':source,'source_sha256':V.sha(source),'incident_walls':INCIDENT,
         'labels':list(labels),'support_sizes':[3,4,5,6],'subcounts':subcounts,'counts':dict(counts),'candidates':records,'best_failed':best,
         'scope':'Finite metric profiles on the complete Lorentzian9-wall span; both norm branches. Stated finite support/alphabet only, all88-wall numerical Coxeter gate. No exact or exhaustive obstruction claim.'}
    (V.ROOT/(PREFIX+'_face_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(records),flush=True)

if __name__=='__main__':run()
