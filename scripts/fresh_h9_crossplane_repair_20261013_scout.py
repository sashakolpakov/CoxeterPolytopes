#!/usr/bin/env python3
"""Two-angle cap proposals in lower-rank remote normal spaces."""
import argparse
from collections import Counter
from itertools import combinations,product
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.linalg import null_space
from scripts import fresh_h9_remote_next_20261013_bridge_verify as B
V=B.V
PREFIX='data/enumerations/fresh_h9_crossplane_repair_20261013'
INCIDENT=['seed1','seed2','seed4','seed6','seed8','wall16','wall17','wall35','nextH4']

def run():
    source=B.PREFIX+'_caps.json';data=V.read(source)
    G,I,h,g,old,names,caps,center,_=B.reconstruct(data);old+=caps;names+=['bridge0','bridge1']
    f=lambda p:sum(p)/(2*V.SCALE)
    C=np.array([[f(p)for p in row]for row in old]);GG=np.array([[f(V.golden(p))for p in row]for row in G])
    hh=np.array([f(V.golden(p))for p in h]);y=np.array([f(V.golden(V.K.parse(p)))for p in data['escape']])
    cov=C@GG;active=[names.index(n)for n in INCIDENT]
    labels=(2,3,4,5,6,8,10,12);profiles=list(product(labels,repeat=2));rhs=-np.array([[2*np.cos(np.pi/m)for m in p]for p in profiles])
    gate=np.array([2*np.cos(np.pi/m)for m in(2,3,4,5,6,7,8,9,10,12,15,16,20,24,30)])
    counts=Counter();accepted=[];allrows=[];seen={};best=None
    for block in combinations(active,7):
        N=null_space(cov[list(block)]);assert N.shape==(10,3)
        H=N.T@GG@N;ev,U=np.linalg.eigh(H);row={'right_walls':[names[i]for i in block],'complement_inertia_numeric':[int((ev>1e-8).sum()),int((ev< -1e-8).sum())],'accepted':0}
        allrows.append(row)
        if ev.min()<1e-8:continue
        N=N@U@np.diag(1/np.sqrt(ev));P=cov@N;hy=hh@N;ty=y@N
        choices=[i for i in range(len(C))if i not in block]
        for a,b in combinations(choices,2):
            A=P[[a,b]];M=A@A.T;det=np.linalg.det(M)
            if det<1e-12*max(1,M[0,0]*M[1,1]):continue
            pinv=np.linalg.solve(M,A);u=rhs@pinv;remainder=2-np.einsum('bi,bi->b',u,u)
            js=np.flatnonzero(remainder>1e-9)
            if not len(js):continue
            v=np.cross(A[0],A[1]);v/=np.linalg.norm(v)
            for branch in(-1,1):
                Z=u[js]+branch*np.sqrt(remainder[js,None])*v
                counts['norm_branches']+=len(js)
                good=(Z@hy>1e-7)&(-Z@ty>1e-7);ks=js[good];Z=Z[good]
                counts['future_target_cutters']+=len(Z)
                if not len(Z):continue
                pair=Z@P.T;acute=pair.max(axis=1)<2e-6;Z=Z[acute];ks=ks[acute];pair=pair[acute]
                counts['acute_target_cutters']+=len(Z)
                if not len(Z):continue
                err=np.abs(pair[:,:,None]+gate[None,None,:]).min(axis=2);err=np.where(pair<=-2+2e-6,0,err)
                maximum=err.max(axis=1)
                j=int(np.argmin(maximum))
                if best is None or maximum[j]<best['maximum_angle_error']:
                    best={'right_walls':row['right_walls'],'anchors':[names[a],names[b]],'labels':list(profiles[ks[j]]),'branch':branch,
                          'maximum_angle_error':float(maximum[j]),'bad_walls':[names[i]for i in np.flatnonzero(err[j]>2e-6)],'coefficients':list(map(float,N@Z[j]))}
                for j in np.flatnonzero(maximum<2e-6):
                    r=N@Z[j];key=tuple(np.round(r,7));counts['compatible_profiles']+=1
                    if key in seen:continue
                    seen[key]=len(accepted);unprescribed=[i for i in range(len(C))if i not in block and i not in(a,b)]
                    other_finite=[names[i]for i in unprescribed if pair[j,i]>-2+2e-6]
                    accepted.append({'right_walls':row['right_walls'],'anchors':[names[a],names[b]],'labels':list(profiles[ks[j]]),
                        'branch':branch,'coefficients':list(map(float,r)),'target_cut':float(-r@y),'strict_center':float(r@hh),
                        'other_finite_relations':other_finite,'all_other_pairs_dotted':not other_finite,
                        'pairings':list(map(float,pair[j]))});row['accepted']+=1
        print({'block':len(allrows),'accepted':len(accepted),'counts':dict(counts),'best_error':None if best is None else best['maximum_angle_error']},flush=True)
    out={'status':'NUMERICAL_RANK7_REMOTE_TWO_ANGLE_CAP_SCOUT','source':source,'source_sha256':V.sha(source),
         'labels':list(labels),'blocks':allrows,'counts':dict(counts),'candidates':accepted,'best_failed':best,
         'scope':'All36 seven-subsets of the stated nine-wall remote set whose complement is numerically positive3D, all retained anchor pairs, stated finite labels and both norm branches. Numerical proposals only. Exact checking and volume gate are required for positive claims.'}
    (V.ROOT/(PREFIX+'_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(accepted),'candidates',flush=True)

if __name__=='__main__':run()
