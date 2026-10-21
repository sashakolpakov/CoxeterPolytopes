#!/usr/bin/env python3
"""Numerical second positive-plane proposal search at the saved 88-wall gap."""
import json
from itertools import combinations
import numpy as np
from scripts import fresh_h9_remote_next_20261013_bridge_verify as B
V=B.V
OUT=B.PREFIX+'_scout.json'

def run():
    source=B.PREFIX+'_caps.json';data=V.read(source)
    G,I,h,g,roots,names,caps,center,_=B.reconstruct(data);roots+=caps;names+=['bridge0','bridge1']
    f=lambda x:sum(x)/(2*V.SCALE)
    C=np.array([[f(x)for x in r]for r in roots]);GG=np.array([[f(V.golden(x))for x in r]for r in G])
    hh=np.array([f(V.golden(x))for x in h]);cc=np.array([f(x)for x in center]);y=np.array([f(V.golden(V.K.parse(x)))for x in data['escape']])
    active=[i for i,x in enumerate(C@y)if x<.002];assert len(active)==9
    labels={2:0.,3:1.,4:np.sqrt(2),5:(1+np.sqrt(5))/2,6:np.sqrt(3),10:np.sqrt((5+np.sqrt(5))/2)}
    gate=np.array([2*np.cos(np.pi/m)for m in(2,3,4,5,6,8,10,12,15,16,20,24,30)])
    rows=[];accepted=[];singles=[]
    for block in combinations(active,8):
        A=C[list(block)];M=A@GG@A.T;ev=np.linalg.eigvalsh(M);sig=[int((ev>1e-7).sum()),int((ev< -1e-7).sum())]
        row={'walls':[names[i]for i in block],'signature_numeric':sig,'branches':0,'target_single_caps':0,'compatible_single_caps':0,'compatible_pairs':0}
        rows.append(row)
        if sig!=[7,1]:continue
        Mi=np.linalg.inv(M)
        def proj(w):return w-(w@GG@A.T)@Mi@A
        for a in range(len(C)):
            if a in block:continue
            u=proj(C[a]);d=u@GG@u
            if d<1e-8:continue
            v=proj(cc);v-=((v@GG@u)/d)*u;e=v@GG@v
            if e<1e-8:continue
            for m,k in labels.items():
                q=2-k*k/d
                if q<1e-8:continue
                for branch in(-1,1):
                    row['branches']+=1
                    r0=-k*u/d+branch*np.sqrt(q/e)*v
                    r1=-np.sqrt(q/d)*u-branch*k*v/np.sqrt(d*e)
                    # Both quarter turns are tested; the strict center fixes orientation.
                    for sign in(-1,1):
                        R=np.array([r0,sign*r1]);P=R@GG@C.T
                        oriented=R@hh>1e-7;cut=-R@y>1e-7
                        compatible=[]
                        for j in range(2):
                            finite=P[j]>-2+1e-7
                            good=oriented[j]and P[j].max()<1e-6 and (not finite.any()or np.abs(P[j,finite,None]+gate).min(axis=1).max()<2e-6)
                            compatible.append(bool(good))
                        row['target_single_caps']+=int(sum(oriented&cut))
                        for j in range(2):
                            if compatible[j]and cut[j]:
                                row['compatible_single_caps']+=1
                                singles.append({'block':row['walls'],'anchor':names[a],'label':m,'branch':branch,'second_sign':sign,'which':j,'cut':float(-R[j]@y),'coefficients':R[j].tolist()})
                        if all(compatible)and any(cut):
                            row['compatible_pairs']+=1
                            accepted.append({'block':row['walls'],'anchor':names[a],'label':m,'branch':branch,'second_sign':sign,'cuts':list(map(float,-R@y)),'coefficients':R.tolist()})
        print(row,flush=True)
    out={'status':'BOUNDED_NUMERICAL_POSITIVE_PLANE_SCOUT','source':source,'source_sha256':V.sha(source),'near_active_walls':[names[i]for i in active],
         'blocks':rows,'pairs':accepted,'single_caps':singles,'scope':'Only the nine 8-subsets of the saved numeric incident 9-set, all existing anchors, labels2,3,4,5,6,10 and both quarter-turn signs. Numerical proposals only; absence does not exclude other incident spans or angles.'}
    (V.ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print('pairs',len(accepted),'singles',len(singles),flush=True)

if __name__=='__main__':run()
