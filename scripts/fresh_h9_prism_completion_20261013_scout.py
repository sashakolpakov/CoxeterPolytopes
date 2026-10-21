#!/usr/bin/env python3
"""Numerical diagnostics for new couplings at the surviving 96-wall gap."""
import json,sys,itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts import fresh_h9_face_bridge_20261013_verify as B
V=B.V;K=B.K
PREFIX='data/enumerations/fresh_h9_prism_completion_20261013'

def setup():
    path=V.ROOT/(PREFIX+'_snapshot.json')
    if path.exists():
        d=json.loads(path.read_text());assert d['source_sha256']==V.sha(B.PREFIX+'_caps.json')
        return d
    data=V.read(B.PREFIX+'_caps.json')
    G,I,h,g,old,names,caps,center,_,_=B.reconstruct(data)
    ff=lambda x:sum(x)/(2*V.SCALE)
    d={'source':B.PREFIX+'_caps.json','source_sha256':V.sha(B.PREFIX+'_caps.json'),
       'gram':[[ff(V.golden(x))for x in row]for row in G],
       'exact_gram':[[K.enc(x)for x in row]for row in G],
       'normals':[[ff(x)for x in r]for r in old+caps], 'names':names+B.NEW,
       'center':[ff(x)for x in center], 'h':[ff(V.golden(x))for x in h],
       'target':[ff(V.golden(K.parse(x)))for x in data['escape']]}
    path.write_text(json.dumps(d,indent=2)+'\n');return d

def run():
    d=setup();G=np.array(d['gram']);R=np.array(d['normals']);h=np.array(d['h']);y=np.array(d['target']);names=d['names']
    cov=R@G;perm=[0,2,1,3,4,6,5,7,8,9]
    print('GRAM',G.tolist(),flush=True)
    print('SWAP_ISOMETRY',np.max(abs(G-G[np.ix_(perm,perm)])),flush=True)
    for i in range(79,len(R)):
        rr=R[i,perm];p=cov@rr;q=rr@G@rr
        valid=(p<-2-1e-7)|np.logical_or.reduce([abs(p+2*np.cos(np.pi/m))<1e-7 for m in(2,3,4,5,6,8,10,12)])
        print('MIRROR',names[i],'norm',q,'cut',-rr@y,'conflicts',[names[j]for j in np.flatnonzero(~valid)],flush=True)
    return d

if __name__=='__main__':run()
