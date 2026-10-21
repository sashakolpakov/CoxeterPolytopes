#!/usr/bin/env python3
"""Bounded numerical search for mutually compatible endpoint covers."""
import argparse,json,time
import numpy as np
from fresh_h9_makarov_joint_20261013_selection import setup,ROOT,PREFIX

def bits(mask):return int.from_bytes(np.packbits(mask,bitorder='little').tobytes(),'little')
def run(seconds=120,simple=False):
    start=time.time();Q,W,c,R0,refs0,null,edges=setup();cv0=R0@Q@null.T>1e-9
    keep=cv0.any(axis=1)
    if simple:
        banks={s:json.loads((ROOT/f'results/{PREFIX}{s}.json').read_text())['caps']for s in ['', '_profiles']}
        for k,(s,i)in enumerate(refs0):
            row=banks[s][i]
            if not s:
                contacts=list(zip(row['old_support'],row['old_labels']))
                labels=row['old_labels']+row['transverse_labels']
            else:
                contacts=row['contacts']+[[row['anchor'],row['anchor_label']]]
                labels=[m for _,m in contacts]
            if any(m not in [2,3,5]for m in labels)or any(w==6 and m!=2 for w,m in contacts):keep[k]=False
    ids=np.flatnonzero(keep);R=R0[ids];cov=cv0[ids];N=len(ids);RQ=R@Q
    edge_masks=[bits(cov[:,j])for j in range(40)];cap_masks=[sum(1<<j for j in np.flatnonzero(row))for row in cov]
    cache={};nodes=0;maxcovered=0;best=[];leaves=0
    def compatible(k):
        if k not in cache:
            ps=RQ@R[k];ok=ps<-2+3e-7
            for m in [2,3,4,5,6,7,8,9,10,12,15,16,17,20,24,30]:ok|=abs(ps+2*np.cos(np.pi/m))<3e-7
            ok[k]=False;cache[k]=bits(ok)
        return cache[k]
    def search(avail,missing,chosen):
        nonlocal nodes,maxcovered,best,leaves
        nodes+=1
        if nodes%100==0 and time.time()-start>seconds:raise TimeoutError
        covered=40-missing.bit_count()
        if covered>maxcovered:
            maxcovered=covered;best=chosen;print('DFS',nodes,'covered',covered,'time',time.time()-start,flush=True)
        if not missing:return chosen
        options=None
        for j in range(40):
            if missing>>j&1:
                a=avail&edge_masks[j]
                if not a:leaves+=1;return None
                if options is None or a.bit_count()<options.bit_count():options=a
        ks=[]
        while options:
            one=options&-options;k=one.bit_length()-1;options-=one;ks.append(k)
        ks.sort(key=lambda k:-(cap_masks[k]&missing).bit_count())
        for k in ks:
            out=search(avail&compatible(k),missing&~cap_masks[k],chosen+[k])
            if out is not None:return out
        return None
    status='NUMERICAL_EXHAUSTED';ans=None
    try:
        ans=search((1<<N)-1,(1<<40)-1,[])
        if ans is not None:status='NUMERICAL_COMPATIBLE_40_COVER'
    except TimeoutError:status='TIME_LIMIT'
    chosen=ans if ans is not None else best
    out={'status':status,'numeric_only':True,'entries':N,'nodes':nodes,'dead_leaves':leaves,'cached_rows':len(cache),
         'max_covered':maxcovered,'selected':[int(ids[k])for k in chosen], 'references':[refs0[ids[k]]for k in chosen],
         'seconds':time.time()-start}
    out['K_rhs_only']=simple
    suffix='_cover_simple'if simple else '_cover_dfs'
    (ROOT/f'results/{PREFIX}{suffix}.json').write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=int,default=120);p.add_argument('--simple',action='store_true');a=p.parse_args();run(a.seconds,a.simple)
