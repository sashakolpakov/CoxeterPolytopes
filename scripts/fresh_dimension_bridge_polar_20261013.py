#!/usr/bin/env python3
"""Exact polar cuts proposed by linear-program vertices of the H20 lift."""
import argparse
import hashlib
import json
from pathlib import Path
import time
from collections import Counter
from fractions import Fraction as F

import numpy as np
from scipy.optimize import linprog
import fresh_dimension_bridge_20261013 as K

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_dimension_bridge_20261013.json'
OUT=ROOT/'data/enumerations/fresh_dimension_bridge_polar_20261013.json'


def div(a,b):
    den=b[0]*b[0]-6*b[1]*b[1];assert den
    return (a[0]*b[0]-6*a[1]*b[1])/den,(a[1]*b[0]-a[0]*b[1])/den


def nullvector(rows):
    a=[list(r)for r in rows];piv=[];level=0
    for j in range(21):
        k=next((i for i in range(level,len(a))if a[i][j]!=K.Z),None)
        if k is None:continue
        a[level],a[k]=a[k],a[level];p=a[level][j]
        a[level]=[div(x,p)for x in a[level]]
        for i in range(len(a)):
            if i==level:continue
            p=a[i][j]
            if p!=K.Z:a[i]=[K.sub(x,K.mul(p,y))for x,y in zip(a[i],a[level])]
        piv.append(j);level+=1
        if level==len(a):break
    assert level==20
    free,=set(range(21))-set(piv);v=[K.Z]*21;v[free]=K.O
    for i,j in enumerate(piv):v[j]=K.neg(a[i][free])
    assert all(K.dot(row,v)==K.Z for row in rows)
    return tuple(v)


def run(limit,attempts):
    started=time.monotonic();base=json.loads((ROOT/SOURCE).read_text())
    Q=tuple(tuple(map(K.dec,r))for r in base['form'])
    roots=list(map(K.parsev,base['normals']));h=K.parsev(base['strict_center'])
    cv=[K.mv(Q,r)for r in roots];norms=[K.dot(r,c)for r,c in zip(roots,cv)]
    qn=np.array([[K.value(x)for x in r]for r in Q]);hn=np.array([K.value(x)for x in h]);hn/=np.sqrt(-hn@qn@hn)
    rng=np.random.default_rng(200071);cuts=[];rejected=[];counts=Counter();remaining=None
    seen=set();old_escape=K.parsev(base['strict_escape']['vector'])
    def numerical():
        rn=np.array([[K.value(x)for x in r]for r in roots]);nn=np.array([K.value(x)for x in norms])
        return rn,np.sqrt(nn)
    R,scales=numerical()
    for step in range(limit):
        accepted=False
        for attempt in range(attempts):
            if step==0 and attempt==0:
                v=np.array([K.value(x)for x in old_escape]);ids=np.argsort(abs(R@qn@v)/scales)[:20]
            else:
                lp=linprog(rng.normal(size=21),A_ub=(R@qn)/scales[:,None],b_ub=np.zeros(len(R)),
                    A_eq=(hn@qn)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
                counts['LP_proposals']+=1
                if not lp.success:counts['LP_failure']+=1;continue
                v=lp.x
                if v@qn@v<1e-8:counts['nonspacelike']+=1;continue
                ids=np.argsort(abs(R@qn@v)/scales)[:20]
            ids=tuple(sorted(map(int,ids)))
            if ids in seen:continue
            seen.add(ids)
            if np.linalg.matrix_rank(R[list(ids)]@qn)<20:counts['dependent_active20']+=1;continue
            r=nullvector([cv[i]for i in ids]);q=K.pair(r,r,Q);ch=K.pair(r,h,Q)
            assert K.sign(q)>0
            if K.sign(ch)>0:r=tuple(map(K.neg,r));ch=K.neg(ch)
            assert K.sign(ch)<0
            # The exact polar vector is also the exact affine vertex ray.
            ps=[K.dot(r,c)for c in cv]
            if any(K.sign(p)>0 for p in ps):counts['wrong_active_set']+=1;continue
            kinds=[];bad=[]
            for i,(p,n)in enumerate(zip(ps,norms)):
                try:kinds.append(K.classify(p,q,n))
                except AssertionError:bad.append(i);kinds.append('forbidden')
            if bad:
                counts['forbidden_polar']+=1
                # Preserve a compact exact obstruction witness for later
                # non-right cap/profile construction, not a no-go claim.
                rejected.append({'active_indices':list(ids),'vector':K.encv(r),
                                 'norm':K.enc(q),'forbidden_indices':bad,
                                 'forbidden_pairings':[K.enc(ps[i])for i in bad]})
                remaining={'vector':K.encv(r),'norm':K.enc(q),'active_indices':list(ids)}
                continue
            target=K.pair(r,old_escape,Q)if step==0 else q
            record={'index':len(roots),'active_indices':list(ids),'vector':K.encv(r),'norm':K.enc(q),
                    'center_pairing':K.enc(ch),'pairings':list(map(K.enc,ps)),
                    'pair_kinds':kinds,'first_old_escape_pairing':K.enc(target)if step==0 else None}
            cuts.append(record);roots.append(r);cv.append(K.mv(Q,r));norms.append(q)
            counts['accepted_polar']+=1;accepted=True;R,scales=numerical()
            print({'step':step,'walls':len(roots),'new_pair_counts':dict(Counter(kinds)),
                   'seconds':round(time.monotonic()-started,2)},flush=True)
            break
        output={'status':'EXACT_POLAR_EXTENSIONS_WITH_NUMERICAL_FRONTIER_SEARCH',
                'source':SOURCE,'source_sha256':hashlib.sha256((ROOT/SOURCE).read_bytes()).hexdigest(),
                'cuts':cuts,'wall_count':len(roots),'counts':dict(counts),
                'rejected_polar_proposals':rejected,'last_escape_ray':remaining,
                'seconds':round(time.monotonic()-started,3),
                'scope':'Every accepted root and its pair relations with all previous walls is exact overQ(sqrt6). LP proposes active sets only. Random objectives do not exhaust vertices; rejected polar angles do not exclude other caps. No finite-volume claim.'}
        OUT.write_text(json.dumps(output,indent=2)+'\n')
        if not accepted:break
    # Use a surviving rejected ray, not a stale witness cut off by a later
    # accepted polar. Perturb toward the unchanged strict interior.
    for row in reversed(rejected):
        v=K.parsev(row['vector'])
        if any(K.sign(K.dot(v,c))>0 for c in cv):continue
        for power in range(1,200):
            w=tuple(K.add(x,K.sc(y,F(1,2**power)))for x,y in zip(v,h))
            q=K.pair(w,w,Q);ch=K.pair(w,h,Q);ps=[K.dot(w,c)for c in cv]
            if K.sign(q)>0 and K.sign(ch)<0 and all(K.sign(x)<0 for x in ps):break
        else:raise AssertionError('No strict perturbation')
        output['strict_escape']={'vector':K.encv(w),'norm':K.enc(q),
            'center_pairing':K.enc(ch),'all_wall_pairings':list(map(K.enc,ps)),
            'perturbation_power':power,'source_active_indices':row['active_indices']}
        break
    else:output['strict_escape']=None
    output['status']='EXACT_POLAR_EXTENSIONS_WITH_STRICT_ESCAPE'if output['strict_escape']else'EXACT_POLAR_EXTENSIONS_VOLUME_UNDECIDED'
    OUT.write_text(json.dumps(output,indent=2)+'\n')
    print({'done':len(cuts),'walls':len(roots),'counts':dict(counts)},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=50);p.add_argument('--attempts',type=int,default=300)
    args=p.parse_args();run(args.limit,args.attempts)
