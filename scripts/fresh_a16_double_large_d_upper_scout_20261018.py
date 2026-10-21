#!/usr/bin/env python3
"""Small complete increasing-family upper for the proposed coupled A16 rule.

Arithmetic scout only: the geometric transfer/retention theorem is separate.
No giant ambient certificate is read and no existing evidence is rewritten.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import fcntl
import json
from pathlib import Path
import signal
import sys

import fresh_d17_node3_large_d_guard_transfer_20261018 as V
import fresh_d17_node3_large_d_guard_upper_20261018 as U
import run_lean_batches_20261018 as R

ROOT=V.ROOT;C=U.C;N=57
TAG='fresh_a16_double_large_d_upper_scout_20261018'


def delta(fs):
    co=Counter(fs)
    if not V.has_guard(fs):return 0
    out=-co['A',16]*(192+(80 if co['A',1] else 0)) if co['A',16]>=2 else 0
    if co['A',16]:
        out+=64*co['A',1]*co['A',13]
        if co['A',2]:out+=64*co['A',1]*co['A',13]+64*co['A',1]*co['D',13]
        if co['A',1]:out+=160*co['A',4]*co['A',12]+128*co['A',2]*co['A',14]
        out+=128*co['A',4]*co['A',13]+128*co['A',4]*co['D',13]+64*co['A',15]
    return out


def value(fs):return V.value(fs)+delta(fs)


def run(saved=None):
    prior={(tuple(r['counts']),r['k5'],r['a12']):r for r in saved['classes']}if saved else{}
    classes=[];checks=preds=options=qoptions=reserves=0
    for cs,used,key in C.fixed_domain():
        U.memory_guard()
        if not cs[7] or not(key[9]or cs[4]or cs[5]or cs[6]):continue
        for k5 in range((N-used)//5+1):
            rem=N-used-5*k5
            if rem<18:continue
            old=prior.get((cs,k5,key[9]))
            if saved:assert old is not None
            fixed=C.fixed_components(cs,key)+[('A',5)]*k5
            fc=C.upper_value(fixed)+64*cs[6]
            weights={r:C.item(r,key,k5)+(128*cs[4]if r==1 else 256*cs[5]if r==2 else 256*cs[4]+160*key[9]if r==4 else 0)for r in C.ITEMS}
            tabs=[]
            for no2 in(False,True):
                ws={r:w for r,w in weights.items()if not no2 or r!=2}
                tab,a,b=U.table(ws,old['tables'][int(no2)]if old else None)
                tab['exclude2']=no2;tabs.append(tab);checks+=a;preds+=b
            rows=[]
            for s in range(18,rem+1):
                n=rem-s
                assert C.item(s,key,k5)==C.actual_base(('D',s))
                base=fc+C.actual_base(('D',s))
                opts=[tabs[0]['values'][n]]
                for a1 in(1,2):
                    og=a1==2 and bool(cs[1]);ng=bool(cs[0]and k5)
                    if n<a1 or not(og or ng)or(a1==2 and not og):continue
                    for a2 in range((n-a1)//2+1):
                        v=tabs[1]['values'][n-a1-2*a2]
                        if v is not None:opts.append(a1*weights[1]+a2*weights[2]+v+(36*a2*(a2-1)if og else 0)+(64*a2*k5 if ng else 0))
                opts=[v for v in opts if v is not None];options+=len(opts)
                delta17=16*k5*key[9]+(8*k5*n if cs[2]else 0)
                ordinary=base+max(opts)+delta17
                quarter=None
                if cs[1]and n>=7:
                    for a1 in range(1,n-5):
                        for a4 in range(1,(n-a1-2)//4+1):
                            for a2 in range(1,(n-a1-4*a4)//2+1):
                                left=n-a1-4*a4-2*a2;v=tabs[1]['values'][left]
                                if v is None:continue
                                z=base+a1*weights[1]+a2*weights[2]+a4*weights[4]+v+delta17+(36*a2*(a2-1)if a1>=2 else 0)+(64*a2*k5 if cs[0]else 0)+64*a2*a4
                                quarter=z if quarter is None else max(quarter,z);qoptions+=1
                upper=max(ordinary,quarter if quarter is not None else -1)
                rows.append(dict(reserved_D_rank=s,residual_rank=n,ordinary_upper=F(ordinary,64),quarter_upper=None if quarter is None else F(quarter,64),order17_upper=F(delta17,64),upper=F(upper,64)))
                reserves+=1
            classes.append(dict(counts=cs,k5=k5,a12=key[9],key=key,remaining_rank=rem,tables=tabs,reserved_rows=rows,upper=max(r['upper']for r in rows)))
    if saved:assert len(prior)==len(classes)
    witnesses=[]
    for fs in([('A',1),('A',16),('A',16),('D',24)],
        [('A',16),('A',16),('D',25)],[('D',16),('D',16),('D',25)]):
        assert sum(r for k,r in fs)==57
        witnesses.append(dict(forest=fs,old=F(V.value(fs),64),proposed=F(value(fs),64)))
    return C.enc(dict(status='PASS_A16_DOUBLE_LARGE_D_UPPER_SCOUT',
        sources={str(p.relative_to(ROOT)):U.digest(p)for p in(Path(__file__),Path(V.__file__),Path(U.__file__),Path(C.__file__),Path(C.C.__file__),Path(C.P.__file__))},
        scale=64,complete_classes=len(classes),reserved_D_options=reserves,
        Bellman_inequalities=checks,attaining_predecessors=preds,ordinary_options=options,quarter_options=qoptions,
        classes=classes,increasing_family_upper=max(r['upper']for r in classes),witnesses=witnesses,
        geometric_theorem_proved=False,new_dimension_bound=None,
        scope='Complete proposed-increase family with actual A16 relaxed to fixedrank16 and reserved D>=18; current D17/oldorder17 positive costs included, both negative transfers ignored in the upper.'))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);R.resources();U.memory_guard()
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    saved=json.loads(path.read_text())if args.verify else None
    out=run(saved);U.memory_guard();R.resources(start=False)
    if args.verify:assert saved==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_classes'],out['increasing_family_upper'],[(r['old'],r['proposed'])for r in out['witnesses']])
