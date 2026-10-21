#!/usr/bin/env python3
"""Complete new finite maximum from three old top classes and four covers."""
import argparse
from collections import Counter
from fractions import Fraction as F
import fcntl
from itertools import combinations_with_replacement
import json
from pathlib import Path
import signal
import fresh_d16_node3_rank16_large_d_guard_scout_20261018 as N

S=N.S;M=N.M;V=N.V;U=N.U;C=N.C;ROOT=N.ROOT
TAG='fresh_d16_node3_guard_complete_finite_20261018'


def value(fs):
    co=Counter(fs);g=V.has_guard(fs);r16=co['A',16]+co['D',16]
    delta=-64*co['D',16]if g and r16>=2 else 0
    if g and r16 and co['A',1]:delta+=256*co['A',2]*co['A',13]
    return M.value(fs)+delta


def selected_array(path,marker):
    with path.open('rb')as stream:
        buf=b''
        while True:
            b=stream.read(65536)
            if not b:raise ValueError('Missing compact class summary')
            buf+=b;i=buf.find(marker)
            if i<0:buf=buf[-len(marker):];continue
            buf=buf[i+len(marker)-1:]
            while True:
                try:return json.JSONDecoder().raw_decode(buf.decode())[0]
                except json.JSONDecodeError:
                    buf+=stream.read(65536);assert len(buf)<2*1024**2


def top_classes():
    oldpath=ROOT/f'results/{C.TAG}.json'
    old=selected_array(oldpath,b'"above_curvature_classes": [')
    assert len(old)==38
    high=[r for r in old if F(r['maximum'])>F(1487,4)]
    expected={(tuple([0]*7+[k,0]),0,0)for k in(0,1,2)}
    assert{(tuple(r['counts']),r['k5'],r['a12'])for r in high}==expected
    assert all(F(r['maximum'])==372 for r in high)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS}
    residual=[t for t in C.TYPES if C.representative(t)in reps and t!=('A',1)]
    fibers=[t for t in C.TYPES if C.representative(t)==C.representative(('D',16))]
    assert all(C.actual_cross(a,b)==0 for a in residual+[('A',1)]for b in residual+[('A',1)])
    out=[];checks=preds=options=0
    for k in(0,1,2):
        rows=[]
        for fixed in combinations_with_replacement(fibers,k):
            rem=57-16*k;tables=[];candidates=[]
            for a1 in range(rem+1):
                anchor=list(fixed)+[('A',1)]*a1;left=rem-a1
                ts=[t for t in residual if t[1]<=left and not(t[0]=='D'and t[1]>=18)]
                base=value(anchor);ws=[value(anchor+[t])-base for t in ts]
                tab,a,b=C.additive(ts,ws,left);checks+=a;preds+=b;tables.append(dict(a1=a1,guard=False,table=tab))
                if tab['values'][left]is not None:candidates.append((base+tab['values'][left],a1,False,None,left,len(tables)-1));options+=1
                if left<18:continue
                anchorg=anchor+[('D',18)];cap=left-18;gb=value(anchorg)
                ts=[t for t in residual if t[1]<=cap];ws=[value(anchorg+[t])-gb for t in ts]
                tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b;tables.append(dict(a1=a1,guard=True,table=tab))
                for s in range(18,left+1):
                    n=left-s;v=tab['values'][n]
                    if v is None:continue
                    chosen=anchor+[('D',s)];bc=value(chosen)
                    for t,w in zip(ts,ws):
                        if t[1]<=n:assert value(chosen+[t])-bc==w
                    candidates.append((bc+v,a1,True,s,n,len(tables)-1));options+=1
            best=max(candidates);z,a1,guard,s,left,ti=best;forest=list(fixed)+[('A',1)]*a1+([('D',s)]if guard else[])
            tab=tables[ti]['table']
            while left:
                t=tab['types'][tab['last'][left]];forest.append(t);left-=t[1]
            assert sum(t[1]for t in forest)==57 and value(forest)==z
            rows.append(dict(fixed=fixed,tables=tables,maximum=F(z,64),attaining_forest=forest,winning_option=best[1:]))
        out.append(dict(fixed_rank16_count=k,rows=rows,maximum=max(r['maximum']for r in rows)))
    assert max(r['maximum']for r in out)<=F(1487,4)
    return dict(inherited_compact_class_maxima=old,inherited_high_classes=high,classes=out,
        Bellman_inequalities=checks,attaining_predecessors=preds,explicit_options=options,
        three_class_maximum=max(r['maximum']for r in out),all_other_original_classes_at_most='1487/4')


def run():
    top=top_classes();scout=V.load(N.TAG);order17=V.load(V.O.TAG);d17=V.load(U.TAG);a16=V.load(M.TAG)
    covers=[F(order17['finite']['increasing_family_upper']),F(d17['increasing_family_upper']),F(a16['increasing_family_upper']),F(scout['upper']['increasing_family_upper'])]
    assert covers==[F(739,2),F(365),F(371),F(11513,32)]and max(covers)<F(1487,4)
    witness=[('A',17),('D',16),('D',24)];assert value(witness)==64*F(1487,4)
    paths={Path(__file__),Path(N.__file__),Path(M.__file__),Path(V.__file__),Path(C.__file__),
        ROOT/f'results/{C.TAG}.json',ROOT/f'results/{N.TAG}.json',ROOT/f'results/{V.O.TAG}.json',ROOT/f'results/{U.TAG}.json',ROOT/f'results/{M.TAG}.json'}
    return V.O.enc(dict(status='PASS_D16_NODE3_GUARD_COMPLETE_FINITE',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(paths)},
        top_classes=top,sequential_increasing_family_uppers=covers,new_exact_global_formal_maximum='1487/4',
        attaining_forest=witness,curvature='364',dimension_bound=57,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);M.R.resources();U.memory_guard()
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();M.R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_exact_global_formal_maximum'],[(r['fixed_rank16_count'],r['maximum'])for r in out['top_classes']['classes']],flush=True)
