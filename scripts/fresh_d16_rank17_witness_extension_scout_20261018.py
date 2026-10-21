#!/usr/bin/env python3
"""Separate R17 witness extension scout; frozen D16 rule is unchanged."""
import argparse
from collections import Counter
from fractions import Fraction as F
import fcntl
from itertools import combinations_with_replacement
import json
from pathlib import Path
import signal
import fresh_d16_node3_guard_complete_finite_20261018 as B

ROOT=B.ROOT;C=B.C;V=B.V;M=B.M;U=B.U
TAG='fresh_d16_rank17_witness_extension_scout_20261018'


def value(fs):
    co=Counter(fs);r16=co['A',16]+co['D',16];r17=co['A',17]+co['D',17]
    if not V.has_guard(fs)or not r17:return B.value(fs)
    debit=-64*co['D',16]if r16<2 else 0
    credit=256*co['A',2]*co['A',13]if co['A',1]and not r16 else 0
    return B.value(fs)+debit+credit


def coupon_upper():
    cs=(0,0,0,0,1,0,0,0,1);key=next(key for counts,used,key in C.fixed_domain()if counts==cs and key[9]==0)
    assert all(C.actual_cross(t,('D',17))==0 for t in C.TYPES if t[1]<=6 or t==('A',13))
    pad=max(0,C.actual_base(('D',17))-C.upper_value([('A',17)]));assert pad==16
    rows=[];tabs=[];checks=preds=0
    for k5 in(0,1):
        fixed=C.fixed_components(cs,key)+[('A',5)]*k5;fc=C.upper_value(fixed)+pad
        ws={r:C.item(r,key,k5)+(256 if r==2 else 128 if r==4 else 0)for r in C.ITEMS}
        tab,a,b=U.table(ws);checks+=a;preds+=b;tabs.append(dict(k5=k5,table=tab))
        for s in range(18,25-5*k5):
            n=27-5*k5-s;opts=[]
            assert C.item(s,key,k5)==C.actual_base(('D',s))
            for a1 in range(1,n-1):
                for a2 in range(1,(n-a1)//2+1):
                    left=n-a1-2*a2;v=tab['values'][left]
                    if v is not None:opts.append((fc+C.actual_base(('D',s))+a1*ws[1]+a2*ws[2]+v,a1,a2,left))
            best=max(opts);rows.append(dict(k5=k5,reserved_D=s,residual=n,upper=F(best[0],64),explicit_counts=best[1:]))
    assert len(rows)==9 and max(r['upper']for r in rows)==F(22917,64)
    return dict(classes=2,tables=tabs,rows=rows,witness_D17_base_pad='1/4',
        Bellman_inequalities=checks,attaining_predecessors=preds,increasing_family_upper='22917/64')


def complete_refinement():
    old=B.selected_array(ROOT/f'results/{C.TAG}.json',b'"above_curvature_classes": [')
    high=[r for r in old if F(r['maximum'])>F(743,2)];assert len(old)==38 and len(high)==5
    extra=[r for r in high if r['counts'][8]]
    assert{tuple(r['counts'])for r in extra}=={tuple([0]*7+[k,1])for k in(0,1)}
    assert all(r['k5']==r['a12']==0 for r in high)
    previous=V.load(B.TAG);assert max(F(r['maximum'])for r in previous['top_classes']['classes'])==371
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS}
    residual=[t for t in C.TYPES if C.representative(t)in reps and t!=('A',1)]
    fibers=[t for t in C.TYPES if C.representative(t)==C.representative(('D',16))]
    assert fibers==[('A',16),('B',16),('D',16)]
    assert C.representative(('A',17))==('A',17)and C.representative(('D',17))==('D',17)
    assert ('A',17)not in residual and ('D',17)in residual
    assert all(C.actual_cross(a,b)==0 for a in residual+[('A',1)]for b in residual+[('A',1)])
    out=[];checks=preds=options=0
    for k in(0,1):
        rows=[]
        for chosen in combinations_with_replacement(fibers,k):
            fixed=[('A',17)]+list(chosen);rem=40-16*k;tables=[];candidates=[]
            for a1 in range(rem+1):
                anchor=fixed+[('A',1)]*a1;left=rem-a1;bc=value(anchor)
                ts=[t for t in residual if t[1]<=left and not(t[0]=='D'and t[1]>=18)];ws=[value(anchor+[t])-bc for t in ts]
                tab,a,b=C.additive(ts,ws,left);checks+=a;preds+=b;tables.append(dict(a1=a1,guard=False,table=tab))
                if tab['values'][left]is not None:candidates.append((bc+tab['values'][left],a1,False,None,left,len(tables)-1));options+=1
                if left<18:continue
                anchorg=anchor+[('D',18)];cap=left-18;bcg=value(anchorg);ts=[t for t in residual if t[1]<=cap]
                ws=[value(anchorg+[t])-bcg for t in ts];tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b;tables.append(dict(a1=a1,guard=True,table=tab))
                for s in range(18,left+1):
                    n=left-s;v=tab['values'][n]
                    if v is None:continue
                    fs=anchor+[('D',s)];base=value(fs)
                    for t,w in zip(ts,ws):
                        if t[1]<=n:assert value(fs+[t])-base==w
                    candidates.append((base+v,a1,True,s,n,len(tables)-1));options+=1
            best=max(candidates);z,a1,guard,s,left,ti=best;fs=fixed+[('A',1)]*a1+([('D',s)]if guard else[]);tab=tables[ti]['table']
            while left:
                t=tab['types'][tab['last'][left]];fs.append(t);left-=t[1]
            assert sum(t[1]for t in fs)==57 and value(fs)==z
            rows.append(dict(fixed=fixed,tables=tables,maximum=F(z,64),attaining_forest=fs,winning_option=best[1:]))
        out.append(dict(fixed_rank16_count=k,rows=rows,maximum=max(r['maximum']for r in rows)))
    assert max(r['maximum']for r in out)<=F(743,2)
    witness=[('A',17),('A',17),('D',23)];assert value(witness)==64*F(743,2)
    return dict(original_compact_class_maxima=old,newly_refined_old_classes=extra,
        fixed_rank17_type='A17_only_D17_is_residual',actual_rank16_fiber=fibers,
        original_residual_pair_crosses_zero=True,
        old_three_top_classes_upper='371',classes=out,Bellman_inequalities=checks,
        attaining_predecessors=preds,explicit_options=options,conditional_complete_new_maximum='743/2',attaining_forest=witness)


def run():
    upper=coupon_upper();finite=complete_refinement();baseline=V.load(B.TAG)
    covers=list(map(F,baseline['sequential_increasing_family_uppers']))+[F(upper['increasing_family_upper'])]
    assert max(covers)<F(743,2)
    sources={Path(__file__),Path(B.__file__),Path(M.__file__),Path(C.__file__),
        ROOT/f'results/{B.TAG}.json',ROOT/f'results/{C.TAG}.json',
        ROOT/'results/fresh_d16_node3_rank16_large_d_guard_transfer_20261018.json'}
    return V.O.enc(dict(status='PASS_D16_RANK17_WITNESS_EXTENSION_SCOUT',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(sources)},
        rule='Extend the existing D16 debit/coupon witness from A16 orD16 to A16 orD16 orA17 orD17; retain the separate whole D>=18 guard.',
        coupon_upper=upper,finite=finite,sequential_increasing_family_uppers=covers,
        inherited_selected_geometry_unchanged=True,inherited_no_consumer_audit_unchanged=True,
        geometric_extension_proof_in_written_report=True,independent_review_complete=False,
        status_of_new_rule='scout_only_not_promoted',new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);M.R.resources();U.memory_guard();p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();M.R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['coupon_upper']['increasing_family_upper'],out['finite']['conditional_complete_new_maximum'],
        [(r['fixed_rank16_count'],r['maximum'])for r in out['finite']['classes']],flush=True)
