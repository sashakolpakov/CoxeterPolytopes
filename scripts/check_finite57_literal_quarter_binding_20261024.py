#!/usr/bin/env python3
"""SOURCE-ONLY prepared quarter-family binding; root owns all execution.

No legacy mathematical worker is imported. The new literal driver supplies
pure coefficient functions and bounded JSON adapters. Frozen upper-table
validity and the prior coefficient/domain batches are explicitly assumed.
Every result is an upper inequality, not a claim of literal attainment by
a relaxed residual reconstruction. No optimizing recurrence is executed.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from pathlib import Path
import signal
import importlib.util
import os

ROOT=Path(__file__).resolve().parents[1]
DRIVER='scripts/check_finite57_literal_binding_batches_20261024.py'
DRIVER_SHA='b3780a3978fd3ac9bf5517b831437f2b304640209641a3c00b83872d79ce4d0f'
BATCHES=('index',)+tuple(f'{kind}-{i}' for kind in ('coarse','fine') for i in range(4))
D=None


def canonical(fs):
    return tuple(sorted(map(D.typ,fs)))


def key(c):
    return (tuple(c['counts']),c['k5'],c['a12'])


def actual_sets(cs,k5,a12):
    """All repeated representative fibers, including mixed assignments."""
    groups=Counter(D.fixed(cs,a12,k5))
    fibers={r:tuple(t for t in D.types() if D.representative(t)==r) for r in groups}
    return {tuple(sorted(sum((list(x) for x in selection),[])))
        for selection in product(*(tuple(combinations_with_replacement(fibers[r],n)) for r,n in groups.items()))}


def expected_domain():
    return {(cs,k5,a12) for cs,k5,a12 in D.domain()
        if cs[1] and 58-sum(t[1] for t in D.fixed(cs,a12,k5))>=7}


def count_options(rem):
    # Complete explicit A1/A4 counts; a2 then has the indicated interval.
    return sum((rem-a1-4*a4)//2 for a1 in range(1,rem-5)
        for a4 in range(1,(rem-a1-2)//4+1))


def index():
    expected=expected_domain()
    assert len(expected)==539
    new=D.document(D.FRONT,('new_context_family',))['new_context_family']
    coarse=new['coarse_rows']
    actual=new['actual_refinements']
    got={(tuple(r['counts']),r['k5'],r['key'][9]) for r in coarse}
    assert got==expected and len(coarse)==539
    total=sum(count_options(58-sum(t[1] for t in D.fixed(cs,a12,k5))) for cs,k5,a12 in expected)
    assert total==new['complete_explicit_options']==37985
    high={(tuple(r['counts']),r['k5'],r['key'][9]) for r in coarse if F(r['upper'])>377}
    assert high=={key(c) for c in actual} and len(actual)==new['coarse_classes_above_curvature']==14
    fibers=0
    options=0
    for c in actual:
        cs,k5,a12=key(c)
        sets={fs for fs in actual_sets(cs,k5,a12) if ('A',10) in fs}
        assert {canonical(r['fixed_actual']) for r in c['rows']}==sets
        assert len(sets)==len(c['rows'])
        rem=58-sum(t[1] for t in D.fixed(cs,a12,k5))
        fibers+=len(sets)
        options+=len(sets)*count_options(rem)
    # Recorded reachable options may omit ranks unreachable in the inherited
    # actual table; this raw-domain counter is a conservative upper count.
    assert fibers==new['actual_fixed_fibers']==77
    assert new['actual_explicit_options']==23376 and options>=23376
    return dict(coarse_classes=539,coarse_options=total,refined_classes=14,actual_fixed_fibers=fibers,
        complete_raw_fine_options=options,recorded_reachable_fine_options=23376,
        scope='complete rank/fiber domains; reachable fine count established by the four fine shards')


def coarse(batch):
    shard=int(batch[-1]);checked=options=0
    new=D.document(D.FRONT,('new_context_family',))['new_context_family']
    rows=new['coarse_rows'];assert len(rows)==539
    tables=D.document(D.COARSE,('finite','tables'))['finite']['tables']
    lookup={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in tables}
    assert len(lookup)==len(tables)==6596
    for i,r in enumerate(rows):
        if i%4!=shard: continue
        cs=r['counts'];k5=r['k5'];a12=r['key'][9]
        fs=D.fixed(cs,a12,k5);rem=58-sum(t[1] for t in fs)
        assert cs[1] and rem==r['residual_rank']>=7
        vs=lookup[tuple(r['key']),k5,True]['values']
        assert len(vs)==59 and vs[0]==0
        ws={rank:D.item(rank,fs) for rank in D.ITEMS}
        fixed_cost=D.bar_fixed(fs)
        bound=64*F(r['upper'])
        row_options=0
        for a1 in range(1,rem-5):
            for a4 in range(1,(rem-a1-2)//4+1):
                for a2 in range(1,(rem-a1-4*a4)//2+1):
                    left=rem-a1-4*a4-2*a2
                    z=fixed_cost+a1*ws[1]+a2*ws[2]+a4*ws[4]+vs[left]
                    if a1>=2: z+=36*a2*(a2-1)
                    if cs[0]: z+=64*a2*k5+40*a1*k5
                    z+=22*cs[-1]+64*a2*a4
                    assert z<=bound,(i,a1,a2,a4,z,bound)
                    row_options+=1
        assert row_options==count_options(rem)
        options+=row_options;checked+=1
    return dict(shard=shard,coarse_rows=checked,complete_explicit_options=options,
        scope='printed barred coefficients bound every coarse quarter option by its frozen row bound',
        premises=['complete domain in index batch','frozen barred no2 upper tables, bound by first-driver coarse batches'])


def residual_alphabet():
    return {t for t in D.types() if t not in (('A',1),('A',2))
        and D.representative(t)[1] in D.ITEMS and D.representative(t)!=('A',17)}


def fine(batch):
    shard=int(batch[-1]);fixed_rows=options=weight_checks=0
    new=D.document(D.FRONT,('new_context_family',))['new_context_family']
    refinements=new['actual_refinements'];assert len(refinements)==14
    alphabet=residual_alphabet()
    with D.View(D.CONDITIONAL) as view:
        ancestors=view.array(('finite','classes'))
        ancestor_index={
            (tuple(ancestors[i]['counts']),ancestors[i]['fixed_A5_count'],ancestors[i]['fixed_A12_count']):i
            for i in range(len(ancestors))}
        assert len(ancestor_index)==len(ancestors)
        for ci,c in enumerate(refinements):
            if ci%4!=shard: continue
            cs,k5,a12=key(c)
            parent=ancestors[ancestor_index[cs,k5,a12]]
            saved={canonical(r['fixed_actual']):r for r in parent['rows']}
            assert set(saved)==actual_sets(cs,k5,a12) and len(saved)==len(parent['rows'])
            selected={fs for fs in saved if ('A',10) in fs}
            assert selected=={canonical(r['fixed_actual']) for r in c['rows']}
            assert len(selected)==len(c['rows'])
            for r in c['rows']:
                fs=list(canonical(r['fixed_actual']));old=saved[tuple(fs)]
                rem=58-sum(t[1] for t in fs)
                assert rem==parent['residual_rank']
                residual=list(map(D.typ,old['residual_types']))
                assert set(residual)==alphabet and len(residual)==len(alphabet)
                for t,w in zip(residual,old['weights'],strict=True):
                    assert D.base(t)+sum(D.cross(t,b) for b in fs)<=w
                    assert all(D.cross(t,a)==0 for a in (('A',1),('A',2),('A',4)))
                    weight_checks+=1
                values=old['values']
                assert len(values)==rem+1 and values[0]==0
                bound=64*F(r['upper']);assert bound<=64*F(2861,8)<377*64
                # The literal context formula uses no other residual type.
                # Its A3 predicate is already true because A1 is present.
                sensitive={('A',n) for n in (1,2,4,9,10,11,13,14,15,16)}|{('D',n) for n in (5,14,15,16)}
                assert not (set(residual)-{('A',3),('A',4)})&sensitive
                test=fs+[('A',1),('A',2),('A',4)]
                assert D.W.context_sum(test+[('A',3)])==D.W.context_sum(test)
                for a1 in range(1,rem-5):
                    for a4 in range(1,(rem-a1-2)//4+1):
                        for a2 in range(1,(rem-a1-4*a4)//2+1):
                            left=rem-a1-4*a4-2*a2
                            if values[left] is None: continue
                            explicit=fs+[('A',1)]*a1+[('A',2)]*a2+[('A',4)]*a4
                            # Every actual A4 was removed into its complete
                            # explicit count. Other residual types leave Psi
                            # fixed by the syntactic support/guard proof above.
                            z=D.value(explicit)+values[left]
                            assert z<=bound,(ci,fs,a1,a2,a4,z,bound)
                            options+=1
                fixed_rows+=1
    return dict(shard=shard,actual_fixed_rows=fixed_rows,reachable_explicit_options=options,
        literal_weight_domination_checks=weight_checks,
        scope='every actual quarter option has a dominating frozen residual upper; no literal attainment claim',
        premises=['index proves all14/77 fibers','saved ancestor upper tables are valid for their full recorded weight alphabet',
            'actual A1/A2/A4 counts are complete; extra A4 in the upper table is a relaxation only'])


def main():
    global D
    assert __debug__,'assertions required'
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',choices=BATCHES,required=True)
    p.add_argument('--parent-lock',action='store_true',required=True)
    args=p.parse_args()
    h=sha256()
    with (ROOT/DRIVER).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''): h.update(block)
    assert h.hexdigest()==DRIVER_SHA,'common driver changed before import'
    spec=importlib.util.spec_from_file_location('finite57_literal_binding_common',ROOT/DRIVER)
    D=importlib.util.module_from_spec(spec);spec.loader.exec_module(D)
    owner=D.parent_lock();signal.alarm(50)
    assert D.digest(DRIVER)==DRIVER_SHA
    for path,pin in D.PINS.items(): assert D.digest(path)==pin,('source/input pin',path)
    D.W=D.module(D.LIB,'finite57_literal_weights')
    D.PROJ=D.module(D.PROJECTOR,'finite57_bounded_projection')
    D.ARR=D.module(D.ARRAY,'finite57_bounded_arrays');D.ARR.PROJ=D.PROJ
    result=index() if args.batch=='index' else coarse(args.batch) if args.batch.startswith('coarse') else fine(args.batch)
    for path,pin in D.PINS.items(): assert D.digest(path)==pin,('changed source/input',path)
    assert D.digest(DRIVER)==DRIVER_SHA and os.getppid()==owner
    print(json.dumps(dict(status='PASS_DECLARED_LITERAL_QUARTER_SCOPE',batch=args.batch,result=result,
        source_sha256=D.digest(Path(__file__).relative_to(ROOT).as_posix()),driver_sha256=DRIVER_SHA,
        library_sha256=D.PINS[D.LIB],inputs_read=sorted(D.READ),native_worker_imports=False,
        full_manuscript_verified=False),sort_keys=True))


if __name__=='__main__': main()
