#!/usr/bin/env python3
"""SOURCE-ONLY: bind ALL actual-fiber prefix/flag options to literal weights.

Root owns execution through the existing bounded runner and shared lock.
No legacy worker is imported, no certificate is generated, and no Bellman
recurrence is run. Frozen residual upper validity and complete fiber-domain
proofs are explicit inputs. A maximum attainer alone is not used as evidence
that other explicit options were charged correctly.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
import signal

ROOT=Path(__file__).resolve().parents[1]
DRIVER='scripts/check_finite57_literal_binding_batches_20261024.py'
DRIVER_SHA='b3780a3978fd3ac9bf5517b831437f2b304640209641a3c00b83872d79ce4d0f'
BATCHES=tuple(f'options-{i}' for i in range(4))
D=None


def options(batch):
    shard=int(batch[-1]);classes=rows=checked=unreachable=0
    alphabet={t for t in D.types() if t not in (('A',1),('A',2))
        and D.representative(t)[1] in D.ITEMS and D.representative(t)!=('A',17)}
    with D.View(D.FRONT) as view:
        groups=view.array(('finite','classes'));assert len(groups)==39
        for ci in range(len(groups)):
            if ci%4!=shard: continue
            group=groups[ci];rem=group['residual_rank']
            for row in group['rows']:
                fixed=list(map(D.typ,row['fixed_actual']));co=Counter(fixed)
                assert rem==58-sum(t[1] for t in fixed)
                tabs={tuple(t['key']):t for t in row['tables']}
                assert len(tabs)==len(row['tables'])
                # Bind each distinct flag table once, outside the count loops.
                for key,tab in tabs.items():
                    has1,no4,no3=key
                    assert all(type(v) is bool for v in key)
                    ts=list(map(D.typ,tab['types']))
                    allowed={t for t in alphabet if not(no4 and t==('A',4)) and not(no3 and t==('A',3))}
                    assert set(ts)==allowed and len(ts)==len(allowed)
                    expected=[D.base(t)+sum(D.cross(t,u) for u in fixed)
                        +(32*co['A',12] if has1 and t==('A',4) else 0) for t in ts]
                    assert expected==tab['weights']
                    assert all(D.cross(t,s)==0 for t in ts for s in (('A',1),('A',2)))
                    assert len(tab['values'])==rem+1 and tab['values'][0]==0
                used=set();bound=64*F(row['maximum'])
                # Prefix costs are derived directly from the entire literal
                # forest, including fixed-fixed interactions, w1/w2 and Psi.
                for a1 in range(rem+1):
                    for a2 in range((rem-a1)//2+1):
                        left=rem-a1-2*a2
                        no4=bool(a1 and a2 and co['A',10])
                        modes=[(False,0)] if a1 or a2 or not co['D',14] else [(True,0),(False,1)]
                        for no3,res3 in modes:
                            if left<3*res3: continue
                            key=(bool(a1),no4,no3);assert key in tabs
                            tab=tabs[key];used.add(key)
                            # This is the same literal item as the first-driver
                            # fiber batch; the equality is independent of native
                            # actual_base/actual_cross/context imports.
                            rest=left-3*res3;values=tab['values']
                            if values[rest] is None:
                                unreachable+=1;continue
                            prefix=fixed+[('A',1)]*a1+[('A',2)]*a2+[('A',3)]*res3
                            # Current Psi on every remaining type is fixed,
                            # except for the documented linear A4/A12 role.
                            # A3 cannot first activate D14 here: either A1/A2
                            # already activates it, it was reserved, or forbidden.
                            assert not (co['D',14] and not a1 and not a2 and not res3 and not no3)
                            assert not (a1 and a2 and co['A',10] and not no4)
                            total=D.value(prefix)+values[rest]
                            assert total<=bound,(ci,fixed,a1,a2,key,res3,total,bound)
                            checked+=1
                assert used==set(tabs),'uncovered or unused stored flag table'
                rows+=1
            classes+=1
    return dict(shard=shard,classes=classes,actual_fixed_fibers=rows,reachable_explicit_options=checked,
        unreachable_options=unreachable,scope='ALL literal fixed/A1/A2/reserved-A3 prefix costs and residual flags bound by each frozen maximum',
        complete_four_shard_totals=dict(actual_fixed_fibers=748,reachable_explicit_options=96528),
        premises=['frozen residual upper validity','complete39/748 domain in first-driver index/fiber batches',
            'symbolic Psi support/removal argument and current residual-zero condition'],
        attainer_only_argument=False)


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
    for path,pin in D.PINS.items(): assert D.digest(path)==pin,('source/input pin',path)
    D.W=D.module(D.LIB,'finite57_literal_weights')
    D.PROJ=D.module(D.PROJECTOR,'finite57_bounded_projection')
    D.ARR=D.module(D.ARRAY,'finite57_bounded_arrays');D.ARR.PROJ=D.PROJ
    result=options(args.batch)
    for path,pin in D.PINS.items(): assert D.digest(path)==pin,('changed source/input',path)
    assert D.digest(DRIVER)==DRIVER_SHA and os.getppid()==owner
    print(json.dumps(dict(status='PASS_DECLARED_LITERAL_ACTUAL_OPTION_SCOPE',batch=args.batch,result=result,
        source_sha256=D.digest(Path(__file__).relative_to(ROOT).as_posix()),driver_sha256=DRIVER_SHA,
        library_sha256=D.PINS[D.LIB],inputs_read=sorted(D.READ),native_worker_imports=False,
        full_manuscript_verified=False),sort_keys=True))


if __name__=='__main__': main()
