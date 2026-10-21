#!/usr/bin/env python3
"""Complete small connected affine search orthogonal to forbidden-angle pairs."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from fresh_h20_null_link_profiles_20261016_review import add,scale,mul,sign

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_persistent_cusp_affine_20261016'
INPUT='results/fresh_h9_persistent_cusp_screen_20261016.json'
ZERO=(F(0),)*4
ONE=(F(1),F(0),F(0),F(0))
VALUES={2:ZERO,3:scale(ONE,-1),4:(F(0),F(-1),F(0),F(0)),
        6:(F(0),F(0),F(-1),F(0)),-1:scale(ONE,-2)}


@lru_cache(None)
def inverse(a):
    if len(a)==1:
        assert a[0]
        return (1/a[0],)
    half=len(a)//2
    x,y=a[:half],a[half:]
    d=2 if half==1 else 3
    den=inverse(add(mul(x,x),scale(mul(y,y),-d)))
    return mul(x,den)+scale(mul(y,den),-1)


@lru_cache(None)
def classify(matrix):
    q=[list(row) for row in matrix]
    positive=zero=0
    for i in range(len(q)):
        s=sign(q[i][i])
        if s<0:
            return 'indefinite',positive,zero
        if s==0:
            if any(q[i][j]!=ZERO for j in range(i+1,len(q))):
                return 'indefinite',positive,zero
            zero+=1
            continue
        positive+=1
        inv=inverse(q[i][i])
        for j in range(i+1,len(q)):
            factor=mul(q[j][i],inv)
            for k in range(j,len(q)):
                q[j][k]=q[k][j]=add(q[j][k],scale(mul(factor,q[i][k]),-1))
    return ('positive_definite' if not zero else 'positive_semidefinite'),positive,zero


def run():
    data=json.loads((ROOT/INPUT).read_text())
    assert data['status']=='PASS_EXACT_49_WALL_PAIR_INVENTORY_FOR_CUSP_SCREEN'
    M=data['pair_matrix']
    results=[]
    global_cache={}
    # The whole old chamber inventory is also needed when a new wall
    # introduces a forbidden-angle pair orthogonal to an old circuit.
    full=dict(walls=[],label=None,common_orthogonal=list(range(49)))
    for bad in data['forbidden_angle_edges']+[full]:
        domain=bad['common_orthogonal']
        queue=[(i,) for i in domain]
        seen=set(queue)
        records=[]
        affine=[]
        while queue:
            support=queue.pop()
            if support not in global_cache:
                gram=tuple(tuple(scale(ONE,2) if i==j else VALUES[M[i][j]] for j in support) for i in support)
                global_cache[support]=classify(gram)
            kind,rank,nullity=global_cache[support]
            records.append(dict(support=list(support),kind=kind,positive_pivots=rank,zero_pivots=nullity))
            if kind=='positive_semidefinite':
                assert nullity==1
                affine.append(list(support))
                continue
            if kind!='positive_definite':
                continue
            # The orthogonal complement of the forbidden-angle positive
            # two-plane has signature (7,1). No PD set can exceed seven;
            # a positive-rank-seven affine circuit is also impossible,
            # but size-eight is retained as a conservative search endpoint.
            assert len(support)<=(7 if bad['walls'] else 9)
            for new in domain:
                if new in support:
                    continue
                labels=[M[new][i] for i in support]
                if any(m not in VALUES for m in labels) or all(m==2 for m in labels):
                    continue
                candidate=tuple(sorted(support+(new,)))
                if candidate not in seen:
                    seen.add(candidate);queue.append(candidate)
        records.sort(key=lambda r:(len(r['support']),r['support']))
        results.append(dict(forbidden_pair=bad['walls'],label=bad['label'],
                            domain=domain,visited=len(records),counts=dict(Counter(r['kind'] for r in records)),
                            affine_circuits=sorted(affine),records=records))
        print('affine domain',bad['walls'],len(records),'circuits',len(affine),flush=True)
    paths=[f'scripts/{TAG}.py',INPUT,'scripts/fresh_h20_null_link_profiles_20261016_review.py']
    entire=results.pop()
    return dict(status='PASS_COMPLETE_OLD49_FORBIDDEN_EDGE_AFFINE_SCREEN',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                forbidden_pairs=len(results),unique_examined_subsets=len(global_cache),
                affine_circuits_total=sum(len(r['affine_circuits']) for r in results),domains=results,
                entire_old_chamber=entire,
                completeness='Every irreducible affine Coxeter diagram has finite labels2,3,4,6 or the rank1 parallel pair; every proper principal subdiagram is positive definite. A connected affine graph admits a connected-growth ordering. The recursion explores every such growth until exact indefiniteness or the first affine circuit. The positive two-plane of the retained forbidden pair bounds connected search size by8; the additional whole-chamber search uses the conservative ambient bound10.',
                scope='Complete inventory of irreducible affine circuits in the fixed49-wall chamber, and their orthogonality to each retained forbidden-angle pair. Absence of a forbidden cusp is not a finite-volume completion theorem or a guarantee that another obstruction is absent.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['affine_circuits_total'],flush=True)
