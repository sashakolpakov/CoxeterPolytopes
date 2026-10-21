#!/usr/bin/env python3
"""Exhaust all actual types in the three remaining representative rows.

The complete frozen Bellman certificate is an input, not rerun. This small
exact refinement removes the mixed-allowance rank4+2rank5+3rank15 row.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product,combinations_with_replacement
from pathlib import Path
import fresh_a14_cd16_ambient59_frontier_20261018 as C
ROOT=C.ROOT;TAG='fresh_ambient59_exact_typed_obstruction_frontier_20261018'
FINAL_MD='257219a7c1e324f44d4c7ed1c39d846a0994da807d2328f07b2387bad9ddb6cf'
enc=C.enc


def partitions(n,least=1):
    if n==0:yield ()
    for k in range(least,n+1):
        for tail in partitions(n-k,k):yield (k,)+tail


def actual_forests(ranks,types):
    return sorted({tuple(sorted(fs))for fs in product(*(types[r]for r in ranks))})


def run():
    old=json.loads((ROOT/f'results/{C.TAG}.json').read_text())
    assert old['status']=='PASS_A14_CD16_AMBIENT59_POSITIVE_IDEAL_FINITE787_OVER2'
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert sha256((ROOT/f'results/{C.TAG}.md').read_bytes()).hexdigest()==FINAL_MD
    ff=old['finite'];curv=16*F(ff['curvature']);assert curv==6264
    expected={(cs,k5)for cs,used,key in C.fixed_domain()for k5 in range((59-used)//5+1)}
    saved={(tuple(r['counts']),r['fixed_A5_count']):r for r in ff['expanded_alternatives']}
    assert len(saved)==len(ff['expanded_alternatives'])==len(expected)==2574 and set(saved)==expected
    above=[r for r in saved.values()if r['sixteenths_upper']>curv]
    assert len(above)==3 and sorted(r['sixteenths_upper']for r in above)==[6272,6276,6296]
    types={n:[t for t in C.C.P.Z.B.B.types()if t[1]==n]for n in range(1,60)}
    records=[];obstructions=[]
    for row in above:
        fixed=tuple(n for n,k in zip(range(9,18),row['counts'])for _ in range(k))+(5,)*row['fixed_A5_count']
        rem=59-sum(fixed);assert rem in (0,4)
        residuals=list(partitions(rem));forests=set()
        for rs in residuals:forests.update(actual_forests(fixed+rs,types))
        rows=[]
        for fs in sorted(forests):
            assert sum(n for k,n in fs)==59
            v=C.value(list(fs));assert v<=row['sixteenths_upper']
            entry=dict(forest=[[k,n]for k,n in fs],sixteenths_actual=v,cost=F(v,16),
                       above_curvature=v>curv,at_curvature=v==curv)
            rows.append(entry)
            if v>curv:obstructions.append(entry)
        maximum=max(r['sixteenths_actual']for r in rows)
        if rem:
            assert fixed==(15,15,15,5,5)and residuals==[(1,1,1,1),(1,1,2),(1,3),(2,2),(4,)]
            assert maximum==6248 # 390.5, below391.5.
        else:assert maximum==row['sixteenths_upper']
        records.append(dict(fixed_ranks=fixed,residual_rank=rem,residual_partitions=residuals,
            counts=row['counts'],fixed_rank5_count=row['fixed_A5_count'],
            prior_representative_upper=F(row['sixteenths_upper'],16),
            actual_maximum=F(maximum,16),actual_type_forests=len(rows),rows=rows))
    assert len(obstructions)==8 and max(r['sixteenths_actual']for r in obstructions)==6296
    assert sorted(Counter(tuple(sorted(n for k,n in r['forest']))for r in obstructions).values())==[3,5]
    assert sorted(r['sixteenths_actual']for r in obstructions)==[6268,6268,6272,6272,6276,6280,6288,6296]
    paths={Path(__file__),Path(C.__file__),ROOT/f'results/{C.TAG}.json',ROOT/f'results/{C.TAG}.md'}
    return enc(dict(status='PASS_AMBIENT59_EIGHT_ACTUAL_TYPED_OBSTRUCTIONS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        complete_saved_expanded_vectors=2574,representative_rows_above=3,typed_refinements=records,
        actual_type_forests_checked=sum(r['actual_type_forests']for r in records),
        removed_representative_pattern='rank4 + 2rank5 + 3rank15',removed_pattern_maximum='781/2',
        exact_finite_maximum='787/2',finite_curvature='783/2',actual_obstructions=obstructions,
        actual_obstruction_count=8,actual_rank_patterns=2,
        all_other_actual_forests_at_most_curvature=True,
        inherited_complete_Bellman_certificate=C.TAG,Bellman_inequalities_replayed=False,
        ideal_rule_unchanged=True,uniform_actual_cusp_surplus='91/2',
        completed_unrestricted_bound=59,new_dimension_bound=None,
        geometric_realizations_claimed=False,
        scope='Exact actual-type refinement of every expanded representative row above curvature; full rank4 residual partitions included. The frozen complete finite certificate bounds every other forest. No geometry or old Bellman computation replay.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['actual_type_forests_checked'],out['actual_obstruction_count'],flush=True)
