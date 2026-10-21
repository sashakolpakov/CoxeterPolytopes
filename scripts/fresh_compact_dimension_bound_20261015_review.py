#!/usr/bin/env python3
"""Independent h-vector and bounded-multiplicity count, not a weight reproof."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb,prod
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact_dimension_bound_20261015_review'

def run():
    source='scripts/fresh_compact_dimension_bound_20261015.py'
    data_path='results/fresh_compact_dimension_bound_20261015.json'
    data=json.loads((ROOT/data_path).read_text())
    for p,h in data['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    assert data['exclusions_proved']==0
    assert F(comb(29,2),1)/F(29,7)==98
    assert 5*98-406==84
    hcoeff=[]
    for i in range(15):
        q=comb(29-i,2)+comb(i,2)-196
        assert q==(14-i)*(15-i)
        hcoeff.append(q)
    assert 5*hcoeff[0]==1050
    assert F(1050,87-84)==350

    # Independent multiplicity-vector enumeration instead of recursive partitions.
    ranks=(4,7,8,9,10);type_counts=(5,4,4,3,3)
    rows=[]
    for mult in product(*(range(29//r+1)for r in ranks)):
        if sum(m*r for m,r in zip(mult,ranks))!=29 or mult[-1]>=3:continue
        ranklist=tuple(r for m,r in zip(mult,ranks)for _ in range(m))
        count=prod(comb(k+m-1,m)for k,m in zip(type_counts,mult))
        no_h4=prod(comb(k+m-1,m)for k,m in zip((4,4,4,3,3),mult))
        rows.append(dict(ranks=list(ranklist),typed_count=count,without_H4_count=no_h4))
    rows.sort(key=lambda x:x['ranks'])
    assert len(rows)==9
    assert sum(x['typed_count']for x in rows)==2126
    assert sum(x['without_H4_count']for x in rows)==1354
    saved={tuple(x['ranks']):x for x in data['extremal_rank_partitions']}
    for row in rows:
        old=saved[tuple(row['ranks'])]
        assert all(old[k]==v for k,v in row.items())
        assert len(old['typed_forests'])==row['typed_count']
    assert [x['coefficient']for x in data['h_vector_slack_coefficients']]==hcoeff
    paths=[f'scripts/{TAG}.py',source,data_path]
    return dict(status='PASS_INDEPENDENT_COMPACT29_REDUCTION_ARITHMETIC_AND_SCOPE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        arithmetic='Standalone Fraction/binomial multiplicity count; no producer imports.',
        h_coefficients=hcoeff,rank_partitions=rows,typed_forests=2126,H4_free=1354,
        extremal_vertex_count_at_least=350,
        inherited='Vinberg connected-component weights and global face inequality; finite spherical classification. No duplicate graph/end-segment census.',
        proof_audit='Exact h-vector identity, deficit inference, all-label scope and arithmetic H4 field restriction pass.',
        scope='Necessary corner types, not an exclusion of compact dimension29. Zeroing weights on >=5-gons preserves face lower bounds and proves an extremal listed vertex incident to a triangle/quadrilateral exists; no star rejection is performed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
