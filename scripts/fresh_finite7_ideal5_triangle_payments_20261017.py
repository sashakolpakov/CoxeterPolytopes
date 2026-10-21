#!/usr/bin/env python3
"""Finite distance seven already suffices for the inherited triangle theorem."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_ideal5_finite8_triangle_payments_20261017 as T

ROOT=T.ROOT
TAG='fresh_finite7_ideal5_triangle_payments_20261017'


def run():
    old=json.loads((ROOT/f'results/{T.TAG}.json').read_text())
    census=json.loads((ROOT/f'results/{T.T.CENSUS}.json').read_text())
    for d in (old,census):
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    rows=[];exceptions=[]
    for diagram in census['diagrams']:
        if diagram['order']<7:continue
        adj,labels=T.T.graph(diagram)
        for marks in combinations(range(diagram['order']),3):
            corners=[]
            for a,b in combinations(marks,2):
                c=T.corner(adj,labels,marks,a,b)
                c['elementary_payment']=bool(c['elementary_payment'] and c['length']<=7)
                corners.append(c)
            paid=sum(c['elementary_payment'] for c in corners)
            rows.append(dict(diagram=diagram['index'],order=diagram['order'],marks=list(marks),
                             corners=corners,elementary_payment_count=paid))
            if paid<2:exceptions.append([diagram['index'],list(marks)])
    assert len(rows)==1025
    assert exceptions==old['elementary_exceptions']
    special=T.exceptional(next(d for d in census['diagrams']if d['index']==25))
    assert special==old['exceptional_resolution']
    assert [len(r['alternate_path'])-1 for r in special['rows']]==[6,6,6,7]
    for row in rows:
        if row['elementary_payment_count']<2:
            e=next(e for e in special['rows']if e['marks']==row['marks'])
            assert row['elementary_payment_count']==1
            assert not next(c for c in row['corners']if c['pair']==e['alternate_pair'])['elementary_payment']
    paths=([Path(__file__),Path(T.__file__),Path(T.T.__file__)]+
        [ROOT/f'{folder}/{s}.{ext}' for s in (T.TAG,T.T.CENSUS,'fresh_cusp_factor_payment_framework_20261017')
         for folder,ext in (('scripts','py'),('results','json'),('results','md'))])
    paths=list(dict.fromkeys(paths))
    return dict(status='PASS_FINITE7_IDEAL5_TRIANGLE_PAYMENTS',
        finite_distance_at_most=7,ideal_distance_at_most=5,
        ideal_feature_pair_distance_at_most=7,ideal_feature_distance_at_most=5,
        whole_ideal_factor_order_at_most=9,triangle_minimum=2,
        counts_by_order=dict(Counter(str(r['order'])for r in rows)),rows=rows,
        elementary_exceptions=exceptions,exceptional_resolution=special,
        new_elementary_exceptions=0,
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        scope='Complete1025 inherited high-core marked triples re-enumerated; inherited all-label low-core Steiner and E8 arbitrary-support proofs now use paths at most7. Shared graph/kernel routines; no independent arithmetic implementation, quadrilateral theorem, finite-forest maximum or ambient dimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['rows']),out['new_elementary_exceptions'],flush=True)
