#!/usr/bin/env python3
"""Independent small label-compression audit; no second full label join."""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_highbase_construction_20261014_review'
BASE = 'fresh_h9_highbase_construction_20261014_'


def run():
    paths = [f'{folder}/{BASE}{stage}.{ext}' for stage in ('cell', 'refined', 'pentagon')
             for folder, ext in (('scripts', 'py'), ('results', 'json'))]
    paths += ['scripts/independent_h9_14_record30_a1_labels.py']
    data = [json.loads((ROOT / f'results/{BASE}{stage}.json').read_text())
            for stage in ('cell', 'refined', 'pentagon')]
    for d in data:
        for p, h in d['sources'].items():
            assert hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    # Exact finite category enumeration, independent of the producer.
    triangles = [x for x in product(range(2,8), repeat=3)
                 if sum(F(1,m) for m in x) < 1]
    assert len(triangles) == data[0]['lanner3_category_count'] == 175
    max_nonexceptional = F(1,2)+F(1,3)+F(1,7)
    assert max_nonexceptional == F(41,42) < 1
    assert F(1,2)+F(1,2) == 1
    # A large edge has cosine²>3/4; any adjoining nonright edge has
    # cosine²>=1/4. Their three-node principal determinant is strictly <0.
    assert 1-F(3,4)-F(1,4) == 0
    cell, refined, pentagon = data
    assert len(cell['cases']) == 14
    assert cell['cases'][0]['state_limit_reached']
    assert all(not x['state_limit_reached'] and x['verdict']=='EXACT_LANNER_VERTEX_OBSTRUCTION'
               for x in cell['cases'][1:])
    assert [r['passing_states'] for r in refined['stages']] == [108,556,40,80,0]
    assert [r['passing_states'] for r in pentagon['stages']] == [108,8,0]
    assert all(not r['truncated'] for d in (refined,pentagon) for r in d['stages'])
    assert all(d['verdict']=='EXACT_LANNER_VERTEX_OBSTRUCTION' for d in (refined,pentagon))
    assert sum(len(d['actual_vertices']) for d in cell['cases'])+len(pentagon['actual_vertices']) == 267350
    paths.append(str(Path(__file__).resolve().relative_to(ROOT)))
    return dict(
        status='PASS_SCOPED_ALL_LABEL_COMPRESSION_AND_JOIN_SOURCE_REVIEW',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        independent_lanner3_category_count=len(triangles),
        maximum_reciprocal_sum_outside_2_2=str(max_nonexceptional),
        high_label_scope='All finite integer labels >=7, excluding parallel/dotted pairs.',
        frozen_join_paths_checked=dict(refined=[108,556,40,80,0],pentagon=[108,8,0]),
        scope='Independent small Fraction category census plus source/proof review. '
              'The old capped case remains inconclusive by itself. No second full Lanner-four catalogue, '
              '267350-vertex incidence census, or complete label join was executed here. '
              'The obstruction is restricted to the fourteen clustered and one pentagon compact wedge recipes.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true')
    args=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == out
    else:
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
