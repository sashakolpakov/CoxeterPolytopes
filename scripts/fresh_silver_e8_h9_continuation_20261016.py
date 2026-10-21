#!/usr/bin/env python3
"""Complete three-norm Vinberg continuation to height100, with prefix binding."""
import argparse,json
from hashlib import sha256
from pathlib import Path
import fresh_silver_e8_h9_20261015 as A
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_silver_e8_h9_continuation_20261016'
SOURCE='results/fresh_silver_e8_h9_20261015_h30.json'
def run():
    old=json.loads((ROOT/SOURCE).read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    out=A.run(A.S(100))
    assert out['levels'][:len(old['levels'])]==old['levels']
    assert out['all_walls'][:old['wall_count']]==old['all_walls']
    out['sources'].update({p:sha256((ROOT/p).read_bytes()).hexdigest()for p in[f'scripts/{TAG}.py',SOURCE]})
    out['matched_prefix']=dict(height=30,time_controls=old['time_controls'],walls=old['wall_count'])
    E=tuple(map(A.S.from_json,old['geometry']['hyperideal_witness']['vector']))
    N=[tuple(map(A.S.from_json,r))for r in out['all_walls']]
    values=[A.H.pair(E,A.GRAM,r)for r in N]
    out['old_escape_test']=dict(point=A.H.enc(E),norm=A.H.pair(E,A.GRAM,E).to_json(),
                              pairings=A.H.enc(values),cutting_walls=[i for i,v in enumerate(values)if v.sign()>0])
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved completeheight100 replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['wall_count'],d['geometry']['status'],d['old_escape_test']['cutting_walls'])
