#!/usr/bin/env python3
"""Enumerate only the new distance levels needed to test two selected H9 corners."""
import argparse,json,time
from hashlib import sha256
from pathlib import Path
import search_fresh_cubic7_h9_20261013 as V

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_cubic_corner_20261014_h141'
def dec(x):return V.C(x['a'],x['b'],x['c'])
def run(verify=False):
    path='results/fresh_cubic7_h9_seed1_h100_20261013.json';old=json.loads((ROOT/path).read_text())
    seed=V.SUPPORTED_E8_CUBE_SEEDS[1];gram=V.gram_for(seed.timelike_coefficient)
    walls=[tuple(dec(x)for x in r)for r in old['all_walls']]
    levels=V._control_levels(seed,141);previous=[l for l in levels if l.height<=100];todo=[l for l in levels if l.height>100]
    assert len(previous)==len(old['levels'])
    for a,b in zip(previous,old['levels']):assert a.height==dec(b['height'])and a.control==dec(b['control'])and a.norm==dec(b['norm'])
    start=time.monotonic();records=[]
    print('New control levels',len(todo),'known walls',len(walls),flush=True)
    for i,l in enumerate(todo):
        vectors,nodes=V.candidates(l.control,l.target,l.norm)
        accepted=[v for v in vectors if all(V.pair(v,gram,w).sign()<=0 for w in walls)]
        assert all(V.pair(v,gram,w).sign()<=0 for j,v in enumerate(accepted)for w in accepted[:j])
        for v in accepted:
            assert V.pair(v,gram,v)==l.norm
            assert all((2*sum((gram[j][k]*v[k]for k in range(10)),V.Z)/l.norm).is_integral for j in range(10))
        walls.extend(accepted);records.append(dict(control=l.control.to_json(),norm=l.norm.to_json(),height=l.height.to_json(),
            candidates=len(vectors),nodes=nodes,new_walls=[V.encode_vector(v)for v in accepted],wall_count=len(walls)))
        if not verify:
            (ROOT/f'results/{TAG}_partial.json').write_text(json.dumps(dict(status='RUNNING',completed_new_levels=len(records),total_new_levels=len(todo),levels=records),indent=2)+'\n')
        print('level',i+1,'/',len(todo),'walls',len(walls),'nodes',nodes,'seconds',round(time.monotonic()-start,1),flush=True)
    paths=[path,'scripts/search_fresh_cubic7_h9_20261013.py',f'scripts/{TAG}.py',
           'src/coxeter_polytopes/exact/cubic7.py','src/coxeter_polytopes/search/cubic7_e8_cubes.py',
           'src/coxeter_polytopes/search/cubic7_root_blocks.py']
    out=dict(status='COMPLETE_NEW_DISTANCE_LEVELS_100_TO_141',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        seed=1,height_bound=141,previous_height_bound=100,complete_new_levels=len(records),levels=records,
        wall_count=len(walls),normals=[V.encode_vector(v)for v in walls],timelike_coefficient=seed.timelike_coefficient.to_json(),
        scope='Targeted height141 continuation chosen from exact finite-corner heights140.3776, not a claim of global reflectivity. Previous levels are inherited and matched; all new levels are enumerated.')
    if verify:
        assert out==json.loads((ROOT/f'results/{TAG}.json').read_text())
    else:
        (ROOT/f'results/{TAG}.json').write_text(json.dumps(out,indent=2)+'\n');(ROOT/f'results/{TAG}_partial.json').unlink()
    print(out['status'],len(walls),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');run(ap.parse_args().verify)
