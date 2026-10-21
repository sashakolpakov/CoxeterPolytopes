#!/usr/bin/env python3
"""Five wedges on an induced pentagon in the 600-cell, not a tetrahedral cluster."""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_h9_highbase_construction_20261014_cell as C

ROOT=C.ROOT;TAG='fresh_h9_highbase_construction_20261014_pentagon'
SOURCE=f'results/{C.TAG}.json';OUT=ROOT/f'results/{TAG}.json'


def run():
    source=json.loads((ROOT/SOURCE).read_text())
    for p,h in source['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    b=source['base'];missing=set(map(tuple,b['minimal_nonfaces']))
    adjacent=lambda i,j:i!=j and tuple(sorted((i,j))) not in missing
    edge=next((i,j) for i,j in combinations(range(120),2) if adjacent(i,j))
    chosen=tuple(k for k in range(120) if all(adjacent(k,i) for i in edge))
    assert len(chosen)==5 and all(sum(adjacent(i,j) for j in chosen)==2 for i in chosen)
    assert sum(tuple(sorted(p)) in missing for p in combinations(chosen,2))==5
    assert not any(all(adjacent(i,j) for i,j in combinations(f,2)) for f in combinations(chosen,4))
    weights,blocks,M,V=C.incidence(b,chosen)
    print('edge-link pentagon',edge,chosen,'vertices',len(V),flush=True)
    trial=C.join(M,V,20000)
    paths=[SOURCE,f'scripts/{TAG}.py',f'scripts/{C.TAG}.py','scripts/independent_h9_14_record30_a1_labels.py']
    return json.loads(json.dumps(dict(status='EXACT_600CELL_PENTAGON_WEDGE_CHECKPOINT',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        edge=edge,doubled_vertices=chosen,weights=weights,blocks=blocks,
        minimal_nonfaces=M,actual_vertices=V,**trial,
        scope='A new five-wedge recipe supported on the pentagonal link of one old 600-cell edge. '
              'No four doubled base vertices form a facet, so it lies outside the tetrahedron-plus-one family. '
              'Necessary complete Lanner-category and finite-vertex joins only; a surviving label state is not a metric realization.')))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['verdict'])
