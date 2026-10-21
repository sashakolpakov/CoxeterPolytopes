#!/usr/bin/env python3
"""Find a cone escape numerically, then certify it entirely in Q(phi)."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from coxeter_polytopes.exact.golden import Golden as G
from coxeter_polytopes.search.og10_golden_h9_vinberg import h9_field_gram
from coxeter_polytopes.search.bugaenko_vinberg import _pair, _mat_vec
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('artifact',type=Path)
a=p.parse_args()
payload=json.loads(a.artifact.read_text())
walls=tuple(tuple(G.from_json(v) for v in r) for r in payload['all_walls'])
gram=h9_field_gram()
if payload.get('a1_norm',2)!=2:
    gram=tuple(tuple(G(payload['a1_norm']) if i==j==9 else gram[i][j] for j in range(10)) for i in range(10))
phi=(1+5**0.5)/2
real=lambda v:float(v.a)+float(v.b)*phi
f=np.array([[real(sum((gram[i][j]*r[j] for j in range(10)),G(0))) for i in range(10)] for r in walls])
rng=np.random.default_rng(20261013)
for k in range(200):
    result=linprog(rng.normal(size=9),A_ub=f[:,1:],b_ub=-f[:,0],bounds=[(None,None)]*9,method='highs')
    if not result.success:
        continue
    active=np.flatnonzero(abs(f[:,1:]@result.x+f[:,0])<1e-7)
    if len(active)!=9:
        continue
    try:
        w=verify_hyperideal_extreme_ray(walls,gram,tuple(int(i) for i in active),9)
    except Exception:
        continue
    # Repeat the key positivity and every wall inequality independently.
    assert _pair(w.ray,gram,w.ray)==w.norm and w.norm.sign()>0
    assert all(_pair(w.ray,gram,r).sign()<=0 for r in walls)
    proof={'source':str(a.artifact),'sample':k,'wall_count':len(walls),'active_facets':list(w.active_facets),'vector':[v.to_json() for v in w.ray],'norm':w.norm.to_json(),'wall_pairings':[v.to_json() for v in w.wall_pairings],'verdict':'EXACT_HYPERIDEAL_ESCAPE_VERIFIED','discovery':'Numerical LP supplies candidate support only; exact kernel, norm and every inequality decide the verdict.'}
    out=a.artifact.with_name(a.artifact.stem+'_escape.json')
    out.write_text(json.dumps(proof,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'output':str(out),'wall_count':len(walls),'norm':str(w.norm),'sample':k}))
    break
else:
    raise SystemExit('No exact witness found; this is undecided, not a finite-volume verdict.')
