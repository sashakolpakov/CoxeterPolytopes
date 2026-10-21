#!/usr/bin/env python3
"""Exact cone census for a fresh golden determinant-neighbor prefix."""
import argparse
import json
from pathlib import Path
from coxeter_polytopes.exact.golden import Golden as G
from coxeter_polytopes.verify.chamber import diagnose_chamber_cone, verify_hyperideal_extreme_ray

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('artifact',type=Path)
p.add_argument('--escape-only',action='store_true')
a=p.parse_args()
d=json.loads(a.artifact.read_text())
matrix=lambda rows:tuple(tuple(G.from_json(x) for x in r) for r in rows)
roots=matrix(d['all_walls']);gram=matrix(d['gram'])
inside=tuple(G.from_json(x) for x in d['strict_interior'])
if a.escape_only:
    import numpy as np
    from scipy.optimize import linprog
    phi=(1+5**.5)/2
    real=lambda x:float(x.a)+float(x.b)*phi
    gf=np.array([[real(x) for x in r] for r in gram])
    f=np.array([[real(sum((gram[i][j]*r[j] for j in range(10)),G(0))) for i in range(10)] for r in roots])
    rng=np.random.default_rng(20261013)
    for _ in range(200):
        lp=linprog(rng.normal(size=9),A_ub=f[:,1:],b_ub=-f[:,0],bounds=[(None,None)]*9,method='highs')
        if not lp.success:continue
        v=np.r_[1.,lp.x]
        if v@gf@v<=1e-10:continue
        active=np.flatnonzero(abs(f@v)<1e-7)
        if len(active)!=9:continue
        simple=tuple(int(i) for i in active)
        try:
            verify_hyperideal_extreme_ray(roots,gram,simple,9)
        except Exception:continue
        break
    else:raise SystemExit('No exact escape discovered; finite volume remains undecided.')
    d['geometry']={"status":"EXACT_HYPERIDEAL_ESCAPE_NO_FULL_CENSUS"}
else:
    g=diagnose_chamber_cone(roots,gram,inside,9)
    d['geometry']={"finite_vertices":len(g.finite_vertices),"cusps":len(g.cusps),"hyperideal_rays":len(g.hyperideal_rays),"opposite_sheet_rays":len(g.opposite_sheet_rays)}
    simple=next((r for r in g.hyperideal_rays if len(r)==9),None)
if simple is not None:
    w=verify_hyperideal_extreme_ray(roots,gram,simple,9)
    d['hyperideal_witness']={"active_facets":list(w.active_facets),"vector":[x.to_json() for x in w.ray],"norm":w.norm.to_json(),"wall_pairings":[x.to_json() for x in w.wall_pairings]}
a.artifact.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
print(json.dumps({'source':str(a.artifact),'wall_count':len(roots),'geometry':d['geometry']}))
