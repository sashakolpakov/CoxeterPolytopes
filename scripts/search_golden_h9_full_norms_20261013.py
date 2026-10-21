#!/usr/bin/env python3
"""Chronological exact mixed norm-two/four search on <-2phi>+E8+<2>."""
from __future__ import annotations
import argparse
from fractions import Fraction
import json
from pathlib import Path
import time

from coxeter_polytopes.exact.golden import Golden as G, PHI, ZERO, ONE
from coxeter_polytopes.search.golden_path_host import bounded_golden_integers
from coxeter_polytopes.search.og10_golden_h9_vinberg import (
    _enumerate_height_level, _positive_definite_enumeration_data,
    h9_cone_roots, h9_field_gram,
)
from coxeter_polytopes.search.bugaenko_vinberg import _pair, _mat_vec
from coxeter_polytopes.verify.chamber import diagnose_chamber_cone, verify_hyperideal_extreme_ray


def run(height,output,skip_geometry=False,a1_norm=2):
    start=time.monotonic()
    controls=[]
    norm_pool=(2,4) if a1_norm==2 else (1,2)
    for q in norm_pool:
        for x in bounded_golden_integers((height*q/(2*PHI),-G(q)/(2*PHI))):
            if x.sign()>0:
                controls.append((2*PHI*x*x/q,q,x))
    controls.sort()
    roots=[]
    levels=[]
    inverse,schur=_positive_definite_enumeration_data(a1_norm)
    gram=h9_field_gram()
    if a1_norm!=2:
        gram=tuple(tuple(G(a1_norm) if i==j==9 else gram[i][j] for j in range(10)) for i in range(10))
    for i,(h,q,x) in enumerate(controls):
        level=_enumerate_height_level(x,roots,inverse,schur,norm=q,a1_norm=a1_norm)
        assert all(_pair(r,gram,s).sign()<=0 for r in level.roots for s in roots)
        roots.extend(level.roots)
        row={"height":h.to_json(),"norm":q,"control":x.to_json(),"new_roots":len(level.roots),"total_positive_walls":len(roots),"nodes":sum(level.visited_by_depth)}
        levels.append(row)
        if level.roots:
            print(json.dumps({"level":i,**row}),flush=True)
    walls=(*h9_cone_roots(),*roots)
    # A tiny positive E8-weight/last-coordinate perturbation of p.
    spatial=_mat_vec(inverse,(ONE,)*8)
    epsilon=G(Fraction(1,2))
    while True:
        interior=(ONE,*(epsilon*v for v in spatial),epsilon)
        if _pair(interior,gram,interior).sign()<0 and all(_pair(interior,gram,w).sign()<0 for w in walls):
            break
        epsilon/=2
    payload={"host":f"<-2phi>+E8+<{a1_norm}>","a1_norm":a1_norm,"dimension":9,"norms":list(norm_pool),"height_bound":height.to_json(),"levels":levels,"positive_wall_count":len(roots),"wall_count":len(walls),"all_walls":[[v.to_json() for v in r] for r in walls],"strict_interior":[v.to_json() for v in interior],"scope":"Exact mixed-norm bounded prefix; no nonreflectivity conclusion from a nonclosed prefix."}
    output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    if skip_geometry:
        print(json.dumps({"phase":"geometry_skipped","walls":len(walls),"seconds":round(time.monotonic()-start,3)}),flush=True)
        return
    print(json.dumps({"phase":"geometry","walls":len(walls)}),flush=True)
    geom=diagnose_chamber_cone(walls,gram,interior,9)
    payload['geometry']={"finite_vertices":len(geom.finite_vertices),"cusps":len(geom.cusps),"hyperideal_rays":len(geom.hyperideal_rays),"opposite_sheet_rays":len(geom.opposite_sheet_rays)}
    simple=next((r for r in geom.hyperideal_rays if len(r)==9),None)
    if simple is not None:
        w=verify_hyperideal_extreme_ray(walls,gram,simple,9)
        payload['hyperideal_witness']={"active_facets":list(w.active_facets),"vector":[v.to_json() for v in w.ray],"norm":w.norm.to_json(),"wall_pairings":[v.to_json() for v in w.wall_pairings]}
    payload['seconds']=round(time.monotonic()-start,3)
    output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({"wall_count":len(walls),"geometry":payload['geometry'],"seconds":payload['seconds']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--height',type=int,default=30)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--skip-geometry',action='store_true')
    p.add_argument('--a1-norm',type=int,choices=(1,2),default=2)
    a=p.parse_args()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    run(G(a.height),a.output,a.skip_geometry,a.a1_norm)
