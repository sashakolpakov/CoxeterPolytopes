#!/usr/bin/env python3
"""Symmetry-reduced three-wall interior sections; exact actual-ridge gate."""
from __future__ import annotations
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import numpy as np
import fresh_h20_global_interior_sections_20261013 as base

ROOT=base.ROOT
SYMMETRIES=ROOT/'data/enumerations/fresh_h20_global_interior_sections_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_triple_central_20261013.json'
ORBIT_OUTPUT=ROOT/'data/enumerations/fresh_h20_triple_orbits_20261013.json'

def triple_orbits(perms):
    def enc(a,b,c):return 44100*a+210*b+c
    def dec(v):a,r=divmod(v,44100);b,c=divmod(r,210);return a,b,c
    remaining={enc(*x) for x in combinations(range(210),3)};result=[]
    while remaining:
        first=min(remaining);seen={first};queue=[first]
        while queue:
            a,b,c=dec(queue.pop())
            for g in perms:
                nxt=enc(*sorted((g[a],g[b],g[c])))
                if nxt not in seen:seen.add(nxt);queue.append(nxt)
        assert seen<=remaining;remaining-=seen
        result.append({'representative':[x+1 for x in dec(first)],'size':len(seen)})
        print('triple orbit',len(result),dec(first),len(seen),'remaining',len(remaining),flush=True)
    assert sum(r['size'] for r in result)==1521520
    return result

def build():
    Q,W,G,h=base.load();sym=json.loads(SYMMETRIES.read_text());perms=[r['wall_permutation'] for r in sym['ambient_symmetry_generators']]
    if ORBIT_OUTPUT.exists():orbits=json.loads(ORBIT_OUTPUT.read_text())['orbits']
    else:
        orbits=triple_orbits(perms)
        ORBIT_OUTPUT.write_text(json.dumps({'schema':'borcherds-triple-orbits-v1','symmetry_source':str(SYMMETRIES.relative_to(ROOT)),'triple_count':1521520,'orbits':orbits},indent=2)+'\n')
    cp=W@Q@h;U=np.array([int(-8//int(cp[i]))*W[i] for i in range(210)],dtype=object)
    assert all(int(v@Q@h)==-8 for v in U)
    rows=[]
    for orb in orbits:
        ids=[x-1 for x in orb['representative']]
        for negative in ids:
            positive=[i for i in ids if i!=negative]
            e=base.primitive(U[positive[0]]+U[positive[1]]-2*U[negative])
            assert int(e@Q@e)>0 and int(e@Q@h)==0
            R,labels,discarded=base.normals(Q,W,h,e);counts,witness=base.actual_bad_ridge(Q,R,h)
            if witness is not None:witness['normal_ambient_wall_labels_one_based']=[[x+1 for x in labels[k]] for k in witness['normal_indices']]
            row={'triple_one_based':orb['representative'],'triple_orbit_size':orb['size'],
                 'negative_wall_one_based':negative+1,'positive_walls_one_based':[x+1 for x in positive],
                 'section_normal':list(map(int,e)),'section_normal_square':int(e@Q@e),
                 'positive_projected_normal_rays':len(R),'pair_gate_counts':counts,
                 'actual_forbidden_ridge':witness,'status':'NON_COXETER_ACTUAL_RIDGE' if witness else 'UNRESOLVED'}
            rows.append(row)
            print('section',orb['representative'],negative+1,row['status'],flush=True)
    out={'schema':'central-three-wall-sections-v1','source':str(base.SOURCE.relative_to(ROOT)),
         'orbit_source':str(ORBIT_OUTPUT.relative_to(ROOT)),'triple_count':1521520,'triple_orbit_count':len(orbits),
         'normal_definition':'u_i+u_j-2u_k, where u_i=8*r_i/(-(h,r_i)); all pass through h',
         'section_count':len(rows),'sections':rows}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    print('DONE',len(rows),Counter(r['status'] for r in rows),flush=True)

if __name__=='__main__':build()
