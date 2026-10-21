#!/usr/bin/env python3
"""Exact capacity-one cliques of incompatible arithmetic cap profiles."""
import argparse
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

import fresh_sqrt5_cross_coface_cap_capacity_20261017_review as V

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_cap_profile_cliques_20261017'
TAGS=dict(old='fresh_sqrt5_cross_orientation_cap_capacity_20261017',
    incidence='fresh_sqrt5_prism_simplex_base_incidence_20261017',
    model='fresh_sqrt5_cross_coface_cap_full_model_20261017')
K=V.K


def integral(z):
    basis=(z.a-z.b).denominator==1 and (2*z.b).denominator==1
    polynomial=(2*z.a).denominator==1 and (z.a*z.a-5*z.b*z.b).denominator==1
    assert basis==polynomial
    return polynomial


def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for obj in (d,d.get('prepared',{})):
        for p,h in obj.get('sources',{}).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d


def run():
    D={k:read(t) for k,t in TAGS.items()};old,I,M=(D[k] for k in ('old','incidence','model'))
    assert old['status']=='PASS_EXACT_CROSS_ORIENTATION_CAP_CAPACITIES'
    dual=next(c for c in M['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={i:F(a) for i,a in dual['face_multipliers']}
    groups=defaultdict(list)
    for c in old['columns']:
        if c['capacity']==1:groups[c['base_type']].append(c)
    inverse={};edges=[];examined=[]
    for bt,cs in sorted(groups.items()):
        for a,b in combinations(cs,2):
            assert a['core']==b['core'] and a['simplex_face']==b['simplex_face']
            assert not set(a['source_bases']) & set(b['source_bases'])
            core=I['cores'][a['core']]
            record=dict(base_type=bt,columns=[a['index'],b['index']],compatible_not_claimed=True)
            if not all(e[2][0] in (3,5) for e in core['edges']):
                record['reason']='outside_selected_physical_Qsqrt5_core_domain';examined.append(record);continue
            if a['core'] not in inverse:
                G,scale=V.physical_core(V.core_graph(core));assert scale is None
                J,det=V.inverse(G);assert det.sign()<0
                n=len(G)
                assert all(sum((G[i][k]*J[k][j] for k in range(n)),K())==int(i==j)
                           for i in range(n) for j in range(n))
                inverse[a['core']]=J
            J=inverse[a['core']];checks=[];reason=None
            for (i,x),(j,y) in product(enumerate(a['orientation_orbit']),enumerate(b['orientation_orbit'])):
                overlap=sorted(set(x[1]) & set(y[1]))
                if not overlap:reason='disjoint_base_subsets_not_tested';break
                if not all(e[0] in (2,3) for e in x[0]+y[0]):reason='outside_selected_contact_domain';break
                u=[i for i,e in enumerate(x[0]) if e[0]==3]
                v=[i for i,e in enumerate(y[0]) if e[0]==3]
                assert sum((J[k][l] for k in u for l in u),K())==K(*a['repeated_wall_product'])
                assert sum((J[k][l] for k in v for l in v),K())==K(*b['repeated_wall_product'])
                t=sum((J[k][l] for k in u for l in v),K());w=t*t
                if integral(w):reason='integral_two_cycle_not_excluded';break
                checks.append(dict(orientations=[i,j],common_base_vertices=overlap,
                    cross_product=t.out(),cross_two_cycle=w.out(),
                    trace=str(2*w.a),norm=str(w.a*w.a-5*w.b*w.b)))
            if reason is None:
                assert len(checks)==len(a['orientation_orbit'])*len(b['orientation_orbit'])
                record['reason']='all_orientation_pairs_excluded_at_same_and_distinct_cofaces'
                record['checks']=checks;edges.append(record)
            else:record['reason']=reason
            examined.append(record)
    edge_set={tuple(e['columns']) for e in edges}
    clique_profiles=[]
    for bt,cs in sorted(groups.items()):
        verts=sorted(c['index'] for c in cs);cliques=[]
        for size in range(2,len(verts)+1):
            for vs in combinations(verts,size):
                if all(p in edge_set for p in combinations(vs,2)):cliques.append(vs)
        clique_profiles.extend((bt,vs) for vs in cliques if not any(set(vs)<set(ws) for ws in cliques))
    columns=deepcopy(old['columns']);required=defaultdict(F)
    for bt,vs in clique_profiles:
        cs=[old['columns'][i] for i in vs];a=cs[0];co=Counter();pr=Counter();bases=[]
        for c in cs:
            co.update(dict(c['positive_face_coefficients']));pr.update(dict(c['provenance_coefficients']))
            bases.extend(c['source_bases'])
        assert len(bases)==len(set(bases))
        columns.append(dict(index=len(columns),kind='distinct_profile_capacity_one_clique',
            base_type=bt,simplex_face=a['simplex_face'],simplex_dimension=a['simplex_dimension'],
            core=a['core'],capacity=1,profile_columns=list(vs),source_bases=sorted(bases),
            positive_face_coefficients=sorted(map(list,co.items())),
            provenance_coefficients=sorted(map(list,pr.items()))))
    for c in columns:
        load=sum((m*alpha.get(fi,F()) for fi,m in c['positive_face_coefficients']),F())
        mass=alpha.get(c['simplex_face'],F())
        c.update(load=str(load),initial_simplex_mass=str(mass),slack=str(c['capacity']*mass-load))
        required[c['simplex_face']]=max(required[c['simplex_face']],load/c['capacity'])
    new=columns[len(old['columns']):]
    assert len(examined)==80 and len(edges)==70 and len(new)==36
    assert Counter(len(c['profile_columns']) for c in new)=={2:22,3:12,4:2}
    assert sum(F(c['slack'])<0 for c in new)==4
    assert len(columns)==851 and sum(F(c['slack'])<0 for c in columns)==46
    paths=[f'scripts/{TAG}.py',str(Path(V.__file__).relative_to(ROOT))]
    paths += [f'results/{t}.json' for t in TAGS.values()]
    return dict(status='PASS_EXACT_DISTINCT_PROFILE_CAP_CLIQUES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=dict(inherited_columns=815,columns=len(columns),examined_profile_pairs=len(examined),
            excluded_profile_pairs=len(edges),new_clique_columns=len(new),
            positive_clique_columns=sum(F(c['load'])>0 for c in new),
            new_violations=4,violations=46),
        pair_checks=examined,incompatibility_edges=edges,
        clique_size_histogram=sorted(Counter(len(c['profile_columns']) for c in new).items()),
        columns=columns,violating_columns=[c['index'] for c in columns if F(c['slack'])<0],
        required_simplex_masses=[[i,str(alpha.get(i,F())),str(a)] for i,a in sorted(required.items())],
        input_bound=dual['bound'],new_repaired_bound=None,optimizer_invocations=0,
        scope='Add36 capacity-one clique columns on an explicitly selected physical-core/contact domain. Each profile already has capacity one; all cross-orientation pairs at a common base type have overlapping bases at the same coface and a nonintegral wall two-cycle at distinct orthogonal A1 cofaces. All815 preceding columns are preserved. No complete compatibility classification, optimum or dimension bound is claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    r=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert r==json.loads(out.read_text())
    else:out.write_text(json.dumps(r,indent=2)+'\n')
    print(r['status'],r['counts'],flush=True)
