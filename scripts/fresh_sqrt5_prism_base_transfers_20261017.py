#!/usr/bin/env python3
"""Export signed zero-budget prism/base transfers; exactly repair a saved dual."""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from math import comb
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
import fresh_sqrt5_rank13_simplex_quad_combined_20261017 as P
import fresh_sqrt5_prism_simplex_base_incidence_20261017 as I

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_prism_base_transfers_20261017'


def run():
    incidence=P.read(I.TAG);d=P.read(I.LP)
    old=next(c for c in d['certificates'] if c['kind']=='EXACT_RESTRICTED_MODEL_DUAL')
    B.dual_bound=D.sparse_dual_bound;B.verify_certificate(d,old)
    alpha={i:F(a) for i,a in old['face_multipliers']}
    byface=defaultdict(list);context_map={}
    for r in incidence['contexts']:
        ctx,fi,ci,k,entries=r
        assert ctx not in context_map;context_map[ctx]=r
        byface[fi].append((ctx,ci,k,dict(entries)))
    raw=defaultdict(list);unknown=Counter();varying=0;all_provenances=0
    for fi,records in sorted(byface.items()):
        face=d['faces'][fi];assert face['kind']=='canonical_dotted_quad_aggregate'
        assert sorted((r[0],r[1]) for r in records)==sorted(map(tuple,face['source_contexts']))
        ks={r[2] for r in records};assert len(ks)==1;k=next(iter(ks))
        alternatives=[r[3] for r in records];keys=set().union(*(set(x) for x in alternatives))
        if any(x!=alternatives[0] for x in alternatives[1:]):varying+=1
        all_provenances+=len(records)
        for t in sorted(keys):
            coefficient=min(x.get(t,0) for x in alternatives)
            if not coefficient:continue
            typ=incidence['base_types'][t]
            if typ['represented_in_lp']:raw[t,k].append([fi,coefficient])
            else:unknown[t,k]+=1
    assert len(byface)==68150 and all_provenances==81183
    columns=[]
    for (t,k),positive in sorted(raw.items()):
        typ=incidence['base_types'][t];sf=typ['old_face'];dimension=typ['simplex_dimension']
        assert 0<=sf<10185
        assert d['faces'][sf]['provenance'][0]['star']==typ['old_star']
        assert d['faces'][sf]['provenance'][0]['simplex_dimension']==dimension
        capacity=comb(dimension+1,k)*(13-dimension);assert capacity>0
        assert all(fi>=10185 and m>0 for fi,m in positive)
        columns.append(dict(index=len(columns),base_type=t,old_star=typ['old_star'],simplex_face=sf,
            prism_dimension=k,simplex_dimension=dimension,vertex_cost=0,
            capacity_per_simplex=capacity,positive_quad_coefficients=positive,
            negative_simplex_coefficient=[sf,-capacity],
            normalization='nonnegative transfer parameter; zero vertex cost; positive minimum quad coefficients and one negative simplex coefficient'))
    # Quad masses do not change. Each transfer has exactly one negative
    # simplex entry, so raising that simplex mass to its largest required
    # load/capacity repairs every transfer without introducing dependencies.
    repaired=dict(alpha);violations=[];requirements=[]
    for col in columns:
        load=sum((m*alpha.get(fi,F(0)) for fi,m in col['positive_quad_coefficients']),F(0))
        sf=col['simplex_face'];capacity=col['capacity_per_simplex'];minimum=load/capacity
        if minimum>repaired.get(sf,F(0)):repaired[sf]=minimum
        old_slack=capacity*alpha.get(sf,F(0))-load
        if old_slack<0:violations.append(col['index'])
        if load:requirements.append([col['index'],str(load),str(minimum),str(old_slack)])
    repaired={i:a for i,a in repaired.items() if a>0}
    assert len(violations)==incidence['counts']['robust_violations']==434
    assert all(repaired.get(fi,F(0))==a for fi,a in alpha.items() if fi>=10185)
    transfer_loads=[]
    for col in columns:
        load=sum((m*repaired.get(fi,F(0)) for fi,m in col['positive_quad_coefficients']),F(0))
        load-=col['capacity_per_simplex']*repaired.get(col['simplex_face'],F(0))
        assert load<=0
        transfer_loads.append(str(load))
    bound,beta,tight,numerator,denominator=D.sparse_dual_bound(d,repaired)
    cert=dict(kind='EXACT_RESTRICTED_MODEL_DUAL',bound=str(bound),objective_numerator=str(numerator),
        normalizer=str(denominator),face_multipliers=[[i,str(a)] for i,a in sorted(repaired.items())],
        vertex_multipliers=[[i,str(a)] for i,a in sorted(beta.items())],
        tight_orbits={str(i):x for i,x in tight.items()},excludes_target12=bound>12)
    B.verify_certificate(d,cert)
    changes=[]
    for sf,a in sorted(repaired.items()):
        previous=alpha.get(sf,F(0))
        if a!=previous:
            assert sf<10185 and a>previous
            changes.append(dict(simplex_face=sf,old_mass=str(previous),new_mass=str(a),increase=str(a-previous),
                                star=d['faces'][sf]['provenance'][0]['star']))
    old_beta={i:F(a) for i,a in old['vertex_multipliers']}
    increases=[[v,str(b-old_beta.get(v,F(0)))] for v,b in beta.items() if b>old_beta.get(v,F(0))]
    extra_demand=numerator-F(old['objective_numerator']);extra_mass=denominator-F(old['normalizer'])
    assert extra_demand>=0 and extra_mass>=0
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'scripts/{P.TAG}.py',
           f'scripts/{I.TAG}.py',f'results/{I.TAG}.json',f'results/{I.LP}.json']
    return dict(status='PASS_EXACT_REPRESENTED_BASE_TRANSFERS_AND_ANALYTIC_DUAL_REPAIR',
        sources={p:P.sha(p) for p in paths},dimension=13,target=12,
        counts=dict(transfer_columns=len(columns),positive_quad_entries=sum(len(c['positive_quad_coefficients']) for c in columns),
            negative_simplex_entries=len(columns),all_quad_rows=len(byface),all_quad_provenances=all_provenances,
            varying_provenance_rows=varying,unknown_base_type_k_pairs=len(unknown),
            old_violated_columns=len(violations),increased_simplex_masses=len(changes),
            newly_positive_simplex_rows=sum(F(x['old_mass'])==0 for x in changes),
            repaired_positive_face_rows=len(repaired),repaired_positive_vertex_masses=len(beta)),
        transfer_columns=columns,provenance_source=f'results/{I.TAG}.json',
        unknown_base_pairs=[[t,k,n] for (t,k),n in sorted(unknown.items())],
        old_bound=old['bound'],old_violated_columns=violations,old_positive_load_requirements=requirements,
        simplex_mass_changes=changes,additional_face_demand=str(extra_demand),additional_vertex_mass=str(extra_mass),
        vertex_mass_increases=increases,repaired_dual=cert,all_transfer_dual_loads=transfer_loads,
        certificate_interpretation='The old columns all have positive vertex costs and are checked by the exact sparse dual verifier. Every added signed transfer has vertex cost0, so its dual load must be nonpositive; all such loads are explicitly recomputed and saved. Quad masses are fixed and each simplex mass is raised independently to the maximum of its old mass and all required base-load/capacity ratios. This is the least componentwise such repair, not an optimized dual.',
        scope='No optimizer. This certificate covers the inherited pair,22clique and6528disjoint-simplex-port columns plus these represented-base conservative transfers only. It does not include the separately studied mixed-three-port columns, whose loads change when simplex masses rise. Missing base simplex types retain unknown mass and have no fabricated zero RHS. A retained bound above12 obstructs only this enlarged partial weighting family; no dimension bound or global realization is claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(B.dump(d)+'\n')
    print(d['status'],d['counts'],d['repaired_dual']['bound'],flush=True)
