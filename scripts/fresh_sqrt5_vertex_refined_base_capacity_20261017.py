#!/usr/bin/env python3
"""Actual prism-base capacity refined by full ambient vertex forest.

Every mass uses the frozen uniform provenance split. No optimizer is used.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import fresh_sqrt5_actual_provenance_dual_20261017 as A
import fresh_sqrt5_rank13_simplex_quad_combined_20261017 as P
import fresh_sqrt5_forced_base_dual_repair_20261017 as R
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_vertex_refined_base_capacity_20261017'
QUADS='fresh_sqrt5_rank13_packing_quad_separation_20261017'
SPECTATORS='fresh_sqrt5_rank13_dotted_b1_spectators_20261017'
THEOREM='results/fresh_vertex_weighted_prism_base_capacity_20261017.md'


def run():
    old=P.read(A.TAG);inc=P.read(R.N.I.TAG);quad=P.read(QUADS);sp=P.read(SPECTATORS)
    d,patch,shared=R.domain();alpha={i:F(a) for i,a in old['face_multipliers']}
    oldmu={v:F(a) for v,a in old['vertex_multipliers']}
    oldconstraints={r['base_type']:r for r in old['actual_shared_base_constraints']}
    templates={r['canonical_core']:r for r in inc['prism_templates']}
    base_lookup={(r['core'],tuple(r['spectator_forest'])):r['index'] for r in inc['base_types']}
    assert len(base_lookup)==len(inc['base_types'])
    vlookup={tuple(v['forest']):v['index'] for v in d['vertices']}
    models=dict(B.C.models());finite={B.C.code(g):name for name,g in models.items()}

    @lru_cache(None)
    def core_corners(ci):
        core=inc['cores'][ci];g={i:{} for i in range(core['nodes'])}
        for i,j,(m,a,b) in core['edges']:g[i][j]=g[j][i]=m
        answer={}
        for p in core['spherical_deletions']:
            h={i:{j:m for j,m in row.items() if j!=p} for i,row in g.items() if i!=p}
            assert all(m>0 for row in h.values() for m in row.values())
            pieces=B.C.pieces(h)
            assert all(sum(map(len,c.values()))==2*len(c)-2 for c in pieces)
            answer[p]=tuple(sorted(finite[B.C.code(c)] for c in pieces))
        assert len(answer)==core['simplex_dimension']+1
        return answer

    profiles={};touched=set()
    for t,oldc in sorted(oldconstraints.items()):
        typ=inc['base_types'][t];dim=typ['simplex_dimension'];corners=[];counts=Counter()
        for p,coreforest in sorted(core_corners(typ['core']).items()):
            forest=tuple(sorted(list(coreforest)+typ['spectator_forest']));v=vlookup[forest]
            assert sum(len(models[n]) for n in forest)==13
            corners.append([p,v]);counts[v]+=1;touched.add(v)
        assert sum(counts.values())==dim+1 and oldc['simplex_dimension']==dim
        profiles[t]=dict(base_type=t,simplex_face=oldc['simplex_face'],simplex_dimension=dim,
            core=typ['core'],spectator_forest=typ['spectator_forest'],corners=corners,
            vertex_counts=[list(x) for x in sorted(counts.items())])
    loads=defaultdict(F);coeff=defaultdict(Counter);literal=[];incidence_count=0;base_count=0
    dim_pairs=Counter();exact_zero_support_checks=0
    for provenance in old['provenance_masses']:
        j=provenance['index'];a=F(provenance['mass']);ci=provenance['canonical_prism_core']
        template=templates[ci];qc=quad['canonical_cores'][ci];k=provenance['prism_dimension']
        assert qc['index']==ci and template['prism_dimension']==k
        spectator=sp['contexts'][provenance['source_context']]['spectator_forest']
        vertices={}
        for rec in qc['aggregate_corner_occurrences']:
            pair=tuple(sorted(rec['vertex_omitted']));forest=tuple(sorted(rec['forest']+spectator))
            if pair in vertices:assert vertices[pair]==forest
            vertices[pair]=forest
        # A star can have additional PD deletion edges outside its unique
        # quadrilateral prism (they can extend a base to a larger simplex).
        # The prism vertices are exactly the union of its complete K2,2
        # partitions, not every PD deletion of the larger wall star.
        prism_vertices={tuple(sorted((u,v))) for left,right in qc['quad_partitions'] for u in left for v in right}
        assert set(vertices)==prism_vertices
        assert set(vertices)<={tuple(sorted(x)) for x in qc['pd_deletions']}
        assert len(vertices)==2*k
        pc=Counter();base_rows=[]
        for base in template['bases']:
            b=base['deleted_normal'];ciS=base['base_core'];extras=base['finite_extra']
            forestS=tuple(sorted(spectator+extras));t=base_lookup[ciS,forestS]
            profile=profiles[t];counts=dict(profile['vertex_counts']);isom=dict(base['core_isomorphism'])
            assert profile['simplex_dimension']==base['simplex_dimension']
            Svertices=dict(profile['corners']);fullS=set(base['full_spherical_deletions'])
            assert {isom[p] for p in fullS}==set(Svertices)
            rows=[]
            for pair,forest in sorted(vertices.items()):
                if b not in pair:continue
                p=pair[1] if pair[0]==b else pair[0]
                assert p in fullS and p in isom
                v=vlookup[forest];pS=isom[p]
                # An actual prism-base vertex is literally a vertex of S,
                # with exactly the same full ambient spherical forest.
                assert Svertices[pS]==v and counts.get(v,0)>0
                exact_zero_support_checks+=1
                rows.append([list(pair),pS,v]);loads[t,v]+=a;coeff[t,v][j]+=1
            assert len(rows)==k
            pc[t]+=k;incidence_count+=len(rows);base_count+=1
            dim_pairs[k,profile['simplex_dimension']]+=1
            base_rows.append(dict(deleted_normal=b,opposite_base_facet=base['opposite_base_facet'],
                                  base_type=t,vertices=rows))
        assert dict(pc)==dict(provenance['base_coefficients'])
        assert len(base_rows)==(4 if k==2 else 2)
        literal.append(dict(provenance=j,source_context=provenance['source_context'],
                            canonical_prism_core=ci,bases=base_rows))
    assert len(literal)==1858
    constraints=[];violations=[];needs={};bybase=defaultdict(list)
    for t,profile in sorted(profiles.items()):
        dim=profile['simplex_dimension'];fi=profile['simplex_face'];need=F(0)
        assert sum((loads[t,v] for v,m in profile['vertex_counts']),F(0))==F(oldconstraints[t]['actual_positive_load'])
        for v,m in profile['vertex_counts']:
            capacity=(13-dim)*m;load=loads[t,v];oldmass=alpha[fi];required=load/capacity
            rec=dict(index=len(constraints),base_type=t,simplex_face=fi,simplex_dimension=dim,
                vertex=v,simplex_vertex_count=m,capacity=capacity,
                provenance_coefficients=[list(x) for x in sorted(coeff[t,v].items())],
                actual_load=str(load),old_simplex_mass=str(oldmass),old_signed_load=str(load-capacity*oldmass),
                required_simplex_mass=str(required))
            constraints.append(rec);bybase[t].append(rec['index'])
            if load>capacity*oldmass:violations.append(rec['index'])
            need=max(need,required)
        needs[fi]=max(alpha[fi],need)
    assert set(loads)<={(t,v) for t,p in profiles.items() for v,m in p['vertex_counts']}
    repaired=dict(alpha);changes=[]
    for fi,a in sorted(needs.items()):
        if a>alpha[fi]:changes.append([fi,str(alpha[fi]),str(a),str(a-alpha[fi])])
        repaired[fi]=a
    assert set(repaired)==set(alpha)
    for rec in constraints:
        rec['repaired_simplex_mass']=str(repaired[rec['simplex_face']])
        rec['repaired_signed_load']=str(F(rec['actual_load'])-rec['capacity']*repaired[rec['simplex_face']])
        assert F(rec['repaired_signed_load'])<=0
    print('VERTEX_REFINED',len(constraints),'VIOLATIONS',len(violations),'RAISED',len(changes),flush=True)

    _,mu,tight,num,_=R.D.sparse_dual_bound(d,repaired)
    features={r['face']:r for r in old['positive_face_mixed_features']}
    assert set(features)==set(repaired)
    columns=old['triple_orbit_columns'];mixed_loads=defaultdict(F)
    for fi,a in repaired.items():
        for j,m in features[fi]['coefficients']:mixed_loads[j]+=a*m
    mixed_changes=[]
    for c in columns:
        j,v,cost=c['index'],c['vertex'],c['vertex_cost'];minimum=mixed_loads[j]/cost
        if minimum>mu.get(v,F(0)):
            mixed_changes.append([j,v,str(mu.get(v,F(0))),str(minimum)]);mu[v]=minimum
    assert set(mu)==set(oldmu)
    positive_loads=defaultdict(F)
    for fi,a in repaired.items():
        for j,m in d['faces'][fi]['coefficients']:positive_loads[j]+=a*m
    for c in d['variables']:assert positive_loads[c['index']]<=c['multiplicity']*mu.get(c['vertex'],F(0))
    for c in columns:assert mixed_loads[c['index']]<=c['vertex_cost']*mu[c['vertex']]
    inherited_signed=R.signed_loads(d,repaired);assert all(a<=0 for a in inherited_signed)
    actual_shared=[]
    for t,rec in sorted(oldconstraints.items()):
        value=F(rec['actual_positive_load'])-rec['capacity']*repaired[rec['simplex_face']]
        assert value<=0;actual_shared.append([t,str(value)])
    for rec in constraints:
        load=sum((m*F(old['provenance_masses'][j]['mass']) for j,m in rec['provenance_coefficients']),F(0))
        assert load==F(rec['actual_load']) and load<=rec['capacity']*repaired[rec['simplex_face']]
    den=sum(mu.values(),F(0));bound=num/den
    ranked=sorted(violations,key=lambda i:F(constraints[i]['required_simplex_mass'])/F(constraints[i]['old_simplex_mass']),reverse=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',f'results/{A.TAG}.md',
        f'scripts/{R.TAG}.py',f'scripts/{R.D.TAG}.py',f'scripts/{B.C.TAG}.py',f'results/{R.N.I.TAG}.json',
        f'results/{QUADS}.json',f'results/{SPECTATORS}.json',THEOREM]
    return dict(status='PASS_EXACT_VERTEX_REFINED_BASE_DUAL_BARRIER' if bound>12 else 'PASS_VERTEX_REFINED_REPAIR_AT_MOST12',
        sources={p:P.sha(p) for p in paths},dimension=13,target=12,excludes_target12=bound>12,
        bound=str(bound),old_uniform_bound=old['bound'],objective_numerator=str(num),normalizer=str(den),
        counts=dict(actual_provenances=len(literal),actual_base_occurrences=base_count,
            actual_base_vertex_occurrences=incidence_count,exact_vertex_type_consistency_checks=exact_zero_support_checks,
            simplex_profiles=len(profiles),labelled_constraints=len(constraints),old_violations=len(violations),
            raised_simplex_rows=len(changes),newly_positive_rows=0,positive_rows=len(repaired),
            positive_forests=len(mu),positive_cost_columns=len(d['variables']),mixed_columns=len(columns),
            mixed_vertex_increases=len({r[1] for r in mixed_changes})),
        full_vertex_forests=[[v,d['vertices'][v]['forest']] for v in sorted(touched)],
        simplex_vertex_profiles=list(profiles.values()),actual_provenance_base_vertices=literal,
        labelled_transfer_columns=constraints,old_violating_columns=violations,
        worst_old_violation_column=ranked[0] if ranked else None,
        old_violation_columns_by_factor=ranked,
        simplex_mass_changes=changes,face_multipliers=[[i,str(a)] for i,a in sorted(repaired.items())],
        vertex_multipliers=[[v,str(a)] for v,a in sorted(mu.items())],
        positive_cost_tight_columns={str(v):j for v,j in tight.items()},
        mixed_column_loads=[[j,str(a)] for j,a in sorted(mixed_loads.items()) if a],mixed_vertex_changes=mixed_changes,
        inherited_signed_loads=list(map(str,inherited_signed)),actual_shared_signed_loads=actual_shared,
        prism_simplex_dimension_pairs=[[k,dim,n] for (k,dim),n in sorted(dim_pairs.items())],
        zero_support_rule='Every actual prism-base vertex is mapped by the saved literal core isomorphism to a spherical deletion of its containing maximal simplex. Both reconstructions must have the identical full ambient finite forest. Thus an absent simplex label has exactly zero incidence coefficient, not an omitted obligation.',
        repair_definition='Keep every frozen actual provenance mass and merged quad total. For each existing simplex row, retain its old mass or increase it to the maximum of all labelled load/capacity ratios. Recompute every positive-cost column, all mixed triples, old shared transfers, actual unlabelled transfers and all new labelled transfers. No optimizer.',
        scope='Full ambient vertex-forest labels refine the proved shared incidence capacity on the same actual positive support. All labels use one consistent frozen provenance distribution. No new face search or common arithmetic-form realization is asserted. Bound>12 is a partial-model obstruction only; a failed analytic repair would not establish feasibility or a dimension bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(B.dump(d)+'\n')
    print(d['status'],d['counts'],d['bound'],flush=True)
