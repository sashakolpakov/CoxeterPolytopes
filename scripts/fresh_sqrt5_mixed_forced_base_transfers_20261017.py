#!/usr/bin/env python3
"""All mixed triple indicators together with the forced-base analytic dual."""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_simplex_quad_combined_20261017 as P
import fresh_sqrt5_forced_base_dual_repair_20261017 as R
import fresh_mixed_face_port_incidence_20261017 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_mixed_forced_base_transfers_20261017'
OLD_MIXED='fresh_sqrt5_mixed_shared_transfers_20261017'


def code(pieces,marks=()):return tuple(sorted(B.C.code(p,marks) for p in pieces))


def run():
    old_mixed=P.read(OLD_MIXED);saved=P.read(R.TAG)
    d,patch,shared=R.domain();incidence=P.read(R.N.I.TAG)
    dual=saved['repaired_dual'];R.verify_dual(d,dual)
    alpha={i:F(a) for i,a in dual['face_multipliers']}
    old_alpha={i:F(a) for i,a in old_mixed['face_multipliers']}
    assert all(alpha[i]==a for i,a in old_alpha.items())
    assert {i for i in alpha if i<patch['inherited_counts']['faces']}==set(old_alpha)
    mu={i:F(a) for i,a in dual['vertex_multipliers']}
    old_mu={i:F(a) for i,a in shared['repaired_dual']['vertex_multipliers']}
    models=dict(B.C.models());finite={B.C.code(g):name for name,g in models.items()}
    vlookup={tuple(v['forest']):v['index'] for v in d['vertices']}
    pairlookup={(v['vertex'],B.dump(v['marked_forest_code'])):v['index']
                for v in d['variables'][:203565]}
    portlookup={(v['vertex'],B.dump(v['marked_forest_code'])):v['index']
                for v in d['variables'] if v.get('kind')=='maximal_simplex_disjoint_port_charge'}
    assert not patch['new_pair_variable_indices']
    forests=set(mu)|{r['vertex'] for r in patch['touched_vertex_budgets']}
    columns=[];orbit={};byvertex={}
    for v in sorted(forests):
        g={}
        for name in d['vertices'][v]['forest']:
            offset=len(g)
            for i,row in models[name].items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==13;pieces=B.C.pieces(g);groups=Counter()
        for U in combinations(range(13),3):groups[code(pieces,U)]+=1
        assert sum(groups.values())==286
        byvertex[v]=[]
        for key,cost in sorted(groups.items()):
            j=len(columns);orbit[v,key]=j;byvertex[v].append(j)
            columns.append(dict(index=j,vertex=v,marked_forest_code=key,vertex_cost=cost))
    # The frozen old feature exporter uses exactly this sorted (vertex,code)
    # order. Reconstruct its complete orbit map, rather than recycling old
    # numeric indices after adding positive forests.
    old_to_new=[]
    for v in sorted(old_mu):old_to_new.extend(byvertex[v])
    assert len(old_to_new)==old_mixed['counts']['triple_orbit_columns']==52736
    features=[];old_provenances=0;old_provenance_dependent=0
    assert {r['face'] for r in old_mixed['positive_rows']}==set(old_alpha)
    for row in old_mixed['positive_rows']:
        fi=row['face'];face=d['faces'][fi]
        if face['kind']=='canonical_dotted_quad_aggregate':
            assert row['provenances']==len(face['source_contexts'])
            old_provenances+=row['provenances']
        else:assert row['provenances']==1 and face['kind']=='maximal_simplex_aggregate'
        coeff=[[old_to_new[j],n] for j,n in row['coefficients']]
        assert all(n>0 for j,n in coeff) and len({j for j,n in coeff})==len(coeff)
        features.append(dict(face=fi,source='inherited_all_provenance_feature',coefficients=coeff,
            provenances=row['provenances'],provenance_dependent=row['provenance_dependent']))
        old_provenance_dependent+=row['provenance_dependent']
    assert old_provenances==1858 and old_provenance_dependent==0
    new_dimensions=Counter();zero_edge_features=0;new_corner_count=0
    for serial,face in enumerate(patch['new_faces']):
        ci=face['core'];core=incidence['cores'][ci]
        g={i:{} for i in range(core['nodes'])}
        for i,j,(m,a,b) in core['edges']:g[i][j]=g[j][i]=m
        offset=core['nodes']
        for name in face['spectator_forest']:
            for i,row in models[name].items():g[offset+i]={offset+j:m for j,m in row.items()}
            offset+=len(models[name])
        assert len(g)==14 and offset==14
        S=set(core['spherical_deletions']);dimension=face['simplex_dimension']
        assert len(S)==dimension+1
        cc=Counter();pp=Counter();uu=Counter();corner_rows=[]
        saved_corners={c['omitted']:c for c in face['corners']};assert set(saved_corners)==S
        for omitted in sorted(S):
            h={i:{j:m for j,m in row.items() if j!=omitted} for i,row in g.items() if i!=omitted}
            assert len(h)==13 and all(m>0 for row in h.values() for m in row.values())
            pieces=B.C.pieces(h);forest=tuple(sorted(finite[B.C.code(p)] for p in pieces))
            v=vlookup[forest];assert v in forests
            block=sorted(S-{omitted});footprint=list(combinations(block,2))
            expected=saved_corners[omitted]
            assert expected['vertex']==v and expected['forest']==list(forest) and expected['literal_port_block']==block
            cp=Counter();cu=Counter();cm=Counter()
            for pair in footprint:cp[pairlookup[v,B.dump(code(pieces,pair))]]+=1
            for port in block:cu[portlookup[v,B.dump(code(pieces,(port,)))]]+=1
            indicators=M.triple_indicators(footprint)
            assert len(indicators)==comb(dimension,3)
            for U in indicators:cm[orbit[v,code(pieces,U)]]+=1
            assert sorted(cp.items())==list(map(tuple,expected['pair_coefficients']))
            assert sorted(cu.items())==list(map(tuple,expected['port_coefficients']))
            pp.update(cp);uu.update(cu);cc.update(cm)
            corner_rows.append([omitted,v,[list(x) for x in sorted(cm.items())]])
        assert sorted(pp.items())==list(map(tuple,face['pair_coefficients']))
        assert sorted(uu.items())==list(map(tuple,face['port_coefficients']))
        assert sum(cc.values())==(dimension+1)*comb(dimension,3)
        if dimension==1:assert not cc;zero_edge_features+=1
        new_dimensions[dimension]+=1;new_corner_count+=len(S)
        features.append(dict(face=face['index'],source='reconstructed_forced_base',base_type=face['base_type'],
            simplex_dimension=dimension,coefficients=[list(x) for x in sorted(cc.items())],
            corner_features=corner_rows,provenances=1,provenance_dependent=False))
        if serial%1000==0:print('MIXED_FORCED_BASE_STARS',serial,'/',len(patch['new_faces']),flush=True)
    assert zero_edge_features==3460 and new_corner_count==9496
    assert len({r['face'] for r in features})==len(features)
    lookup={r['face']:r for r in features};assert set(alpha)<=set(lookup)
    loads=defaultdict(F);origins=defaultdict(list)
    for fi,a in sorted(alpha.items()):
        for j,n in lookup[fi]['coefficients']:
            loads[j]+=a*n;origins[j].append([fi,n,str(a)])
    assert all(columns[j]['vertex'] in mu for j in loads)
    best=defaultdict(F);violations=[]
    for c in columns:
        j,v,cost=c['index'],c['vertex'],c['vertex_cost'];minimum=loads[j]/cost
        best[v]=max(best[v],minimum)
        if minimum>mu.get(v,F(0)):
            violations.append(dict(column=j,vertex=v,forest=d['vertices'][v]['forest'],
                marked_forest_code=c['marked_forest_code'],vertex_cost=cost,
                mixed_load=str(loads[j]),old_vertex_mass=str(mu.get(v,F(0))),required_mass=str(minimum),
                positive_rows=origins[j]))
    new_mu={v:max(mu.get(v,F(0)),best[v]) for v in sorted(forests)}
    new_mu={v:a for v,a in new_mu.items() if a>0}
    # Recount every inherited/new positive-cost column with the fixed masses.
    positive_loads=defaultdict(F)
    for fi,a in alpha.items():
        for j,n in d['faces'][fi]['coefficients']:positive_loads[j]+=a*n
    for v in d['variables']:
        assert positive_loads[v['index']]<=v['multiplicity']*new_mu.get(v['vertex'],F(0))
    for c in columns:assert loads[c['index']]<=c['vertex_cost']*new_mu.get(c['vertex'],F(0))
    signed=R.signed_loads(d,alpha);assert all(a<=0 for a in signed)
    numerator=sum((d['faces'][i]['required']*a for i,a in alpha.items()),F(0))
    assert numerator==F(dual['objective_numerator'])
    old_den=sum(mu.values(),F(0));assert old_den==F(dual['normalizer'])
    denominator=sum(new_mu.values(),F(0));bound=numerator/denominator
    changes=[[v,str(mu.get(v,F(0))),str(a),str(a-mu.get(v,F(0)))]
             for v,a in new_mu.items() if a>mu.get(v,F(0))]
    paths=[f'scripts/{TAG}.py',f'scripts/{B.C.TAG}.py',f'scripts/{P.TAG}.py',f'scripts/{M.TAG}.py',
        f'results/{M.TAG}.md',f'scripts/{R.TAG}.py',f'results/{R.TAG}.json',f'results/{R.N.TAG}.json',
        f'results/{R.N.I.TAG}.json',f'results/{R.N.I.LP}.json',f'results/{R.S.TAG}.json',
        f'scripts/{OLD_MIXED}.py',f'results/{OLD_MIXED}.json']
    return dict(status='PASS_EXACT_MIXED_FORCED_BASE_TRANSFER_DUAL_BARRIER' if bound>12 else
        'PASS_EXACT_MIXED_FORCED_BASE_REPAIR_AT_MOST12',sources={p:P.sha(p) for p in paths},
        dimension=13,target=12,excludes_target12=bound>12,
        counts=dict(positive_cost_columns=len(d['variables']),signed_zero_cost_columns=len(signed),
            positive_face_rows=len(alpha),positive_vertex_forests=len(mu),
            exported_triple_forests=len(forests),extra_zero_load_forests=len(forests-set(mu)),
            triple_orbit_columns=len(columns),positive_forest_triple_orbits=sum(len(byvertex[v]) for v in mu),
            inherited_mixed_columns=len(old_to_new),inherited_positive_row_features=len(old_mixed['positive_rows']),
            new_full_star_features=len(patch['new_faces']),new_corner_occurrences=new_corner_count,
            zero_mixed_maximal_edges=zero_edge_features,positive_quad_provenances=old_provenances,
            old_provenance_dependent_rows=old_provenance_dependent,violating_mixed_columns=len(violations),
            increased_vertex_masses=len(changes),augmented_positive_vertices=len(new_mu)),
        old_forced_base_bound=dual['bound'],bound=str(bound),objective_numerator=str(numerator),normalizer=str(denominator),
        additional_vertex_mass=str(denominator-old_den),face_multipliers=dual['face_multipliers'],
        augmented_vertex_multipliers=[[v,str(a)] for v,a in sorted(new_mu.items())],
        all_signed_transfer_loads=list(map(str,signed)),vertex_mass_changes=changes,
        triple_orbit_columns=columns,old_to_new_mixed_column=old_to_new,
        face_features=features,new_simplex_dimension_counts=dict(sorted(new_dimensions.items())),
        mixed_column_loads=[[j,str(a)] for j,a in sorted(loads.items()) if a],violations=violations,
        provenance_policy='The frozen old feature certificate reconstructs every merged quad provenance, retaining coefficientwise minima; its complete old orbit map is reconstructed and remapped. All old face masses are unchanged. Each new full base star is separate and has its full literal corner footprint rebuilt. No representative is chosen for a merged quad, and no common actual provenance splitting for mixed and signed-transfer minima is asserted.',
        edge_policy='A maximal edge has one port at each endpoint, hence an empty occupied pair footprint and zero mixed-three-port coefficient. Its one-port cost and zero-demand transfer row remain active and are rechecked. A triangle also has zero mixed-three-port indicator at each corner; tetrahedra and four-simplexes contribute their full three-port cliques.',
        feature_scope='Every one of the4195 new base stars is exported, including zero-mass stars, and every old positive row uses the complete frozen all-provenance feature. Old zero-mass row features need not be exported for a dual certificate, but would need full reconstruction before an augmented primal optimization. All positive-forest triple orbits are tested; the27 extra exported forests have zero load.',
        scope='Exact fixed-alpha dual repair for the union of all current positive-cost columns, all mixed-three-port indicators, and all represented/new shared base transfers. Other absent face/base/field classes remain outside. A bound>12 is only a partial weighting-family obstruction. No optimizer, dimension theorem, compact realization, or simultaneous provenance splitting is claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(B.dump(d)+'\n')
    print(d['status'],d['counts'],d['bound'],flush=True)
