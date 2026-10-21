#!/usr/bin/env python3
"""Conservative mixed three-port separation retaining every merged provenance."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C
import fresh_mixed_face_port_incidence_20261017 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_mixed_three_port_dual_20261017'
MODEL='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'
SEP='fresh_sqrt5_rank13_packing_quad_separation_20261017'
CORE='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017'
SP='fresh_sqrt5_rank13_dotted_b1_spectators_20261017'
TRI='fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017'


def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for p,h in d.get('sources',{}).items():assert sha(p)==h,(tag,p)
    return d
def freeze(x):return tuple(map(freeze,x)) if isinstance(x,list) else x
def code(pieces,marks):return tuple(sorted(C.code(p,marks) for p in pieces))
def pair(a,b):return tuple(sorted((a,b)))


def run():
    model=read(MODEL);sep=read(SEP);cores=read(CORE)['records']
    spectators=read(SP)['contexts'];stars=read(TRI)['stars']
    dual=next(c for c in model['certificates'] if c['kind']=='EXACT_RESTRICTED_MODEL_DUAL')
    alpha={i:F(a) for i,a in dual['face_multipliers']}
    mu={i:F(a) for i,a in dual['vertex_multipliers']}
    oldloads=defaultdict(F)
    for fi,a in alpha.items():
        for j,n in model['faces'][fi]['coefficients']:oldloads[j]+=a*n
    oldminimum=defaultdict(F)
    for j,value in oldloads.items():
        var=model['variables'][j];v=var['vertex']
        oldminimum[v]=max(oldminimum[v],value/var['multiplicity'])
    assert dict(oldminimum)==mu
    models=dict(C.models());finite={C.code(g):name for name,g in models.items()}
    vlookup={tuple(v['forest']):v['index'] for v in model['vertices']}
    pairlookup={(v['vertex'],freeze(v['marked_forest_code'])):v['index']
                for v in model['variables'][:model['counts']['pair_variables']]}
    orbit={};columns=[]
    for v in sorted(mu):
        g={}
        for name in model['vertices'][v]['forest']:
            offset=len(g)
            for i,row in models[name].items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==13
        pieces=C.pieces(g);groups=defaultdict(list)
        for U in combinations(range(13),3):groups[code(pieces,U)].append(U)
        assert sum(map(len,groups.values()))==286
        for key,Us in sorted(groups.items()):
            j=len(columns);orbit[v,key]=j
            columns.append(dict(index=j,vertex=v,forest=model['vertices'][v]['forest'],
                marked_forest_code=key,vertex_cost=len(Us),literal_triple_orbit=[list(U) for U in Us]))

    def corner(g,footprint,cc,pp):
        assert len(g)==13 and all(m>0 for row in g.values() for m in row.values())
        pieces=C.pieces(g);forest=tuple(sorted(finite[C.code(p)] for p in pieces))
        v=vlookup[forest];assert v in mu
        footprint=set(map(tuple,footprint))
        for a,b in footprint:
            assert a in g and b in g
            pp[pairlookup[v,code(pieces,(a,b))]]+=1
        for U in M.triple_indicators(footprint):cc[orbit[v,code(pieces,U)]]+=1

    def simplex(face):
        star=stars[face['provenance'][0]['star']]
        g=C.graph(range(14),{(i,j):m for i,j,m in star['edges']})
        S={r['omitted'] for r in star['spherical_deletions']}
        cc=Counter();pp=Counter()
        for o in S:
            h={i:{j:m for j,m in row.items() if j!=o} for i,row in g.items() if i!=o}
            corner(h,list(combinations(sorted(S-{o}),2)),cc,pp)
        return cc,pp

    templates={}
    for core in sep['canonical_cores']:
        source=cores[core['source_core']];n=source['nodes']
        g=C.graph(range(n),dict(zip(combinations(range(n),2),source['row']['labels'])))
        # C.graph drops right pairs but retains dotted label0.
        at=defaultdict(set)
        for row in core['aggregate_corner_occurrences']:
            om=tuple(sorted(row['vertex_omitted']));e=tuple(sorted(row['marked_corner_roots']))
            assert e not in at[om];at[om].add(e)
        assert sum(map(len,at.values()))==4*core['required']
        templates[core['index']]=n,g,at

    def quad(ctxid,ci):
        ctx=spectators[ctxid];assert ctx['index']==ctxid
        assert sep['source_core_maps'][ctx['core']]['canonical_core']==ci
        n,g,at=templates[ci]
        full={i:dict(row) for i,row in g.items()};offset=n
        for name in ctx['spectator_forest']:
            for i,row in models[name].items():full[offset+i]={offset+j:m for j,m in row.items()}
            offset+=len(models[name])
        assert len(full)==15
        cc=Counter();pp=Counter()
        for omitted,footprint in sorted(at.items()):
            h={i:{j:m for j,m in row.items() if j not in omitted}
               for i,row in full.items() if i not in omitted}
            corner(h,footprint,cc,pp)
        return cc,pp

    loads=defaultdict(F);origins=defaultdict(list)
    records=[];counts=Counter();trace=hashlib.sha256()
    for serial,(fi,a) in enumerate(sorted(alpha.items())):
        face=model['faces'][fi];assert face['index']==fi
        expected={j:n for j,n in face['coefficients'] if j<model['counts']['pair_variables']}
        if face['kind']=='maximal_simplex_aggregate':
            cc,pp=simplex(face);assert dict(pp)==expected
            candidates=[cc];counts['positive_simplex_rows']+=1
        else:
            assert face['kind']=='canonical_dotted_quad_aggregate'
            candidates=[]
            for ctx,ci in face['source_contexts']:
                cc,pp=quad(ctx,ci);assert dict(pp)==expected,(fi,ctx)
                candidates.append(cc)
            counts['positive_quad_rows']+=1
            counts['positive_quad_provenances']+=len(candidates)
        keys=set().union(*(set(c) for c in candidates))
        lower={j:min(c.get(j,0) for c in candidates) for j in keys}
        upper={j:max(c.get(j,0) for c in candidates) for j in keys}
        different=lower!=upper
        if different:counts['rows_with_provenance_dependent_mixed_coefficients']+=1
        lower={j:n for j,n in lower.items() if n}
        counts['positive_row_mixed_nonzeros']+=len(lower)
        for j,n in lower.items():
            loads[j]+=a*n;origins[j].append([fi,n,str(a)])
        rec=dict(face=fi,provenances=len(candidates),coefficients=[[j,n] for j,n in sorted(lower.items())],
                 provenance_dependent=different)
        records.append(rec)
        trace.update((json.dumps([fi,sorted(lower.items()),sorted(upper.items())],separators=(',',':'))+'\n').encode())
        if serial%300==0:print('MIXED_POSITIVE_ROWS',serial,'/',len(alpha),flush=True)
    violations=[];best={}
    for col in columns:
        j,v,cost=col['index'],col['vertex'],col['vertex_cost']
        necessary=loads[j]/cost
        if necessary>best.get(v,(F(),None))[0]:best[v]=necessary,j
        if necessary>mu[v]:
            violations.append(dict(**col,total_incidence=str(loads[j]),vertex_mass=str(mu[v]),
                required_vertex_mass=str(necessary),ratio=str(necessary/mu[v]),positive_rows=origins[j]))
    extra=sum((max(F(),value-mu[v]) for v,(value,j) in best.items()),F())
    newmu={v:max(mu[v],best.get(v,(F(),None))[0]) for v in mu}
    assert sum(newmu.values(),F())==sum(mu.values(),F())+extra
    for col in columns:assert loads[col['index']]<=newmu[col['vertex']]*col['vertex_cost']
    # All positive merged rows happen to have identical complete mixed vectors.
    # Hence the resulting dual also lifts to the model with every full star split
    # into its own row, by assigning the old mass to one of its provenances.
    assert counts.get('rows_with_provenance_dependent_mixed_coefficients',0)==0
    numerator=sum((a*model['faces'][i]['required'] for i,a in alpha.items()),F())
    normalizer=sum(mu.values(),F())
    assert numerator==F(dual['objective_numerator']) and normalizer==F(dual['normalizer'])
    assert numerator/normalizer==F(dual['bound'])
    counts.update(positive_rows=len(alpha),positive_forests=len(mu),triple_orbit_columns=len(columns),
                  violating_columns=len(violations),violating_forests=len({v['vertex'] for v in violations}),
                  rows_with_provenance_dependent_mixed_coefficients=0)
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'scripts/{M.TAG}.py',f'results/{M.TAG}.md']
    paths += [f'results/{tag}.json' for tag in [MODEL,SEP,CORE,SP,TRI]]
    assert numerator/(normalizer+extra)>12
    return dict(status='PASS_EXACT_ALL_THREE_PORT_MIXED_CHARGE_OBSTRUCTION',
        sources={p:sha(p) for p in paths},counts=dict(counts),
        provenance_policy='For each inherited merged quad row, independently reconstruct EVERY source context and its complete occupied port-pair footprints. Verify its full pair coefficients against that row. Each new mixed column uses the minimum integer coefficient over all those provenances. This weakens the row and remains valid for each actual star; no first-provenance assumption is made.',
        vertex_cost_policy='Sum the budget-one mixed indicator over all triples in a labelled forest orbit. The integer column cost is the orbit size, and the face coefficient is the number of selected triples at its actual aggregate vertices.',
        old_dual_bound=str(numerator/normalizer),additional_vertex_mass_with_fixed_face_multipliers=str(extra),
        old_demand_over_necessary_new_vertex_mass=str(numerator/(normalizer+extra)),
        augmented_dual_vertex_multipliers=[[v,str(value)] for v,value in sorted(newmu.items())],
        split_positive_quad_mass=[[r['face'],model['faces'][r['face']]['source_contexts'][0][0],str(alpha[r['face']])]
            for r in records if model['faces'][r['face']]['kind']=='canonical_dotted_quad_aggregate'],
        positive_row_feature_trace_sha256=trace.hexdigest(),positive_rows=records,
        violations=violations,
        strongest_violation=max(violations,key=lambda r:F(r['ratio'])) if violations else None,
        scope='Exact dual lower bound for the saved combined domain even after adding ALL shared three-port indicator charges. All positive-forest triple orbits are checked; every other forest has zero face load. The positive merged rows have identical mixed coefficient vectors across EVERY provenance, so the dual also lifts to distinct full-star rows. No optimizer, primal rule, compact realization or dimension theorem; full face coverage remains absent.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['old_demand_over_necessary_new_vertex_mass'],flush=True)
