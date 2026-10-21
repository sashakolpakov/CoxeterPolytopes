#!/usr/bin/env python3
"""Add raw76 class42946 at three forests disjoint from the two-class patch.

The old six patched graphs are inherited unchanged. Only this class's new
stars are tested; no global optimization or old-pair retest is performed.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import json
import time
import fresh_new5073_raw76_high_mass_local_covers_20261018 as B

Q=B.Q;L=B.L;M=B.M;Z=B.Z;H=B.H;W=B.W;ROOT=B.ROOT
TAG='fresh_new5073_raw76_third_class_local_covers_20261018'
STAGE=Path('/private/tmp')/(TAG+'_stage.json')
FACE=42946


def run(saved=None):
    started=time.monotonic();inp=Q.Inputs();base=inp.read(B.TAG);local=inp.read(L.TAG)
    assert base['represented_classes']==326 and base['represented_forests']==600
    assert base['face_multipliers']==local['face_multipliers']
    alpha={f:F(x) for f,x in base['face_multipliers']};mu={v:F(x) for v,x in base['vertex_multipliers']}
    previous_mu={v:F(x) for v,x in local['vertex_multipliers']}
    selected=next(r for r in base['screen']['rows'] if r['face']==FACE)
    assert F(selected['mass'])==alpha[FACE] and selected['source_context']==59136
    assert selected['spectator_forest']==['A1','A1','A3','H3']
    quads=inp.read(M.QUAD)
    raw=B.normal(M.enc(Q.geometry(inp,[{k:v for k,v in selected.items() if k!='spectator_forest'}],quads,alpha)))
    # The inherited embedding body loops over the passed classes. Its former
    # two-class metadata constants are replaced by direct occurrence counts.
    assert len(raw['selected_faces'])==1 and len(raw['selected_faces'][0]['corners'])==4
    raw['counts']['faces']=1
    raw['counts']['actual_corner_occurrences']=sum(len(r['corners']) for r in raw['selected_faces'])
    targets={r['vertex'] for r in raw['forest_records']};assert targets=={534,788,1665}
    oldpatch={r['vertex']:r for r in base['forest_records']};assert len(oldpatch)==6 and not targets&set(oldpatch)
    assert all(mu.get(v,F())==previous_mu.get(v,F()) for v in targets)
    domains,represented=L.current_domains();p=Z.prepare();patches={r['vertex']:r for r in M.read(Z.TAG)['forest_records']}
    assert FACE not in represented
    assert FACE not in {r['face'] for r in base['selected_faces']}
    localrows={r['vertex']:r for r in local['forest_records']}
    prepared=[B.merge(r,p,localrows,domains,patches,alpha) for r in raw['forest_records']]
    assert sum(u['new_pair_positions'] for u,*_ in prepared)==26910
    assert sum(u['new_literal_stars'] for u,*_ in prepared)==68
    print('RAW76_THIRD_PREPARED',len(prepared),'forests',26910,'newpairs','seconds',round(time.monotonic()-started,2),flush=True)
    paths=set(inp.paths)
    for mod in (B,Q,L,M,Z,H,W,B.A,B.R):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    paths.add(str(Path(__file__).relative_to(ROOT)))
    sources={f:M.digest(f) for f in sorted(paths)}
    if saved:assert sources==saved['sources']
    cache={r['vertex']:r for r in saved['forest_records']} if saved else {}
    if saved:assert set(cache)==targets and len(saved['forest_records'])==3
    out=[];newmu=dict(mu);columns=[];counts=Counter()
    for u,g,owner,pairs,oldedges,demand in prepared:
        v=u['vertex'];physical=Z.physical(u,g,pairs)
        assert {tuple(sorted(e)) for e in g.edges() if max(e)<u['old_literal_stars']}==oldedges
        enum,error=H.vectors_for(g,owner,len(demand));assert error is None,(v,error)
        vectors,reps,ncliques,npatterns=enum;vectors=list(map(list,vectors));reps=list(map(list,reps))
        c=cache[v]['certificate'] if saved else H.exact_profile(vectors,demand)
        value=W.check_certificate(vectors,demand,c);assert value>=F(u['baseline']['exact_local_optimum'])
        newmu[v]=max(mu.get(v,F()),value);coeff=Counter()
        for w,z in zip(map(F,c['weights']),u['corner_orbits']):
            for f,m in z['positive_face_coefficients']:coeff[f]+=w*m
        coeffs=[[f,str(x)] for f,x in sorted(coeff.items()) if x]
        assert sum(alpha.get(f,F())*F(x) for f,x in coeffs)==value
        row=B.normal(M.enc(dict(vertex=v,domain=u,physical=physical,
            graph_edges_sha256=L.jhash(sorted(tuple(sorted(e)) for e in g.edges())),
            maximal_cliques=ncliques,count_patterns=npatterns,maximal_count_vectors=vectors,literal_representatives=reps,
            certificate=c,exact_local_optimum=value,baseline_vertex_mass=mu.get(v,F()),repaired_vertex_mass=newmu[v],
            multiplier_increase=newmu[v]-mu.get(v,F()),positive_face_coefficients=coeffs)))
        if saved:assert row==cache[v]
        out.append(row)
        if value>mu.get(v,F()):columns.append(dict(vertex=v,vertex_cost='1',positive_face_coefficients=coeffs,exact_load=str(value)))
        counts.update(new_pairs=len(pairs),new_exclusions=physical['new_exclusions'],new_unknowns=physical['new_unknown_pairs'],
            new_stars=u['new_literal_stars'],actual_corners=len(u['added_occurrences']),maximal_cliques=ncliques,
            count_vectors=len(vectors),mixture_terms=len(c['mixture']))
        print('RAW76_THIRD_DONE',v,'excluded',physical['new_exclusions'],'unknown',physical['new_unknown_pairs'],
            'vectors',len(vectors),'mu_increase',str(newmu[v]-mu.get(v,F())),'seconds',round(time.monotonic()-started,2),flush=True)
    latest={v:dict(source=L.TAG,row=r) for v,r in localrows.items()}
    for v,r in oldpatch.items():latest[v]=dict(source=B.TAG,row=r)
    inherited=[dict(vertex=v,source=z['source'],record_sha256=L.jhash(z['row'])) for v,z in sorted(latest.items()) if v not in targets]
    assert len(inherited)==597 and all(v in latest for v in targets)
    model=M.read(L.MODEL);loads=[];two=0
    for col in model['prepared']['inherited_several_vertex_columns']+model['prepared']['new_columns']:
        load=sum(alpha.get(f,F())*m for f,m in col['integer_face_coefficients'])
        before=sum(mu.get(v,F())*m for v,m in col['integer_vertex_capacities'])
        after=sum(newmu.get(v,F())*m for v,m in col['integer_vertex_capacities']);assert load<=before<=after
        two+=len(col['integer_vertex_capacities'])>1
        loads.append(dict(variable=col['variable'],all_vertex_costs=col['integer_vertex_capacities'],load=str(load),old_cost=str(before),new_cost=str(after)))
    assert len(loads)==412 and two==14
    num,den=F(base['objective_numerator']),sum(newmu.values(),F())
    result=B.normal(M.enc(dict(status='PASS_NEW5073_RAW76_THIRD_CLASS_EXACT_LOCAL_COVERS',sources=sources,
        selected_face=FACE,raw_domain=raw,forest_records=out,inherited_forests=inherited,
        face_multipliers=base['face_multipliers'],vertex_multipliers=sorted((v,x) for v,x in newmu.items() if x),
        baseline_ratio=base['bound'],objective_numerator=num,normalizer=den,bound=num/den,columns=columns,counts=dict(counts),
        represented_classes=327,represented_forests=600,changed_forests=3,
        previous_patch_vertices=sorted(oldpatch),previous_unknown_pairs_preserved=base['counts']['new_unknowns'],
        checked_recent_model_loads=loads,preserved_multi_owner_columns=14,global_optimizer_invocations=0,new_dimension_bound=None,
        scope='One complete raw76 class added at three forests disjoint from the preceding six. All old graph outcomes,312 preceding unknowns and14 multi-owner costs remain fixed. Complete new-pair tests and three exact graph-cone optima;597 inherited domains/covers. No old-pair retest or global solve. A repaired fixed witness does not decide target feasibility or give a dimension theorem.')))
    print(result['status'],result['counts'],'ratio',float(num/den),flush=True)
    return result


def finalize():
    out=json.loads(STAGE.read_text())
    # Only report metadata may have finished review since arithmetic staging.
    # All computational inputs and the producer itself must be byte-identical.
    for p,h in out['sources'].items():
        now=M.digest(p)
        if p.endswith('.md'):out['sources'][p]=now
        else:assert now==h,p
    path=ROOT/f'results/{TAG}.json';assert not path.exists()
    path.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS_FINALIZED_UNCHANGED_COMPUTATION_WITH_FROZEN_REPORT_BINDINGS',out['bound'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--stage',action='store_true');ap.add_argument('--finalize',action='store_true');a=ap.parse_args()
    assert sum((a.verify,a.stage,a.finalize))<=1
    if a.finalize:finalize()
    else:
        p=ROOT/f'results/{TAG}.json';saved=json.loads(p.read_text()) if a.verify else None
        if not a.verify:assert not (STAGE if a.stage else p).exists()
        out=run(saved)
        if a.verify:assert out==saved
        else:(STAGE if a.stage else p).write_text(json.dumps(out,indent=2)+'\n')
