#!/usr/bin/env python3
"""Fresh exact local covers at the5121 dual masses on571 frozen graph cones.

Complete count-vector domains are inherited; only bounded local LPs propose
new rational certificates. Replay invokes no solver. No global optimization.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import time
import fresh_cross_vertex_support_migration_20261018 as M
import fresh_current_historical_weighted_corners_20261017 as H

ROOT=M.ROOT;W=H.W
TAG='fresh_cross_vertex_current290_local_covers_20261018'
PROGRESS=Path('/private/tmp')/(TAG+'_progress.jsonl')


def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def run(saved=None):
    start=time.monotonic();cert=M.dual(M.MODEL)
    alpha={i:F(x) for i,x in cert['face_multipliers']};mu={v:F(x) for v,x in cert['vertex_multipliers']}
    assert len(alpha)==5121
    domains=M.current_domains();assert len(domains)==571
    stems=tuple(dict.fromkeys((M.MODEL,*M.DOMAIN_STEMS)))
    sources=M.bindings(stems)
    for mod in (M,H,W):
        p=str(Path(mod.__file__).relative_to(ROOT));sources[p]=M.digest(p)
    sources[f'scripts/{TAG}.py']=M.digest(f'scripts/{TAG}.py')
    key=jhash([sources,cert['face_multipliers']]);cached={}
    if saved is not None:
        assert saved['sources']==sources and saved['face_multipliers']==cert['face_multipliers']
        cached={r['vertex']:r for r in saved['forest_records']}
    elif PROGRESS.exists():
        for line in PROGRESS.read_text().splitlines():
            r=json.loads(line);assert r['key']==key and r['row']['vertex'] not in cached
            cached[r['row']['vertex']]=r['row']
    order=sorted(domains,key=lambda v:(-mu.get(v,F()),v))
    records=[];columns=[];newmu=dict(mu);zero=0;solves=0
    for step,v in enumerate(order,1):
        d=domains[v];old=d['row'];groups=d['groups'];vectors=old['maximal_count_vectors']
        demand=[sum((alpha.get(f,F())*F(m) for f,m in g['positive_face_coefficients']),F()) for g in groups]
        assert len(groups)==len(demand)==len(vectors[0])
        if v in cached:c=cached[v]['certificate']
        else:
            assert saved is None
            c=H.exact_profile(vectors,demand);solves+=bool(any(demand))
        value=W.check_certificate(vectors,demand,c)
        weights=list(map(F,c['weights']));terms=Counter()
        for a,g in zip(weights,groups):
            for f,m in g['positive_face_coefficients']:terms[f]+=a*F(m)
        positive=[[f,str(a)] for f,a in sorted(terms.items()) if a]
        assert sum((alpha.get(f,F())*F(a) for f,a in positive),F())==value
        newmu[v]=max(mu.get(v,F()),value)
        # The bound input proves these are complete maximal count vectors.
        # Verify every saved literal representative's orbit count afresh;
        # physical compatibility and clique-domain completeness are inherited.
        owner={i:g['index'] for g in groups for i in g['literal_stars']}
        assert sum(len(g['literal_stars']) for g in groups)==len(owner)
        for vec,clique in zip(vectors,old['literal_representatives']):
            assert len(clique)==len(set(clique)) and all(i in owner for i in clique)
            assert vec==[Counter(owner[i] for i in clique).get(j,0) for j in range(len(groups))]
        row=M.enc(dict(vertex=v,source=d['source'],group_source=d['group_source'],
            orbit_coefficients=[g['positive_face_coefficients'] for g in groups],
            complete_count_vectors=len(vectors),count_vectors_sha256=jhash(vectors),
            literal_representatives_sha256=jhash(old['literal_representatives']),
            demands=demand,certificate=c,baseline_vertex_mass=mu.get(v,F()),
            exact_local_optimum=value,repaired_vertex_mass=newmu[v],
            required_increase=newmu[v]-mu.get(v,F())))
        if v in cached:assert row==cached[v]
        records.append(row);zero+=not any(demand)
        if value>mu.get(v,F()):
            columns.append(dict(index=len(columns),vertex=v,vertex_cost='1',
                positive_face_coefficients=positive,orbit_weights=c['weights'],
                exact_load=str(value),baseline_vertex_mass=str(mu.get(v,F())),source_geometry=d['source']))
        if saved is None and v not in cached:
            with PROGRESS.open('a') as f:f.write(json.dumps(dict(key=key,row=row),separators=(',',':'))+'\n')
        if step%50==0 or value>mu.get(v,F()):
            print('NEW_ALPHA_LOCAL_COVERS',step,'/571','violations',len(columns),'seconds',round(time.monotonic()-start,2),flush=True)
    assert len(records)==571 and len({r['vertex'] for r in records})==571
    assert all(newmu.get(v,F())>=x for v,x in mu.items())
    model=M.read(M.MODEL);multi=[]
    for col in model['prepared']['new_columns']:
        load=sum((alpha.get(i,F())*m for i,m in col['integer_face_coefficients']),F())
        old_cost=sum((mu.get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        new_cost=sum((newmu.get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        assert load<=old_cost<=new_cost
        if len(col['integer_vertex_capacities'])>1:
            multi.append(dict(variable=col['variable'],all_vertex_costs=col['integer_vertex_capacities'],
                              load=str(load),old_cost=str(old_cost),new_cost=str(new_cost)))
    assert len(multi)==14
    num=F(cert['objective_numerator']);den=sum(newmu.values(),F());ratio=num/den
    for row in records:
        assert F(row['exact_local_optimum'])<=newmu.get(row['vertex'],F())
    return dict(status='PASS_NEW5121_COMPLETE_CURRENT290_EXACT_LOCAL_COVERS',sources=sources,
        model_sha256=model['prepared']['model_sha256'],face_multipliers=cert['face_multipliers'],
        initial_vertex_multipliers=cert['vertex_multipliers'],vertex_multipliers=[[v,str(x)] for v,x in sorted(newmu.items()) if x],
        initial_bound=cert['bound'],objective_numerator=str(num),normalizer=str(den),bound=str(ratio),
        forest_records=records,columns=columns,preserved_multi_owner_columns=multi,
        counts=dict(forests=571,exact_new_demand_local_optima=571,zero_demand_forests=zero,
                    nonzero_demand_forests=571-zero,new_charge_columns=len(columns),raised_vertex_masses=len(columns),
                    mixture_terms=sum(len(r['certificate']['mixture']) for r in records),
                    complete_count_vectors=sum(r['complete_count_vectors'] for r in records)),
        newly_solved_local_LPs=571-zero,global_optimizer_invocations=0,replay_optimizer_invocations=0,
        clique_enumerations=0,old_alpha_optima_reused=0,
        preservation='Face masses are exactly the new5121 dual. Vertex allowances only increase, preserving every inherited full-budget column. All14 added multi-owner columns explicitly retain both costs. Every current290 actual-corner count-vector cone receives a fresh exact charge/cover certificate at the new demands.',
        scope='Complete current290-class/571-forest graph cones only, inheriting frozen necessary graphs and complete clique-count domains including all owned unknowns. New positive historical classes outside290, missing S/simplex and Q coverage remain outside. No wall-pair recensus, global optimization, target12 primal, complete hyperbolic realizability or dimension theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';saved=json.loads(p.read_text()) if args.verify else None
    result=run(saved)
    if args.verify:assert result==saved
    else:
        assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],'bound',float(F(result['bound'])),flush=True)
