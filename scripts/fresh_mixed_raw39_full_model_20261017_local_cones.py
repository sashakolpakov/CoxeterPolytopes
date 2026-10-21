#!/usr/bin/env python3
"""Complete local simplex packing at the new full dual masses.

Only bounded per-forest packing solves; no second global optimization.
All covers and literal marginal features are reconstructed at the new alpha.
"""
import argparse
from collections import defaultdict,Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import time
import fresh_mixed_raw39_full_model_20261017 as E
import fresh_sqrt5_restricted_flag_local_cones_20261017 as R

ROOT=E.ROOT
TAG='fresh_mixed_raw39_full_model_20261017_local_cones'
PROGRESS=Path('/private/tmp')/(TAG+'_progress.jsonl')
A,L=R.A,R.L


def verify_model(d,alpha,mu):
    loads=defaultdict(F)
    for fi,a in alpha.items():
        assert a>0
        for j,m in d['faces'][fi]['coefficients']:loads[j]+=m*a
    for c in d['variables']:
        if c['vertex']is None:assert loads[c['index']]<=0
        else:assert loads[c['index']]<=c['cost']*mu.get(c['vertex'],F())
    return dict(all_columns=len(d['variables']),zero_budget_signed=len(d['zeros']),
        vertex_cost_signed=sum(c['kind']=='forced_quad_incidence_with_vertex_cost'for c in d['variables']))


def verify_mixture(columns,weights,terms):
    byid={c['index']:c for c in columns};cover={j:F()for j in byid};mass=F()
    for term in terms:
        a=F(term['mass']);assert a>0;used=0
        for j,mask in term['blocks']:
            assert j in byid and mask in byid[j]['literal_masks'] and not used&mask
            used|=mask;cover[j]+=a
        mass+=a
    assert all(cover[c['index']]>=w for c,w in zip(columns,weights))
    return mass


def run(saved=None):
    model=E.read(E.TAG);d,prepared=E.build();assert json.loads(E.dump(prepared))==model['prepared']
    cert=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL');E.verify(d,cert)
    alpha={i:F(a)for i,a in cert['face_multipliers']};mu={v:F(a)for v,a in cert['vertex_multipliers']}
    prior=E.read('fresh_mixed_actual_corner_full_model_20261017_local_cones')
    oldalpha={i:F(a)for i,a in prior['face_multipliers']}
    columns,features,dimensions,corners=A.corner_features(d,alpha)
    support=dict(positive_faces=len(alpha),previous_positive_faces=len(oldalpha),
        newly_positive_faces=len(set(alpha)-set(oldalpha)),dropped_faces=len(set(oldalpha)-set(alpha)),
        changed_common_faces=sum(alpha[i]!=oldalpha[i]for i in set(alpha)&set(oldalpha)),
        positive_simplex_faces=len(features),positive_quad_faces=len(alpha)-len(features),
        new_simplex_feature_rows=len({f['face']for f in features}-{f['face']for f in prior['positive_simplex_features']}))
    print('NEW_RAW39_DUAL_SUPPORT',support,flush=True)
    assert len(alpha)==5059 and len(features)==2886

    columns=json.loads(E.dump(columns));features=json.loads(E.dump(features))
    literal=R.literal_orbits(d,columns)
    print('NEW_FULL_DUAL_FEATURES',len(features),len(columns),corners,dimensions,flush=True)
    weights=defaultdict(F);byvertex=defaultdict(list);features_by_v=defaultdict(list)
    colvertex={c['index']:c['vertex']for c in columns}
    for f in features:
        touched=set()
        for j,m in f['block_orbits']:weights[j]+=m*alpha[f['face']];touched.add(colvertex[j])
        for v in touched:features_by_v[v].append(f)
    assert set(weights)==set(colvertex) and all(w>0 for w in weights.values())
    for c in literal:byvertex[c['vertex']].append(c)
    for cs in byvertex.values():cs.sort(key=lambda x:x['index'])
    order=sorted(byvertex,key=lambda v:(-mu.get(v,F()),-len(byvertex[v]),v))
    input_key=sha256(E.dump([prepared['model_sha256'],cert['face_multipliers']]).encode()).hexdigest()
    previous={}
    if saved is None and PROGRESS.exists():
        for line in PROGRESS.read_text().splitlines():
            r=json.loads(line);assert r['key']==input_key and r['vertex']not in previous
            previous[r['vertex']]=r['certificate']
    saved_rows={}if saved is None else{r['vertex']:r for r in saved['local_cones']}
    repaired=dict(mu);records=[];cuts=[];open_opt=[];started=time.monotonic()
    for serial,v in enumerate(order):
        cs=byvertex[v];ws=[weights[c['index']]for c in cs];pack=None;stop=None
        try:pack=L.Packings(cs)
        except L.ResourceStop as exc:stop=str(exc)
        if saved is not None:c=saved_rows[v]['certificate']
        elif v in previous:c=previous[v]
        elif pack is not None:c=L.solve_and_recover(pack,ws,max(mu.get(v,F()),max(ws)))
        else:
            c=dict(status='OPEN_ENUMERATION_WITH_EXACT_SINGLETON_COVER',reason=stop,charges=[],
                exact_lower_bound='0',exact_cover_mass=str(sum(ws,F())),
                packing_mixture=[dict(mass=str(w),blocks=[[x['index'],x['literal_masks'][0]]])for x,w in zip(cs,ws)])
        if pack is not None:
            lower,upper,slacks=L.verify_certificate(pack,ws,c)
            count_vectors=len(pack.vectors)
            vector_digest=sha256(E.dump(pack.vectors).encode()).hexdigest()
        else:
            assert c['status']=='OPEN_ENUMERATION_WITH_EXACT_SINGLETON_COVER' and c['reason']==stop
            lower=F();upper=F(c['exact_cover_mass']);slacks=[];count_vectors=None;vector_digest=None
        assert verify_mixture(cs,ws,c['packing_mixture'])==upper
        repaired[v]=max(mu.get(v,F()),upper)
        if c['status'].startswith('OPEN'):open_opt.append(v)
        if lower>mu.get(v,F()):
            q={j:F(a)for j,a in c['charges']};terms=[]
            for f in features_by_v[v]:
                a=sum((m*q.get(j,F())for j,m in f['block_orbits']),F())
                if a:terms.append([f['face'],str(a)])
            assert sum((alpha[fi]*F(a)for fi,a in terms),F())==lower
            cuts.append(dict(index=len(cuts),vertex=v,vertex_cost=1,positive_simplex_coefficients=terms,
                exact_load=str(lower),old_vertex_mass=str(mu.get(v,F())),orbit_charges=c['charges']))
        rec=dict(vertex=v,forest=d['budgets'][v]['forest'],status=c['status'],certificate=c,
            literal_block_orbits=cs,required_marginals=[[x['index'],str(w)]for x,w in zip(cs,ws)],
            old_vertex_mass=str(mu.get(v,F())),repaired_vertex_mass=str(repaired[v]),
            packing_count_vectors=count_vectors,count_vector_digest=vector_digest,cover_slacks=slacks)
        records.append(rec)
        if saved is None and v not in previous:
            with PROGRESS.open('a')as stream:stream.write(E.dump(dict(key=input_key,vertex=v,certificate=c))+'\n')
        if pack is not None:pack.dp.cache_clear()
        if serial%40==0 or lower>mu.get(v,F()):
            print('NEW_LOCAL_CONES',serial+1,'/',len(order),'cuts',len(cuts),'openopt',len(open_opt),
                'seconds',round(time.monotonic()-started,1),flush=True)
    checked=verify_model(d,alpha,repaired)
    assert all(repaired.get(v,F())>=a for v,a in mu.items())
    num=sum((a*d['faces'][fi]['required']for fi,a in alpha.items()),F());assert num==F(cert['objective_numerator'])
    den=sum(repaired.values(),F());ratio=num/den
    paths=[Path(__file__),Path(E.__file__),ROOT/f'results/{E.TAG}.json',Path(R.__file__),Path(A.__file__),Path(L.__file__),Path(L.H.__file__),
           ROOT/'results/fresh_simplex_disjoint_ports_20261017.md',ROOT/'results/fresh_mixed_actual_corner_full_model_20261017_local_cones.json']
    result=dict(status='PASS_NEW_MIXED_RAW39_DUAL_COMPLETE_LOCAL_PACKING_COVERS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        inherited_model_sha256=prepared['model_sha256'],input_support=support,initial_bound=cert['bound'],bound=str(ratio),
        objective_numerator=str(num),normalizer=str(den),face_multipliers=cert['face_multipliers'],
        initial_vertex_multipliers=cert['vertex_multipliers'],vertex_multipliers=[[v,str(a)]for v,a in sorted(repaired.items())if a],
        positive_simplex_features=features,block_orbit_columns=columns,literal_block_orbits=literal,
        local_cones=records,new_exact_charge_columns=cuts,complete_positive_simplex_forest_coverage=True,
        open_local_optima=open_opt,all_exported_model_columns_checked=checked,global_optimizer_invocations=0,
        replay_optimizer_invocations=0,reused_old_fixed_alpha_covers=0,
        counts=dict(positive_faces=len(alpha),simplex_faces=len(features),quad_faces=len(alpha)-len(features),
            forests=len(records),block_orbits=len(columns),corner_occurrences=corners,dimensions=dimensions,
            exact_local_optima=sum(r['status']=='EXACT_LOCAL_OPTIMUM'for r in records),
            other_exact_covers=sum(r['status']!='EXACT_LOCAL_OPTIMUM'for r in records),
            packing_terms=sum(len(r['certificate']['packing_mixture'])for r in records),
            new_cut_forests=len(cuts),raised_vertex_masses=sum(repaired[v]>mu.get(v,F())for v in repaired)),
        resource_limits=dict(count_vector_entries_per_forest=L.ENTRY_LIMIT,enumeration_seconds=L.ENUM_SECONDS,local_solver_seconds=L.SOLVE_SECONDS),
        theorem='Each saved rational mixture consists of disjoint literal simplex port blocks and covers every actual new-alpha orbit marginal. Simultaneous whole-packing automorphism averaging gives uniform literal marginals, so mixture mass bounds every nonnegative orbit-local simplex charge at once. All finite-model columns, including52 signed vertex-cost inequalities, are rechecked after increasing mu.',
        scope='Complete local simplex-packing coverage at the newly optimized full dual alpha, not reused old covers. Q provenance and joint simplex/prism families are not resolved by this local cone. No second global solve, target12 primal feasibility or dimension theorem.')
    print(result['status'],result['counts'],float(ratio),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=ROOT/f'results/{TAG}.json';saved=json.loads(out.read_text())if a.verify else None
    d=json.loads(E.dump(run(saved)))
    if a.verify:assert d==saved
    else:out.write_text(E.dump(d)+'\n')
