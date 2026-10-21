#!/usr/bin/env python3
"""One changed global solve with the complete raw76 marker-family exports.

All439856 inherited columns remain literal. No physical pair test, ordinary
clique census or local optimization is repeated. The242 structural charges
are valid at arbitrary masses; old fixed-mass repairs are not new optima.
"""
import argparse,json,shutil,time
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import fresh_current327_full_model_20261018 as B

ROOT=B.ROOT;P=B.P;M=B.M;L=B.L
TAG='fresh_raw76_complete_family_full_model_20261018'
FAMILY='fresh_raw76_extra_vertex_complete_family_20261018'
OUT=ROOT/f'results/{TAG}.json'
jhash=B.jhash;digest=B.digest;normal=B.normal


def family_columns():
    doc=M.read(FAMILY);old=M.read(B.TAG)
    assert doc['status']=='PASS_RAW76_COMPLETE_FAMILY_COUPLED_CHARGES'
    for p,h in doc['sources'].items():assert M.digest(p)==h,p
    dual=next(c for c in old['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    rep=doc['exported']['repair'];assert rep['face_multipliers']==dual['face_multipliers']
    alpha={f:F(x)for f,x in dual['face_multipliers']}
    oldmu={v:F(x)for v,x in dual['vertex_multipliers']}
    mu={v:F(x)for v,x in rep['vertex_multipliers']}
    assert all(mu.get(v,F())>=x for v,x in oldmu.items())
    assert F(rep['objective_numerator'])/sum(mu.values(),F())==F(rep['bound'])>12
    domains,represented=L.current_domains();assert len(domains)==600 and len(represented)==324
    vectors={v:r['domain']['maximal_count_vectors']for v,r in domains.items()}
    for tag,_,_ in B.INPUTS[1:]:
        for r in M.read(tag)['forest_records']:vectors[r['vertex']]=r['maximal_count_vectors']
    records={r['vertex']:r for r in doc['forest_records']}
    assert len(records)==len(doc['forest_records'])==235
    assert [(r['face'],r['vertex'])for r in doc['forest_records']]==[(r['face'],r['target_vertex'])for r in doc['family']]
    columns={};inequalities=0
    assert len(doc['exported']['columns'])==242
    for index,col in enumerate(doc['exported']['columns']):
        v=col['vertex'];r=records[v];groups=r['old_groups'];owner={i:g['index']for g in groups for i in g['literal_stars']}
        assert set(owner)==set(range(len(r['old_nodes'])))
        assert all(Counter(o['face']for o in g['origins'])==Counter(dict(g['positive_face_coefficients']))for g in groups)
        w0=list(map(F,col['old_orbit_weights']));w=list(map(F,col['retained_orbit_weights']))
        assert len(w0)==len(w)==len(groups)and all(x>=0 for x in w0)
        vv=vectors[v]if groups else [[]]
        assert len(vv)==col['old_vector_count']and jhash(vv)==col['old_vectors_sha256']
        for a in vv:
            assert len(a)==len(w0)and all(isinstance(n,int)and n>=0 for n in a)
            assert sum((x*n for x,n in zip(w0,a)),F())<=1;inequalities+=1
        omitted=sorted({owner[i]for i,k in r['marker_ordinary_edges']})
        assert omitted==col['omitted_orbits']and w==[F()if i in omitted else x for i,x in enumerate(w0)]
        assert not omitted or(v==1304 and omitted==[0]and groups[0]['positive_face_coefficients']==[[60808,1]])
        types=Counter(r['forest']);cap=min(types['A1']-1,types['H3'])
        assert cap==r['marker_capacity']and F(col['marker_weight'])==F(1,cap)
        assert len(r['markers'])==types['A1']*(types['A1']-1)*types['H3']
        face=Counter()
        for x,g in zip(w,groups):
            for f,n in g['positive_face_coefficients']:face[f]+=x*n
        face[r['face']]+=F(1,cap)
        assert [[f,str(x)]for f,x in sorted(face.items())if x]==col['positive_face_coefficients']
        assert col['source_face']==r['face']and F(col['vertex_cost'])==1
        ff,cc,scale=P.canonical(col);key=ff,cc
        load=sum((alpha.get(f,F())*x for f,x in ff),F())
        before=sum((oldmu.get(v,F())*x for v,x in cc),F())
        fixed=sum((mu.get(v,F())*x for v,x in cc),F())
        factor=F(scale['denominator'],scale['integer_gcd'])
        assert load==F(col['load'])*factor and before==F(col['baseline_cost'])*factor and load<=fixed
        assert load-before==F(col['violation'])*factor
        columns.setdefault(key,[]).append(dict(package=FAMILY,column_index=index,vertex=v,
            record_sha256=jhash(r),column_sha256=jhash(col),scaling=scale,
            old_dual_load=str(load),old_dual_cost=str(before),old_dual_violation=str(load-before)))
    assert len(columns)==242 and sum(len(a)for a in columns.values())==242 and inequalities==1048
    assert sum(F(a[0]['old_dual_violation'])>0 for a in columns.values())==10
    info=dict(family_classes=235,marker_orientations=3204,exported_charges=242,distinct_charges=242,
        old_charge_vector_inequalities=inequalities,old_dual_violations=10,violated_target_types=8,
        preserved_cross_unknowns=6,preserved_marker_unknown_edges=750,
        fixed_mass_exported_only_ratio=rep['bound'],physical_tests_repeated=0,clique_census_repeated=False)
    return columns,alpha,mu,info


def build():
    assert shutil.disk_usage(ROOT).free>200*1024**2
    start=time.monotonic();columns,alpha,mu,info=family_columns()
    print('RAW76_COMPLETE_FAMILY_CHARGES',info,flush=True)
    d,parent=B.build();assert parent==M.read(B.TAG)['prepared']
    nv,nf,nb,nz=map(len,(d['variables'],d['faces'],d['budgets'],d['zeros']))
    assert (nv,nf,nb,nz)==(439856,82806,6558,96176)
    oldvars=digest(d['variables']);oldzeros=digest(d['zeros']);oldactive=set(d['active'])
    oldfaces=[list(map(list,f['coefficients']))for f in d['faces']]
    oldbudgets=[list(map(list,b['coefficients']))for b in d['budgets']]
    oldmeta=[(f['index'],f['kind'],f['required'])for f in d['faces']]
    costs=P.cost_rows(d);faces=[[]for _ in range(nv)]
    for f in d['faces']:
        for j,a in f['coefficients']:faces[j].append((f['index'],a))
    matches={}
    for j,(ff,cc)in enumerate(zip(faces,costs)):
        if not ff or not cc or any(a<=0 for _,a in ff):continue
        assert all(isinstance(a,int)for _,a in ff+cc)
        div=gcd(*(a for _,a in ff+cc));assert div>0
        key=tuple((i,a//div)for i,a in sorted(ff)),tuple((i,a//div)for i,a in sorted(cc))
        if key in columns:matches.setdefault(key,[]).append(j)
    del faces,costs
    added=[];reused=[]
    for(ff,caps),origins in sorted(columns.items()):
        assert len(caps)==1 and all(0<=v<nb and a>0 for v,a in caps)
        assert all(0<=f<nf and a>0 for f,a in ff)
        if(ff,caps)in matches:
            assert all(F(o['old_dual_violation'])<=0 for o in origins)
            reused.append(dict(integer_face_coefficients=[list(x)for x in ff],
                integer_vertex_capacities=[list(x)for x in caps],matching_variables=matches[ff,caps],origins=origins))
            continue
        j=len(d['variables']);v,cap=caps[0]
        d['variables'].append(dict(index=j,kind='complete_raw76_marker_family_export',vertex=v,cost=cap,source_index=len(added)))
        for f,a in ff:d['faces'][f]['coefficients'].append([j,a])
        for v,a in caps:d['budgets'][v]['coefficients'].append([j,a])
        added.append(dict(variable=j,integer_face_coefficients=[list(x)for x in ff],
            integer_vertex_capacities=[list(x)for x in caps],origins=origins))
    assert len(added)+len(reused)==242
    assert sum(any(F(o['old_dual_violation'])>0 for o in c['origins'])for c in added)==10
    assert digest(d['variables'][:nv])==oldvars and digest(d['zeros'])==oldzeros
    assert(len(d['faces']),len(d['budgets']),len(d['zeros']))==(nf,nb,nz)
    for i,f in enumerate(d['faces']):
        assert(f['index'],f['kind'],f['required'])==oldmeta[i]
        assert [[j,a]for j,a in f['coefficients']if j<nv]==oldfaces[i]
        assert all(a>0 for j,a in f['coefficients']if j>=nv)
    for i,b in enumerate(d['budgets']):
        assert [[j,a]for j,a in b['coefficients']if j<nv]==oldbudgets[i]
        assert all(a>0 for j,a in b['coefficients']if j>=nv)
    d['active']=sorted({j for f in d['faces']for j,a in f['coefficients']if a>0})
    assert oldactive<=set(d['active'])and all(j<nv for j in d['repair_indices'])
    multi=parent['inherited_multi_vertex_columns'];assert len(multi)==14 and all(len(c['integer_vertex_capacities'])>1 for c in multi)
    loads=P.loads(d,alpha)
    for j,cc in enumerate(P.cost_rows(d)):
        assert loads[j]<=sum((mu.get(v,F())*a for v,a in cc),F()),j
    h=sha256()
    for group in(d['variables'],d['faces'],d['budgets'],d['zeros']):h.update(digest(group).encode());h.update(b'\n')
    paths={str(Path(mod.__file__).relative_to(ROOT))for mod in(B,B.B,P,P.BASE,P.ENGINE,M,L)}
    paths.add(f'scripts/{TAG}.py')
    for tag in(B.TAG,FAMILY):paths.update(f'{folder}/{tag}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    prepared=dict(sources={p:M.digest(p)for p in sorted(paths)},dimension=13,ground_field='Q(sqrt5)',target=12,
        inherited_model_sha256=parent['model_sha256'],model_sha256=h.hexdigest(),
        provenance_to_face=parent['provenance_to_face'],new_columns=added,reused_columns=reused,
        family_audit=info,inherited_multi_vertex_columns=multi,
        counts=dict(inherited_variables=nv,variables=len(d['variables']),faces=nf,vertex_forests=nb,
            active_variables=len(d['active']),zero_budget_signed_variables=nz,inherited_vertex_cost_signed_variables=52,
            exported_family_columns=242,appended_family_columns=len(added),reused_family_columns=len(reused),
            old_dual_violated_new_columns=10,preserved_multi_owner_columns=14,
            face_matrix_nonzeros=sum(len(f['coefficients'])for f in d['faces'])),
        preservation='All439856 inherited variable records, requirements, signed face coefficients, full budget coefficients, zero records and provenance mappings remain literal;52 positive-cost signed variables and14 multi-owner columns are retained.',
        scope='Finite exported compact arithmetic H13/Qsqrt5 matrix plus all242 complete raw76 family charges. Full family relative to fixed327 ordinary domain, not complete enlarged local cone facets or all actual classes. One global solve and exact full-matrix primal/dual recovery; no physical pair or clique census. Changed masses do not inherit old attaining covers. No chamber or dimension theorem.')
    print('PREPARED_RAW76_COMPLETE_FAMILY_MODEL',prepared['counts'],prepared['model_sha256'],'seconds',round(time.monotonic()-start,2),flush=True)
    return d,normal(prepared)


if __name__=='__main__':
    ap=argparse.ArgumentParser();mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--solve',action='store_true');mode.add_argument('--verify',action='store_true');a=ap.parse_args()
    d,prepared=build()
    if a.verify:
        saved=json.loads(OUT.read_text());assert saved['prepared']==prepared and saved['solver_invocations']==1
        assert saved['certificates']
        for c in saved['certificates']:P.verify(d,c)
        print('PASS_RAW76_COMPLETE_FAMILY_FULL_MATRIX_REPLAY',saved['status'],flush=True)
    else:
        if OUT.exists():
            assert a.solve;saved=json.loads(OUT.read_text());assert saved['prepared']==prepared and saved['solver_invocations']==0
        else:saved=dict(status='PREPARED_RAW76_COMPLETE_FAMILY_MODEL',prepared=prepared,solver_invocations=0,solver=None,certificates=[])
        if a.prepare:OUT.write_text(P.BASE.dump(saved)+'\n')
        else:
            saved.update(status='ONE_CHANGED_RAW76_COMPLETE_FAMILY_SOLVE_STARTED',solver_invocations=1);OUT.write_text(P.BASE.dump(saved)+'\n')
            solver,certificates=P.optimize(d)
            status='PASS_EXACT_RAW76_COMPLETE_FAMILY_MODEL_CERTIFICATES'if certificates else 'UNRESOLVED_SINGLE_RAW76_COMPLETE_FAMILY_SOLVE'
            if any(c.get('meets_target12')for c in certificates):status='PASS_EXACT_RAW76_COMPLETE_FAMILY_PARTIAL_PRIMAL_TARGET12'
            saved.update(status=status,solver=solver,certificates=certificates);OUT.write_text(P.BASE.dump(saved)+'\n')
            print(status,[(c['kind'],c.get('maximum',c.get('bound')))for c in certificates],flush=True)
