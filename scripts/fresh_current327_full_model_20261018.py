#!/usr/bin/env python3
"""One changed global solve with34 actual new5073/raw76 charge exports.

All439822 inherited columns, including14 multi-owner charges, stay literal.
No native pair/clique census or local optimization is repeated. Old attaining
covers certify only their original masses; charge inequalities are universal
on their bound conservative graph domains.
"""
import argparse,json,shutil,time
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import fresh_current324_full_model_20261018 as B
import fresh_current324_new5073_local_covers_20261018 as L

P=B.P;M=B.M;H=B.H;ROOT=B.ROOT
TAG='fresh_current327_full_model_20261018';OUT=ROOT/f'results/{TAG}.json'
INPUTS=(('fresh_current324_new5073_local_covers_20261018',600,25),
        ('fresh_new5073_raw76_high_mass_local_covers_20261018',6,6),
        ('fresh_new5073_raw76_third_class_local_covers_20261018',3,3))


def jhash(x):return B.jhash(x)
def digest(x):return B.digest(x)
def normal(x):return B.normal(x)


def charge_columns():
    domains,represented=L.current_domains();assert len(domains)==600 and len(represented)==324
    final=M.read(INPUTS[-1][0]);assert final['status']=='PASS_NEW5073_RAW76_THIRD_CLASS_EXACT_LOCAL_COVERS'
    alpha={f:F(x) for f,x in final['face_multipliers']};mu={v:F(x) for v,x in final['vertex_multipliers']}
    previous=M.read(B.TAG)
    olddual=next(c for c in previous['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert olddual['face_multipliers']==final['face_multipliers']
    oldmu={v:F(x) for v,x in olddual['vertex_multipliers']}
    assert all(mu.get(v,F())>=x for v,x in oldmu.items()) and F(final['bound'])<12
    assert final['represented_classes']==327 and final['represented_forests']==600
    columns={};sources={};input_counts=[];inequalities=terms=0
    for tag,nrows,nexports in INPUTS:
        doc=M.read(tag)
        assert doc['status'].startswith('PASS') and doc['face_multipliers']==final['face_multipliers']
        rows={r['vertex']:r for r in doc['forest_records']}
        assert len(rows)==len(doc['forest_records'])==nrows and len(doc['columns'])==nexports
        for path,h in doc['sources'].items():assert M.digest(path)==h,(tag,path)
        seen=set()
        for index,col in enumerate(doc['columns']):
            v=col['vertex'];assert v not in seen;seen.add(v);r=rows[v];c=r['certificate']
            if tag==INPUTS[0][0]:
                dd=domains[v];vectors=dd['domain']['maximal_count_vectors'];coeffs=r['orbit_coefficients'];demands=r['demands']
                assert r['source_record_sha256']==jhash(dd['record']) and r['source']==dd['source']
                assert coeffs==[g['positive_face_coefficients'] for g in dd['groups']]
                assert len(vectors)==r['complete_count_vectors']
                assert jhash(vectors)==r['count_vectors_sha256']
            else:
                vectors=r['maximal_count_vectors'];coeffs=[g['positive_face_coefficients'] for g in r['domain']['corner_orbits']]
                demands=r['domain']['demands']
            weights=list(map(F,c['weights']));assert len(weights)==len(coeffs) and all(w>=0 for w in weights)
            dem=[sum((alpha.get(f,F())*F(a) for f,a in co),F()) for co in coeffs]
            assert list(map(str,dem))==demands
            value=H.W.check_certificate(vectors,dem,c)
            assert value==F(r['exact_local_optimum'])<=mu.get(v,F())
            inequalities+=len(vectors);terms+=len(c['mixture']);face=Counter()
            for w,co in zip(weights,coeffs):
                for f,a in co:face[f]+=w*F(a)
            expected=[[f,str(a)] for f,a in sorted(face.items()) if a]
            assert expected==col['positive_face_coefficients'] and F(col['vertex_cost'])==1
            ff,caps,scale=P.canonical(col);key=ff,caps
            load=sum((alpha.get(f,F())*a for f,a in ff),F())
            oldcost=sum((oldmu.get(v,F())*a for v,a in caps),F())
            fixedcost=sum((mu.get(v,F())*a for v,a in caps),F())
            assert load>oldcost and load<=fixedcost
            columns.setdefault(key,[]).append(dict(package=tag,column_index=index,vertex=v,
                record_sha256=jhash(r),scaling=scale,old_dual_load=str(load),old_dual_cost=str(oldcost)))
        input_counts.append(dict(package=tag,available_records=nrows,exported_charges=nexports))
        for folder,ext in (('scripts','py'),('results','json'),('results','md')):
            p=f'{folder}/{tag}.{ext}';sources[p]=M.digest(p)
    assert sum(len(v) for v in columns.values())==34 and len(columns)==34
    return columns,alpha,mu,dict(sources=sources,input_counts=input_counts,
        exported_charge_certificates=34,charge_inequalities=inequalities,mixture_terms=terms,
        distinct_charges=34,old_dual_violations=34,fixed_mass_ratio=final['bound'])


def build():
    assert shutil.disk_usage(ROOT).free>200*1024**2
    start=time.monotonic();columns,alpha,mu,info=charge_columns()
    print('CURRENT327_NEW_CHARGES',{k:v for k,v in info.items() if k!='sources'},flush=True)
    d,parent=B.build();assert parent==M.read(B.TAG)['prepared']
    nv,nf,nb,nz=map(len,(d['variables'],d['faces'],d['budgets'],d['zeros']))
    assert (nv,nf,nb,nz)==(439822,82806,6558,96176)
    oldvars=digest(d['variables']);oldzeros=digest(d['zeros'])
    oldfaces=[list(map(list,f['coefficients'])) for f in d['faces']]
    oldbudgets=[list(map(list,b['coefficients'])) for b in d['budgets']]
    oldmeta=[(f['index'],f['kind'],f['required']) for f in d['faces']];oldactive=set(d['active'])
    costs=P.cost_rows(d);faces=[[] for _ in range(nv)]
    for f in d['faces']:
        for j,a in f['coefficients']:faces[j].append((f['index'],a))
    matches={}
    for j,(ff,cc) in enumerate(zip(faces,costs)):
        if not ff or not cc or any(a<=0 for _,a in ff):continue
        assert all(isinstance(a,int) for _,a in ff+cc)
        div=gcd(*(a for _,a in ff+cc));assert div>0
        key=(tuple((i,a//div) for i,a in sorted(ff)),tuple((i,a//div) for i,a in sorted(cc)))
        if key in columns:matches.setdefault(key,[]).append(j)
    assert not matches,'Every exported charge strictly violated the saved exact old dual.'
    del faces,costs
    added=[]
    for (ff,caps),origins in sorted(columns.items()):
        assert len(caps)==1 and caps[0][1]>0
        assert all(0<=f<nf and a>0 for f,a in ff) and all(0<=v<nb and a>0 for v,a in caps)
        j=len(d['variables']);v,cap=caps[0]
        d['variables'].append(dict(index=j,kind='current327_saved_local_export',vertex=v,cost=cap,source_index=len(added)))
        for f,a in ff:d['faces'][f]['coefficients'].append([j,a])
        for v,a in caps:d['budgets'][v]['coefficients'].append([j,a])
        added.append(dict(variable=j,integer_face_coefficients=[list(x) for x in ff],
            integer_vertex_capacities=[list(x) for x in caps],origins=origins))
    assert len(added)==34 and digest(d['variables'][:nv])==oldvars and digest(d['zeros'])==oldzeros
    assert (len(d['faces']),len(d['budgets']),len(d['zeros']))==(nf,nb,nz)
    for i,f in enumerate(d['faces']):
        assert (f['index'],f['kind'],f['required'])==oldmeta[i]
        assert [[j,a] for j,a in f['coefficients'] if j<nv]==oldfaces[i]
        assert all(a>0 for j,a in f['coefficients'] if j>=nv)
    for i,b in enumerate(d['budgets']):
        assert [[j,a] for j,a in b['coefficients'] if j<nv]==oldbudgets[i]
        assert all(a>0 for j,a in b['coefficients'] if j>=nv)
    d['active']=sorted({j for f in d['faces'] for j,a in f['coefficients'] if a>0});assert oldactive<=set(d['active'])
    assert all(j<nv for j in d['repair_indices'])
    multi=parent['inherited_several_vertex_columns'];assert len(multi)==209
    assert sum(len(c['integer_vertex_capacities'])>1 for c in multi)==14
    # Check the final fixed-mass repair against the entire new matrix, before
    # proposing a new global optimum. Signed face loads and zero costs included.
    load=P.loads(d,alpha)
    for j,cost in enumerate(P.cost_rows(d)):
        assert load[j]<=sum((mu.get(v,F())*a for v,a in cost),F()),j
    modelhash=sha256()
    for group in (d['variables'],d['faces'],d['budgets'],d['zeros']):modelhash.update(digest(group).encode());modelhash.update(b'\n')
    paths={str(Path(mod.__file__).relative_to(ROOT)) for mod in (B,P,P.BASE,P.ENGINE,M,H,H.W)}
    paths.update((f'scripts/{TAG}.py',f'results/{B.TAG}.json',f'results/{B.TAG}.md'))
    sources={**info.pop('sources'),**{p:M.digest(p) for p in paths}}
    prepared=dict(sources=sources,dimension=13,ground_field='Q(sqrt5)',target=12,
        inherited_model_sha256=parent['model_sha256'],model_sha256=modelhash.hexdigest(),
        provenance_to_face=parent['provenance_to_face'],new_columns=added,local_audit=info,
        inherited_multi_vertex_columns=[c for c in multi if len(c['integer_vertex_capacities'])>1],
        counts=dict(inherited_variables=nv,variables=len(d['variables']),faces=nf,vertex_forests=nb,
            active_variables=len(d['active']),zero_budget_signed_variables=nz,inherited_vertex_cost_signed_variables=52,
            appended_local_columns=len(added),reused_local_columns=0,preserved_multi_owner_columns=14,
            old_dual_violated_new_columns=34,face_matrix_nonzeros=sum(len(f['coefficients']) for f in d['faces'])),
        preservation='All439822 inherited variables, face requirements and coefficients, budget coefficients, zero records and provenance mappings are literal-preserved, including52 positive-cost signed columns and14 two-owner costs. Every34 exported charge is independently primitive-normalized and absent from the entire old matrix.',
        scope='Changed finite exported compact arithmetic H13/Qsqrt5 matrix. One global solve with exact full-matrix primal/dual recovery and separate replay. Local charge families for327 represented classes are not exhausted; old attaining covers are not claimed valid at new masses. No target feasibility outside this finite matrix, chamber or dimension theorem follows.')
    print('PREPARED_CURRENT327_MODEL',prepared['counts'],prepared['model_sha256'],'seconds',round(time.monotonic()-start,2),flush=True)
    return d,normal(prepared)


if __name__=='__main__':
    ap=argparse.ArgumentParser();mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--solve',action='store_true');mode.add_argument('--verify',action='store_true');a=ap.parse_args()
    d,prepared=build()
    if a.verify:
        saved=json.loads(OUT.read_text());assert saved['prepared']==prepared and saved['solver_invocations']==1
        for c in saved['certificates']:P.verify(d,c)
        assert saved['certificates']
        print('PASS_CURRENT327_FULL_MATRIX_REPLAY',saved['status'],flush=True)
    else:
        if OUT.exists():
            assert a.solve;saved=json.loads(OUT.read_text());assert saved['prepared']==prepared and saved['solver_invocations']==0
        else:saved=dict(status='PREPARED_CURRENT327_FULL_MODEL',prepared=prepared,solver_invocations=0,solver=None,certificates=[])
        if a.prepare:OUT.write_text(P.BASE.dump(saved)+'\n')
        else:
            saved.update(status='ONE_CHANGED_CURRENT327_SOLVE_STARTED',solver_invocations=1);OUT.write_text(P.BASE.dump(saved)+'\n')
            solver,certificates=P.optimize(d)
            status='PASS_EXACT_CURRENT327_MODEL_CERTIFICATES' if certificates else 'UNRESOLVED_SINGLE_CURRENT327_SOLVE'
            if any(c.get('meets_target12') for c in certificates):status='PASS_EXACT_CURRENT327_PARTIAL_PRIMAL_TARGET12'
            saved.update(status=status,solver=solver,certificates=certificates);OUT.write_text(P.BASE.dump(saved)+'\n')
            print(status,[(c['kind'],c.get('maximum',c.get('bound'))) for c in certificates],flush=True)
