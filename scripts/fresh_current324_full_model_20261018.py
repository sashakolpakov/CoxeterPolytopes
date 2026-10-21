#!/usr/bin/env python3
"""One changed full-model solve after the current324 local charge certificates.

Every inherited column/cost is preserved. Saved local charges are universally
valid on their bound conservative graph; their old attaining mixtures are not
claimed optimal or covering at newly optimized face masses.
"""
import argparse,json,shutil,time
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import fresh_cross_vertex_current_full_model_20261018 as P
import fresh_cross_vertex_support_migration_20261018 as M
import fresh_current_historical_weighted_corners_20261017 as H

ROOT=P.ROOT;TAG='fresh_current324_full_model_20261018';OUT=ROOT/f'results/{TAG}.json'
INPUTS=(('fresh_cross_vertex_current290_local_covers_20261018',571,16),
        ('fresh_cross_vertex_historical32_local_covers_20261018',161,27),
        ('fresh_new5121_raw4_two_face_local_covers_20261018',6,4))


def digest(x):return P.BASE.digest(x)
def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def normal(x):return json.loads(P.BASE.dump(x))


def local_columns():
    """Read all738 saved charge certificates; no graph census or local solve."""
    domains=M.current_domains();assert len(domains)==571
    final=M.read(INPUTS[-1][0]);alpha={f:F(x) for f,x in final['face_multipliers']};mu={v:F(x) for v,x in final['vertex_multipliers']}
    olddual=next(c for c in M.read(P.TAG)['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert olddual['face_multipliers']==final['face_multipliers']
    oldmu={v:F(x) for v,x in olddual['vertex_multipliers']}
    assert all(mu.get(v,F())>=x for v,x in oldmu.items())
    docs={};sources={};columns={};counts=[];export_keys=set();certificates=inequalities=zero=0
    for tag,expected,expected_exports in INPUTS:
        doc=M.read(tag);docs[tag]=doc
        assert doc['status'].startswith('PASS') and len(doc['forest_records'])==expected and len(doc['columns'])==expected_exports
        assert doc['face_multipliers']==final['face_multipliers']
        for path,h in doc['sources'].items():assert M.digest(path)==h,(tag,path)
        npositive=0
        for index,row in enumerate(doc['forest_records']):
            v=row['vertex'];c=row['certificate']
            if tag==INPUTS[0][0]:
                dd=domains[v];vectors=dd['row']['maximal_count_vectors'];groups=dd['groups'];reps=dd['row']['literal_representatives']
                assert row['count_vectors_sha256']==jhash(vectors) and row['literal_representatives_sha256']==jhash(reps)
                assert row['orbit_coefficients']==[g['positive_face_coefficients'] for g in groups]
            else:
                vectors=row['maximal_count_vectors'];groups=row['domain']['corner_orbits'];reps=row['literal_representatives']
            weights=list(map(F,c['weights']));assert len(weights)==len(groups) and all(x>=0 for x in weights)
            demand=[sum((alpha.get(f,F())*F(m) for f,m in g['positive_face_coefficients']),F()) for g in groups]
            assert list(map(str,demand))==(row['demands'] if tag==INPUTS[0][0] else row['domain']['demands'])
            value=H.W.check_certificate(vectors,demand,c)
            assert value==F(row['exact_local_optimum'])<=mu.get(v,F())
            owner={s:g['index'] for g in groups for s in g['literal_stars']}
            assert len(owner)==sum(len(g['literal_stars']) for g in groups) and len(vectors)==len(reps)
            for vector,cl in zip(vectors,reps):
                assert len(cl)==len(set(cl)) and vector==[Counter(owner[i] for i in cl).get(j,0) for j in range(len(groups))]
            terms=Counter()
            for w,g in zip(weights,groups):
                for f,m in g['positive_face_coefficients']:terms[f]+=w*F(m)
            faces=[[f,str(x)] for f,x in sorted(terms.items()) if x]
            certificates+=1;inequalities+=len(vectors)
            if not faces:zero+=1;continue
            assert sum((alpha.get(f,F())*F(x) for f,x in faces),F())==value
            raw=dict(vertex=v,vertex_cost='1',positive_face_coefficients=faces)
            ff,caps,scale=P.canonical(raw);key=(ff,caps)
            columns.setdefault(key,[]).append(dict(package=tag,record_index=index,vertex=v,
                record_sha256=jhash(row),scaling=scale))
            npositive+=1
        for index,col in enumerate(doc['columns']):
            ff,caps,scale=P.canonical(col);key=(ff,caps);assert key in columns
            assert any(x['package']==tag and x['vertex']==col['vertex'] for x in columns[key])
            # All47 exported violated old-alpha columns genuinely strengthen
            # the previous optimized matrix, not just its209 later additions.
            load=sum((alpha.get(f,F())*a for f,a in ff),F());oldcost=sum((oldmu.get(v,F())*a for v,a in caps),F())
            assert load>oldcost,(tag,index)
            export_keys.add(key)
        counts.append(dict(package=tag,records=expected,nonzero_charges=npositive,exports=expected_exports))
        for folder,ext in (('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=M.digest(path)
    assert certificates==738 and len(export_keys)==47
    for tag in M.DOMAIN_STEMS:
        path=f'results/{tag}.json';sources[path]=M.digest(path)
    return columns,export_keys,dict(sources=sources,input_counts=counts,saved_charge_certificates=certificates,
        saved_charge_inequalities=inequalities,zero_charge_records=zero,distinct_local_charges=len(columns),
        old_dual_violated_exports=47,baseline_fixed_mass_bound=final['bound'],alpha=alpha,mu=mu)


def build():
    assert shutil.disk_usage(ROOT).free>200*1024**2
    start=time.monotonic();columns,exports,info=local_columns()
    print('CURRENT324_CHARGE_AUDIT', {k:v for k,v in info.items() if k not in ('sources','alpha','mu')},flush=True)
    d,parent=P.build();assert parent==M.read(P.TAG)['prepared']
    nv,nf,nb,nz=map(len,(d['variables'],d['faces'],d['budgets'],d['zeros']))
    assert (nv,nf,nb,nz)==(439619,82806,6558,96176)
    old_variables=digest(d['variables']);old_zeros=digest(d['zeros'])
    oldfaces=[list(map(list,f['coefficients'])) for f in d['faces']]
    oldbudgets=[list(map(list,b['coefficients'])) for b in d['budgets']]
    oldmeta=[(f['index'],f['kind'],f['required']) for f in d['faces']];oldactive=set(d['active'])
    costs=P.cost_rows(d);terms=[[] for _ in range(nv)]
    for f in d['faces']:
        for j,a in f['coefficients']:terms[j].append((f['index'],a))
    inherited={}
    for j,(ff,cc) in enumerate(zip(terms,costs)):
        if not ff or not cc or any(a<=0 for f,a in ff):continue
        assert all(isinstance(a,int) for i,a in ff+cc)
        divisor=gcd(*(a for i,a in ff+cc));assert divisor>0
        key=(tuple((i,a//divisor) for i,a in sorted(ff)),tuple((i,a//divisor) for i,a in sorted(cc)))
        if key in columns:inherited.setdefault(key,[]).append(j)
    assert not exports&set(inherited)
    del terms,costs
    added=[];reused=[]
    for (ff,caps),origins in sorted(columns.items()):
        assert len(caps)==1 and caps[0][1]>0
        load=sum((info['alpha'].get(f,F())*a for f,a in ff),F());cost=sum((info['mu'].get(v,F())*a for v,a in caps),F())
        assert load<=cost
        record=dict(integer_face_coefficients=[list(x) for x in ff],integer_vertex_capacities=[list(x) for x in caps],
            origins=origins,violates_previous_optimized_dual=(ff,caps) in exports)
        if (ff,caps) in inherited:
            reused.append(dict(**record,inherited_variables=inherited[(ff,caps)]));continue
        j=len(d['variables']);v,cap=caps[0]
        d['variables'].append(dict(index=j,kind='current324_saved_local_charge',vertex=v,cost=cap,source_index=len(added)))
        for f,a in ff:d['faces'][f]['coefficients'].append([j,a])
        for v,a in caps:d['budgets'][v]['coefficients'].append([j,a])
        added.append(dict(**record,variable=j))
    assert digest(d['variables'][:nv])==old_variables and digest(d['zeros'])==old_zeros
    assert (len(d['faces']),len(d['budgets']),len(d['zeros']))==(nf,nb,nz)
    for i,f in enumerate(d['faces']):
        assert (f['index'],f['kind'],f['required'])==oldmeta[i]
        assert [[j,a] for j,a in f['coefficients'] if j<nv]==oldfaces[i]
        assert all(a>0 for j,a in f['coefficients'] if j>=nv)
    for i,b in enumerate(d['budgets']):
        assert [[j,a] for j,a in b['coefficients'] if j<nv]==oldbudgets[i]
        assert all(a>0 for j,a in b['coefficients'] if j>=nv)
    assert sum(x['violates_previous_optimized_dual'] for x in added)==47
    d['active']=sorted({j for f in d['faces'] for j,a in f['coefficients'] if a>0});assert oldactive<=set(d['active'])
    assert all(j<nv for j in d['repair_indices'])
    assert sum(len(c['integer_vertex_capacities'])>1 for c in parent['new_columns'])==14
    modelhash=sha256()
    for group in (d['variables'],d['faces'],d['budgets'],d['zeros']):modelhash.update(digest(group).encode());modelhash.update(b'\n')
    paths={str(Path(mod.__file__).relative_to(ROOT)) for mod in (P,P.BASE,P.ENGINE,M,H,H.W)}
    paths.update((f'scripts/{TAG}.py',f'results/{P.TAG}.json',f'results/{P.TAG}.md'))
    sources={**info.pop('sources'),**{path:M.digest(path) for path in paths}}
    del info['alpha'],info['mu']
    prepared=dict(sources=sources,dimension=13,ground_field='Q(sqrt5)',target=12,
        inherited_model_sha256=parent['model_sha256'],model_sha256=modelhash.hexdigest(),
        provenance_to_face=parent['provenance_to_face'],new_columns=added,reused_local_columns=reused,local_audit=info,
        inherited_several_vertex_columns=parent['new_columns'],
        counts=dict(inherited_variables=nv,variables=len(d['variables']),faces=nf,vertex_forests=nb,
            active_variables=len(d['active']),zero_budget_signed_variables=nz,inherited_vertex_cost_signed_variables=52,
            appended_local_columns=len(added),reused_local_columns=len(reused),preserved_multi_owner_columns=14,
            old_dual_violated_new_columns=47,face_matrix_nonzeros=sum(len(f['coefficients']) for f in d['faces'])),
        preservation='Every439619 old variable, face requirement/coefficient, budget coefficient and zero record is literal-preserved, including52 positive-cost signed variables and14 two-owner columns. All nonzero saved local charges from738 records are appended or matched to identical positive-scaled inherited columns; no attaining old-mass cover is assigned to a new mass.',
        scope='Changed finite exported compact-arithmetic H13/Qsqrt5 matrix. One bounded solve plus exact full-matrix recovery/replay. Complete324 graph-cone charge families are not all exported, and missing face/field/global geometric coverage remains. No dimension theorem or realization follows from this finite solve.')
    print('PREPARED_CURRENT324_MODEL',prepared['counts'],prepared['model_sha256'],'seconds',round(time.monotonic()-start,2),flush=True)
    return d,normal(prepared)


if __name__=='__main__':
    ap=argparse.ArgumentParser();mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--solve',action='store_true');mode.add_argument('--verify',action='store_true')
    args=ap.parse_args();d,prepared=build()
    if args.prepare:
        assert not OUT.exists();OUT.write_text(P.BASE.dump(dict(status='PREPARED_CURRENT324_FULL_MODEL',prepared=prepared,solver_invocations=0,solver=None,certificates=[]))+'\n')
    elif args.verify:
        saved=json.loads(OUT.read_text());assert saved['prepared']==prepared
        for certificate in saved['certificates']:P.verify(d,certificate)
        print('PASS_CURRENT324_FULL_MATRIX_REPLAY',saved['status'],flush=True)
    else:
        saved=json.loads(OUT.read_text());assert saved['prepared']==prepared and saved['solver_invocations']==0
        saved.update(status='ONE_CHANGED_CURRENT324_SOLVE_STARTED',solver_invocations=1);OUT.write_text(P.BASE.dump(saved)+'\n')
        solver,certificates=P.optimize(d)
        status='PASS_EXACT_CURRENT324_MODEL_CERTIFICATES' if certificates else 'UNRESOLVED_SINGLE_CURRENT324_SOLVE'
        if any(c.get('meets_target12') for c in certificates):status='PASS_EXACT_CURRENT324_PARTIAL_PRIMAL_TARGET12'
        saved.update(status=status,solver=solver,certificates=certificates);OUT.write_text(P.BASE.dump(saved)+'\n')
        print(status,[(c['kind'],c.get('maximum',c.get('bound'))) for c in certificates],flush=True)
