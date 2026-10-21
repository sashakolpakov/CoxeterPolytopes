#!/usr/bin/env python3
"""One changed full-model solve with all five frozen marker export packets.

No geometry or clique search is repeated. Exact charge reconstruction and
joint face/budget normalization precede one direct HiGHS invocation; exact
full-matrix recovery and a separate saved replay decide feasibility.
"""
import argparse,gc,json,time
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import fresh_current327_full_model_20261018 as B

ROOT=B.ROOT;P=B.P;M=B.M
TAG='fresh_all_marker_families_full_model_20261018'
OUT=ROOT/f'results/{TAG}.json'
INPUTS=(
 ('fresh_raw76_extra_vertex_complete_family_20261018',242,235),
 ('fresh_raw56_raw57_extra_vertex_family_20261018',1212,1212),
 ('fresh_raw56_raw57_raw76_joint_marker_family_20261018',1430,1394),
 ('fresh_raw108_complete_marker_family_20261018',415,383),
 ('fresh_raw112_complete_marker_family_20261018',529,383))
REPORTS={
 'fresh_raw76_extra_vertex_complete_family_20261018':'9c3faef2ba41b7b3a38e4d42edc744e9388257c279c2036fa8e2320b19163f4f',
 'fresh_raw56_raw57_extra_vertex_family_20261018':'b4dfbe0a4cf0ecb2344f4c400644d9a65f75050edbf343229808cbabe5bd51be',
 'fresh_raw56_raw57_raw76_joint_marker_family_20261018':'2f4170240c76dad4ac739181f25c2736f1824c2726938bd71839d439ee0d383f',
 'fresh_raw108_complete_marker_family_20261018':'1ccf6c8cb0e9e8c28f139972d14f520dec9dcfc6a67921ae57b6936976ca048b',
 'fresh_raw112_complete_marker_family_20261018':'6d55a0961d27d6fcf6cf37cd071b54c4d440c05543be5ef247adfbd8313e7273',
 B.TAG:'69e0e115e01a82f2e5080a0703e76dc36710e6d5203faff49a5f9d3d2a5c36b8'}
jhash=B.jhash;digest=B.digest;normal=B.normal


def saved_vectors():
    """Chronological current600 vectors; no native reconstruction or census."""
    vectors={};cache={};paths=set()
    for tag,_,_ in B.INPUTS:
        doc=M.read(tag);paths.add(f'results/{tag}.json')
        for r in doc['forest_records']:
            v=r['vertex']
            if 'maximal_count_vectors'in r:vv=r['maximal_count_vectors']
            else:
                source=r['geometry_source'];paths.add(f'results/{source}.json')
                if source not in cache:cache[source]={z['vertex']:z for z in M.read(source)['forest_records']}
                vv=cache[source][v]['maximal_count_vectors']
                assert len(vv)==r['complete_count_vectors']and jhash(vv)==r['count_vectors_sha256']
            vectors[v]=vv
    assert len(vectors)==600
    return vectors,paths


def columns():
    vectors,paths=saved_vectors();old=M.read(B.TAG)
    dual=next(c for c in old['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};oldmu={v:F(x)for v,x in dual['vertex_multipliers']}
    docs={s:M.read(s)for s,_,_ in INPUTS};bindings={}
    for s,h in REPORTS.items():assert M.digest(f'results/{s}.md')==h,(s,'report')
    for s,d in docs.items():
        assert d['status'].startswith('PASS')and d['optimizer_invocations']==0
        if s in(INPUTS[1][0],INPUTS[2][0]):assert d['current_graphs_changed']is False
        elif s!=INPUTS[0][0]:assert d['old_outcomes_changed']is False
        for p,h in d['sources'].items():
            if p not in bindings:bindings[p]=M.digest(p)
            assert bindings[p]==h,(s,p)
    rows={s:{r['vertex']:r for r in d['forest_records']}for s,d in docs.items()}
    quotas={s:{q['vertex']:q for q in d['exported'].get('target_allowances',[])}for s,d in docs.items()}
    rep=docs[INPUTS[-1][0]]['exported']['repair'];assert rep['face_multipliers']==dual['face_multipliers']
    mu={v:F(x)for v,x in rep['vertex_multipliers']};assert all(mu.get(v,F())>=x for v,x in oldmu.items())
    assert F(rep['bound'])<12
    result={};counts=[];checks=boundary_checks=0
    for stage,(s,ncols,nrows)in enumerate(INPUTS):
        doc=docs[s];rr=rows[s]
        assert len(rr)==len(doc['forest_records'])==nrows and len(doc['exported']['columns'])==ncols
        assert doc['exported']['repair']['face_multipliers']==dual['face_multipliers']
        for index,c in enumerate(doc['exported']['columns']):
            v=c['vertex'];r=rr[v]
            if stage==2:
                original=rows[INPUTS[1][0]].get(v,rows[INPUTS[0][0]].get(v));assert original is not None
                assert all(jhash(original[k])==r[h]for k,h in(('old_groups','old_groups_sha256'),('old_nodes','old_nodes_sha256'),('old_graph_edges','old_graph_sha256')))
            else:original=r
            groups=original['old_groups'];owner={i:g['index']for g in groups for i in g['literal_stars']}
            assert set(owner)==set(range(len(original['old_nodes'])))
            assert all(Counter(o['face']for o in g['origins'])==Counter(dict(g['positive_face_coefficients']))for g in groups)
            ww=list(map(F,c['old_orbit_weights']));w=list(map(F,c['retained_orbit_weights']))
            assert len(ww)==len(w)==len(groups)and all(z>=0 for z in ww)
            vv=vectors[v]if groups else[[]]
            assert len(vv)==c['old_vector_count']and jhash(vv)==c['old_vectors_sha256']
            for q in vv:
                assert len(q)==len(ww)and all(isinstance(k,int)and k>=0 for k in q)
                assert sum(z*k for z,k in zip(ww,q))<=1;checks+=1
            if stage==0:
                cross=r['marker_ordinary_edges'];faces={76:r['face']};terms=[[76,c['source_face'],c['marker_weight']]]
                order=[76];boundary=[[r['marker_capacity']]]
            elif stage==1:
                cross=[p for z in r['marker_families']for p in z['marker_ordinary_edges']]
                faces={z['raw_family']:z['face']for z in r['marker_families']};terms=c['marker_terms'];order=sorted(faces)
                cap={z['raw_family']:z['marker_capacity']for z in r['marker_families']}
                boundary=[[cap[k]if k==j else 0 for k in order]for j in order]
            elif stage==2:
                cross=[]
                for prior in INPUTS[:2]:
                    z=rows[prior[0]].get(v)
                    if z:cross.extend(z['marker_ordinary_edges']if prior==INPUTS[0]else[p for f in z['marker_families']for p in f['marker_ordinary_edges']])
                faces=dict(r['source_faces']);terms=c['marker_terms'];order=r['envelope']['raw_order'];boundary=r['envelope']['boundary_counts']
            else:
                cross=[p for z in r['marker_families']for p in z['marker_ordinary_edges']];terms=c['marker_terms'];faces={}
                for prior,_,_ in INPUTS[:stage+1]:
                    z=rows[prior].get(v)
                    if not z:continue
                    if prior==INPUTS[0][0]:faces[76]=z['face']
                    elif prior==INPUTS[2][0]:faces.update(z['source_faces'])
                    else:faces.update({f.get('kind',f.get('raw_family')):f['face']for f in z['marker_families']})
                env=quotas[s][v]['envelope'];order=env['kind_order'];boundary=env['boundary_counts']
            omitted={owner[i]for i,k in cross}
            if stage>=3:
                priorquota=quotas[INPUTS[3][0]].get(v)if stage==4 else None
                priorjoint=rows[INPUTS[2][0]].get(v)
                omitted.update(priorquota['omitted_orbits']if priorquota else priorjoint['ordinary_omitted_orbits']if priorjoint else[])
            assert c['omitted_orbits']==sorted(omitted)and w==[F()if i in omitted else z for i,z in enumerate(ww)]
            terms0={k:(f,F(z))for k,f,z in terms};assert len(terms0)==len(terms)
            assert all(k in order and faces[k]==f and z>0 for k,(f,z)in terms0.items())
            weights=[terms0.get(k,(None,F()))[1]for k in order]
            for q in boundary:
                assert len(q)==len(order)and sum(a*b for a,b in zip(q,weights))<=1;boundary_checks+=1
            face=Counter()
            for z,g in zip(w,groups):
                for f,k in g['positive_face_coefficients']:face[f]+=z*k
            for f,z in terms0.values():face[f]+=z
            assert [[f,str(z)]for f,z in sorted(face.items())if z]==c['positive_face_coefficients']and c['vertex_cost']=='1'
            ff,cc,scale=P.canonical(c);key=ff,cc
            load=sum(alpha.get(f,F())*z for f,z in ff);before=sum(oldmu.get(v,F())*z for v,z in cc)
            fixed=sum(mu.get(v,F())*z for v,z in cc);assert load<=fixed
            assert sum(alpha.get(f,F())*F(z)for f,z in c['positive_face_coefficients'])==F(c['load'])
            result.setdefault(key,[]).append(dict(package=s,column_index=index,column_sha256=jhash(c),scaling=scale,
                old_dual_load=str(load),old_dual_cost=str(before),old_dual_violation=str(load-before)))
        counts.append(dict(package=s,exported_columns=ncols,target_forests=nrows))
        paths.update(f'{folder}/{s}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    assert sum(len(z)for z in result.values())==3828
    info=dict(input_counts=counts,export_origins=3828,distinct_joint_primitive_charges=len(result),
        ordinary_vector_inequalities=checks,marker_boundary_inequalities=boundary_checks,
        old_dual_violated_distinct_charges=sum(F(z[0]['old_dual_violation'])>0 for z in result.values()),
        fixed_mass_ratio=rep['bound'],physical_tests_repeated=0,clique_census_repeated=False,
        charge_validity='Frozen complete family geometry and marker envelopes inherited; saved ordinary-vector, omission, marker-boundary and source-coefficient arithmetic reconstructed.')
    return result,alpha,mu,info,paths


def build():
    start=time.monotonic();cols,alpha,mu,info,paths=columns()
    print('COMPLETE_FAMILY_EXPORTS',info,flush=True);gc.collect()
    d,parent=B.build();assert parent==M.read(B.TAG)['prepared']
    nv,nf,nb,nz=map(len,(d['variables'],d['faces'],d['budgets'],d['zeros']))
    assert(nv,nf,nb,nz)==(439856,82806,6558,96176)
    oldvars=digest(d['variables']);oldzeros=digest(d['zeros']);oldactive=set(d['active'])
    oldfaces=[list(map(list,f['coefficients']))for f in d['faces']];oldbudgets=[list(map(list,b['coefficients']))for b in d['budgets']]
    oldmeta=[(f['index'],f['kind'],f['required'])for f in d['faces']]
    costs=P.cost_rows(d);facecols=[[]for _ in range(nv)]
    for f in d['faces']:
        for j,a in f['coefficients']:facecols[j].append((f['index'],a))
    matches={}
    for j,(ff,cc)in enumerate(zip(facecols,costs)):
        if not ff or not cc or any(a<=0 for _,a in ff):continue
        div=gcd(*(a for _,a in ff+cc));key=tuple((i,a//div)for i,a in sorted(ff)),tuple((i,a//div)for i,a in sorted(cc))
        if key in cols:matches.setdefault(key,[]).append(j)
    del costs,facecols
    added=[];reused=[]
    for(ff,cc),origins in sorted(cols.items()):
        assert len(cc)==1 and all(0<=v<nb and a>0 for v,a in cc)and all(0<=f<nf and a>0 for f,a in ff)
        record=dict(integer_face_coefficients=list(map(list,ff)),integer_vertex_capacities=list(map(list,cc)),origins=origins)
        if(ff,cc)in matches:
            assert all(F(o['old_dual_violation'])<=0 for o in origins)
            reused.append(dict(record,matching_variables=matches[ff,cc]));continue
        j=len(d['variables']);v,cap=cc[0]
        d['variables'].append(dict(index=j,kind='complete_marker_families_export',vertex=v,cost=cap,source_index=len(added)))
        for f,a in ff:d['faces'][f]['coefficients'].append([j,a])
        for v,a in cc:d['budgets'][v]['coefficients'].append([j,a])
        added.append(dict(record,variable=j))
    assert len(added)+len(reused)==len(cols)
    assert digest(d['variables'][:nv])==oldvars and digest(d['zeros'])==oldzeros
    assert(len(d['faces']),len(d['budgets']),len(d['zeros']))==(nf,nb,nz)
    for i,f in enumerate(d['faces']):
        assert(f['index'],f['kind'],f['required'])==oldmeta[i]
        assert [[j,a]for j,a in f['coefficients']if j<nv]==oldfaces[i]
        assert all(a>0 for j,a in f['coefficients']if j>=nv)
    for i,b in enumerate(d['budgets']):
        assert [[j,a]for j,a in b['coefficients']if j<nv]==oldbudgets[i]
        assert all(a>0 for j,a in b['coefficients']if j>=nv)
    d['active']=sorted({j for f in d['faces']for j,a in f['coefficients']if a>0});assert oldactive<=set(d['active'])
    assert all(j<nv for j in d['repair_indices'])
    multi=parent['inherited_multi_vertex_columns'];assert len(multi)==14
    loads=P.loads(d,alpha)
    for j,cc in enumerate(P.cost_rows(d)):assert loads[j]<=sum(mu.get(v,F())*a for v,a in cc),j
    modelhash=sha256()
    for group in(d['variables'],d['faces'],d['budgets'],d['zeros']):modelhash.update(digest(group).encode());modelhash.update(b'\n')
    for mod in(B,B.B,P,P.BASE,P.ENGINE,M):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    paths.add(str(Path(__file__).relative_to(ROOT)))
    paths.update(f'results/{B.TAG}.{ext}'for ext in('json','md'))
    prepared=dict(sources={p:M.digest(p)for p in sorted(paths)},dimension=13,ground_field='Q(sqrt5)',target=12,
        inherited_model_sha256=parent['model_sha256'],model_sha256=modelhash.hexdigest(),provenance_to_face=parent['provenance_to_face'],
        family_audit=info,new_columns=added,reused_columns=reused,inherited_multi_vertex_columns=multi,
        counts=dict(inherited_variables=nv,variables=len(d['variables']),faces=nf,vertex_forests=nb,active_variables=len(d['active']),
            zero_budget_signed_variables=nz,inherited_vertex_cost_signed_variables=52,export_origins=3828,
            distinct_family_charges=len(cols),appended_family_columns=len(added),reused_family_columns=len(reused),
            preserved_multi_owner_columns=14,face_matrix_nonzeros=sum(len(f['coefficients'])for f in d['faces'])),
        preservation='Every439856 current327 variable, face coefficient and requirement, budget coefficient and zero record is literal-preserved. All52 positive-budget signed variables and14 multi-owner columns remain.',
        scope='Finite compact arithmetic H13/Qsqrt5 exported matrix with all frozen raw56/57/76/108/112 family exports. No complete enlarged ordinary cone, missing-class or field coverage is inferred. A target12 primal would certify this finite matrix only, not an H13 exclusion.')
    print('PREPARED_ALL_MARKER_FAMILIES',prepared['counts'],prepared['model_sha256'],'seconds',round(time.monotonic()-start,2),flush=True)
    return d,normal(prepared)


def optimize(d):
    import highspy,numpy as np,scipy
    from scipy.sparse import coo_matrix
    active=d['active'];indices={j:i for i,j in enumerate(active)};t=len(active)
    rr=[];cc=[];vv=[];rhs=[]
    for f in d['faces']:
        row=len(rhs);rhs.append(-f['required'])
        for j,a in f['coefficients']:
            if j in indices:rr.append(row);cc.append(indices[j]);vv.append(-a)
    for b in d['budgets']:
        row=len(rhs);rhs.append(0)
        for j,a in b['coefficients']:
            if j in indices:rr.append(row);cc.append(indices[j]);vv.append(a)
        rr.append(row);cc.append(t);vv.append(-1)
    matrix=coo_matrix((np.asarray(vv,dtype=float),(rr,cc)),shape=(len(rhs),t+1)).tocsc()
    lp=highspy.HighsLp();lp.num_col_=t+1;lp.num_row_=len(rhs)
    lp.col_cost_=np.r_[np.zeros(t),1.];lp.col_lower_=np.zeros(t+1);lp.col_upper_=np.full(t+1,highspy.kHighsInf)
    lp.row_lower_=np.full(len(rhs),-highspy.kHighsInf);lp.row_upper_=np.asarray(rhs,dtype=float)
    lp.a_matrix_.format_=highspy.MatrixFormat.kColwise;lp.a_matrix_.start_=matrix.indptr;lp.a_matrix_.index_=matrix.indices;lp.a_matrix_.value_=matrix.data
    h=highspy.Highs()
    for k,v in dict(solver='simplex',simplex_strategy=1,time_limit=900.,threads=4,presolve='on',random_seed=0).items():
        assert h.setOptionValue(k,v)==highspy.HighsStatus.kOk
    assert h.passModel(lp)==highspy.HighsStatus.kOk
    old=next(c for c in M.read(B.TAG)['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    P.verify(d,old) # Literal zero-extension remains feasible in the new model.
    x={j:F(n,old['integer_denominator'])for j,n in old['rounded_variable_numerators']}
    for j in old['repair_variable_indices']:x[j]=x.get(j,F())+F(old['repair_delta'])
    seed=highspy.HighsSolution();seed.col_value=[float(x.get(j,F()))for j in active]+[float(F(old['maximum']))];seed.value_valid=True
    assert h.setSolution(seed)==highspy.HighsStatus.kOk
    print('ONE_CHANGED_DIRECT_HIGHS_BEGIN',matrix.shape,matrix.nnz,'warm_primal',old['maximum'],'limit900',flush=True)
    start=time.monotonic();runstatus=h.run();elapsed=time.monotonic()-start
    sol=h.getSolution();info=h.getInfo();status=h.getModelStatus()
    solver=dict(interface='highspy',highs_version=h.version(),scipy_version=scipy.__version__,numpy_version=np.__version__,
        method='dual simplex',simplex_strategy=1,presolve='on',threads=4,time_limit_seconds=900,
        run_status=str(runstatus),model_status=str(status),model_status_string=h.modelStatusToString(status),
        success=status==highspy.HighsModelStatus.kOptimal,elapsed_seconds=elapsed,iterations=info.simplex_iteration_count,
        objective=info.objective_function_value,primal_solution_status=info.primal_solution_status,dual_solution_status=info.dual_solution_status,
        value_valid=sol.value_valid,dual_valid=sol.dual_valid,warm_start='Inherited exact primal with every new variable zero, converted to floats.')
    rawpath=ROOT/f'results/{TAG}_numerical.json'
    assert not rawpath.exists();rawpath.write_text(json.dumps(dict(solver=solver,col_value=list(sol.col_value),row_dual=list(sol.row_dual)))+'\n')
    print('ONE_CHANGED_DIRECT_HIGHS_END',solver,flush=True)
    certificates=[]
    if sol.value_valid:certificates.append(P.ENGINE.recover_primal(d,sol.col_value))
    if sol.dual_valid:
        try:certificates.append(P.recover_dual(d,sol.row_dual))
        except AssertionError as error:
            solver['dual_exact_recovery_failure']=str(error)
    if not any(c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL'for c in certificates):certificates.append(old)
    for c in certificates:P.verify(d,c)
    return solver,certificates,str(rawpath.relative_to(ROOT)),M.digest(str(rawpath.relative_to(ROOT)))


if __name__=='__main__':
    ap=argparse.ArgumentParser();mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--solve',action='store_true');mode.add_argument('--verify',action='store_true');args=ap.parse_args()
    d,prepared=build()
    if args.verify:
        saved=json.loads(OUT.read_text());assert saved['prepared']==prepared and saved['solver_invocations']==1
        assert saved['certificates']and M.digest(saved['numerical_snapshot'])==saved['numerical_snapshot_sha256']
        for c in saved['certificates']:P.verify(d,c)
        print('PASS_ALL_MARKER_FAMILIES_FULL_MATRIX_REPLAY',saved['status'],flush=True)
    else:
        if OUT.exists():
            assert args.solve;saved=json.loads(OUT.read_text());assert saved['prepared']==prepared and saved['solver_invocations']==0
        else:saved=dict(status='PREPARED_ALL_MARKER_FAMILIES_MODEL',prepared=prepared,solver_invocations=0,solver=None,certificates=[])
        if args.prepare:OUT.write_text(P.BASE.dump(saved)+'\n')
        else:
            saved.update(status='ONE_CHANGED_ALL_MARKER_FAMILIES_SOLVE_STARTED',solver_invocations=1);OUT.write_text(P.BASE.dump(saved)+'\n')
            solver,certificates,path,h=optimize(d)
            kinds={c['kind']for c in certificates}
            status='PASS_EXACT_ALL_MARKER_FAMILIES_INTERVAL'if len(kinds)==2 else'PARTIAL_ALL_MARKER_FAMILIES_CERTIFICATE'
            if any(c.get('meets_target12')for c in certificates):status='PASS_EXACT_ALL_MARKER_FAMILIES_PARTIAL_PRIMAL_TARGET12'
            saved.update(status=status,solver=solver,certificates=certificates,numerical_snapshot=path,numerical_snapshot_sha256=h)
            OUT.write_text(P.BASE.dump(saved)+'\n');print(status,[(c['kind'],c.get('maximum',c.get('bound')))for c in certificates],flush=True)
