#!/usr/bin/env python3
"""Install sixteen frozen coface transfers; test an exact cheap lower bound.

All face masses are variables in this changed finite matrix. This packet tests
one feasible assignment before spending any numerical optimizer budget.
The31 absent simplex rows remain uninstalled. No geometry is repeated.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_all_marker_families_full_model_20261018 as A

ROOT=A.ROOT;P=A.P;M=A.M
TAG='fresh_raw357_coface_full_model_20261018'
COFACE='fresh_raw357_simplex_coface_transfer_20261018'
MODEL='fresh_all_marker_families_dual_continuation_20261018'
PINS={
 f'scripts/{COFACE}.py':'e68e8e31cd60b21c41ea546a88f4bbd82b59a1330daa7cf65d620c68b879358f',
 f'results/{COFACE}.json':'82e49ed6794e759614dfd3de30829f67a93bc083e12194600b8802f2dbf6d94b',
 f'results/{COFACE}.md':'b62128fefc397cc3faf465dc23ff0a6503ede5a786a534313b032f7f193090dd',
 'scripts/check_raw357_simplex_coface_saved_20261018.py':'83c28a342b19759cbc24651ccb3d04974cb1fe19bc2c1549ec4ff125d356008f'}


def build():
    for p,h in PINS.items():assert M.digest(p)==h,p
    doc=M.read(COFACE);old=M.read(MODEL)
    assert doc['status']=='PASS_RAW357_COMPLETE_TYPED_SIMPLEX_COFACE_TRANSFER'
    for source in (doc,old):
        for p,h in source['sources'].items():assert M.digest(p)==h,p
    d,parent=A.build();assert parent==old['prepared']
    nv,nf,nb,nz=map(len,(d['variables'],d['faces'],d['budgets'],d['zeros']))
    assert(nv,nf,nb,nz)==(442210,82806,6558,96176)
    hashes={k:A.digest(d[k])for k in('variables','faces','budgets','zeros')}
    oldfaces=[list(map(list,f['coefficients']))for f in d['faces']]
    oldmeta=[(f['index'],f['kind'],f['required'])for f in d['faces']]
    oldactive=set(d['active']);sourcefaces={r['source_face']for r in doc['family']}
    assert len(sourcefaces)==47 and len(doc['family'])==47
    olddual=next(c for c in old['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    oldprimal=next(c for c in old['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    for c in(olddual,oldprimal):P.verify(d,c)
    alpha={f:F(x)for f,x in olddual['face_multipliers']}
    deleted=[[f,str(alpha[f])]for f in sorted(sourcefaces)if f in alpha]
    assert len(deleted)==9 and sum(F(x)for _,x in deleted)==F(83284683,10**10)
    # Inspect, rather than assume, the sign of every old coefficient at a
    # changed face. Full feasibility below also checks all unrelated rows.
    negative_terms=[[f,j,a]for f,_ in deleted for j,a in d['faces'][f]['coefficients']if a<0]
    records=[]
    for k,c in enumerate(doc['existing_target_transfer_columns']):
        q,s,cap=c['source_face'],c['simplex_face'],c['capacity']
        assert q in sourcefaces and 0<=q<nf and 0<=s<nf and cap in(4,6,8)
        assert d['faces'][q]['kind']=='canonical_dotted_quad_aggregate'and d['faces'][q]['required']==1
        assert d['faces'][s]['kind']=='maximal_simplex_aggregate'and d['faces'][s]['required']==20
        assert c['integer_face_coefficients']==[[q,1],[s,-cap]]
        assert c['integer_vertex_capacities']==[]and c['vertex_cost']==0
        j=len(d['variables'])
        d['variables'].append(dict(index=j,kind='raw357_typed_simplex_coface_transfer',vertex=None,cost=0,source_index=k))
        d['faces'][q]['coefficients'].append([j,1]);d['faces'][s]['coefficients'].append([j,-cap])
        d['zeros'].append(dict(variable=j,simplex_face=s,capacity=cap,positive_coefficients=[[q,1]]))
        records.append(dict(variable=j,source_column=k,integer_face_coefficients=[[q,1],[s,-cap]],integer_vertex_capacities=[]))
    assert len(records)==16 and len({r['source_column']for r in records})==16
    assert A.digest(d['variables'][:nv])==hashes['variables']and A.digest(d['zeros'][:nz])==hashes['zeros']
    assert A.digest(d['budgets'])==hashes['budgets']
    assert(len(d['faces']),len(d['budgets']))==(nf,nb)
    for i,f in enumerate(d['faces']):
        assert(f['index'],f['kind'],f['required'])==oldmeta[i]
        assert [[j,a]for j,a in f['coefficients']if j<nv]==oldfaces[i]
    d['active']=sorted({j for f in d['faces']for j,a in f['coefficients']if a>0})
    assert set(d['active'])==oldactive|set(range(nv,nv+16))
    assert all(j<nv for j in d['repair_indices'])
    multi=parent['inherited_multi_vertex_columns'];assert len(multi)==14
    for r in multi:
        j=r['variable'];assert j<nv
        assert len(r['integer_vertex_capacities'])>1
    h=sha256()
    for k in('variables','faces','budgets','zeros'):h.update(A.digest(d[k]).encode());h.update(b'\n')
    prepared=dict(inherited_model_sha256=parent['model_sha256'],model_sha256=h.hexdigest(),
        inherited_counts=parent['counts'],new_columns=records,inherited_multi_vertex_columns=multi,
        counts=dict(inherited_variables=nv,variables=len(d['variables']),faces=nf,vertex_forests=nb,
            active_variables=len(d['active']),zero_budget_signed_variables=len(d['zeros']),
            inherited_positive_budget_signed_variables=52,preserved_multi_owner_columns=14,
            added_zero_budget_transfers=16,absent_simplex_rows_uninstalled=31,
            face_matrix_nonzeros=sum(len(f['coefficients'])for f in d['faces'])),
        changed_alpha_old_negative_terms=negative_terms,
        preservation='Every inherited variable, signed coefficient, face requirement, vertex budget and zero record is literal-preserved; only16 signed zero-budget columns are appended. No new face or vertex row.')
    return d,prepared,doc,old,oldprimal,olddual,deleted


def run():
    d,prepared,doc,old,primal,prior,deleted=build()
    gone={f for f,_ in deleted};alpha={f:F(x)for f,x in prior['face_multipliers']if f not in gone}
    mu={v:F(x)for v,x in prior['vertex_multipliers']}
    values=P.loads(d,alpha);costs=P.cost_rows(d)
    violations=[[j,str(values[j]),str(sum((a*mu.get(v,F())for v,a in cc),F()))]
        for j,cc in enumerate(costs)if values[j]>sum((a*mu.get(v,F())for v,a in cc),F())]
    numerator=sum((d['faces'][f]['required']*a for f,a in alpha.items()),F());denominator=sum(mu.values(),F())
    assert numerator==F(prior['objective_numerator'])-sum(F(x)for _,x in deleted)
    cert=dict(kind='EXACT_PARTIAL_MODEL_DUAL',face_multipliers=[[f,str(x)]for f,x in sorted(alpha.items())],
        vertex_multipliers=prior['vertex_multipliers'],objective_numerator=str(numerator),normalizer=str(denominator),
        bound=str(numerator/denominator),excludes_target12=numerator>12*denominator,
        recovery='Delete exactly nine positive raw357 source masses, retain every other face and vertex multiplier; check the full changed exact matrix.')
    P.verify(d,primal)
    certificates=[primal]
    if not violations:P.verify(d,cert);certificates.append(cert)
    paths=set(PINS)|{f'scripts/{TAG}.py',f'results/{MODEL}.json',f'scripts/{MODEL}.py',
        f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',f'results/{A.TAG}.md',str(Path(P.__file__).relative_to(ROOT))}
    out=dict(status='PASS_RAW357_COFACE_EXTENSION_EXACT_INTERVAL'if not violations else'PASS_RAW357_COFACE_EXTENSION_INFEASIBLE_DELETION_DIAGNOSTIC',
        sources={p:M.digest(p)for p in sorted(paths)},prepared=prepared,certificates=certificates,
        zeroed_source_multipliers=deleted,deleted_required_mass=str(sum(F(x)for _,x in deleted)),
        candidate=cert,deletion_violations=violations,full_matrix_column_checks=len(values),
        exact_interval=dict(lower=str(numerator/denominator),upper=primal['maximum'],
            gap=str(F(primal['maximum'])-numerator/denominator),lower_exceeds12=numerator>12*denominator)
            if not violations else None,
        optimizer_invocations=0,face_masses_free_in_model=True,
        scope='Changed finite exported model with16 complete typed coface transfers. Exact feasible lower and inherited upper certificates if deletion passes; no optimum, completed geometric coverage, H13 exclusion or dimension theorem.')
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['prepared']['counts'],out['exact_interval'],len(out['deletion_violations']),flush=True)
