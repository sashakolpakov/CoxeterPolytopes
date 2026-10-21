#!/usr/bin/env python3
"""Prepared sparse append of 795 current-union and raw33/295 corner capacities.

The inherited 52 signed vertex-cost columns remain ordinary budget columns;
new columns are all nonnegative face payments with a positive vertex cost.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import lcm
from pathlib import Path
import shutil

import fresh_raw23_mixed_full_model_20261017 as BASE

ROOT=BASE.ROOT
TAG='fresh_raw33_295_mixed_full_model_20261017'
OUT=ROOT/f'results/{TAG}.json'
INPUTS=(
    ('fresh_current220_actual_corner_union_20261017','columns','positive_face_coefficients',443),
    ('fresh_raw33_295_mixed_actual_corners_20261017','columns','positive_face_coefficients',352),
)
ENGINE=BASE.ENGINE
dump,read,digest=BASE.dump,BASE.read,BASE.digest


def build():
    assert shutil.disk_usage(ROOT).free>200*1024**2
    d,inherited=BASE.build()
    assert json.loads(dump(inherited))==read(BASE.TAG)['prepared']
    nv,nf,nb,nz=len(d['variables']),len(d['faces']),len(d['budgets']),len(d['zeros'])
    assert (nv,nf,nb,nz)==(438266,82806,6558,96176)
    hashes=digest(d['variables']),digest(d['zeros'])
    old_faces=[list(map(list,f['coefficients']))for f in d['faces']]
    old_budgets=[list(map(list,b['coefficients']))for b in d['budgets']]
    old_meta=[(f['index'],f['kind'],f['required'])for f in d['faces']]
    old_active=set(d['active'])
    records={};counts=[];packages=[];bound_hashes={}
    for tag,field,coeff_field,expected in INPUTS:
        p=read(tag);packages.append(p);counts.append([tag,len(p[field])])
        for path,expected_hash in p['sources'].items():
            if path not in bound_hashes:
                h=sha256()
                with (ROOT/path).open('rb')as f:
                    for block in iter(lambda:f.read(1<<20),b''):h.update(block)
                bound_hashes[path]=h.hexdigest()
            assert bound_hashes[path]==expected_hash,(tag,path)
        assert len(p[field])==expected
        if 'fixed_witness_repair' in p:assert p['fixed_witness_repair']['face_masses_unchanged']
        for c in p[field]:
            v=c['vertex'];cost=F(c['vertex_cost'])
            terms=tuple(sorted((fi,F(m))for fi,m in c[coeff_field]))
            assert cost>0 and terms and all(m>0 for fi,m in terms)
            assert len({fi for fi,m in terms})==len(terms)
            assert 0<=v<nb
            assert all(0<=fi<nf for fi,m in terms)
            key=v,cost,terms
            records.setdefault(key,[]).append([tag,c['index']])
    assert counts==[[tag,n]for tag,field,coeff_field,n in INPUTS]
    parent_dual=next(c for c in read(BASE.TAG)['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    previous=parent_dual
    for package in packages:
        assert package['face_multipliers']==parent_dual['face_multipliers']
        assert package['initial_vertex_multipliers']==previous['vertex_multipliers']
        prior_mu={i:F(a)for i,a in previous['vertex_multipliers']}
        next_mu={i:F(a)for i,a in package['vertex_multipliers']}
        assert all(next_mu.get(i,F())>=a for i,a in prior_mu.items())
        assert all(a>=0 for a in next_mu.values())
        assert F(package['fixed_witness_repair']['objective_numerator'])/sum(F(a)for i,a in package['vertex_multipliers'])==F(package['fixed_witness_repair']['bound'])
        previous=package
    assert packages[-1]['fixed_witness_repair']['bound_below12']
    columns=[]
    for ci,((v,cost,terms),origins)in enumerate(sorted(records.items())):
        den=lcm(cost.denominator,*(m.denominator for fi,m in terms))
        j=len(d['variables']);sc=int(cost*den)
        assert F(sc,den)==cost and sc>0
        d['variables'].append(dict(index=j,kind='raw33_295_mixed_actual_corner_capacity',vertex=v,cost=sc,source_index=ci))
        d['budgets'][v]['coefficients'].append([j,sc])
        for fi,m in terms:
            sm=int(m*den);assert F(sm,den)==m and sm>0
            d['faces'][fi]['coefficients'].append([j,sm])
        columns.append(dict(variable=j,vertex=v,vertex_cost=str(cost),positive_face_coefficients=[[fi,str(m)]for fi,m in terms],
            positive_denominator=den,scaled_vertex_cost=sc,sources=origins))
    assert len(d['faces'])==nf and len(d['budgets'])==nb and len(d['zeros'])==nz
    assert digest(d['variables'][:nv])==hashes[0] and digest(d['zeros'])==hashes[1]
    for i,f in enumerate(d['faces']):
        assert (f['index'],f['kind'],f['required'])==old_meta[i]
        assert [[j,m]for j,m in f['coefficients']if j<nv]==old_faces[i]
        assert all(m>0 for j,m in f['coefficients']if j>=nv)
    for i,b in enumerate(d['budgets']):
        assert [[j,m]for j,m in b['coefficients']if j<nv]==old_budgets[i]
    d['active']=sorted({j for f in d['faces']for j,m in f['coefficients']if m>0})
    assert old_active<=set(d['active'])
    eligible=set(d['repair_indices'])
    assert all(j<nv for j in eligible)
    assert all(sum(m for j,m in f['coefficients']if j in eligible)>0 for f in d['faces'])
    assert all(m>0 for f in d['faces']for j,m in f['coefficients']if j in eligible)
    h=sha256()
    for group in(d['variables'],d['faces'],d['budgets'],d['zeros']):h.update(digest(group).encode());h.update(b'\n')
    paths=[Path(__file__),Path(BASE.__file__),ROOT/f'results/{BASE.TAG}.json',ROOT/f'results/{BASE.TAG}.md']
    for tag,field,coeff_field,expected in INPUTS:
        paths +=[ROOT/f'{folder}/{tag}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    prepared=dict(sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        dimension=13,target=12,inherited_model_sha256=inherited['model_sha256'],model_sha256=h.hexdigest(),
        provenance_to_face=inherited['provenance_to_face'],new_columns=columns,
        counts=dict(faces=nf,vertex_forests=nb,variables=len(d['variables']),inherited_variables=nv,
            active_variables=len(d['active']),zero_budget_signed_variables=nz,vertex_cost_signed_variables=52,
            exported_input_columns=sum(n for tag,field,coeff_field,n in INPUTS),distinct_added_columns=len(columns),exact_duplicate_exports=sum(n for tag,field,coeff_field,n in INPUTS)-len(columns),
            face_matrix_nonzeros=sum(len(f['coefficients'])for f in d['faces'])),input_column_counts=counts,
        preservation='Every inherited438266 variable,82806 face obligation,6558 vertex budget and96176 zero-budget signed column retains its exact coefficients. The52 signed vertex-cost incidence columns are preserved and never converted to zero-budget transfers. All795 exported current220-union and raw33/295 mixed actual-corner capacities are represented; only exactly identical new (vertex,cost,face-coefficient) tuples are merged. Every variable remains nonnegative, including those with signed face coefficients.',
        recovery='The inherited exact primal repair changes only strictly positive old pair/port coefficients. The inherited dual recovery raises simplex masses to meet zero-budget transfers and recomputes every vertex-cost load, including the52 signed columns. New columns have only positive entries and positive vertex costs.',
        finite_family_scope='The finite exported matrix is enlarged. Complete local packing cones and the millions of implicit joint-family constraints are not all columns of this matrix. Their fixed-alpha certificates cannot be reused at newly optimized face masses without reconstruction or a valid marginal comparison.',
        scope='One bounded finite partial-model optimization. A target12 primal still requires missing face, forest, field and charge coverage before any dimension theorem. An exact dual obstructs this finite charge family; it is not simultaneous ambient geometric realization.')
    return d,prepared


def verify(d,c):
    BASE.verify(d,c)
    if c['kind']=='EXACT_PARTIAL_MODEL_DUAL':
        alpha={i:F(a)for i,a in c['face_multipliers']};mu={v:F(a)for v,a in c['vertex_multipliers']}
        variables={x['index']:x for x in d['variables']if x['kind']=='raw33_295_mixed_actual_corner_capacity'}
        loads={j:F()for j in variables}
        for fi,a in alpha.items():
            for j,m in d['faces'][fi]['coefficients']:
                if j in loads:loads[j]+=m*a
        assert all(loads[j]<=x['cost']*mu.get(x['vertex'],F())for j,x in variables.items())


if __name__=='__main__':
    p=argparse.ArgumentParser();m=p.add_mutually_exclusive_group(required=True)
    m.add_argument('--prepare',action='store_true');m.add_argument('--solve',action='store_true');m.add_argument('--verify',action='store_true')
    args=p.parse_args();d,prepared=build();prepared=json.loads(dump(prepared))
    print('PREPARED_RAW33_295_MIXED_ACTUAL_CORNER_FULL_MODEL',prepared['counts'],prepared['model_sha256'],flush=True)
    if args.verify:
        saved=json.loads(OUT.read_text());assert saved['prepared']==prepared
        for c in saved.get('certificates',[]):verify(d,c)
        print('PASS_SAVED_RAW33_295_MIXED_ACTUAL_CORNER_FULL_MODEL',saved['status'],flush=True)
    elif args.prepare:
        if OUT.exists():assert json.loads(OUT.read_text())['solver_invocations']==0
        OUT.write_text(dump(dict(status='PREPARED_RAW33_295_MIXED_ACTUAL_CORNER_FULL_MODEL',prepared=prepared,
            solver_invocations=0,solver=None,certificates=[]))+'\n')
    else:
        saved=json.loads(OUT.read_text());assert saved['prepared']==prepared and saved['solver_invocations']==0
        assert shutil.disk_usage(ROOT).free>200*1024**2
        saved.update(status='ONE_BOUNDED_SOLVE_STARTED',solver_invocations=1);OUT.write_text(dump(saved)+'\n')
        solver,certificates=ENGINE.optimize(d)
        for c in certificates:verify(d,c)
        status='PASS_EXACT_PARTIAL_PRIMAL_TARGET12'if any(c.get('meets_target12')for c in certificates)else \
            'PASS_EXACT_RAW33_295_MIXED_ACTUAL_CORNER_MODEL_CERTIFICATES'if certificates else'UNRESOLVED_SINGLE_BOUNDED_SOLVE'
        saved.update(status=status,solver=solver,certificates=certificates);OUT.write_text(dump(saved)+'\n')
        print(status,[(c['kind'],c.get('maximum',c.get('bound')))for c in certificates],flush=True)
