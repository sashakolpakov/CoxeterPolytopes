#!/usr/bin/env python3
"""One bounded finite-model solve with all exported caps and actual corners.

The forced-quad incidence columns spend a vertex budget despite having a
negative simplex coefficient. They are not zero-budget transfer columns.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import lcm
from pathlib import Path
import json
import shutil

import fresh_sqrt5_cross_coface_cap_full_model_20261017 as BASE

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_actual_corner_enlarged_model_20261017'
OUT=ROOT/f'results/{TAG}.json'
WITNESS='fresh_sqrt5_scaled_core_cap_edge_repair_20261017'
INCIDENCE='fresh_sqrt5_core23_25_spectator_incidence_20261017'
FAMILY='fresh_dominant_cap_compatibility_family_20261017'
JOINT='fresh_sqrt5_profile_clique_joint_ports_20261017'
ENGINE=BASE.N.M
digest=BASE.N.digest


def dump(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def read(tag):return BASE.read(tag)


def cap_key(c):
    positive=tuple(sorted((fi,F(m))for fi,m in c['positive_face_coefficients']))
    assert len({fi for fi,m in positive})==len(positive)
    assert positive and all(m>0 for fi,m in positive)
    assert F(c['capacity'])>0
    return c['simplex_face'],F(c['capacity']),positive


def build():
    assert shutil.disk_usage(ROOT).free>200*1024**2,'Need200MiB free for compact certificate and solver headroom'
    d,inherited=BASE.build()
    saved=read(BASE.TAG)
    assert json.loads(dump(inherited))==saved['prepared']
    witness,incidence,family,joint=map(read,(WITNESS,INCIDENCE,FAMILY,JOINT))
    oldvars,oldfaces,oldzeros=len(d['variables']),len(d['faces']),len(d['zeros'])
    assert (oldvars,oldfaces,oldzeros)==(428768,82806,94722)
    oldhash=(digest(d['variables']),digest(d['zeros']))
    oldfacescopy=[list(map(list,f['coefficients']))for f in d['faces']]
    oldbudgets=[list(map(list,b['coefficients']))for b in d['budgets']]
    oldmeta=[(f['index'],f['kind'],f['required'])for f in d['faces']]
    oldactive=set(d['active']);pmap=dict(inherited['provenance_to_face'])
    assert len(pmap)==1858 and len(set(pmap.values()))==1858
    actual_rows=set(pmap.values());sources=[];cap_records={};raw_cap_count=0
    for tag in witness['cap_packages']:
        pack=read(tag);sources.append(tag)
        cc=pack.get('columns',pack.get('witnesses'))
        assert cc is not None
        for ci,c in enumerate(cc):
            k=cap_key(c);raw_cap_count+=1
            if 'provenance_coefficients'in c:
                assert sorted((pmap[j],F(m))for j,m in c['provenance_coefficients'])==list(k[2])
            else:
                assert 'provenance'in c and k[2]==((pmap[c['provenance']],F(1)),)
            cap_records.setdefault(k,[]).append([tag,c.get('index',ci)])
    assert raw_cap_count==2071 and len(cap_records)==1695
    oldcap=read(BASE.CAP);oldkeys={cap_key(c)for c in oldcap['columns']}
    assert len(oldkeys)==241 and oldkeys<=set(cap_records)
    newkeys=sorted(set(cap_records)-oldkeys)
    assert len(newkeys)==1454
    blocks=[];scalings=[];cap_bindings=[];vertex_signed=[]

    def add_positive(kind,source,ci,vertex,cost,positive,negative=()):
        """One nonnegative payment variable; its vertex cost is positive.

        Negative S coefficients are permitted for the proved forced-face
        incidence family. No such variable is appended to d['zeros'].
        """
        cost=F(cost);positive=[(fi,F(m))for fi,m in positive]
        negative=[(fi,F(m))for fi,m in negative]
        assert cost>0 and positive and all(m>0 for fi,m in positive+negative)
        assert len({fi for fi,m in positive+negative})==len(positive)+len(negative)
        den=lcm(cost.denominator,*(m.denominator for fi,m in positive+negative))
        j=len(d['variables']);scaled_cost=int(cost*den)
        assert F(scaled_cost,den)==cost and 0<=vertex<len(d['budgets'])
        d['variables'].append(dict(index=j,kind=kind,vertex=vertex,cost=scaled_cost,source_index=ci))
        d['budgets'][vertex]['coefficients'].append([j,scaled_cost])
        for fi,m in positive:
            assert 0<=fi<oldfaces
            d['faces'][fi]['coefficients'].append([j,int(m*den)])
        for fi,m in negative:
            assert d['faces'][fi]['kind']in('maximal_simplex_aggregate','maximal_edge_incidence')
            d['faces'][fi]['coefficients'].append([j,-int(m*den)])
        scalings.append(dict(variable=j,source=source,source_index=ci,positive_denominator=den,
            original_vertex_cost=str(cost),scaled_vertex_cost=scaled_cost,
            positive_face_count=len(positive),negative_face_count=len(negative)))
        if negative:vertex_signed.append(dict(variable=j,vertex=vertex,cost=scaled_cost,
            positive_coefficients=[[fi,int(m*den)]for fi,m in positive],
            negative_coefficients=[[fi,int(m*den)]for fi,m in negative]))
        return j

    start=len(d['variables'])
    for ci,k in enumerate(newkeys):
        fi,capacity,positive=k
        assert d['faces'][fi]['kind']in('maximal_simplex_aggregate','maximal_edge_incidence')
        assert all(q in actual_rows and d['faces'][q]['kind']=='canonical_dotted_quad_aggregate'for q,m in positive)
        den=lcm(capacity.denominator,*(m.denominator for q,m in positive))
        cap=int(capacity*den);terms=[[q,int(m*den)]for q,m in positive]
        assert cap>0 and F(cap,den)==capacity
        j=len(d['variables'])
        d['variables'].append(dict(index=j,kind='exported_geometric_cap_transfer',vertex=None,cost=0,source_index=ci))
        d['faces'][fi]['coefficients'].append([j,-cap])
        for q,m in terms:d['faces'][q]['coefficients'].append([j,m])
        d['zeros'].append(dict(variable=j,simplex_face=fi,capacity=cap,positive_coefficients=terms))
        cap_bindings.append(dict(variable=j,sources=cap_records[k],positive_denominator=den,
            simplex_face=fi,capacity=str(capacity),positive_face_coefficients=[[q,str(m)]for q,m in positive]))
    blocks.append(dict(kind='exported_geometric_cap_transfer',start=start,count=len(newkeys)))

    start=len(d['variables']);assert len(incidence['constraints'])==52
    for c in incidence['constraints']:
        assert len(c['vertex_capacity'])==1
        v,cost=c['vertex_capacity'][0]
        for q,m in c['positive_faces']:assert q in actual_rows
        add_positive('forced_quad_incidence_with_vertex_cost',INCIDENCE,c['index'],v,cost,
                     c['positive_faces'],c['negative_faces'])
    blocks.append(dict(kind='forced_quad_incidence_with_vertex_cost',source=INCIDENCE,start=start,count=52))

    start=len(d['variables']);assert len(family['columns'])==294
    for c in family['columns']:
        assert c['vertex_cost']==1
        assert all(q in actual_rows for q,m in c['positive_face_coefficients'])
        add_positive('actual_raw23_corner_clique',FAMILY,c['index'],c['vertex'],1,c['positive_face_coefficients'])
    blocks.append(dict(kind='actual_raw23_corner_clique',source=FAMILY,start=start,count=294))

    start=len(d['variables']);assert len(joint['critical_columns'])==98
    for c in joint['critical_columns']:
        add_positive('saved_profile_clique_joint_port_cut',JOINT,c['index'],c['vertex'],
                     c['vertex_cost'],c['positive_face_coefficients'])
    blocks.append(dict(kind='saved_profile_clique_joint_port_cut',source=JOINT,start=start,count=98))

    assert (len(d['variables']),len(d['faces']),len(d['budgets']),len(d['zeros']))==(430666,82806,6558,96176)
    assert digest(d['variables'][:oldvars])==oldhash[0] and digest(d['zeros'][:oldzeros])==oldhash[1]
    signed_vertices={c['variable']for c in vertex_signed}
    assert len(signed_vertices)==52 and signed_vertices.isdisjoint(c['variable']for c in d['zeros'])
    for i,f in enumerate(d['faces']):
        assert (f['index'],f['kind'],f['required'])==oldmeta[i]
        co=Counter()
        for j,m in f['coefficients']:co[j]+=m
        f['coefficients']=[[j,m]for j,m in sorted(co.items())if m]
        assert [[j,m]for j,m in f['coefficients']if j<oldvars]==oldfacescopy[i]
        assert all(m>0 or j in signed_vertices for j,m in f['coefficients']if d['variables'][j]['vertex']is not None)
    for v,b in enumerate(d['budgets']):
        assert [[j,m]for j,m in b['coefficients']if j<oldvars]==oldbudgets[v]
    d['active']=sorted({j for f in d['faces']for j,m in f['coefficients']if m>0})
    assert oldactive<=set(d['active'])
    eligible=set(d['repair_indices']);assert not eligible&signed_vertices
    assert all(sum(m for j,m in f['coefficients']if j in eligible)>0 for f in d['faces'])
    assert all(m>0 for f in d['faces']for j,m in f['coefficients']if j in eligible)
    for c in d['zeros']:
        assert d['variables'][c['variable']]['vertex']is None
        assert d['faces'][c['simplex_face']]['kind']in('maximal_simplex_aggregate','maximal_edge_incidence')
        assert all(d['faces'][fi]['kind']=='canonical_dotted_quad_aggregate'for fi,m in c['positive_coefficients'])

    h=sha256()
    for group in(d['variables'],d['faces'],d['budgets'],d['zeros']):h.update(digest(group).encode());h.update(b'\n')
    paths=[Path(__file__),Path(BASE.__file__),ROOT/f'results/{BASE.TAG}.json']
    for source in [WITNESS,INCIDENCE,FAMILY,JOINT]+sources:
        paths +=[ROOT/f'{folder}/{source}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths=sorted(set(paths))
    prepared=dict(sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        dimension=13,target=12,inherited_model_sha256=inherited['model_sha256'],model_sha256=h.hexdigest(),
        provenance_to_face=inherited['provenance_to_face'],counts=dict(faces=len(d['faces']),vertex_forests=len(d['budgets']),
            variables=len(d['variables']),inherited_variables=oldvars,active_variables=len(d['active']),
            zero_budget_signed_variables=len(d['zeros']),vertex_cost_signed_variables=len(vertex_signed),
            input_cap_records=raw_cap_count,unique_input_caps=len(cap_records),cap_additions=len(newkeys),
            forced_incidence_columns=52,actual_corner_columns=294,saved_joint_cuts=98,
            face_matrix_nonzeros=sum(len(f['coefficients'])for f in d['faces'])),
        added_column_blocks=blocks,new_cap_bindings=cap_bindings,exact_vertex_column_scalings=scalings,
        vertex_cost_signed_columns=vertex_signed,
        preservation='All82806 old face obligations,6558 budgets,428768 variables and94722 zero-budget transfers retain their coefficients. All2071 exported cap records are present through1695 exact distinct columns, including241 already inherited. No frozen source is changed.',
        recovery='The52 forced-incidence columns have a negative S term and positive single-vertex cost. They are handled by ordinary per-vertex dual loads, never as zero-budget constraints. Zero-budget Q-positive/S-negative columns still admit monotone S-mass repair. Recomputed mu handles every vertex-cost column, including signed ones. Exact primal repair changes only inherited strictly positive pair/port columns.',
        finite_family_scope='All exported historical constraints and new selected actual-corner charges are retained. The implicit complete local simplex-packing cones and all2.7M joint-family inequalities are not claimed to be fully represented by this finite matrix. Their old fixed-alpha certificates cannot be imported for new optimized masses.',
        scope='One enlarged finite partial model. Any exact target12 primal still needs missing face/forest and field coverage and appropriate full charge verification; no dimension theorem follows directly.')
    return d,prepared


def verify(d,c):
    BASE.N.verify(d,c)
    if c['kind']=='EXACT_PARTIAL_MODEL_DUAL':
        alpha={fi:F(a)for fi,a in c['face_multipliers']};mu={v:F(a)for v,a in c['vertex_multipliers']}
        # Recount all nonzero-cost signed columns independently of the
        # monotone zero-budget repair mechanism.
        ids={x['index']:x for x in d['variables']if x['kind']=='forced_quad_incidence_with_vertex_cost'}
        load={j:F()for j in ids}
        for fi,a in alpha.items():
            for j,m in d['faces'][fi]['coefficients']:
                if j in load:load[j]+=m*a
        assert len(ids)==52
        assert all(load[j]<=x['cost']*mu.get(x['vertex'],F())for j,x in ids.items())


if __name__=='__main__':
    p=argparse.ArgumentParser();mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--solve',action='store_true');mode.add_argument('--verify',action='store_true')
    args=p.parse_args();d,prepared=build();prepared=json.loads(dump(prepared))
    print('PREPARED_ACTUAL_CORNER_ENLARGED_MODEL',prepared['counts'],prepared['model_sha256'],flush=True)
    if args.verify:
        saved=json.loads(OUT.read_text());assert saved['prepared']==prepared
        for c in saved.get('certificates',[]):verify(d,c)
        print('PASS_SAVED_ACTUAL_CORNER_ENLARGED_MODEL',saved['status'],flush=True)
    elif args.prepare:
        if OUT.exists():assert json.loads(OUT.read_text())['solver_invocations']==0
        OUT.write_text(dump(dict(status='PREPARED_ACTUAL_CORNER_ENLARGED_MODEL',prepared=prepared,
            solver_invocations=0,solver=None,certificates=[]))+'\n')
    else:
        saved=json.loads(OUT.read_text());assert saved['prepared']==prepared and saved['solver_invocations']==0
        assert shutil.disk_usage(ROOT).free>200*1024**2
        saved.update(status='ONE_BOUNDED_SOLVE_STARTED',solver_invocations=1);OUT.write_text(dump(saved)+'\n')
        solver,cert=ENGINE.optimize(d)
        for c in cert:verify(d,c)
        status='PASS_EXACT_PARTIAL_PRIMAL_TARGET12'if any(c.get('meets_target12')for c in cert)else \
            'PASS_EXACT_ACTUAL_CORNER_MODEL_CERTIFICATES'if cert else 'UNRESOLVED_SINGLE_BOUNDED_SOLVE'
        saved.update(status=status,solver=solver,certificates=cert);OUT.write_text(dump(saved)+'\n')
        print(status,[(c['kind'],c.get('maximum',c.get('bound')))for c in cert],flush=True)
