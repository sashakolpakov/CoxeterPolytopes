#!/usr/bin/env python3
"""Exact saved-only support inventory for the raw76 family fallback witness.

No geometry, local optimization or new global solve. Marker-family coverage
is kept separate from complete ordinary corner domains.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw76_complete_family_support_20261018'
MODEL='fresh_raw76_complete_family_full_model_20261018'
PARENT='fresh_current327_full_model_20261018'
FAMILY='fresh_raw76_extra_vertex_complete_family_20261018'
LOCAL='fresh_current324_new5073_local_covers_20261018'
PATCH='fresh_new5073_raw76_high_mass_local_covers_20261018'
THIRD='fresh_new5073_raw76_third_class_local_covers_20261018'
FORESTS='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'
FORCED='fresh_sqrt5_forced_base_simplex_domain_20261017'
ACTUAL='fresh_sqrt5_actual_provenance_dual_20261017'
QUAD='fresh_sqrt5_rank13_packing_quad_separation_20261017'
TAGS=(MODEL,PARENT,FAMILY,LOCAL,PATCH,THIRD,FORESTS,FORCED,ACTUAL,QUAD)
def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    return x
def dual(d):return next(c for c in d['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')

def run():
    ds={t:read(t)for t in TAGS};d=ds[MODEL];old=ds[PARENT]
    assert d['status'].startswith('PASS_EXACT')and d['solver_invocations']==1
    cert=dual(d);previous=dual(old)
    a={f:F(x)for f,x in cert['face_multipliers']};oa={f:F(x)for f,x in previous['face_multipliers']}
    assert all(x>0 for x in a.values())
    assert d['optimization_resolved']is False and d['changed_face_masses']is False
    assert a==oa
    pm={f:pi for pi,f in d['prepared']['provenance_to_face']}
    assert len(pm)==len(d['prepared']['provenance_to_face'])
    assert d['prepared']['provenance_to_face']==old['prepared']['provenance_to_face']
    faces=list(ds[FORESTS]['faces']);assert len(faces)==78335
    for r in ds[FORCED]['new_faces']:
        assert r['index']==len(faces);faces.append(r)
    assert len(faces)==82530
    pr=ds[ACTUAL]['provenance_masses'];cores=ds[QUAD]['canonical_cores']
    ordinary=set(ds[LOCAL]['represented_classes'])
    ordinary.update(r['face']for r in ds[PATCH]['selected_faces']);ordinary.add(ds[THIRD]['selected_face'])
    marker={r['face']for r in ds[FAMILY]['family']};assert len(ordinary)==327 and len(marker)==235
    rows=[]
    for fi,mass in sorted(a.items()):
        if fi in pm:
            pi=pm[fi];r=pr[pi];ci=r['canonical_prism_core'];ctx=r['source_context']
            z=dict(face=fi,kind='Q',provenance_kind='explicit',provenance=pi,
                   source_context=ctx,canonical_core=ci,raw_family=cores[ci]['source_core'],required=cores[ci]['required'])
        else:
            assert fi<len(faces);ff=faces[fi]
            if ff['kind']=='canonical_dotted_quad_aggregate':
                contexts=ff['source_contexts'];z=dict(face=fi,kind='Q',required=ff['required'],source_contexts=contexts)
                if len(contexts)==1:
                    ctx,ci=contexts[0];z.update(provenance_kind='unique_saved_context',provenance=None,
                        source_context=ctx,canonical_core=ci,raw_family=cores[ci]['source_core'])
                else:z.update(provenance_kind='merged_unexpanded',raw_families=sorted({cores[ci]['source_core']for ctx,ci in contexts}))
            else:
                assert 'simplex'in ff['kind']or'edge'in ff['kind']
                z=dict(face=fi,kind=ff['kind'],required=ff['required'])
        z.update(mass=mass,old_mass=oa.get(fi,F()),change=mass-oa.get(fi,F()),
                 demand_mass=mass*z['required'],ordinary327=fi in ordinary,marker235=fi in marker)
        if fi in marker:assert z.get('raw_family')==76
        rows.append(z)
    assert sum((r['demand_mass']for r in rows),F())==F(cert['objective_numerator'])
    assert F(cert['objective_numerator'])/F(cert['normalizer'])==F(cert['bound'])
    families={}
    for r in rows:
        if 'raw_family'not in r:continue
        k=r['raw_family'];z=families.setdefault(k,dict(raw_family=k,classes=0,mass=F(),demand_mass=F(),
            ordinary_classes=0,marker_classes=0,missing_classes=0,missing_mass=F(),missing_demand=F(),faces=[]))
        z['classes']+=1;z['mass']+=r['mass'];z['demand_mass']+=r['demand_mass'];z['faces'].append(r['face'])
        z['ordinary_classes']+=r['ordinary327'];z['marker_classes']+=r['marker235']
        if not r['ordinary327']and not r['marker235']:
            z['missing_classes']+=1;z['missing_mass']+=r['mass'];z['missing_demand']+=r['demand_mass']
    missing=[r for r in rows if r['kind']=='Q'and'raw_family'in r and not r['ordinary327']and not r['marker235']]
    missing.sort(key=lambda r:(-r['demand_mass'],-r['mass'],r['face']))
    ranked=sorted(families.values(),key=lambda r:(-r['missing_demand'],r['raw_family']))
    marker_before=sum((oa.get(f,F())for f in marker),F());marker_after=sum((a.get(f,F())for f in marker),F())
    counts=dict(old_positive=len(oa),positive=len(a),newly_positive=len(set(a)-set(oa)),newly_zero=len(set(oa)-set(a)),
        newly_positive_outside327=len(set(a)-set(oa)-ordinary),ordinary327_positive=len(set(a)&ordinary),
        marker235_positive=len(set(a)&marker),positive_Q=sum(r['kind']=='Q'for r in rows),
        positive_merged_unexpanded=sum(r.get('provenance_kind')=='merged_unexpanded'for r in rows),
        unambiguous_missing_classes=len(missing),unambiguous_missing_families=sum(r['missing_classes']>0 for r in ranked))
    sources={f'scripts/{TAG}.py':digest(f'scripts/{TAG}.py')}
    sources.update({f'results/{t}.json':digest(f'results/{t}.json')for t in TAGS})
    return enc(dict(status='PASS_SAVED_RAW76_FAMILY_MODEL_SUPPORT',sources=sources,counts=counts,
        bound=cert['bound'],objective_numerator=cert['objective_numerator'],normalizer=cert['normalizer'],
        marker_mass_before=marker_before,marker_mass_after=marker_after,positive_faces=rows,
        ranked_unambiguous_missing_families=ranked,ranked_unambiguous_missing_classes=missing,
        scope='Saved exact fallback dual/provenance inventory only. Face masses are unchanged after the timed-out proposal; this is not a new optimized support. Ordinary327 and marker235 cover different constraints; neither is a complete actual face domain. Merged contexts are not assigned a unique raw family. Rankings measure weighted demand, not a possible improvement or a geometric obstruction. No geometry, local LP or global solve.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(out.read_text())==result
    else:assert not out.exists();out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],result['counts'],flush=True)
    print('MARKER_MASS',result['marker_mass_before'],result['marker_mass_after'],flush=True)
    print('TOP_MISSING_FAMILIES',[(r['raw_family'],r['missing_classes'],r['missing_demand'])for r in result['ranked_unambiguous_missing_families'][:12]],flush=True)
    print('TOP_MISSING_CLASSES',[(r['face'],r['raw_family'],r['source_context'],r['canonical_core'],r['mass'],r['required'])for r in result['ranked_unambiguous_missing_classes'][:12]],flush=True)
