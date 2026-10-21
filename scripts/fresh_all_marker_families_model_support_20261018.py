#!/usr/bin/env python3
"""Saved exact dual support after the five-packet changed global solve.

Provenance and weighted demand only; no geometric incidence or gain claim.
Ordinary corner coverage and marker-family coverage remain distinct.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_all_marker_families_model_support_20261018'
MODEL='fresh_all_marker_families_dual_continuation_20261018';PARENT='fresh_current327_full_model_20261018'
FAMILIES=('fresh_raw76_extra_vertex_complete_family_20261018','fresh_raw56_raw57_extra_vertex_family_20261018',
 'fresh_raw108_complete_marker_family_20261018','fresh_raw112_complete_marker_family_20261018')
LOCAL='fresh_current324_new5073_local_covers_20261018';PATCH='fresh_new5073_raw76_high_mass_local_covers_20261018';THIRD='fresh_new5073_raw76_third_class_local_covers_20261018'
FORESTS='fresh_sqrt5_rank13_combined_disjoint_ports_20261017';FORCED='fresh_sqrt5_forced_base_simplex_domain_20261017'
ACTUAL='fresh_sqrt5_actual_provenance_dual_20261017';QUAD='fresh_sqrt5_rank13_packing_quad_separation_20261017'
TAGS=(MODEL,PARENT,*FAMILIES,LOCAL,PATCH,THIRD,FORESTS,FORCED,ACTUAL,QUAD)


def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    if isinstance(x,(tuple,list)):return list(map(enc,x))
    return x
def dual(d):return next(c for c in d['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')


def run():
    ds={s:read(s)for s in TAGS};d=ds[MODEL];old=ds[PARENT]
    assert d['status'].startswith('PASS_EXACT')and d['solver_invocations']==2
    cert=dual(d);prior=dual(old);a={f:F(x)for f,x in cert['face_multipliers']};oa={f:F(x)for f,x in prior['face_multipliers']}
    assert all(x>0 for x in a.values())
    pm={f:pi for pi,f in d['prepared']['provenance_to_face']};assert len(pm)==len(d['prepared']['provenance_to_face'])
    assert d['prepared']['provenance_to_face']==old['prepared']['provenance_to_face']
    faces=list(ds[FORESTS]['faces']);assert len(faces)==78335
    for r in ds[FORCED]['new_faces']:assert r['index']==len(faces);faces.append(r)
    assert len(faces)==82530
    provenance=ds[ACTUAL]['provenance_masses'];cores=ds[QUAD]['canonical_cores']
    ordinary=set(ds[LOCAL]['represented_classes']);ordinary.update(r['face']for r in ds[PATCH]['selected_faces']);ordinary.add(ds[THIRD]['selected_face']);assert len(ordinary)==327
    markers={}
    for stem in FAMILIES:
        for r in ds[stem]['family']:
            raw=r.get('raw_family',76 if stem==FAMILIES[0]else None);assert raw in(56,57,76,108,112)
            f=r['face'];assert f not in markers;markers[f]=raw
    assert len(markers)==2035
    rows=[]
    for f,mass in sorted(a.items()):
        if f in pm:
            pi=pm[f];r=provenance[pi];ci=r['canonical_prism_core']
            z=dict(face=f,kind='Q',provenance_kind='explicit',provenance=pi,source_context=r['source_context'],
                canonical_core=ci,raw_family=cores[ci]['source_core'],required=cores[ci]['required'])
        else:
            assert f<len(faces);ff=faces[f]
            if ff['kind']=='canonical_dotted_quad_aggregate':
                contexts=ff['source_contexts'];z=dict(face=f,kind='Q',required=ff['required'],source_contexts=contexts)
                if len(contexts)==1:
                    ctx,ci=contexts[0];z.update(provenance_kind='unique_saved_context',provenance=None,
                        source_context=ctx,canonical_core=ci,raw_family=cores[ci]['source_core'])
                else:z.update(provenance_kind='merged_unexpanded',raw_families=sorted({cores[ci]['source_core']for ctx,ci in contexts}))
            else:
                assert 'simplex'in ff['kind']or'edge'in ff['kind'];z=dict(face=f,kind=ff['kind'],required=ff['required'])
        if f in markers:assert z.get('raw_family')==markers[f]
        z.update(mass=mass,old_mass=oa.get(f,F()),change=mass-oa.get(f,F()),demand_mass=mass*z['required'],
            ordinary327=f in ordinary,marker_family=markers.get(f))
        rows.append(z)
    assert sum(r['demand_mass']for r in rows)==F(cert['objective_numerator'])
    assert F(cert['objective_numerator'])/F(cert['normalizer'])==F(cert['bound'])
    families={}
    for r in rows:
        if 'raw_family'not in r:continue
        k=r['raw_family'];z=families.setdefault(k,dict(raw_family=k,classes=0,mass=F(),demand_mass=F(),ordinary_classes=0,
            marker_classes=0,missing_classes=0,missing_mass=F(),missing_demand=F(),faces=[]))
        z['classes']+=1;z['mass']+=r['mass'];z['demand_mass']+=r['demand_mass'];z['faces'].append(r['face'])
        z['ordinary_classes']+=r['ordinary327'];z['marker_classes']+=r['marker_family']is not None
        if not r['ordinary327']and r['marker_family']is None:
            z['missing_classes']+=1;z['missing_mass']+=r['mass'];z['missing_demand']+=r['demand_mass']
    missing=sorted((r for r in rows if r['kind']=='Q'and'raw_family'in r and not r['ordinary327']and r['marker_family']is None),key=lambda r:(-r['demand_mass'],-r['mass'],r['face']))
    ranked=sorted(families.values(),key=lambda r:(-r['missing_demand'],r['raw_family']))
    marker_mass=[dict(raw_family=k,classes=sum(v==k for v in markers.values()),old_mass=sum(oa.get(f,F())for f,v in markers.items()if v==k),
        new_mass=sum(a.get(f,F())for f,v in markers.items()if v==k))for k in(56,57,76,108,112)]
    counts=dict(old_positive=len(oa),positive=len(a),newly_positive=len(set(a)-set(oa)),newly_zero=len(set(oa)-set(a)),
        changed_existing=sum(a[f]!=oa[f]for f in set(a)&set(oa)),ordinary327_positive=len(set(a)&ordinary),
        marker2035_positive=len(set(a)&set(markers)),newly_positive_outside_ordinary_and_markers=len(set(a)-set(oa)-ordinary-set(markers)),
        positive_merged_unexpanded=sum(r.get('provenance_kind')=='merged_unexpanded'for r in rows),
        unambiguous_missing_classes=len(missing),unambiguous_missing_families=sum(r['missing_classes']>0 for r in ranked))
    sources={f'scripts/{TAG}.py':digest(f'scripts/{TAG}.py')};sources.update({f'results/{s}.json':digest(f'results/{s}.json')for s in TAGS})
    return enc(dict(status='PASS_SAVED_ALL_MARKER_FAMILIES_MODEL_SUPPORT',sources=sources,counts=counts,
        bound=cert['bound'],objective_numerator=cert['objective_numerator'],normalizer=cert['normalizer'],
        solver_success=d['solver']['success'],face_masses_changed=a!=oa,marker_family_mass=marker_mass,
        positive_faces=rows,ranked_unambiguous_missing_families=ranked,ranked_unambiguous_missing_classes=missing,
        scope='Saved exact dual and provenance only. Ordinary327 and marker2035 are different incomplete constraints, not complete geometric coverage. Merged contexts remain ambiguous. Weighted-demand ranking is not an incidence theorem or a bound on possible gain. No geometry, local LP, new global solve or realization claim.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
    print('TOP_MISSING_FAMILIES',[(r['raw_family'],r['missing_classes'],r['missing_demand'])for r in out['ranked_unambiguous_missing_families'][:12]],flush=True)
