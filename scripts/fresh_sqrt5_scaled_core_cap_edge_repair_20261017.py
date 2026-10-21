#!/usr/bin/env python3
"""Exact no-optimizer repair; only maximal-edge face masses increase."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_sqrt5_cross_coface_cap_full_model_20261017 as N

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_scaled_core_cap_edge_repair_20261017'
BASE='fresh_sqrt5_core23_25_spectator_incidence_20261017'
COVER='fresh_sqrt5_a3_d5_cap_local_cover_update_20261017'
JOINT='fresh_sqrt5_a3_d5_cap_joint_ports_20261017'
CAPS='fresh_sqrt5_scaled_core_cap_incompatibilities_20261017'


def dump(x):return json.dumps(x,sort_keys=True,separators=(',',':'))


def run():
    base,old,joint,newcaps=map(N.read,(BASE,COVER,JOINT,CAPS))
    assert old['complete_positive_simplex_forest_coverage']
    assert base['face_multipliers']==old['face_multipliers']==joint['face_multipliers']
    olda={i:F(a)for i,a in base['face_multipliers']};alpha=dict(olda)
    oldmu={v:F(a)for v,a in base['vertex_multipliers']}
    assert all(oldmu.get(v,F())>=F(a)for v,a in old['vertex_multipliers'])
    assert all(oldmu.get(v,F())>=F(a)for v,a in joint['vertex_multipliers'])
    for c in newcaps['columns']:
        load=sum((m*olda.get(fi,F())for fi,m in c['positive_face_coefficients']),F())
        fi=c['simplex_face'];alpha[fi]=max(alpha.get(fi,F()),load/c['capacity'])
    alpha={i:a for i,a in alpha.items()if a}
    changes={i:a-olda.get(i,F())for i,a in alpha.items()if a>olda.get(i,F())}
    assert changes=={80010:F(9,2500000000),82048:F(2991,200000000)}
    assert set(alpha)==set(olda)
    columns=old['block_orbit_columns'];features=old['positive_simplex_features']
    literals={c['index']:c for c in old['literal_block_orbits']}
    byface={f['face']:f for f in features};assert len(byface)==len(features)
    for fi in changes:
        f=byface[fi];assert f['simplex_dimension']==1
        assert sum(m for j,m in f['block_orbits'])==2
        assert all(columns[j]['block_size']==1 for j,m in f['block_orbits'])
    # A one-port block contributes zero to every selected joint S feature.
    # Check all literal k-subsets once; Q features and all other S masses stay fixed.
    joint_zero_checks=[]
    for k in(3,4,5):
        count=0
        for U in combinations(range(13),k):
            hit=int(0 in U)
            hist=[int(hit==r)for r in range(2,k+1)]
            assert hist==[0]*(k-1)
            assert int(hit>=2)==0 and max(hit-1,0)==0 and int(hit==3)==0
            count+=1
        joint_zero_checks.append([k,count])
    weights=defaultdict(F);oldweights=defaultdict(F)
    for f in features:
        for j,m in f['block_orbits']:
            weights[j]+=m*alpha[f['face']];oldweights[j]+=m*olda[f['face']]
    assert {j for j in weights if weights[j]!=oldweights[j]}=={2643,3137,4043,4044}
    mu=dict(oldmu);records=[];terms_count=0;additions=0
    for r in old['local_covers']:
        v=r['vertex'];js=[c['index']for c in r['literal_block_orbits']]
        covered={j:F()for j in js};mass=F();terms=[]
        for term in r['inherited_remapped_packing_mixture']+r['singleton_packing_increments']:
            a=F(term['mass']);assert a>0;used=0
            for j,mask in term['blocks']:
                assert j in covered and mask in literals[j]['literal_masks'] and not used&mask
                used|=mask;covered[j]+=a
            mass+=a;terms_count+=1;terms.append(term)
        assert mass==F(r['new_cover_mass'])
        assert all(F(a)==oldweights[j]for j,a in r['required_marginals'])
        assert all(covered[j]>=oldweights[j]for j in js)
        before=mass;increments=[]
        for j in js:
            deficit=max(F(),weights[j]-covered[j])
            if deficit:
                assert columns[j]['block_size']==1
                mask=literals[j]['literal_masks'][0];assert mask.bit_count()==1
                increments.append(dict(mass=str(deficit),blocks=[[j,mask]]))
                covered[j]+=deficit;mass+=deficit;additions+=1
        assert all(covered[j]>=weights[j]for j in js)
        mu[v]=max(mu.get(v,F()),mass)
        records.append(dict(vertex=v,forest=r['forest'],old_cover_mass=str(before),new_cover_mass=str(mass),
            literal_block_orbits=r['literal_block_orbits'],required_marginals=[[j,str(weights[j])]for j in js],
            covered_marginals=[[j,str(covered[j])]for j in js],
            inherited_remapped_packing_mixture=terms,singleton_packing_increments=increments))
    assert terms_count==3712 and additions<=4
    print('EDGE_LOCAL_COVER',terms_count,additions,changes,flush=True)
    d,prepared=N.build();assert json.loads(json.dumps(prepared))==N.read(N.TAG)['prepared']
    for fi in changes:
        assert d['faces'][fi]['required']==0
        assert d['faces'][fi]['kind']!='canonical_dotted_quad_aggregate'
    qrows={f['index']for f in d['faces']if f['kind']=='canonical_dotted_quad_aggregate'}
    assert all(alpha.get(i,F())==olda.get(i,F())for i in qrows)
    assert sum(i in qrows for i in alpha)==336
    pmap=dict(prepared['provenance_to_face'])
    for c in newcaps['columns']:
        assert sorted([pmap[p],m]for p,m in c['provenance_coefficients'])==c['positive_face_coefficients']
    loads=defaultdict(F)
    for fi,a in alpha.items():
        for j,m in d['faces'][fi]['coefficients']:loads[j]+=a*m
    # The complete local cover already provides all singleton-port increases.
    # Retain a direct inherited-matrix check rather than assuming that fact.
    for c in d['variables']:
        if c['vertex']is None:assert loads[c['index']]<=0
        else:assert loads[c['index']]<=c['cost']*mu.get(c['vertex'],F()),(c,loads[c['index']])
    cap_inputs=[]
    for tag in old['cap_packages']:
        x=N.read(tag)
        if 'columns'not in x:
            assert x['status']=='PASS_TWO_EXACT_A3_D5_CAP_WITNESSES'
            x=dict(columns=[dict(index=i,simplex_face=w['simplex_face'],capacity=w['capacity'],
                positive_face_coefficients=w['positive_face_coefficients'])for i,w in enumerate(x['witnesses'])])
        cap_inputs.append((tag,x))
    cap_inputs.append((CAPS,newcaps));assert sum(len(x['columns'])for tag,x in cap_inputs)==2071
    cap_checks=[]
    for tag,x in cap_inputs:
        for c in x['columns']:
            load=sum((m*alpha.get(fi,F())for fi,m in c['positive_face_coefficients']),F())
            rhs=c['capacity']*alpha.get(c['simplex_face'],F());assert load<=rhs
            cap_checks.append([tag,c['index'],str(load),str(rhs)])
    incidence_checks=[]
    for c in base['constraints']:
        load=sum((m*alpha.get(fi,F())for fi,m in c['positive_faces']),F())
        load-=sum((m*alpha.get(fi,F())for fi,m in c['negative_faces']),F())
        rhs=sum((m*mu.get(v,F())for v,m in c['vertex_capacity']),F());assert load<=rhs
        incidence_checks.append([c['index'],str(load),str(rhs)])
    assert len(incidence_checks)==52
    num=sum((d['faces'][fi]['required']*a for fi,a in alpha.items()),F())
    assert num==F(base['objective_numerator'])
    den=sum(mu.values(),F());bound=num/den
    paths=[Path(__file__),Path(N.__file__),ROOT/f'results/{N.TAG}.json']
    paths +=[ROOT/f'results/{t}.json'for t in(BASE,COVER,JOINT,CAPS)]
    paths +=[ROOT/f'results/{t}.json'for t in old['cap_packages']]
    return dict(status='PASS_EXACT_SCALED_CORE_EDGE_CAP_REPAIR_WITH_UNCHANGED_JOINT_LOADS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        inherited_model_sha256=prepared['model_sha256'],cap_package=CAPS,
        cap_packages=[t for t,x in cap_inputs],initial_bound=base['bound'],bound=str(bound),
        objective_numerator=str(num),normalizer=str(den),bound_above_12=bound>12,
        face_multipliers=[[i,str(a)]for i,a in sorted(alpha.items())],
        vertex_multipliers=[[v,str(a)]for v,a in sorted(mu.items())if a],
        simplex_mass_raises=[[i,str(olda[i]),str(alpha[i])]for i in sorted(changes)],
        vertex_mass_raises=[[v,str(oldmu.get(v,F())),str(a)]for v,a in sorted(mu.items())if a>oldmu.get(v,F())],
        positive_simplex_features=features,block_orbit_columns=columns,
        literal_block_orbits=old['literal_block_orbits'],local_covers=records,
        complete_positive_simplex_forest_coverage=True,cap_checks=cap_checks,
        forced_incidence_checks=incidence_checks,joint_singleton_zero_checks=joint_zero_checks,
        counts=dict(positive_faces=len(alpha),simplex_faces=len(features),actual_quad_faces=336,
            forests=len(records),block_orbits=len(columns),old_packing_terms=terms_count,
            added_singleton_terms=additions,raised_simplex_masses=len(changes),
            increased_vertex_masses=sum(a>oldmu.get(v,F())for v,a in mu.items()),
            inherited_signed_columns=len(d['zeros']),cap_checks=len(cap_checks),
            forced_incidence_checks=len(incidence_checks),joint_families_preserved=17),
        optimizer_invocations=0,repeated_large_joint_scans=0,
        inherited_joint_status='All changed S rows are maximal edges: every changed local block has size1, hence every3/4/5-port S histogram used by the17 inherited families is exactly zero. All Q masses and all other S masses are fixed, and mu only increases. Thus every inherited joint load is unchanged.',
        scope='Exact partial-model witness after15 new cap columns. Complete saved local packing mixtures and inherited matrix/cap/incidence rows are rechecked; all17 joint families persist by an explicit zero-feature identity. No optimization, whole-geometric-coverage claim or dimension theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(dump(run()));f=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(f.read_text())
    else:f.write_text(dump(d)+'\n')
    print(d['status'],d['counts'],float(F(d['bound'])),d['vertex_mass_raises'],flush=True)
