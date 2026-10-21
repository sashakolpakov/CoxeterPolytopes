#!/usr/bin/env python3
"""Exact cap multiplicity obstruction across orthogonal A1 cofaces.

The geometry is in the report: fixed core contacts determine a common
projection v; 4<tau<6 gives -1<v.v<0 for caps at distinct A1 copies.
"""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_cross_coface_cap_capacity_20261017'
INPUTS=dict(slices='fresh_sqrt5_labelled_coface_base_packing_20261017',
    incidence='fresh_sqrt5_prism_simplex_base_incidence_20261017',
    actual='fresh_sqrt5_actual_provenance_dual_20261017',
    quad='fresh_sqrt5_rank13_packing_quad_separation_20261017',
    raw='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017')


def sign(a,b):
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if (a>0)==(b>0):return (a>0)-(a<0)
    return ((a>0)-(a<0))*((a*a>5*b*b)-(a*a<5*b*b))


def run():
    inputs={k:json.loads((ROOT/f'results/{tag}.json').read_text()) for k,tag in INPUTS.items()}
    src,inc,act,quad,raw=(inputs[k] for k in ('slices','incidence','actual','quad','raw'))
    assert src['status']=='PASS_EXACT_LABELLED_COFACE_SUBSET_PACKING_DIAGNOSTIC'
    alpha={i:F(a) for i,a in src['face_multipliers']}
    cores={c['core']:c for c in src['core_symmetry_and_subset_orbits']}
    profiles={r['base_type']:r for r in src['coface_slice_demands']}
    a1s={r['base_type']:r for r in src['coface_slice_demands'] if r['coface_label']==['finite','A1',0]}
    templates={t['canonical_core']:t for t in inc['prism_templates']}
    groups=defaultdict(list)
    occurrences=[]
    counts=Counter()
    for B in src['actual_bases']:
        t=B['base_type'];p=profiles[t];c=cores[p['core']]
        if B['coface_location']!=['finite','A1',0]:continue
        counts['A1_base_occurrences']+=1
        if len(c['automorphisms'])!=1:continue
        assert c['automorphisms'][0]==list(range(c['nodes']))
        counts['trivial_core_automorphism']+=1
        ci=act['provenance_masses'][B['provenance']]['canonical_prism_core']
        si=quad['canonical_cores'][ci]['source_core']
        rec=raw['records'][si]
        tpl=next(x for x in templates[ci]['bases'] if x['deleted_normal']==B['deleted'])
        b=B['deleted'];q=B['coface_root'];mp=dict(tpl['core_isomorphism'])
        assert tpl['base_core']==p['core'] and tpl['opposite_base_facet']==q
        if set(range(rec['nodes']))!={b,q}|set(mp):continue
        assert b!=q and b not in mp and q not in mp
        assert sorted(mp.values())==list(range(c['nodes']))
        assert tpl['finite_extra']==['A1']
        counts['sole_A1_extra']+=1
        es=dict(zip(combinations(range(rec['nodes']),2),rec['row']['labels']))
        def label(i,j):return es[tuple(sorted((i,j)))]
        assert all(label(q,i)==2 for i in mp)
        if label(b,q)!=0:continue
        tau=tuple(map(F,rec['parameters'][','.join(map(str,sorted((b,q))))]))
        if not (sign(tau[0]-4,tau[1])>0 and sign(6-tau[0],-tau[1])>0):continue
        counts['forbidden_distinct_copy_pair']+=1
        projection_norm=(F(2)-tau[0]/2,-tau[1]/2)
        assert sign(*projection_norm)<0 and sign(projection_norm[0]+1,projection_norm[1])>0
        contacts=[]
        for old,new in sorted(mp.items(),key=lambda x:x[1]):
            m=label(b,old)
            parameters=rec['parameters'][','.join(map(str,sorted((b,old))))] if m==0 else None
            contacts.append((new,m,None if parameters is None else tuple(map(str,map(F,parameters)))))
        key=(t,tuple(contacts),tuple(B['canonical_vertices']),tuple(map(str,tau)))
        occurrence=dict(base_index=B['index'],base_type=t,provenance=B['provenance'],
            face=B['face'],canonical_prism_core=ci,source_core=si,deleted_normal=b,
            original_coface_root=q,core=p['core'],core_isomorphism=sorted(map(list,mp.items())),
            base_vertices=B['canonical_vertices'],core_contacts=contacts,
            dotted_squared_norm2_product=list(map(str,tau)),
            common_projection_norm=list(map(str,projection_norm)),
            A1_copies=a1s[t]['coface_multiplicity'])
        occurrences.append(occurrence)
        groups[key].append(occurrence)
    columns=[]
    required=defaultdict(F)
    for key,bs in sorted(groups.items()):
        t=key[0];p=profiles[t]
        co=Counter(b['provenance'] for b in bs)
        face_co=Counter(b['face'] for b in bs)
        load=sum((m*alpha.get(fi,F()) for fi,m in face_co.items()),F())
        old=alpha.get(p['simplex_face'],F())
        required[p['simplex_face']]=max(required[p['simplex_face']],load)
        columns.append(dict(index=len(columns),base_type=t,core=p['core'],
            simplex_face=p['simplex_face'],simplex_dimension=p['simplex_dimension'],
            capacity=1,A1_copies=a1s[t]['coface_multiplicity'],
            source_bases=[b['base_index'] for b in bs],
            provenance_coefficients=sorted(map(list,co.items())),
            positive_face_coefficients=sorted(map(list,face_co.items())),
            base_vertices=bs[0]['base_vertices'],core_contacts=bs[0]['core_contacts'],
            dotted_squared_norm2_product=bs[0]['dotted_squared_norm2_product'],
            common_projection_norm=bs[0]['common_projection_norm'],
            load=str(load),initial_simplex_mass=str(old),slack=str(old-load)))
    violations=[c['index'] for c in columns if F(c['slack'])<0]
    assert len(occurrences)==len(columns)==241
    assert sum(F(c['load'])>0 for c in columns)==70 and len(violations)==53
    witness=next(c for c in columns if c['base_type']==54999)
    assert witness['positive_face_coefficients']==[[26259,1]]
    assert witness['simplex_face']==80143 and witness['A1_copies']==9
    assert F(witness['load'])==9*F(witness['initial_simplex_mass'])
    paths=[f'scripts/{TAG}.py']+[f'results/{tag}.json' for tag in INPUTS.values()]
    paths += ['results/fresh_sqrt5_flag_refined_base_capacity_20261017.md',
              'results/fresh_prism_simplex_shared_capacity_20261017.md']
    return dict(status='PASS_EXACT_CROSS_COFACE_A1_CAP_CAPACITIES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=dict(counts,columns=len(columns),positive_columns=70,violations=len(violations)),
        occurrences=occurrences,columns=columns,violating_columns=violations,
        necessary_simplex_mass_lower_bounds=[[fi,str(alpha.get(fi,F())),str(v)] for fi,v in sorted(required.items())],
        explicit_nine_copy_witness=witness,input_bound=src['input_bound'],
        new_repaired_bound=None,optimizer_invocations=0,
        scope='Valid additional signed capacities coupling identical isolated A1 cofaces of one maximal simplex. Selected source stars have trivial Lorentz-core automorphism, no other internal finite component after base deletion, and4<tau<6. All omitted types have zero positive coefficient. The unchanged old witness violates53columns; no repaired dual, primal feasibility, complete ambient model or dimension bound is asserted here.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],'old witness invalid; repair separate',flush=True)
