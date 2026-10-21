#!/usr/bin/env python3
"""Check saved contact-sum certificates with no producer imports."""
import json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import check_raw23_four_face_saved_20261018 as H

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw23_spacelike_contact_sum_20261018'
OLD='fresh_raw23_union_polar_completion_preflight_20261018'

def run():
    read=lambda tag:json.loads((ROOT/f'results/{tag}.json').read_text())
    data=read(TAG);old=read(OLD);H.bind(data);H.bind(old)
    assert data['status']=='PASS_RAW23_SPACELIKE_FACE_CONTACT_SUM'and len(data['cases'])==len(old['cases'])==3
    for case,o in zip(data['cases'],old['cases'],strict=True):
        assert case['kinds']==o['kinds'];rays=o['initial_cone']['rays'];xyz={r['index']:H.vec(r['coordinates'])for r in rays}
        zeros={r['index']:set(r['zero_facets'])for r in rays};contacts={r['index']:[-p for p in H.vec(r['products'])]for r in rays}
        rows=H.mat(o['initial_cone']['normal_rows']);metric=H.mat(o['metric']);sc=H.vec(o['original_squared_scales'])
        space=tuple(r['index']for r in rays if r['norm_sign']>0);common=tuple(sorted(set.intersection(*(zeros[i]for i in space))))
        assert len(space)==32 and common==(10,11)and case['spacelike_rays']==list(space)and case['common_active_walls']==list(common)
        assert tuple(i for i in xyz if set(common)<=zeros[i])==space
        assert all(p.sign()>=0 for ps in contacts.values()for p in ps)
        remaining=tuple(i for i in range(12)if i not in common)
        lower=[2 if(s-H.Q(4)).sign()>=0 else 1 for s in sc]
        assert all((s-H.Q(l*l)).sign()>=0 for s,l in zip(sc,lower))and lower==case['contact_lower_bounds']
        cache={}
        def rank(ids):
            if ids not in cache:cache[ids]=H.rank([xyz[i]for i in ids])
            return cache[ids]
        bs=tuple(case['full_face_basis']);assert len(bs)==len(set(bs))==8 and set(bs)<=set(space)and rank(bs)==8
        gram=lambda ids:[[H.bil(xyz[i],metric,xyz[j])for j in ids]for i in ids]
        assert H.mat(case['full_face_gram'])==gram(bs)and H.pd(gram(bs))
        domain=[s for size in range(11)for s in combinations(remaining,size)]
        assert len(domain)==len(case['active_decisions'])==1024;grouped={};fulls={}
        for subset,r in zip(domain,case['active_decisions'],strict=True):
            active=common+subset;full=tuple(i for i in space if set(active)<=zeros[i]);dim=rank(full)
            assert r['additional_active_walls']==list(subset)and r['active_rank']==H.rank([rows[i]for i in active])
            assert r['all_extreme_rays']==list(full)and r['dimension']==dim and r['origin_only']==(not full)
            if not full:continue
            closed=tuple(sorted(set.intersection(*(zeros[i]for i in full))))
            assert H.rank([rows[i]for i in closed])+dim==10 and r['closed_active_walls']==list(closed)
            grouped.setdefault(closed,[]).append(list(subset));fulls[closed]=full
        assert [f['active_walls']for f in case['faces']]==[list(k)for k in sorted(grouped)]
        dims=Counter();survive=[];excluded=0
        for f in case['faces']:
            active=tuple(f['active_walls']);full=fulls[active];bs=tuple(f['basis_rays']);dim=rank(full)
            assert f['all_extreme_rays']==list(full)and f['generating_subsets']==grouped[active]and f['dimension']==dim
            assert len(bs)==len(set(bs))==dim and set(bs)<=set(full)and rank(bs)==dim
            g=gram(bs);assert H.pd(g)and H.mat(f['plane_gram'])==g
            nonzero=[j for j in range(12)if any(contacts[i][j]!=0 for i in full)]
            assert nonzero==sorted(set(range(12))-set(active))and f['nonzero_wall_forms']==nonzero
            assert f['individual_contact_lower_bounds']==[lower[j]for j in nonzero]
            p=[sum((contacts[i][j]for j in nonzero),H.Q())for i in bs];z=H.vec(f['projection_coefficients'])
            assert H.vec(f['summed_contact_on_basis'])==p and len(z)==dim
            assert all(sum((g[i][j]*z[j]for j in range(dim)),H.Q())==p[i]for i in range(dim))
            q=sum((a*b for a,b in zip(p,z)),H.Q());L=sum(lower[j]for j in nonzero);margin=H.Q(L*L)-2*q
            assert q.sign()>0 and H.Q(*f['summed_form_squared_norm'])==q
            assert f['required_contact_sum']==L and H.Q(*f['Cauchy_margin'])==margin
            assert f['relative_interior_excluded']==(margin.sign()>0)
            dims[dim]+=1;excluded+=margin.sign()>0
            if margin.sign()<=0:survive.append(list(full))
        assert dict(dims)=={1:32,2:112,3:184,4:180,5:112,6:44,7:10,8:1}
        assert sorted(survive)==[[6],[31]]and excluded==673
        assert case['counts']==dict(active_subsets=1024,nonzero_faces=675,excluded_faces=673,
            excluded_ray_interiors=30,excluded_higher_face_interiors=643,unexcluded_ray_interiors=2)
        assert case['remaining_projected_ray_directions']==[6,31]
        print('CHECKED_CONTACT_SUM',case['kinds'],dict(dims),excluded,flush=True)
    for k,v in dict(active_subsets_checked=3072,classified_nonzero_faces=2025,excluded_higher_face_interiors=1929,
        excluded_ray_interiors=90,remaining_ray_directions=6,full_spacelike_face_dimension=8).items():assert data[k]==v
    assert data['arbitrary_passive_components_allowed']and data['all_nonzero_spacelike_ray_combinations_reduced_to_six_ray_directions']
    assert all(data[k]is False for k in('extension_ground_field_restricted','directions_with_timelike_generators_excluded',
        'arbitrary_completion_excluded','capacity_or_model_changed','zero_projection_excluded'))and data['new_dimension_bound']is None
    print('PASS_RAW23_SPACELIKE_CONTACT_SUM_SAVED',flush=True)

if __name__=='__main__':run()
