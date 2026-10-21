#!/usr/bin/env python3
"""Independent saved final-ray bounds, using standalone Q(sqrt5) arithmetic."""
import json
from pathlib import Path
import check_raw23_four_face_saved_20261018 as H

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw23_spacelike_pole_classification_20261018'
SUM='fresh_raw23_spacelike_contact_sum_20261018'
OLD='fresh_raw23_union_polar_completion_preflight_20261018'

def run():
    read=lambda tag:json.loads((ROOT/f'results/{tag}.json').read_text())
    data,summed,old=read(TAG),read(SUM),read(OLD)
    for d in(data,summed,old):H.bind(d)
    assert data['status']=='PASS_RAW23_NONZERO_SPACELIKE_FACE_POLE_CLASSIFICATION'
    assert len(data['cases'])==len(summed['cases'])==len(old['cases'])==3
    strict=forced=0
    for c,s,o in zip(data['cases'],summed['cases'],old['cases'],strict=True):
        assert c['kinds']==s['kinds']==o['kinds']and s['remaining_projected_ray_directions']==[6,31]
        assert [r['ray']for r in c['ray_tests']]==[6,31]
        by={r['index']:r for r in o['initial_cone']['rays']};metric=H.mat(o['metric']);sc=H.vec(o['original_squared_scales']);kept=[]
        for r in c['ray_tests']:
            rr=by[r['ray']];x=H.vec(rr['coordinates']);p=[-a for a in H.vec(rr['products'])];norm=H.bil(x,metric,x)
            assert norm.sign()>0 and all(a.sign()>=0 for a in p)
            assert H.vec(r['coordinates'])==x and H.Q(*r['squared_norm'])==norm and H.vec(r['contact_magnitudes'])==p
            bound=[2*a*a/norm/s for a,s in zip(p,sc)];assert H.vec(r['physical_contact_square_upper'])==bound
            below=[j for j,a in enumerate(p)if a!=0 and(bound[j]-1).sign()<0]
            assert r['strict_excluding_walls']==below
            if below:
                strict+=1;assert r['status']=='strict_contact_exclusion'and r['forced_pole']is None;continue
            forced+=1;kept.append(r['ray']);a=r['anchor_wall']
            assert p[a]!=0 and bound[a]==sc[a]==1
            scale=H.Q(1)/p[a];pole=[scale*z for z in x];squares=[scale*scale*z*z/t for z,t in zip(p,sc)]
            assert H.bil(pole,metric,pole)==2
            labels=[2 if z==0 else 3 if q==1 else 5 for z,q in zip(p,squares)]
            assert all(q in(H.Q(),H.Q(1),H.Q('3/2','1/2'))for q in squares)
            pp=r['forced_pole'];assert H.Q(*pp['scale'])==scale and H.vec(pp['coordinates'])==pole
            assert H.vec(pp['physical_contact_squares'])==squares and pp['labels']==labels
            assert labels.count(3)==2 and labels.count(5)==1
            assert r['status']=='unique_pure_pole'and r['anchor_physical_upper']=='1'
            assert r['forced_projection_squared_norm']=='2'and r['forced_passive_squared_norm']=='0'
        expect=[]if o['kinds']==[1,6]else[31]if o['kinds']==[1,1]else[6]
        assert kept==expect==c['compatible_poles_in_nonzero_spacelike_face']
    assert strict==data['additional_strict_ray_exclusions']==4 and forced==data['forced_original_poles']==2
    assert data['nonzero_spacelike_generated_domain_complete']and data['forced_passive_part_zero']
    assert all(data[k]is False for k in('extension_field_restricted','zero_projection_excluded',
        'timelike_generator_directions_excluded','arbitrary_completion_excluded','model_or_capacity_changed'))
    assert data['new_dimension_bound']is None
    print('PASS_RAW23_SPACELIKE_POLES_SAVED',strict,forced)

if __name__=='__main__':run()
