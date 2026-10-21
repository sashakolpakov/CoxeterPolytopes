#!/usr/bin/env python3
"""Order12 B64/C63 with an A2/D11 short-leaf coupon upgrade.

Complete inherited selected support is reconstructed. The arbitrary-support
sealing proof is inherited and reviewed explicitly in the written report.
No finite-forest maximum or ambient theorem is asserted here.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_order12_affine_b_half_transfer_20261017 as B
import fresh_order12_affine_c_endpoint_reduction_20261017 as C

ROOT=B.ROOT
TAG='fresh_order12_bc_unit_coupon_transfer_20261017'
CDROP={(0,10),(1,11)}
STEMS=(B.TAG,C.TAG,B.DIAG,B.OVER,B.D13,
       'fresh_affine_d11_leaf_pair_transfer_20261017',
       'fresh_cycle12_ideal5_transfer_20261017',
       'fresh_finite_d17_long_pair_support_20261017',
       'fresh_ideal5_weighted_dimension66_20261017',
       'fresh_order14_affine_b_fork_projection_transfer_20261017',
       'fresh_large_finite_d_leaf_half_transfer_20261017')


def improved_corners(payment):
    """Recount actual finite components, with no ideal B contribution."""
    rows=[]
    for c in payment['corners']:
        if c['kind']!='finite':
            assert c['kind']=='affineB11_order12' and F(c['ideal_weight'])==F(1,2)
            rows.append(dict(marks=c['marks'],kind=c['kind'],new_weight='0'))
            continue
        marks=c['marks'];cc=c['components']
        a,b=[next(r for r in cc if m in r['nodes']) for m in marks]
        assert c['signature'][1:]==[0,0]
        unit=int(a is b);assert unit==c['finite_unit']
        assert not unit or len(a['nodes'])<=13
        coupon=F(c['finite_half']);upgrade=None
        for small,large,lm in ((a,b,marks[1]),(b,a,marks[0])):
            if small['type'] in ('A1','A2') and large['type']=='D11' and large['roles'].get(str(lm))=='short_fork_leaf':
                assert a is not b and coupon==F(1,2)
                coupon=F(1)
                upgrade=dict(small_type=small['type'],small_nodes=small['nodes'],
                    D11_nodes=large['nodes'],short_leaf_mark=lm,
                    source='existing_A1_upgrade' if small['type']=='A1' else 'new_A2_upgrade')
        rows.append(dict(marks=marks,kind='finite',components=cc,
            finite_unit=unit,old_half=c['finite_half'],new_coupon=str(coupon),
            upgrade=upgrade,new_weight=str(unit+coupon)))
    total=sum(F(c['new_weight']) for c in rows)
    assert total>=1
    return dict(corners=rows,new_finite_total=str(total),new_ideal_total='0')


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    over=data[B.OVER];support=data[B.DIAG];patterns=support['high_patterns']
    assert len(patterns)==758
    oldb=data[B.TAG];oldc=data[C.TAG]
    brec=next(r for r in support['records'] if r['kind']=='B')
    crec=next(r for r in support['records'] if r['kind']=='C')
    # Re-enumerate all eligible patterns, including every failed pattern.
    all_support={};affected={}
    for kind,rec,drop,mapper in [('B',brec,B.DROP,B.P.embeddings),('C',crec,CDROP,C.path_maps)]:
        saved={r['case']:r for r in rec['high_case_embeddings']};rows=[];hits=[]
        for p in patterns:
            maps=mapper(p['labels'],rec['target'])
            pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
            if maps:
                assert sorted(maps)==sorted(saved[p['index']]['embeddings'])
                assert pairs==sorted(map(tuple,saved[p['index']]['marked_images']))
                rows.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
                if set(pairs)&drop:hits.append(rows[-1])
            else:assert p['index'] not in saved
        assert len(rows)==len(saved)
        all_support[kind]=rows;affected[kind]=hits
    assert [r['index'] for r in affected['B']]==B.EXPECTED
    assert [r['index'] for r in affected['C']]==[164,1452]
    brows=[];reconstructed={};outcomes=Counter();totals=Counter();upgrades=Counter()
    for saved in affected['B']:
        c=over['isolation']['cases'][saved['index']];core=over['isolation']['cores'][c['core_key']]
        row=B.inspect_case(c,core)
        inherited=next(r for r in oldb['cases'] if r['index']==c['index'])
        assert all(inherited[k]==v for k,v in row.items())
        reconstructed[c['index']]=row;variants=[]
        for key in ('sealed_payment','no_extra_root_payment','one_extra_root_payment'):
            if key in row:
                value=improved_corners(row[key]);value['variant']=key
                variants.append(value);totals[value['new_finite_total']]+=1
                for corner in value['corners']:
                    if corner.get('upgrade'):upgrades[corner['upgrade']['source']]+=1
        if not variants:
            assert row['outcome']=='whole_finite_B11_before_any_affine_completion'
            assert row['finite_unit']==1
        outcomes[row['outcome']]+=1
        brows.append(dict(index=c['index'],embeddings=saved['embeddings'],
            inherited_geometry=row,payment_variants=variants))
    assert outcomes==Counter(sealed_radical_core=8,
        whole_finite_B11_before_any_affine_completion=2,zero_or_one_sealed_unit_contact_root=4)
    assert sum(totals.values())==16 and min(map(F,totals))>=1
    crows=[]
    for saved in affected['C']:
        c=over['isolation']['cases'][saved['index']];core=over['isolation']['cores'][c['core_key']]
        checked=reconstructed[c['index']]
        assert checked['outcome']=='whole_finite_B11_before_any_affine_completion'
        ids=c['marked_component'];inv=[[F(x) for x in row] for row in core['inverse']]
        scale=list(map(F,core['gauge']));target=crec['target'];maps=[]
        for mp in saved['embeddings']:
            image={t:ids[i] for i,t in enumerate(mp)};missing=set(range(12))-set(image)
            assert len(missing)==1
            entries=[]
            for new in missing:
                for old,node in image.items():
                    if target[new][old]:
                        label=target[new][old];assert label==4
                        bound=2*scale[node]*inv[node][node]
                        assert node in c['preferred_corner'] and node in c['possible_contacts']
                        assert bound==16
                        entries.append(dict(new_target_node=new,old_target_node=old,
                            old_core_node=node,native_label=label,norm_two_contact_square='2',
                            old_gauge_square=str(scale[node]),inverse_diagonal=str(inv[node][node]),
                            forced_projection_lower_bound=str(bound)))
            assert len(entries)==1
            maps.append(dict(map=mp,forced_entries=entries))
        crows.append(dict(index=c['index'],maps=maps,
            whole_finite_component=checked['whole_finite_component'],retained_finite_unit=1,
            end_projection_certificate=checked))
    assert sum(len(r['maps']) for r in crows)==4
    bpaid=set(combinations(range(12),2))-B.DROP
    cpaid=set(map(tuple,oldc['new_paid_pairs']))-CDROP
    assert len(bpaid)==64 and len(cpaid)==63
    assert bpaid==set(map(tuple,oldb['basic_pairs_retained_at_weight_one']))
    assert cpaid==set(map(tuple,oldc['basic_pairs']))
    for paid,rec,drop in ((bpaid,brec,B.DROP),(cpaid,crec,CDROP)):
        assert set(map(tuple,rec['short8_pairs']))<=paid
        assert all({tuple(sorted(a[i] for i in p)) for p in drop}==drop
                   for a in B.P.embeddings(rec['target'],rec['target']))
        assert all(not(drop&set(map(tuple,r['marked_images']))) for r in rec['low_case_embeddings'])
    assert len(support['low_corner_domain']['cases'])==66
    assert sum(len(r['variants']) for r in support['low_corner_domain']['cases'])==69
    assert not crec['immediate_pairs']
    gaps=B.immediate_gap(data[B.D13]);assert gaps==oldb['immediate_D13_strict_gap']
    costs=dict(data['fresh_ideal5_weighted_dimension66_20261017']['small_order_costs']['12'])
    assert costs==dict(A=60,B=65,C=65,D=62)
    costs.update(B=64,C=63);assert max(costs.values())==64
    paths=[Path(__file__),Path(B.__file__),Path(C.__file__),Path(B.L.__file__),Path(B.Q.__file__),Path(B.P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER12_B64_C63_WITH_A2_D11_UNIT_UPGRADE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        complete_eligible_high_patterns=758,reconstructed_high_support=all_support,
        B=dict(old_cost=65,new_cost=64,removed_half_pairs=sorted(B.DROP),
            retained_unit_pairs=sorted(bpaid),affected_maps=sum(len(r['embeddings']) for r in affected['B']),
            cases=brows,outcomes=dict(outcomes),new_finite_total_histogram=dict(totals)),
        C=dict(old_cost=65,new_cost=63,removed_unit_pairs=sorted(CDROP),
            inherited_endpoint_pair_stays_absent=[0,11],retained_unit_pairs=sorted(cpaid),
            affected_maps=4,cases=crows),
        new_order12_costs=costs,new_rank11_envelope='16/3',
        coupon_upgrade='Whole A2 against a short fork leaf of a separate whole finite D11: half becomes unit; every A2 mark is eligible.',
        new_coupon_increment_per_A2_D11_component_pair='2',
        doubled_finite_budget_increment='4*a2*d11',used_upgrades=dict(upgrades),
        existing_A1_D11_full_upgrade_required=True,
        low_actual_cases=66,low_actual_variants=69,immediate_D13_strict_gap=gaps,
        large_D_D10_coupon_used=False,
        coupling='Independent of the large-D/D10 transfer geometrically; its new A2/D11 increment must be included jointly in any finite-forest budget.',
        triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
        scope='Complete selected B/C support and exact native projection/whole finite corner reconstruction, including both zero-residual variants. Inherited all-label core completeness and arbitrary-support radical/residual sealing are used, not rerun as a new geometric census. No global finite maximum or ambient theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_order12_costs'],out['B']['new_finite_total_histogram'],out['used_upgrades'],flush=True)
