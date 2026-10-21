#!/usr/bin/env python3
"""Whole-face retention with two quarter orbits and unused D15 center zero.

No new native search or forest DP. Arbitrary-support proof is in the report,
using the bound complete selected/whole-face and chosen-output theorems.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
import fresh_d14_small_d15_ambient58_frontier_20261018 as C
import fresh_d14_a16_d15_joint_face_composition_20261018 as G
import fresh_a16_restoration_alternative_cost_diagnostic_20261018 as N

ROOT=C.ROOT
TAG='fresh_a15_d15_quarter_a16_zero_center_retention_20261018'
STEMS=(C.TAG,G.TAG,N.TAG,N.S,N.A,N.D,C.P.TAG,C.P.I.TAG,G.S.AW)
context=C.context

def actual_base(t):
    return C.actual_base(t)+(16 if t in(('A',15),('D',15))else 0)-(16 if t==('A',16)else 0)
def upper_base(t):
    return C.upper_base(t)+(16 if t==('D',15)else 0)-(16 if t==('D',16)else 0)
def center_debit(a,b):
    return sum(12*s[1]for s,t in((a,b),(b,a))if 3<=s[1]<=5 and t==('D',15))
def actual_cross(a,b):return C.actual_cross(a,b)-center_debit(a,b)
def upper_cross(a,b):return C.upper_cross(a,b)-center_debit(a,b)
def value(fs,upper=False):
    base,cross=(upper_base,upper_cross)if upper else(actual_base,actual_cross)
    return sum(base(t)for t in fs)+sum(cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs,upper)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    joint,selected,center,inventory=(docs[s]for s in(G.TAG,N.S,C.P.TAG,C.P.I.TAG))
    assert joint['baseline']==selected['baseline']==center['baseline']=='fresh_d12_half_face_composition_20261018'
    assert joint['triangle_minimum']=='2'and joint['quadrilateral_minimum']=='1'
    assert joint['selected_contexts']==17 and joint['selected_maps']==36 and joint['whole_variants']==24
    assert joint['proper_completion_exclusions']==4
    assert center['changed_rule']['new']=='3/8'
    assert center['output_routes']==inventory['universal_output_inventory']
    assert inventory['literal_inventory']['small_rank3_to5_D15_corners']==0
    assert inventory['later_finite_inventory']['D15_cross_corners']==0
    # Complete current same-component tables, with automorphism-invariant orbits.
    tables={}
    for stem,typ,orbit in((N.A,'A15',N.AA),(N.D,'D15',N.DA)):
        w={tuple(r['pair']):F(r['weight'])for r in docs[stem]['weighted_pairs']}
        assert sum(w.values())==F(C.actual_base((typ[0],15)),32)==100
        assert all(w[p]==0 for p in orbit)
        for p in orbit:w[p]=F(1,4)
        assert sum(w.values())==F(actual_base((typ[0],15)),32)==F(201,2)
        if typ=='A15':assert {tuple(sorted((14-i,14-j)))for i,j in orbit}==orbit
        else:assert {tuple(sorted((27-i if i in(13,14)else i,27-j if j in(13,14)else j)))for i,j in orbit}==orbit
        tables[typ]=[dict(pair=p,weight=w[p])for p in sorted(w)]
    aw={tuple(r['pair']):F(r['weight'])for r in selected['A16_weighted_pairs']}
    olda={tuple(r['pair']):F(r['weight'])for r in docs[G.S.AW]['A16_weighted_pairs']}
    assert aw[1,14]==F(1,2)and olda[1,14]==0
    aw[1,14]=0;assert aw==olda
    assert sum(aw.values())==F(actual_base(('A',16)),32)==F(441,4)
    tables['A16']=[dict(pair=p,weight=aw[p])for p in sorted(aw)]
    assert {tuple(sorted((15-i,15-j)))for i,j in {(1,14)}}=={(1,14)}
    variants=[];replaced=[]
    for p in selected['selected_proofs']:
        for v in p.get('variants',[]):
            corners=[]
            for c in v['corners']:
                role=N.role(c);old=F(c['payment'])
                new={'D15':F(1,4),'A15':F(1,4),'A16':F()}.get(role,old)
                if role in('D15','A15'):assert old==0
                if role=='A16':assert old==F(1,2)
                corners.append(dict(marks=c['marks'],components=c['components'],images=c['images'],
                    previous=old,new=new,quarter_role=role))
            assert len({tuple(r['marks'])for r in corners})==4
            total=sum(r['new']for r in corners);assert total>=1
            row=dict(index=p['index'],extra_contact=v['extra_contact'],corners=corners,total=total)
            variants.append(row)
            if any(r['quarter_role']=='A16'for r in corners):
                assert sorted(r['new']for r in corners)==[F(),F(1,4),F(1,4),F(1,2)]
                replaced.append((p['index'],v['extra_contact']))
    assert len(variants)==24 and replaced==[(673,None),(1558,None),(1646,13),(1933,3)]
    # The historical lower table deliberately did not credit restored A16.
    literal=selected['literal_interfaces']
    assert literal['credit_for_restored_A16_counted']is False
    assert len(literal['records'])==216 and all(F(r['total'])>=1 for r in literal['records'])
    assert sum(len(r['corners'])for r in selected['later_finite_interfaces'])==564
    assert all(c['new_debit']=='0'for r in selected['later_finite_interfaces']for c in r['corners'])
    assert len(selected['affine_D12_two_A1_interfaces'])==16
    # The center is intrinsic; leaf halves and all smaller-rank roles survive.
    marked=[]
    for r in(3,4,5):
        ws=[(i,j,F()if j==12 else F(1,2))for i in range(r)for j in(12,13,14)]
        assert sum(w for i,j,w in ws)==r
        assert sorted((i,27-j if j in(13,14)else j,w)for i,j,w in ws)==sorted(ws)
        marked.append(dict(rank=r,weights=ws,new_component_cost=r,saving_from_three_eighth=F(3*r,8)))
    for a in C.TYPES:
        assert actual_base(a)>=0 and upper_base(C.representative(a))>=actual_base(a)
        for b in C.TYPES:
            assert 0<=actual_cross(a,b)==actual_cross(b,a)<=upper_cross(C.representative(a),C.representative(b))
    assert actual_cross(('A',1),('D',15))==50
    assert all(actual_cross(t,('D',15))==96 for t in C.TYPES if t[1]==2)
    assert all(actual_cross(t,('D',15))==32*t[1]for t in C.TYPES if 3<=t[1]<=5)
    assert upper_base(('D',15))==32*F(201,2)and upper_base(('D',16))==32*F(441,4)
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return N.enc(dict(status='PASS_A15_D15_QUARTER_A16_ZERO_CENTER_FACE_RETENTION',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
        baseline=G.TAG,scale=32,whole_variants=variants,replaced_A16_variants=replaced,
        A15_restored_pairs=sorted(N.AA),D15_restored_pairs=sorted(N.DA),restored_weight='1/4',
        A16_lowered_pair=[1,14],A16_new_weight='0',weighted_tables=tables,
        new_actual_costs={'A15':'201/2','D15':'201/2','A16':'441/4','D14':'179/2','D16':'439/4'},
        new_representative_costs={'rank15':'201/2','rank16':'441/4'},
        small_D15_center=dict(ranks=[3,4,5],old='3/8',new='0',marked_rules=marked),
        A1_D15_total='25/16',rank2_D15_total='3',D15_short_leaf_weights='1/2',A1_D15_node11='1/16',
        unchanged_historical_lower_totals=216,restored_A16_credit_in_historical_table=False,
        retained_later_terminal_corners=564,retained_affine_D12_interfaces=16,
        complete_selected_patterns=600,complete_selected_maps=1104,affected_selected_maps=36,
        proper_completion_exclusions_retained=4,whole_face_proof_output_replay=True,
        simultaneous_type_pairs=len(C.TYPES)**2,ideal_weights_changed=False,contextual_rules_changed=False,
        new_cross_coupons=[],triangle_minimum='2',quadrilateral_minimum='1',
        forest_DP_performed=False,new_maximum_claimed=False,completed_unrestricted_bound=58,new_dimension_bound=None,
        scope='Complete current face-proof retention: positive A15/D15 restorations replace the only newly introduced A16 obligation, all original baseline outputs survive, and the inherited explicit chosen-output inventory permits zero at the unused rank3–5/D15 center. Saved checks supplement arbitrary-support theorems; no new census, forest optimum or dimension result.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_actual_costs'],flush=True)
