#!/usr/bin/env python3
"""Saved whole-corner and witness diagnostic; no new face theorem or DP."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import fresh_d14_small_d15_ambient58_frontier_20261018 as C

ROOT=C.ROOT;enc=C.enc
TAG='fresh_a16_restoration_alternative_cost_diagnostic_20261018'
S=C.D14
A='fresh_finite_a15_three_coupon_zero_transfer_20261018'
D='fresh_finite_d15_a4_a11_half_zero_transfer_20261018'
STEMS=(C.TAG,S,A,D,C.P.TAG)
DA={(1,13),(1,14)};AA={(0,13),(1,14)};A16={(1,14)}

def role(c):
    if not c['images']:return None
    whole=[g for g in c['components']if all(m in g['nodes']for m in c['marks'])]
    assert len(whole)==1
    images={tuple(p)for p in c['images']};t=whole[0]['type']
    if t=='D15'and images==DA:return 'D15'
    if t=='A15'and images==AA:return 'A15'
    if t=='A16'and images==A16:return 'A16'

def delta(fs,x,y,z,center):
    co=Counter(fs)
    small=sum(n*k for (t,n),k in co.items()if 3<=n<=5)
    return 2*x*co['D',15]+2*y*co['A',15]+(z-F(1,2))*co['A',16]-(F(3,8)-center)*small*co['D',15]

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    for stem,orbit in((A,AA),(D,DA)):
        w={tuple(r['pair']):F(r['weight'])for r in docs[stem]['weighted_pairs']}
        assert all(w[p]==0 for p in orbit)and sum(w.values())==100
    variants=[];required=[]
    for p in docs[S]['selected_proofs']:
        for v in p.get('variants',[]):
            corners=[]
            for c in v['corners']:
                rr=role(c)
                if rr in('A15','D15'):assert c['payment']=='0'
                if rr=='A16':assert c['payment']=='1/2'
                corners.append(dict(marks=c['marks'],images=c['images'],components=[g['type']for g in c['components']],
                    current_payment=c['payment'],alternative_role=rr))
            r=dict(index=p['index'],extra_contact=v['extra_contact'],corners=corners,total=v['total'])
            variants.append(r)
            if any(c['alternative_role']=='A16'for c in corners):
                assert Counter(c['alternative_role']for c in corners)==Counter({None:1,'D15':1,'A15':1,'A16':1})
                assert sum(F(c['current_payment'])for c in corners if c['alternative_role']is None)==F(1,2)
                required.append(r)
    assert len(variants)==24
    assert [(r['index'],r['extra_contact'])for r in required]==[(673,None),(1558,None),(1646,13),(1933,3)]
    classes=docs[C.TAG]['actual_typed_refinement']['classes'];witnesses=[];classmax=[]
    for ci,c in enumerate(classes):
        for ri,r in enumerate(c['rows']):
            fs=tuple(map(tuple,r['attaining_forest']));v=F(C.value(fs),32)
            assert sum(n for k,n in fs)==58 and v==F(r['thirtyseconds_actual_maximum'],32)
            out=dict(class_index=ci,row_index=ri,forest=fs,current_cost=v,
                class_current_maximum=c['actual_maximum'],current_class_above=F(c['actual_maximum'])>377)
            witnesses.append(out)
            if v==F(c['actual_maximum'])and out['current_class_above']:classmax.append(len(witnesses)-1)
    assert len(classes)==20 and len(witnesses)==443 and sum(F(c['actual_maximum'])>377 for c in classes)==16
    scenarios=[]
    for center in(F(0),F(1,8),F(1,4),F(3,8)):
        for ix in range(9):
            for iy in range(9-ix):
                x,y,z=F(ix,16),F(iy,16),F(8-ix-iy,16)
                face_values=[]
                for v in variants:
                    total=sum({'D15':x,'A15':y,'A16':z}.get(c['alternative_role'],F(c['current_payment']))for c in v['corners'])
                    assert total>=1
                    face_values.append(total)
                values=[r['current_cost']+delta(r['forest'],x,y,z,center)for r in witnesses]
                byclass={ci:max(v for v,r in zip(values,witnesses)if r['class_index']==ci)for ci in range(len(classes))}
                mx=max(values)
                scenarios.append(dict(D15_pair_weight=x,A15_pair_weight=y,A16_singleton_weight=z,small_D15_center_weight=center,
                    actual_base_deltas={'D15':2*x,'A15':2*y,'A16':z-F(1,2)},
                    additional_center_decrease=F(3,8)-center,
                    minimum_saved_whole_variant_payment=min(face_values),
                    maximum_over_saved_witnesses=mx,
                    attaining_saved_witness_indices=[i for i,v in enumerate(values)if v==mx],
                    saved_classes_with_known_obstruction=sum(v>377 for v in byclass.values()),
                    changed_saved_class_maxima=byclass,
                    original_16_class_maximizer_values=[values[i]for i in classmax]))
    best=min(r['maximum_over_saved_witnesses']for r in scenarios)
    bests=[r for r in scenarios if r['maximum_over_saved_witnesses']==best]
    blocker=[('A',1)]+[('A',2)]*9+[('A',13)]*3
    assert F(C.value(blocker),32)==F(6093,16)
    assert all(delta(blocker,r['D15_pair_weight'],r['A15_pair_weight'],r['A16_singleton_weight'],r['small_D15_center_weight'])==0 for r in scenarios)
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_SAVED_A16_RESTORATION_ALTERNATIVE_DIAGNOSTIC',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
        baseline=C.TAG,whole_variants_checked=24,required_four_variants=required,
        complete_required_role_identity='D14half + D15(1,13/14) x + A15(0,13)/(1,14) y + A16(1,14) z; x+y+z>=1/2',
        all_saved_variants=variants,witnesses=witnesses,original_above_classes=16,
        complete_saved_actual_fiber_witnesses=443,grid_step='1/16',center_probes=['0','1/8','1/4','3/8'],
        scenarios=scenarios,best_saved_witness_maximum=best,best_grid_scenarios=bests,
        invariant_obstructing_forest=blocker,invariant_obstruction_cost='6093/16',
        saved_witness_lower_bound_only=True,all_unrecorded_forests_reoptimized=False,
        universal_retention_of_decreased_A16_weight_proved=False,
        further_center_zero_theorem_proved=False,new_DP_performed=False,new_geometry_enumerated=False,
        completed_unrestricted_bound=58,new_dimension_bound=None,
        scope='Exact canonical-role extraction from24 frozen whole variants and cost shifts of443 saved actual-fiber witnesses. Every proposal passes these saved whole variants; arbitrary-support local theorems and global output retention are separate. Values are lower obstructions, not reoptimized upper bounds or complete actual forest enumeration.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['best_saved_witness_maximum'],len(out['best_grid_scenarios']),flush=True)
