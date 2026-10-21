#!/usr/bin/env python3
"""Require an unmarked whole A1 for the newly restored rank15 quarters.

Complete selected-output retention plus a monotone refinement of all three
affected above-curvature classes. No global dynamic program is repeated.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement,product
import fresh_a16_d16_witness_context_ambient58_frontier_20261018 as C

ROOT=C.ROOT;L=C.L;W=C.W;Q=W.C.Q;R=W.C.R;T=R.T
TAG='fresh_rank15_a1_context_composed_ambient58_frontier_20261018'
FINAL_REPORTS={C.TAG:'37ead4131f4676475d44bb92fd621eba68a4859c9715383d37cbba2cfaac2123',
 L.TAG:'2701281f8d5d3660a7ff792c5829c59c8e30b5f168ddea1aa33d509def66271f',
 W.TAG:'c8acf530f290cd3ea1ebae1bcf88af06b82f17423952ca0cb0d3a0819ca0e004'}
STEMS=(C.TAG,L.TAG,W.TAG,Q.TAG,R.TAG)

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def rank15(t):return t in(('A',15),('D',15))
def actual_base(t):return C.actual_base(t)-(32 if rank15(t)else 0)
actual_cross=C.actual_cross
def context(fs):
    co=Counter(fs)
    return C.context(fs)+(32*(co['A',15]+co['D',15])if co['A',1]>0 else 0)
def value(fs):
    return sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)

def local(docs):
    variants=[];required=[]
    for p,d in zip(docs[Q.TAG]['whole_variants'],docs[L.TAG]['current_interfaces']['retained_D14_variants'],strict=True):
        assert(p['index'],p['extra_contact'])==(d['index'],d['extra_contact'])
        cs=[]
        for c,delta in zip(p['corners'],d['corners'],strict=True):
            assert c['marks']==delta['marks']and F(delta['delta'])>=0
            role=c['quarter_role'];assert role==Q.N.role(c)
            has_a1=any(b['type']=='A1'for b in c['components'])
            before=F(c['new'])
            if role in('A15','D15'):
                assert before==F(1,4)
                shared=next(b for b in c['components']if set(c['marks'])<=set(b['nodes']))
                assert shared['type']==role
                assert all(not set(c['marks'])&set(b['nodes'])for b in c['components']if b['type']=='A1')
            retained=F()if role in('A15','D15')and not has_a1 else before
            cs.append(dict(marks=c['marks'],role=role,has_unmarked_whole_A1=has_a1,
                quarter_lower=before,guarded_quarter_lower=retained,
                rank16_transfer_increment=delta['delta'],joint_lower=retained+F(delta['delta'])))
        total=sum(c['joint_lower']for c in cs);assert total>=1
        if any(c['quarter_role']=='A16'for c in p['corners']):
            paying=[c for c in cs if c['role']in('A15','D15')]
            assert len(paying)==2 and all(c['has_unmarked_whole_A1']for c in paying)
            required.append([p['index'],p['extra_contact']])
        variants.append(dict(index=p['index'],extra_contact=p['extra_contact'],corners=cs,total=total))
    assert len(variants)==24 and required==[[673,None],[1558,None],[1646,13],[1933,3]]
    # The new rank16 proof spends only its selected rank16 weight and A4/A12
    # coupons. No positive terminal output is a restored rank15 quarter.
    positive=[]
    for p in docs[L.TAG]['local_proofs']['whole_proofs']:
        for c in p['corners']:
            if F(c['payment'])<=0:continue
            ms=c['marks'];blocks=c['components']
            aa,bb=[next(b for b in blocks if m in b['nodes'])for m in ms]
            assert not(aa is bb and aa['type']in('A15','D15'))
            positive.append(dict(index=p['index'],marks=ms,payment=c['payment'],marked_types=[aa['type'],bb['type']]))
    assert len(positive)==12
    return dict(whole_variants=variants,required_guarded_variants=required,
        new_rank16_positive_outputs=positive,triangle_minimum=2,quadrilateral_minimum=1)

def refined_classes(docs):
    old=docs[C.TAG];ref=old['actual_typed_refinement'];ff=old['finite']
    assert ref['exact_actual_maximum']=='380'and ref['classes_still_above']==9
    above=[r for r in ref['classes']if r['actual_obstruction_retained']]
    affected=[r for r in above if r['counts'][6]>0]
    assert len(above)==9 and len(affected)==3 and all(r['residual_rank']==0 for r in affected)
    fibers={}
    for t in C.TYPES:fibers.setdefault(C.representative(t),[]).append(t)
    rows=[];assignments=0
    for r in affected:
        key=(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count'])
        u=next(z for z in ff['above_curvature']if(tuple(z['counts']),z['fixed_A5_count'],z['fixed_A12_count'])==key)
        rep=C.fixed_components(tuple(r['counts']),u['key'])+[('A',5)]*r['fixed_A5_count']
        groups=sorted(Counter(rep).items())
        expected=sorted(tuple(sorted(t for group in choice for t in group))for choice in product(*(combinations_with_replacement(fibers[t],n)for t,n in groups)))
        assert expected==sorted(tuple(map(tuple,z['fixed_actual']))for z in r['rows'])
        costs=[]
        for z in r['rows']:
            fs=list(map(tuple,z['fixed_actual']));assert sum(n for k,n in fs)==58
            assert z['actual_A1_count']==z['actual_A2_count']==0 and z['values']==[0]
            assert ('A',1)not in fs
            oldcost=C.value(fs);newcost=value(fs)
            assert oldcost==z['sixtyfourths_actual_maximum']
            assert newcost==oldcost-32*sum(rank15(t)for t in fs)
            costs.append(dict(forest=fs,previous_cost=F(oldcost,64),new_cost=F(newcost,64)));assignments+=1
        maximum=max(z['new_cost']for z in costs)
        rows.append(dict(counts=r['counts'],fixed_A5_count=r['fixed_A5_count'],fixed_A12_count=r['fixed_A12_count'],
            previous_maximum=r['actual_maximum'],new_maximum=maximum,complete_actual_assignments=costs,
            maximizers=[z['forest']for z in costs if z['new_cost']==maximum]))
    assert sorted(r['new_maximum']for r in rows)==[F(377),F(378),F(379)]
    unchanged=[r for r in above if r not in affected]
    for r in unchanged:
        assert r['counts'][6]==0
        assert any(value(list(map(tuple,fs)))==64*F(r['actual_maximum'])for fs in r['maximizers'])
    assert len(unchanged)==6
    maximum=max([F(r['actual_maximum'])for r in unchanged]+[r['new_maximum']for r in rows])
    witness=[('D',17)]*2+[('D',24)]
    assert maximum==F(1519,4)and value(witness)==C.value(witness)==64*maximum
    return dict(complete_predecessor_upper_classes=39,previous_above_curvature_classes=9,
        changed_complete_classes=rows,checked_fixed_assignments=assignments,
        unchanged_above_curvature_classes=6,new_above_curvature_classes=8,
        global_exact_maximum=maximum,global_attaining_forest=witness,
        pointwise_upper_from_complete_predecessor=True,global_DP_repeated=False,
        all_other_predecessor_classes_remain_at_most_curvature=True)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    for s,h in FINAL_REPORTS.items():assert digest(f'results/{s}.md')==h
    assert C.L is L and C.W is W and L.context is W.context
    assert docs[L.TAG]['baseline']==W.TAG
    for t in C.TYPES:
        assert C.actual_base(t)==L.actual_base(t)
        assert actual_base(t)>=0
        for u in C.TYPES:assert C.actual_cross(t,u)==L.actual_cross(t,u)
    geometry=local(docs);finite=refined_classes(docs)
    ii=docs[C.TAG]['ideal']
    assert ii['complete_profiles']==2552 and ii['all_positive']
    assert ii['absolute_surplus']=='36'and ii['normalized_minimum']=='87/10240'
    assert ii['weights_and_envelopes_unchanged']and ii['inherited_positive_low8_concentration_and_factor_fiber']
    d=57;g0=sum(min(j,d-j)for j in range(d+1));curvature=F(5*g0-d*(d+1),2)
    strict=5*sum(2*j-d for j in range(d+1)if 2*j>d)
    excess=F(finite['global_exact_maximum'])-curvature
    assert(g0,curvature,strict,excess)==(812,F(377),4205,F(11,4))
    v0=F(strict)/excess;minimum=-(-v0.numerator//v0.denominator);assert minimum==1530
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return W.enc(dict(status='PASS_COMPOSED_RANK15_A1_CONTEXT_AMBIENT58_FRONTIER',
        sources={p:digest(p)for p in sorted(paths)},baseline=L.TAG,
        local_proof=geometry,finite=finite,
        new_rule='Restore A15(0,13)/(1,14) and D15(1,13)/(1,14) quarters only with a whole unmarked A1 component',
        ordinary_A15_D15_base='100',scale64_context='32*(count(A15)+count(D15)) if count(A1)>0; otherwise0',
        contextual_witness_multiplicity_not_counted=True,pointwise_cost_nonincreasing=True,
        finite_base_identities=len(C.TYPES),finite_cross_identities=len(C.TYPES)**2,
        conditional_predecessor_geometry_discharged=True,ideal_weights_changed=False,
        positive_ideal_profiles=2552,actual_cusp_surplus='36',finite_curvature=str(curvature),
        finite_excess=str(excess),strict_global_constant=str(strict),
        necessary_ambient58_inequality='(11/4)*v_f >= 4205 + 36*v_i',
        integer_vertex_requirement='v_f >= ceil((16820+144*v_i)/11)',minimum_finite_vertices=minimum,
        completed_dimension_bound=58,new_dimension_bound=None,geometric_example_claimed=False,
        scope='Whole-A1 witnessed rank15 repair under complete current chosen-output replay, simultaneous A16/D16 composition, and monotone exact refinement of all affected above-curvature actual classes. The frozen complete count supplies all remaining upper classes and ideal polynomial checks. No global DP, convolution or new dimension theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['global_exact_maximum'],out['finite']['new_above_curvature_classes'],out['minimum_finite_vertices'],flush=True)
