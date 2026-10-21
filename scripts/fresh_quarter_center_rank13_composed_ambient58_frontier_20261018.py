#!/usr/bin/env python3
"""Compose two local rules with the complete scale64 actual forest count.

The representative upper is deliberately NOT substituted for the exact
actual-type refinement. This gives a necessary inequality, not bound57.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
import fresh_a15_d15_quarter_a16_zero_center_retention_20261018 as Q
import fresh_rank13_variable_context_retention_20261018 as R
import fresh_quarter_center_rank13_ambient58_frontier_20261018 as C

ROOT=C.ROOT
TAG='fresh_quarter_center_rank13_composed_ambient58_frontier_20261018'
IDEAL=C.C.C.L.TAG
TAIL='fresh_a15_d15_dimension58_20261018'
PREVIOUS='fresh_d14_small_d15_composed_ambient58_frontier_20261018'
FINAL_REPORTS={Q.TAG:'e4542f034976b30547c07ee98896b5837264c17d38ee5532fcfadde86544cea1',
 R.TAG:'e139bf96e5d9d4210f3fc26b2ce98fa652ff4d75925ed4914f2b598eed4c4807',
 C.TAG:'746ad72af5625b28398cf81cfbdc0de04cd9132baf4d2241fec365a5bb381b77'}

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    assert all(h!='PENDING'for h in FINAL_REPORTS.values())
    tags=(Q.TAG,R.TAG,C.TAG,IDEAL,TAIL,PREVIOUS)
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in tags}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    for s,h in FINAL_REPORTS.items():assert digest(f'results/{s}.md')==h,s
    q,r,c,idealbase,tail,prior=(docs[s]for s in tags)
    assert q['baseline']==r['baseline']==R.G.TAG
    assert q['triangle_minimum']==r['triangles']=='2'
    assert q['quadrilateral_minimum']==r['quadrilaterals']=='1'
    assert not q['ideal_weights_changed']and not r['ideal_weights_changed']
    assert not q['contextual_rules_changed']and r['other_context_rule_unchanged']
    assert q['small_D15_center']['new']=='0'and r['chosen_parameter']=='9/64'
    epsilon=R.TARGET-R.OLD;assert epsilon==F(5,64)
    assert C.Q is Q and C.R is R and C.C.context is R.C.context
    assert 256*epsilon==20 and 64*4*R.TARGET==36
    variants=[]
    rv=r['later_finite_interfaces']['new_D14_variants']
    for a,b in zip(q['whole_variants'],rv,strict=True):
        assert(a['index'],a['extra_contact'])==(b['index'],b['extra_contact'])
        corners=[]
        for ac,bc in zip(a['corners'],b['corners'],strict=True):
            assert ac['marks']==bc['marks']and bc['coefficient']>=0
            after=F(ac['new'])+epsilon*bc['coefficient'];assert after>=0
            corners.append(dict(marks=ac['marks'],quarter_lower=ac['new'],
                parameter_coefficient=bc['coefficient'],joint_lower=str(after)))
        total=sum(F(x['joint_lower'])for x in corners);assert total>=1
        variants.append(dict(index=a['index'],extra_contact=a['extra_contact'],corners=corners,total=str(total)))
    assert len(variants)==24
    # No positive contribution in the four rank13 seals spends either
    # negative role of the quarter/zero-center change. Positive restorations
    # can be ignored in a lower bound.
    seals=[]
    for s in r['local_seals']:
        g=R.A.matrix(s['gram']);rows=[]
        for x in s['corners']:
            lower=F(x['parameter_values']['9/64'])
            if lower==0:continue
            marks=x['marks'];blocks=x['components']
            aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
            if aa is bb and aa['type']=='A16':
                assert (1,14)not in R.T.mark_pairs(g,aa['nodes'],marks,16,'A')
            if aa is not bb:
                for small,large,lm in((aa,bb,marks[1]),(bb,aa,marks[0])):
                    if 3<=len(small['nodes'])<=5 and large['type']=='D15':
                        assert (12,)not in R.T.mark_pairs(g,large['nodes'],[lm],15,'D')
            rows.append(dict(marks=marks,retained_lower=str(lower)))
        assert sum(F(x['retained_lower'])for x in rows)>=1
        seals.append(dict(name=s['name'],positive_contributions=rows))
    assert len(seals)==4
    for a in C.TYPES:
        assert C.actual_base(a)==2*Q.actual_base(a)
        for b in C.TYPES:
            assert F(C.actual_cross(a,b),64)==R.actual_cross(a,b)-F(Q.center_debit(a,b),32)
    ff,ref,ii=c['finite'],c['actual_typed_refinement'],c['ideal']
    assert ff['certified_upper']=='3573/8'and ref['exact_actual_maximum']=='761/2'
    assert ff['exact_actual_maximum']==ref['exact_actual_maximum']
    assert ff['actual_maximum_claimed']and not ff['actual_attainment_of_global_upper_claimed']
    classes=ref['classes'];assert len(classes)==ref['complete_expanded_upper_classes']==40
    expected={(tuple(x['counts']),x['fixed_A5_count'],x['fixed_A12_count'])for x in ff['above_curvature']}
    observed={(tuple(x['counts']),x['fixed_A5_count'],x['fixed_A12_count'])for x in classes}
    assert observed==expected and len(observed)==len(classes)
    assert ref['all_other_classes_at_most_curvature']and ref['both_context_guards_exact']
    assert sum(F(x['actual_maximum'])>377 for x in classes)==ref['classes_still_above']==12
    checked_contexts=0
    for row in classes:
        assert max(x['sixtyfourths_actual_maximum']for x in row['rows'])==64*F(row['actual_maximum'])
        assert row['actual_obstruction_retained']==(F(row['actual_maximum'])>377)
        for forest in row['maximizers']:
            fs=list(map(tuple,forest))
            assert sum(n for k,n in fs)==58 and C.value(fs)==64*F(row['actual_maximum'])
            assert C.context(fs)==64*R.context(fs);checked_contexts+=1
    for row in ff['expanded_alternatives']:
        fs=[(s[0],int(s[1:]))for s in row['forest']]
        assert C.context(fs,True)==64*R.context(fs,upper=True);checked_contexts+=1
    witness=[('A',16)]*2+[('D',26)]
    assert F(C.value(witness),64)==F(761,2)
    assert prior['finite_exact_maximum']=='763/2'
    assert ii['complete_profiles']==2552 and ii['all_positive']
    assert ii['weights_and_envelopes_unchanged']and ii['inherited_positive_low8_concentration_and_factor_fiber']
    assert ii['absolute_surplus']==idealbase['ideal']['absolute_surplus']=='36'
    assert ii['normalized_minimum']=='87/10240'
    d=57;g0=sum(min(j,d-j)for j in range(d+1))
    curvature=F(5*g0-d*(d+1),2)
    strict=5*sum(2*j-d for j in range(d+1)if 2*j>d)
    excess=F(ref['exact_actual_maximum'])-curvature
    assert(g0,curvature,strict,excess)==(812,F(377),4205,F(7,2))
    assert c['finite_curvature']==str(curvature)and c['strict_global_constant']==strict
    assert c['conditional_on_joint_quarter_center_and_rank13_whole_face_composition']
    lower=F(strict)/excess;ceil_lower=-(-lower.numerator//lower.denominator)
    assert ceil_lower==1202 and tail['dimension_bound']==58
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in tags for folder,ext in
        [('scripts','py'),('results','json'),('results','md')]}
    return dict(status='PASS_COMPOSED_QUARTER_CENTER_RANK13_AMBIENT58_FRONTIER',
        sources={p:digest(p)for p in sorted(paths)},hypotheses=tail['hypotheses'],
        ambient_dimension=58,completed_dimension_bound=58,new_dimension_bound=None,
        conditional_count_geometric_hypothesis_discharged=True,
        triangle_minimum=2,quadrilateral_minimum=1,quarter_variants=variants,rank13_seals=seals,
        finite_base_identities=len(C.TYPES),finite_cross_identities=len(C.TYPES)**2,
        context_rule_identical=True,checked_context_witnesses=checked_contexts,ideal_rule_unchanged=True,
        finite_exact_actual_maximum='761/2',representative_upper='3573/8',previous_exact_actual_maximum='763/2',
        finite_curvature=str(curvature),finite_excess=str(excess),
        finite_attaining_formal_forest=[list(t)for t in witness],
        above_curvature_upper_classes=40,above_curvature_actual_classes=12,
        every_obstructing_forest_enumerated=False,geometric_realization_asserted=False,
        positive_ideal_profiles=2552,actual_cusp_surplus='36',strict_global_constant=str(strict),
        necessary_ambient58_inequality='(7/2)*v_f >= 4205 + 36*v_i',
        integer_vertex_requirement='v_f >= ceil((8410 + 72*v_i)/7)',minimum_finite_vertices=ceil_lower,
        contradiction_obtained=False,DP_replayed=False,geometric_census_repeated=False,
        scope='Joint chosen-output local composition, complete actual-type refinement of every obstructing representative class, unchanged positive-cusp theorem and exact final arithmetic. The coarse representative upper is not used as the final actual bound. This gives a necessary vertex inequality, not a new dimension theorem or polytope.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite_exact_actual_maximum'],out['above_curvature_actual_classes'],out['minimum_finite_vertices'],flush=True)
