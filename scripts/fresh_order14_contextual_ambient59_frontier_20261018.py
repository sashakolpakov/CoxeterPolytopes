#!/usr/bin/env python3
"""Compose the proved order14 envelope and contextual finite certificate.

No optimization, historical census, or convolution is rerun. The complete
terminal domain and every saved polynomial/margin identity are checked.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order14_contextual_ambient59_frontier_20261018'
FRONT='fresh_ambient59_payment_frontier_20261018'
FINITE='fresh_rank13_contextual_coupon_rebalance_20261018'
B='fresh_order14_b_basic_distance9_transfer_20261018'
C='fresh_order14_c_basic_distance9_transfer_20261018'
D='fresh_order14_d_near_fork_half_transfer_20261018'
OLD='fresh_depth6_dimension59_20261018'
STEMS=(FRONT,FINITE,B,C,D,OLD)
FROZEN_D_MD='26dc0b458e6d84b7553572c0e051f172f63279bc315781fdb80c936b0b313e61'
FROZEN_FINITE_MD='2b5d7ac343174bc4e6df4ad1c74c7c0edbb03a92c322867d38e7a5cc3e4b6de7'


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    return x


def high_parts(n,low=9):
    if n==0:yield ();return
    for r in range(low,n+1):
        for tail in high_parts(n-r,r):yield (r,)+tail


def summary(rows):
    order=sorted(rows,key=lambda r:(r['normalized_margin'],r['ranks']))
    minimum=order[0]['normalized_margin'];absolute=min(r['absolute_margin'] for r in rows)
    return dict(count=len(rows),negative_count=sum(r['normalized_margin']<0 for r in rows),
        zero_count=sum(r['normalized_margin']==0 for r in rows),normalized_minimum=minimum,
        normalized_minimizers=[r['ranks'] for r in rows if r['normalized_margin']==minimum],
        absolute_minimum=absolute,absolute_minimizers=[r['ranks'] for r in rows if r['absolute_margin']==absolute],
        negative_rows=[dict(ranks=r['ranks'],normalized_margin=r['normalized_margin'],absolute_margin=r['absolute_margin']) for r in order if r['normalized_margin']<0])


def conditional(rows,savings):
    updated=[]
    for r in rows:
        m=r['normalized_margin']+sum(F(r['ranks'].count(k),k+1)*v for k,v in savings.items())
        updated.append(dict(ranks=r['ranks'],link_vertices=r['link_vertices'],normalized_margin=m,absolute_margin=m*r['link_vertices']))
    return dict(total_cost_savings_by_affine_rank=savings,summary=summary(updated),
        geometry_proved=False,additional_finite_coupons_counted=False,
        terminal_domain_only=True,low9_concentration_requires_revalidation=any(k<=9 for k in savings))


def threshold(rows,rank,fixed=None):
    fixed=fixed or {};values=[];unaffected=[]
    for r in rows:
        margin=r['normalized_margin']+sum(F(r['ranks'].count(k),k+1)*v for k,v in fixed.items())
        n=r['ranks'].count(rank)
        if n:values.append((-margin*F(rank+1,n),r['ranks']))
        elif margin<=0:unaffected.append(dict(ranks=r['ranks'],margin=margin))
    bound=max(x for x,_ in values)
    return dict(varied_affine_rank=rank,affine_order=rank+1,fixed_total_cost_savings=fixed,
        strict_total_cost_saving_threshold=bound,critical_profiles=[p for x,p in values if x==bound],
        unaffected_nonpositive=unaffected,sufficient_for_strict_terminal_positivity=not unaffected,
        geometry_proved=False,low9_concentration_requires_revalidation=rank<=9 or any(k<=9 for k in fixed))


def ideal(old):
    domain=sorted(((low,)+p if low else p) for low in range(9) for p in high_parts(58-low))
    prior={tuple(r['ranks']):r for r in old['rows']}
    assert len(domain)==len(set(domain))==len(prior)==len(old['rows'])==1890
    assert set(domain)==set(prior)
    env={int(k):F(v) for k,v in old['current_small_rank_envelopes'].items()}
    env.update({r['order']-1:F(r['envelope']) for r in old['exact_family_tail']})
    assert env[13]==F(41,7)
    newenv={**env,13:F(40,7)};rows=[]
    for rs in domain:
        r=prior[rs];h=r['h_polynomial'];v=prod(x+1 for x in rs)
        assert sum(rs)==58 and len(h)==59 and h==h[::-1] and sum(h)==v==r['link_vertices']
        assert all(isinstance(x,int) and x>0 for x in h)
        g=sum(min(j,58-j)*x for j,x in enumerate(h));curvature=F(5*g-58*v,2)
        assert g==r['G'] and curvature==F(r['curvature'])
        assert 2*g==58*v-sum(abs(2*j-58)*x for j,x in enumerate(h))
        oldcost=sum(env[x] for x in rs);cost=sum(newenv[x] for x in rs)
        assert oldcost==F(r['normalized_cost']) and oldcost-cost==F(rs.count(13),7)
        before=curvature/v-oldcost;after=curvature/v-cost
        assert before==F(r['normalized_margin']) and v*before==F(r['absolute_margin'])
        assert after==before+F(rs.count(13),7)
        rows.append(dict(ranks=rs,link_vertices=v,G=g,curvature=curvature,old_normalized_cost=oldcost,
            normalized_cost=cost,order14_saving=F(rs.count(13),7),old_normalized_margin=before,
            normalized_margin=after,absolute_margin=v*after))
    out=summary(rows);assert out['negative_count']==7 and out['zero_count']==0
    assert out['normalized_minimum']==-F(5241,57460)
    assert old['summary']['negative_count']==10 and old['summary']['zero_count']==0
    fixed={12:F(1),14:F(1,4),15:F(1,2),16:F(3,4)}
    scenarios=[conditional(rows,{9:F(1,2)}),conditional(rows,fixed),
        conditional(rows,{12:F(1),14:F(1,4),16:F(5,4)})]
    assert all(x['summary']['negative_count']==x['summary']['zero_count']==0 for x in scenarios)
    return dict(complete_profiles=1890,complete_low9_domain=True,saved_polynomials_checked=True,
        fresh_convolution_performed=False,old_negative_count=10,rows=rows,summary=out,
        current_envelopes=newenv,unchanged_exact_family_tail=old['exact_family_tail'],
        single_rank_thresholds=[threshold(rows,r) for r in (9,10,12,14,15,16)],
        coordinated_thresholds=[threshold(rows,r,{k:v for k,v in fixed.items() if k!=r}) for r in fixed],
        conditional_scenarios=scenarios,negative_absolute_margin_not_extended_to_arbitrary_profiles=True)


def finite(d,local):
    z=d['conditional_neighbor_preflight'];assert z['neighbor_half_increase_is_conditional_preflight']
    assert local['new_coupon']['canonical_large_node']==9
    assert local['new_coupon']['component_pair_increment']=='1/2'
    assert local['new_coupon']['old_actual_component_pair_cost']=='3/2'
    assert local['new_coupon']['new_actual_component_pair_cost']==z['A1_D13_component_cost']=='2'
    assert d['ideal_weights_unchanged'] and d['delta']=='1/16'
    assert z['exact_maximum']=='795/2' and z['curvature']=='783/2'
    assert z['fixed_vectors']==len(z['alternatives'])==961 and z['interaction_keys']==960
    assert len(z['tables'])==1920 and z['simultaneous_type_pairs']==32400
    assert max(r['sixteenths_upper'] for r in z['alternatives'])==6360
    assert max(r['sixteenths_actual'] for r in z['alternatives'])==6360
    above=[r for r in z['alternatives'] if r['sixteenths_upper']>6264]
    assert above==z['above_curvature'] and len(above)==14
    assert z['max_with_rank13']=='1563/4' and z['max_with_rank10']=='386'
    assert z['A2_A13_component_cost']==z['A2_D13_component_cost']=='39/8'
    assert z['A1_A13_component_cost']=='25/16' and z['context_pair_payment']=='1/8'
    assert z['context_total']=='k*(k-1)/4 if whole A10 and at least two whole A1 are present; otherwise0'
    assert z['complete_context_count'] and z['actual_D15_cost']=='201/2'
    keep=[k for k in z if k not in ('tables','alternatives','neighbor_half_increase_is_conditional_preflight')]
    return dict(source_certificate=FINITE,source_field='conditional_neighbor_preflight',
        source_preflight_geometric_hypothesis_discharged=True,**{k:z[k] for k in keep},
        finite_deficit='6',new_DP_performed=False,saved_Bellman_replay_performed=False,
        inherited_fixed_rank17='A17',inherited_residual_rank17='D17',
        complete_certificate_scope='The bound source includes all961 vectors,1920 ordinary/no-rank2 tables, type domination, last-item inequalities and actual attainment; this binder checks its exact interfaces without rerunning it.')


def literal(docs):
    r=docs[FINITE]['combined_literal_interface'];dd=docs[D]['current_literal_interface']
    assert r['configurations']==dd['configurations']==216 and dd['changed']==[]
    assert docs[B]['current_literal_interface']['changed_corners']==docs[C]['current_literal_interface']['changed_corners']==[]
    key=lambda row:(row['source'],json.dumps(row['case'],sort_keys=True))
    before={key(x):x for x in dd['rows']};assert len(before)==216
    rows=[]
    for row in r['rows']:
        base=before[key(row)];assert len(base['corners'])==len(row['corners'])
        total=F()
        for c,p in zip(row['corners'],base['corners']):
            assert c['marks']==p['marks'] and F(c['previous_lower'])==F(p['retained_payment_lower_bound'])
            assert F(p['debit'])==F(p['increment'])==0
            assert F(c['retained_payment_lower_bound'])==F(c['previous_lower'])+F(c['delta'])>=0
            total+=F(c['retained_payment_lower_bound'])
        assert total==F(row['total_lower_bound'])>=1
        rows.append(dict(source=row['source'],case=row['case'],total_lower_bound=total))
    assert r['changed_roles']=={'reduce_A2_D13_middle_or_center':2,'context_A2_A2_with_A10_and_two_A1':1}
    return dict(configurations=216,root_role_changes=r['changed_roles'],new_order14_literal_changes=0,
        rows=rows,triangle_minimum=2,quadrilateral_minimum=1,
        universal_actual_corner_proof_required_beyond_literals=True)


def run():
    assert FROZEN_D_MD!='PENDING','Wait for the final local report.'
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    assert sha256((ROOT/f'results/{D}.md').read_bytes()).hexdigest()==FROZEN_D_MD
    assert sha256((ROOT/f'results/{FINITE}.md').read_bytes()).hexdigest()==FROZEN_FINITE_MD
    assert docs[OLD]['dimension_bound']==59 and docs[D]['current_order14_family_costs']==dict(A=70,B=80,C=80,D=80)
    assert docs[B]['new_cost']==docs[C]['new_cost']==docs[D]['new_cost']=='80'
    assert docs[FRONT]['finite_curvature']=='783/2' and docs[FRONT]['strict_global_constant']=='4350'
    paths={Path(__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_ORDER14_CONTEXTUAL_AMBIENT59_FRONTIER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        ambient_dimension=59,link_dimension=58,completed_unrestricted_dimension_bound=59,new_dimension_bound=None,
        current_order14_family_costs=docs[D]['current_order14_family_costs'],order14_envelope='40/7',
        finite=finite(docs[FINITE],docs[D]),ideal=ideal(docs[FRONT]['ideal']),
        combined_literal_interface=literal(docs),finite_curvature='783/2',strict_global_constant='4350',
        constant_counted_once=True,triangle_minimum=2,quadrilateral_minimum=1,
        current_rule_does_not_close_ambient59=True,
        scope='Compatible proved B/C/D14 and contextual finite rule. The previous optional A1/D13-neighbor count now has its local geometric proof. Complete1890 ideal saved-row updates and exact finite certificate interfaces; no new DP, convolution, historical census, negative-surplus extension or dimension58 theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'ideal negatives',out['ideal']['summary']['negative_count'],flush=True)
