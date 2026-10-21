#!/usr/bin/env python3
"""Final dimension57 composition with independent ideal arithmetic checks.

No finite Bellman table, geometric census or local matrix computation is
rerun. The complete preceding finite proof and the rank58 local repair
are explicit pinned mathematical inputs, not inferred from sampled faces.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_d13_context_dimension57_20261018'
LOCAL='fresh_d13_three_rank15_context_transfer_20261018'
B13='fresh_b13_endpoint_retention_20261018'
COUNT='fresh_d17_a16_coupled_ambient58_frontier_20261018'
IDEAL='fresh_order10_low8_ambient58_frontier_20261018'
TAIL='fresh_a15_d15_dimension58_20261018'
PINS={
 LOCAL:'91e93c3f7ff5eb4e64b435eea2892bf20af08711aed131110f40e898a9322014',
 B13:'fedf5747f39c26e006084a7af623402205d5114841e31a941d29ef5bd23a384d',
 COUNT:'6f510a745724c64197d4b202068d6e29d2f3ebe381821db0c5f67db07fadec97',
 IDEAL:'77149e2f7e6d6a96f720bbca5bbef5ba7647ca1039be6386bac2c9f0216cf983',
 TAIL:'c153463b2b9b043eb3f076db1e6504548e0fc001072dd1455c7bd06f2014f845'}

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
@lru_cache(None)
def high_parts(n,least=8):
    if n==0:return((),)
    return tuple((r,)+p for r in range(least,n+1)for p in high_parts(n-r,r))
def polynomial(rs):
    h=[1]
    for r in rs:
        out=[0]*(len(h)+r)
        for i,a in enumerate(h):
            for j in range(r+1):out[i+j]+=a
        h=out
    return h

def check_ideal(ii):
    rows=ii['rows'];domain=sorted((low,)+p if low else p for low in range(8)for p in high_parts(57-low))
    assert len(domain)==len(set(domain))==len(rows)==ii['complete_profiles']==2552
    assert sorted(tuple(r['ranks'])for r in rows)==domain
    env={int(k):F(v)for k,v in ii['current_envelopes'].items()}
    assert all(env[r]==F(r,2)for r in range(1,9))and env[9]==F(71,16)
    checked=[]
    for row in rows:
        rs=row['ranks'];h=polynomial(rs);v=prod(r+1 for r in rs)
        assert sum(rs)==57 and h==row['h_polynomial']and len(h)==58 and h==h[::-1]
        assert sum(h)==v==row['link_vertices']and all(x>0 for x in h)
        g=sum(min(j,57-j)*a for j,a in enumerate(h));assert g==row['G']
        assert 2*g==57*v-sum(abs(2*j-57)*a for j,a in enumerate(h))
        curv=F(5*g-57*v,2);cost=sum(env[r]for r in rs);margin=curv/v-cost
        assert curv==F(row['curvature'])and cost==F(row['normalized_cost'])
        assert margin==F(row['normalized_margin'])>0 and margin*v==F(row['absolute_margin'])>=36
        checked.append((margin,margin*v))
    assert min(a for a,b in checked)==F(87,10240)and min(b for a,b in checked)==36
    cc=ii['concentration'];assert cc['cutoff']==8 and cc['changed_rank9_never_transferred']
    transfers=cc['complete_low_rank_transfers'];assert len(transfers)==28
    assert {(r['a'],r['b'])for r in transfers}=={(a,b)for a in range(1,8)for b in range(a,8)}
    for row in transfers:
        a,b=row['a'],row['b']
        def law(x,y):return[F(sum(i+j==k for i in range(x+1)for j in range(y+1)),(x+1)*(y+1))for k in range(x+y+1)]
        old,new=law(a,b),law(a-1,b+1)
        assert sum(old)==sum(new)==1 and sum(i*(y-x)for i,(x,y)in enumerate(zip(old,new)))==0
        gaps=[sum(max(i-t,0)*(y-x)for i,(x,y)in enumerate(zip(old,new)))for t in range(a+b+1)]
        assert list(map(str,gaps))==row['stop_loss_gaps']and min(gaps)>=0
        assert (a+1)*(b+1)==row['old_product']>row['new_product']==a*(b+2)
        assert F(a,2)+F(b,2)==F(a-1,2)+F(b+1,2)
    assert ii['actual_charge_bounded_above_by_envelope']and ii['absolute_surplus_valid_for_every_actual_cusp']
    return dict(complete_profiles=2552,polynomials_reconstructed=2552,low_rank_stop_loss_checks=28,
        normalized_minimum='87/10240',absolute_surplus='36',actual_charge_at_most_envelope=True,
        extension='Forward low8 transfers do not increase normalized margin or V, preserve normalized cost, and end at a positive terminal margin. Therefore original absolute surplus is at least36.',
        no_negative_margin_extension=True)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in PINS}
    bindings=0
    for s,d in docs.items():
        assert d['status'].startswith('PASS')and digest(f'results/{s}.md')==PINS[s]
        for p,h in d['sources'].items():assert digest(p)==h,(s,p);bindings+=1
    local,b13,count,ideal,tail=(docs[s]for s in(LOCAL,B13,COUNT,IDEAL,TAIL))
    for s in(B13,COUNT):
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            p=f'{folder}/{s}.{ext}';assert local['sources'][p]==digest(p)
    for folder,ext in(('scripts','py'),('results','json'),('results','md')):
        p=f'{folder}/{IDEAL}.{ext}';assert count['sources'][p]==digest(p)
    assert local['ambient_dimension']==58 and local['changed_pairs']==[[0,11],[0,12]]
    assert b13['changed_pair']==[0,12]and b13['old_B13_cost']=='78'and b13['new_B13_cost']=='77'
    assert not b13['ideal_weights_changed']and b13['added_coupons']==[]
    selected=local['selected_support'];assert(selected['patterns'],selected['maps'],selected['affected_contexts'],selected['affected_maps'])==(801,862,33,76)
    assert Counter(r['method']for r in local['selected_proofs'])=={
        'radical_full_core_seal':4,'genuine_new_fork_leaf_projection_exclusion':5,
        'untouched_small_factor':12,'one_contact_inverse_rank_bound':12}
    assert len(local['quiet_repairs'])==24 and len(local['one_sided_whole_models'])==2
    for r in local['one_sided_whole_models']:assert F(r['quiet_half'])+F(r['new_half'])==F(r['total'])==1
    ci=local['current_interfaces'];assert ci['literal_configurations']==216 and len(ci['later_corners'])==564
    assert len(ci['D14_variants'])==24 and ci['D12_positive_domain']==332 and ci['no_changed_same_component_role_in_these_domains']
    assert count['exact_actual_finite_maximum']=='378'and count['finite']['all_other_38_classes_at_most_curvature']
    assert count['representative_cover']['inherited_expanded_classes']==3298
    assert count['finite']['actual_fixed_fibers']==748
    old={tuple(sorted(tuple(t)for t in fs))for fs in count['finite']['exact_obstructing_forests']}
    wanted={tuple(sorted(((k,13),)+p))for k in('B','D')for p in combinations_with_replacement((('A',15),('D',15)),3)}
    assert old==wanted and len(old)==8
    ff=local['finite'];assert ff['maximum']=='377'and ff['ideal_weights_unchanged']
    got=set()
    for r in ff['all_old_costly_forests']:
        fs=tuple(sorted(tuple(t)for t in r['forest']));got.add(fs)
        assert sum(n for k,n in fs)==58 and r['old']=='378'and r['new']=='377'
        assert 77+3*100==F(r['new'])
    assert got==old and len(ff['all_old_costly_forests'])==8
    expected={tuple(sorted(tuple(fixed)+p))for fixed in((('D',28),),(('A',12),('D',16)))
        for p in combinations_with_replacement((('A',15),('D',15)),2)}
    got=set()
    for r in ff['complete_increasing_families']:
        fs=tuple(sorted(tuple(t)for t in r['forest']));got.add(fs);co=Counter(fs)
        assert sum(n for k,n in fs)==58 and co['A',15]+co['D',15]==2
        new=F(170+200)+F(1,2)if co['D',28]else F(66)+F(217,2)+200+1
        assert F(r['new'])==new==F(r['previous'])+F(r['increment'])<=377
    assert got==expected and len(got)==len(ff['complete_increasing_families'])==6
    ii=check_ideal(ideal['ideal'])
    assert count['positive_ideal_profiles']==2552 and count['actual_cusp_surplus']=='36'
    d=57;g0=sum(min(j,d-j)for j in range(d+1));beta=F(5*g0-d*(d+1),2)
    correction=sum(2*j-d for j in range(d+1)if 2*j>d);constant=5*correction
    assert(g0,beta,correction,constant)==(812,F(377),841,4205)
    assert beta==F(count['finite_curvature'])==F(ideal['finite_curvature'])==F(ff['curvature'])
    assert constant==int(count['strict_global_constant'])==int(ideal['strict_global_constant'])==int(ff['strict_global_constant'])
    assert ideal['necessary_ambient58_inequality']=='(47/8)*v_f >= 4205 + 36*v_i'
    assert ideal['constant_counted_once']and tail['dimension_bound']==58
    assert tail['hypotheses']==['finite-sided real hyperbolic Coxeter polytope','finite volume']
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in PINS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return dict(status='PASS_COMPLETE_D13_CONTEXT_DIMENSION57_COMPOSITION',sources={p:digest(p)for p in sorted(paths)},
        hypotheses=tail['hypotheses'],excluded_dimension=58,inherited_dimension_bound=58,dimension_bound=57,
        arithmeticity_required=False,compactness_required=False,simple_cusps_required=False,facet_count_bound_required=False,
        triangle_minimum=2,quadrilateral_minimum=1,finite_actual_maximum='377',finite_curvature='377',finite_slack='0',
        complete_predecessor_actual_obstructions=8,complete_new_increasing_forests=6,all_other_actual_forests_at_most_curvature=True,
        ideal=ii,strict_global_constant='4205',constant_counted_once=True,
        contradiction='0 >= 4205 + 36*v_i > 0',
        local_scope='Dimension58 only: three whole unmarked A15/D15 witnesses for D13 halves; two such witnesses for the distinct D28 or A12/D16 half. Complete arbitrary-support selected proof and universal output repairs; ideal and actual-face triangle outputs retained.',
        finite_scope='Frozen complete3298 representative cover/748 actual fibers plus both increasing families; the new local count checks all eight prior obstructions and every possible increasing rank58 family. No new Bellman replay.',
        inherited_binding_checks=bindings,higher_dimension_tail='The completed bound58 excludes every dimension at least59; the present contradiction excludes58.',
        new_finite_volume_example=False,external_refereeing=False,formal_proof_assistant_verified=False,
        verification_scope='Hash/interface/finite-closure composition; independent complete2552 ideal polynomial products, exact margins and28 stop-loss transfers; independent curvature/extrema arithmetic. No historical geometry, native census or finite DP replay.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'dimension',out['dimension_bound'],'finite',out['finite_actual_maximum'],'surplus',out['ideal']['absolute_surplus'],'constant',out['strict_global_constant'],flush=True)
