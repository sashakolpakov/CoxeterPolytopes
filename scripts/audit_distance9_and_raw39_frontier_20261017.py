#!/usr/bin/env python3
"""Incremental archive of distance9, new-mass packing and raw39 actual corners.

The recursive archive walker is inherited from the frozen previous auditor.
Every inherited file is rehashed before its earlier dependency graph is reused.
New small checks are independent Fraction compositions, not geometry proofs.
"""
import argparse
import ast
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import re
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/distance9_and_raw39_frontier_20261017.md'
PREVIOUS='results/distance10_and_mixed_corner_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_distance10_and_mixed_corner_frontier_20261017.py'
OUT=ROOT/'results/distance9_and_raw39_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
STEMS=(
 'fresh_selective_distance9_dimension108_20261017',
 'fresh_selective_distance9_dimension108_20261017_review',
 'fresh_distance9_low_core_20261017',
 'fresh_distance9_supporting_extensions_review_20261017',
 'fresh_mixed_actual_corner_full_model_20261017_local_cones',
 'fresh_mixed_actual_corner_full_model_20261017_joint_ports',
 'fresh_raw39_actual_corner_compatibility_20261017',
 'fresh_zero_arc_cap_high_r_20261017',
)
EXTRA_REPORTS=()

@lru_cache(None)
def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())

def key(r):return(r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['match'])))
def masses(rows):
    d={i:F(a)for i,a in rows};assert len(d)==len(rows)and all(a>0 for a in d.values());return d

def arithmetic_checks(local,joint,raw):
    model=read('fresh_mixed_actual_corner_full_model_20261017')
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    assert local['face_multipliers']==joint['face_multipliers']==raw['face_multipliers']==dual['face_multipliers']
    assert local['inherited_model_sha256']==model['prepared']['model_sha256']
    alpha=masses(local['face_multipliers']);assert len(alpha)==4661
    initial=masses(local['initial_vertex_multipliers']);mu=masses(local['vertex_multipliers'])
    assert initial==masses(dual['vertex_multipliers'])and len(initial)==780
    assert local['global_optimizer_invocations']==local['replay_optimizer_invocations']==0
    assert local['complete_positive_simplex_forest_coverage']and not local['reused_old_fixed_alpha_covers']
    columns=local['block_orbit_columns'];literal=local['literal_block_orbits']
    assert len(columns)==len(literal)==5173
    marginals=Counter();features=local['positive_simplex_features'];feature_faces=set();corners=0
    for f in features:
        fi=f['face'];assert fi in alpha and fi not in feature_faces;feature_faces.add(fi)
        d=f['simplex_dimension'];assert 1<=d<=4 and sum(m for j,m in f['block_orbits'])==d+1
        for j,m in f['block_orbits']:
            assert m>0 and columns[j]['block_size']==d;marginals[j]+=m*alpha[fi];corners+=m
    assert len(features)==2753 and corners==8483
    for i,(c,l)in enumerate(zip(columns,literal)):
        assert c['index']==l['index']==i and c['vertex']==l['vertex']and c['block_size']==l['block_size']
        assert len(l['literal_masks'])==len(set(l['literal_masks']))>0
        assert all(0<mask<1<<13 and mask.bit_count()==c['block_size']for mask in l['literal_masks'])
    rebuilt=dict(initial);covered=set();terms=0;literal_blocks=0;statuses=Counter()
    for row in local['local_cones']:
        v=row['vertex'];assert v not in covered;covered.add(v);statuses[row['status']]+=1
        want=dict((j,F(a))for j,a in row['required_marginals'])
        assert want=={j:a for j,a in marginals.items()if columns[j]['vertex']==v}
        assert row['literal_block_orbits']==[literal[j]for j in want]
        got=Counter();total=F();cert=row['certificate']
        for term in cert['packing_mixture']:
            mass=F(term['mass']);assert mass>0;used=0;terms+=1
            for j,mask in term['blocks']:
                assert j in want and mask in literal[j]['literal_masks']and not used&mask
                used|=mask;got[j]+=mass;literal_blocks+=1
            total+=mass
        assert total==F(cert['exact_cover_mass'])and all(got[j]>=a for j,a in want.items())
        assert F(row['old_vertex_mass'])==initial.get(v,F())
        rebuilt[v]=max(initial.get(v,F()),total)
        assert F(row['repaired_vertex_mass'])==rebuilt[v]
    assert terms==4363 and len(covered)==760 and covered=={c['vertex']for c in columns}
    assert {i:a for i,a in rebuilt.items()if a}==mu
    assert sum(mu.get(i,F())>initial.get(i,F())for i in set(mu)|set(initial))==205
    assert statuses=={'EXACT_LOCAL_OPTIMUM':675,'EXACT_COVER_BRACKET':85}
    for c in local['new_exact_charge_columns']:
        load=sum((alpha[i]*F(a)for i,a in c['positive_simplex_coefficients']),F())
        assert load==F(c['exact_load'])>F(c['old_vertex_mass'])and load<=mu[c['vertex']]
    assert len(local['new_exact_charge_columns'])==132
    assert local['all_exported_model_columns_checked']==dict(all_columns=432021,vertex_cost_signed=52,zero_budget_signed=96176)
    numerator=F(dual['objective_numerator']);assert F(local['objective_numerator'])==numerator
    assert sum(mu.values(),F())==F(local['normalizer'])
    assert numerator/sum(mu.values(),F())==F(local['bound'])==F(247898037929594052000000000000000,19983545698211739397486458908069)>12
    assert masses(joint['initial_vertex_multipliers'])==mu and joint['optimizer_invocations']==0
    jm=masses(joint['vertex_multipliers']);rebuilt=dict(mu);nc=0
    for c in joint['critical_columns']:
        load=sum((alpha[i]*F(a)for i,a in c['positive_face_coefficients']),F());cost=F(c['vertex_cost']);v=c['vertex']
        assert load==F(c['load'])and cost>0 and load>cost*mu.get(v,F())
        assert F(c['required_vertex_mass'])==load/cost and F(c['old_vertex_mass'])==mu.get(v,F())
        assert load<=cost*jm[v];rebuilt[v]=max(rebuilt.get(v,F()),load/cost);nc+=1
    assert rebuilt==jm and nc==183 and sum(jm.get(i,F())>mu.get(i,F())for i in set(jm)|set(mu))==33
    stages=joint['complete_family_stages'];assert [(r['subset_size'],r['inequalities'],r['violated_columns'])for r in stages]==[(3,58468,25),(4,120593,24),(5,2881140,134)]
    assert sum(r['inequalities']for r in stages)==3060201
    assert joint['counts']==dict(actual_quad_faces=473,actual_quad_vertex_occurrences=3232,block_orbits=5173,explicit_provenance_Q_rows=433,increased_vertex_masses=33,positive_faces=4661,simplex_faces=2753,tested_forests=764,total_inequalities=3060201,unique_saved_raw23_raw25_Q_rows=40,unselected_Q_rows_with_zero_new_coefficients=1435,violated_columns=183)
    assert len(joint['zero_coefficient_unselected_Q_rows'])==1435
    assert F(joint['objective_numerator'])==numerator and sum(jm.values(),F())==F(joint['normalizer'])
    assert numerator/sum(jm.values(),F())==F(joint['bound'])==F(247898037929594052000000000000000,20013747715459140911449997327399)>12
    assert masses(raw['initial_vertex_multipliers'])==jm
    rm=masses(raw['vertex_multipliers']);rebuilt=dict(jm);violations=0
    for c in raw['columns']:
        load=sum((alpha[i]*F(a)for i,a in c['positive_face_coefficients']),F());cost=F(c['vertex_cost']);v=c['vertex']
        assert load==F(c['load'])and F(c['baseline_vertex_mass'])==jm.get(v,F())and cost>0
        assert c['violated']==(load>cost*jm.get(v,F()));violations+=c['violated']
        assert load<=cost*rm[v];rebuilt[v]=max(rebuilt.get(v,F()),load/cost)
    assert rebuilt==rm and len(raw['columns'])==129 and violations==4
    changed=[[i,str(jm.get(i,F())),str(rm[i])]for i in sorted(rm)if rm[i]>jm.get(i,F())]
    rep=raw['fixed_witness_repair'];assert changed==rep['vertex_updates']and len(changed)==3
    assert [r[0]for r in changed]==[141,277,5587]
    assert rep['face_masses_unchanged']and rep['optimizer_invocations']==0
    assert F(rep['objective_numerator'])==numerator and F(rep['baseline_bound'])==F(joint['bound'])
    assert sum(rm.values(),F())==F(rep['normalizer'])and numerator/sum(rm.values(),F())==F(rep['bound'])==F(247898037929594052000000000000000,20020725665111814911449997327399)>12
    assert not rep['bound_below12']and not raw['omitted_ambiguous_raw39_rows']
    assert raw['raw_source']==39 and raw['canonical_core']==29 and raw['dotted_square']==['41/2','17/2']
    assert len(raw['selected_faces'])==31 and sum(len(f['corners'])for f in raw['selected_faces'])==124
    assert sum(f['provenance_kind']=='existing_explicit_map'for f in raw['selected_faces'])==1
    assert sum(f['provenance_kind']=='unique_saved_context'for f in raw['selected_faces'])==30
    assert len(raw['forest_records'])==88 and sum(len(r['literal_stars'])for r in raw['forest_records'])==1246
    assert sum(r['tested_literal_pairs']for r in raw['forest_records'])==21115
    assert sum(len(r['incompatibility_certificates'])for r in raw['forest_records'])==21047
    return dict(unchanged_positive_faces=4661,local_cover_forests=760,local_cover_orbits=5173,literal_packing_terms=terms,literal_packing_blocks=literal_blocks,new_local_cuts=132,local_bound=local['bound'],selected_joint_Q=473,joint_zero_coefficient_Q=1435,joint_inequalities=3060201,joint_cuts=183,joint_bound=joint['bound'],raw39_columns=129,raw39_violations=4,raw39_vertex_increases=3,final_arithmetic_bound=rep['bound'])

def theorem_checks(p,review,low,extensions):
    assert p['dimension_bound']==108 and p['inherited_dimension_bound']==116
    assert p['unrestricted_finite_volume']and not p['arithmetic_hypothesis']and p['cusp_restriction']is None
    assert p['excluded_dimensions']==list(range(109,117))
    assert review['status']=='PASS_INDEPENDENT_FRACTION_SELECTIVE_DISTANCE9_DIMENSION108_REVIEW'and review['dimension_bound']==108
    rows=p['reduced_profiles'];expected=set()
    for d in range(108,116):
        for k in range(d//18+1):
            for low_rank in range(18):
                high=d-18*k-low_rank
                if high==0 or high>=19:expected.add(tuple(sorted([18]*k+([low_rank]if low_rank else[])+([high]if high else[]))))
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8);assert beta>9*n
    assert len(rows)==len(expected)==756 and {tuple(sorted(r['ranks']))for r in rows}==expected
    def scalar(r):
        rs=r['ranks'];d=sum(rs);v=prod(t+1 for t in rs);h=[1]
        for t in rs:
            q=[0]*(len(h)+t)
            for i,a in enumerate(h):
                for j in range(t+1):q[i+j]+=a
            h=q
        assert h==h[::-1]and sum(h)==v
        G=sum(a*min(i,d-i)for i,a in enumerate(h));curvature=F(5*G-d*v,2*v)
        cost=sum((min(F(9),F(t,2))for t in rs),F())
        if 'G'in r:assert r['G']==G and r['link_vertices']==v and r['link_dimension']==d
        assert F(r['normalized_curvature'])==curvature and F(r['normalized_cost'])==cost
        assert F(r['margin'])==curvature-cost
        return curvature-cost
    assert min(map(scalar,rows))==F(45873,2476099)>0
    assert scalar(p['boundary_at_n108'])==F(-4659,260642)
    br=review['budgets'];assert len(br['profiles'])==756
    byshape={tuple(sorted(r['ranks'])):(F(r['normalized_curvature']),F(r['normalized_cost']),F(r['margin']))for r in rows}
    assert byshape=={tuple(sorted(r['ranks'])):(F(r['normalized_curvature']),F(r['normalized_cost']),F(r['margin']))for r in br['profiles']}
    assert scalar(br['boundary'])==F(-4659,260642)
    facet=[]
    for start in(79,80):
        vals=[]
        for t in range(4):
            d=start+2*t;beta=F(d*(d-6),8)if d%2==0 else F((d-1)*(d-5),8);vals.append(8*(beta-9*(d+1)))
        c2=(vals[2]-2*vals[1]+vals[0])/2;co=[vals[0],vals[1]-vals[0]-c2,c2]
        assert all(c>0 for c in co)and all(vals[t]==sum(c*t**i for i,c in enumerate(co))for t in range(4))
        facet.append(list(map(int,co)))
    assert facet==p['compact_facet_parity_slack_times8']==br['compact_facet_slack_coefficients']==[[12,160,4],[88,164,4]]
    assert p['compact_facet_ambient_bound']==review['compact_facet_ambient_bound']==79
    c,j,iso=p['census'],p['disjoint_obstructions'],p['isolation'];rc,ri=review['classification'],review['isolation']
    assert c['marked_contexts']==rc['marked_contexts']==384 and c['allowed_ports']==rc['allowed_ports']==378
    assert len(c['ports'])==len(rc['ports'])==1262 and c['ports']==rc['ports']
    assert c['overlap_generated']==rc['candidate_matchings']==312316
    assert c['consistent_overlaps']==rc['consistent_matchings']==20718
    cc=c['overlap_counts'];assert cc['inverse_rejection']==18721 and cc['edge_pd_rejection']==1048
    assert rc['overlap_counts']=={('edge_rejection'if k=='edge_pd_rejection'else k):v for k,v in cc.items()}
    left={key(r):r for r in c['overlap_records']};right={key(r):r for r in rc['overlap_records']}
    assert len(left)==len(right)==20718 and left.keys()==right.keys()
    for k,a in left.items():
        kind='inverse_rejection'if'inverse_rejection'in a else'edge_rejection'if'edge_pd_rejection'in a else a['surviving_type']
        assert right[k]['kind']==kind
    assert j['candidate_joins']==len(j['rows'])==len(rc['disjoint_joins'])==4442
    assert j['white_rejections']==sum('white_rejection'in r for r in j['rows'])==1103
    assert j['remaining_joins']==3339 and j['remaining_context_pair_multiplicities']=={'1':2007,'2':348,'3':180,'4':24}
    for a,b in zip(j['rows'],rc['disjoint_joins']):
        assert a['left']==b['left']and a['right']==b['right']
        assert('white_rejection'in a)==('white_obstruction'in b)
    assert len(j['underlying_inertia_types'])==len(rc['surviving_disjoint_types'])==85
    assert all(r['signature'][1]==2 for r in j['underlying_inertia_types'])
    def mk(r):return(r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['edges'])))
    assert len(j['multiple_edge_obstructions'])==len(rc['multiple_edge_obstructions'])==1332
    assert {mk(r)for r in j['multiple_edge_obstructions']}=={mk(r)for r in rc['multiple_edge_obstructions']}
    assert len(iso['cases'])==len(ri['preferred_corners'])==949 and iso['component_order_cap']==17
    pc={key(r):r for r in iso['cases']};rr={key(r):r for r in ri['preferred_corners']};assert pc.keys()==rr.keys()
    short=extra=0;methods=Counter()
    for k,a in pc.items():
        b=rr[k];methods[a['method']]+=1
        for field in('preferred_corner','omitted_black','blocks','possible_contacts','marked_component','method','component_order_bound'):
            assert a[field]==b[field],(k,field)
        assert a['component_order_bound']<=17
        if a['method']=='radical':assert a['shortcut_bound']==b['shortcut_bound']<=9
        elif a['method']=='relevant_component':
            for field in('relevant_nodes','contact_domains','end_sets'):assert a[field]==b[field],(k,field)
            bx={r['node']:r for r in b['extension_cases']};assert len(bx)==len(a['extension_cases'])
            for r in a['extension_cases']:
                z=bx[r['node']]
                for field in('projection','rho','alone_tail_bound'):assert r[field]==z[field]
                aa,bb=r['coexisting_outside'],z['joint'];assert bool(aa)==bool(bb)
                if aa:
                    assert aa['outside_node']==bb['node']
                    for field in('determinant','kind','tail'):assert aa[field]==bb[field]
                    if 'inverseGG'in aa:assert aa['inverseGG']==bb['inverseGG']
                own={q['node']:q for q in z['pair_checks']}
                for t in r['pair_exclusions']:
                    q=own[t['node']];assert q['reason']==('forced_orthogonal'if t['kind']=='orthogonal_by_gap'else t['kind'])
                    if 'projection'in t:assert t['projection']==q['projection']and t['residual_product']==q['residual_product']
        else:
            assert a['end_sets']==b['end_sets']
            assert [{k:v for k,v in r.items()if k!='contact'}for r in a['single_contact_candidates']]==b['single_contact_candidates']
        cross=a['extra_cross_bond_checks'];assert cross['shortcut_rejections']==b['extra_cross_shortcuts']
        assert len(cross['edge_pd_z_rejections'])==len(b['extra_cross_pd_rejections'])
        short+=cross['shortcut_rejections'];extra+=len(cross['edge_pd_z_rejections'])
    assert methods['relevant_component']==69 and sum(methods.values())==949
    assert short==iso['extra_cross_shortcuts']==ri['extra_cross_shortcuts']==60831
    assert extra==iso['extra_cross_pd_rejections']==ri['extra_cross_pd_rejections']==888
    assert review['saved_comparison']['saved_nonnegative_witnesses']==1103+1332+1048+888==4371
    ec={key(r):r for r in extensions['contexts']};assert len(ec)==69
    assert ec.keys()=={k for k,r in pc.items()if r['method']=='relevant_component'}
    for k,e in ec.items():
        a=pc[k]
        for field in('preferred_corner','omitted_black','possible_contacts','marked_component','relevant_nodes','contact_domains','end_sets'):assert a[field]==e[field],(k,field)
        assert a['component_order_bound']==e['component_bound']<=16
        bynode={z['node']:z for z in e['extensions']};assert len(bynode)==len(a['extension_cases'])
        for r in a['extension_cases']:
            z=bynode[r['node']]
            for field in('projection','rho','alone_tail_bound'):assert r[field]==z[field]
            x,y=r['coexisting_outside'],z['coexisting_outside'];assert bool(x)==bool(y)
            if x:
                assert x['outside_node']==y['outside_node']and x['determinant']==y['determinant']and x['kind']==y['kind']
                assert x['tail']==y['tail_bound']
                if 'inverseGG'in x:assert F(x['inverseGG'])==F(y['inverse'][0][0])
    assert extensions['counts']==dict(contexts=69,core_types=9,relevant_root_cases=24,no_relevant_root_contexts=46,max_marked_component_order=16)
    assert low['status']=='PASS_COMPLETE_DISTANCE9_SMALL_CORE_DOMAIN'and p['low_core_lemma']==STEMS[2]
    for order,total,types in((4,392,32),(5,1050,14)):
        d=low['classification'][str(order)];assert d['order']==order and d['all_labelled_minimal_indefinite']==total and d['canonical_types']==len(d['diagrams'])==types
        assert sum(r['labelled_multiplicity']for r in d['diagrams'])==total
    lj=low['disjoint_and_overlap'];assert len(lj['joins'])==104 and len(lj['overlaps'])==30 and len(lj['multiple_edge_obstructions'])==6
    assert lj['low_first_simple_joins']==95 and lj['low_first_possible_double_joins']==9 and lj['low_first_two_negative']==93 and lj['low_first_white_rejections']==11
    assert lj['low_first_overlap_candidates']==30 and lj['low_first_surviving_overlaps']==6 and lj['ordered_surviving_overlaps']==12
    assert len(low['extensions']['cases'])==6 and low['extensions']['component_order_cap']==12
    assert all(r['selected_corner']['component_order_bound']<=12 for r in low['extensions']['cases'])
    return dict(scalar_profile_checks=756,independent_profile_matches=756,unrestricted_dimension_bound=108,compact_facet_ambient_bound=79,independent_port_matches=1262,independent_overlap_matches=20718,independent_corner_matches=949,high_disjoint_joins=4442,multiple_edge_obstructions=1332,reviewed_nonnegative_witnesses=4371,extra_cross_shortcuts=60831,extra_cross_pd_obstructions=888,independent_supporting_extension_contexts=69,small_core_types=46,small_core_join_cases=104,small_core_surviving_ordered_overlaps=12)


def construction_checks(high):
    assert high['status']=='PASS_EXACT_CONTINUOUS_ZERO_ARC_HIGH_R_BRIDGE_GATE'
    assert len(high['records'])==high['counts']['R_families']==10
    assert {r['R_index']for r in high['records']}==set(range(10))
    tc=Counter();maxima={};recipes=roots=0;left_endpoints=set()
    def interval(pair):
        lo,hi=map(F,pair);assert lo<=hi;return lo,hi
    def negative(pair):assert interval(pair)[1]<0
    def gap(proof):
        lo,hi=interval(proof['product'])
        if proof['kind']=='strict_obtuse':assert lo>0
        else:
            assert proof['kind']=='strict_consecutive_angle_gap'
            assert interval(proof['lower_threshold'])[1]<lo<=hi<interval(proof['upper_threshold'])[0]
    def leaf(proof):
        if proof.get('reason')=='assigned_cosine_gap':gap(proof);tc['assigned_gap_leaves']+=1;return
        if proof.get('reason')=='R_no_real_norm_branch':negative(proof['remainder']);return
        assert {b['sign']for b in proof['branches']}=={-1,1}and len(proof['branches'])==2
        for branch in proof['branches']:
            if 'native_obstruction'in branch:gap(branch['native_obstruction']);continue
            q=branch['proof'];reason=q.get('reason')
            if reason=='all_X_profiles_fail_projected_finite_AX_corner':
                assert interval(q['delta'])[0]>=0 and interval(q['K'])[0]>0 and F(q['w_lower']['relaxed_lower'])>0
                for field in('F_over_w','t_derivative_over_w','w_derivative_upper','A0_old44'):negative(q[field])
                upper=interval(q['maximum_AX_divided_by_sqrt_2q'])[1];assert upper<-1
                maxima['regular']=max(upper,maxima.get('regular',upper))
            elif reason=='X63_orientation_forces_X30_nonnegative':assert interval(q['N30'])[0]>=0
            elif reason=='negative_corner_residual':negative(q['d'])
            else:
                assert {r['sign']for r in q['branches']}=={-1,1}and len(q['branches'])==2
                for r in q['branches']:
                    if r['reason']=='strict_positive_native_quadratic':
                        a,c,b=map(interval,r['witness']['coefficients']);assert a[0]>0 and b[0]>=0 and c[0]>=0
                    else:
                        for pair in r['F_shifted_coefficients']+r['t_derivative_shifted_coefficients']:negative(pair)
                        for field in('s_derivative_leading','s_derivative_at_boundary'):negative(r[field])
                        upper=interval(r['maximum_AX'])[1];assert upper<-2
                        maxima['ordinary']=max(upper,maxima.get('ordinary',upper))
    def walk(node,depth=0):
        tc['nodes']+=1
        if 'children'in node:
            assert len(node['children'])==2 and 'split'in node and 'mid'in node;tc['splits']+=1
            a,b=[walk(c,depth+1)for c in node['children']];k=node['split'];mid=F(node['mid'])
            assert set(a)==set(b)and k in a and a[k][1]==b[k][0]==mid and a[k][0]<mid<b[k][1]
            assert all(a[j]==b[j]for j in a if j!=k)
            a[k]=(a[k][0],b[k][1]);return a
        tc['leaves']+=1;tc['max_depth']=max(tc['max_depth'],depth);leaf(node['proof'])
        return {k:interval(pair)for k,pair in node['box'].items()}
    labels=[]
    for row in high['records']:
        labels.append(row['R_recipe']['labels'])
        assert sum(len(c['X_recipes'])for c in row['cases'])==row['X_recipe_count'];recipes+=row['X_recipe_count']
        rp={'R'+str(j)for j in row['R_recipe']['high_positions']}
        for case in row['cases']:
            roots+=1;domain=walk(case['certificate'])
            expected=rp|{k for k,pos in(('X8',1),('XR',2))if any(a[pos]==6 for a in case['X_recipes'])}
            assert set(domain)==expected
            for lo,hi in domain.values():assert lo>0 and lo*lo<=3 and hi==2;left_endpoints.add(lo)
    assert labels==[[2,2,6],[2,6,2],[3,2,6],[4,2,6],[5,2,6],[6,2,2],[6,2,3],[6,2,4],[6,2,5],[6,2,6]]
    assert len(left_endpoints)==1 and roots==120 and recipes==150
    assert high['counts']==dict(nodes=18240,leaves=9180,splits=9060,assigned_gap_leaves=6211,max_depth=25,R_families=10,X_recipes=150)
    assert all(tc[k]==high['counts'][k]for k in('nodes','leaves','splits','assigned_gap_leaves','max_depth'))
    assert tc['nodes']==2*tc['splits']+roots and tc['leaves']==tc['splits']+roots
    return dict(construction_high_R_families=10,construction_high_R_X_recipes=recipes,construction_high_R_nodes=tc['nodes'],construction_high_R_leaves=tc['leaves'],construction_high_R_assigned_gap_leaves=tc['assigned_gap_leaves'],construction_high_R_maxdepth=tc['max_depth'],construction_saved_upper_bounds={k:str(v)for k,v in maxima.items()})


def small_checks():
    p,review,low,extensions,local,joint,raw,high=map(read,STEMS)
    a=theorem_checks(p,review,low,extensions);b=arithmetic_checks(local,joint,raw);c=construction_checks(high)
    return dict(**a,**b,**c,scope='Independent Fraction scalar, saved-mixture/mass/column composition, matched primary and independent classifications, and scoped construction partition-tree, saved-sign and count checks. No optimizer, historical census or whole geometry replay; no new arithmetic dimension or complete joint-domain claim.')


def run():
    prior=json.loads((ROOT/PREVIOUS).read_text());frozen={};preserved={}
    for section in('source_hashes','top_artifact_hashes','report_hashes'):
        for name,h in prior[section].items():
            assert digest(ROOT/name)==h,('frozen predecessor',name)
            assert name not in frozen or frozen[name]==h
            frozen[name]=h
        preserved[section]=len(prior[section])
    assert digest(ROOT/PREVIOUS_AUDITOR)==prior['auditor_sha256']
    frozen[PREVIOUS_AUDITOR]=prior['auditor_sha256'];frozen[PREVIOUS]=digest(ROOT/PREVIOUS)
    print('ARCHIVE_FROZEN_BYTES_PASS',len(frozen),flush=True)
    sources=dict(frozen);sources[HELPER]=digest(ROOT/HELPER)
    top=[ROOT/f'results/{s}.json'for s in STEMS]
    pending=list(top);seen=set();new_manifests=[];bibliography=list(prior['inherited_bibliographic_entries'])
    references=records_count=0;skipped=set();new_jsons=[]
    while pending:
        path=pending.pop();name=str(path.relative_to(ROOT))
        if name in seen:continue
        seen.add(name)
        if name in frozen:
            assert digest(path)==frozen[name];skipped.add(name);continue
        doc=json.loads(path.read_text());new_jsons.append(name)
        if not isinstance(doc,dict):continue
        bound=[]
        if 'database_sha256'in doc and name not in A.HISTORICAL_DATABASE_METADATA:
            assert isinstance(doc.get('database'),str);bound.append((doc['database'],doc['database_sha256']))
        mdocs=[doc]+([doc['prepared']]if isinstance(doc.get('prepared'),dict)else[])
        for md in mdocs:
            for field in('sources','source_hashes'):
                for source,value in md.get(field,{}).items():
                    if isinstance(value,dict):
                        assert field=='sources'and'title'in value
                        bibliography.append(dict(document=name,citation=source))
                        for key,h in value.items():
                            if key.endswith('_sha256'):
                                q=value.get(key[:-7])
                                if isinstance(q,str)and(ROOT/q).is_file():bound.append((q,h))
                    else:bound.append((source,value))
        if 'source_sha256'in doc:
            self_source='scripts/'+path.stem+'.py'
            assert(ROOT/self_source).is_file();bound.append((self_source,doc['source_sha256']))
        records=list(A.record_refs(doc));records_count+=len(records)
        for source,h in bound+records:
            assert re.fullmatch('[0-9a-f]{64}',h),(name,source)
            source=A.EXTERNAL_SOURCE_ALIASES.get(source,source);target=ROOT/source
            assert target.is_file()and digest(target)==h,(name,source)
            assert source not in sources or sources[source]==h
            sources[source]=h;references+=1
            if target.suffix=='.json':pending.append(target)
        if bound or records:new_manifests.append(dict(path=name,status=doc.get('status'),source_references=len(bound),record_references=len(records)))
    for path in top:sources[str(path.relative_to(ROOT))]=digest(path)
    print('ARCHIVE_NEW_GRAPH_PASS',len(new_manifests),'inherited_json_skips',len(skipped),flush=True)
    checks=small_checks();print('ARCHIVE_FRACTION_CHECKS_PASS',checks['scalar_profile_checks'],flush=True)
    reports=[ROOT/CHECKPOINT]+[ROOT/f'results/{s}.md'for s in STEMS]+[ROOT/p for p in EXTRA_REPORTS]
    links=[]
    for report in reports:
        assert report.is_file(),report
        for target in re.findall(r'\]\(([^)]+)\)',report.read_text()):
            if '://'in target or target.startswith('#'):continue
            resolved=(report.parent/target.split('#')[0]).resolve()
            assert resolved.is_file()or resolved==OUT,(report,target)
            links.append(dict(report=str(report.relative_to(ROOT)),target=target))
    syntax=[]
    for name in sorted(set(sources)-set(frozen)):
        if name.endswith('.py'):ast.parse((ROOT/name).read_text());syntax.append(name)
    ast.parse(Path(__file__).read_text());syntax.append(str(Path(__file__).relative_to(ROOT)))
    all_manifests={r['path']:r for r in prior['manifests']}
    for r in new_manifests:assert r['path']not in all_manifests;all_manifests[r['path']]=r
    return dict(status='PASS_INCREMENTAL_FROZEN_ARCHIVE_DISTANCE9_AND_RAW39',
                scope='Archive integrity and exact composition of distance9 dimension108, new-mass complete local packing covers and selected joint-port repairs, raw39 actual-corner constraints, and the scoped high-R continuous gate; not formalization, a new arithmetic dimension theorem or repeated historical verification.',
                checkpoint=CHECKPOINT,previous_archive=PREVIOUS,previous_archive_sha256=digest(ROOT/PREVIOUS),
                previous_frozen_sections_preserved=preserved,previous_auditor_preserved=True,
                frozen_json_graph_policy='Rehash all inherited files first; do not reparse their already verified dependency graph.',
                inherited_json_documents_not_reparsed=sorted(skipped),new_json_documents_parsed=sorted(new_jsons),
                historical_nonbinding_database_metadata=sorted(A.HISTORICAL_DATABASE_METADATA),
                recovered_historical_source_aliases=A.EXTERNAL_SOURCE_ALIASES,
                recovered_provenance_sha256=digest(ROOT/A.RECOVERED_PROVENANCE),
                inherited_bibliographic_entries=bibliography,small_certificate_checks=checks,
                report_count=len(reports),local_link_count=len(links),
                source_manifest_count=len(all_manifests),new_source_manifest_count=len(new_manifests),
                source_reference_count=prior['source_reference_count']+references,
                new_source_reference_count=references,
                nested_record_reference_count=prior['nested_record_reference_count']+records_count,
                unique_source_count=len(sources),manifests=sorted(all_manifests.values(),key=lambda r:r['path']),
                source_hashes=dict(sorted(sources.items())),
                source_hash_scope='Cumulative predecessor files and new recursive sources; entrypoint snapshots are historical.',
                top_artifact_hashes={str(p.relative_to(ROOT)):digest(p)for p in top},
                report_hashes={str(p.relative_to(ROOT)):digest(p)for p in reports},local_links=links,
                new_python_syntax_checks=syntax,
                entrypoint_hashes={p:digest(ROOT/p)for p in('README.md','RESEARCH_PLAN.md')},
                entrypoint_scope='Current entrypoints only; prior snapshots remain historical.',auditor_sha256=digest(Path(__file__)))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--small-only',action='store_true');args=ap.parse_args()
    if args.small_only:
        print(json.dumps(small_checks(),indent=2));raise SystemExit(0)
    d=run()
    if args.verify:assert d==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],{k:d[k]for k in('report_count','local_link_count','source_manifest_count','new_source_manifest_count','unique_source_count')},flush=True)
