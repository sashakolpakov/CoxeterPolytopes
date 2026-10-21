#!/usr/bin/env python3
"""Incremental frozen archive; saved crosschecks do not repeat geometric censuses."""
import argparse,ast,json,re
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from itertools import combinations
import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A
ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/spectator_spectrum_and_c3_frontier_20261017.md'
PREVIOUS='results/ideal5_local_and_raw23_model_frontier_20261017_audit.json'
PREVIOUS_AUDITOR='scripts/audit_ideal5_local_and_raw23_model_frontier_20261017.py'
OUT=ROOT/'results/spectator_spectrum_and_c3_frontier_20261017_audit.json'
HELPER='scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
STEMS=(
 'fresh_spectator_energy_below_one_20261017',
 'fresh_spectator_energy_one_boundary_20261017',
 'fresh_raw33_295_mixed_actual_corners_20261017',
 'fresh_raw33_295_mixed_full_model_20261017',
 'fresh_raw33_295_mixed_full_model_20261017_support',
 'fresh_weighted_actual_corner_capacities_20261017',
 'fresh_c3_companion_caps_20261017',
 'fresh_c3_companion_continuation_20261017',
 'fresh_c3_companion_region_closure_20261017',
 'fresh_distance8_adaptive_leaf_overlap_20261017',
 'fresh_distance8_high_overlap_20261017',
 'fresh_distance8_ideal5_small_factor_20261017',
 'fresh_finite8_ideal5_combined_payments_20261017',
)
REPORT_STEMS=STEMS
EXTRA_REPORTS=()
PENDING_INPUTS=()
COMPLETED_UNRESTRICTED_BOUND=77
@lru_cache(None)
def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())
def inverse_column_check(g,col,k):
    assert len(col)==len(g)
    assert all(sum((a*b for a,b in zip(row,col)),F())==int(i==k)for i,row in enumerate(g))
def bd_matrix(kind,n):
    g=[[F(2 if i==j else 0)for j in range(n)]for i in range(n)]
    if kind=='B':
        g[-1][-1]=4;edges=[(i,i+1,2 if i==n-2 else 1)for i in range(n-1)]
    else:edges=[(i,i+1,1)for i in range(n-3)]+[(n-3,n-2,1),(n-3,n-1,1)]
    for i,j,w in edges:g[i][j]=g[j][i]=-F(w)
    return g

def spectrum_checks():
    d=read('fresh_spectator_energy_below_one_20261017')
    assert d['counts']==dict(A_ranks=64,A_adjacent_pairs=2016,BD_templates=124)
    assert len(d['A_identity_checks'])==64 and len(d['BD_obstruction_checks'])==124
    pairs=0
    for r in d['A_identity_checks']:
        n=r['rank'];assert F(r['energy'])==F(n,n+1)
        assert list(map(F,r['inverse_leaf_column']))==[F(n-i,n+1)for i in range(n)]
        assert [x['marks']for x in r['adjacent_contact_pairs']]==[[i,i+1]for i in range(n-1)]
        for p in r['adjacent_contact_pairs']:
            assert p['separate_edge_extensions_PD']==all(F((j+1)*(n-j),n+1)<2 for j in p['marks'])
        pairs+=len(r['adjacent_contact_pairs'])
    assert pairs==2016
    for r in d['BD_obstruction_checks']:
        g=bd_matrix(r['kind'],r['rank']);v=list(map(F,r['first_inverse_column']))
        inverse_column_check(g,v,0);assert v[0]==1 and all(t>0 for t in v)
        if r['second_inverse_column']is not None:
            w=list(map(F,r['second_inverse_column']));inverse_column_check(g,w,1);assert w[1]==2
        assert all(F(t)>0 for t in r['positive_pivots'])
    assert F(d['nonunit_A2_excess_positive_square_difference'])==1-2*F(2,3)**2>0
    for r in d['rational_energy_applications']:
        e=F(r['rational_energy']);n=e/(1-e)
        assert n==F(r['forced_rank']) and r['allowed']==(n.denominator==1 and n>=1)
    b=read('fresh_spectator_energy_one_boundary_20261017')
    assert b['counts']==dict(one_contact_BD_identities=125,twin_identities=62,isolated_identities=3)
    for r in b['one_nonempty_support']:
        g=bd_matrix('D'if r['kind']=='A3_middle'else r['kind'],r['rank'])
        v=list(map(F,r['positive_inverse_contact_column']));inverse_column_check(g,v,0)
        assert v[0]==F(r['energy'])==1 and F(r['strict_extension_schur'])==1 and all(t>0 for t in v)
    for r in b['opposite_twin_contacts']:
        n=r['rank'];assert F(r['separate_projection_norm'])==F(n,4)
        assert F(r['separate_extension_schur'])==2-F(n,4)
        assert r['strictly_allowed']==(n<=7)
        z=[F(0)]*(n-2)+[F(1),F(-1)];v=list(map(F,r['inverse_difference']))
        assert v==[t/2 for t in z]
        assert all(sum((a*c for a,c in zip(row,v)),F())==z[i]for i,row in enumerate(bd_matrix('D',n)))
        assert sum((a*c for a,c in zip(z,v)),F())==F(r['energy'])==1
    for r in b['isolated_cases']:
        diag,x,y=map(lambda key:list(map(F,r[key])),('gram','x','y'))
        xx=sum((t*t/g for t,g in zip(x,diag)),F());yy=sum((t*t/g for t,g in zip(y,diag)),F())
        assert [xx,yy]==list(map(F,r['separate_projection_norms']))and xx<2 and yy<2
        assert sum(((a-b)**2/g for a,b,g in zip(x,y,diag)),F())==1
    return dict(A_ranks=64,A_adjacent_pairs=2016,BD_columns=124,boundary_one_contact=125,
      boundary_twins=62,boundary_isolated=3,scope='Independent Fraction substitution in saved identities and rational consequences. All-rank completeness remains the separately reviewed written proof.')


def cone(record):
    bd,ic=record['boundary'],record['PD_incidence']
    assert not bd['warnings'] and ic['status']=='numerically_separated_acute_incidence'
    finite={tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}
    assert len(finite)==bd['counts']['finite']==ic['finite']
    assert finite=={tuple(sorted(r))for r in ic['finite_supports']}
    assert len(bd['rays'])==sum(bd['counts'].values()) and not bd['counts'].get('ideal',0)
    active_edges={tuple(sorted(e['sides']))for e in bd['edges']if e['kind']!='outside'}
    assert active_edges=={tuple(sorted(e))for e in ic['edge_supports']}
    assert len(active_edges)==ic['geodesic_edges']
    assert bd['open_ends']==ic['open_ends']==2*len(active_edges)-3*len(finite)
    return bd

def construction_checks():
    first=read('fresh_c3_companion_caps_20261017');cont=read('fresh_c3_companion_continuation_20261017');end=read('fresh_c3_companion_region_closure_20261017')
    for d in(first,cont,end):assert d['retained_control_walls']==63 and d['deleted_walls']==[]and not d['unresolved']
    assert first['old_finite_vertices_preserved_by_exact_C3']==199
    bd=cone(first['baseline']);assert bd['counts']==dict(finite=27,exterior=7)and bd['open_ends']==11
    assert first['counts']==dict(norm_equations=145,norm_branches=274,coherent_companions=54,companion_pairs_tested=1431)
    assert len(first['companions'])==54 and not first['coherent_pairs']
    assert [len(d['recipes'])for d in first['finite_diagram_domains']]==[24,51,40]
    assert [sum(bool(r['high_positions'])for r in d['recipes'])for d in first['finite_diagram_domains']]==[3,5,2]
    profiles=Counter();triangles=0
    for r in first['companions']:
        bd=cone(r);assert bd['counts']==dict(finite=29,exterior=7)and bd['open_ends']==11
        profiles[r['profile']]+=1
        assert r['center_product']<0 and abs(r['norm_residual'])<1e-7
        if r['frame']==[7,30,67]:
            rays=[v for v in bd['rays']if 67 in v['active']]
            assert len(rays)==3 and all(v['kind']=='finite'for v in rays);triangles+=1
    assert profiles=={0:47,1:3,2:1,4:3}and triangles==39
    assert cont['counts']==dict(norm_equations=2837,norm_branches=5674,coherent_continuations=552,continuation_pairs_tested=5196,coherent_continuation_pairs=100)
    assert cont['parent_companions']==54 and cont['C3_retained']and not cont['improvements']
    hist=Counter();pairs=0;prefixes={}
    for r in cont['records']:
        p=first['companions'][r['parent']];assert r['parent_profile']==p['profile']and r['parent_labels']==p['labels']
        for c in r['continuations']:
            bd=cone(c);hist[bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
            assert c['center_product']<0 and abs(c['norm_residual'])<1e-7
        for c in r['coherent_pairs']:
            assert c['relation']=='dotted'and c['native_product']<-2
            bd=cone(c);assert bd['counts']==dict(finite=33,exterior=7)and bd['open_ends']==11
            prefixes[r['parent'],tuple(c['indices'])]=c;pairs+=1
    assert hist=={(31,7,11):551,(31,7,13):1}and pairs==100
    assert end['counts']==dict(norm_equations=1010,norm_branches=2020,native_candidates=340)
    assert end['prefixes']==len(end['records'])==100 and end['boundary_cycle']==[30,52,70,68,69]
    assert not end['hits']and not end['improvements']
    candidates=recipes=0;maxA=maxL=float('-inf')
    for r in end['records']:
        assert(r['parent'],tuple(r['pair']))in prefixes
        assert all(0<=q[0]<7 and all(2<=m<=5 for m in q[1:])for q in r['complete_finite_domain'])
        recipes+=len(r['complete_finite_domain'])
        for c in r['native_candidates']:
            bd=cone(c);assert bd['counts']==dict(finite=35,exterior=7)and bd['open_ends']==11
            assert c['five_finite_corners']==[True,True,False,False,False]
            assert c['center_product']<0
            native=dict(c['native_products']);maxA=max(maxA,native[68]);maxL=max(maxL,native[69]);candidates+=1
    assert candidates==340 and recipes==1010 and maxA<-4.77 and maxL<-30.17
    return dict(first_companions=54,whole_C_triangles=39,continuations=552,coherent_pairs=100,
      final_native_candidates=340,final_required_corner_pattern=[True,True,False,False,False],
      best_scored_open_ends=11,exact_control_open_ends=10,old_wall_deletions=0,
      scope='Saved numerical complete-cone/PD support and graph-count comparisons only. No numerical recomputation, interval re-enclosure, exact new chamber or finite-volume claim.')


def raw33_checks():
    d=read('fresh_raw33_295_mixed_actual_corners_20261017');counts=d['counts']
    assert(counts['actual_Q_classes'],counts['corner_occurrences'],counts['literal_orientations'],counts['forests'])==(243,2024,35299,527)
    assert counts['columns']==len(d['columns'])==352 and counts['violated_columns']==119 and counts['vertex_mass_increases']==72
    assert counts['mixed_columns']==208 and counts['violated_mixed_columns']==93 and counts['mixed_additional_vertex_raises']==55
    old=read('fresh_current220_actual_corner_union_20261017')
    assert d['face_multipliers']==old['face_multipliers']and d['initial_vertex_multipliers']==old['vertex_multipliers']
    del old
    faces={r['face']:r for r in d['selected_faces']};origins={}
    for r in faces.values():
        for a in r['corners']:
            key=(a['face'],tuple(a['omitted']),a['vertex']);assert key not in origins;origins[key]=a
    assert len(origins)==2024
    total=reused=excluded=0;seen_cols=set()
    for row in d['forest_records']:
        previous=[];owned_total=owned_excluded=0
        for dom in row['inherited_pair_domains']:
            mp=dict(dom['new_to_old_star_indices']);assert len(mp)==len(set(mp.values()))
            cur=set(mp);owned=len(cur)*(len(cur)-1)//2
            for k in range(1,len(previous)+1):
                for sub in combinations(previous,k):
                    m=len(cur.intersection(*sub));owned+=(-1)**k*m*(m-1)//2
            assert owned==dom['pairs_reused']and 0<=dom['exclusions_reused']<=owned
            previous.append(cur);owned_total+=owned;owned_excluded+=dom['exclusions_reused']
        n=len(row['literal_stars']);pairs=n*(n-1)//2
        assert pairs==row['tested_literal_pairs']==owned_total+row['newly_tested_pairs']
        assert pairs-row['unknown_literal_pairs']==owned_excluded+len(row['incompatibility_certificates'])
        groups=row['corner_orbits'];assert sorted(i for g in groups for i in g['literal_stars'])==list(range(n))
        for g in groups:
            coeff=Counter()
            for a in g['origins']:
                assert origins[a['face'],tuple(a['omitted']),a['vertex']]==a;coeff[a['face']]+=1
            assert dict(coeff)==dict(g['positive_face_coefficients'])
        for ci in row['columns']:
            assert ci not in seen_cols;seen_cols.add(ci);c=d['columns'][ci]
            assert c['vertex']==row['vertex']
            indices=c['orbit_clique']if c['kind']=='capacity_one_orbit_clique'else c['orbit_indices']
            coeff=Counter()
            for j in indices:coeff.update(dict(groups[j]['positive_face_coefficients']))
            assert sorted(coeff.items())==list(map(tuple,c['positive_face_coefficients']))
            if c['kind']=='capacity_one_orbit_clique':assert c['vertex_cost']==1
            else:
                assert c['vertex_cost']==len(c['literal_coloring'])
                assert sorted(i for part in c['literal_coloring']for i in part)==sorted(i for j in indices for i in groups[j]['literal_stars'])
        total+=pairs;reused+=owned_total;excluded+=owned_excluded
    assert seen_cols==set(range(352))and total==counts['literal_pairs_tested']==4988983
    assert reused==counts['pairs_reused']==4785958 and excluded==counts['exclusions_reused']==4686927
    alpha={i:F(a)for i,a in d['face_multipliers']};initial={i:F(a)for i,a in d['initial_vertex_multipliers']}
    expected=dict(initial);single=dict(initial);violations=mixed=0
    for c in d['columns']:
        load=sum((F(m)*alpha[i]for i,m in c['positive_face_coefficients']),F());v=c['vertex'];cost=F(c['vertex_cost'])
        assert load==F(c['load']) and cost>0
        violation=load>cost*initial.get(v,F());assert violation==c['violated'];violations+=violation
        expected[v]=max(expected.get(v,F()),load/cost)
        if c['cross_cohort']:mixed+=violation
        else:single[v]=max(single.get(v,F()),load/cost)
    final={i:F(a)for i,a in d['vertex_multipliers']};cert=d['fixed_witness_repair']
    assert expected==final and violations==119 and mixed==93
    assert sum(expected.values())==F(cert['normalizer'])
    objective=F(cert['objective_numerator']);full=objective/sum(final.values());without=objective/sum(single.values())
    assert full==F(cert['bound'])<12<without==F(cert['without_cross_cohort_bound'])
    assert sum(final.get(v,F())>initial.get(v,F())for v in final)==72
    return dict(classes=243,corners=2024,columns=352,full_ratio=str(full),without_cross_ratio=str(without),
      saved_pairs=total,first_owned_pairs=reused,scope='Saved first-owner inclusion-exclusion, occurrence coefficients, coloring partitions and exact loads. No repeated wall-pair tests or independent graph-clique arithmetic.')


def model_support_checks():
    d=read('fresh_raw33_295_mixed_full_model_20261017');p=d['prepared'];n=p['counts']
    assert d['solver_invocations']==1 and d['solver']['status']==0 and d['solver']['success']
    assert(n['inherited_variables'],n['variables'],n['faces'],n['vertex_forests'],n['zero_budget_signed_variables'],n['vertex_cost_signed_variables'])==(438266,439061,82806,6558,96176,52)
    assert len(p['new_columns'])==n['distinct_added_columns']==795
    assert {x['variable']for x in p['new_columns']}==set(range(438266,439061))
    q=next(c for c in d['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    a={i:F(v)for i,v in q['face_multipliers']};m={i:F(v)for i,v in q['vertex_multipliers']}
    assert len(a)==len(q['face_multipliers'])==4713 and len(m)==len(q['vertex_multipliers'])==986
    assert all(0<=i<82806 and v>0 for i,v in a.items())and all(0<=i<6558 and v>0 for i,v in m.items())
    assert sum(m.values())==F(q['normalizer'])and F(q['objective_numerator'])/sum(m.values())==F(q['bound'])>12
    assert q['excludes_target12']
    for c in p['new_columns']:
        assert c['positive_denominator']>0 and c['scaled_vertex_cost']==F(c['vertex_cost'])*c['positive_denominator']
        assert sum((F(v)*a.get(i,F())for i,v in c['positive_face_coefficients']),F())<=F(c['vertex_cost'])*m.get(c['vertex'],F())
    primal=next(c for c in d['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    assert primal['integer_denominator']>0 and F(primal['repair_delta'])>=0
    assert len({i for i,v in primal['rounded_variable_numerators']})==len(primal['rounded_variable_numerators'])
    assert all(0<=i<439061 and isinstance(v,int)and v>=0 for i,v in primal['rounded_variable_numerators'])
    budgets={i:F(v)for i,v in primal['vertex_budgets']}
    assert len(budgets)==len(primal['vertex_budgets'])==6557 and set(budgets)<=set(range(6558))
    assert all(v>=0 for v in budgets.values())
    maximum=max(budgets.values());assert maximum==F(primal['maximum'])and primal['meets_target12']==(maximum<=12)
    assert F(q['bound'])<=maximum and maximum==F(6191130293139,500000000000)
    support=read('fresh_raw33_295_mixed_full_model_20261017_support')
    assert support['face_multipliers']==q['face_multipliers']and support['vertex_multipliers']==q['vertex_multipliers']
    assert support['optimizer_invocations']==0
    qr={r['face']:r for r in support['quad_rows']};Q=set(qr);S=set(a)-Q
    assert len(Q)==1993 and len(S)==2720
    dom=support['explicit_geometry_domains'];old=set(dom['current220_rows']);new=set(dom['raw33_295_rows'])
    assert len(old)==220 and len(new)==243 and old<=new and len(Q&new)==184
    assert Counter(r['provenance_kind']for r in qr.values())==dict(explicit=377,unique_saved_context=1214,merged_unexpanded=402)
    for i,r in qr.items():assert F(r['mass'])==a[i]and F(r['demand_mass'])==a[i]*r['required']
    for block,tag,ncols in zip(support['appended_export_loads'],('fresh_current220_actual_corner_union_20261017','fresh_raw33_295_mixed_actual_corners_20261017'),(443,352)):
        export=read(tag);assert block['source']==tag and len(block['loads'])==len(export['columns'])==ncols
        for c,r in zip(export['columns'],block['loads']):
            load=sum((F(v)*a.get(i,F())for i,v in c['positive_face_coefficients']),F());capacity=F(c['vertex_cost'])*m.get(c['vertex'],F())
            assert r['index']==c['index']and r['vertex']==c['vertex']
            assert F(r['load'])==load<=capacity==F(r['capacity'])and F(r['slack'])==capacity-load
        assert block['positive_loads']==sum(F(r['load'])>0 for r in block['loads'])
        assert block['exactly_tight_positive']==sum(F(r['load'])>0 and F(r['slack'])==0 for r in block['loads'])
        del export
    ratios=[]
    for r,keep in zip(support['exact_restrictions'],(Q&new,{i for i,v in qr.items()if v['provenance_kind']!='merged_unexpanded'})):
        cert=r['certificate'];aa={i:F(v)for i,v in cert['face_multipliers']};mm={i:F(v)for i,v in cert['vertex_multipliers']}
        assert aa=={i:a[i]for i in S|keep}and all(0<v<=m[i]for i,v in mm.items())
        assert r['kept_S']==2720 and r['kept_Q']==len(keep)and r['dropped_Q']==len(Q-keep)
        assert r['alpha_unchanged_on_retained_rows']and not r['outside_S_restored']
        objective=F(q['objective_numerator'])-sum((a[i]*qr[i]['required']for i in Q-keep),F())
        assert objective==F(cert['objective_numerator'])and sum(mm.values())==F(cert['normalizer'])
        ratio=objective/sum(mm.values());assert ratio==F(cert['bound'])>12 and cert['excludes_target12'];ratios.append(str(ratio))
    return dict(variables=439061,faces=82806,vertex_budgets=6558,appended_loads=795,signed_cost_columns_preserved=52,
      finite_model_dual=q['bound'],finite_model_primal=str(maximum),restriction_ratios=ratios,
      scope='Saved append loads, budget maxima and restriction accounting only; no LP rerun, full matrix rebuild or new-alpha implicit packing coverage.')


def marked_key(r):return(r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['match'])))

def profile_margin(ranks):
    h=[1]
    for r in ranks:
        out=[0]*(len(h)+r)
        for i,v in enumerate(h):
            for j in range(r+1):out[i+j]+=v
        h=out
    d=sum(ranks);V=sum(h);G=sum(v*min(i,d-i)for i,v in enumerate(h))
    def cost(r):return F(r,2)if r<=12 else{13:F(87,14),14:F(20,3),15:F(111,16),16:F(127,17)}.get(r,F(163,23))
    return h,V,G,F(5*G-d*V,2*V)-sum(map(cost,ranks),F())

def finite8_checks():
    adaptive=read('fresh_distance8_adaptive_leaf_overlap_20261017')
    assert len(adaptive['cases'])==9 and adaptive['method_counts']==dict(adaptive_untouched_marked_leaf=3,sealed_marked_leaf=4,adaptive_leaf_or_bounded_relevant_tail=2)
    assert Counter(r['method']for r in adaptive['cases'])==adaptive['method_counts']
    corners=0
    for r in adaptive['cases']:
        g=[list(map(F,row))for row in r['gram']];J=[list(map(F,row))for row in r['inverse']]
        for k in range(len(g)):inverse_column_check(g,[row[k]for row in J],k)
        corners+=len(r['all_corners'])
    assert corners==36
    ov=read('fresh_distance8_high_overlap_20261017');c=ov['census'];iso=ov['isolation']
    assert(c['marked_contexts'],c['candidate_matchings'],c['consistent_matchings'])==(384,6390952,82197)
    assert len(c['overlap_records'])==82197 and len(iso['cores'])==154 and len(iso['cases'])==2614
    assert len({marked_key(r)for r in c['overlap_records']})==82197
    kinds=Counter(r['kind']for r in c['overlap_records'])
    assert kinds['inverse_rejection']==75947 and kinds['edge_rejection']==3636
    assert sum(kinds.values())-75947-3636==2614
    assert kinds==c['overlap_counts']
    assert Counter(r['method']for r in iso['cases'])==iso['method_counts']
    assert sum(iso['method_counts'].values())==2614
    assert {marked_key(r)for r in iso['cases']}=={marked_key(r)for r in c['overlap_records']if r['kind']not in('inverse_rejection','edge_rejection')}
    assert len({r['core_key']for r in iso['cases']})==154
    adaptive_keys={marked_key(r)for r in adaptive['cases']}
    assert adaptive_keys=={marked_key(r)for r in iso['cases']if r['method']=='nine_context_adaptive_finite_payment'}
    shortcuts=witnesses=0
    cases={r['index']:r for r in iso['cases']}
    assert set(cases)==set(range(2614))
    for r in iso['cases']:
        assert len(r['all_corners'])==4 and r['core_key']in iso['cores']
        shortcuts+=r['extra_cross_shortcuts']
        for w in r['extra_cross_pd_rejections']:
            assert F(w['square'])<=0 and len(w['nodes'])==len(w['vector'])
            v=list(map(F,w['vector']));assert all(t>=0 for t in v)and any(v)
            witnesses+=1
    assert shortcuts==iso['extra_cross_shortcuts']==154502 and witnesses==iso['extra_cross_pd_rejections']==5528
    small=read('fresh_distance8_ideal5_small_factor_20261017');oldsmall=read('fresh_ideal5_small_factor_payments_20261017')
    patterns={r['index']:r for r in small['patterns']}
    assert len(patterns)==len(small['patterns'])==1153 and len(small['records'])==16
    assert Counter(str(r['order_bound'])for r in patterns.values())==small['order_bound_counts']=={'14':313,'15':354,'16':346,'17':140}
    assert set(patterns)=={i for i,r in cases.items()if r.get('component_order_bound',0)>=14}
    for i,r in patterns.items():
        parent=cases[i];assert marked_key(r)==marked_key(parent)
        assert r['nodes']==parent['marked_component']and r['method']==parent['method']and r['order_bound']==parent['component_order_bound']
        assert [r['nodes'][j]for j in r['marks']]==parent['preferred_corner']
    embeddings=0;sealed={}
    for row,old in zip(small['records'],oldsmall['records']):
        assert(row['order'],row['kind'])==(old['order'],old['kind'])
        assert row['target']==old['target']and row['retained_paid_pairs']==old['paid_pairs']
        paid=set(map(tuple,row['retained_paid_pairs']));drop=set(map(tuple,row['dropped_pairs']))
        assert len(paid)==row['retained_cost']==old['cost']and drop==set(map(tuple,old['dropped_pairs']))
        union=set();affected=set();seen=set()
        for rec in row['case_embeddings']:
            i=rec['case'];assert i not in seen;seen.add(i);pat=patterns[i]
            assert row['order']<=pat['order_bound']and len(pat['nodes'])<=row['order']
            images=set();maps=set()
            for mp in rec['embeddings']:
                assert len(mp)==len(pat['nodes'])==len(set(mp))and all(0<=v<row['order']for v in mp)
                assert tuple(mp)not in maps;maps.add(tuple(mp))
                assert all(pat['labels'][a][b]==row['target'][mp[a]][mp[b]]for a in range(len(mp))for b in range(len(mp)))
                images.add(tuple(sorted(mp[j]for j in pat['marks'])));embeddings+=1
            assert images==set(map(tuple,rec['marked_images']));union|=images
            if images&drop:affected.add(i)
        assert union==set(map(tuple,row['preferred_pairs']))and union<=paid|drop
        assert affected==set(row['dropped_affected_cases'])
        corr=row['old_sealed_case_correspondence']
        if row['kind']=='D'and row['order']in(16,17):
            assert len(corr)==9 and {r['new_case']for r in corr}==affected
            assert sorted(r['old_case']for r in corr)==old['dropped_affected_cases']
            for r in corr:assert marked_key(r)==marked_key(patterns[r['new_case']])
            sealed[row['order']]=9
        else:assert not corr
    assert sealed=={16:9,17:9}
    binder=read('fresh_finite8_ideal5_combined_payments_20261017');oldbinder=read('fresh_ideal5_combined_payments_20261017')
    assert binder['finite_vertex_budget']=='8n'and binder['finite_budget_dimension_maximum']==80
    assert F(binder['triangle_minimum'])==2 and F(binder['quadrilateral_minimum'])==1
    assert(binder['high_overlap_contexts'],binder['low_overlap_contexts'],binder['low_component_order_bound'],binder['selected_small_patterns'])==(2614,66,12,1153)
    assert binder['full_new_overlap_domain']and binder['old_order16_17_affected_keys_matched']
    assert binder['ideal_rank_envelope']==oldbinder['ideal_rank_envelope']and binder['small_order_costs']==oldbinder['small_order_costs']
    assert len(binder['triangle_pair_retention'])==len(binder['short8_reduction'])==16
    for tr,short,rec in zip(binder['triangle_pair_retention'],binder['short8_reduction'],small['records']):
        assert(tr['order'],tr['kind'])==(short['order'],short['kind'])==(rec['order'],rec['kind'])
        target=rec['target'];n=len(target);ds=[[0 if i==j else 1 if target[i][j]else n+1 for j in range(n)]for i in range(n)]
        features={i for i,row in enumerate(target)if sum(bool(t)for t in row)>=3 or any(t>=4 for t in row)}
        for k in range(n):
            for i in range(n):
                for j in range(n):ds[i][j]=min(ds[i][j],ds[i][k]+ds[k][j])
        near={i for i in range(n)if any(ds[i][f]<=5 for f in features)}
        required={p for p in combinations(range(n),2)if ds[p[0]][p[1]]<=5 or(ds[p[0]][p[1]]<=8 and set(p)&near)}
        paid=set(map(tuple,rec['retained_paid_pairs']))
        assert tr['retained']and required==set(map(tuple,tr['required_pairs']))and required<=paid
        missing={p for p in combinations(range(n),2)if ds[p[0]][p[1]]<=8 and p not in paid}
        assert missing==set(map(tuple,short['unpaid_short8_pairs']))and(not missing or rec['kind']=='A')
    frontier=binder['fixed_envelope_frontier']
    assert [r['ambient_dimension']for r in frontier['diagnostic_minima']]==list(range(70,78))
    for r in frontier['diagnostic_minima']:
        n=r['ambient_dimension'];assert sum(r['minimizing_ranks'])==n-1
        assert profile_margin(r['minimizing_ranks'])[3]==F(r['minimum_margin'])<0
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta-8*n==F(r['finite8_slack'])>=0
    boundary=frontier['boundary_at_n77'];h,V,G,margin=profile_margin([12,16,16,16,16])
    assert h==boundary['h_polynomial']and V==boundary['link_vertices']==1085773 and G==boundary['G']
    assert margin==F(boundary['margin'])==F(-16664,1085773)
    assert -margin*V==frontier['unnormalized_cusp_deficit']==16664
    assert F(frontier['finite8_slack_at_n77'])==68 and frontier['strict_global_shape_constant_at_n77']==380
    assert F(frontier['required_normalized_single_rank12_cost_decrease'])==-margin
    assert F(frontier['required_each_rank16_normalized_decrease_if_shared_equally'])==-margin/4
    assert F(frontier['required_order13_pair_cost_decrease'])==-13*margin
    assert F(frontier['required_order17_pair_cost_decrease_each'])==-17*margin/4
    return dict(adaptive_cores=9,adaptive_corners=36,consistent_overlap_records=82197,overlap_contexts=2614,
      overlap_cores=154,extra_cross_witness_signs=5528,preferred_patterns=1153,saved_induced_maps_checked=embeddings,
      ideal_targets=16,sealed_replacements=18,triangle_minimum='2',quadrilateral_minimum='1',finite_budget='8n',
      ideal_envelope_unchanged=True,boundary77_margin=str(margin),dimension_improvement=False,
      scope='Saved census partitions, witness signs, literal induced maps, paid-set identities and selected scalar profiles. No candidate matching/extension census, proof of embedding exhaustiveness, or new global bound is reconstructed by this archive.')


def weighted_checks():
    d=read('fresh_weighted_actual_corner_capacities_20261017');support=read('fresh_raw33_295_mixed_full_model_20261017_support')
    assert not d['open_forests']and len(d['forest_records'])==486
    assert d['face_multipliers']==support['face_multipliers']
    alpha={i:F(a)for i,a in d['face_multipliers']};baseline={v:F(a)for v,a in support['vertex_multipliers']}
    restricted=support['exact_restrictions'][0]['certificate'];rbase={v:F(a)for v,a in restricted['vertex_multipliers']}
    assert d['restricted_face_multipliers']==restricted['face_multipliers']
    new=dict(baseline);rnew=dict(rbase);orbits=nodes=vectors=terms=cliques=0;rows={}
    for row in d['forest_records']:
        v=row['vertex'];assert v not in rows;rows[v]=row
        assert F(row['baseline_full_mu'])==baseline.get(v,F())and F(row['baseline_restricted_mu'])==rbase.get(v,F())
        demand=[sum((F(m)*alpha.get(i,F())for i,m in cc),F())for cc in row['orbit_coefficients']]
        assert demand==list(map(F,row['demands']))and len(demand)==len(row['orbit_indices'])<=14
        assert all(x>0 for x in demand)
        cert=row['certificate'];weights=list(map(F,cert['weights']));assert len(weights)==len(demand)and all(x>=0 for x in weights)
        vs=row['maximal_count_vectors'];assert len(vs)==len(row['literal_representatives'])==row['distinct_count_vectors']
        for vec in vs:
            assert len(vec)==len(demand)and all(isinstance(x,int)and x>=0 for x in vec)
            assert sum((w*c for w,c in zip(weights,vec)),F())<=1
        cover=[F()for _ in demand];mass=F()
        for j,value in cert['mixture']:
            value=F(value);assert value>0 and 0<=j<len(vs)
            literal=row['literal_representatives'][j]
            assert len(literal)==len(set(literal))==sum(vs[j])and set(literal)<=set(row['literal_nodes'])
            for k in range(len(demand)):cover[k]+=value*vs[j][k]
            mass+=value;terms+=1
        assert all(x>=y for x,y in zip(cover,demand))
        assert mass==F(cert['optimal_mass'])==sum((x*y for x,y in zip(weights,demand)),F())
        new[v]=max(new.get(v,F()),mass);rnew[v]=max(rnew.get(v,F()),mass)
        orbits+=len(demand);nodes+=len(row['literal_nodes']);vectors+=len(vs);cliques+=row['maximal_literal_cliques']
    assert(orbits,nodes,vectors,terms,cliques)==(1582,29459,1477,1299,95989)
    assert d['counts']==dict(distinct_count_vectors=1477,forests=486,full_violations=84,
      literal_stars=29459,maximal_count_vectors=1477,maximal_literal_cliques=95989,
      mixture_terms=1299,orbits=1582,restricted_violations=86)
    ratios=[]
    for repair,mu,base in zip(d['repairs'],(new,rnew),(baseline,rbase)):
        assert {v:F(x)for v,x in repair['vertex_multipliers']}=={v:x for v,x in mu.items()if x}
        assert F(repair['normalizer'])==sum(mu.values())
        ratio=F(repair['objective_numerator'])/sum(mu.values());assert ratio==F(repair['bound']);ratios.append(ratio)
        assert repair['vertex_updates']==[[v,str(base.get(v,F())),str(x)]for v,x in sorted(mu.items())if x>base.get(v,F())]
    assert len(d['repairs'][0]['vertex_updates'])==84 and len(d['repairs'][1]['vertex_updates'])==86
    assert ratios[0]>12>ratios[1]and len(d['columns'])==86
    for col in d['columns']:
        row=rows[col['vertex']];coef=Counter()
        for w,cc in zip(row['certificate']['weights'],row['orbit_coefficients']):
            for i,m in cc:coef[i]+=F(w)*m
        assert col['positive_face_coefficients']==[[i,str(x)]for i,x in sorted(coef.items())if x]
        assert F(col['vertex_cost'])==1 and F(col['load'])==sum((x*alpha.get(i,F())for i,x in coef.items()),F())
        assert col['violates_full']==(F(col['load'])>baseline.get(col['vertex'],F()))
        assert col['violates_restricted']==(F(col['load'])>rbase.get(col['vertex'],F()))
    return dict(forests=486,orbits=orbits,literal_stars=nodes,saved_maximal_cliques=cliques,count_vectors=vectors,
      exact_mixture_terms=terms,columns=86,full_ratio=str(ratios[0]),restricted_ratio=str(ratios[1]),
      scope='Saved nonnegative charges, vector upper bounds, equal-mass covering mixtures, occurrence coefficients and monotone repairs. Separate no-producer checker establishes literal graph compatibility/charge maxima; this archive does not enumerate cliques or assert physical realizability.')

def small_checks():
    return dict(completed_unrestricted_bound=77,dimension_improvement=False,
      spectator_spectrum=spectrum_checks(),construction=construction_checks(),raw33_geometry=raw33_checks(),
      model_support=model_support_checks(),finite8_local=finite8_checks(),weighted_capacities=weighted_checks(),
      review_scopes=dict(C3='Root full source/math reads of all three numerical packages; no root/cone recomputation.',
        raw33_geometry='Root review plus frontier full source/math review after producer report freeze; no repeated pair census.',
        full_model_and_support='Root narrow source/math reads; independent saved append-load and restriction accounting.',
        spectrum_below_one='Exact_hub full source/math read, no independent finite identity replay.',
        spectrum_boundary='Construction_audit full source/math read, no duplicate 125/62 producer run.',
        finite8_binder='Root full composition/source read; exact_hub independently read final source/report after freeze, without upstream census/scalar replay.',
        weighted_capacities='Root full source/report read and saved Fraction accounting; separate no-producer integer clique and literal-mixture checker, no native pair census.'),
      finite_volume_example=False,scope='Incremental preservation and saved certificate crosschecks. Exact local coverage does not itself lower77; numerical C3 extensions do not constitute an exact new chamber.')

def run(allow_pending=False):
    assert allow_pending or not PENDING_INPUTS, ('archive inputs not final', PENDING_INPUTS)
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
    checks=small_checks();print('ARCHIVE_FRACTION_CHECKS_PASS',len(STEMS),flush=True)
    reports=[ROOT/CHECKPOINT]+[ROOT/f'results/{s}.md'for s in REPORT_STEMS]+[ROOT/p for p in EXTRA_REPORTS]
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
    return dict(status='PASS_INCREMENTAL_SPECTATOR_SPECTRUM_AND_C3_ARCHIVE',
                scope='Preserves the completed unrestricted77 proof and binds local spectator/finite8 coverage, raw33/295 arithmetic capacities and finite-model support, and numerical deletion-free C3 construction diagnostics. No automatic dimension improvement, optimized full-coverage feasibility or new finite-volume example. Saved crosschecks do not replace geometric proofs, contact arithmetic or matrix certificate replays.',
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
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--small-only',action='store_true');ap.add_argument('--inventory',action='store_true');args=ap.parse_args()
    if args.small_only:
        print(json.dumps(small_checks(),indent=2));raise SystemExit(0)
    d=run(allow_pending=args.inventory)
    if args.inventory:
        print('READ_ONLY_ARCHIVE_INVENTORY',{k:d[k]for k in('report_count','local_link_count','source_manifest_count','new_source_manifest_count','source_reference_count','unique_source_count')},flush=True)
        raise SystemExit(0)
    if args.verify:assert d==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],{k:d[k]for k in('report_count','local_link_count','source_manifest_count','new_source_manifest_count','unique_source_count')},flush=True)
