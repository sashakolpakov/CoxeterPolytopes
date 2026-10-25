#!/usr/bin/env python3
"""Pure transports of existing old-affine source proofs and CURRENT receipts.

The enclosing final consumer supplies exact byte readers, the unchanged
3a1d legacy helper with its original HIGH projection, current W and ONLY
F.actual-validated units. No worker is imported or prepared here. Missing
source applicability remains a refusal, not an intrinsic-unit promotion.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
import re

SOURCE_PINNER=None
SOURCE_READER=None
SOURCE_ACCEPTANCE=None
APPLICATION_SOURCE_SHA=None
LEGACY_HELPER=None
L=None
CACHE={}
OWN='scripts/finite57_selected_old_affine_source_applications_20261025.py'
LEGACY='scripts/finite57_selected_legacy_receipt_applications_20261025.py'
LEGACY_SHA='3a1d80666a78b1876dc476ae8d44394cad5204b4127118102301ed830c52440c'
D16_SOURCE='results/fresh_d16_partial_weight_transfer_20261017.json'
D16_NINE='results/fresh_order16_d_distance9_half_transfer_20261018.json'
D17_SOURCE='results/fresh_d17_partial_weight_transfer_20261017.json'
D17_REPORT='results/fresh_d17_partial_weight_transfer_20261017.md'
D17_SINGLE='results/fresh_d17_singleton_pair_coupon_20261017.json'
D17_E6='results/fresh_order17_e6_unit_coupon_zero_transfer_20261018.json'
D17_PROPER='results/fresh_order17_affine_d_three_quarter_transfer_20261017.json'
AFFINE_REVIEW='results/audit_finite57_selected_affine_completion_chain_fresh_20261024.md'
AFFINE_ROOT='results/review_finite57_selected_affine_completions_root_20261025.md'
RESIDUAL='results/fresh_distance9_supporting_extensions_review_20261017.md'
DIRECT='scripts/check_finite57_terminal_collectors_direct_receivers_20261025.py'
DIRECT_SHA='a2a9d61a7637e45d5a4d7af1c3283e120ca7685d07e9f686ae7b6ee585520eca'
CONTRACT='results/finite57_selected_old_affine_source_application_contract_20261025.md'
PINS={
 CONTRACT:'2756281c94fae7de3ed6edd89264018467e66b42cc9fec78005bf58a21acce39',
 D16_SOURCE:'6d3bbabeba06ea257f4d1e783ae359b1cbf9448287ec8b058969152c827f3ff8',
 D16_NINE:'90a9ebd1b975415eda67a1259e52ee9ea6207cb83cc069a2bf3d6f06e1046100',
 D17_SOURCE:'116838eb27f56e5bcf0ad2350cf3d7deb2834787dc9aab5fa0cee544dadb9b71',
 D17_REPORT:'46af0e7682a4fa9b6cda516679cb690abd05e42f5a9629b5a23bf76a3b244e48',
 D17_SINGLE:'04b875fb77babed891db12b97796fa619843f0908f9217b3c7707a5044b94621',
 D17_E6:'e6c9bb5d5820880f413cfadf7998330f91040077968d822d498471be0bda4d0b',
 D17_PROPER:'8356c0b60c8d28da69a884609371e87e4206ce14976ca20b0fa1ec3737e412c6',
 AFFINE_REVIEW:'e740b5ccf7ed66546ceca59507a3c3fdedd6179069732d12f73ca6dbb37c3555',
 AFFINE_ROOT:'6fe04d21f048a55997b583bc7bc54689d0f79ac994f008526795fb981cbf719f',
 RESIDUAL:'4a074a5c70634d64581c90c1884c142bd6b560c05c1760c416db712f03062677',
 DIRECT:DIRECT_SHA,
}
PREMISES={
 'finite57.oldAffine.exact_source_maps':
 'The complete source-native map, exact original physical high/core/base/marks and source-case selectors apply to THIS actual whole target. Counts or unmarked-pattern membership do not substitute for the full map.',
 'finite57.oldAffine.complete_contact_completion':
 'The complete all-label contact/end/native-gap and residual/no-bridge statements apply at the retained exact numerical predicates. Missing-root newness and uniqueness/exterior-root exclusion precede every forced completion zero. Recorded inverse/signature/inertia propositions remain explicit accepted computed inputs unless their exact identity is checked in this call.',
 'finite57.oldAffine.current_whole_receipts':
 'The exact positive-affine-kernel/radical or no-contact seal makes the specified actual paying factors whole at the original corners, including indirect supporting-white tails. Fully validated CURRENT receipt values or explicit current scalar evaluations then apply at those SAME factors, without obsolete finite units or selected half credits.'}


def need(ok,message):
    if not ok:raise ValueError(message)


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=True)


def rational(value):
    need(type(value) in (str,int) and len(str(value))<=128,'bounded exact source rational')
    return F(value)


def pin(path,expected=None):
    digest=PINS[path] if expected is None else expected
    need(callable(SOURCE_PINNER) and type(digest) is str and re.fullmatch('[0-9a-f]{64}',digest),
         'mandatory exact source byte pin')
    SOURCE_PINNER(path,digest)
    return dict(path=path,sha256=digest)


def source(path):
    pin(path);need(callable(SOURCE_READER),'mandatory complete source reader')
    if path not in CACHE:CACHE[path]=SOURCE_READER(path,PINS[path])
    need(type(CACHE[path]) is dict,'complete pinned source object')
    return CACHE[path]


def one(rows,predicate,message):
    need(type(rows) is list,'complete source list')
    matched=[p for p in rows if predicate(p)]
    need(len(matched)==1,message)
    return matched[0]


def acceptance(ids):
    ref=SOURCE_ACCEPTANCE
    need(type(ref) is dict and set(ref)=={'path','sha256'} and
         type(ref['path']) is str and ref['path'].startswith('results/') and
         type(APPLICATION_SOURCE_SHA) is str and re.fullmatch('[0-9a-f]{64}',APPLICATION_SOURCE_SHA),
         'separate exact supplemental root acceptance and helper-source identity')
    pin(OWN,APPLICATION_SOURCE_SHA);pin(ref['path'],ref['sha256'])
    doc=SOURCE_READER(ref['path'],ref['sha256'])
    need(doc['status']=='CURRENT_SOURCE_APPLICATION_PREMISES_ROOT_ACCEPTED' and
         doc['accepted_application_sources'].get(OWN)==APPLICATION_SOURCE_SHA and
         type(doc['accepted_supplemental_premise_ids']) is list and
         set(ids)<=set(doc['accepted_supplemental_premise_ids']),
         'root accepts THIS new source and its explicit supplemental statements')
    return dict(ref)


def native_full_map(key,W,g):
    graph=W.classical_graph(key['target_kind'],key['target_order'],affine=True)
    base,mp=key['base_nodes'],key['induced_map']
    for ai,i in enumerate(base):
        for aj,j in enumerate(base):
            label=graph[mp[ai]].get(mp[aj],0) if ai!=aj else 0
            need(rational(g[i][j])==(2 if ai==aj else -1 if label==3 else 0) and
                 label in (0,3),'THIS entire simply-laced full map, every edge and nonedge')
    return graph


def original_current(key,units,gram,white):
    """CURRENT actual receipt at this exact physical star; no old units used."""
    p,_=L.authoritative(key);context=L.context(gram,white)
    actual_corners={tuple(q['preferred_corner']) for q in p['all_corners']}
    matches=[]
    for shard in range(4):
        doc,evidence=L.unit(units,DIRECT,f'original-{shard}',DIRECT_SHA,
            'PASS_NAMED_CURRENT_TERMINAL_ADAPTER','consumed_inputs',(L.TEX,L.LITERAL))
        for ordinal,row in enumerate(doc['result']['rows']):
            if row['context_sha256']!=context['context_sha256']:continue
            if {tuple(c['marks']) for c in row['corners']}!=actual_corners:continue
            selected=[c for c in row['corners'] if c['marks']==key['selected_marks']]
            if (len(selected)!=1 or selected[0]['kind']!='ideal' or
                    selected[0]['forest']!=[['D',key['target_order']]]):continue
            chosen=selected[0]
            need(key['marked_image'] in chosen['marked_images'],
                 'complete actual selected affine marked-image orbit contains THIS full-key image')
            corners=row['corners'];marks=[tuple(c['marks']) for c in corners]
            need(len(corners)==len(set(marks))==4 and row['direct_current_receiver_condition'] is True and
                 row['legacy_pointwise_bounds_used'] is False and
                 len(row['final_collector_joins'])==4 and all(rational(c['weight'])>=0 for c in corners) and
                 sum((rational(c['weight']) for c in corners),F())==rational(row['total'])>=1 and
                 all(not c.get('guard_possible',False) or rational(c['weight'])<=F(1,2) for c in corners),
                 'complete actual CURRENT distinct-corner sum with conservative unknown negative guards')
            matches.append((row,evidence,ordinal))
    if not matches:return None
    row,evidence,ordinal=min(matches,key=lambda q:canonical([q[0]['source'],q[0]['index'],q[1]['batch']]))
    return dict(physical=context,paying_corners=[dict(marks=c['marks'],current_lower=c['weight'])
        for c in row['corners']],actual_receipts=[evidence],
        direct_current_source_case=[row['source'],row['index']],actual_result_selector=['result','rows',ordinal],
        current_total=row['total'],current_corner_rows=row['corners'],
        final_collector_joins=row['final_collector_joins'],
        equivalent_current_receipts=[dict(source=r['source'],case=r['index'],receipt=e['path'])
                                    for r,e,_ in matches],alternative_receipt_totals_not_added=True)


def e6_current(key,W,core,white):
    doc=source(D17_E6);rows=doc['changed_cases']
    need([p['case'] for p in rows]==[1981,2374,2549,2551],'complete exact four E6 unit source cases')
    matches=[(i,p) for i,p in enumerate(rows) if p['case']==key['high_index']]
    if not matches:return None
    need(len(matches)==1,'one exact source E6 whole-factor occurrence')
    ri,p=matches[0];g=L.matrix_values(core['gram'])
    need(p['original_marks']==key['selected_marks'] and key['marked_image'] in doc['changed_pairs'] and
         W.affine_weight('D',17,*key['marked_image'])==0 and
         p['distinct_unit_marks']!=key['selected_marks'],'THIS whole affine zero uses a distinct E6/D11 unit')
    p_high,_=L.authoritative(key)
    need(p['distinct_unit_marks'] in [r['preferred_corner'] for r in p_high['all_corners']],
         'paying E6/D11 marks are an ORIGINAL adjacent corner')
    roles={}
    for typ,order,lengths,distance in (('E6',6,[1,2,2],2),('D11',11,[1,1,8],1)):
        role=p['coupon_role'][typ];nodes=role['nodes'];mark=role['actual_mark']
        need(len(nodes)==len(set(nodes))==order and mark in nodes,'complete physical source factor/mark')
        adj={i:[j for j in nodes if j!=i and g[i][j]] for i in nodes}
        need(all(g[i][i]==2 and all(g[i][j]==-1 for j in adj[i]) for i in nodes) and
             sum(map(len,adj.values()))==2*(order-1),
             'literal simply-laced finite source tree')
        forks=[i for i in nodes if len(adj[i])==3]
        need(forks==[role['fork']] and all(len(adj[i]) in (1,2,3) for i in nodes),
             'unique physical source fork')
        ds={forks[0]:0};todo=[forks[0]]
        for i in todo:
            for j in adj[i]:
                if j not in ds:ds[j]=ds[i]+1;todo.append(j)
        leaves=[i for i in nodes if len(adj[i])==1]
        eligible=sorted(i for i in leaves if ds[i]==distance)
        need(set(ds)==set(nodes) and sorted(ds[i] for i in leaves)==lengths and
             role['leaf_arm_lengths']==lengths and eligible==role['eligible_marks'] and
             len(eligible)==2 and mark in eligible,'exact complete source long-E6-end/short-D11-leaf role')
        roles[typ]=role
    small,large=roles['E6'],roles['D11']
    need(set(small['nodes']).isdisjoint(large['nodes']) and
         set(small['nodes']+large['nodes'])==set(white)|set(p['distinct_unit_marks']) and
         {small['actual_mark'],large['actual_mark']}==set(p['distinct_unit_marks']) and
         all(g[i][j]==0 for i in small['nodes'] for j in large['nodes']),
         'complete actual separated whole paying factors at THIS corner')
    evaluations=[dict(function='finite_cross_weight',small_type=['E',6],small_mark=m,
        large_type=['D',11],large_mark=n,result=str(W.finite_cross_weight(('E',6),m,('D',11),n)),
        execution_scope='THIS_ACTUAL_F_CALL_NOT_A_PRIOR_RECEIPT')
        for m in sorted(W.E6_LONG_ENDS) for n in (9,10)]
    need(len(evaluations)==4 and min(rational(r['result']) for r in evaluations)>=1,
         'every complete canonical E6/D11 eligible mark pair has CURRENT unit')
    return dict(physical=L.context(core['gram'],white),
        paying_corners=[dict(marks=p['distinct_unit_marks'],current_lower='1',
                             source_selector=['changed_cases',ri,'coupon_role'])],
        actual_receipts=[],current_total='1',consumer_literal_evaluations=evaluations,
        exact_source_selector=[D17_E6,'changed_cases',ri],whole_factor_physical_roles=roles)


def d17_whole(key,W,units):
    if key['target_order']!=17 or len(key['base_nodes'])!=17:return None
    p,core=L.authoritative(key);g=L.matrix_values(core['gram']);n=len(g)
    if n!=19 or p['method']!='radical' or core['signature']!=[17,1,1]:return None
    u,v=L.opposite_pairs(p,n);black=set(u+v);omitted=sorted(black-set(key['selected_marks']))
    white=sorted(set(range(n))-black)
    need(set(key['base_nodes'])==set(range(n))-set(omitted) and p['omitted_black']==omitted and
         p['possible_contacts']==omitted and len(omitted)==2 and
         any(a in u and b in v or b in u and a in v for a,b in [omitted]),
         'whole17 affine block and exactly two ORIGINAL adjacent omitted black contacts')
    L.source_physical(key,core['gram'],white,True);native_full_map(key,W,core['gram'])
    block=one(p['blocks'],lambda q:q['nodes']==key['base_nodes'],
              'the actual preferred whole affine block')
    kernel=[rational(x) for x in block['positive_kernel']]
    need(len(kernel)==17 and all(x>0 for x in kernel) and
         all(sum(g[i][j]*kernel[t] for t,j in enumerate(key['base_nodes']))==0
             for i in key['base_nodes']),'exact positive affine kernel identity checked in THIS call')
    need(type(core['radical']) is list and len(core['radical'])==1,'exhibited full-core radical relation')
    z=[rational(x) for x in core['radical'][0]]
    need(len(z)==n and any(z) and all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n)) and
         z[omitted[0]]*z[omitted[1]]<0,
         'exact full-core zero relation with opposite nonzero omitted coefficients')
    current=original_current(key,units,core['gram'],white)
    singleton=False;e6=False
    if current is None:
        current=e6_current(key,W,core,white);e6=current is not None
    if current is None:
        doc=source(D17_SINGLE);q=doc['literal_core']
        need(doc['affected_overlap_cases']==[key['high_index']] and key['marked_image']==doc['removed_pair'] and
             q['gram']==core['gram'] and q['white']==white and q['affine_nodes']==key['base_nodes'] and
             [rational(x) for x in q['radical']]==z,
             'no actual original receipt: require exact independent singleton whole-core source')
        paid=[];evaluations=[]
        for ci,corner in enumerate(q['corners']):
            if corner.get('half_coupon')!='1/2':continue
            parts=corner['components']
            small=one(parts,lambda c:c['type']=='A4' and c['role']=='endpoint',
                      'exact singleton whole A4 endpoint')
            large=one(parts,lambda c:c['type']=='D13' and c['role']=='long_arm_neighbor_of_fork',
                      'exact singleton whole D13 long-arm neighbor')
            need(len(small['nodes'])==4 and len(large['nodes'])==13 and
                 set(small['nodes']).isdisjoint(large['nodes']) and
                 set(small['nodes']+large['nodes'])==set(white)|set(corner['marks']) and
                 {small['mark'],large['mark']}==set(corner['marks']) and
                 all(g[i][j]==0 for i in small['nodes'] for j in large['nodes']) and
                 sum(bool(g[small['mark']][j]) for j in small['nodes'] if j!=small['mark'])==1,
                 'same actual separated whole coupon factors and physical endpoint')
            branch=one(large['nodes'],lambda i:sum(bool(g[i][j]) for j in large['nodes'] if j!=i)==3,
                       'unique physical D13 fork')
            need(g[branch][large['mark']]==-1 and
                 sum(bool(g[large['mark']][j]) for j in large['nodes'] if j!=large['mark'])==2,
                 'large physical mark is the unique long-arm neighbor, not a short leaf or center')
            values=[W.finite_cross_weight(('A',4),m,('D',13),9) for m in (0,3)]
            need(min(values)>=F(1,2),'CURRENT A4 endpoint/D13 long-neighbor halves')
            paid.append(dict(marks=corner['marks'],current_lower=str(min(values)),
                             source_selector=['literal_core','corners',ci]))
            evaluations.extend(dict(function='finite_cross_weight',small_type=['A',4],small_mark=m,
                large_type=['D',13],large_mark=9,result=str(value),
                execution_scope='THIS_ACTUAL_F_CALL_NOT_A_PRIOR_RECEIPT')
                for m,value in zip((0,3),values,strict=True))
        need(len(paid)==2 and len({tuple(c['marks']) for c in paid})==2 and
             all(c['marks']!=key['selected_marks'] for c in paid) and
             sum((rational(c['current_lower']) for c in paid),F())>=1,
             'two distinct actual singleton CURRENT halves; no ideal unit is credited')
        current=dict(physical=L.context(core['gram'],white),paying_corners=paid,actual_receipts=[],
                     current_total='1',consumer_literal_evaluations=evaluations)
        singleton=True
    named=[i for i,q in enumerate(source(D17_SOURCE)['cases']) if q['index']==key['high_index']]
    if named:
        need(len(named)==1,'one exact named P17 source occurrence')
        q=source(D17_SOURCE)['cases'][named[0]]
        need(q['gram']==core['gram'] and q['white']==white and q['affine_nodes']==key['base_nodes'] and
             [rational(x) for x in q['radical']]==z and
             [rational(x) for x in q['positive_affine_kernel']]==kernel,
             'same named P17 physical core/kernel/radical/whole factor')
    proof_source=D17_E6 if e6 else D17_SINGLE if singleton else D17_SOURCE if named else D17_REPORT
    selector=current['exact_source_selector'][1:] if e6 else ['literal_core'] if singleton else ['cases',named[0]] if named else dict(
        format='markdown',section='Actual-kernel sealing applies uniformly')
    return dict(application_kind='EXACT_CURRENT_RECEIVER_SOURCE_BINDING',source=proof_source,
        proof_selector=selector,
        map_source=L.HIGH,map_selector=['isolation','cases',key['high_index']],
        physical=current['physical'],source_sealing=dict(kind='EXACT_WHOLE_D17_KERNEL_RADICAL_SEAL',
            complete_actual_positive_kernel=[str(x) for x in kernel],full_radical=[str(x) for x in z],
            omitted_original_black=omitted,complete_possible_contacts=p['possible_contacts'],
            exact_kernel_and_radical_linear_identities_checked_in_this_call=True,
            inherited_signature=core['signature'],signature_and_actual_Lorentz_geometry_regenerated=False,
            geometry='PSD actual corner plus positive affine kernel seals every additional white against the whole17factor. Pairing the actual zero relation leaves neither omitted contact or BOTH; BOTH is the actual adjacent omitted-black/white/black path2. Under NO-S8 every further white is orthogonal to the ENTIRE old core, sealing every paying factor and indirect tail.'),
        paying_corners=current['paying_corners'],actual_receipts=current['actual_receipts'],
        current_receiver_binding=current,consumer_literal_evaluations=current.get('consumer_literal_evaluations',[]),
        proof_references=[pin(AFFINE_REVIEW),pin(AFFINE_ROOT),pin(RESIDUAL)])


def d16_control(key,p,core,control):
    g=L.matrix_values(core['gram']);J=L.matrix_values(core['inverse']);n=len(g)
    need(core['gauge']==[1]*n and core['signature']==[n-1,1,0] and
         control['allowed_contacts']==p['possible_contacts'] and control['end_sets']==p['end_sets'],
         'same exact norm2 Lorentz core and complete original contact/end domain')
    allowed=control['allowed_contacts'];ends=control['end_sets']
    need(len(ends)==2 and not set(ends[0])&set(ends[1]) and set(ends[0]+ends[1])==set(allowed) and
         all(J[i][j]>0 for end in ends for i in end for j in end) and
         all(J[i][i]+J[j][j]+2*J[i][j]>2 for end in ends for i,j in combinations(end,2)) and
         control['projection_diagonals']=={str(i):str(J[i][i]) for i in allowed},
         'complete positive same-end and two-contact predicates precede forced zeros')
    candidates=[i for i in allowed if J[i][i]<=2]
    need(control['candidate_nodes']==candidates,'complete candidate domain from THIS saved inverse')
    u,v=L.opposite_pairs(p,n);du=L.graph_distances(g,u);dv=L.graph_distances(g,v)
    need(max((du[i] for i in ends[0]),default=0)+max((dv[i] for i in ends[1]),default=0)+2<=8,
         'actual original opposite-end short8 path bound')
    if p['method']=='lorentz_residual':
        need(control['kind']=='at_most_one_core_contact_root' and all(J[i][i]>1 for i in candidates) and
             candidates==[r['contact_node'] for r in p['single_contact_candidates']] and
             control['pair_tests']==p['pair_tests'],'complete ordinary one-contact-root source domain')
        expected=[]
        for i,j in product(candidates,repeat=2):
            a=J[i][j];residual=(2-J[i][i])*(2-J[j][j]);short=min(du[i]+dv[j],du[j]+dv[i])+3
            need(a*a>residual and (a>0 or short<=8),'every ordered contact pair, including repetitions')
            expected.append(dict(nodes=[i,j],projection_square=str(a*a),residual_product=str(residual),
                                 sign=1 if a>0 else -1,shortcut=short if a<0 else None))
        need(expected==control['pair_tests'],'complete source pair-table identity')
    else:
        need(p['method']=='relevant_component' and control['kind']=='one_relevant_root_with_optional_exterior_root' and
             p['relevant_nodes']==[control['relevant_node']] and candidates==p['contact_domains'] and
             len(p['extension_cases'])==1,'exact complete relevant/exterior source domain')
        q=p['extension_cases'][0];good=control['relevant_node'];other=control['exterior_node']
        S=[[2-J[good][good],-J[good][other]],[-J[good][other],2-J[other][other]]]
        determinant=S[0][0]*S[1][1]-S[0][1]**2
        need(q['node']==good and good in key['base_nodes'] and other not in key['base_nodes'] and
             1<J[good][good]<2 and control['pair_checks']==q['pair_checks'] and
             control['double_checks']==q['double_checks'] and
             control['residual_gram']==[[str(x) for x in r] for r in S] and
             rational(q['joint']['determinant'])==determinant and
             ((determinant==0 and S[0][0]>0 and S[1][1]>0 and S[0][1]<0 and
               q['joint']['kind']=='positive_residual_relation' and q['joint']['tail']==0) or
              (determinant>0 and S[1][1]/determinant==2 and
               q['joint']['kind']=='saturated_single_tail' and q['joint']['tail']==1)) and
             control['joint_tail_bound']==q['joint']['tail'] in (0,1) and
             control['joint_component_order_bound']==len(key['base_nodes'])+1+q['joint']['tail']<16,
             'complete exterior multiplicity/residual/no-bridge bound before excluding exterior entry')
    return g,J


def d16_saved(key,W,units):
    if key['target_order']!=16:return None
    p,core=L.authoritative(key);matches=[]
    for path,index_field in ((D16_SOURCE,'index'),(D16_NINE,'case')):
        doc=source(path)
        for ri,row in enumerate(doc['cases']):
            if row[index_field]!=key['high_index']:continue
            for ci,q in enumerate(row['completions']):
                embeddings=q.get('embeddings',[dict(map=q.get('embedding'),marked_image=q.get('marked_image'))])
                for mi,mp in enumerate(embeddings):
                    if mp['map']==key['induced_map'] and mp['marked_image']==key['marked_image']:
                        matches.append((path,ri,row,ci,q,mi))
    if not matches:return None
    need(len(matches)==1,'unique exact disjoint old/new D16 source-map occurrence')
    path,ri,row,ci,q,mi=matches[0];selector=['cases',ri,'completions',ci]
    need(row['core_key']==key['core_key'] and row['core_gram']==core['gram'] and
         row['core_inverse']==core['inverse'] and row['preferred_corner']==key['selected_marks'] and
         row.get('base_nodes',row.get('marked_component'))==key['base_nodes'] and
         row['method']==p['method'],'same entire physical source D16 occurrence')
    g,J=d16_control(key,p,core,row['constraints']);u,v=L.opposite_pairs(p,len(g))
    white=sorted(set(range(len(g)))-set(u+v));L.source_physical(key,core['gram'],white,True)
    graph=native_full_map(key,W,core['gram']);n=len(g);base=key['base_nodes']
    image=dict(zip(key['induced_map'],base));missing=sorted(set(range(16))-set(image))
    image.update({t:n+i for i,t in enumerate(missing)});nn=n+len(missing)
    complete=[r+[F()]*len(missing) for r in g]+[[F()]*nn for _ in missing]
    for t in missing:
        complete[image[t]][image[t]]=2
        for old,label in graph[t].items():
            need(label==3,'D completion uses exact required simple target bonds')
            complete[image[t]][image[old]]=complete[image[old]][image[t]]=F(-1)
    required=[dict(new_node=i,core_nodes=[j for j in range(n) if complete[i][j]])
              for i in range(n,nn) if any(complete[i][j] for j in range(n))]
    outcome=q['outcome'];excluded=outcome!='sealed_completion'
    if outcome=='excluded_forced_entry_projection':
        need(path==D16_NINE and q['forced_entries'] and
             all(e['native_label']==3 and e['physical_contact_square']==1 and
                 e['new_target_node'] in missing and
                 key['induced_map'][base.index(e['old_literal_node'])]==e['old_target_node'] and
                 graph[e['new_target_node']].get(e['old_target_node'])==3 and
                 rational(e['projection_lower_bound'])==J[e['old_literal_node']][e['old_literal_node']]>2
                 for e in q['forced_entries']),
             'exact source special first-entry witnesses belong to THIS full native map')
    else:
        need(q['new_roots']==list(range(n,nn)) and q['required_core_contacts']==required and
             all(len(r['core_nodes'])==1 for r in required),
             'same source forced first roots and contact incidences')
    if outcome=='excluded_multiple_core_contact_roots':
        need(len(required)>1 and row['constraints']['kind']=='at_most_one_core_contact_root',
             'two required genuine roots contradict complete ordinary uniqueness')
    elif outcome in ('excluded_projection_above_norm','excluded_forced_entry_projection'):
        need(required and any(J[i][i]>2 for r in required for i in r['core_nodes']),
             'the given full D16 map forces projection above actual norm2')
    elif outcome=='excluded_lorentz_index':
        need(len(required)==1,'exactly one allowed entry precedes forced principal completion')
        if row['constraints']['kind']=='one_relevant_root_with_optional_exterior_root':
            need(q['exterior_root_exclusion']['required_factor_order']==16 and
                 q['exterior_root_exclusion']['maximum_order_if_exterior_root_present']==
                    row['constraints']['joint_component_order_bound']<16,
                 'the whole target rules out the exterior root BEFORE completing zeros')
        need(L.matrix_values(q['gram'])==complete and type(q['signature']) is list and
             len(q['signature'])==3 and q['signature'][1]>1 and sum(q['signature'])==nn and
             type(q['forced_completion_reason']) is str and q['forced_completion_reason'],
             'same forced actual principal and its explicitly inherited computed index-two predicate')
    elif outcome=='sealed_completion':
        need(not missing and not required and row['constraints']['candidate_nodes']==[],
             'whole D16 and no possible further old-core contact')
        current=original_current(key,units,q['gram'],q['white'])
        if current is None:
            need(path==D16_NINE,'no actual CURRENT original receipt for this old whole D16 source star')
            role=q['current_payment']['existing_whole_component_coupon']
            need(role['small_type']=='D5' and role['small_role']=='any_node' and
                 role['large_type']=='D11' and role['large_role']=='fork_center' and
                 role['marks']!=key['selected_marks'],'exact new-distance9 distinct whole D5/D11 center coupon')
            corner=one(q['corners'],lambda c:c['marks']==role['marks'],'same source finite coupon corner')
            for typ,nodes,mark in (('D5',role['small_nodes'],role['small_mark']),
                                   ('D11',role['large_nodes'],role['large_mark'])):
                part=one(corner['components'],lambda c:c['type']==typ and c['nodes']==nodes,
                         'same source whole typed coupon component')
                need(mark in nodes,'actual coupon mark belongs to this whole component')
                if typ=='D11':
                    need(part['roles'].get(str(mark))=='fork_center' and
                         sum(bool(rational(q['gram'][mark][j])) for j in nodes if j!=mark)==3,
                         'actual D11 marked node is its unique fork center')
            selected=W.affine_weight('D',16,*key['marked_image'])
            values=[W.finite_cross_weight(('D',5),m,('D',11),8) for m in range(5)]
            lower=min(values)
            need(selected==F(1,2) and lower>=F(1,2),'CURRENT exact selected half plus conservative any-D5/center half')
            current=dict(physical=L.context(q['gram'],q['white']),
                paying_corners=[dict(marks=key['selected_marks'],current_lower=str(selected)),
                                dict(marks=role['marks'],current_lower=str(lower))],
                actual_receipts=[],current_total=str(selected+lower),
                consumer_literal_evaluations=[dict(function='finite_cross_weight',small_type=['D',5],
                    small_mark=m,large_type=['D',11],large_mark=8,result=str(value),
                    execution_scope='THIS_ACTUAL_F_CALL_NOT_A_PRIOR_RECEIPT')
                    for m,value in enumerate(values)])
    else:raise ValueError('unsupported exact D16 source outcome')
    return dict(application_kind='EXCLUDED_COMPLETION_UNDER_NO_S8' if excluded else
                'EXACT_CURRENT_RECEIVER_SOURCE_BINDING',source=path,proof_selector=selector,
        map_source=path,map_selector=selector+(['embeddings',mi] if path==D16_SOURCE else ['embedding']),
        physical=L.context(core['gram'],white),source_sealing=dict(
            kind='ACCEPTED_COMPLETE_D16_CONTACT_AND_NO_BRIDGE_SOURCE_IMPLICATION',
            complete_constraints=row['constraints'],exact_outcome=outcome,
            source_completion=q,requested_target_to_literal=[image[t] for t in range(16)],
            numerical_signature_inverse_and_inertia_regenerated=False,
            forced_zero_scope='Only after complete one-root or relevant/exterior/no-bridge argument; no old exterior component enters by an indirect tail.'),
        paying_corners=[] if excluded else current['paying_corners'],
        actual_receipts=[] if excluded else current['actual_receipts'],
        current_receiver_binding=None if excluded else current,
        consumer_literal_evaluations=[] if excluded else current.get('consumer_literal_evaluations',[]),
        source_branch_conclusion=dict(kind='THIS_EXACT_ACTUAL_COMPLETION_IMPOSSIBLE_UNDER_NO_S8' if excluded
            else 'CURRENT_DISTINCT_CORNER_TOTAL_AT_LEAST_ONE_UNDER_NO_S8',source_outcome=outcome),
        proof_references=[pin(AFFINE_REVIEW),pin(AFFINE_ROOT),pin(RESIDUAL)])


def d17_proper(key,W,units):
    if key['target_order']!=17 or len(key['base_nodes'])>=17:return None
    doc=source(D17_PROPER)
    matches=[(i,p) for i,p in enumerate(doc['cases']) if p['index']==key['high_index'] and type(p['proof']) is dict]
    if not matches:return None
    need(len(matches)==1,'unique complete proper-D17 source occurrence')
    ri,row=matches[0];proof=row['proof'];p,core=L.authoritative(key)
    maps=[(i,q) for i,q in enumerate(proof['maps']) if q['map']==key['induced_map']]
    if not maps:return None
    need(len(maps)==1 and key['marked_image'] in doc['changed_pairs'] and
         W.affine_weight('D',17,*key['marked_image'])==0,
         'exact full proper map in the current zero orbit; no historical three-quarter credit')
    mi,q=maps[0]
    need(row['core_key']==key['core_key'] and row['method']==p['method']=='lorentz_residual' and
         proof['base_nodes']==key['base_nodes'] and proof['core_gram']==core['gram'] and
         proof['core_inverse']==core['inverse'],'same original physical core and complete source map')
    control=proof['complete_residual_contact_control']
    g,J=d16_control(key,p,core,control);n=len(g);u,v=L.opposite_pairs(p,n)
    white=sorted(set(range(n))-set(u+v));L.source_physical(key,core['gram'],white,True)
    graph=native_full_map(key,W,core['gram']);image=dict(zip(key['induced_map'],key['base_nodes']))
    missing=sorted(set(range(17))-set(image));entries=[]
    for new in missing:
        neighbors=[(old,node) for old,node in image.items() if graph[new].get(old)]
        need(len(neighbors)<=1,'every source first-entry target node has at most one old-base neighbor')
        if neighbors:
            old,node=neighbors[0]
            need(graph[new][old]==3 and node in control['allowed_contacts'],
                 'genuine first supporting-white contact uses THIS complete native domain')
            entries.append(dict(new_target_node=new,core_node=node,native_label=3,
                                projection_square=str(J[node][node])))
    need(entries and q['missing_target_nodes']==missing and q['forced_entries']==entries,
         'complete exact first-entry requirements of THIS full native map')
    outcome=q['outcome']
    if outcome=='excluded_multiple_core_contact_roots':
        need(len(entries)>1,'distinct required first roots contradict at-most-one contact root')
    elif outcome=='excluded_projection_above_norm':
        need(len(entries)==1 and rational(entries[0]['projection_square'])>2 and
             rational(q['projection_square'])==rational(entries[0]['projection_square']),
             'actual mandatory simple entry has projection above physical norm2')
    elif outcome=='excluded_forced_lorentz_index':
        need(len(entries)==1 and entries[0]['core_node'] in control['candidate_nodes'],
             'the unique allowed entry exhausts all old-core contacts BEFORE forced zeros')
        image.update({t:n+i for i,t in enumerate(missing)});nn=n+len(missing)
        h=[r+[F()]*len(missing) for r in g]+[[F()]*nn for _ in missing]
        for t in missing:
            h[image[t]][image[t]]=2
            for old,label in graph[t].items():
                need(label==3,'exact proper affine-D target bond')
                h[image[t]][image[old]]=h[image[old]][image[t]]=F(-1)
        need(q['target_to_literal']==[image[t] for t in range(17)] and
             L.matrix_values(q['gram'])==h and len(q['signature'])==3 and
             q['signature'][1]>1 and sum(q['signature'])==nn,
             'same forced actual principal and explicit inherited computed index-two signature')
    else:raise ValueError('unsupported exact proper-D17 source outcome')
    selector=['cases',ri,'proof','maps',mi]
    return dict(application_kind='EXCLUDED_COMPLETION_UNDER_NO_S8',source=D17_PROPER,
        proof_selector=selector,map_source=D17_PROPER,map_selector=selector+['map'],
        physical=L.context(core['gram'],white),source_sealing=dict(
            kind='COMPLETE_PROPER_D17_ONE_CONTACT_NO_BRIDGE_IMPLICATION',
            complete_residual_contact_control=control,full_source_map_record=q,
            missing_root_newness='Other original corner blocks are orthogonal to the base and omitted blacks are absent. The unique first root excludes any second old-core contact or indirect old-block bridge; every later missing root is genuinely new and old-core-orthogonal.',
            numerical_inertia_and_actual_geometric_precursors_regenerated=False),
        paying_corners=[],actual_receipts=[],source_branch_conclusion=dict(
            kind='THIS_EXACT_ACTUAL_COMPLETION_IMPOSSIBLE_UNDER_NO_S8',source_outcome=outcome),
        proof_references=[pin(AFFINE_REVIEW),pin(AFFINE_ROOT),pin(RESIDUAL)])


# Source-only affine-B fragment for the unified old-affine application helper.
# Shared parent names: L=LEGACY_HELPER3a1d, F=Fraction, need, one,
# source(own pinned path), pin(own pinned path), canonical, rational.
# The parent merges B_SOURCE_PINS into its OWN registry; do not mutate L.PINS.
B14_SOURCE = 'results/fresh_order14_b_terminal_residual_transfer_20261018.json'
B14_PRODUCER = 'scripts/fresh_order14_b_terminal_residual_transfer_20261018.py'
B14_REPORT = 'results/fresh_order14_b_terminal_residual_transfer_20261018.md'
B15_SOURCE = 'results/fresh_order15_bd_terminal_coupon_transfer_20261018.json'
B15_PRODUCER = 'scripts/fresh_order15_bd_terminal_coupon_transfer_20261018.py'
B15_REPORT = 'results/fresh_order15_bd_terminal_coupon_transfer_20261018.md'
B_CONTACT_PRODUCER = 'scripts/fresh_d16_partial_weight_transfer_20261017.py'
B_RESIDUAL_REPORT = 'results/fresh_distance9_supporting_extensions_review_20261017.md'
B_AFFINE_PROOF = 'results/audit_finite57_selected_affine_completion_chain_fresh_20261024.md'
B_AFFINE_ROOT = 'results/review_finite57_selected_affine_completions_root_20261025.md'
B_GENERAL_REPORT = 'results/finite57_old_affine_b_given_map_source_contract_20261025.md'
B13_OLD = 'results/fresh_order13_affine_b_singleton_transfer_20261017.json'
B13_FORK = 'results/fresh_order13_affine_b_fork_projection_transfer_20261017.json'
B13_OLD_PRODUCER = 'scripts/fresh_order13_affine_b_singleton_transfer_20261017.py'
B13_FORK_PRODUCER = 'scripts/fresh_order13_affine_b_fork_projection_transfer_20261017.py'
B13_OLD_REPORT = 'results/fresh_order13_affine_b_singleton_transfer_20261017.md'
B13_FORK_REPORT = 'results/fresh_order13_affine_b_fork_projection_transfer_20261017.md'
B13_QUIET = 'scripts/check_finite57_quiet_b13_current_priority_bridge_20261025.py'
B13_QUIET_SHA = '5e5c1b9736701dbb48ff3f3dc6bba6b7f87a314758d5fe0b16802607044263e4'
B13_QUIET_ROOT = 'results/review_finite57_quiet_b13_current_priority_bridge_root_20261025.md'
B13_QUIET_CONTRACT = 'results/finite57_quiet_b13_current_priority_bridge_contract_20261025.md'
B_SOURCE_PINS = {
    B_GENERAL_REPORT:'c51b37dd9da254ed7e638dc80a18cae9726016e39ea8b16ce43d2be23be6ed79',
    'results/fresh_distance8_high_overlap_20261017.json':'f8a2ed3bdef6064febf7336bbe4ad75d71c7d35b5747d9f57b0ebf26b03d0cc9',
    B13_OLD:'d21a7a2daa5f215dae29f48730d1f0b6fc6c5ee0632715dcb3db4882a79b372c',
    B13_FORK:'e9b32647aa9f29bc3da18388bba719718b6ac1ea520d01a1d9029f2b7f270d55',
    B13_OLD_PRODUCER:'35749c95282c04909e1861706ff8115cdd45709e1933748ebae231f92a2ae0f5',
    B13_FORK_PRODUCER:'ba82121f61158352d8d9b8c9d07901f4d74097fc608a0593fe6d74002469cd9e',
    B13_OLD_REPORT:'84262d376910422fe21b252af4269509b3dd319db5b4ae54da2273c64ccca2f0',
    B13_FORK_REPORT:'33144882bd974c15902d7c734c33f9a4f1586344e71a2cfbf91f8036ed3f13d0',
    B13_QUIET:B13_QUIET_SHA,
    B13_QUIET_ROOT:'2d1ae426c9c9eed1b81fa8fd88ef55b7c9f21439af8e042c4af1da09f9fac4fd',
    B13_QUIET_CONTRACT:'9667c9028656440af5123f151c1ceccc0f1346eda538fa6c198f61332f38f59f',
    B14_SOURCE:'54f2e36ba4c3e5815e8664504738009267ee7c9880c8510480f666c09d9c951c',
    B14_PRODUCER:'8f8f25effcb6a2d5088865fdf9d6bd197f1af114b097815a715e1a92b2834936',
    B14_REPORT:'fe70262bedc1dd220aa4e91a703b6eba89e8053a3042d25e108833daee7aa435',
    B15_SOURCE:'dda5e90c98c99b99bdf4b667421adfb2dc1637862e7e70a2b906bc1dd5a178ec',
    B15_PRODUCER:'f909ab330ed53c0a0faf219b171e05e010dc9c8f79b00fd37d4fd8bcf55a5486',
    B15_REPORT:'8a8b2667fd05e6d318192918655e06c7b4663c605f956862a111ddfd428a16a3',
    B_CONTACT_PRODUCER:'4f60793bce5c1c5a91d59eda946bbc71cb3956bbb35fa6bc68b61bca28cc55ee',
    B_RESIDUAL_REPORT:'4a074a5c70634d64581c90c1884c142bd6b560c05c1760c416db712f03062677',
    B_AFFINE_PROOF:'e740b5ccf7ed66546ceca59507a3c3fdedd6179069732d12f73ca6dbb37c3555',
    B_AFFINE_ROOT:'6fe04d21f048a55997b583bc7bc54689d0f79ac994f008526795fb981cbf719f',
}


def b13_cross_current(corner,W):
    parts = corner['components'];marks = corner['marks']
    chosen = [one(parts,lambda b:m in b['nodes'],'same whole source marked finite factor') for m in marks]
    if chosen[0] == chosen[1]:
        return F(),[]  # No old same-component unit is borrowed after later deletions.
    alternatives = []
    for si,li in ((0,1),(1,0)):
        small,large = chosen[si],chosen[li]
        st,lt = small['type'],large['type'];sm,lm = marks[si],marks[li]
        if not (st[0] in 'ABD' and st[1:].isdigit() and int(st[1:]) <= 5):
            continue
        role = large['roles'].get(str(lm));large_images = None
        if lt == 'D12' and role in ('fork_center','short_fork_leaf'):
            large_images = [9] if role == 'fork_center' else [10,11]
        elif lt == 'D11' and role == 'short_fork_leaf' and st in ('A1','A2'):
            large_images = [9,10]
        elif lt == 'B10' and role == 'doubled_terminal' and st in ('A1','A2'):
            large_images = [9]
        if large_images is None:
            continue
        sr = int(st[1:]);lr = int(lt[1:])
        evaluations = [dict(function='finite_cross_weight',small_type=[st[0],sr],small_mark=i,
            large_type=[lt[0],lr],large_mark=j,
            result=str(W.finite_cross_weight((st[0],sr),i,(lt[0],lr),j)))
            for i in range(sr) for j in large_images]
        value = min(rational(q['result']) for q in evaluations)
        alternatives.append((value,evaluations))
    if not alternatives:
        return F(),[]
    return max(alternatives,key=lambda t:t[0])


def b13_full_application(key,W,units):
    if not key['affine'] or key['target_kind'] != 'B' or key['target_order'] != 13:
        return None
    if len(key['base_nodes']) != 13 or key['marked_image'] not in ([0,10],[9,11],[9,12],[10,11],[10,12]):
        return None
    p,core = L.authoritative(key);old = source(B13_OLD);fork = source(B13_FORK)
    support = one(fork['reconstructed_high_support'],lambda r:r['index'] == key['high_index'],
                  'complete original B13 support occurrence')
    need(support['embeddings'].count(key['induced_map']) == 1
         and key['marked_image'] in support['marked_images'],
         'THIS complete B13 native map in the complete original source census')
    support_ordinal = fork['reconstructed_high_support'].index(support)
    if key['marked_image'] == [0,10]:
        collection = old['selected_support']['cases']
        geometry = one(collection,lambda r:r['index'] == key['high_index'],'same source B13 singleton seal')
        src = B13_OLD;selector = ['selected_support','cases',collection.index(geometry)]
    elif key['marked_image'] in ([10,11],[10,12]):
        collection = old['conditional_long_orbit']['cases']
        geometry = one(collection,lambda r:r['index'] == key['high_index'],'same original B13 long-orbit seal')
        src = B13_OLD;selector = ['conditional_long_orbit','cases',collection.index(geometry)]
    else:
        collection = fork['cases']
        case = one(collection,lambda r:r['index'] == key['high_index'],'same source B13 fork seal')
        need(case['method'] == 'radical' and case['outcome'] == 'sealed_A2_D11_unit'
             and case['embeddings'].count(key['induced_map']) == 1,'same full affine B13 source branch')
        geometry = case['proof']['geometry'];src = B13_FORK
        selector = ['cases',collection.index(case),'proof','geometry']
    need(geometry['index'] == key['high_index'] and geometry['core_key'] == key['core_key']
         and geometry['affine_nodes'] == key['base_nodes']
         and geometry['gram'] == core['gram'] and geometry['gauge_squared'] == core['gauge']
         and geometry['signature'] == core['signature'] == [13,1,1]
         and geometry['omitted'] == p['omitted_black']
         and p['possible_contacts'] == p['omitted_black'] and p['shortcut_bound'] == 2
         and len(geometry['positive_affine_kernel']) == 13
         and all(rational(v) > 0 for v in geometry['positive_affine_kernel'])
         and len(geometry['omitted_radical_coefficients']) == 2
         and rational(geometry['omitted_radical_coefficients'][0])*
             rational(geometry['omitted_radical_coefficients'][1]) < 0,
         'same complete actual positive-affine-kernel/full-core-radical B13 seal')
    L.source_physical(key,geometry['gram'],geometry['white'])
    g = L.matrix_values(core['gram']);scales = [rational(v) for v in core['gauge']]
    target = W.classical_graph('B',13,affine=True)
    for ai,i in enumerate(key['base_nodes']):
        for aj,j in enumerate(key['base_nodes'][:ai]):
            label = target[key['induced_map'][ai]].get(key['induced_map'][aj],0)
            need((label == 0 and g[i][j] == 0) or
                 (label in (3,4) and g[i][j] < 0 and
                  g[i][j]**2 == (1 if label == 3 else 2)*scales[i]*scales[j]),
                 'same complete physical B13 native labels and nonedges')
    need(W.affine_weight('B',13,*key['marked_image']) == 0,'selected B13 intrinsic CURRENT zero remains zero')
    physical = L.context(geometry['gram'],geometry['white']);evidence = [];evaluations = []
    if key['marked_image'] in ([10,11],[10,12]):
        doc,receipt_evidence = L.unit(units,B13_QUIET,'quiet-b13',B13_QUIET_SHA,
            'PASS_NAMED_QUIET_B13_CURRENT_PRIORITY_BRIDGE','consumed_inputs',(L.LITERAL,L.TEX))
        need(doc['consumed_inputs'].get(B13_OLD) == B_SOURCE_PINS[B13_OLD],
             'actual current long-orbit receipt uses THIS exact original source')
        promotion = one(doc['current_promotion_receipts'],lambda r:r['case'] == key['high_index'],
                        'one exact actual CURRENT B13 long-orbit promotion')
        actual = promotion['current_receipt']
        need(actual['source'] == B13_OLD[8:-5] and actual['index'] == key['high_index']
             and actual['context_sha256'] == physical['context_sha256']
             and actual['direct_current_receiver_condition'] is True
             and actual['legacy_pointwise_bounds_used'] is False
             and promotion['saved_old_promotion_payment_used'] is False
             and [r['marks'] for r in actual['corners']] == [r['marks'] for r in geometry['corners']],
             'actual current promotion is THIS exact source physical four-corner model')
        paid = [dict(marks=r['marks'],current_lower=r['weight']) for r in actual['corners']]
        need(sum((rational(r['current_lower']) for r in paid),F()) ==
             rational(actual['total']) >= 1,'actual CURRENT distinct-corner B13 total')
        evidence.append(receipt_evidence)
    else:
        paid = []
        for ci,corner in enumerate(geometry['corners']):
            if corner['kind'] != 'finite':
                continue
            value,scalars = b13_cross_current(corner,W)
            if value > 0:
                paid.append(dict(marks=corner['marks'],current_lower=str(value),
                    source_selector=selector+['corners',ci]))
                evaluations += scalars
        need(sum((rational(r['current_lower']) for r in paid),F()) >= 1,
             'CURRENT distinct whole B13 source coupons pay; old same-component values are unused')
    need(paid and len({tuple(sorted(r['marks'])) for r in paid}) == len(paid),
         'same actual B13 face corners are distinct, never count a native map as another corner')
    return dict(source=src,proof_selector=selector,map_source=B13_FORK,
        map_selector=['reconstructed_high_support',support_ordinal,'embeddings',
            support['embeddings'].index(key['induced_map'])],physical=physical,
        source_sealing=dict(kind='accepted_complete_positive_affine_kernel_and_actual_zero_core_radical',
            source_signature=core['signature'],source_gauge=core['gauge'],
            positive_affine_kernel=geometry['positive_affine_kernel'],radical=geometry['radical'],
            original_omitted_black=geometry['omitted'],
            opposite_omitted_radical_coefficients=geometry['omitted_radical_coefficients'],
            complete_original_contacts=p['possible_contacts'],independent_contact_path_length=2,
            every_finite_factor_used_is_whole_under_NO_S8=True,
            source_finite_matrix_kernel_and_native_role_predicates_runtime_regenerated=False),
        paying_corners=paid,actual_receipts=evidence,consumer_literal_evaluations=evaluations,
        proof_references=[pin(B_GENERAL_REPORT),pin(B13_OLD_PRODUCER),pin(B13_FORK_PRODUCER),
            pin(B13_OLD_REPORT),pin(B13_FORK_REPORT),pin(B13_QUIET_CONTRACT),pin(B13_QUIET_ROOT),
            pin(B_AFFINE_PROOF),pin(B_AFFINE_ROOT)])


def b_source_frame(key,W):
    """Bind one GIVEN full native map to the original physical source."""
    L.full_key(key)
    if not key['affine'] or key['target_kind'] != 'B' or not 14 <= key['target_order'] <= 17:
        return None
    p,core = L.authoritative(key)
    if len(key['base_nodes']) >= key['target_order']:
        return None
    n = len(core['gram'])
    if core.get('inverse') is None or core['signature'] != [n-1,1,0] or core['gauge'] != [1]*n:
        return None
    if p['method'] not in ('lorentz_residual','relevant_component','marked_component_isolated_by_end_inverses'):
        return None
    g = L.matrix_values(core['gram']);J = L.matrix_values(core['inverse'])
    if len(J) != n:
        return None
    need(all(g[i][i] == 2 and all(g[i][j] == g[j][i] and g[i][j] in (-1,0)
         and J[i][j] == J[j][i] for j in range(i)) for i in range(n)),
         'same complete simply-laced norm2 physical core and symmetric inverse')
    u,v = L.opposite_pairs(p,n);white = sorted(set(range(n))-set(u+v))
    L.source_physical(key,core['gram'],white,True)
    target = W.classical_graph('B',key['target_order'],affine=True)
    for ai,i in enumerate(key['base_nodes']):
        for aj,j in enumerate(key['base_nodes'][:ai]):
            label = target[key['induced_map'][ai]].get(key['induced_map'][aj],0)
            need(label in (0,3) and g[i][j] == (-1 if label == 3 else 0),
                 'THIS complete native map preserves every old simple label and nonedge')
    allowed = p['possible_contacts'];ends = p['end_sets']
    need(type(allowed) is list and len(allowed) == len(set(allowed))
         and all(type(i) is int and 0 <= i < n for i in allowed)
         and type(ends) is list and len(ends) == 2
         and all(type(e) is list and len(e) == len(set(e)) for e in ends)
         and not set(ends[0]) & set(ends[1])
         and set(ends[0]+ends[1]) == set(allowed),
         'complete source necessary contacts and disjoint full two-end cover')
    du,dv = L.graph_distances(g,u),L.graph_distances(g,v)
    need(len(du) == len(dv) == n and min(du[i] for i in v) > 8
         and any(all(du[i] <= k for i in ends[0]) and all(dv[i] <= 6-k for i in ends[1])
                 for k in range(7))
         and all(J[i][j] > 0 for end in ends for i in end for j in end),
         'complete strict positive inverse end blocks and actual original short8 distance split')
    corner = set(white) | set(key['selected_marks'])
    need(set(key['base_nodes']) <= corner
         and p['omitted_black'] == sorted(set(u+v)-set(key['selected_marks']))
         and all(g[i][j] == 0 for i in corner-set(key['base_nodes']) for j in key['base_nodes']),
         'genuinely new first entry: omitted original blacks absent and other old corner blocks orthogonal')
    image = dict(zip(key['induced_map'],key['base_nodes']))
    missing = sorted(set(range(key['target_order']))-set(image))
    entries = [dict(new_target_node=new,old_target_node=old,core_node=image[old],
                    native_label=label)
               for new in missing for old,label in target[new].items() if old in image]
    need(entries and all(e['native_label'] in (3,4) for e in entries),
         'nonempty exact mandatory native first-entry domain')
    return dict(p=p,core=core,g=g,J=J,target=target,white=white,allowed=allowed,
                ends=ends,du=du,dv=dv,missing=missing,entries=entries,
                physical=L.context(core['gram'],white))


def b_case(doc,field,key):
    matches = [(i,p) for i,p in enumerate(doc[field]) if p['index'] == key['high_index']]
    if not matches:
        return None
    return one(matches,lambda t:True,'unique exact affine-B saved source occurrence')


def b_map(rows,key):
    return one([[i,p] for i,p in enumerate(rows)],lambda t:t[1]['map'] == key['induced_map'],
               'THIS complete affine-B native map in the exact original proof domain')


def b_bound(key,frame,src,selector,map_selector,seal,reason,producer,report):
    seal = dict(seal,original_high_source=L.HIGH,
        original_high_selector=['isolation','cases',key['high_index']],
        exact_core_selector=['isolation','cores',key['core_key']],
        source_signature=frame['core']['signature'],source_gauge=frame['core']['gauge'])
    return dict(application_kind='EXCLUDED_COMPLETION_UNDER_NO_S8',
        source=src,proof_selector=selector,map_source=src,map_selector=map_selector,
        physical=frame['physical'],source_sealing=seal,paying_corners=[],actual_receipts=[],
        source_branch_conclusion=dict(kind='THIS_EXACT_SELECTED_COMPLETION_IMPOSSIBLE_UNDER_NO_S8',
            reason=reason,source_numerical_predicates_runtime_regenerated=False),
        proof_references=[pin(B_GENERAL_REPORT),pin(producer),pin(report),pin(B_CONTACT_PRODUCER),
                          pin(B_RESIDUAL_REPORT),pin(B_AFFINE_PROOF),pin(B_AFFINE_ROOT)])


def b14_saved_application(key,W,frame):
    if key['target_order'] != 14 or key['marked_image'] not in ([9,12],[9,13]):
        return None
    doc = source(B14_SOURCE);found = b_case(doc,'cases',key)
    if found is None:
        return None
    pi,p = found;mi,m = b_map(p['maps'],key)
    support = one(doc['reconstructed_support'],lambda r:r['index'] == key['high_index'],
                  'same complete B14 support occurrence')
    need(support['embeddings'].count(key['induced_map']) == 1
         and key['marked_image'] in support['marked_images']
         and p['core_key'] == key['core_key'] and p['method'] == frame['p']['method']
         and p['base_nodes'] == key['base_nodes'] and p['preferred_corner'] == key['selected_marks']
         and p['core_gram'] == frame['core']['gram'] and p['core_inverse'] == frame['core']['inverse'],
         'same B14 complete support/full map/core/physical inverse/ordered marks')
    need(p['end_projection_proof']['allowed_contacts'] == frame['allowed']
         and p['end_projection_proof']['end_sets'] == frame['ends'],
         'same B14 complete original end/contact proof')
    entries = m['forced_entries']
    need([ {k:e[k] for k in ('new_target_node','old_target_node','core_node','native_label')}
           for e in entries] == frame['entries'],
         'same exact source B14 forced first entries')
    seal = dict(kind='accepted_complete_B14_saved_map_exclusion',
        source_selector=['cases',pi,'maps',mi],complete_end_proof=p['end_projection_proof'],
        forced_entries=entries,saved_numeric_predicates_runtime_regenerated=False)
    if m['outcome'] == 'first_entry_projection_above_norm':
        need(m['witness'] in entries and rational(m['witness']['projection_lower_bound']) > 2,
             'named physical B14 first-entry projection obstruction')
        seal['forced_projection_witness'] = m['witness']
        reason = 'Required physical norm2 first entry exceeds its projection norm bound, unless this same face has the independent short8 escape.'
    else:
        need(m['outcome'] == 'doubled_terminal_residual_cauchy_contradiction'
             and frame['missing'] == [10,11] and len(entries) == 1
             and entries[0]['new_target_node'] == m['first_new_target_node'] == 10
             and entries[0]['old_target_node'] == 9 and entries[0]['native_label'] == 3
             and entries[0]['core_node'] == m['old_entry_node']
             and m['doubled_terminal_target_node'] == 11
             and m['native_terminal_bond'] == frame['target'][10][11] == 4,
             'exact source B14 one simple entry followed by a genuinely new doubled terminal')
        i = m['old_entry_node'];rho = rational(m['first_entry_residual_square'])
        need(0 < rho < 1 and rho == 2-frame['J'][i][i]
             and rational(m['terminal_physical_norm']) == 2
             and rational(m['required_residual_inner_product_square_at_least']) == 2
             and rational(m['Cauchy_upper_bound']) == 2*rho < 2
             and rational(m['strict_gap']) == 2-2*rho,
             'same source strict residual Cauchy contradiction')
        signs = m['terminal_complete_contact_signs']
        need([r['core_node'] for r in signs] == frame['allowed']
             and all(rational(r['projection_coefficient']) == frame['J'][i][r['core_node']]
                     and (r['disposition'] == 'nonnegative_projection' and frame['J'][i][r['core_node']] >= 0
                          or r['disposition'] == 'opposite_black_short_path'
                          and frame['J'][i][r['core_node']] < 0
                          and type(r['path_length']) is int and 0 < r['path_length'] <= 8)
                     for r in signs),
             'COMPLETE terminal contact coefficients: every negative term has the same-face short8 escape')
        seal.update(first_entry_residual_square=m['first_entry_residual_square'],
            same_end_two_contact_exclusions=m['same_end_two_contact_exclusions'],
            terminal_complete_contact_signs=signs,terminal_is_new_reason=m['terminal_is_new_reason'],
            required_residual_inner_product_square_at_least='2',
            Cauchy_upper_bound=m['Cauchy_upper_bound'],strict_gap=m['strict_gap'])
        reason = 'Under NO-S8 the arbitrary-contact terminal residual has squared product at least2, but positive-complement Cauchy bounds it by2rho<2.'
    return b_bound(key,frame,B14_SOURCE,['cases',pi],['cases',pi,'maps',mi],
                   seal,reason,B14_PRODUCER,B14_REPORT)


def b15_saved_application(key,W,frame):
    if key['target_order'] != 15 or key['marked_image'] not in ([9,13],[9,14]):
        return None
    doc = source(B15_SOURCE);found = b_case(doc,'B_cases',key)
    if found is None:
        return None
    pi,case = found;p = case['proof'];mi,m = b_map(p['maps'],key)
    support = one(doc['B_reconstructed_support'],lambda r:r['index'] == key['high_index'],
                  'same complete B15 reconstructed support occurrence')
    need(support['embeddings'].count(key['induced_map']) == 1
         and key['marked_image'] in support['marked_images']
         and case['core_key'] == key['core_key'] and case['method'] == frame['p']['method']
         and p['base_nodes'] == key['base_nodes']
         and p['core_gram'] == frame['core']['gram'] and p['core_inverse'] == frame['core']['inverse'],
         'same B15 complete support/full map/core/physical inverse/ordered actual high row')
    need(p['end_projection_proof']['allowed_contacts'] == frame['allowed']
         and p['end_projection_proof']['end_sets'] == frame['ends']
         and m['missing_target_nodes'] == frame['missing'],
         'same B15 complete old-core contacts/ends and native missing-node domain')
    entries = m['forced_entries']
    need([{k:e[k] for k in ('new_target_node','old_target_node','core_node','native_label')}
          for e in entries] == frame['entries'],'same exact B15 mandatory physical native entries')
    seal = dict(kind='accepted_complete_B15_saved_map_exclusion',
        source_selector=['B_cases',pi,'proof','maps',mi],complete_end_proof=p['end_projection_proof'],
        forced_entries=entries,saved_numeric_predicates_runtime_regenerated=False)
    if m['outcome'] == 'forced_entry_projection_exclusion':
        need(m['witness'] in entries and rational(m['witness']['projection_lower_bound']) > 2,
             'named physical B15 first-entry projection obstruction')
        seal['forced_projection_witness'] = m['witness']
        reason = 'A mandatory genuinely new physical norm2 first entry exceeds its projection bound under NO-S8.'
    else:
        need(m['outcome'] == 'forced_completed_index_two' and len(entries) == 1
             and entries[0]['native_label'] == 3
             and frame['missing'] in ([11,12],[10,11,12])
             and entries[0]['new_target_node'] == min(frame['missing']),
             'exact source B15 residual/no-bridge forced-tail branch')
        control = m['residual_contact_control']
        need(control['allowed_contacts'] == frame['allowed'] and control['end_sets'] == frame['ends']
             and set(control['candidate_nodes']) ==
                 {i for i in frame['allowed'] if frame['J'][i][i] <= 2}
             and all(rational(control['projection_diagonals'][str(i)]) == frame['J'][i][i]
                     for i in frame['allowed']),
             'same complete necessary contact and native single-contact candidate domain')
        if frame['p']['method'] == 'lorentz_residual':
            need(control['kind'] == 'at_most_one_core_contact_root'
                 and [r['nodes'] for r in control['pair_tests']] ==
                     [[i,j] for i in control['candidate_nodes'] for j in control['candidate_nodes']]
                 and m['exterior_root_order_exclusion'] is None,
                 'complete ordered residual pair tests, including repeated coordinates')
        else:
            exterior = m['exterior_root_order_exclusion']
            need(frame['p']['method'] == 'relevant_component'
                 and control['kind'] == 'one_relevant_root_with_optional_exterior_root'
                 and control['relevant_node'] == entries[0]['core_node']
                 and exterior['actual_factor_order'] == 15
                 and exterior['whole_order_with_exterior_root'] == control['joint_component_order_bound'] < 15,
                 'optional exterior-root NO-BRIDGE order bound applies to THIS whole order15 target')
        congruence = m['terminal_congruence']
        need(congruence == dict(target_node=12,physical_norm='2',scale='1/sqrt2',
             scaled_norm='1',scaled_native4_bond='-1')
             and m['signature'] in ([16,2,0],[17,2,0])
             and len(m['gram']) == len(frame['g'])+len(frame['missing'])
             and len(m['target_to_literal']) == len(set(m['target_to_literal'])) == 15
             and len(m['positive_affine_kernel']) == 15
             and all(rational(v) > 0 for v in m['positive_affine_kernel']),
             'exact saved physical-terminal congruence and two-negative completed principal')
        need([m['target_to_literal'][t] for t in key['induced_map']] == key['base_nodes'],
             'same full native image in the literal completed principal')
        seal.update(residual_contact_control=control,
            exterior_root_order_exclusion=m['exterior_root_order_exclusion'],
            terminal_congruence=congruence,target_to_literal=m['target_to_literal'],
            completed_gram_source_selector=['B_cases',pi,'proof','maps',mi,'gram'],
            accepted_saved_completed_signature=m['signature'],
            accepted_saved_positive_affine_kernel=m['positive_affine_kernel'],
            completion_zero_geometry='Only after complete at-most-one-contact or strict exterior-root/no-bridge exclusion: every later required target white is genuinely new and old-core-orthogonal.')
        reason = 'The accepted complete residual contact/no-bridge proof forces THIS saved principal Gram before zeros are written; its declared two-negative inertia is impossible in the ambient Lorentz space.'
    return b_bound(key,frame,B15_SOURCE,['B_cases',pi],['B_cases',pi,'proof','maps',mi],
                   seal,reason,B15_PRODUCER,B15_REPORT)


B_EXACT_CORE_CACHE = set()
B_GJ_PRODUCTS = 0
B_MAX_GJ_PRODUCTS = 2_000_000


def b_checked_exact_core(key,frame):
    """Exact given-core identity only; no native-map or new-root enumeration."""
    global B_GJ_PRODUCTS
    core = frame['core'];g,J = frame['g'],frame['J'];n = len(g)
    digest = L.sha256(canonical([core['gram'],core['inverse'],core['gauge'],core['signature']]).encode()).hexdigest()
    cache_key = (key['core_key'],digest)
    if cache_key not in B_EXACT_CORE_CACHE:
        if cache_key in L.CORE_CACHE:
            old_g,old_J,scales = L.CORE_CACHE[cache_key]
            need(old_g == g and old_J == J and scales == [F(1)]*n,
                 'same exact core identity already checked by frozen legacy helper in THIS process')
        else:
            for i in range(n):
                for j in range(n):
                    value = F()
                    for k in range(n):
                        if g[i][k]:
                            B_GJ_PRODUCTS += 1
                            need(B_GJ_PRODUCTS <= B_MAX_GJ_PRODUCTS,
                                 'complete exact affine-B inverse identity budget; refuse, never truncate')
                            value += g[i][k]*J[k][j]
                    need(value == int(i == j),'exact physical GJ identity on THIS complete core input')
        B_EXACT_CORE_CACHE.add(cache_key)
    return digest


def b_orthogonal_terminal_cauchy(key,W,frame):
    """The successor's native zeros and complete contacts prove them physically."""
    if set(key['base_nodes']) != set(frame['white']) | set(key['selected_marks']):
        return None  # Otherwise an unassigned original corner factor needs separate control.
    J = frame['J'];target = frame['target']
    outside = [i for i in frame['allowed'] if i not in key['base_nodes']]
    if not all(J[i][i] > 2 for i in outside):
        return None
    for terminal in frame['missing']:
        if len(target[terminal]) != 1:
            continue
        first,label = next(iter(target[terminal].items()))
        if first not in frame['missing'] or label not in (3,4):
            continue
        for entry in frame['entries']:
            if entry['new_target_node'] != first or entry['native_label'] != 3:
                continue
            i = entry['core_node']
            if i not in frame['allowed'] or not 1 < J[i][i] <= 2:
                continue
            rho = 2-J[i][i];bond_square = F(1 if label == 3 else 2)
            if not bond_square > 2*rho:
                continue
            digest = b_checked_exact_core(key,frame)
            return dict(application_kind='EXCLUDED_COMPLETION_UNDER_NO_S8',
                source=L.HIGH,proof_selector=['isolation','cases',key['high_index']],
                map_source=B_GENERAL_REPORT,map_selector=dict(
                    format='GIVEN_COMPLETE_NATIVE_MAP_AND_STRICT_OLD_CORNER_TERMINAL_CERTIFICATE',
                    induced_map=key['induced_map'],first=first,terminal=terminal),
                physical=frame['physical'],source_sealing=dict(
                    kind='strict_old_corner_covered_and_new_terminal_orthogonal_Cauchy',
                    core_input_sha256=digest,checked_in_this_actual_process_on_exact_core=True,
                    source_signature=frame['core']['signature'],source_gauge=frame['core']['gauge'],
                    complete_original_actual_corner=key['base_nodes'],
                    omitted_original_black=frame['p']['omitted_black'],
                    complete_allowed_contacts=frame['allowed'],complete_end_sets=frame['ends'],
                    terminal_target_node=terminal,terminal_only_native_neighbor=first,
                    terminal_has_zero_product_with_every_mapped_old_corner_root=True,
                    all_possible_contacts_outside_the_mapped_corner=[
                        dict(core_node=j,inverse_diagonal=str(J[j][j])) for j in outside],
                    terminal_old_core_orthogonality_reason='All original actual corner roots are already mapped, so the distinct missing terminal is new. Its native zeros exclude all mapped roots. Any nonzero remaining old contact has inverse diagonal>2 and native magnitude>=1 on a positive inverse end; both-end contacts give independent S8. Therefore every old-core product is zero.',
                    first_target_node=first,forced_old_target_neighbor=entry['old_target_node'],
                    forced_old_core_contact=i,first_native_label=3,
                    first_projection_square_lower=str(J[i][i]),
                    first_residual_square_upper=str(rho),
                    arbitrary_additional_first_root_contacts_not_set_zero=True,
                    terminal_residual_square='2',successor_native_label=label,
                    required_residual_product_square=str(bond_square),
                    Cauchy_upper_bound=str(2*rho),strict_gap=str(bond_square-2*rho),
                    source_contact_and_signature_completeness_runtime_regenerated=False),
                paying_corners=[],actual_receipts=[],
                source_branch_conclusion=dict(kind='THIS_EXACT_SELECTED_COMPLETION_IMPOSSIBLE_UNDER_NO_S8',
                    reason='The genuinely new native terminal is forced old-core-orthogonal before any product is set zero. Its required bond to the first root violates positive-complement Cauchy despite arbitrary other contacts of that first root.'),
                proof_references=[pin(B_GENERAL_REPORT),pin(B14_PRODUCER),pin(B14_REPORT),
                    pin(B_AFFINE_PROOF),pin(B_AFFINE_ROOT)])
    return None


def b_native_successor_cauchy(key,W,frame):
    """Strict typed extension of the existing two-root/saturated-entry proof."""
    J = frame['J'];target = frame['target'];allowed = frame['allowed'];ends = frame['ends']
    for first in frame['missing']:
        entries = [e for e in frame['entries'] if e['new_target_node'] == first]
        if len(entries) != 1 or entries[0]['native_label'] != 3:
            continue
        i = entries[0]['core_node']
        if i not in allowed or not 1 < J[i][i] <= 2:
            continue
        end = next(e for e in ends if i in e)
        # A forced simple entry already has magnitude1. Every additional same-
        # end native contact has magnitude>=1 and would exceed norm2 here.
        if not all(J[i][i]+J[j][j]+2*J[i][j] > 2 for j in end if j != i):
            continue
        signs = []
        for j in allowed:
            a = J[i][j]
            length = min(frame['du'][i]+frame['dv'][j],frame['du'][j]+frame['dv'][i])+3
            if a < 0 and length > 8:
                break
            signs.append(dict(core_node=j,projection_coefficient=str(a),
                disposition='same_face_white_interior_short8' if a < 0 else 'nonnegative_projection',
                path_length=length if a < 0 else None))
        else:
            rho = 2-J[i][i]
            for second,label in target[first].items():
                if second not in frame['missing'] or label not in (3,4):
                    continue
                bond_square = F(1 if label == 3 else 2)
                if not bond_square > 2*rho:
                    continue
                digest = b_checked_exact_core(key,frame)
                seal = dict(kind='strict_native_first_entry_and_successor_residual_Cauchy',
                    original_high_selector=['isolation','cases',key['high_index']],
                    complete_core_selector=['isolation','cores',key['core_key']],
                    core_input_sha256=digest,checked_in_this_actual_process_on_exact_core=True,
                    source_signature=frame['core']['signature'],source_gauge=frame['core']['gauge'],
                    complete_allowed_contacts=allowed,complete_end_sets=ends,
                    first_target_node=first,old_target_node=entries[0]['old_target_node'],
                    forced_old_core_contact=i,first_native_label=3,first_physical_contact='1',
                    inverse_diagonal=str(J[i][i]),first_residual_square=str(rho),
                    first_complete_same_end_extra_contact_lowers=[
                        dict(other_core_node=j,lower=str(J[i][i]+J[j][j]+2*J[i][j]))
                        for j in end if j != i],
                    successor_target_node=second,successor_native_label=label,
                    required_successor_bond_square=str(bond_square),
                    successor_is_genuinely_new='The first root has exactly one old-core contact at an already mapped base node. An old missing successor in another original corner block would be an additional old-core contact; omitted original blacks are absent.',
                    complete_successor_contact_coefficients=signs,
                    successor_projection_square_nonnegative='Every allowed one-end inverse block is entrywise positive; a core-orthogonal successor has projection zero.',
                    successor_residual_square_upper='2',
                    required_residual_product_square_lower=str(bond_square),
                    Cauchy_upper_bound=str(2*rho),strict_gap=str(bond_square-2*rho),
                    independent_same_face_S8_premise_explicit=True,
                    source_contact_and_signature_completeness_runtime_regenerated=False)
                return dict(application_kind='EXCLUDED_COMPLETION_UNDER_NO_S8',
                    source=L.HIGH,proof_selector=['isolation','cases',key['high_index']],
                    map_source=B14_REPORT,map_selector=dict(
                        format='GIVEN_COMPLETE_NATIVE_MAP_AND_VERIFIED_STRICT_RESIDUAL_CERTIFICATE',
                        induced_map=key['induced_map'],first=first,successor=second),
                    physical=frame['physical'],source_sealing=seal,
                    paying_corners=[],actual_receipts=[],
                    source_branch_conclusion=dict(
                        kind='THIS_EXACT_SELECTED_COMPLETION_IMPOSSIBLE_UNDER_NO_S8',
                        reason='NO-S8 removes every negative old-projection term. The required native successor bond forces residual product square above positive-complement Cauchy, including saturated rho=0.'),
                    proof_references=[pin(B_GENERAL_REPORT),pin(B14_PRODUCER),pin(B14_REPORT),
                        pin(B_CONTACT_PRODUCER),pin(B_RESIDUAL_REPORT),
                        pin(B_AFFINE_PROOF),pin(B_AFFINE_ROOT)])
    return None


def b_apply(key,W,units):
    L.full_key(key)
    if key['affine'] and key['target_kind'] == 'B' and key['target_order'] == 13:
        return b13_full_application(key,W,units)
    frame = b_source_frame(key,W)
    if frame is None:
        return None
    need(W.affine_weight('B',key['target_order'],*key['marked_image']) == 0,
         'CURRENT selected affine-B intrinsic zero remains zero')
    bound = None
    if key['target_order'] == 14:
        bound = b14_saved_application(key,W,frame)
    elif key['target_order'] == 15:
        bound = b15_saved_application(key,W,frame)
    if bound is None:
        bound = b_orthogonal_terminal_cauchy(key,W,frame)
    return bound if bound is not None else b_native_successor_cauchy(key,W,frame)


# Additive terminal-tail successor. Existing B/D functions and contracts stay intact.
B_TERMINAL_TAIL_NOTE = 'results/finite57_affine_b_longer_terminal_tail_source_proposal_20261025.md'
PINS[B_TERMINAL_TAIL_NOTE] = 'c2cc35cb9e0739845ac0ead4b4b9ac6f74f227bb3c578f47a65bc8e42acda14f'


def b_terminal_tail_contact_profile(key,frame):
    """An old CONTACT theorem is transported; no old B15 map is relabelled."""
    p,core = frame['p'],frame['core'];matches = []
    for path,field in ((D16_SOURCE,'cases'),(B15_SOURCE,'B_cases')):
        doc = source(path)
        need(doc['sources'].get(L.HIGH) == L.PINS[L.HIGH],
             'contact theorem was produced from THIS exact complete original HIGH source')
        for ri,row in enumerate(doc[field]):
            if row['index'] != key['high_index']:
                continue
            if path == D16_SOURCE:
                proof = row
                need(row['preferred_corner'] == key['selected_marks'],
                     'contact profile uses THIS original ordered actual corner')
                controls = [(['constraints'],row['constraints'])]
            else:
                proof = row['proof']
                controls = [(['proof','maps',mi,'residual_contact_control'],m['residual_contact_control'])
                            for mi,m in enumerate(proof['maps']) if 'residual_contact_control' in m]
            need(row['core_key'] == key['core_key'] and row['method'] == p['method']
                 and proof['core_gram'] == core['gram'] and proof['core_inverse'] == core['inverse']
                 and proof['base_nodes'] == key['base_nodes'],
                 'same complete original high/core/base/physical inverse CONTACT theorem')
            for suffix,control in controls:
                matches.append(dict(source=path,selector=[field,ri]+suffix,control=control))
    if not matches:
        return None
    need(all(canonical(q['control']) == canonical(matches[0]['control']) for q in matches),
         'all saved complete contact profiles for this original occurrence agree exactly')
    chosen = min(matches,key=lambda q:canonical([q['source'],q['selector']]))
    # This validates the COMPLETE ordinary/relevant numerical hypotheses,
    # including all contact pairs, repeated coordinates and exterior no-bridge.
    g,J = d16_control(key,p,core,chosen['control'])
    need(g == frame['g'] and J == frame['J'],'same exact physical core in both contact frames')
    return chosen,matches


def b_complete_terminal_tail_application(key,W,units):
    if not key['affine'] or key['target_kind'] != 'B' or key['target_order'] not in (16,17):
        return None
    frame = b_source_frame(key,W)
    if frame is None or W.affine_weight('B',key['target_order'],*key['marked_image']) != 0:
        return None
    N = key['target_order'];missing = frame['missing']
    if len(missing) < 2:
        return None
    first = missing[0]
    if (missing != list(range(first,N-2)) or
            set(key['induced_map']) != set(range(first)) | {N-2,N-1}):
        return None
    graph = frame['target']
    image = dict(zip(key['induced_map'],key['base_nodes']))
    if first == 0:
        return None
    entry = dict(new_target_node=first,old_target_node=first-1,
                 core_node=image[first-1],native_label=3)
    need(frame['entries'] == [entry],
         'THIS given map has exactly its required simple first native spine entry')
    for t in missing:
        expected = {}
        if t == first:
            expected[first-1] = 3
        else:
            expected[t-1] = 4 if t == N-3 else 3
        if t < N-3:
            expected[t+1] = 4 if t+1 == N-3 else 3
        need(graph[t] == expected,
             'complete missing suffix is an unbranched native path with ONLY its final doubled bond')
    profile = b_terminal_tail_contact_profile(key,frame)
    if profile is None:
        return None
    chosen,all_profiles = profile;control = chosen['control'];i = entry['core_node']
    if i not in control['candidate_nodes'] or not 1 < frame['J'][i][i] <= 2:
        return None
    exterior = None
    if control['kind'] == 'one_relevant_root_with_optional_exterior_root':
        need(control['relevant_node'] == i and
             control['joint_component_order_bound'] < N,
             'THIS forced first root is relevant and the actual whole target excludes its exterior root')
        exterior = dict(actual_whole_target_order=N,
            complete_source_whole_order_with_exterior_root=control['joint_component_order_bound'],
            exterior_node=control['exterior_node'],
            no_bridge_scope='Full inherited residual positive-relation/saturated-tail theorem, including exterior multiplicity and all indirect supporting-white bridges; not immediate-neighbor counting.')
    else:
        need(control['kind'] == 'at_most_one_core_contact_root' and
             frame['p']['method'] == 'lorentz_residual',
             'complete ordinary uniqueness theorem controls every additional old-core contact')
    core_digest = b_checked_exact_core(key,frame)
    rho = 2-frame['J'][i][i]
    need(0 <= rho < 1,'strict positive-complement first-root residual bound')
    tail = missing[1:]
    # A terminal-only positive congruence records the EXACT vector identity:
    # old physical norms stay2; terminal coefficient1/sqrt2 gives diagonal1.
    scaled_tail = [[F(2 if a == b else -1 if abs(a-b) == 1 else 0)
                    for b in range(len(tail))] for a in range(len(tail))]
    scaled_tail[-1][-1] = F(1)
    need(sum((x for row in scaled_tail for x in row),F()) == 1,
         'exact terminal-scaled tail vector norm1')
    coefficients = ['1']*(len(tail)-1)+['1/sqrt(2)']
    return dict(application_kind='EXCLUDED_COMPLETION_UNDER_NO_S8',
        source=chosen['source'],proof_selector=chosen['selector'],
        map_source=L.HIGH,map_selector=dict(
            format='GIVEN_COMPLETE_CURRENT_NATIVE_MAP_ON_ORIGINAL_HIGH',
            original_high_selector=['isolation','cases',key['high_index']],
            target=dict(affine=True,kind='B',order=N),induced_map=key['induced_map'],
            old_contact_theorem_is_not_current_target_map_membership=True),
        physical=frame['physical'],source_sealing=dict(
            kind='COMPLETE_SAVED_CONTACT_THEOREM_AND_EXACT_NATIVE_B_TERMINAL_TAIL',
            all_matching_original_contact_profile_selectors=[
                dict(source=q['source'],selector=q['selector']) for q in all_profiles],
            complete_contact_control=control,first_native_entry=entry,
            whole_target_order=N,strict_exterior_exclusion=exterior,
            complete_original_possible_contacts=frame['allowed'],complete_original_end_sets=frame['ends'],
            source_signature=frame['core']['signature'],source_gauge=frame['core']['gauge'],
            core_input_sha256=core_digest,checked_in_this_actual_process_on_exact_core=True,
            required_complete_missing_suffix=missing,
            later_root_newness_and_orthogonality='Only AFTER complete per-root contact/ordinary uniqueness or actual-whole-order exterior/no-bridge exclusion: every later missing target root is new and old-core-orthogonal. An old root in another original corner block would create a second old-core contact at the first bridge.',
            first_projection_square=str(frame['J'][i][i]),first_residual_square=str(rho),
            physical_tail_target_nodes=tail,physical_tail_root_norms=['2']*len(tail),
            vector_coefficients=coefficients,
            terminal_only_congruent_tail_gram=[[str(x) for x in row] for row in scaled_tail],
            vector_norm_square='1',first_root_vector_inner_product='-1',
            Cauchy_required_square='1',Cauchy_upper_bound=str(rho),strict_gap=str(1-rho),
            new_inertia_computation_used=False,
            source_signature_contact_completeness_and_actual_Lorentz_geometry_regenerated=False),
        paying_corners=[],actual_receipts=[],
        source_branch_conclusion=dict(kind='THIS_EXACT_ACTUAL_COMPLETION_IMPOSSIBLE_UNDER_NO_S8',
            reason='The full given-map tail is proved new and old-core-orthogonal before zeros. Its physical terminal-scaled vector has norm1 and product-1 with the first residual of norm square<1, contradicting positive-space Cauchy.'),
        proof_references=[pin(B_TERMINAL_TAIL_NOTE),pin(AFFINE_REVIEW),pin(AFFINE_ROOT),
                          pin(RESIDUAL),pin(B_CONTACT_PRODUCER)])


_B_PRE_TERMINAL_TAIL_APPLY = b_apply


def b_apply(key,W,units):
    previous = _B_PRE_TERMINAL_TAIL_APPLY(key,W,units)
    return previous if previous is not None else b_complete_terminal_tail_application(key,W,units)


PINS.update(B_SOURCE_PINS)


def apply_selected(key,W,units):
    global L
    L=LEGACY_HELPER
    need(L is not None and L.APPLICATION_SOURCE==LEGACY,'injected unchanged legacy helper definitions')
    pin(LEGACY,LEGACY_SHA);L.full_key(key)
    pin(L.TEX,L.TEX_SHA);pin(L.LITERAL,L.LITERAL_SHA)
    need(getattr(W,'TEX_SHA256',None)==L.TEX_SHA,'the frozen current literal/canonical manuscript')
    if not key['affine'] or key['target_kind'] not in ('B','D'):return None
    bound=None
    if key['target_kind']=='D':
        if key['target_order']==17:
            bound=d17_whole(key,W,units) if len(key['base_nodes'])==17 else d17_proper(key,W,units)
        else:bound=d16_saved(key,W,units)
    elif key['target_kind']=='B':bound=b_apply(key,W,units)
    if bound is None:return None
    application_kind=bound.get('application_kind','EXACT_CURRENT_RECEIVER_SOURCE_BINDING')
    excluded=application_kind=='EXCLUDED_COMPLETION_UNDER_NO_S8'
    ids=['finite57.oldAffine.exact_source_maps','finite57.oldAffine.complete_contact_completion']
    if not excluded:ids.append('finite57.oldAffine.current_whole_receipts')
    accepted=acceptance(ids);frozen=json.loads(canonical(key))
    return dict(status='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION',
        application_kind=application_kind,
        theorem_id='finite57.old_affine_source.'+key['target_kind']+str(key['target_order'])+'.current.v1',
        application_key=dict(selected_full_key=frozen,source=bound['source'],
            source_proof_selector=bound['proof_selector'],source_map_selector=bound['map_selector']),
        physical_context=bound['physical'],source_sealing=bound['source_sealing'],
        actual_terminal_receipts=bound['actual_receipts'],current_paying_corners=bound['paying_corners'],
        current_receiver_binding=bound.get('current_receiver_binding'),
        consumer_literal_evaluations=bound.get('consumer_literal_evaluations',[]),
        source_branch_conclusion=bound.get('source_branch_conclusion',
            dict(kind='CURRENT_DISTINCT_CORNER_TOTAL_AT_LEAST_ONE_UNDER_NO_S8')),
        source_references=[pin(bound['source']),pin(bound['map_source'],
            L.PINS[bound['map_source']] if bound['map_source'] not in PINS else None)],
        proof_references=bound['proof_references']+[pin(CONTRACT)],supplemental_source_premise_acceptance=accepted,
        supplemental_source_premises=[dict(id=i,statement=PREMISES[i],runtime_regenerated=False,
            source_selector=bound['proof_selector'],same_full_key=frozen,accepted_by=accepted) for i in ids],
        assumptions=[dict(kind='same_actual_source_face_and_whole_target',selected_full_key=frozen),
            dict(kind='NO_SAME_FACE_ADJACENT_WHITE_INTERIOR_PATH_AT_MOST_8',runtime_verified=False),
            dict(kind='complete_source_numerical_and_geometric_statements',runtime_regenerated=False,
                 source=bound['source'],selector=bound['proof_selector'])],
        provider_composition_premise=dict(status='EXACT_SOURCE_NO_S8_IMPLICATION_BOUND_TO_THIS_FULL_KEY',
            same_full_key=frozen,implication='The accepted exhaustive source implication at its complete retained predicates excludes this completion or pays this SAME actual face under NO-S8.',
            independent_S8_composition='If the same face were unpaid, independent current S8 forbids every such path; the exact source implication then gives a contradiction.',
            requiresS8_flag_alone_is_not_a_disjunction=True),
        conclusion=dict(quantity='total_current_charge_of_same_actual_face',lower='1'),
        intrinsic_corner_unit=False,actual_short_path_asserted=False,theorem_execution=False,
        geometric_formalization=False,full_transitive_regeneration=False,global_finite57_verified=False)


if __name__=='__main__':
    raise SystemExit('Pure application definitions only; root must load exact reviewed source into F')
