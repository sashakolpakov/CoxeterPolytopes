#!/usr/bin/env python3
"""Pure full-domain affine-D applications at accepted current premises.

No repository import, file I/O, native-map enumeration or standalone job.
F supplies byte-pinning readers, its actual-validated inventory, current W
and the already checked independent-S8 runtime record. Shared expanded
contracts MUST be retained once by the enclosing actual selected receipt.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
import re

SOURCE_PINNER = None
SOURCE_READER = None
SOURCE_ACCEPTANCE = None
APPLICATION_SOURCE_SHA = None
OWN = 'scripts/finite57_full_universal_affine_d_applications_20261025.py'
PREMISE_IDS = ('finite57.fullD.complete_profiles','finite57.fullD.current_boundary_roles',
               'finite57.fullD.whole_face_scope')
CACHE = {}
CONTRACTS = {}
FIELDS = {'high_index','core_key','base_nodes','selected_marks','target_kind',
          'target_order','affine','induced_map','marked_image'}
PREFIX = 'finite57.full_universal_affine_D.current.20261025'
LEAF = 'results/fresh_affine_d12_leaf_pair_transfer_20261017.json'
D14 = 'results/fresh_affine_d13_long_pair_transfer_20261017.json'
D15 = 'results/fresh_affine_d14_long_pair_transfer_20261017.json'
BOUNDARY = 'results/fresh_d12_universal_ideal_output_repairs_20261018.json'
MODEL = 'results/fresh_d13_three_rank15_context_transfer_20261018.json'
LOCAL_D13 = 'results/fresh_order13_d_center_leaf_half_transfer_20261018.json'
DIRECT = 'scripts/check_finite57_terminal_collectors_direct_receivers_20261025.py'
FINITE = 'scripts/check_finite57_terminal_collectors_syntax_safe_20261025.py'
UNIVERSAL = 'scripts/check_finite57_current_universal_receivers_20261025.py'
MICRO = 'scripts/check_finite57_final_quiet_cases_20261025.py'
QUIET = 'scripts/check_finite57_quiet_b13_current_priority_bridge_20261025.py'
ONE = 'scripts/check_finite57_final_guard_one_sided_current_receivers_20261025.py'
S8 = 'scripts/check_finite57_current_s8_applications_20261025.py'
G4 = 'scripts/check_finite57_g4_lower_cycle_terminal_bindings_20261025.py'
ZERO = 'scripts/check_finite57_current_d13_center_zero_outputs_20261025.py'
CASE_INDICES = tuple(3+6*i+j for i in range(8) for j in range(3))
PINS = {
 LEAF:'8f24d38cca6526a25652835c04923f542338590509c1199bcba2707520f62eb5',
 D14:'79989c70ace2cc1591d88e4167de39dc243977b748a6d369a495c313945c8a2a',
 D15:'42fea3ab62fd837ed3652a289d889107f6933d323096c7961c971ab28cd023fd',
 BOUNDARY:'82850cfdfe1ff1a9d56f4f76f5466d02d7539d74670e4172aa3f7c7643a0916d',
 MODEL:'2c7d615cd38123bef20cffad68de8ab30baa91826c147a6bb084892590e8ca86',
 LOCAL_D13:'02c74a83e76be9b903e1013228cba6de72beeda4bbb743e2f6687a548d314465',
 DIRECT:'a2a9d61a7637e45d5a4d7af1c3283e120ca7685d07e9f686ae7b6ee585520eca',
 FINITE:'930ffff5601e6096359d3b2f8846fa6d7d050f95fb2aabe449447e0d3d0bf9ca',
 UNIVERSAL:'29ddbcf65cae2a7c7ff7703a1b4a81bad3ca1fdb0ade4bdbffda26fb9c8a1631',
 MICRO:'9e4e56ce9b290e0b660f8b3c7c87442a986c9d01a1da93ed40bd6b8f36d335b4',
 QUIET:'5e5c1b9736701dbb48ff3f3dc6bba6b7f87a314758d5fe0b16802607044263e4',
 ONE:'0892b27a233d761b9568067a3630c73b63f25e551e20b4fd9b69ae6ba99425a3',
 S8:'e9a8074d79b6861f04678d10fd9abe20259a36d7c7e04ed21f815836b81c5755',
 G4:'f17014594587a68ef766ae245db372550fe89adc02bedf9e44c477a8d79eb242',
 ZERO:'08b011fb62f9c02f5eea1b9008f9fc865b4f3e4adf1ad84b8ceb22b73722cf3f',
 'scripts/finite57_literal_weights_audit_20261024.py':'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0',
 'manuscripts/finite_volume57.tex':'c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0',
 'results/audit_finite57_nine_fork_profiles_and_twelve_d13_repairs_20261024.md':'0ecf41083477787cfd20af8dc899138e1480db0aec7bf5ce6115838ef6631146',
 'results/review_finite57_remaining_universal_profiles_root_20261025.md':'a02e9daa3e3f42a3d3b832bc561702b87e658ad95693644bf7469fe91c1b6201',
 'results/finite57_current_universal_receiver_retention_bridge_20261025.md':'309bbd85a766466bcb3bb4437d7a5d50ffe8c614c4892f40f4ccc06d842893eb',
 'results/finite57_current_universal_receiver_bridge_root_acceptance_20261025.json':'6bc292fb08a50c15cfa1cff2ef9e43def44d4febd942904f0cdd1033eac4838e',
 'results/fresh_d14_small_component_context_retention_20261018.md':'10fb43eb7646aeb36b1f14c776d3defc76a6f63923d8d96c2f29bbbe01170f02',
 'results/review_finite57_current_S8_conditional_bridge_fresh_20261025.md':'2c97525c6005627fb78e627dbca31895271112fd02556f72a6de7034f9d23f48',
 'results/finite57_current_S8_conditional_bridge_root_acceptance_20261025.json':'7cfc187fe9a25e488e0952d063ab7771702a2732baeafaea6077c8ef628cff4f',
 'results/finite57_full_universal_affine_d_current_dependency_note_20261025.md':'cf91aba49029289e0034f87ef00e29be84b0daf3eb8c2b234e02b2b930e870a6',
}


def need(ok,message):
    if not ok:raise ValueError(message)


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=True)


def digest(value):
    return sha256(canonical(value).encode()).hexdigest()


def pin(path):
    need(path in PINS and callable(SOURCE_PINNER),'mandatory byte pinner and fixed source')
    SOURCE_PINNER(path,PINS[path])
    return dict(path=path,sha256=PINS[path])


def source(path):
    pin(path)
    need(callable(SOURCE_READER),'mandatory complete source reader')
    if path not in CACHE:CACHE[path]=SOURCE_READER(path,PINS[path])
    need(type(CACHE[path]) is dict,'complete parsed source document')
    return CACHE[path]


def source_acceptance():
    ref=SOURCE_ACCEPTANCE
    need(type(ref) is dict and set(ref)=={'path','sha256'} and
         type(ref['path']) is str and ref['path'].startswith('results/') and
         '..' not in ref['path'].split('/') and type(ref['sha256']) is str and
         re.fullmatch('[0-9a-f]{64}',ref['sha256']) and
         type(APPLICATION_SOURCE_SHA) is str and re.fullmatch('[0-9a-f]{64}',APPLICATION_SOURCE_SHA) and
         callable(SOURCE_PINNER) and callable(SOURCE_READER),
         'mandatory injected exact helper-source and supplemental root acceptance')
    SOURCE_PINNER(OWN,APPLICATION_SOURCE_SHA);SOURCE_PINNER(ref['path'],ref['sha256'])
    doc=SOURCE_READER(ref['path'],ref['sha256'])
    need(type(doc) is dict and doc['status']=='CURRENT_SOURCE_APPLICATION_PREMISES_ROOT_ACCEPTED' and
         doc['accepted_application_sources'][OWN]==APPLICATION_SOURCE_SHA and
         type(doc['accepted_supplemental_premise_ids']) is list and
         set(PREMISE_IDS)<=set(doc['accepted_supplemental_premise_ids']),
         'root explicitly accepts THIS helper and all additional complete full-D premises')
    return dict(path=ref['path'],sha256=ref['sha256'],application_source_sha256=APPLICATION_SOURCE_SHA,
                accepted_supplemental_premise_ids=list(PREMISE_IDS))


def indices(values,n,size=None):
    need(type(values) is list and len(values)==len(set(values)) and
         all(type(v) is int and 0<=v<n for v in values) and
         (size is None or len(values)==size),'distinct bounded physical indices')


def full_key(key):
    need(type(key) is dict and set(key)==FIELDS,'complete physical nine-field key')
    need(type(key['affine']) is bool and type(key['high_index']) is int and
         0<=key['high_index']<2614 and type(key['core_key']) is str and
         0<len(key['core_key'])<=256 and type(key['target_order']) is int,
         'strict full-key scalar fields')
    n=key['target_order'];need(3<=n<=59,'bounded target order')
    indices(key['base_nodes'],64);indices(key['selected_marks'],64,2)
    indices(key['induced_map'],n,len(key['base_nodes']));indices(key['marked_image'],n,2)
    need(1<=len(key['base_nodes'])<=17 and set(key['selected_marks'])<=set(key['base_nodes']) and
         key['marked_image']==sorted(key['induced_map'][key['base_nodes'].index(v)]
                                     for v in key['selected_marks']),
         'complete map and exact two marked images')
    return key


def require(units,worker,batch,evidence):
    need((worker,batch) in units,'missing ACTUAL full-universal-D dependency '+worker+' '+batch)
    entry,doc=units[worker,batch]
    need(type(entry) is dict and type(doc) is dict and entry['worker']==worker and
         entry['batch']==batch and doc['batch']==batch and
         entry['inner_status']==doc['status'],'same F.actual-validated entry/inner dependency')
    if worker in PINS:
        pin(worker);need(doc['source_sha256']==PINS[worker],'exact dependency worker source')
    for name in ('sha256','stdout_sha256'):
        need(type(entry[name]) is str and re.fullmatch('[0-9a-f]{64}',entry[name]),
             'actual ledger and stdout identities retained')
    identity=(worker,batch)
    if identity in evidence:need(evidence[identity]==entry,'conflicting actual dependency entries')
    evidence[identity]=entry
    return doc


def simple_gram(n,edges):
    g=[['2' if i==j else '0' for j in range(n)] for i in range(n)]
    for a,b in edges:
        need(a!=b and g[a][b]=='0','simple fixed source construction has unique edges')
        g[a][b]=g[b][a]='-1'
    return g


def affine_d(n):
    return [(i,i+1) for i in range(n-5)]+[(0,n-4),(0,n-3),(n-5,n-2),(n-5,n-1)]


def coefficient(value):
    need(type(value) in (str,int),'literal simple source coefficient')
    return {'-(1)':'-1','-(2)':'-2'}.get(str(value),str(value))


def component_map(g,c):
    """Deterministic A/D path orientation, not a native-map search."""
    nodes=c['nodes'];indices(nodes,len(g));n=len(nodes)
    match=re.fullmatch('([AD])([1-9][0-9]*)',c['type'])
    need(match is not None and int(match[2])==n,'saved finite A/D type/order')
    kind=match[1];adj={v:[w for w in nodes if w!=v and coefficient(g[v][w])=='-1'] for v in nodes}
    need(all(coefficient(g[v][w])==('2' if v==w else '-1' if w in adj[v] else '0')
             for v in nodes for w in nodes),'complete literal finite component Gram')
    if kind=='A':
        ends=[v for v in nodes if len(adj[v])<=1]
        need((n==1 and len(ends)==1) or (n>1 and len(ends)==2 and
             all(len(adj[v]) in (1,2) for v in nodes)),'saved A path degree pattern')
        order=[min(ends)]
        while len(order)<n:
            nxt=[v for v in adj[order[-1]] if v not in order]
            need(len(nxt)==1,'single complete saved A path');order.append(nxt[0])
    else:
        branch=[v for v in nodes if len(adj[v])==3]
        need(n>=5 and len(branch)==1 and all(len(adj[v]) in (1,2,3) for v in nodes),
             'unique saved D fork')
        b=branch[0];short=sorted(v for v in adj[b] if len(adj[v])==1)
        ends=[v for v in nodes if len(adj[v])==1 and v not in short]
        need(len(short)==2 and len(ends)==1,'two short D leaves and unique long end')
        order=[ends[0]]
        while order[-1]!=b:
            nxt=[v for v in adj[order[-1]] if v not in order]
            need(len(nxt)==1,'single complete D long arm');order.append(nxt[0])
        order+=short
    need(len(order)==n and set(order)==set(nodes),'complete deterministic finite map')
    return (kind,n),dict(zip(order,range(n)))


def current_role(g,white,corner,W,whole_factor_premise=None):
    marks=corner['marks'];indices(marks,len(g),2)
    indices(white,len(g));need(not set(white)&set(marks),'actual black marks are absent from white')
    cs=corner['components'];need(type(cs) is list and cs,'saved finite component partition')
    nodes=[v for c in cs for v in c['nodes']]
    need(len(nodes)==len(set(nodes)) and set(nodes)==set(white)|set(marks),
         'every physical actual corner node belongs to exactly one saved component')
    forest=[];maps=[];locations=[]
    for i,c in enumerate(cs):
        typ,mp=component_map(g,c);forest.append(typ);maps.append(mp)
        need(c['marks']==[v for v in marks if v in mp] or
             set(c['marks'])==set(marks)&set(mp),'complete saved physical marked membership')
    for a,c in enumerate(cs):
        for b,d in enumerate(cs):
            if a!=b:need(all(coefficient(g[v][w])=='0' for v in c['nodes'] for w in d['nodes']),
                         'literal whole components are pairwise orthogonal')
    for m in marks:
        hits=[(i,mp[m]) for i,mp in enumerate(maps) if m in mp]
        need(len(hits)==1,'unique actual mark component');locations.append(hits[0])
    def flip(i,v):
        kind,n=forest[i]
        return n-1-v if kind=='A' else 2*n-3-v if v>=n-2 else v
    (ca,a),(cb,b)=locations
    if ca==cb:
        orbit={(a,b),(flip(ca,a),flip(cb,b))}
        alternatives=[[(ca,x),(cb,y)] for x,y in sorted(orbit)]
    else:
        alternatives=[[(ca,x),(cb,y)] for x in sorted({a,flip(ca,a)})
                      for y in sorted({b,flip(cb,b)})]
    values={W.finite_corner_weight(forest,*pair) for pair in alternatives}
    need(len(values)==1 and min(values)>=0,'CURRENT lower invariant under complete fixed A/D role orbit')
    value=next(iter(values))
    return dict(actual_marks=marks,forest=[list(t) for t in forest],
        canonical_marks=[list(v) for v in locations],
        all_canonical_mark_images=[[[i,v] for i,v in pair] for pair in alternatives],current_lower=str(value),
        component_membership=[dict(nodes=c['nodes'],actual_marks=c['marks']) for c in cs],
        whole_factor_premise=whole_factor_premise or 'complete source zero-energy or energy-one seal; all other actual spectator components are orthogonal to the entire source core',
        evaluation_origin='THIS_ACTUAL_F_APPLICATION_NOT_A_HISTORICAL_RECEIPT')


def saved_roles(W):
    repairs=source(BOUNDARY)['affine_D13']['rows']
    expected={(side,a,b,leaf) for side in ('left','right')
              for a,b in ([(0,v) for v in (9,12,13)] if side=='left' else [(v,9) for v in (0,10,11)])
              for leaf in ((12,13) if side=='left' else (10,11))}
    got=[];out=[];single=0;double=0
    for i,p in enumerate(repairs):
        a,b=p['old_marks'];side=p['side'];leaf=p['attachment'];got.append((side,a,b,leaf))
        edges=affine_d(14)+([(14,1),(15,leaf)] if side=='left' else [(14,leaf),(15,8)])
        g=simple_gram(16,edges);white=sorted(set(range(14))-{a,b})
        slots=p['corners'];unit_only=p['old_marks']==[0,9]
        need(type(slots) is list and len(slots)==2 and type(slots[0]) is dict and
             ((unit_only and slots[1] is None) or
              (not unit_only and type(slots[1]) is dict and slots[1]['marks']==[14,15])),
             'exact boundary branch schema: four unit-only null second slots or eight distinct repairs')
        single+=int(unit_only);double+=int(not unit_only)
        roles=[current_role(g,white,c,W) for c in slots if c is not None]
        need(len(roles) in (1,2) and len({tuple(sorted(c['actual_marks'])) for c in roles})==len(roles),
             'distinct exact boundary receiver corners')
        need(Q(roles[0]['current_lower'])==(1 if unit_only else Q(1,2)),
             'unit-only D12 center pair or changed D12 leaf half is the first actual receiver')
        total=sum((Q(c['current_lower']) for c in roles),Q())
        need(total>=1,'CURRENT twelve-boundary source-role sum')
        out.append(dict(source_selector=[BOUNDARY,'affine_D13','rows',i],
            complete_key=[p['old_marks'],side,leaf],physical_context_sha256=digest([g,white]),
            source_corner_slot_status=['CORNER','STRUCTURAL_UNUSED_NONE' if unit_only else 'CORNER'],
            current_roles=roles,distinct_current_total=str(total)))
    need(len(got)==len(set(got))==12 and set(got)==expected and single==4 and double==8,
         'complete exact twelve boundary keys and four/eight receiver branch partition')
    rows=source(LEAF)['two_A1_spectators']['rows'];got=[];two=[]
    for i,p in enumerate(rows):
        a,b=p['marks'];c,d=p['leaf_attachments'];got.append((a,b,c,d))
        need(p['white']==sorted(set(range(13))-{a,b})+[15,16] and p['signature']==[15,1,1],
             'complete two-A1 source white and declared inertia')
        expected_g=simple_gram(17,affine_d(13)+[(13,c),(14,d),(13,15),(14,16)])
        need([[coefficient(v) for v in r] for r in p['gram']]==expected_g,
             'complete literal two-A1 physical source construction')
        need(len(p['corners'])==3,'complete three nonselected source corners')
        roles=[current_role(p['gram'],p['white'],p['corners'][j],W) for j in (0,2)]
        need(len({tuple(sorted(r['actual_marks'])) for r in roles})==2,
             'two ORIGINAL distinct side receivers')
        total=sum((Q(r['current_lower']) for r in roles),Q())
        need(total>=1,'CURRENT two-A1 side receivers pay without the middle A15 corner')
        two.append(dict(source_selector=[LEAF,'two_A1_spectators','rows',i],
            complete_key=[p['marks'],p['leaf_attachments']],physical_context_sha256=digest([p['gram'],p['white']]),
            current_roles=roles,distinct_current_total=str(total),middle_corner_used=False))
    expected={(a,b,c,d) for a in (9,10) for b in (11,12) for c in (9,10) for d in (11,12)}
    need(len(got)==len(set(got))==16 and set(got)==expected,'complete exact16 two-A1 source keys')
    return dict(boundary12=out,two_A1_16=two)


def order13_source_roles(W):
    doc=source(LEAF);nonextreme=doc['nonextreme_zero_cases'];got=[];out=[]
    for i,p in enumerate(nonextreme):
        a,b=p['marks'];side=p['nonextreme_side'];extreme=p['extreme_profile'];profile=p['nonextreme_profile']
        got.append((a,b,side,extreme,profile))
        need(p['white']==sorted(set(range(13))-{a,b}),'complete original nonextreme white')
        role=current_role(p['gram'],p['white'],p['chosen'],W)
        unit=p['chosen']['full_A1_D12_upgrade']
        need(type(unit) is bool,'exact source coupon-vs-D13 branch')
        if unit:
            need(Q(role['current_lower'])==1 and extreme=='isolated',
                 'CURRENT whole-A1/D12-shortleaf unit for exact isolated extreme branch')
        else:
            need(Q(role['current_lower'])==1 and extreme=='mate3' and i in CASE_INDICES,
                 'exact quiet D13 ordinary unit has separately joined negative-guard repair')
            role['ordinary_value_before_negative_guard']=role['current_lower']
            role['current_lower']='1/2'
        out.append(dict(source_selector=[LEAF,'nonextreme_zero_cases',i],
            complete_key=[p['marks'],side,extreme,profile],
            physical_context_sha256=digest([p['gram'],p['white']]),current_role=role,
            conditional_negative_D13_guard=not unit,
            required_guard_interface=None if unit else dict(quiet_source_index=i,
                ordered_witness_domain=[list(k) for k in product(('A','D'),repeat=3)])))
    expected={(a,b,s,e,p) for a in (9,10) for b in (11,12) for s in ('left','right')
              for e in ('isolated','mate3') for p in ('isolated4','mate4','center3')}
    need(len(got)==len(set(got))==48 and set(got)==expected,'complete48 original nonextreme source keys')
    one=[];got=[]
    for i,p in enumerate(doc['arbitrary_one_sided_cases']):
        a,b=p['marks'];c,d=p['leaf_attachments'];contact=p['only_spectator_contacting_black']
        got.append((a,b,c,d,contact))
        g=simple_gram(15,affine_d(13)+[(13,c),(14,d)])
        white=sorted(set(range(13))-{a,b});chosen=p['chosen']
        need(contact in (13,14) and contact not in chosen['marks'],
             'whole paying corner excludes the spectator-contacting opposite black')
        role=current_role(g,white,chosen,W,whole_factor_premise=
            'All actual spectators are orthogonal to the whole affine factor and the noncontacting new black. The PAYING corner excludes the spectator-contacting opposite black, so precisely its displayed marked factors are whole. The excluded black may have a nonzero spectator column; orthogonality to the entire four-black core is not assumed.')
        unit=chosen['full_A1_D12_upgrade']
        need(type(unit) is bool and Q(role['current_lower'])==1,
             'CURRENT one-sided whole coupon or ordinary D13 receiver')
        if not unit:
            role['ordinary_value_before_negative_guard']=role['current_lower']
            role['current_lower']='1/2'
        one.append(dict(source_selector=[LEAF,'arbitrary_one_sided_cases',i],
            complete_key=[p['marks'],p['leaf_attachments'],contact],
            physical_context_sha256=digest([g,white]),current_role=role,
            conditional_negative_D13_guard=not unit,
            required_guard_interface=None if unit else 'complete actual ordered one-sided-eight receivers at ambient58 energy-one classification'))
    expected={(a,b,c,d,k) for a in (9,10) for b in (11,12) for c in (9,10)
              for d in (11,12) for k in (13,14)}
    need(len(got)==len(set(got))==32 and set(got)==expected,'complete32 arbitrary one-sided source keys')
    return dict(nonextreme48=out,arbitrary_one_sided32=one)


def current_sum(row):
    corners=row['current_receivers'];marks=[tuple(sorted(c['current']['marks'])) for c in corners]
    need(corners and len(marks)==len(set(marks)),'distinct ORIGINAL current receiver corners')
    total=sum((Q(c['current']['weight']) for c in corners),Q())
    need(total==Q(row['current_total']) and total>=1,'actual distinct current receiver sum')


def lower_priorities(units,evidence):
    d=require(units,UNIVERSAL,'order12-d',evidence)['result']
    pairs=[(a,b) for a in (8,9) for b in (10,11)]
    expected1={canonical([list(p),s,b]) for p in pairs for s in ('left','right')
               for b in ('isolated4','mate3_1','D5','A1_A3')}
    expected2={canonical([list(p),a,b,s]) for p in pairs for a in (8,9)
               for b in (10,11) for s in (12,13)}
    for field,expected in (('current32',expected1),('separate32',expected2)):
        rows=d[field];keys=[canonical(r['application_key']) for r in rows]
        need(len(keys)==len(set(keys)) and set(keys)==expected,'exact actual lower-order '+field)
        for row in rows:current_sum(row)
    for a,b in pairs:
        apps=[p for p in d['theorem_applications'] if p['application_key']==['D',12,a,b] and
              p['theorem']=='finite57.current.affine_order12_D_opposite_leaves.whole_face.v1']
        need(len(apps)==1 and apps[0]['status']=='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION',
             'exact accepted order12-D theorem domain')
    local=require(units,UNIVERSAL,'five-local',evidence)
    theorem='finite57.current.D13_center_leaf.whole_face.v1'
    for a,b in ((0,11),(0,12),(8,9),(8,10)):
        apps=[p for p in local['result']['theorem_applications'] if p['theorem']==theorem and
              p['application_key']==['D',13,a,b]]
        need(len(apps)==1 and apps[0]['status']=='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION',
             'exact current order13 center/leaf priority application')
    # These are the actual finite interfaces named by the frozen local producer.
    for batch in ('universal-bc13','universal-d13','universal-d14-near','universal-d15-basic'):
        require(units,FINITE,batch,evidence)
    for batch in ('order11-bc','order12-near'):require(units,G4,batch,evidence)
    finite=require(units,FINITE,'universal-d13',evidence)['result']['outputs']
    originals=source(LOCAL_D13)['finite_corner_representatives']
    expected={canonical([p['side'],p['attached_leaf'],int(p.get('spectator_A1',False))]):p for p in originals}
    need(len(finite)==len(originals)==len(expected)==6 and
         {canonical(p['index']) for p in finite}==set(expected),
         'complete exact source-scoped six center/leaf finite keys')
    for row in finite:
        p=expected[canonical(row['index'])]
        need(row['source']==LOCAL_D13.removeprefix('results/').removesuffix('.json') and
             row['corner']['marks']==p['marks'] and Q(row['required'])==Q(p['existing_payment']) and
             Q(row['corner']['weight'])>=Q(row['required']),
             'same source/marks and actual current center/leaf finite lower')
    return dict(order12_current32=[p['application_key'] for p in d['current32']],
        order12_separate32=[p['application_key'] for p in d['separate32']],
        order13_center_leaf_keys=[['D',13,a,b] for a,b in ((0,11),(0,12),(8,9),(8,10))],
        no_higher_rank_back_edge=True)


def original_equalities(units,evidence):
    rows=[];domain=None
    for shard in range(4):
        part=require(units,DIRECT,f'original-{shard}',evidence)['result'];keys=part['complete_domain_keys']
        need(len(keys)==len({canonical(k) for k in keys})==216,'exact original216 source/case domain')
        if domain is None:domain=keys
        need(keys==domain and part['shard']==shard,'same complete original domain and actual shard')
        current=part['rows'];got=[]
        for row in current:
            cs=row['corners'];marks=[tuple(sorted(c['marks'])) for c in cs]
            need(cs and len(marks)==len(set(marks)) and row['direct_current_receiver_condition'] is True and
                 row['legacy_pointwise_bounds_used'] is False and
                 sum((Q(c['weight']) for c in cs),Q())==Q(row['total']) and Q(row['total'])>=1 and
                 len(row['final_collector_joins'])==len(cs),'actual complete CURRENT original face receipts')
            need(all(not c.get('guard_possible',False) or Q(c['weight'])<=Q(1,2) for c in cs),
                 'possible final D13 guard handled conservatively')
            got.append([row['source'],canonical(row['index'])]);rows.append(row)
        need(got==domain[shard::4],'exact original receipt rows, no skipped or repeated keys')
    out=[]
    for path,selector in ((D14,('extreme_equality','rows')),(D15,('equality','rows'))):
        data=source(path)
        for part in selector:data=data[part]
        need(len(data)==72,'complete accepted nine-pair equality product')
        keys=[]
        for ordinal,p in enumerate(data):
            stem=path.removeprefix('results/').removesuffix('.json')
            matched=[r for r in rows if r['source']==stem and r['index']==ordinal]
            need(len(matched)==1,'each exact whole-source equality has an actual CURRENT receipt')
            row=matched[0]
            if path==D14:
                left,right=p['leaf_attachments'];attach=p['A3_attached_to']
                g=simple_gram(19,[(i,i+1) for i in range(9)]+[(0,10),(0,11),(9,12),(9,13),
                    (14,left),(15,right),(16,attach),(16,17),(17,18)])
                # Frozen DIRECT.gram hashes these raw literal strings.
                # Preserve that historical representation identity exactly.
                g=[['-(1)' if v=='-1' else v for v in row] for row in g]
            else:g=p['gram']
            need(row['context_sha256']==digest([g,p['white']]) and
                 [c['marks'] for c in row['corners']]==[c['marks'] for c in p['corners']],
                 'same complete source Gram/white/ordered actual corners as CURRENT receipt')
            key=[p['marks'],p['leaf_attachments'],p['A3_attached_to'] if path==D14 else p['spectator_attached_to']]
            keys.append(canonical(key))
            out.append(dict(source_selector=[path,*selector,ordinal],complete_key=key,
                actual_source_case=[stem,ordinal],physical_context_sha256=row['context_sha256'],
                current_total=row['total'],ordered_actual_marks=[c['marks'] for c in row['corners']],
                current_receiver_rows_sha256=digest(row),negative_D13_guard_conservative=True))
        n=14 if path==D14 else 15
        expected={canonical([[a,b],[c,d],side]) for a in (0,n-4,n-3)
                  for b in (n-5,n-2,n-1) for c in (n-4,n-3)
                  for d in (n-2,n-1) for side in (n,n+1)}
        need(len(keys)==len(set(keys))==72 and set(keys)==expected,
             'each exact nine-pair/attachment/spectator-side equality key occurs once')
    return out


def guard_receivers(units,evidence):
    model=source(MODEL);one=require(units,ONE,'one-sided-eight',evidence)['receivers']
    expected={(a,b,c) for a in (9,10) for b,c in product(('A','D'),repeat=2)}
    need(len(one)==8 and {tuple(r['key']) for r in one}==expected,'all eight ordered one-sided keys')
    for row in one:
        need(Q(row['quiet_current']['weight'])==Q(row['opposite_current']['weight'])==Q(1,2) and
             Q(row['distinct_current_total'])==1 and set(row['quiet_marks'])!=set(row['opposite_marks']) and
             all(sum(rank for _,rank in row[k]['forest'])==58 and row[k]['complete'] is True
                 for k in ('quiet_current','opposite_current')),'same actual rank58 distinct one-sided halves')
        ordinal=row['source_selector'][1];p=model['one_sided_whole_models'][ordinal]
        need(row['source_selector']==['one_sided_whole_models',ordinal] and
             p['active_attachment']==row['key'][0] and p['quiet_corner']==row['quiet_marks'] and
             p['new_paying_corner']==row['opposite_marks'] and
             row['physical_model_sha256']==digest([p['gram'],p['white']]),
             'one-sided source model, physical identity and exact paying corners')
    legacy=(UNIVERSAL,'final-quiet');micro={(MICRO,f'quiet-case-{i:02d}') for i in CASE_INDICES}
    actual_micro={k for k in units if k[0]==MICRO}
    need(actual_micro<=micro and not (legacy in units and actual_micro),
         'no unknown quiet case or mixed complete/partial quiet mode')
    quiet=[];consumed=None
    if legacy in units:quiet=require(units,*legacy,evidence)['result']['theorem_applications']
    else:
        need(actual_micro==micro,'all24 actual quiet cases required')
        for i in CASE_INDICES:
            doc=require(units,MICRO,f'quiet-case-{i:02d}',evidence);part=doc['result']
            need(doc['base_universal_source_sha256']==PINS[UNIVERSAL] and
                 type(part['case_index']) is int and part['case_index']==i and
                 part['original24_case_domain']==list(CASE_INDICES) and
                 part['case_complete_at_declared_premises'] is True and
                 part['quiet24_complete_at_declared_premises'] is False,'THIS exact actual quiet case')
            if consumed is None:consumed=doc['consumed_inputs']
            need(doc['consumed_inputs']==consumed,'same quiet source/rule inputs for all24 actual cases')
            apps=part['theorem_applications'];keys=[tuple(p['application_key']) for p in apps]
            need(len(keys)==len(set(keys))==8 and set(keys)=={(i,*k) for k in product(('A','D'),repeat=3)},
                 'all eight ordered witness-copy keys of THIS physical quiet case')
            quiet.extend(apps)
    expected={(i,*k) for i in CASE_INDICES for k in product(('A','D'),repeat=3)}
    keys=[tuple(p['application_key']) for p in quiet]
    need(len(keys)==len(set(keys))==192 and set(keys)==expected,'complete192 exact quiet contexts')
    priority=require(units,QUIET,'quiet-b13',evidence)['exact_quiet_consumer_joins']
    providers={canonical(p['consumer_key']):p for p in priority};used=set()
    need(len(providers)==len(priority)==64,'exact64 distinct quiet-B13 consumer bindings')
    for app in quiet:
        need(app['same_three_witness_copies_at_both_corners'] is True and
             app['intrinsic_zero_credited_as_unit'] is False,'same actual whole witness copies; no zero promotion')
        if app['distinct_current_total'] is not None:
            need(Q(app['distinct_current_total'])>=1,'actual distinct quiet receiver total')
        else:
            token=canonical(app['application_key']);need(token in providers,'missing exact quiet-B13 priority')
            p=providers[token]
            need(p['intrinsic_current_weight']=='0' and
                 app['original_source_context_sha256']==p['physical_context_sha256'] and
                 app['original_source_selector'][1]=='nonextreme_zero_cases' and
                 app['original_source_selector'][2]==p['provider_key']['quiet_source_index'] and
                 app['paying_receiver']['actual_corner']==p['original_paying_marks'] and
                 app['paying_receiver']['current']['marked_images']==p['provider_key']['marked_images'] and
                 app['paying_receiver']['current']['weight']=='0','same original context/index/marks/complete B13 images')
            used.add(token)
    need(used==set(providers),'all and only exact quiet-B13 consumers are joined')
    return dict(one_sided_keys=sorted([list(k) for k in expected_one()]),
        quiet_keys=[p['application_key'] for p in quiet],quiet_B13_keys=[p['consumer_key'] for p in priority],
        receiver_tables_sha256=digest([one,quiet,priority]))


def expected_one():
    return {(a,b,c) for a in (9,10) for b,c in product(('A','D'),repeat=2)}


def support(units,W,spine):
    acceptance=source_acceptance()
    need(type(spine) is dict and spine.get('theorem')=='finite57.current.S8.whole_face.v1' and
         spine.get('status')=='RUNTIME_APPLICATION_INTERFACES_AT_ACCEPTED_COMPLETE_SOURCE_PREMISES' and
         spine.get('dag')==['S5','F6','S6','F7','S7','F8','S8'] and
         all(spine.get(k) is False for k in ('protected74_used_as_premise',
             'adaptive9_used_as_premise','largeD_half_used_as_premise','order12_D_used_as_premise')),
         'explicit already checked independent acyclic S8 runtime interfaces')
    g4=('order10-cycle','cycle11','order10-b','order10-d','order11-bc','order11-bcd',
        'cycle12','cycle14','cycle16','cycle18-tree')
    expected={(G4,b) for b in g4}|{(ZERO,'d13-center-zero')}|{
        (S8,f'deep-{i:02d}') for i in range(16)}|{(S8,'tree-immediate')}
    need({tuple(p) for p in spine['actual_units']}==expected and
         len(spine['actual_units'])==len(expected),'complete exact independent S8 interface set')
    fingerprint=digest([sorted((list(k),e) for k,(e,_) in units.items()),spine,acceptance])
    cachekey=('support',id(W),fingerprint)
    if cachekey in CACHE:return CACHE[cachekey]
    evidence={}
    for worker,batch in sorted(expected):require(units,worker,batch,evidence)
    require(units,S8,'association-controls',evidence)
    refs=[pin(path) for path in PINS]
    roles=saved_roles(W);roles.update(order13_source_roles(W));lower=lower_priorities(units,evidence)
    equalities=original_equalities(units,evidence);guards=guard_receivers(units,evidence)
    cid=PREFIX+'.contract.v1'
    contract=dict(status='COMPLETE_CURRENT_SOURCE_APPLICATION_SUPPORT_AT_EXPLICIT_PREMISES',
        contract_id=cid,source_inputs=refs,actual_receipts=[evidence[k] for k in sorted(evidence)],
        supplemental_source_acceptance=acceptance,
        independent_S8=spine,lower_priority_bindings=lower,original_equality_receivers=equalities,
        current_literal_role_evaluations=roles,final_guard_bindings=guards,
        assumptions=[
            dict(id='actual-star',kind='accepted_geometric_statement',runtime_regenerated=False,
                 statement='The complete original quadrilateral star, strict white/edge positivity, PSD actual corners, Lorentz index and actual opposite-pair interpretation hold.'),
            dict(id='native-classification',kind='accepted_complete_classification',runtime_regenerated=False,
                 statement='The native alphabet, finite/affine one-root profiles and full positive-definite spectator energy classifications hold at their all-root scopes.'),
            dict(id='finite57.fullD.complete_profiles',kind='SUPPLEMENTAL_ACCEPTED_COMPLETE_SOURCE_PREMISE',runtime_regenerated=False,
                 statement='The complete all-nine end profiles, full-metric inverse/Schur/height statements, native gaps, unrestricted opposite products and exact computed numerical identities hold at the pinned 0ecf report/a02 root scope; these are additional complete premises, not an enlargement of an ordinary manifest record.',
                 accepted_by=acceptance),
            dict(id='finite57.fullD.current_boundary_roles',kind='SUPPLEMENTAL_ACCEPTED_COMPLETE_SOURCE_PREMISE',runtime_regenerated=False,
                 statement='The energy3/4 and1/2 whole A3/A1 equality seals, exact original144 current output interfaces, complete twelve boundary receivers, order13 energy-one/zero-energy exhaustiveness and ambient58 final rank15 guard classification hold at the pinned full-source/current-bridge scopes. THIS F call evaluates saved typed current roles and binds external actual receipts; it does not regenerate the geometric/numerical precursor proofs.',
                 accepted_by=acceptance),
            dict(id='finite57.fullD.whole_face_scope',kind='SUPPLEMENTAL_ACCEPTED_COMPLETE_SOURCE_PREMISE',runtime_regenerated=False,
                 statement='The complete order13 opposite-leaf and order14/15 nine-pair theorems apply to every actual whole factor with all supporting roots and all native labels. Lower immediate whole-face priorities and the independent S8 theorem have their stated actual same-face meanings. This full scope is separately accepted, not inferred from the narrower five-local ordinary premise.',
                 accepted_by=acceptance),
            dict(kind='same_actual_whole_affine_target',runtime_regenerated=False,
                 statement='THIS selected full key is realized by the whole varying actual affine factor; all other actual factors at the ideal corner are finite and orthogonal to it.')],
        proof_logic='Assume the SAME actual face has total current charge<1. Independent S8 forbids any adjacent white-interior path<=8. The complete pinned universal end/energy case split then excludes the realized factor or supplies the distinct CURRENT receivers/lower whole-face priority recorded here, a contradiction.',
        historical_finite_order_cutoff_used=False,historical_selected_half_restored=False,
        current_role_evaluation_scope='Executed by THIS actual enclosing F application; no standalone numerical receipt is asserted.',
        geometric_premises_proved_by_code=False,full_transitive_regeneration=False,
        global_finite57_verified=False)
    if cid in CONTRACTS:need(canonical(CONTRACTS[cid])==canonical(contract),'shared contract conflict')
    CONTRACTS[cid]=contract;CACHE[cachekey]=cid
    return cid


def apply_selected(key,W,units,typed_predicate):
    full_key(key)
    if not key['affine'] or key['target_kind']!='D' or key['target_order'] not in (13,14,15):return None
    n=key['target_order'];s=n-5;pair=key['marked_image']
    left={n-4,n-3} if n==13 else {0,n-4,n-3};right={n-2,n-1} if n==13 else {s,n-2,n-1}
    if not (len(set(pair)&left)==len(set(pair)&right)==1):return None
    need(type(typed_predicate) is dict and typed_predicate.get('same_full_key')==key and
         typed_predicate.get('whole_actual_affine_factor') is True,
         'THIS exact selected key and explicit whole-actual-target hypothesis')
    cid=support(units,W,typed_predicate['actual_runtime_S8'])
    graph=W.classical_graph('D',n,affine=True);distance=W._distances(graph,pair[0])[pair[1]]
    return dict(status='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION',
        theorem=f'finite57.current.full_affine_D_order{n}.whole_face.v1',contract_id=cid,
        selected_full_key=key,typed_hypotheses=dict(whole_actual_affine_factor=True,
            affine_kind='D',affine_order=n,affine_rank=n-1,canonical_marks=pair,
            complete_left_mark_domain=sorted(left),complete_right_mark_domain=sorted(right),
            marked_distance=distance,current_intrinsic_weight=str(W.affine_weight('D',n,*pair))),
        conclusion=dict(quantity='total_current_charge_of_same_actual_face',lower='1'),
        intrinsic_corner_unit=False,actual_short_path_asserted=False,
        geometric_premises_proved_here=False)


def application_contracts():
    """Caller MUST serialize these expanded contracts, not only their hashes."""
    return dict(CONTRACTS)
