#!/usr/bin/env python3
"""Pure physical-key applications of reviewed affine-D source arguments.

No repository import, file I/O, native-map search, or standalone execution.
The final reader supplies its PINNED original high/core objects, frozen W,
F.actual-validated units, and byte-verifying SOURCE_PINNER/SOURCE_READER.
Predicates concern physical profiles, never a high-index/shard allowlist.
Unsupported profiles return None; geometry is retained as an assumption.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

from fractions import Fraction as F
from hashlib import sha256
import json
import re

SOURCE_PINNER = None
SOURCE_READER = None
CACHE = {}
FIELDS = {'high_index','core_key','base_nodes','selected_marks','target_kind',
          'target_order','affine','induced_map','marked_image'}
HIGH = 'results/fresh_distance8_high_overlap_20261017.json'
LEAF = 'results/fresh_affine_d12_leaf_pair_transfer_20261017.json'
D14 = 'results/fresh_affine_d13_long_pair_transfer_20261017.json'
D15 = 'results/fresh_affine_d14_long_pair_transfer_20261017.json'
MODEL = 'results/fresh_d13_three_rank15_context_transfer_20261018.json'
TEX = 'manuscripts/finite_volume57.tex'
LITERAL = 'scripts/finite57_literal_weights_audit_20261024.py'
LOCAL = 'results/audit_finite57_local_routing_whole_components_fresh_20261024.md'
PROFILES = 'results/audit_finite57_nine_fork_profiles_and_twelve_d13_repairs_20261024.md'
PROFILE_ROOT = 'results/review_finite57_remaining_universal_profiles_root_20261025.md'
BRIDGE = 'results/finite57_current_universal_receiver_retention_bridge_20261025.md'
BRIDGE_ROOT = 'results/finite57_current_universal_receiver_bridge_root_acceptance_20261025.json'
ONE = 'scripts/check_finite57_final_guard_one_sided_current_receivers_20261025.py'
ONE_SHA = '0892b27a233d761b9568067a3630c73b63f25e551e20b4fd9b69ae6ba99425a3'
ONE_ROOT = 'results/finite57_final_guard_one_sided_source_clearance_root_20261025.json'
PINS = {
    HIGH:'f8a2ed3bdef6064febf7336bbe4ad75d71c7d35b5747d9f57b0ebf26b03d0cc9',
    LEAF:'8f24d38cca6526a25652835c04923f542338590509c1199bcba2707520f62eb5',
    D14:'79989c70ace2cc1591d88e4167de39dc243977b748a6d369a495c313945c8a2a',
    D15:'42fea3ab62fd837ed3652a289d889107f6933d323096c7961c971ab28cd023fd',
    MODEL:'2c7d615cd38123bef20cffad68de8ab30baa91826c147a6bb084892590e8ca86',
    TEX:'c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0',
    LITERAL:'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0',
    LOCAL:'12f865b1373e043d5d1b9bbebb13b3e93a5e9b3a1519e3073f50506d8c4a4785',
    PROFILES:'0ecf41083477787cfd20af8dc899138e1480db0aec7bf5ce6115838ef6631146',
    PROFILE_ROOT:'a02e9daa3e3f42a3d3b832bc561702b87e658ad95693644bf7469fe91c1b6201',
    BRIDGE:'309bbd85a766466bcb3bb4437d7a5d50ffe8c614c4892f40f4ccc06d842893eb',
    BRIDGE_ROOT:'6bc292fb08a50c15cfa1cff2ef9e43def44d4febd942904f0cdd1033eac4838e',
    ONE:ONE_SHA,
    ONE_ROOT:'dab08dc5abf52a07cc57dff538ebf7ed0891ce36f1bc4b961d847ca4f6f585da',
}


def need(ok,message):
    if not ok:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=True)


def pin(name,expected=None):
    expected = PINS[name] if expected is None else expected
    need(type(expected) is str and re.fullmatch('[0-9a-f]{64}',expected)
         and callable(SOURCE_PINNER), 'mandatory exact source/evidence byte pinner')
    SOURCE_PINNER(name,expected)
    return dict(path=name,sha256=expected)


def source(name):
    pin(name)
    need(callable(SOURCE_READER),'mandatory complete parsed-source reader')
    if name not in CACHE:
        CACHE[name] = SOURCE_READER(name,PINS[name])
    need(type(CACHE[name]) is dict,'complete named source document')
    return CACHE[name]


def literal(value):
    need(type(value) in (str,int) and len(str(value)) <= 128,'bounded literal Gram coefficient')
    return {'-(1)':'-1','-(2)':'-2'}.get(str(value),str(value))


def indices(values,n,size=None):
    need(type(values) is list and len(values) == len(set(values))
         and all(type(v) is int and 0 <= v < n for v in values)
         and (size is None or len(values) == size),'ordered distinct physical indices')


def physical(key,p,core):
    need(type(key) is dict and set(key) == FIELDS,'complete nine-field physical selected key')
    need(type(key['affine']) is bool and type(key['high_index']) is int
         and 0 <= key['high_index'] < 2614 and type(key['target_order']) is int
         and type(key['core_key']) is str,'selected key scalar types')
    if not key['affine'] or key['target_kind'] != 'D' or key['target_order'] not in (13,14,15):
        return None
    N = key['target_order'];s = N-5
    need(type(p) is dict and type(core) is dict and p['index'] == key['high_index']
         and p['core_key'] == key['core_key'] and p['marked_component'] == key['base_nodes']
         and p['preferred_corner'] == key['selected_marks'],
         'caller PINNED high/core object is THIS original selected occurrence')
    raw = core['gram'];n = len(raw)
    need(type(raw) is list and 4 <= n <= 64 and all(type(row) is list and len(row) == n for row in raw),
         'complete bounded original physical Gram')
    g = [[literal(v) for v in row] for row in raw]
    need(all(g[i][j] == g[j][i] for i in range(n) for j in range(n)), 'symmetric physical Gram')
    nodes,marks,mp = (key[name] for name in ('base_nodes','selected_marks','induced_map'))
    indices(nodes,n);indices(marks,n,2);indices(mp,N,len(nodes));indices(key['marked_image'],N,2)
    need(10 <= len(nodes) <= N and set(marks) <= set(nodes)
         and key['marked_image'] == sorted(mp[nodes.index(v)] for v in marks),
         'complete ordered native map/physical marked association')
    to = dict(zip(nodes,mp));back = {v:k for k,v in to.items()}
    edge = lambda i,j: abs(i-j) == 1 and max(i,j) <= s or {i,j} in (
        {0,N-4},{0,N-3},{s,N-2},{s,N-1})
    if any(g[i][j] != ('2' if i == j else '-1' if edge(to[i],to[j]) else '0')
           for i in nodes for j in nodes):
        return None
    left,right = {0,N-4,N-3},{s,N-2,N-1}
    if not (len(set(key['marked_image']) & left) == len(set(key['marked_image']) & right) == 1):
        return None
    a = next(m for m in marks if to[m] in left);b = next(m for m in marks if to[m] in right)
    oldleft = p['left'][1];oldright = p['right'][1];copy = p['right_copy_map']
    need(type(copy) in (list,dict),'saved old-to-new copy-map schema')
    def copied(i):
        value = copy[i] if type(copy) is list else copy[str(i)] if str(i) in copy else copy[i]
        need(type(value) is int and 0 <= value < n,'integer copied original root')
        return value
    pairs = [oldleft,[copied(i) for i in oldright]]
    for pair in pairs:indices(pair,n,2)
    black = set(pairs[0]) | set(pairs[1])
    need(len(black) == 4 and all(len(set(pair) & set(marks)) == 1 for pair in pairs),
         'two disjoint ORIGINAL opposite pairs with one selected mark each')
    opposite = {v:next(w for w in pair if w != v) for pair in pairs for v in pair}
    c,d = opposite[a],opposite[b]
    need(set(p['omitted_black']) == {c,d} and not ({c,d} & set(nodes)), 'same original omitted black roots')
    corners = p['all_corners'];expected = {tuple(sorted((i,j))) for i in pairs[0] for j in pairs[1]}
    need(type(corners) is list and len(corners) == 4
         and {tuple(sorted(q['preferred_corner'])) for q in corners} == expected
         and all(set(q['omitted_black']) == black-set(q['preferred_corner']) for q in corners),
         'complete four ORIGINAL adjacent corners, not arbitrary graph pairs')
    if any(g[i][j] != '0' for i,j in ((a,d),(b,c),(c,d))):return None
    white = sorted(set(range(n))-black)
    other = sorted(set(white)-set(nodes))
    if any(g[i][j] != '0' for i in other for j in nodes):return None
    to.update({c:N,d:N+1});back.update({N:c,N+1:d})
    return dict(N=N,s=s,g=g,n=n,p=p,nodes=nodes,marks=marks,a=a,b=b,c=c,d=d,
        to=to,back=back,white=white,other=other,missing=sorted(set(range(N))-set(mp)),
        corner_keys=expected,physical_context_sha256=sha256(canonical([raw,white]).encode()).hexdigest())


def original_corner(x,marks):
    need(tuple(sorted(marks)) in x['corner_keys'],'paying pair is an ORIGINAL actual adjacent corner')
    return marks


def coupon(x,W):
    if x['N'] != 13 or x['missing']:return None
    a,b = x['to'][x['a']],x['to'][x['b']]
    if a not in (9,10) or b not in (11,12):return None
    rows = source(LEAF)['nonextreme_zero_cases']
    need(type(rows) is list and len(rows) == 48 and len({(tuple(q['marks']),q['nonextreme_side'],
        q['extreme_profile'],q['nonextreme_profile']) for q in rows}) == 48,'complete original source48 profile keys')
    back = x['back'];g = x['g'];matches = []
    for i,q in enumerate(rows):
        if q['marks'] != [a,b] or q['chosen']['full_A1_D12_upgrade'] is not True:continue
        if all(g[back[j]][back[k]] == literal(q['gram'][j][k]) for j in range(15) for k in range(15)):
            matches.append((i,q))
    if not matches:return None
    need(len(matches) == 1,'unique complete source48 principal Gram identity')
    i,q = matches[0]
    if any(g[j][k] != '0' for j in x['other'] for k in (x['c'],x['d'])):return None
    need(q['white'] == sorted(set(range(13))-{a,b}) and q['chosen']['signature'] == [13,0,0]
         and W.finite_cross_weight(('A',1),0,('D',12),10) == 1
         and W.finite_cross_weight(('A',1),0,('D',12),11) == 1,
         'exact current whole-A1/whole-D12 short-leaf unit, not old finite<=17')
    paid = original_corner(x,[back[m] for m in q['chosen']['marks']])
    return dict(kind='CURRENT_COUPON',source=LEAF,selector=['nonextreme_zero_cases',i],
        source_context=dict(gram=q['gram'],white=q['white']),
        original_to_source=sorted(x['to'].items()),paying_corners=[dict(marks=paid,current_lower='1')],
        source_seal='Nonextreme8+4=R12 and full positive spectator metric force both ENTIRE disjoint columns zero; every spectator remains orthogonal to the full star.',
        proof_references=[pin(LOCAL),pin(BRIDGE),pin(BRIDGE_ROOT)],actual_receipts=[])


def strict_profile(x):
    N,s,g,to,back = (x[k] for k in ('N','s','g','to','back'))
    if N not in (14,15) or len(x['missing']) > 1:return None
    pair = [to[x['a']],to[x['b']]]
    if (pair[0] == 0) == (pair[1] == s):return None
    ci = 0 if pair[0] == 0 else 1;li = 1-ci
    old = [x['a'],x['b']];new = [x['c'],x['d']]
    if g[new[0]][new[0]] != '2' or g[new[1]][new[1]] != '2':return None
    center = 0 if ci == 0 else s;outer = {N-4,N-3} if ci == 0 else {N-2,N-1}
    entries = [1,2] if ci == 0 else [s-1,s-2]
    leafcenter = 0 if li == 0 else s
    leaf_known = g[new[li]][old[li]] == '0' and leafcenter in back and g[new[li]][back[leafcenter]] == '-1'
    leaf_known = leaf_known and all(g[new[li]][v] == ('-1' if to[v] == leafcenter else '0')
                                   for v in x['nodes'] if v != old[li])
    center_outer = g[new[ci]][old[ci]] == '-1' and all(v in back and g[new[ci]][back[v]] == '0' for v in entries)
    center_outer = center_outer and all(g[new[ci]][v] == '0' for v in x['nodes']
                                        if to[v] not in outer | {center})
    if not leaf_known and not center_outer:return None
    stem = D14 if N == 14 else D15;doc = source(stem)
    prefix = 'formulas' if N == 14 else 'profile_formulas';rows = doc[prefix]['end_rows']
    need(type(rows) is list and len(rows) == 96,'complete original96 end-identity domain')
    profile,side = ('center3','left' if li == 0 else 'right') if leaf_known else ('outer_general','left' if ci == 0 else 'right')
    found = [(i,q) for i,q in enumerate(rows) if q['marks'] == pair and q['side'] == side and q['profile'] == profile]
    need(len(found) == 1,'one exact original center/leaf profile identity')
    i,q = found[0];R = N-3
    need(q['R'] == R,'exact source resistance for THIS marked pair')
    bound = 10
    gap = str(F(R-bound,R))
    need(F(gap) > 0,'strict positive source determinant margin')
    return dict(kind='EXCLUDED_COMPLETION_UNDER_NO_S8',source=stem,selector=[prefix,'end_rows',i],
        profile_predicate='leaf_center3_t0' if leaf_known else 'center_outer_t1',
        observed_native_map=sorted(to.items()),missing_target_nodes=x['missing'],
        missing_outer_contacts_retained=sorted(set(x['missing']) & outer),
        full_contact_classification='Complete native one-root profiles; observed inner-entry zeros exclude inner3/inner4/inner_next3. Missing outer contacts remain arbitrary with u^2+v^2<4. A witnessed leaf center3 excludes every extra contact by that same complete classification.',
        determinant=dict(resistance=R,cost_sum_upper=bound,positive_lower=gap,
            bounds='leaf_center3<=4 and every center<=6' if leaf_known else 'outer-only center with opposite t=1 <=2 and every leaf<=8',
            full_spectator_correction='nonnegative in the ORIGINAL positive spectator metric'),
        paying_corners=[],actual_receipts=[],proof_references=[pin(PROFILES),pin(PROFILE_ROOT)])


def one_sided_receipts(units,attachment):
    batch = 'one-sided-eight';status = 'PASS_NAMED_FINAL_GUARD_EIGHT_ONE_SIDED_CURRENT_RECEIVERS'
    pin(ONE);pin(ONE_ROOT)
    need((ONE,batch) in units,'required ACTUAL eight one-sided receiver receipt is missing')
    entry,doc = units[ONE,batch]
    need(entry['worker'] == ONE and entry['batch'] == doc['batch'] == batch
         and entry['inner_status'] == doc['status'] == status and doc['source_sha256'] == ONE_SHA,
         'exact actual0892 worker/source/batch/status association')
    for name in (MODEL,TEX,LITERAL):
        need(doc['consumed_inputs'].get(name) == PINS[name],'actual receiver consumed THIS exact model/current rule')
    evidence = dict(worker=ONE,batch=batch,source_sha256=ONE_SHA)
    for path,sha in (('path','sha256'),('stdout_path','stdout_sha256')):
        pin(entry[path],entry[sha]);evidence[path] = entry[path];evidence[sha] = entry[sha]
    models = source(MODEL)['one_sided_whole_models']
    rows = doc['receivers'];keys = [[a,b,c] for a in (9,10) for b in ('A','D') for c in ('A','D')]
    need(type(rows) is list and [q['key'] for q in rows] == keys,'complete exact8 ordered one-sided model keys')
    selected = []
    for q in rows:
        a = q['key'][0];selector = ['one_sided_whole_models',a-9];model = models[a-9]
        need(q['source_selector'] == selector and q['quiet_marks'] == [14,9] and q['opposite_marks'] == [13,14]
             and q['inactive_whole_factors'] == [[kind,15] for kind in q['key'][1:]]
             and q['active_whole_D15_nodes'] == list(range(15,30))
             and q['physical_model_sha256'] == sha256(canonical([model['gram'],model['white']]).encode()).hexdigest()
             and q['quiet_current']['marks'] == model['quiet_corner']
             and q['opposite_current']['marks'] == model['new_paying_corner']
             and q['quiet_current']['sealed_witness_components'] == [list(range(15,30))]
             and q['opposite_current']['sealed_witness_components'] == [],
             'same physical source model/ordered whole witnesses/actual corner records')
        for current in (q['quiet_current'],q['opposite_current']):
            need(current['kind'] == 'finite' and current['complete'] is True
                 and current['ignored_unsealed_components'] == [] and F(current['weight']) == F(1,2)
                 and sum(t[1] for t in current['forest']) == 58,'actual complete rank58 CURRENT half receiver')
        need(F(q['distinct_current_total']) == 1,'actual distinct-corner current total is1')
        if a == attachment:selected.append(q)
    need(len(selected) == 4 and doc['global_finite57_verified'] is False
         and doc['current_terminal_bridge_complete'] is False,'scoped four model applications, no global promotion')
    return selected,evidence


def one_sided(x,W,units):
    if x['N'] != 13 or x['missing'] or not 4 <= len(x['other']) <= 13:return None
    g,to,back = x['g'],x['to'],x['back'];a,b = to[x['a']],to[x['b']]
    if a not in (9,10) or b not in (11,12):return None
    contacts = [[to[v] for v in x['nodes'] if g[new][v] != '0'] for new in (x['c'],x['d'])]
    if any(len(v) != 1 for v in contacts) or contacts[0][0] not in (9,10) or contacts[1][0] not in (11,12):return None
    if any(g[new][back[v[0]]] != '-1' for new,v in zip((x['c'],x['d']),contacts)):return None
    support = [[v for v in x['other'] if g[new][v] != '0'] for new in (x['c'],x['d'])]
    if sorted(map(bool,support)) != [False,True]:return None
    active = 0 if support[0] else 1
    if len(support[active]) != 1:return None
    root = support[active][0];new = (x['c'],x['d'])[active]
    if g[new][root] != '-1':return None
    # A single rooted simple path, walked deterministically; no embedding search.
    path = [];previous = None;current = root
    while current is not None:
        if current in path:return None
        path.append(current)
        neighbors = [v for v in x['other'] if v != current and g[current][v] != '0']
        if any(g[current][v] != '-1' for v in neighbors):return None
        onward = [v for v in neighbors if v != previous]
        if len(onward) > 1:return None
        previous,current = current,onward[0] if onward else None
    if set(path) != set(x['other']) or any(g[i][i] != '2' for i in range(x['n'])):return None
    # Reflect the spine when active is right; swap fork leaves to put marks at9,11.
    reflect = (lambda j:8-j if j <= 8 else j+2 if j <= 10 else j-2) if active else (lambda j:j)
    image = {v:reflect(j) for v,j in to.items() if j < 13}
    oldactive,oldquiet = (x['a'],x['b']) if not active else (x['b'],x['a'])
    if image[oldactive] == 10:image = {v:(19-j if j in (9,10) else j) for v,j in image.items()}
    if image[oldquiet] == 12:image = {v:(23-j if j in (11,12) else j) for v,j in image.items()}
    image[new] = 13;inactive_new = (x['d'],x['c'])[active];image[inactive_new] = 14
    image.update({v:15+i for i,v in enumerate(path)})
    attachment = image[back[contacts[active][0]]]
    if image[back[contacts[1-active][0]]] != 12:return None
    model_index = attachment-9
    models = source(MODEL)['one_sided_whole_models']
    need(type(models) is list and [p['active_attachment'] for p in models] == [9,10],'complete two-model domain')
    model = models[model_index]
    if any(g[i][j] != literal(model['gram'][image[i]][image[j]]) for i in range(x['n']) for j in range(x['n'])):return None
    need(model['quiet_corner'] == [14,9] and model['new_paying_corner'] == [13,14], 'exact source model corner association')
    inverse = {j:v for v,j in image.items()}
    quiet = original_corner(x,[inverse[j] for j in model['quiet_corner']])
    pay = original_corner(x,[inverse[j] for j in model['new_paying_corner']])
    rows = source(LEAF)['arbitrary_one_sided_cases']
    selected = [(i,q) for i,q in enumerate(rows) if q['marks'] == [a,b]
        and q['leaf_attachments'] == [v[0] for v in contacts]
        and q['only_spectator_contacting_black'] == 13+active]
    need(len(rows) == 32 and len(selected) == 1,'one exact original32 conditional support key')
    si,q = selected[0]
    need(set(q['chosen']['marks']) == {to[v] for v in quiet}
         and W.finite_same_weight(('D',13),0,11) == W.finite_same_weight(('D',13),0,12) == 1,
         'same source quiet D13 role; ordinary guard-FALSE value is1')
    actual,evidence = one_sided_receipts(units,attachment)
    return dict(kind='CURRENT_ONE_SIDED_GUARD_REPAIR',source=LEAF,selector=['arbitrary_one_sided_cases',si],
        model_source=MODEL,model_selector=['one_sided_whole_models',model_index],
        original_to_model=sorted(image.items()),observed_active_path=path,
        paying_corners=[dict(marks=quiet,guard_false_current_lower='1',guard_true_current_lower='1/2'),
                       dict(marks=pay,guard_true_current_lower='1/2')],
        actual_receiver_keys=[q['key'] for q in actual],actual_receipts=[evidence],
        source_seal='Under NO-S8 the energy-one both-support classification cannot retain a contacted component containing this path of >=4 nodes: A1/A3 are too small, and D4..D7 give a white-interior path4. Thus the full opposite support column vanishes. Under the ACTUAL rank58 guard exactly one whole D15 is active; two inactive whole A15/D15 factors remain orthogonal to the ENTIRE star at BOTH corners.',
        proof_references=[pin(LOCAL),pin(BRIDGE),pin(BRIDGE_ROOT),pin(ONE_ROOT)])


def apply_selected(key,highrow,core,W,units):
    """None outside supported profiles; caller retains exact pinned HIGH origin."""
    x = physical(key,highrow,core)
    if x is None:return None
    bound = strict_profile(x) or coupon(x,W) or one_sided(x,W,units)
    if bound is None:return None
    frozen = json.loads(canonical(key));excluded = bound['kind'] == 'EXCLUDED_COMPLETION_UNDER_NO_S8'
    return dict(status='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION',
        theorem_id='finite57.affine_D_source.'+bound['kind']+'.v1',
        application_kind=bound['kind'],application_key=dict(selected_full_key=frozen,
            source=bound['source'],source_selector=bound['selector']),
        original_source=pin(HIGH),original_core_key=key['core_key'],
        original_physical_context_sha256=x['physical_context_sha256'],
        rule_identity=dict(manuscript=pin(TEX),literal_source=pin(LITERAL)),
        exact_source_application=bound,
        assumptions=[
            dict(kind='actual_same_face_and_WHOLE_selected_affine_factor',same_full_key=frozen,
                statement='Caller supplies the exact original high/core objects from the pinned route domain. THIS full map represents the entire actual varying affine factor, with every omitted target root actual, rather than a displayed component.'),
            dict(kind='complete_native_profiles_and_source_computed_identities',runtime_regenerated=False,
                statement='Accepted actual-star/factor-fiber, strict actual edge positivity, complete native finite/affine classification, positive FULL spectator metric, source inverse/determinant/corner/image and equality-sealing propositions remain premises.'),
            dict(kind='NO_ADJACENT_WHITE_INTERIOR_PATH_AT_MOST_8_IN_THIS_SAME_FACE',runtime_verified=False,
                statement='Independent current S8 is a separate whole-face theorem. This helper binds the exhaustive source implication under NO-S8; no short path is asserted from a flag.'),
            dict(kind='actual_rank58_guard_and_whole_rank15_classification',required=bound['kind']=='CURRENT_ONE_SIDED_GUARD_REPAIR',
                runtime_geometry_verified=False,
                statement='The one-sided repair splits on the ACTUAL guard. Three genuine whole unmarked rank15 A/D witnesses and the accepted complete energy-one classification supply exactly one active D15 and two inactive factors. The displayed path is only a principal prefix; no roots are padded.')],
        provider_composition_premise=dict(status='EXACT_SOURCE_NO_S8_IMPLICATION_BOUND',same_full_key=frozen,
            implication='NO same-face adjacent white-interior path<=8 excludes this exact completion.' if excluded else
                'NO same-face adjacent white-interior path<=8 implies current DISTINCT-corner total>=1 for this exact physical source occurrence.',
            independent_S8_composition='An unpaid face cannot have an adjacent white-interior path<=8 by independent current S8; the bound source implication then contradicts unpaidness.'),
        conclusion=dict(quantity='total_current_charge_of_same_actual_face',lower='1',
            route='impossible exact completion under NO-S8' if excluded else 'current distinct actual-corner receipts under NO-S8'),
        intrinsic_corner_unit=False,theorem_execution=False,selected_provider_union_complete=False,
        current_terminal_bridge_complete=False,full_transitive_regeneration=False,
        geometric_formalization=False,global_finite57_verified=False)


if __name__ == '__main__':
    raise SystemExit('Pure source-application helper: root must load pinned definitions in the final reader')
