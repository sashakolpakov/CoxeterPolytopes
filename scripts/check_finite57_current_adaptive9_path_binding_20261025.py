#!/usr/bin/env python3
"""SOURCE ONLY: exact adaptive9 prefixes at the CURRENT ambient58 rule.

The all-native profile/whole-bound proof is an explicit accepted premise.
Independent current S8 remains a structured OPEN prerequisite, never an
intrinsic unit. No legacy mathematical worker is imported by this adapter.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import types

ROOT = Path(__file__).resolve().parents[1]
ROUTE = 'scripts/index_finite57_current_terminal_routes_syntax_safe_20261025.py'
ROUTE_SHA = '43a0c4e986cbdb92d112771bf0edc1d128ce2d0feced30fffd9d9825dda7c29d'
EXTRA = 'results/finite57_current_adaptive9_extra_inputs_20261025.json'
EXTRA_SHA = 'a79e1b3d662dccbbf46ef16d977a4123f909be8b8b3aa9bde604ddf1c54f2a29'
ADAPT = 'fresh_distance8_adaptive_leaf_overlap_20261017'
INTERFACE = 'fresh_d12_low_protected_adaptive_interfaces_20261018'
MAX_OUTPUT = 500 * 1024
CASES = (
 (22,(1,2),25,(1,2),((5,9),(6,8),(7,7),(8,6))),
 (25,(0,1),25,(1,2),((3,9),)),
 (25,(0,1),25,(1,2),((5,9),(6,8),(7,7),(8,6),(9,5))),
 (25,(0,1),25,(1,4),((9,3),)),
 (25,(1,2),22,(1,2),((6,8),(7,7),(8,6),(9,5))),
 (25,(1,2),25,(0,1),((9,3),)),
 (25,(1,2),25,(0,1),((5,9),(6,8),(7,7),(8,6),(9,5))),
 (25,(1,2),25,(1,2),((4,9),(5,8),(6,7),(7,6),(8,5),(9,4))),
 (25,(1,4),25,(0,1),((3,9),)),
)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative, cap=4*1024**2):
    p = (ROOT/relative).resolve()
    need(p.is_relative_to(ROOT) and p.is_file(), 'repository source path')
    with p.open('rb') as stream:
        b = stream.read(cap+1)
    need(0 < len(b) <= cap, 'bounded source read')
    return b


def rational(x):
    need(type(x) in (str, int) and len(str(x)) <= 128, 'bounded rational')
    return F(x)


def profile_key(p):
    return (tuple(p['support']), p['scale_square'])


def source_key(p):
    return (p['left'][0], tuple(p['left'][1]), p['right'][0],
            tuple(p['right'][1]), tuple(map(tuple, p['match'])))


def prepare():
    encoded = raw(ROUTE)
    need(sha256(encoded).hexdigest() == ROUTE_SHA, 'exact route controller')
    r = types.ModuleType('finite57_current_adaptive_route')
    r.__file__ = str(ROOT/ROUTE)
    exec(compile(encoded, r.__file__, 'exec'), r.__dict__)
    c = r.prepare()  # Requires the already-held shared parent lock.
    extra = raw(EXTRA)
    need(sha256(extra).hexdigest() == EXTRA_SHA, 'extra manifest identity')
    entries = json.loads(extra, object_pairs_hook=c.unique)['inputs']
    expected = {f'{folder}/{stem}.{ext}' for stem in (ADAPT, INTERFACE)
                for folder, ext in (('scripts','py'), ('results','json'), ('results','md'))}
    need(len(entries) == 6 and {p['path'] for p in entries} == expected,
         'exact six additional adaptive source inputs')
    for p in entries:
        if p['path'] in c.PINS:
            need(c.PINS[p['path']] == p['sha256'], 'base/extra collision')
        c.PINS[p['path']] = p['sha256']
        need(c.digest(p['path']) == p['sha256'], 'extra source/data identity')
        c.READ.add(p['path'])
    return r, c, extra


def typed_current_interface(c):
    """Complete finite domains, ignoring positive guards conservatively."""
    counts = {'sealed_marked_simple_leaf9': 0, 'whole_order_at_most16_path9': 0}
    digest = sha256()
    for kind in ('A','B','D'):
        for order in range(10,59):
            g = c.W.classical_graph(kind, order)
            ds = [c.W._distances(g, i) for i in range(order)]
            for a, b in combinations(range(order), 2):
                c.E.charge(c.FUEL)
                if ds[a][b] != 9:
                    continue
                sealed = any(len(g[m]) == 1 and
                    g[m][next(iter(g[m]))] == 3 and
                    len(g[next(iter(g[m]))]) == 2 for m in (a,b))
                bounded = order <= 16
                if not (sealed or bounded):
                    continue
                value = c.W.finite_same_weight((kind,order), a,b)
                # The final negative guard is included even though its
                # D13 marked distance11 cannot enter this distance9 domain.
                if (kind,order) == ('D',13) and (a,b) in ((0,11),(0,12)):
                    value = min(value,F(1,2))
                need(value == 1, 'CURRENT adaptive typed output is not unit')
                for name, applies in (('sealed_marked_simple_leaf9',sealed),
                                      ('whole_order_at_most16_path9',bounded)):
                    if applies:
                        counts[name] += 1
                        digest.update(c.canonical([name,kind,order,a,b,str(value)]).encode()+b'\n')
    need(all(counts.values()), 'nonempty exact typed interfaces')
    return dict(counts=counts, ordered_domain_sha256=digest.hexdigest(),
        ambient_dimension=58, current_literal_only=True,
        domain='All finite A/B/D whole ranks10..58 at graph distance9; stated sealed-leaf or whole-order<=16 predicate. Exceptional ranks<=8 cannot have distance9.',
        positive_context_witnesses_credited=False, negative_D13_guard_accounted=True)


def adjacency(c, gram, nodes):
    labels = c.E.labels(c.E.matrix(gram), nodes)
    return {v:{nodes[j]:label for j,label in enumerate(labels[i]) if label}
            for i,v in enumerate(nodes)}


def path(adj, marks):
    parents={marks[0]:None};todo=[marks[0]]
    for v in todo:
        for w in sorted(adj[v]):
            if w not in parents:
                parents[w]=v;todo.append(w)
    need(marks[1] in parents, 'connected physical marked component')
    p=[marks[1]]
    while p[-1] != marks[0]:
        p.append(parents[p[-1]])
    return list(reversed(p))


def leaf(record, adj, marks, untouched):
    m=record['leaf'];v=record['neighbor']
    need(record['corner'] == marks and record['distance'] == 9 and
         record['base_order'] == len(adj) >= 10, 'exact leaf/corner association')
    need(m in marks and set(adj[m]) == {v} and adj[m][v] == 3 and
         len(adj[v]) == 2 and m in untouched and v in untouched,
         'sealed marked simple leaf AND sealed degree-two neighbor')
    return dict(leaf=m,neighbor=v,basis='finite57.current.finite_leaf9_degree2.v1')


def distances(adj, starts):
    result={v:0 for v in starts};todo=list(starts)
    for v in todo:
        for w in adj[v]:
            if w not in result:
                result[w]=result[v]+1;todo.append(w)
    need(len(result)==len(adj), 'connected complete physical core')
    return result


def correlated_pairs(rows, left, right, inverse, du, dv):
    need(type(rows) is list and len(rows)<=4096, 'bounded correlated pair list')
    expected={(profile_key(a),profile_key(b)) for a,b in product(left,right)}
    lmap={profile_key(p):p for p in left};rmap={profile_key(p):p for p in right}
    observed=set()
    for row in rows:
        key=(profile_key(row['left']),profile_key(row['right']))
        need(key in expected and key not in observed, 'exact correlated profile-pair domain')
        observed.add(key)
        a,b=row['left'],row['right']
        need(a==lmap[key[0]] and b==rmap[key[1]], 'complete correlated profile association')
        s=sum((inverse[i][j] for i in a['support'] for j in b['support']),F())
        sq=a['scale_square']*b['scale_square']*s*s
        rp=rational(a['residual'])*rational(b['residual'])
        need(s == rational(row['projection_base']) < 0 and
             sq == rational(row['projection_square']) > rp == rational(row['residual_product']),
             'inverse/profile/correlated residual association')
        shortcut=min(min(du[i]+dv[j],du[j]+dv[i])+3
                     for i in a['support'] for j in b['support'])
        need(type(row['shortcut']) is int and 0 < row['shortcut']==shortcut <= 8,
             'accepted pair obstruction has an independent short8 escape')
    need(observed == expected, 'missing correlated profile pair')


def adaptive(c,r):
    high,cores,_=r.domain(c)
    consumers=[p for p in high if p['method']=='nine_context_adaptive_finite_payment']
    cases=c.field(ADAPT,('cases',))
    need(type(cases) is list and len(cases)==9 and
         [p['index'] for p in cases]==list(range(9)) and
         tuple(source_key(p) for p in cases)==CASES, 'complete fixed nine source cases')
    need(len(consumers)==9 and {source_key(p) for p in consumers}==set(CASES),
         'exact high-prefix/source-case bijection')
    interface=c.field(INTERFACE,('adaptive',))
    need(len(interface)==9 and {p['index'] for p in interface}=={p['index'] for p in consumers},
         'exact nine later interface keys')
    indexed={source_key(p):p for p in consumers};interfaces={p['index']:p for p in interface}
    out=[]
    for case in cases:
        c.E.charge(c.FUEL)
        highrow=indexed[source_key(case)];core=cores[highrow['core_key']]
        need(highrow['local_case']==case['index'] and highrow['payment']==case['payment'] and
             core['gram']==case['gram'] and core['inverse']==case['inverse'] and
             core['gauge']==[1]*len(core['gram']), 'exact physical high/adaptive core association')
        marks=case['payment']['corner']
        corners=[p for p in highrow['all_corners'] if p['preferred_corner']==marks]
        need(len(corners)==1, 'one exact chosen physical corner')
        nodes=corners[0]['marked_component'];adj=adjacency(c,core['gram'],nodes)
        pp=path(adj,marks)
        q=interfaces[highrow['index']]
        need(len(pp)==10 and q['marks']==marks and q['retained_path']==pp and q['distance']==9,
             'complete retained path matches later current interface')
        n=len(core['gram']);iv=[[rational(x) for x in row] for row in case['inverse']]
        need(len(iv)==n and all(len(row)==n for row in iv), 'bounded full inverse domain')
        full=adjacency(c,core['gram'],list(range(n)))
        copied={int(k):v for k,v in case['right_copy_map'].items()}
        du=distances(full,case['left'][1]);dv=distances(full,[copied[i] for i in case['right'][1]])
        profiles=case['profiles'];need(0<len(profiles)<=64, 'bounded source profile domain')
        keys=set()
        for p in profiles:
            support=p['support'];k=profile_key(p)
            need(k not in keys and support==sorted(set(support)) and
                 all(type(i) is int and 0<=i<n for i in support) and
                 type(p['scale_square']) is int and p['scale_square'] in (1,2),
                 'unique complete native profile key')
            keys.add(k)
            projection=p['scale_square']*sum((iv[i][j] for i in support for j in support),F())
            need(projection==rational(p['projection']) and
                 2-projection==rational(p['residual'])>=0, 'physical profile/inverse association')
        touched={i for p in profiles for i in p['support']};untouched=set(range(n))-touched
        payment=case['payment'];branches=[];method=case['method']
        if case['index'] in (1,3,5,8):
            need(method=='sealed_marked_leaf', 'sealed four-case dispatch')
            branches.append(leaf(payment,adj,marks,untouched))
        elif case['index'] in (0,4,7):
            need(method=='adaptive_untouched_marked_leaf', 'correlated three-case dispatch')
            alternatives=payment['leaf_alternatives']
            need(len(alternatives)==2 and {p['leaf'] for p in alternatives}==set(marks),
                 'both complete marked-leaf alternatives')
            need(all(not set(marks)<=set(p['support']) for p in profiles), 'one root cannot meet both marks')
            for p in alternatives:
                item=leaf(p,adj,marks,untouched|{p['leaf']})
                item['branch_guard']='No extra supporting white meets this mark. The complete correlated pair exclusions force at least one such alternative unless the whole face has an independent S8 escape.'
                branches.append(item)
            correlated_pairs(payment['pair_exclusions'],
                [p for p in profiles if marks[0] in p['support']],
                [p for p in profiles if marks[1] in p['support']],iv,du,dv)
        else:
            need(case['index'] in (2,6) and method=='adaptive_leaf_or_bounded_relevant_tail',
                 'two-case far/whole-bound dispatch')
            far=[1,2,3] if case['index']==2 else [11,12,13]
            good,bad=(11,13) if case['index']==2 else (1,3)
            need(payment['far_support_nodes']==far and payment['good_node']==good and
                 payment['outside_node']==bad, 'exact far/good/outside node association')
            gs=[p for p in profiles if good in p['support']]
            need(len(gs)==1 and profile_key(gs[0])==((good,),1), 'unique good unit contact')
            item=leaf(payment['if_far_present'],adj,marks,untouched|{good})
            item['branch_guard']='A far-side extra root occurs; complete pair-shortcut exclusions then keep the good marked leaf untouched unless independent S8 applies.'
            branches.append(item)
            correlated_pairs(payment['pair_exclusions'],
                [p for p in profiles if set(p['support'])&set(far)],gs,iv,du,dv)
            active=[p for p in profiles if not set(p['support'])&set(far)]
            need([p for p in active if set(p['support'])&set(nodes)]==gs,
                 'exact remaining marked-component contact domain')
            sealed=set(range(n))-{i for p in active for i in p['support']}
            fixed=[v for v in nodes if len(adj[v])==1 and v in sealed and
                   adj[v][next(iter(adj[v]))]==3 and
                   next(iter(adj[v])) in sealed and len(adj[next(iter(adj[v]))])==2]
            need(fixed==payment['if_far_absent_fixed_simple_leaves'] and bool(fixed) and len(nodes)==12,
                 'whole finite branch has an actual sealed simple leaf')
            need((iv[good][good],iv[bad][bad],iv[good][bad])==(F(5,4),F(1),F(1,2)) and
                 payment['residual_gram']==[['3/4','-1/2'],['-1/2','1']] and
                 payment['residual_inverse']==[['2','1'],['1','3/2']],
                 'exact residual inverse/coordinate association')
            remaining={profile_key(p):p for p in active
                       if p!=gs[0] and profile_key(p)!=((bad,),1)}
            seen=set()
            for check in payment['other_profile_exclusions']:
                p=check['profile'];key=profile_key(p)
                need(key in remaining and key not in seen and p==remaining[key],
                     'exact remaining-profile exclusion domain')
                seen.add(key)
                s=sum((iv[i][good] for i in p['support']),F())
                sq=p['scale_square']*s*s;rp=rational(p['residual'])*rational(gs[0]['residual'])
                need(s>0 and sq==rational(check['projection_square'])>
                     rp==rational(check['residual_product']),
                     'exact remaining-profile/residual incompatibility association')
            need(seen==set(remaining), 'missing remaining-profile incompatibility')
            need(payment['alone_tail_rank_bound']==3 and payment['joint_tail_rank_bound']==1 and
                 payment['whole_component_bound']==16, 'accepted complete whole-bound interface')
            branches.append(dict(whole_finite_order_bound=16,
                basis='finite57.current.finite_whole_order16_path9.v1',
                branch_guard='No extra supporting white contacts the specified far support nodes.',
                geometry_premise='Accepted complete active-profile/first-extra-root/tail proof; not inferred from the displayed12 nodes alone.'))
        out.append(dict(consumer_key=dict(high_index=highrow['index'],method=highrow['method'],
            core_key=highrow['core_key'],local_case=case['index'],selected_marks=marks),
            physical_context_sha256=sha256(c.canonical([core['gram'],nodes,marks]).encode()).hexdigest(),
            retained_path=pp,method=method,current_alternative_branches=branches,
            required_independent_theorem='finite57.current.S8.whole_face.v1',
            independent_S8_verified=False,all_additional_support_geometry='Accepted all-native/profile/whole-bound theorem at its explicitly bound finite data premises.'))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=('adaptive9',),required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args=parser.parse_args(); own=str(Path(__file__).resolve().relative_to(ROOT));source=raw(own)
    r,c,extra=prepare();basis=typed_current_interface(c)
    rows=adaptive(c,r)
    for p in c.READ:need(c.digest(p)==c.PINS[p], 'consumed source/input changed')
    need(sha256(raw(ROUTE)).hexdigest()==ROUTE_SHA,
         'route controller changed during batch')
    need(sha256(raw(r.COLLECTOR)).hexdigest()==r.COLLECTOR_SHA,
         'underlying immutable collector changed during batch')
    need(raw(EXTRA)==extra and raw(own)==source, 'adapter/extra manifest changed during batch')
    need(sha256(c.raw(c.MANIFEST)).hexdigest()==c.MANIFEST_SHA, 'base manifest changed during batch')
    out=dict(status='PASS_CURRENT_ADAPTIVE9_LITERAL_BINDING_WITH_EXPLICIT_S8_PREREQUISITE',
        batch=args.batch,source_sha256=sha256(source).hexdigest(),typed_current_basis=basis,
        current_prefix_providers=rows,consumed_inputs={p:c.PINS[p] for p in sorted(c.READ)},
        geometric_completeness_and_inverse_validity='Explicit accepted source/computed premises; not regenerated.',
        independent_S8_verified=False,current_terminal_bridge_complete=False,
        global_finite57_verified=False,fuel_used=c.FUEL[0])
    payload=json.dumps(out,sort_keys=True)
    need(len(payload.encode())+1<=MAX_OUTPUT, 'complete output cap; never truncate provider keys')
    print(payload)


if __name__=='__main__':main()
