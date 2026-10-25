#!/usr/bin/env python3
"""SOURCE ONLY: twenty-two adaptive leaf/pair/key/whole-bound controls.

All matrices and profiles here are representation fixtures, not geometric
examples or Gram-inverse certificates. One fixed187KiB saved producer input
supplies all nine branch fixtures, with synthetic high/interface wrappers.
No legacy worker is read or run. Root exclusively launches after review.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from copy import deepcopy
import fcntl
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import types

ROOT=Path(__file__).resolve().parents[1]
SOURCE='scripts/check_finite57_current_adaptive9_copy_map_safe_20261025.py'
SOURCE_SHA='a62d18107a54c512f7abcb1ceb9b02dd187eed38805928b5f0d0bf096659a0e2'
FIXTURE='results/fresh_distance8_adaptive_leaf_overlap_20261017.json'
FIXTURE_SHA='c738c6af3f910e74b24c17a10c3f38cab2ed3835f4a6661c745046abb98ec8ab'


def need(ok,message):
    if not ok:raise ValueError(message)


def raw(path):
    p=(ROOT/path).resolve()
    need(p.is_relative_to(ROOT) and p.is_file(),'bounded repository source')
    with p.open('rb') as stream:b=stream.read(65537)
    need(0<len(b)<=65536,'bounded source bytes')
    return b


def fixture_bytes():
    p=ROOT/FIXTURE
    with p.open('rb') as stream:b=stream.read(256*1024+1)
    need(len(b)==187026 and sha256(b).hexdigest()==FIXTURE_SHA,'exact bounded nine-case fixture')
    return b


def unique(pairs):
    d={}
    for k,v in pairs:
        need(k not in d,'duplicate fixture JSON key');d[k]=v
    return d


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    parser.parse_args()
    lock=ROOT/'formal/.proof-batches/exclusive.lock'
    with lock.open('r+') as stream:
        try:fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:pass
        else:
            fcntl.flock(stream.fileno(),fcntl.LOCK_UN)
            raise ValueError('the shared parent proof lock must already be held')
    encoded=raw(SOURCE);need(sha256(encoded).hexdigest()==SOURCE_SHA,'exact reviewed adapter source')
    own=str(Path(__file__).resolve().relative_to(ROOT));own_bytes=raw(own)
    sut=types.ModuleType('finite57_adaptive9_association_fixture')
    sut.__file__=str(ROOT/SOURCE)
    exec(compile(encoded,sut.__file__,'exec'),sut.__dict__)
    passed=[]
    def accept(name,fn):
        fn();passed.append(name)
    def refuse(name,fn,message=None):
        try:fn()
        except ValueError as error:
            if message is not None:need(message in str(error),'wrong refusal stage')
            passed.append(name)
        else:raise ValueError('required refusal did not occur: '+name)

    adj={i:{} for i in range(10)}
    for i in range(9):adj[i][i+1]=adj[i+1][i]=3
    leaf=dict(corner=[0,9],leaf=0,neighbor=1,base_order=10,distance=9)
    call=lambda rec=leaf,g=adj,u=set(range(10)):sut.leaf(rec,g,[0,9],u)
    accept('simple_marked_leaf_with_degree_two_neighbor',call)
    wrong=deepcopy(leaf);wrong['neighbor']=2
    refuse('wrong_neighbor_association',lambda:call(wrong))
    doubled=deepcopy(adj);doubled[0][1]=doubled[1][0]=4
    refuse('nonsimple_leaf_is_not_this_interface',lambda:call(leaf,doubled))
    fork=deepcopy(adj);fork[10]={1:3};fork[1][10]=3
    fork_record=deepcopy(leaf);fork_record['base_order']=11
    refuse('degree_three_neighbor_refused',lambda:call(fork_record,fork))
    refuse('unsealed_neighbor_refused',lambda:call(leaf,adj,set(range(10))-{1}))
    reversed_corner=deepcopy(leaf);reversed_corner['corner']=[9,0]
    refuse('ordered_corner_association',lambda:call(reversed_corner))

    left=dict(support=[0],scale_square=1,projection='5/4',residual='3/4')
    right=dict(support=[1],scale_square=1,projection='5/4',residual='3/4')
    iv=[[F(5,4),F(-1)],[F(-1),F(5,4)]]
    pair=dict(left=left,right=right,projection_base='-1',projection_square='1',
              residual_product='9/16',shortcut=3)
    du,dv={0:0,1:10},{0:10,1:0}
    pairs=lambda rows=[pair],inverse=iv,a=du,b=dv:sut.correlated_pairs(rows,[left],[right],inverse,a,b)
    accept('correlated_pair_inverse_and_shortcut_associated',pairs)
    refuse('missing_correlated_pair',lambda:pairs([]))
    refuse('duplicate_correlated_pair',lambda:pairs([pair,pair]))
    altered=deepcopy(pair);altered['left']['residual']='1/4'
    refuse('same_key_with_wrong_profile_values',lambda:pairs([altered]))
    wrong_iv=deepcopy(iv);wrong_iv[0][1]=wrong_iv[1][0]=F(-1,2)
    refuse('same_pair_at_wrong_inverse_entry',lambda:pairs([pair],wrong_iv))
    shortcut=deepcopy(pair);shortcut['shortcut']=4
    refuse('wrong_physical_shortcut_association',lambda:pairs([shortcut]))
    other=deepcopy(pair);other['left']['support']=[1]
    refuse('unsupported_ordered_profile_key',lambda:pairs([other]))
    changed=deepcopy(pair);changed['projection_square']='2'
    refuse('wrong_squared_projection',lambda:pairs([changed]))
    long_pair=deepcopy(pair);long_pair['shortcut']=11
    refuse('actual_shortcut_exceeds_eight',lambda:pairs([long_pair],iv,{0:4,1:10},{0:10,1:4}))

    # These fixtures stop at the new finite-domain join before any
    # geometry/matrix field is requested. They do not simulate a proof.
    cases=[]
    for i,(di,u,dj,v,match) in enumerate(sut.CASES):
        cases.append(dict(index=i,left=[di,list(u)],right=[dj,list(v)],
                          match=[list(x) for x in match]))
    high=[dict(**p,method='nine_context_adaptive_finite_payment') for p in cases]
    interfaces=[dict(index=i) for i in range(9)]
    def domain_call(case_rows,high_rows,interface_rows):
        c=types.SimpleNamespace(field=lambda stem,selector:case_rows if stem==sut.ADAPT else interface_rows)
        r=types.SimpleNamespace(domain=lambda _:(high_rows,{},{}))
        sut.adaptive(c,r)
    refuse('missing_ninth_literal_case',lambda:domain_call(cases[:-1],high,interfaces),
           'complete fixed nine source cases')
    repeated=deepcopy(high);repeated[-1]=deepcopy(repeated[0])
    refuse('high_prefix_bijection_not_just_count',lambda:domain_call(cases,repeated,interfaces),
           'exact high-prefix/source-case bijection')
    wrong_interface=deepcopy(interfaces);wrong_interface[-1]['index']=99
    refuse('later_interface_key_mismatch',lambda:domain_call(cases,high,wrong_interface),
           'exact nine later interface keys')

    # Actual saved matrices/profile values exercise every source branch.
    # The high-index and later-interface wrappers below are synthetic:
    # this tests the adapter's positive schema, not the2614-domain theorem.
    fixed=fixture_bytes()
    saved=json.loads(fixed,object_pairs_hook=unique,
        parse_constant=lambda _:(_ for _ in ()).throw(ValueError('nonfinite fixture constant')))
    saved_cases=saved['cases']
    def complete_fixture(case_rows):
        high_rows=[];cores={};interface_rows=[]
        def fixture_labels(g,nodes):
            out=[]
            for i in nodes:
                row=[]
                for j in nodes:
                    value=F(g[i][j])
                    need(value in (F(2),F(0),F(-1)),'fixture uses only its actual simply-laced entries')
                    row.append(0 if i==j or value==0 else 3)
                out.append(row)
            return out
        for case in case_rows:
            core_key='synthetic-wrapper-'+str(case['index'])
            cores[core_key]=dict(gram=case['gram'],inverse=case['inverse'],gauge=[1]*len(case['gram']))
            high_rows.append(dict(index=case['index'],left=case['left'],right=case['right'],
                match=case['match'],local_case=case['index'],payment=case['payment'],
                all_corners=case['all_corners'],core_key=core_key,
                method='nine_context_adaptive_finite_payment'))
            marks=case['payment']['corner']
            corner=next(p for p in case['all_corners'] if p['preferred_corner']==marks)
            nodes=corner['marked_component'];lab=fixture_labels(case['gram'],nodes)
            adj={v:{nodes[j]:label for j,label in enumerate(lab[i]) if label} for i,v in enumerate(nodes)}
            interface_rows.append(dict(index=case['index'],marks=marks,distance=9,
                                       retained_path=sut.path(adj,marks)))
        fuel=[0]
        def charge(_):
            fuel[0]+=1;need(fuel[0]<=1000,'tiny fixture fuel')
        c=types.SimpleNamespace(field=lambda stem,selector:case_rows if stem==sut.ADAPT else interface_rows,
            canonical=lambda v:json.dumps(v,sort_keys=True,separators=(',',':')),
            E=types.SimpleNamespace(matrix=lambda g:g,labels=fixture_labels,charge=charge),FUEL=fuel)
        r=types.SimpleNamespace(domain=lambda _:(high_rows,cores,{}))
        out=sut.adaptive(c,r)
        need(len(out)==9 and {p['method'] for p in out}=={
            'sealed_marked_leaf','adaptive_untouched_marked_leaf','adaptive_leaf_or_bounded_relevant_tail'},
            'complete nine source branches reached')
        need(all(p['independent_S8_verified'] is False for p in out),'S8 open prerequisite survives every provider')
        need(sum('whole_finite_order_bound' in b for p in out for b in p['current_alternative_branches'])==2,
             'both complete whole-bound branches reached')
    accept('complete_nine_case_positive_fixture_all_three_methods',lambda:complete_fixture(deepcopy(saved_cases)))
    wrong_bound=deepcopy(saved_cases);wrong_bound[2]['payment']['whole_component_bound']=17
    refuse('whole_bound_seventeen_refused_at_bound_interface',lambda:complete_fixture(wrong_bound),
           'accepted complete whole-bound interface')
    def dictionary_copy_fixture():
        rows=deepcopy(saved_cases)
        for case in rows:
            indexed=case['right_copy_map']
            need(type(indexed) is list,'saved producer uses an indexed copy-map list')
            mapping={str(i):v for i,v in enumerate(indexed)}
            old_order=sut.NATIVE_ORDERS[case['right'][0]]
            n=len(case['gram'])
            need(sut.copy_map(indexed,old_order,n)==sut.copy_map(mapping,old_order,n),
                 'list/dictionary old-node copy-map equivalence')
            case['right_copy_map']=mapping
        complete_fixture(rows)
    accept('complete_nine_case_dictionary_copy_map_schema_equivalence',dictionary_copy_fixture)
    refuse('incomplete_integer_copy_map_domain_refused',
           lambda:sut.copy_map({'0':0,'1':1,'3':3},3,4),
           'complete integer old-node copy-map domain')
    need(len(passed)==22,'complete twenty-two necessary association controls')
    need(raw(SOURCE)==encoded and raw(own)==own_bytes,'test/adapter source changed')
    need(fixture_bytes()==fixed,'saved fixture changed')
    print(json.dumps(dict(status='PASS_ADAPTIVE9_ASSOCIATION_CONTROLS_ONLY',tests=passed,
        count=len(passed),source_sha256=SOURCE_SHA,test_sha256=sha256(own_bytes).hexdigest(),
        saved_fixture_sha256=FIXTURE_SHA,synthetic_high_and_later_wrappers=True,
        authoritative_high_domain_read=False,geometric_theorem_tested=False,
        independent_S8_verified=False,current_terminal_bridge_complete=False,
        scope='Twenty-two schema/association controls: the original twenty plus complete saved nine-case list/dictionary copy-map equivalence and an incomplete integer-domain refusal. No current typed-rule enumeration, inverse-validity or whole-support theorem is tested.'),sort_keys=True))


if __name__=='__main__':main()
