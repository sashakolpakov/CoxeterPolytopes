#!/usr/bin/env python3
"""SOURCE ONLY: twelve protected order/dependency/structure controls.

The positive and altered records are tiny representation fixtures. They
test exact32+32 receiver association and the rank11/order12 interface;
they do not prove the geometric whole-face theorem or current funding.
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
from itertools import product
import json
from pathlib import Path
import types

ROOT=Path(__file__).resolve().parents[1]
SOURCE='scripts/index_finite57_current_protected74_order_safe_applications_20261025.py'
SOURCE_SHA='957b1ef8c332de9aa5af142b9574d3ec625887bc6b57eb56c7e5c80ffacc4ab5'
MODULES=(
 ('a','scripts/check_finite57_current_adaptive9_path_binding_20261025.py',
  '76022d95fcc5ec6e61bf1ccaea65635adda1ab8593d0bf921d75754f14bd139c'),
 ('r','scripts/index_finite57_current_terminal_routes_syntax_safe_20261025.py',
  '43a0c4e986cbdb92d112771bf0edc1d128ce2d0feced30fffd9d9825dda7c29d'),
 ('w','scripts/finite57_literal_weights_audit_20261024.py',
  'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0'),
 ('e','scripts/finite57_terminal_literal_engine_small_affine_20261025.py',
  '39dd55d2b2abef3db045e646edcc95bc1c8363f483d9f853d8c1d4e1b65d0e1f'),
)
PRODUCER='fresh_affine_d11_leaf_pair_transfer_20261017'
REPAIR='fresh_d12_universal_ideal_output_repairs_20261018'
PROVIDER='finite57.current.affine_order12_D_opposite_leaves.whole_face.v1'


def need(ok,message):
    if not ok:raise ValueError(message)


def raw(path):
    p=(ROOT/path).resolve()
    need(p.is_relative_to(ROOT) and p.is_file(),'bounded repository source')
    with p.open('rb') as stream:b=stream.read(65537)
    need(0<len(b)<=65536,'bounded source bytes')
    return b


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    parser.parse_args()
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('r+') as stream:
        try:fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:pass
        else:
            fcntl.flock(stream.fileno(),fcntl.LOCK_UN)
            raise ValueError('the shared parent proof lock must already be held')
    encoded=raw(SOURCE)
    need(sha256(encoded).hexdigest()==SOURCE_SHA,'exact import-repaired protected source')
    own=str(Path(__file__).resolve().relative_to(ROOT));own_bytes=raw(own)
    sut=types.ModuleType('finite57_protected74_order_fixture')
    sut.__file__=str(ROOT/SOURCE)
    exec(compile(encoded,sut.__file__,'exec'),sut.__dict__)

    marks=list(product((8,9),(10,11)))
    original={
        (PRODUCER,('removed_pairs',)):[list(p) for p in marks],
        (REPAIR,('affine_D11','rows')):[
            dict(old_marks=list(p),spectator_free_extreme_side=side,branch=branch)
            for p,side,branch in product(marks,('left','right'),('isolated4','mate3_1','D5','A1_A3'))],
        (PRODUCER,('two_A1_A3_spectators','rows')):[
            dict(marks=list(p),leaf_attachments=[i,j],A1_contacting_black=side)
            for p in marks for i in (8,9) for j in (10,11) for side in (12,13)],
    }
    def fixture(data=None,value=F(0),distance=9):
        data=original if data is None else data
        def field(stem,selector):
            need((stem,selector) in data,'unexpected producer/selector association')
            return data[stem,selector]
        def graph(kind,order,affine=False):
            need((kind,order,affine)==('D',12,True),'wrong rank/order graph interface')
            return None
        def weight(kind,order,*pair):
            need((kind,order)==('D',12) and pair in marks,'wrong marked zero domain')
            return value
        def distances(_,start):
            need(start in (8,9),'wrong first marked vertex')
            return {10:distance,11:distance}
        return types.SimpleNamespace(field=field,W=types.SimpleNamespace(
            affine_weight=weight,classical_graph=graph,_distances=distances))
    baseline=sut.order12_dependencies(fixture())
    need(baseline['producer']==PRODUCER and baseline['whole_affine_order']==12 and
         baseline['affine_rank']==11 and baseline['theorem']==PROVIDER and
         baseline['current_distinct_receipts_verified'] is False,
         'positive rank/order metadata and open funding premise')
    need(len(baseline['exact_current_repair_keys'])==32 and
         len(baseline['exact_both_nonempty_keys'])==32,
         'both complete receiver key domains are retained')
    need(sut.PRIORITY_SOURCES[PROVIDER]==PRODUCER and
         all(sut.PRIORITIES['D',12,*p]==PROVIDER for p in marks),
         'exact protected request-to-source association')

    rejected=[]
    def refuse(name,fn,message):
        try:fn()
        except ValueError as error:
            need(message in str(error),'wrong refusal stage: '+name)
            rejected.append(name)
        else:raise ValueError('required refusal did not occur: '+name)

    wrong=deepcopy(original)
    wrong[PRODUCER,('removed_pairs',)]=[list(p) for p in product((9,10),(11,12))]
    refuse('rank12_order13_pairs_are_not_rank11_order12',
           lambda:sut.order12_dependencies(fixture(wrong)),'exact rank11/order12/distance9 zero domain')
    refuse('nonzero_current_value_is_not_the_zero_interface',
           lambda:sut.order12_dependencies(fixture(value=F(1,2))),
           'exact rank11/order12/distance9 zero domain')
    refuse('distance_ten_is_not_the_order12_distance9_interface',
           lambda:sut.order12_dependencies(fixture(distance=10)),
           'exact rank11/order12/distance9 zero domain')
    wrong=deepcopy(original)
    wrong[REPAIR,('affine_D11','rows')][0]['old_marks']=[8,12]
    refuse('current_repair_key_associated_with_wrong_marked_pair',
           lambda:sut.order12_dependencies(fixture(wrong)),
           'exact current32 D11 repair keys')
    wrong=deepcopy(original)
    wrong[REPAIR,('affine_D11','rows')]=wrong[REPAIR,('affine_D11','rows')][:24]
    refuse('quiet24_cannot_replace_current32_repair_domain',
           lambda:sut.order12_dependencies(fixture(wrong)),
           'exact current32 D11 repair keys')
    wrong=deepcopy(original)
    wrong[PRODUCER,('two_A1_A3_spectators','rows')][0]['leaf_attachments']=[8,12]
    refuse('both_nonempty_receiver_has_wrong_attachment_association',
           lambda:sut.order12_dependencies(fixture(wrong)),
           'exact separate32 both-nonempty A1+A3 receivers')
    need(len(rejected)==6,'six necessary order-association controls')

    # Read exact helper source bytes; no preparation/main entry is called.
    # Fixtures below certify program rejection semantics, not geometric
    # realizability or the complete74 data domain.
    modules={};module_bytes={}
    for name,path,digest in MODULES:
        code=raw(path);need(sha256(code).hexdigest()==digest,'exact structural helper source')
        module=types.ModuleType('finite57_protected_fixture_'+name)
        module.__file__=str(ROOT/path)
        exec(compile(code,module.__file__,'exec'),module.__dict__)
        modules[name]=module;module_bytes[path]=code
    a,r,w,e=(modules[name] for name in ('a','r','w','e'))
    def context(double=False):
        gram=[['2' if i==j else '0' for j in range(10)] for i in range(10)]
        for i in range(9):gram[i][i+1]=gram[i+1][i]='-sqrt(2)' if double and i==8 else '-1'
        inverse=[['1' if i==j else '0' for j in range(10)] for i in range(10)]
        forbidden={9} if double else {0,1}
        row=dict(index=0,method='finite_leaf9_or_ideal_feature5',core_key='fixture',
            marked_component=list(range(10)),preferred_corner=[0,9],distance=9,base_order=10,
            possible_contacts=[i for i in range(10) if i not in forbidden],
            protected_marked_leaves=[9] if double else [0],
            protected_simple_leaves=[] if double else [0],
            feature_distance=0 if double else 999,forced_finite=not double,ideal_basic_pair=double)
        core=dict(gram=gram,inverse=inverse,gauge=['1']*10)
        return row,core
    def protected_fixture(row,core,*,duplicate=False,sparse=False):
        # Exactly one synthetic row lies in shard0; the73 additional unique
        # wrappers lie in shard1 and supply no graph or mathematical claim.
        rows=[row]+[dict(index=4*i+1,method=row['method']) for i in range(73)]
        interfaces=[dict(index=p['index'],marks=[0,9],retained_path=list(range(10)),distance=9) for p in rows]
        target=w.classical_graph('B',13,affine=True) if sparse else w.classical_graph('A',10)
        mapping=list(range(1,11)) if sparse else list(range(10))
        def maps(pattern,candidate,fuel,engine):
            if candidate==target:
                yield mapping[:]
                if duplicate:yield mapping[:]
        route=types.SimpleNamespace(domain=lambda _:(rows,{'fixture':core},{}),embeddings=maps)
        c=types.SimpleNamespace(field=lambda stem,selector:interfaces,W=w,E=e,FUEL=[0],
            canonical=lambda value:json.dumps(value,sort_keys=True,separators=(',',':')))
        return sut.protected(0,a,route,c)
    row,core=context()
    baseline=protected_fixture(deepcopy(row),deepcopy(core))
    need(len(baseline)==1 and baseline[0]['direct_unit_full_maps']==1 and
         baseline[0]['exact_small_current_priority_requests']==[] and
         baseline[0]['independent_S8_verified'] is False,
         'positive sealed-path synthetic wrapper retains one current unit and open S8')
    wrong=deepcopy(row);wrong['possible_contacts'].append(0)
    refuse('falsely_claimed_sealed_marked_leaf',
           lambda:protected_fixture(wrong,deepcopy(core)),
           'exact protected leaf/neighbor/feature seals')
    wrong=deepcopy(row);wrong['possible_contacts'].append(1)
    refuse('degree_two_neighbor_not_sealed',
           lambda:protected_fixture(wrong,deepcopy(core)),
           'exact protected leaf/neighbor/feature seals')

    # Exercise the ACTUAL induced-injection routine rather than a stub.
    # Path3 has a nonedge between its endpoints: a triangle has none.
    pattern=[[0,3,0],[3,0,3],[0,3,0]]
    triangle=[{1:3,2:3},{0:3,2:3},{0:3,1:3}]
    path4=w.classical_graph('A',4)
    need(len(list(r.embeddings(pattern,path4,[0],e)))==4,
         'positive strict partial-injection fixture retains all four full maps')
    need(list(r.embeddings(pattern,triangle,[0],e))==[],
         'a full map cannot send the required endpoint nonedge to an edge')
    wrong_label=[[0,4,0],[4,0,3],[0,3,0]]
    need(list(r.embeddings(wrong_label,path4,[0],e))==[],
         'a full map cannot alter an actual native edge label')
    rejected.append('wrong_full_map_nonedge_and_native_label')
    refuse('duplicate_complete_protected_map_key',
           lambda:protected_fixture(deepcopy(row),deepcopy(core),duplicate=True),
           'unique complete protected target map')
    sparse_row,sparse_core=context(double=True)
    need(w.affine_weight('B',13,1,10)<1 and ('B',13,1,10) not in sut.PRIORITIES,
         'sparseB13 is not the universal inner-terminal domain')
    refuse('sparseB13_is_not_misrouted_to_universal_BC13_inner',
           lambda:protected_fixture(sparse_row,sparse_core,sparse=True),
           'unsupported CURRENT protected output')

    # Limit only this SYNTHETIC interface test to whole order20. The
    # production large_interface source and its full18..58 domain stay
    # unchanged. This fixture still uses actual order20 classical graphs.
    native_range=range
    def fixture_range(*args):
        return (20,) if args==(18,59) else native_range(*args)
    sut.range=fixture_range
    def large(weight):
        proxy=types.SimpleNamespace(**vars(w));proxy.finite_same_weight=weight
        c=types.SimpleNamespace(W=proxy,E=e,FUEL=[0],
            canonical=lambda value:json.dumps(value,sort_keys=True,separators=(',',':')))
        return sut.large_interface(c)
    ordinary=w.finite_same_weight
    large_baseline=large(ordinary)
    need(large_baseline['counts']['finite_D_short_leaf_halves']==2 and
         large_baseline['largeD_current_provider_verified'] is False,
         'both current D20 short-leaf roles retain only half and an open whole-face provider')
    def promoted(typ,p,q):
        if typ==('D',20) and {p,q}=={9,18}:return F(1)
        return ordinary(typ,p,q)
    refuse('largeD_short_leaf_half_cannot_be_promoted_to_unit',
           lambda:large(promoted),'CURRENT largeD short-leaf9 must retain its half')
    del sut.range
    need(len(rejected)==12,'six order and six structural controls')
    for path,code in module_bytes.items():need(raw(path)==code,'structural helper changed')
    need(raw(SOURCE)==encoded and raw(own)==own_bytes,'test/adapter source changed')
    print(json.dumps(dict(status='PASS_PROTECTED74_ORDER_AND_STRUCTURE_CONTROLS_ONLY',
        positive_representation_fixtures=4,negative_controls=rejected,
        source_sha256=SOURCE_SHA,test_sha256=sha256(own_bytes).hexdigest(),
        authoritative_data_read=False,geometric_theorem_tested=False,
        current_distinct_receipts_verified=False,independent_S8_verified=False,
        current_terminal_bridge_complete=False,
        helper_source_sha256={path:digest for _,path,digest in MODULES},
        scope='Synthetic complete32+32 key fixture plus six order and six structural controls. Actual pinned injection/graph helpers are used on tiny fixtures; no full data enumeration, whole-factor seal or current receiver funding is certified.'),sort_keys=True))


if __name__=='__main__':main()
