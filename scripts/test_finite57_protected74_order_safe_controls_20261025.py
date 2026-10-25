#!/usr/bin/env python3
"""SOURCE ONLY: six protected-order/dependency rejection controls.

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
SOURCE_SHA='e48fb5b1ccd7d2f0d87831c811cbdcdb3c4b02f928ba95be21c056d4cdb9402e'
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
    need(len(rejected)==6,'six necessary negative controls')
    need(raw(SOURCE)==encoded and raw(own)==own_bytes,'test/adapter source changed')
    print(json.dumps(dict(status='PASS_PROTECTED74_ORDER_ASSOCIATION_CONTROLS_ONLY',
        positive_representation_fixtures=1,negative_controls=rejected,
        source_sha256=SOURCE_SHA,test_sha256=sha256(own_bytes).hexdigest(),
        authoritative_data_read=False,geometric_theorem_tested=False,
        current_distinct_receipts_verified=False,independent_S8_verified=False,
        current_terminal_bridge_complete=False,
        scope='One synthetic complete32+32 key fixture and six required rejections; no enumeration, whole-factor seal or current receiver funding is certified.'),sort_keys=True))


if __name__=='__main__':main()
