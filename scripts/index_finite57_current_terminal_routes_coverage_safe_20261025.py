#!/usr/bin/env python3
"""SOURCE ONLY: complete current marked-completion obligation index (G6).

This source emits exact unbound obligations. It does NOT discharge them and
does NOT emit a full-terminal PASS. Root alone may launch a named shard under
the existing bounded runner after independent source clearance.
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
import json
from pathlib import Path
import types

ROOT=Path(__file__).resolve().parents[1]
COLLECTOR='scripts/check_finite57_terminal_collectors_coverage_safe_20261025.py'
COLLECTOR_SHA='68f7a441f4e9701090e30d9a049430a23a3e5cae660e86d1f248a63a26a87dee'
SHARDS=64
MAX_OUTPUT=500*1024


def need(ok,message):
    if not ok:raise ValueError(message)


def raw(relative,cap=2*1024**2):
    path=(ROOT/relative).resolve()
    need(path.is_relative_to(ROOT)and path.is_file(),'source path')
    with path.open('rb')as stream:data=stream.read(cap+1)
    need(0<len(data)<=cap,'bounded source read')
    return data


def distance(labels,a,b):
    dist={a:0};queue=[a]
    for i in queue:
        for j,label in enumerate(labels[i]):
            if label and j not in dist:dist[j]=dist[i]+1;queue.append(j)
    need(b in dist,'connected marked component required')
    return dist[b]


def embeddings(pattern,target,fuel,engine):
    """Every induced labelled injection; degree is a necessary lower filter."""
    n,m=len(pattern),len(target)
    order=sorted(range(n),key=lambda i:(-sum(bool(x)for x in pattern[i]),i))
    used,assignment=set(),{}
    def visit(depth):
        engine.charge(fuel)
        if depth==n:
            yield [assignment[i]for i in range(n)]
            return
        i=order[depth];degree=sum(bool(x)for x in pattern[i])
        for j in range(m):
            if j in used or degree>len(target[j]):continue
            if any(pattern[i][a]!=target[j].get(b,0)for a,b in assignment.items()):continue
            used.add(j);assignment[i]=j
            yield from visit(depth+1)
            del assignment[i];used.remove(j)
    yield from visit(0)


def prepare():
    encoded=raw(COLLECTOR);need(sha256(encoded).hexdigest()==COLLECTOR_SHA,'collector source identity')
    c=types.ModuleType('finite57_route_collector');c.__file__=str(ROOT/COLLECTOR)
    exec(compile(encoded,str(ROOT/COLLECTOR),'exec'),c.__dict__)
    c.parent_lock()
    manifest=c.raw(c.MANIFEST)
    need(sha256(manifest).hexdigest()==c.MANIFEST_SHA,'collector manifest identity')
    doc=json.loads(manifest,object_pairs_hook=c.unique)
    c.PINS={r['path']:r['sha256']for r in doc['inputs']}
    need(len(c.PINS)==len(doc['inputs']),'unique input paths')
    c.W=c.module(c.LIB,'finite57_route_literal')
    c.E=c.module(c.ENGINE,'finite57_route_engine')
    c.J=c.module(c.PROJECTOR,'finite57_route_projection')
    c.J.json=types.SimpleNamespace(loads=lambda x:json.loads(x,object_pairs_hook=c.unique,
        parse_constant=lambda x:(_ for _ in()).throw(ValueError('nonfinite JSON number'))))
    primary='manuscripts/finite_volume57.tex'
    need(c.digest(primary)==c.PINS[primary]==c.W.TEX_SHA256,'literal/primary source binding')
    c.READ.add(primary)
    return c


def domain(c):
    doc=c.document(c.HIGH,('isolation','cases'),('isolation','cores'))['isolation']
    cases=doc['cases'];need(type(cases)is list and len(cases)==2614,'complete2614 high domain')
    need([p['index']for p in cases]==list(range(2614)),'exact authoritative high row indices')
    counts={}
    for p in cases:counts[p['method']]=counts.get(p['method'],0)+1
    need(counts=={'radical':1329,'lorentz_residual':823,'relevant_component':259,
                 'marked_component_isolated_by_end_inverses':120,'finite_leaf9_or_ideal_feature5':74,
                 'nine_context_adaptive_finite_payment':9},'exact six-method high partition')
    return cases,doc['cores'],counts


def selected(shard,c):
    cases,cores,counts=domain(c);out=[];direct=sha256();maps=0;paid=0;prefix=[]
    for p in cases:
        if p['index']%SHARDS!=shard:continue
        if p['method']in('finite_leaf9_or_ideal_feature5','nine_context_adaptive_finite_payment'):
            prefix.append(dict(index=p['index'],method=p['method'],
                required_provider='G4/G6 independent protected/adaptive current-output binding'))
            continue
        h=p['component_order_bound'];nodes=p['marked_component'];marks=p['preferred_corner']
        need(type(h)is int and 0<=h<=17,'whole-completion bound type/range')
        need(len(marks)==2 and set(marks)<=set(nodes),'both selected marks lie in the bounded component')
        core=cores[p['core_key']];g=c.E.matrix(core['gram']);lab=c.E.labels(g,nodes)
        c.E.indices(nodes,len(g));c.E.indices(marks,len(g),2)
        local=[nodes.index(m)for m in marks]
        if distance(lab,*local)<=8:
            prefix.append(dict(index=p['index'],method='independent_S8',
                required_provider='Current S8 whole-face theorem; not a same-corner unit assertion'))
            continue
        need(len(nodes)>=10,'no-short8 connected component has at least ten nodes')
        if h<len(nodes):
            prefix.append(dict(index=p['index'],method='display_exceeds_proved_whole_bound',
                required_provider='Exact generic extension bound certificate; impossible actual component'))
            continue
        for affine,kinds in((False,('A','B','D')),(True,('A','B','C','D'))):
            for order in range(len(nodes),h+1):
                for kind in kinds:
                    target=c.W.classical_graph(kind,order,affine=affine)
                    for mp in embeddings(lab,target,c.FUEL,c.E):
                        pair=tuple(sorted(mp[i]for i in local));maps+=1
                        weight=c.W.affine_weight(kind,order,*pair)if affine else c.W.finite_same_weight((kind,order),*pair)
                        conditional=False
                        if not affine and(kind,order)==('D',13)and pair in((0,11),(0,12)):
                            # This is the guard-TRUE obligation. Guard-FALSE has
                            # ordinary unit1. No displayed spectator capacity
                            # is used to rule out the guard.
                            weight=F(1,2);conditional=True
                        key=dict(high_index=p['index'],core_key=p['core_key'],base_nodes=nodes,
                            selected_marks=marks,target_kind=kind,target_order=order,
                            affine=affine,induced_map=mp,marked_image=list(pair))
                        if weight==1:
                            paid+=1;direct.update(c.canonical(key).encode()+b'\n')
                        else:
                            need(0<=weight<1,'literal target weight domain')
                            out.append(dict(key=key,intrinsic_current_lower=str(weight),
                                conditional_negative_D13_guard=conditional,
                                required_link='G5 guarded selected exclusion'if conditional else
                                    'exact current G1/G2/G3/G4 whole receipt or exact proper-completion exclusion',
                                binding_status='UNBOUND'))
    return dict(shard=shard,shards=SHARDS,complete_authoritative_partition=counts,
        all_induced_maps_visited=maps,direct_whole_factor_unit_maps=paid,
        direct_map_key_sha256=direct.hexdigest(),prefix_obligations=prefix,
        unbound_current_completion_obligations=out,
        scope='Complete named-shard obligation INDEX. It does not discharge a replacement, its geometry or its computed premise.')


def index(c):
    cases,cores,counts=domain(c)
    return dict(authoritative_rows=2614,exact_method_partition=counts,
        shards=SHARDS,shard_keys=[list(range(i,2614,SHARDS))for i in range(SHARDS)],
        complete_finite_target_families=['A','B','D'],complete_affine_target_families=['A','B','C','D'],
        whole_order_bound=17,exceptional_exclusion='Connected marked path length>=9 needs>=10 nodes; exceptional finite/affine types are smaller.',
        proof_link_contract=dict(full_key=['high_index','core_key','base_nodes','selected_marks','target_kind','target_order','affine','induced_map','marked_image'],
            permitted_links=['current whole-case receipt with exact complete variant domain',
                             'exact proper-map exclusion outcome at the same physical core and target map',
                             'named independently proved universal priority at this exact affine pair',
                             'final guard-TRUE selected exclusion with guard-FALSE ordinary unit'],
            forbidden=['lookup by count alone','historical paid-pair inclusion','arbitrary A1 padding',
                       'an inherited lower called current equality','a proper exclusion taken from a guard-dependent theorem without its guard']),
        status='INDEX_ONLY_FINAL_LINK_TABLE_NOT_PREPARED')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=('index',)+tuple(f'selected-{i:02d}'for i in range(SHARDS)),required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args=parser.parse_args()
    relative=str(Path(__file__).resolve().relative_to(ROOT));source=raw(relative)
    c=prepare()
    out=index(c)if args.batch=='index'else selected(int(args.batch[-2:]),c)
    for path in c.READ:need(c.digest(path)==c.PINS[path],'consumed pinned source/input changed')
    need(sha256(raw(COLLECTOR)).hexdigest()==COLLECTOR_SHA,'collector changed during route index')
    need(raw(relative)==source,'route-index source changed during batch')
    need(sha256(c.raw(c.MANIFEST)).hexdigest()==c.MANIFEST_SHA,'collector manifest changed during route index')
    payload=json.dumps(dict(status='CURRENT_TERMINAL_OBLIGATION_INDEX_ONLY',batch=args.batch,
        source_sha256=sha256(source).hexdigest(),collector_sha256=COLLECTOR_SHA,
        consumed_inputs={p:c.PINS[p]for p in sorted(c.READ)},fuel_used=c.FUEL[0],result=out),sort_keys=True)
    need(len(payload.encode())<=MAX_OUTPUT,'output cap: refuse rather than truncate obligation keys')
    print(payload)


if __name__=='__main__':main()
