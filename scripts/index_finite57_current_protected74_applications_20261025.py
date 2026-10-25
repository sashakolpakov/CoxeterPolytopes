#!/usr/bin/env python3
"""SOURCE ONLY: exact protected74 prefix applications and priority requests.

Current units are bound directly. A reduced intrinsic value produces an
exact current whole-face priority request, never a fictitious corner unit.
Independent S8 and those current providers remain explicitly OPEN.
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
from itertools import combinations
import json
from pathlib import Path
import types

ROOT=Path(__file__).resolve().parents[1]
LOADER='scripts/check_finite57_current_adaptive9_path_binding_20261025.py'
LOADER_SHA='76022d95fcc5ec6e61bf1ccaea65635adda1ab8593d0bf921d75754f14bd139c'
PRIORITY_INPUTS='results/finite57_current_protected74_priority_inputs_20261025.json'
PRIORITY_INPUTS_SHA='21e30762c875cdfd24060fb767b85f6d0e73103c540088358011d767b28b3401'
MAX_OUTPUT=500*1024
SHARDS=4
PRIORITY_SOURCES={
 'finite57.current.C10_endpoint_zero.whole_face.v1':'fresh_order10_c_endpoint_zero_transfer_20261018',
 'finite57.current.BC11_zero.whole_face.v1':'fresh_order11_bc_zero_transfer_20261018',
 'finite57.current.B12_center_terminal.whole_face.v1':'fresh_order12_b_center_terminal_transfer_20261018',
 'finite57.current.BC12_near_terminal.whole_face.v1':'fresh_order12_bc_near_terminal_half_transfer_20261018',
 'finite57.current.D12_opposite_leaves.whole_face.v1':'fresh_affine_d12_leaf_pair_transfer_20261017',
 'finite57.current.BC13_inner_terminal.whole_face.v1':'fresh_order13_bc_inner_terminal_transfer_20261018',
 'finite57.current.D13_center_leaf.whole_face.v1':'fresh_order13_d_center_leaf_half_transfer_20261018',
 'finite57.current.D14_near_fork.whole_face.v1':'fresh_order14_d_near_fork_half_transfer_20261018',
 'finite57.current.D15_basic_fractional.whole_face.v1':'fresh_order15_bd_basic_fractional_transfer_20261018',
 'finite57.current.largeD_leaf9_distinct_half.whole_face.v1':'fresh_large_finite_d_leaf_half_transfer_20261017',
}

# Exact all-native whole-face domains, NOT intrinsic units or full-map
# claims about selected-only sparse/projection packets.
PRIORITIES={
 ('C',10,0,9):'finite57.current.C10_endpoint_zero.whole_face.v1',
 **{('B',11,8,j):'finite57.current.BC11_zero.whole_face.v1' for j in (9,10)},
 **{('C',11,*p):'finite57.current.BC11_zero.whole_face.v1' for p in ((0,9),(0,10),(1,10))},
 ('B',12,0,9):'finite57.current.B12_center_terminal.whole_face.v1',
 **{('B',12,8,j):'finite57.current.BC12_near_terminal.whole_face.v1' for j in (10,11)},
 **{('D',12,i,j):'finite57.current.D12_opposite_leaves.whole_face.v1'
    for i in (8,9) for j in (10,11)},
 **{('B',13,8,j):'finite57.current.BC13_inner_terminal.whole_face.v1' for j in (11,12)},
 **{('D',13,*p):'finite57.current.D13_center_leaf.whole_face.v1'
    for p in ((0,11),(0,12),(8,9),(8,10))},
 **{('D',14,*p):'finite57.current.D14_near_fork.whole_face.v1'
    for p in ((1,12),(1,13),(8,10),(8,11))},
 **{('D',15,*p):'finite57.current.D15_basic_fractional.whole_face.v1'
    for p in ((2,13),(2,14),(8,11),(8,12))},
}


def need(ok,message):
    if not ok:raise ValueError(message)


def raw(name):
    p=(ROOT/name).resolve()
    need(p.is_relative_to(ROOT) and p.is_file(),'bounded repository source')
    with p.open('rb') as stream:b=stream.read(4*1024**2+1)
    need(0<len(b)<=4*1024**2,'bounded source bytes')
    return b


def prepare():
    encoded=raw(LOADER)
    need(sha256(encoded).hexdigest()==LOADER_SHA,'exact shared adaptive/current loader')
    a=types.ModuleType('finite57_protected_current_shared');a.__file__=str(ROOT/LOADER)
    exec(compile(encoded,a.__file__,'exec'),a.__dict__)
    r,c,extra=a.prepare()
    priority=raw(PRIORITY_INPUTS)
    need(sha256(priority).hexdigest()==PRIORITY_INPUTS_SHA,'exact priority source register')
    entries=json.loads(priority,object_pairs_hook=c.unique)['inputs']
    expected={f'{folder}/{stem}.{ext}' for stem in PRIORITY_SOURCES.values()
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))}
    need(len(entries)==len(expected)==30 and {p['path'] for p in entries}==expected,
         'exact thirty named priority proof/source/data identities')
    for p in entries:
        if p['path'] in c.PINS:need(c.PINS[p['path']]==p['sha256'],'priority/base identity collision')
        c.PINS[p['path']]=p['sha256']
        need(c.digest(p['path'])==p['sha256'],'bound current-priority proof/source premise')
        c.READ.add(p['path'])
    return a,r,c,extra,priority


def large_interface(c):
    """Complete ambient58 typed domain; no artificial spectator credits."""
    digest=sha256();counts={'finite_units':0,'finite_D_short_leaf_halves':0,'ideal_units':0}
    for ideal in (False,True):
        for kind in (('B','C','D') if ideal else ('A','B','D')):
            for order in range(18,59):
                g=c.W.classical_graph(kind,order,affine=ideal)
                ds=[c.W._distances(g,i) for i in range(order)]
                features=c.W._features(g,finite=not ideal)
                for p,q in combinations(range(order),2):
                    c.E.charge(c.FUEL)
                    if ds[p][q]!=9 or not any(len(g[m])==1 for m in (p,q)):continue
                    if ideal:
                        if not features or min(ds[m][f] for m in (p,q) for f in features)>5:continue
                        value=c.W.affine_weight(kind,order,p,q)
                        need(value==1,'large protected ideal feature5 CURRENT unit')
                        name='ideal_units'
                    else:
                        value=c.W.finite_same_weight((kind,order),p,q)
                        if value==F(1,2):
                            need(kind=='D' and order>=20 and
                                 {p,q} in ({order-11,order-2},{order-11,order-1}),
                                 'only the exact largeD short-leaf9 reduction is admitted')
                            name='finite_D_short_leaf_halves'
                        else:
                            need(value==1,'large protected finite leaf9 CURRENT value')
                            name='finite_units'
                    counts[name]+=1
                    digest.update(c.canonical([name,kind,order,p,q,str(value)]).encode()+b'\n')
    need(all(counts.values()),'complete nonempty three large-current domains')
    return dict(counts=counts,ordered_domain_sha256=digest.hexdigest(),actual_ambient=58,
        finite_D_half_requires='finite57.current.largeD_leaf9_distinct_half.whole_face.v1',
        largeD_current_provider_verified=False,
        domain='Every finite rank18..58 and affine B/C/D order18..58 with sealed marked leaf9; ideal cases additionally have a persistent feature within5. AffineA is excluded by its lack of leaves.')


def protected(shard,a,r,c):
    cases,cores,_=r.domain(c)
    chosen=[p for p in cases if p['method']=='finite_leaf9_or_ideal_feature5']
    interfaces=c.field(a.INTERFACE,('protected',))
    need(len(chosen)==74 and len(interfaces)==74 and
         {p['index'] for p in chosen}=={p['index'] for p in interfaces},
         'exact complete74 high/later interface keys')
    later={p['index']:p for p in interfaces};out=[]
    for p in chosen:
        if p['index']%SHARDS!=shard:continue
        c.E.charge(c.FUEL)
        core=cores[p['core_key']];g=c.E.matrix(core['gram'])
        nodes=p['marked_component'];marks=p['preferred_corner'];n=len(g)
        c.E.indices(nodes,n);c.E.indices(marks,n,2)
        adj=a.adjacency(c,core['gram'],nodes);pp=a.path(adj,marks)
        record=later[p['index']]
        need(len(pp)==10 and record['marks']==marks and record['retained_path']==pp and
             record['distance']==p['distance']==9 and p['base_order']==len(nodes),
             'exact protected physical path/later association')
        allowed=set(p['possible_contacts'])
        need(allowed<=set(range(n)) and core['inverse'] is not None and len(core['gauge'])==n,
             'complete physical allowed-contact/inverse domain')
        forbidden={v for v in nodes if v not in allowed or
                   a.rational(core['gauge'][v])*a.rational(core['inverse'][v][v])>2}
        leaves=[v for v in nodes if len(adj[v])==1 and v in forbidden]
        simple=[v for v in leaves if next(iter(adj[v])) in forbidden and
                len(adj[next(iter(adj[v]))])==2 and adj[v][next(iter(adj[v]))]==3]
        marked=[v for v in marks if v in leaves]
        distances={m:a.distances(adj,[m]) for m in marks}
        features=[v for v in nodes if len(adj[v])>=3 or any(label!=3 for label in adj[v].values())]
        near=min((distances[m][f] for m in marks for f in features),default=999)
        need(type(p['forced_finite']) is bool and type(p['ideal_basic_pair']) is bool and
             marked==p['protected_marked_leaves'] and bool(marked) and
             simple==p['protected_simple_leaves'] and near==p['feature_distance'] and
             p['forced_finite']==bool(simple and len(nodes)>=10) and
             p['ideal_basic_pair']==(near<=5) and (p['forced_finite'] or near<=5),
             'exact protected leaf/neighbor/feature seals from SAME inverse and contact set')
        labels=c.E.labels(g,nodes);local=[nodes.index(m) for m in marks]
        direct=sha256();mapped=0;direct_count=0;requests=[];seen=set()
        for ideal,kinds in ((False,('A','B','D')),(True,('A','B','C','D'))):
            if ideal and p['forced_finite']:continue
            for order in range(len(nodes),18):
                for kind in kinds:
                    target=c.W.classical_graph(kind,order,affine=ideal)
                    for mp in r.embeddings(labels,target,c.FUEL,c.E):
                        # Every first extra target neighbor is genuinely a
                        # new supporting root: other displayed components
                        # are orthogonal to the base component.
                        if any(len(target[mp[nodes.index(v)]])!=len(adj[v]) for v in forbidden):continue
                        pair=tuple(sorted(mp[i] for i in local))
                        need(c.W._distances(target,pair[0])[pair[1]]==9 and
                             any(len(target[mp[nodes.index(v)]])==1 for v in marked),
                             'sealed marked leaf preserves actual marked distance9')
                        key=dict(high_index=p['index'],core_key=p['core_key'],base_nodes=nodes,
                            selected_marks=marks,target_kind=kind,target_order=order,
                            affine=ideal,induced_map=mp,marked_image=list(pair))
                        encoded=c.canonical(key);need(encoded not in seen,'unique complete protected target map')
                        seen.add(encoded);mapped+=1
                        value=c.W.affine_weight(kind,order,*pair) if ideal else c.W.finite_same_weight((kind,order),*pair)
                        if not ideal and (kind,order)==('D',13) and pair in ((0,11),(0,12)):
                            value=min(value,F(1,2))
                        if value==1:
                            direct_count+=1;direct.update(encoded.encode()+b'\n')
                        else:
                            need(ideal and (kind,order,*pair) in PRIORITIES,
                                 'unsupported CURRENT protected output: '+c.canonical([key,str(value)]))
                            requests.append(dict(key=key,intrinsic_current=str(value),
                                required_current_whole_face_provider=PRIORITIES[kind,order,*pair],
                                exact_source_theorem=PRIORITY_SOURCES[PRIORITIES[kind,order,*pair]],
                                current_priority_verified=False,
                                source_application='Whole actual affine marked factor, all-native opposite contacts and arbitrary additional supporting forest; finite source receipts must still be joined.'))
        out.append(dict(consumer_key=dict(high_index=p['index'],method=p['method'],core_key=p['core_key'],selected_marks=marks),
            physical_context_sha256=sha256(c.canonical([core['gram'],nodes,marks]).encode()).hexdigest(),
            retained_path=pp,sealed_core_nodes=sorted(forbidden),protected_marked_leaves=marked,
            forced_finite=p['forced_finite'],persistent_feature_distance=near,
            complete_small_full_maps_after_seals=mapped,direct_unit_full_maps=direct_count,
            direct_unit_map_sha256=direct.hexdigest(),exact_small_current_priority_requests=requests,
            higher_order_interface='Complete typed family predicate checked separately for all18..58; no arbitrary spectator filling or finite-order<=17 assertion.',
            required_independent_S8='finite57.current.S8.whole_face.v1',independent_S8_verified=False))
    need(len(out)==sum(p['index']%SHARDS==shard for p in chosen),'complete named protected shard')
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=tuple(f'protected-{i}' for i in range(SHARDS)),required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args=parser.parse_args();own=str(Path(__file__).resolve().relative_to(ROOT));source=raw(own)
    a,r,c,extra,priority=prepare();large=large_interface(c);rows=protected(int(args.batch[-1]),a,r,c)
    for p in c.READ:need(c.digest(p)==c.PINS[p],'consumed protected input changed')
    need(sha256(raw(LOADER)).hexdigest()==LOADER_SHA,
         'shared protected loader changed')
    need(sha256(raw(a.ROUTE)).hexdigest()==a.ROUTE_SHA and
         sha256(raw(r.COLLECTOR)).hexdigest()==r.COLLECTOR_SHA,
         'route/immutable collector changed')
    need(raw(a.EXTRA)==extra and raw(PRIORITY_INPUTS)==priority and raw(own)==source and
         sha256(c.raw(c.MANIFEST)).hexdigest()==c.MANIFEST_SHA,'protected controller/manifest changed')
    payload=json.dumps(dict(status='CURRENT_PROTECTED74_APPLICATION_INDEX_WITH_EXPLICIT_OPEN_PRIORITIES',
        batch=args.batch,source_sha256=sha256(source).hexdigest(),large_current_interface=large,
        current_prefix_applications=rows,consumed_inputs={p:c.PINS[p] for p in sorted(c.READ)},
        independent_S8_verified=False,current_terminal_bridge_complete=False,global_finite57_verified=False,
        universal_current_priority_receipts_verified=False,fuel_used=c.FUEL[0],
        premise='Accepted full-core allowed-contact/end-block inverse seal and current actual-factor geometry. Saved inverse validity/all-native domain completeness are not regenerated.'),sort_keys=True)
    need(len(payload.encode())+1<=MAX_OUTPUT,'complete protected output exceeds cap; refuse rather than truncate maps')
    print(payload)


if __name__=='__main__':main()
