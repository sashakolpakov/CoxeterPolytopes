#!/usr/bin/env python3
"""Raw357 conservative capacities and complete ordinary-anchor transport.

A structural pilot, not a new column export or global optimization.
"""
import argparse,json,time
from collections import Counter,deque
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_raw357_relative_marker_pilot_20261018 as P
G=P.G;L=P.L;U=P.U;H=P.H;R=P.R;T=P.T;D=L.D;ROOT=P.ROOT
TAG='fresh_raw357_capacity_ordinary_pilot_20261018';normal=P.normal


def generators(target):
    graph,parts=U.graph_components(target);out=[];seen=set();identity=dict(zip(range(13),range(13)))
    def add(p):
        assert set(p)==set(p.values())==set(range(13))
        assert all(target[i].get(j,2)==target[p[i]].get(p[j],2)for i,j in combinations(range(13),2))
        key=tuple(p[i]for i in range(13))
        if key!=tuple(range(13))and key not in seen:seen.add(key);out.append(p)
    groups=[]
    for part in parts:
        g=graph.subgraph(part)
        for mp in nx.algorithms.isomorphism.GraphMatcher(g,g,edge_match=lambda a,b:a['label']==b['label']).isomorphisms_iter():add(identity|mp)
        for group in groups:
            iso=nx.algorithms.isomorphism.GraphMatcher(graph.subgraph(group[-1]),g,edge_match=lambda a,b:a['label']==b['label'])
            if iso.is_isomorphic():
                mp=next(iso.isomorphisms_iter());add(identity|mp|{b:a for a,b in mp.items()});group.append(part);break
        else:groups.append([part])
    return out


def role_blocks(target,m):
    cs=P.components(target,m)
    # Components are ordered by the original raw357 node order.
    k=m['kind']
    assert [len(c)for c in cs]==([2,4,2]if k==3574 else[2,3,1,2]if k==3575 else[2,2,2,2])
    return cs


def cycle_key(u,v):
    return min((u,v),(v,tuple(reversed(u))),(tuple(reversed(u)),tuple(reversed(v))),(tuple(reversed(v)),u))


def rule(target,a,b):
    ka,kb=a['kind'],b['kind']
    if ka>kb:return rule(target,b,a)
    aa=role_blocks(target,a);bb=role_blocks(target,b)
    paira=(aa[0],aa[-1]);pairb=(bb[0],bb[-1]);same=cycle_key(*paira)==cycle_key(*pairb)
    disjoint=not(set(sum(paira,()))&set(sum(pairb,())))
    if ka==kb==3574:return same or(disjoint and aa[1]==bb[1])
    if ka==kb==3575:return aa[1:3]==bb[1:3]and(same or disjoint)
    if ka==kb==3576:return aa[1:3]==bb[1:3]and(same or disjoint)
    if(ka,kb)==(3574,3575):return same
    return False


def capacity(target,forest,family,cases):
    co=Counter(forest);b,h=co['A2'],co['H4'];out=[]
    for f in family:
        kind=f['kind'];markers=f['markers'];cap=4*max(h,b//2)if kind==3574 else 4*(b//2)if kind==3575 else 4*((b-1)//2)
        chosen=[]
        for marker in markers:
            if all(rule(target,marker,other)for other in chosen):chosen.append(marker)
        assert len(chosen)==cap
        for a,z in combinations(chosen,2):
            key,pat,payload=P.signature(target,kind,a,kind,z)
            assert key in cases and cases[key]['result']['allowed']and cases[key]['payload']==normal(payload)
        out.append(dict(kind=kind,capacity=cap,witness=[m['index']for m in chosen]))
    # The only mixed new-kind branch is A/B. A clique containing both uses
    # one common oriented A2 pair; B also fixes its A1 and H3 components.
    caps={r['kind']:r['capacity']for r in out};joint=max(caps.values())
    if 3574 in caps and 3575 in caps:
        assert b==2 and co['H4']==co['H3']==1
        joint=max(joint,caps[3574]+caps[3575])
    assert joint<=8
    return dict(individual=out,joint_new_marker_upper=joint,
        interpretation='Conservative pair-graph capacity only. Different ordinary/old marker coefficients still require one shared allowance; no charge is exported.')


def run():
    start=time.monotonic();data=json.loads((ROOT/f'results/{P.TAG}.json').read_text())
    assert data['status']=='PASS_RAW357_COMPLETE_RELATIVE_MARKER_PILOT'
    for p,h in data['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    cases={H.frozen([c['kinds'][0],c['kinds'][1],c['pattern']]):c for c in data['relative_cases']}
    assert len(cases)==578
    for c in data['relative_cases']:
        if c['kinds'][1]in P.KINDS:
            target={int(i):{int(j):m for j,m in r.items()}for i,r in c['finite_frame'].items()}
            assert rule(target,c['anchor'],c['other'])==c['result']['allowed']
    inp=D.Q.Inputs();domains,_=D.L.current_domains();prepared=D.Z.prepare()
    local={r['vertex']:r for r in inp.read(D.L.TAG)['forest_records']}
    patches={r['vertex']:r for r in inp.read(D.Z.TAG)['forest_records']}
    records=[];anchor_tests=unknowns=transported_pairs=0
    for record in data['forest_records']:
        v=record['vertex'];forest=record['forest'];target={int(i):{int(j):m for j,m in r.items()}for i,r in record['finite_frame'].items()}
        row={'target_forest':forest};nodes,groups,graph,owner,got,history=U.complete_current(v,row,prepared,domains,local,patches)
        assert got==forest
        canonical=[U.transport_native(s,dict(zip(range(13),range(13))))for s in nodes]
        keys={H.frozen(s['contacts']):i for i,s in enumerate(canonical)};assert len(keys)==len(nodes)
        gg=generators(target);saved_generators=[]
        for perm in gg:
            action=[]
            for i,s in enumerate(nodes):
                moved=U.transport_native(s,perm);j=keys[H.frozen(moved['contacts'])]
                assert all(normal(moved[k])==normal(canonical[j][k])for k in H.NATIVE_FIELDS)and owner[i]==owner[j]
                action.append(j)
            assert len(set(action))==len(nodes)
            saved_generators.append(dict(frame_permutation=sorted(perm.items()),ordinary_permutation=action))
        J,ac,blocks,tens=T.projection_data(target);families=[]
        for f in record['marker_families']:
            markers=f['markers'];mk={H.frozen(m['native']['contacts']):m['index']for m in markers};assert len(mk)==len(markers)
            acts=[]
            for perm in gg:
                action=[]
                for m in markers:
                    moved=U.transport_native(m['native'],perm);j=mk[H.frozen(moved['contacts'])]
                    assert all(normal(moved[k])==normal(markers[j]['native'][k])for k in H.NATIVE_FIELDS)
                    ports=sorted(sorted(perm[i]for i in p)for p in m['triangle_port_pairs'])
                    assert ports==markers[j]['triangle_port_pairs'];action.append(j)
                assert len(set(action))==len(markers);acts.append(action)
            reached={0};todo=deque([0])
            while todo:
                i=todo.popleft()
                for action in acts:
                    j=action[i]
                    if j not in reached:reached.add(j);todo.append(j)
            assert reached==set(range(len(markers)))
            anchor=markers[0];decisions=[]
            for i,s in enumerate(nodes):
                ports=[p for p in anchor['triangle_port_pairs']if p in s['occupied_pairs']]
                result=dict(allowed=False,triangle_port_exclusion=ports,expressions=[])if ports else L.pair_test(J,tens,s,anchor['native'])
                decisions.append(dict(ordinary=i,result=result))
            count=sum(r['result']['allowed']for r in decisions)
            anchor_tests+=len(nodes);unknowns+=count;transported_pairs+=count*len(markers)
            families.append(dict(kind=f['kind'],source_face=f['source_face'],actual_extra_occurrences=1,
                marker_count=len(markers),generator_marker_permutations=acts,anchor_decisions=decisions,
                anchor_unknowns=count,retained_ordinary_pairs=count*len(markers)))
        records.append(dict(vertex=v,forest=forest,finite_frame=target,ordinary_nodes=nodes,ordinary_groups=groups,
            old_graph_edges=sorted(tuple(sorted(e))for e in graph.edges()),old_owner=owner,prior_history=history,
            generators=saved_generators,marker_families=families,
            capacity=capacity(target,forest,record['marker_families'],cases)))
        if len(records)%20==0 or any(f['anchor_unknowns']for f in families):print('FOREST',len(records),v,'anchors',anchor_tests,'unknowns',unknowns,'seconds',round(time.monotonic()-start,1),flush=True)
    assert len(records)==138 and anchor_tests==7282
    paths=set(inp.paths)|set(prepared['sources'])|set(data['sources'])
    paths.update({str(Path(__file__).relative_to(ROOT)),f'scripts/{P.TAG}.py',f'results/{P.TAG}.json'})
    for mod in(D,D.L,D.Z,U,H):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    return normal(dict(status='PASS_RAW357_COMPLETE_CAPACITY_ORDINARY_PILOT',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
        forest_records=records,counts=dict(source_classes=47,pointed_domains=141,target_forests=138,role_orientations=10624,
            relative_cases=578,retained_relative_cases=27,ordinary_anchor_tests=anchor_tests,
            retained_anchor_unknowns=unknowns,retained_ordinary_pairs=transported_pairs),
        old_outcomes_changed=False,source_quad_demand=1,extra_occurrences_per_source=3,
        exported_columns=[],optimizer_invocations=0,
        scope='Complete conservative raw357 marker capacities and ordinary compatibility by native/owner-preserving generators. Every old Q-Q outcome and all retained unknowns remain unchanged. No joint old/new one-budget charge or global optimization is asserted.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists(),path
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
