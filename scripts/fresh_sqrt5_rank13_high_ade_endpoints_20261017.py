#!/usr/bin/env python3
"""Complete one-wall finite transitions across the rank13 high-baseline ADE set."""
import argparse,hashlib,json
from collections import Counter,deque
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank13_high_ade_endpoints_20261017'
DIMENSION=13;THRESHOLD=F(12)

def family(g):
    n=len(g)
    if n==1:return 'A1'
    degrees=[len(g[i])for i in g]
    if sum(degrees)!=2*(n-1)or max(degrees)>3:return None
    branch=[i for i in g if len(g[i])==3]
    if not branch:return 'A'+str(n)
    if len(branch)!=1:return None
    center=branch[0];arms=[]
    for v in g[center]:
        previous=center;length=1
        while len(g[v])==2:
            nxt=next(x for x in g[v]if x!=previous);previous,v=v,nxt;length+=1
        arms.append(length)
    arms.sort()
    if arms[:2]==[1,1]:return 'D'+str(n)
    return {(1,2,2):'E6',(1,2,3):'E7',(1,2,4):'E8'}.get(tuple(arms))

def classify(g):
    unseen=set(g);names=[]
    while unseen:
        first=min(unseen);seen={first};stack=[first]
        for i in stack:
            for j in g[i]:
                if j not in seen:seen.add(j);stack.append(j)
        unseen-=seen;name=family({i:g[i]&seen for i in seen})
        if name is None:return None
        names.append(name)
    return tuple(sorted(names))

def component(name,weighted):
    g={i:set(x)for i,x in weighted.items()};assert family(g)==name
    pairs=[]
    for i,j in combinations(g,2):
        distances={i:0};queue=deque([i])
        while queue:
            a=queue.popleft()
            for b in g[a]:
                if b not in distances:distances[b]=distances[a]+1;queue.append(b)
        weight=int(distances[j]<=2);assert weight==C.base(weighted,i,j)
        pairs.append(dict(pair=[i,j],distance=distances[j],baseline=weight))
    return dict(name=name,rank=len(g),edges=[list(p)for p in combinations(g,2)if p[1]in g[p[0]]],
        pairs=pairs,baseline=sum(x['baseline']for x in pairs)),g

def run():
    models=[(n,g)for n,g in C.models()if n[0]in'ADE'and len(g)<=DIMENSION]
    comp=[];graphs={}
    for name,weighted in models:
        c,g=component(name,weighted);comp.append(c);graphs[name]=g
    def enumerate_forests(rank,start=0):
        if not rank:yield ();return
        for j in range(start,len(comp)):
            if comp[j]['rank']<=rank:
                for tail in enumerate_forests(rank-comp[j]['rank'],j):yield(j,)+tail
    census=[]
    for inds in enumerate_forests(DIMENSION):
        names=tuple(sorted(comp[j]['name']for j in inds));base=sum(comp[j]['baseline']for j in inds)
        census.append(dict(forest=list(names),baseline=base,margin=str(F(base)-THRESHOLD),selected=base>=THRESHOLD))
    census.sort(key=lambda x:(-x['baseline'],x['forest']))
    selected=[x for x in census if x['selected']];lookup={tuple(x['forest']):i for i,x in enumerate(selected)}
    assert len(census)==276 and len(selected)==215
    targets=[];all_edges=set();counts=Counter()
    for idx,rec in enumerate(selected):
        g={};components=[]
        for name in rec['forest']:
            h=graphs[name];offset=len(g)
            for i in h:g[offset+i]={offset+j for j in h[i]}
            components.append(dict(type=name,nodes=list(range(offset,offset+len(h)))))
        assert len(g)==DIMENSION and classify(g)==tuple(rec['forest'])
        ports=[];destinations=set()
        for port in range(DIMENSION):
            remaining=[i for i in g if i!=port];old={i:g[i]-{port}for i in remaining};attachments=[]
            for size in range(4):
                for contacts in combinations(remaining,size):
                    h={i:set(x)for i,x in old.items()};h[DIMENSION]=set(contacts)
                    for i in contacts:h[i].add(DIMENSION)
                    result=classify(h);counts['tested']+=1
                    if result is None:counts['nonspherical']+=1;continue
                    destination=lookup.get(result)
                    if destination is None:counts['outside_high_set']+=1;continue
                    attachments.append(dict(contacts=list(contacts),destination=destination,
                        copied_contacts=set(contacts)==g[port]))
                    destinations.add(destination);all_edges.add(tuple(sorted((idx,destination))));counts['accepted']+=1
            assert any(a['destination']==idx and a['copied_contacts']for a in attachments)
            ports.append(dict(port=port,attachments=attachments))
        targets.append(dict(index=idx,**rec,components=components,
            edges=[list(p)for p in combinations(g,2)if p[1]in g[p[0]]],ports=ports,
            endpoint_attachments=sum(len(x['attachments'])for x in ports),destination_types=sorted(destinations)))
        print('ADE13_ENDPOINTS',idx,len(selected),targets[-1]['endpoint_attachments'],len(destinations),flush=True)
    # Connected components of the type-level graph, counting self-edges only
    # as local compatibility and never as a simultaneous realization claim.
    adjacency={i:set()for i in range(len(targets))}
    for i,j in all_edges:adjacency[i].add(j);adjacency[j].add(i)
    unseen=set(adjacency);type_components=[]
    while unseen:
        seen={min(unseen)};queue=list(seen)
        for i in queue:
            for j in adjacency[i]:
                if j not in seen:seen.add(j);queue.append(j)
        unseen-=seen;type_components.append(sorted(seen))
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
        'results/fresh_sqrt5_discount_structural_floor_20261016.json']
    return dict(status='PASS_COMPLETE_RANK13_HIGH_BASELINE_ADE_FINITE_ENDPOINTS',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dimension=DIMENSION,
        strict_h_threshold=str(THRESHOLD),components=comp,forest_census=census,selected_types=len(selected),
        all_ADE_types=len(census),counts=dict(counts),targets=targets,type_adjacency_edges=[list(p)for p in sorted(all_edges)],
        type_adjacency_components=type_components,
        completeness='All spherical simply-laced components are A,D,E trees. The finite destination root has degree0,1,2,or3 and all its nonright contacts have label3. Every such contact subset of all12 retained old roots is enumerated. The full resulting forest, not just the component at the new root, must belong to the selected set.',
        all_special_triangle_exclusion='If all three triangle vertices belong to this ADE set, every entry of its full normalized star Gram is rational. Its physical negative direction contradicts arithmetic conjugate positivity over Q(sqrt5), independently of whether the three forests agree.',
        scope='Complete finite one-wall transitions and baseline margins for a selected SET of215 rank13 ADE forests, not a framing-parameter classification. All non-ADE vertex types and61 lower-baseline ADE forests are outside this set. No compact realization, short-face coverage, joint incidence inequality or dimension bound is claimed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['selected_types'],data['counts'],list(map(len,data['type_adjacency_components'])),flush=True)
