#!/usr/bin/env python3
"""Jointly test ten same-quartic, no-removal cutters of the selected gap."""
import json
import time
from itertools import combinations
import mpmath as mp
from scripts.fresh_h9_split_endpoint_walk_exact_20261013 import (
    source,read,OPTIONAL_SCOUT,SCOUT,PREFIX,ROOT,BASIS,enc,parse,sha,angle_pool,
    ANGLE_SOURCE,target_proposal,make,relation,finite_vertex,affine_census,volume_gate)

CHOSEN=(2720,2721,2722,2723,2735,2747,6746,6747,6748,6758)

def cliques(graph):
    result=[]
    def bk(r,p,x):
        if not p and not x:result.append(sorted(r));return
        u=max(p|x,key=lambda j:len(p&graph[j]))if p or x else None
        for j in sorted(p-(graph[u]if u is not None else set())):
            bk(r|{j},p&graph[j],x&graph[j]);p.remove(j);x.add(j)
    bk(set(),set(range(len(graph))),set())
    return sorted(result,key=lambda x:(-len(x),x))

def run():
    mp.mp.dps=90;start=time.monotonic();G,I,h,old,scout,banks=source();optional=read(OPTIONAL_SCOUT)
    assert optional['source_snapshot']==scout['source_snapshot']
    proof=target_proposal(scout,I,h,old);target=tuple(map(parse,proof['y']));pool=angle_pool(ANGLE_SOURCE)
    caps=[];new=[]
    archive={'status':'EXACT_TEN_NO_REMOVAL_EDGE_CAPS','source':SCOUT,'proposal_source':OPTIONAL_SCOUT,
             'source_sha256':{p:sha(p)for p in [SCOUT,OPTIONAL_SCOUT,*banks]},
             'independent_seed_basis':list(BASIS),'gram':[[enc(x)for x in r]for r in G],
             'inverse':[[enc(x)for x in r]for r in I],'strict_center':list(map(enc,h)),
             'old_walls':list(old),'target':proof,'caps':caps,
             'scope':'Ten no-removal cutters in the same-quartic block/anchor subfamily. Actual jointly retained configurations are certified separately; no exhaustion of all edge caps is claimed.'}
    for index in CHOSEN:
        row,record=make(optional['candidates'][index],index,I,h,old,target,pool)
        new.append(row);caps.append(record)
        (ROOT/(PREFIX+'_family.json')).write_text(json.dumps(archive,indent=2)+'\n')
        print({'cap':index,'seconds':round(time.monotonic()-start,1)},flush=True)
    relations={};graph=[set()for _ in new]
    for i,j in combinations(range(len(new)),2):
        rr=relation(new[i],new[j],pool);relations[f'{i}|{j}']=rr
        if rr['kind']in('finite','dotted'):graph[i].add(j);graph[j].add(i)
    maximal=cliques(graph)
    print({'maximal_selections':maximal},flush=True)
    pairs={}
    for a,b in combinations(old.values(),2):
        rr=relation(a,b,pool);assert rr['kind']in('finite','dotted')
        pairs['|'.join(sorted((a['name'],b['name'])))]=rr
    cases=[]
    result={'source':PREFIX+'_family.json','new_pair_relations':relations,
            'maximal_accepted_cliques':maximal,'cases':cases,
            'scope':'Every maximal selection in the archived accepted compatibility graph is tested. Exact escapes prove only these literal configurations infinite volume; smaller subcollections are their enlargements.'}
    for ids in maximal:
        selected=old|{new[i]['name']:new[i]for i in ids};checks=dict(pairs)
        for i in ids:checks.update({'|'.join(sorted((new[i]['name'],n))):rr for n,rr in caps[i]['pair_checks'].items()})
        for i,j in combinations(ids,2):checks['|'.join(sorted((new[i]['name'],new[j]['name'])))]=relations[f'{min(i,j)}|{max(i,j)}']
        census=affine_census(list(selected),checks)
        vertices=[finite_vertex(new[i],caps[i],I,h,selected)for i in ids]
        geo=volume_gate(I,h,list(selected.values()))
        cases.append({'new_caps':ids,'walls':list(selected),'facet_count':len(selected),
                      'pair_checks':checks,'affine_census':census,'actual_finite_vertices':vertices,
                      'geometry':geo,'status':geo['status']})
        (ROOT/(PREFIX+'_family_geometry.json')).write_text(json.dumps(result,indent=2)+'\n')
        print({'selection':ids,'facets':len(selected),'affine':len(census['affine_circuits']),
               'status':geo['status'],'max_norm':geo['maximum_numeric_vertex_norm'],
               'seconds':round(time.monotonic()-start,1)},flush=True)

if __name__=='__main__':run()
