#!/usr/bin/env python3
"""Exact local endpoint and affine gates for the three new remote vertices."""
import argparse
from itertools import combinations
import json
from scripts import fresh_h9_remote_caps_20261013_verify as V
from scripts import fresh_h9_remote_caps_20261013_coupled_verify as W
K=V.K
OUT=V.PREFIX+'_frontier.json'

def psd(matrix):
    m=[r[:]for r in matrix];rank=0
    for j in range(len(m)):
        if any(K.sign(m[i][i])<0 for i in range(j,len(m))):return 'indefinite'
        pivot=next((i for i in range(j,len(m))if K.sign(m[i][i])>0),None)
        if pivot is None:
            if any(m[a][b]!=K.Z for a in range(j,len(m))for b in range(j,len(m))):return'indefinite'
            return'positive_semidefinite'
        if pivot!=j:
            m[j],m[pivot]=m[pivot],m[j]
            for row in m:row[j],row[pivot]=row[pivot],row[j]
        q=m[j][j];rank+=1
        for a in range(j+1,len(m)):
            for b in range(j+1,len(m)):m[a][b]=K.sub(m[a][b],K.div(K.mul(m[a][j],m[j][b]),q))
    return'positive_definite'

def run():
    d0=V.read(V.PREFIX+'_caps.json');d1=V.read(W.PREFIX+'_caps.json');V.check(d0);W.check(d1)
    old=V.read(V.OLD);geometry=V.read(old['source']);names=[r['name']for r in old['roots']]
    pairs=dict(geometry['pair_checks']);assert not geometry['affine_census']['affine_circuits']
    def put(a,b,kind,label=None):
        pairs['|'.join(sorted((a,b)))]={'kind':kind}|({'label':label}if label is not None else{})
    for i in range(4):
        n='remote'+str(i)
        for oldn in names:put(n,oldn,'finite',2)if oldn in d0['orthogonal_old_walls']or oldn.startswith('remote')else put(n,oldn,'dotted')
        names.append(n)
    for i in range(2):
        n='coupled'+str(i)
        for oldn in names:
            if oldn in W.S_NAMES or oldn=='coupled0':put(n,oldn,'finite',2)
            elif i==0 and oldn=='wall17':put(n,oldn,'finite',3)
            else:put(n,oldn,'dotted')
        names.append(n)
    assert len(names)==85 and len(pairs)==3570
    def pair(a,b):return {'kind':'finite','label':1}if a==b else pairs['|'.join(sorted((a,b)))]
    def kind(ss):
        if any(pair(a,b)['kind']=='dotted'for a,b in combinations(ss,2)):return'dotted_pair'
        labels={2:K.Z,3:K.neg(K.O),5:K.neg(K.P)}
        M=[[K.scale(K.O,2)if a==b else labels[pair(a,b)['label']]for b in ss]for a in ss]
        return psd(M)
    # Every nondotted subset containing a new facet is contained in one
    # of these positive-definite maximal blocks. Hence no new affine circuit.
    blocks=[d0['orthogonal_old_walls']+['remote'+str(i)for i in range(4)],
            W.S_NAMES+['coupled0','coupled1'],W.S_NAMES+['wall17','coupled0']]
    assert all(kind(ss)=='positive_definite'for ss in blocks)
    for name in ['remote'+str(i)for i in range(4)]+['coupled0','coupled1']:
        neighbors=[n for n in names if n==name or pair(n,name)['kind']=='finite']
        available=[set(ss)for ss in blocks if name in ss]
        assert len(neighbors)<=10
        for mask in range(1<<len(neighbors)):
            ss=[n for i,n in enumerate(neighbors)if mask>>i&1]
            if name not in ss or any(pair(a,b)['kind']=='dotted'for a,b in combinations(ss,2)):continue
            assert any(set(ss)<=b for b in available)
    vertices={'remote9A1':blocks[0],'attachedH3':blocks[2],'coupled9A1':blocks[1]}
    records=[]
    for vertex,incident in vertices.items():
        assert kind(incident)=='positive_definite'
        for drop in incident:
            block=[n for n in incident if n!=drop];assert kind(block)=='positive_definite'
            tests=[];ends=[]
            for n in names:
                if n in incident:continue
                result=kind(block+[n]);tests.append({'ninth':n,'result':result})
                if result=='positive_definite':ends.append(n)
                assert result!='positive_semidefinite','Unexpected ideal endpoint; inspect'
            assert len(ends)<=1
            records.append({'vertex':vertex,'removed_facet':drop,'edge_walls':block,
                            'status':'finite_endpoint'if ends else'open_to_light_cone',
                            'other_endpoint_facet':ends[0]if ends else None,'ninth_wall_tests':tests})
    counts={n:{'finite':sum(r['status']=='finite_endpoint'for r in records if r['vertex']==n),
               'open':sum(r['status']=='open_to_light_cone'for r in records if r['vertex']==n)}for n in vertices}
    assert counts=={'remote9A1':{'finite':0,'open':9},'attachedH3':{'finite':2,'open':7},'coupled9A1':{'finite':1,'open':8}}
    # The old spherical edge had no second finite or ideal endpoint before
    # the two new caps: the only PSD ninth old facet is seed4.
    old_edge=W.S_NAMES+['wall17'];old_tests=[{'ninth':n,'result':kind(old_edge+[n])}for n in names[:-2]if n not in old_edge]
    assert [r['ninth']for r in old_tests if r['result']in('positive_definite','positive_semidefinite')]==['seed4']
    return {'status':'PASS_EXACT_THREE_VERTEX_FRONTIER','source_sha256':{p:V.sha(p)for p in(V.PREFIX+'_caps.json',W.PREFIX+'_caps.json',V.OLD,old['source'])},
            'total_facets':85,'pair_count':3570,'new_affine_circuits':0,'inherited_affine_circuits':0,
            'maximal_nondotted_blocks_containing_new_facets':blocks,'selected_vertex_counts':counts,
            'directed_edges_checked':len(records),'ninth_wall_tests':sum(len(r['ninth_wall_tests'])for r in records),
            'old_open_edge':{'walls':old_edge,'old_ninth_wall_tests':old_tests,'new_endpoint':'coupled0'},
            'edges':records,'scope':'Complete endpoint tests only for the27 outgoing rays at the three specified new vertices, plus the old edge now closed. This is not a census of every finite vertex or component in the85-wall chamber. All85 remain infinite-volume by the independently replayed strict escape.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert V.read(OUT)==out
    else:(V.ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['selected_vertex_counts'])
