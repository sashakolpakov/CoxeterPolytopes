#!/usr/bin/env python3
"""Exact new H4 vertex frontier and affine gate, including the pi/10 pair."""
import argparse
from itertools import combinations
import json
from scripts import fresh_h9_remote_next_20261013_verify as N
from scripts.fresh_h9_remote_caps_20261013_frontier import psd
V=N.V;W=N.W;K=V.K
OUT=N.PREFIX+'_frontier.json'

def pair_table():
    source=V.read(V.OLD);old=V.read(source['source']);pairs=dict(old['pair_checks']);names=[r['name']for r in source['roots']]
    first=V.read(V.PREFIX+'_caps.json')
    def put(a,b,m):pairs['|'.join(sorted((a,b)))]={'kind':'dotted'}if m is None else{'kind':'finite','label':m}
    for i in range(4):
        n='remote'+str(i)
        for b in names:put(n,b,2 if b in first['orthogonal_old_walls']or b.startswith('remote')else None)
        names.append(n)
    for i in range(2):
        n='coupled'+str(i)
        for b in names:put(n,b,2 if b in W.S_NAMES or b=='coupled0'else 3 if i==0 and b=='wall17'else None)
        names.append(n)
    for b in names:put('nextH4',b,2 if b in N.S_NAMES+['wall17']else 3 if b=='coupled0'else 10 if b=='coupled1'else None)
    names.append('nextH4');assert len(names)==86 and len(pairs)==3655
    return names,pairs

def run():
    data=V.read(N.PREFIX+'_caps.json');N.check(data)
    earlier=V.read(V.PREFIX+'_frontier.json');assert earlier['new_affine_circuits']==earlier['inherited_affine_circuits']==0
    names,pairs=pair_table()
    def pair(a,b):return pairs['|'.join(sorted((a,b)))]
    def kind(ss):
        if any(pair(a,b)['kind']=='dotted'for a,b in combinations(ss,2)):return'dotted_pair'
        if {'nextH4','coupled0','coupled1'}<=set(ss):
            # Diagonal2, offdiag0,-1,-L, L²=2+phi: determinant2-2phi<0.
            assert K.sign(K.sub(K.scale(K.O,2),K.scale(K.P,2)))<0
            return'indefinite_pi3_pi10_triple'
        if {'nextH4','coupled1'}<=set(ss):
            rest=[n for n in ss if n not in('nextH4','coupled1')]
            assert all(pair(n,a)=={'kind':'finite','label':2}for n in rest for a in('nextH4','coupled1'))
            assert K.sign(K.sub(K.scale(K.O,2),K.P))>0 # det I2(10)=4-L².
            return kind(rest)if rest else'positive_definite'
        labels={2:K.Z,3:K.neg(K.O),5:K.neg(K.P)}
        M=[[K.scale(K.O,2)if a==b else labels[pair(a,b)['label']]for b in ss]for a in ss]
        return psd(M)
    A=N.S_NAMES+['wall17','coupled0','nextH4'];B=N.S_NAMES+['coupled1','nextH4']
    assert kind(A)==kind(B)=='positive_definite'
    neighbors=[n for n in names if n=='nextH4'or pair(n,'nextH4')['kind']=='finite'];assert len(neighbors)==10
    for mask in range(1<<len(neighbors)):
        ss=[n for i,n in enumerate(neighbors)if mask>>i&1]
        if'nextH4'not in ss:continue
        result=kind(ss)
        if result in('dotted_pair','indefinite_pi3_pi10_triple'):continue
        assert result=='positive_definite'and(set(ss)<=set(A)or set(ss)<=set(B))
    edges=[]
    for removed in A:
        block=[n for n in A if n!=removed];assert kind(block)=='positive_definite'
        tests=[];ends=[]
        for n in names:
            if n in A:continue
            result=kind(block+[n]);tests.append({'ninth':n,'result':result})
            if result in('positive_definite','positive_semidefinite'):ends.append((n,result))
        assert len(ends)<=1 and all(t=='positive_definite'for n,t in ends)
        edges.append({'removed':removed,'block':block,'status':'finite_endpoint'if ends else'open_to_light_cone',
                      'other_endpoint':ends[0][0]if ends else None,'ninth_wall_tests':tests})
    assert [(r['removed'],r['other_endpoint'])for r in edges if r['other_endpoint']]==[('nextH4','seed0')]
    isolated=[]
    for n in names:
        if n in B:continue
        result=kind(B+[n]);assert result not in('positive_definite','positive_semidefinite');isolated.append({'ninth':n,'result':result})
    # Before this insertion, dropping seed0 at the H3 vertex was open.
    previous_edge=next(r for r in earlier['edges']if r['vertex']=='attachedH3'and r['removed_facet']=='seed0')
    assert previous_edge['status']=='open_to_light_cone'and set(previous_edge['edge_walls'])==set(A)-{'nextH4'}
    return {'status':'PASS_EXACT_H4_ENDPOINT_FRONTIER','source_sha256':{p:V.sha(p)for p in(N.PREFIX+'_caps.json',V.PREFIX+'_frontier.json',V.OLD)},
            'facets':86,'all_pair_count':3655,'new_affine_circuits':0,'inherited_affine_circuits':0,
            'maximal_positive_definite_blocks_containing_new_facet':[A,B],
            'indefinite_triangle':{'walls':['nextH4','coupled0','coupled1'],'determinant':K.enc(K.sub(K.scale(K.O,2),K.scale(K.P,2)))},
            'new_finite_vertex_count':1,'new_vertex_type':'H4+5A1','new_vertex_finite_edges':1,'new_vertex_open_rays':8,
            'edges':edges,'pi10_spherical_line':{'walls':B,'rank':8,'ninth_wall_tests':isolated,'finite_or_ideal_endpoints':0},
            'scope':'Every positive-semidefinite subset containing the new wall lies in one of the two stated positive-definite blocks. Thus the cap creates exactly one finite vertex and no affine circuits. All9 outgoing rays of that vertex are checked, and the additional rank8 I2(10)+6A1 face is an actual complete H1 line with no finite or ideal wall endpoint in this86-wall chamber. Further compatible walls may change this frontier; no universal completion obstruction is claimed.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert V.read(OUT)==out
    else:(V.ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_vertex_finite_edges'],'finite',out['new_vertex_open_rays'],'open')
