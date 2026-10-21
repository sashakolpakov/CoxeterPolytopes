#!/usr/bin/env python3
"""Exact new-vertex and affine census for the positive-plane bridge pair."""
import argparse
import json
from itertools import combinations
from scripts import fresh_h9_remote_next_20261013_bridge_verify as B
from scripts import fresh_h9_remote_next_20261013_frontier as F
from scripts.fresh_h9_remote_caps_20261013_frontier import psd
V=B.V;K=V.K
OUT=B.PREFIX+'_frontier.json'

def pair_table():
    names,pairs=F.pair_table()
    for i in range(2):
        n='bridge'+str(i)
        for old in names:
            m=2 if old in B.N.S_NAMES+['wall22','nextH4','bridge0'] else 3 if i==0 and old=='wall38' else None
            pairs['|'.join(sorted((n,old))) ]={'kind':'dotted'}if m is None else{'kind':'finite','label':m}
        names.append(n)
    assert len(names)==88 and len(pairs)==3828
    return names,pairs

def run():
    data=V.read(B.PREFIX+'_caps.json');B.check(data)
    previous=V.read(B.N.PREFIX+'_frontier.json');assert previous==F.run()
    assert previous['new_affine_circuits']==previous['inherited_affine_circuits']==0
    names,pairs=pair_table()
    def pair(a,b):return pairs['|'.join(sorted((a,b)))]
    def kind(ss):
        if any(pair(a,b)['kind']=='dotted'for a,b in combinations(ss,2)):return'dotted_pair'
        vals={2:K.Z,3:K.neg(K.O),5:K.neg(K.P)}
        M=[[K.scale(K.O,2)if a==b else vals[pair(a,b)['label']]for b in ss]for a in ss]
        return psd(M)
    A=B.N.S_NAMES+['wall22','bridge0','bridge1'];C=B.N.S_NAMES+['nextH4','bridge0','bridge1']
    small=[['bridge0','wall38',n]for n in B.N.S_NAMES]
    maximal=[A,C]+small
    for row in maximal:assert kind(row)=='positive_definite'
    examined=0
    for new in('bridge0','bridge1'):
        neighbors=[n for n in names if n==new or pair(n,new)['kind']=='finite']
        for mask in range(1<<len(neighbors)):
            ss=[n for j,n in enumerate(neighbors)if mask>>j&1]
            if new not in ss:continue
            examined+=1;result=kind(ss)
            if result in('positive_definite','positive_semidefinite'):
                assert result=='positive_definite'and any(set(ss)<=set(m)for m in maximal)
    edges=[]
    for vertex in(A,C):
        for removed in vertex:
            block=[n for n in vertex if n!=removed];assert kind(block)=='positive_definite'
            tests=[];ends=[]
            for n in names:
                if n in vertex:continue
                result=kind(block+[n]);tests.append({'ninth':n,'result':result})
                if result in('positive_definite','positive_semidefinite'):ends.append((n,result))
            expected=['nextH4']if removed=='wall22'else['wall22']if removed=='nextH4'else[]
            assert [n for n,t in ends]==expected and all(t=='positive_definite'for n,t in ends)
            edges.append({'vertex_endpoint':vertex[-3],'removed':removed,'block':block,
                          'status':'finite_endpoint'if ends else'open_to_light_cone','other_endpoint':ends[0][0]if ends else None,'ninth_wall_tests':tests})
    return {'status':'PASS_EXACT_POSITIVE_PLANE_BRIDGE_FRONTIER','facets':88,'all_pair_count':3828,
            'source_sha256':{p:V.sha(p)for p in(B.PREFIX+'_caps.json',B.N.PREFIX+'_frontier.json')},
            'inherited_affine_circuits':0,'new_affine_circuits':0,'new_psd_subsets_checked':examined,
            'maximal_positive_definite_blocks_containing_new_facets':maximal,
            'new_finite_vertices':2,'new_vertex_types':['9A1','9A1'],'new_component_finite_edges':1,'new_component_open_rays':16,
            'edges':edges,'scope':'Complete positive-semidefinite-subset check for sets containing either of the two new walls. The two9A1 blocks give exactly two new finite vertices; the six H3 blocks have rank3. No new affine circuits occur. All18 outgoing rays of the new vertices are checked; one finite edge joins them and sixteen rays remain open. Existing finite vertices persist by the exact acute/common-center projection argument.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert V.read(OUT)==out
    else:(V.ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_finite_vertices'],'new vertices',out['new_component_open_rays'],'open rays')
