#!/usr/bin/env python3
"""Separate symbolic H4-prism Gram checks and nested 160-digit replay."""
import argparse,json,sys
from pathlib import Path
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_face_bridge_20261013_verify as C
V=C.V
OUT=V.ROOT/'results/fresh_h9_face_bridge_20261013_review.json'

def run():
    phi=(1+sp.sqrt(5))/2
    L=sp.Matrix([[2,-1,0,0],[-1,2,-phi,0],[0,-phi,2,-1],[0,0,-1,2]])
    p=sp.Matrix([0,0,0,-1]);qw=sp.simplify((p.T*L.inv()*p)[0])
    assert qw<0
    assert sp.simplify(qw-(4-6*phi)/11)==0
    assert sp.simplify(L.inv()[1,1]+(6+24*phi)/11)==0
    assert sp.simplify(L.inv()[3,1]+(8+10*phi)/11)==0
    assert sp.simplify(qw+2*(9+3*phi)/11)==2
    assert sp.simplify(L.inv()[1,1]+1+2*(17+24*phi)/22)==2
    assert sp.simplify(L.inv()[3,1])<-2
    cross=-sp.sqrt(4-2*qw);assert cross<-2
    G=sp.zeros(6)
    G[:4,:4]=L;G[4,4]=G[5,5]=2
    G[:4,5]=p;G[5,:4]=p.T;G[4,5]=G[5,4]=cross
    assert sp.simplify(G.det())==0
    assert sp.simplify(G[:5,:5].det())<0
    for omit in range(4):
        side=[i for i in range(4)if i!=omit]
        for base in(4,5):
            ids=side+[base];M=G.extract(ids,ids)
            assert all(sp.simplify(M[:k,:k].det())>0 for k in range(1,5))
    assert all(sp.simplify(x)<0 for x in L.inv())
    source=C.PREFIX+'_caps.json';saved=V.read(C.PREFIX+'_verify.json')
    V.SCALE=10**160
    V.PHI=V.div(V.add(V.rat(1),V.sqrt(V.rat(5))),V.rat(2));V.radical.cache_clear()
    replay=C.check(V.read(source));factor=10**80;count=0
    def compare(a,b,path='root'):
        nonlocal count
        if a==b:return
        if isinstance(a,dict):
            assert isinstance(b,dict)and a.keys()==b.keys(),path
            for k in a:compare(a[k],b[k],path+'.'+k)
        elif isinstance(a,list):
            assert isinstance(b,list)and len(a)==len(b),path
            if len(a)==2 and all(isinstance(x,str)for x in a):
                try:lo,hi=map(int,a);ll,hh=map(int,b)
                except ValueError:pass
                else:
                    assert lo*factor<=ll<=hh<=hi*factor,path
                    count+=1;return
            for i,(x,y)in enumerate(zip(a,b,strict=True)):compare(x,y,path+f'[{i}]')
        else:assert a==b,path
    compare(saved,replay)
    assert replay['prism_vertices']==8 and len(replay['actual_finite_vertices'])==6
    assert replay['vertex_component_frontier']['vertices']==10
    assert replay['vertex_component_frontier']['finite_edges']==19
    assert replay['vertex_component_frontier']['open_rays']==52
    # Correct the historical textual type using the actual exact Gram.
    common=replay['prism_common_facets']
    vertex=next(v for v in replay['actual_finite_vertices']if set(common)<=set(v['incident']))
    indices=[vertex['incident'].index(n)for n in common]
    M=[[C.K.parse(vertex['gram'][i][j])for j in indices]for i in indices]
    assert M==[[C.K.scale(C.K.O,2)if i==j else C.K.Z for j in range(5)]for i in range(5)]
    return {'status':'PASS_SYMBOLIC_PRISM_AND_160_DIGIT_REPLAY',
        'source':source,'source_sha256':V.sha(source),
        'verifier_sha256':V.sha('scripts/fresh_h9_face_bridge_20261013_verify.py'),
        'facets':96,'new_actual_vertices':6,'compact_H4_prism_vertices':8,
        'prism_common_type':'5A1','new_pair_counts':replay['new_pair_counts'],
        'nested_changed_enclosures':count,'new_affine_circuits':0,
        'strict_escape_norm':replay['remaining_escape']['norm'],
        'scope':'Independent six-facet prism Gram/norm/vertex checks and source review, followed by rigorous nested160-digit replay of the existing interval implementation. Whole H9 remains infinite-volume.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['nested_changed_enclosures'])
