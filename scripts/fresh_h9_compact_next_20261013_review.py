#!/usr/bin/env python3
"""Separate four-cap identities, tetrahedral Gram check and 160-digit replay."""
import argparse,json
from pathlib import Path
import sys
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_compact_next_20261013_verify as C
V=C.V
OUT=V.ROOT/'results/fresh_h9_compact_next_20261013_review.json'

def run():
    phi=(1+sp.sqrt(5))/2
    L=sp.Matrix([[2,-1,0,0],[-1,2,-phi,0],[0,-phi,2,-1],[0,0,-1,2]])
    assert sp.simplify(L.det()-(5-4*phi))==0
    assert L.det()<0
    for omit in range(4):
        ids=[i for i in range(4)if i!=omit];M=L.extract(ids,ids)
        assert all(sp.simplify(M[:k,:k].det())>0 for k in range(1,4))
    assert all(sp.simplify(x)<0 for x in L.inv())
    # Euclidean coordinates in the boosted spacelike three-frame.
    Z=sp.Matrix([[sp.sqrt(sp.Rational(2,3)),sp.sqrt(sp.Rational(4,3)),0],
                 [sp.sqrt(sp.Rational(2,3)),-1/sp.sqrt(3),1],
                 [sp.sqrt(sp.Rational(2,3)),-1/sp.sqrt(3),-1]])
    assert sp.simplify(Z*Z.T)==2*sp.eye(3)
    qa=(4-2*phi)/3;nu=1/qa-2
    assert qa>0 and sp.simplify(nu)>0 and sp.simplify(1/qa-nu)==2
    k=sp.Rational(13,4);ct=(k+1/k)/2;st=(k-1/k)/2
    boost=sp.Matrix([[ct,st],[st,ct]])
    assert boost.T*sp.diag(-1,1)*boost==sp.diag(-1,1)
    source=C.PREFIX+'_caps.json';saved=V.read(C.PREFIX+'_verify.json')
    assert saved.pop('source_sha256')==V.sha(source)
    V.SCALE=10**160
    V.PHI=V.div(V.add(V.rat(1),V.sqrt(V.rat(5))),V.rat(2));V.radical.cache_clear()
    replay=C.check(V.read(source));count=0;factor=10**80
    def compare(old,new,path='root'):
        nonlocal count
        if old==new:return
        if isinstance(old,dict):
            assert isinstance(new,dict)and old.keys()==new.keys(),path
            for key in old:compare(old[key],new[key],path+'.'+key)
        elif isinstance(old,list):
            assert isinstance(new,list)and len(old)==len(new),path
            if len(old)==2 and all(isinstance(x,str)for x in old):
                try:a,b=map(int,old);c,d=map(int,new)
                except ValueError:pass
                else:
                    assert a*factor<=c<=d<=b*factor,path
                    count+=1;return
            for i,(a,b)in enumerate(zip(old,new,strict=True)):compare(a,b,path+f'[{i}]')
        else:assert old==new,path
    compare(saved,replay)
    assert len(replay['actual_finite_vertices'])==4
    assert replay['new_vertex_edge_frontier']['finite_edges']==6
    assert replay['new_vertex_edge_frontier']['open_rays']==24
    return {'status':'PASS_FOUR_CAP_IDENTITIES_TETRAHEDRON_AND_160_DIGIT_REPLAY',
      'source':source,'source_sha256':V.sha(source),'verifier_sha256':V.sha('scripts/fresh_h9_compact_next_20261013_verify.py'),
      'facets':94,'compact_face':'H3 tetrahedron with diagram3-5-3',
      'new_actual_vertices':4,'closed_face_edges':6,'open_new_vertex_directions':24,
      'new_psd_subsets_checked':960,'new_affine_circuits':0,
      'nested_changed_enclosures':count,'strict_escape_norm':replay['remaining_escape']['norm'],
      'scope':'Separate symbolic norm/orthogonality/boost and compact tetrahedral Gram checks; source review and160-digit replay of the existing interval arithmetic. All changed enclosures nest and exact metadata agree. Not a separate interval implementation. The strict escape still proves infinite volume.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['nested_changed_enclosures'])
