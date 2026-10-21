#!/usr/bin/env python3
"""Symbolic bridge review and 160-digit replay of the same interval engine."""
import argparse
import json
from pathlib import Path
import sys
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_remote_next_20261013_bridge_verify as C

V=C.V
OUT=V.ROOT/'results/fresh_h9_remote_next_20261013_bridge_review.json'


def run():
    d,e,q=sp.symbols('d e q',positive=True)
    r0=sp.Matrix([-1/d,-sp.sqrt(q/e)])
    r1=sp.Matrix([-sp.sqrt(q/d),1/sp.sqrt(d*e)])
    gram=sp.diag(d,e)
    assert sp.simplify((r0.T*gram*r1)[0])==0
    assert sp.simplify((r0.T*gram*r0)[0]-(q+1/d))==0
    assert sp.simplify((r1.T*gram*r1)[0]-(q+1/d))==0
    assert sp.simplify(d*r0[0])==-1
    p=sp.symbols('p');A=sp.Matrix([[2,p],[p,2]])
    assert sp.simplify(A*sp.Matrix([[2,-p],[-p,2]])/(4-p*p))==sp.eye(2)
    source=C.PREFIX+'_caps.json';saved=V.read(C.PREFIX+'_verify.json')
    assert saved.pop('source_sha256')==V.sha(source)
    V.SCALE=10**160
    V.PHI=V.div(V.add(V.rat(1),V.sqrt(V.rat(5))),V.rat(2));V.radical.cache_clear()
    replay=C.check(V.read(source));count=0;factor=10**80
    def compare(old,new,path='root'):
        nonlocal count
        if isinstance(old,dict):
            assert isinstance(new,dict)and old.keys()==new.keys(),path
            for key in old:compare(old[key],new[key],path+'.'+key)
        elif isinstance(old,list):
            assert isinstance(new,list)and len(old)==len(new),path
            if len(old)==2 and all(isinstance(x,str)for x in old):
                try:a,b=map(int,old);c,d=map(int,new)
                except ValueError:pass
                else:
                    assert a*factor<=c<=d<=b*factor,(path,'higher precision is not nested')
                    count+=1;return
            for i,(a,b)in enumerate(zip(old,new,strict=True)):compare(a,b,path+f'[{i}]')
        else:assert old==new,path
    compare(saved,replay)
    assert count==448
    return {'status':'PASS_SYMBOLIC_BRIDGE_REVIEW_AND_160_DIGIT_REPLAY',
            'source':source,'source_sha256':V.sha(source),
            'verifier_source_sha256':V.sha('scripts/fresh_h9_remote_next_20261013_bridge_verify.py'),
            'facets':88,'new_norms':[2,2],'new_mutual_label':2,
            'new_old_pair_counts':replay['new_old_pair_counts'],
            'finite_vertices':2,'finite_vertex_types':['9A1','9A1'],
            'finite_edge':replay['finite_edge'],'nested_interval_enclosures':count,
            'remaining_escape_norm':replay['remaining_escape']['norm'],
            'strict_escape_inequalities':len(replay['remaining_escape']['strict_pairings']),
            'scope':'Separate symbolic proof of the projection inverse,norms,mutual orthogonality and pi3 anchor; source review and160-digit execution of the existing interval engine. All448 saved80-digit enclosures contain the new enclosures, with exact metadata equality. This is not a separate interval arithmetic implementation. The strict escape still proves infinite volume.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['nested_interval_enclosures'])
