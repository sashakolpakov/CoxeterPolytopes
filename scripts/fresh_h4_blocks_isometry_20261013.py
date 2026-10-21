#!/usr/bin/env python3
"""Recover exact integral H4/E8 transports; do not count new bases as hosts."""
import itertools
import json
from pathlib import Path
import networkx as nx
from coxeter_polytopes.exact.golden import Golden as G, PHI
from coxeter_polytopes.exact.matrix import determinant
from scripts.fresh_h4_blocks_h9_search_20261013 import inverse, mv, pair, Host

ROOT=Path(__file__).resolve().parents[1]
def matrix(rows):return tuple(tuple(G.from_json(x) for x in r) for r in rows)
def mm(a,b):return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))),G(0))
                             for j in range(len(b[0]))) for i in range(len(a)))
def transpose(a):return tuple(zip(*a))
def encode(a):return [[x.to_json() for x in row] for row in a]
def graph(m):
    g=nx.Graph();g.add_nodes_from(range(len(m)))
    for i in range(len(m)):
        for j in range(i):
            if m[i][j]:g.add_edge(i,j,weight=m[i][j])
    return g


def run():
    new=json.loads((ROOT/'results/fresh_h4_blocks_h8_2phi_h30_20261013.json').read_text())
    old=json.loads((ROOT/'data/known_examples/bugaenko_compact_h8.json').read_text())
    b=matrix(new['gram']);r=matrix(new['all_walls'])
    a=matrix(old['lattice_gram']);s=matrix(old['simple_roots'])
    n=tuple(tuple(pair(x,b,y) for y in r) for x in r)
    iso=nx.algorithms.isomorphism.GraphMatcher(graph(n),graph(matrix(old['wall_gram'])),
                                              edge_match=lambda x,y:x['weight']==y['weight'])
    assert iso.is_isomorphic();mapping=iso.mapping
    indices=next(row for row in itertools.combinations(range(11),9)
                 if determinant(tuple(r[i] for i in row)))
    rn=tuple(tuple(r[j][i] for j in indices) for i in range(9))
    so=tuple(tuple(s[mapping[j]][i] for j in indices) for i in range(9))
    t=mm(so,inverse(rn))
    assert mm(mm(transpose(t),a),t)==b
    assert all(mv(t,r[i])==s[mapping[i]] for i in range(11))
    assert all(x.is_integral for row in (*t,*inverse(t)) for x in row)
    record={'verdict':'EXACT_INTEGRAL_BUGAENKO_H8_ISOMETRY','wall_map':mapping,
            'basis_walls':indices,'transform':encode(t),'determinant':determinant(t).to_json(),
            'integral':True,'inverse_integral':True,'gram_congruence_verified':True}
    (ROOT/'results/fresh_h4_blocks_h8_bugaenko_match_20261013.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    extended=tuple(tuple(t[i][j] if i<9 and j<9 else G(i==j) for j in range(10)) for i in range(10))
    summaries=[]
    for name,alpha,binary,reverse in (
        ('2phi',2*PHI,None,False),
        ('phi',PHI,((-PHI/2,(PHI-1)/2),(-G(1)/2,PHI/2)),False),
        ('sqrt5',2*PHI-1,((PHI,G(1)),(PHI,2*PHI-1)),True),
    ):
        f=[[G(i==j) for j in range(10)] for i in range(10)]
        if binary is not None:
            for i in range(2):
                for j in range(2):f[(0,9)[i]][(0,9)[j]]=binary[i][j]
        transport=mm(mm(f,extended),inverse(f)) if reverse else mm(mm(inverse(f),extended),f)
        source=Host(alpha).gram
        target=tuple(tuple(-alpha if i==j==0 else G(1) if i==j==9 else
                           a[i][j] if i<9 and j<9 else G(0) for j in range(10)) for i in range(10))
        assert mm(mm(transpose(transport),target),transport)==source
        assert all(x.is_integral for row in (*transport,*inverse(transport)) for x in row)
        rec={'verdict':'EXACT_INTEGRAL_ISOMETRY_TO_PREVIOUS_E8_HOST','alpha':alpha.to_json(),
             'source':'<-alpha>+H4+H4+<1>', 'target':'<-alpha>+E8+<1>',
             'transform':encode(transport),'target_gram':encode(target),
             'determinant':determinant(transport).to_json(),'integral':True,'inverse_integral':True,
             'gram_congruence_verified':True,
             'scope':'Identical integral orthogonal and reflection groups under this transport; different Vinberg basepoint.'}
        (ROOT/f'results/fresh_h4_blocks_h9_{name}_isometry_20261013.json').write_text(json.dumps(rec,indent=2,sort_keys=True)+'\n')
        summaries.append({'alpha':name,'verdict':rec['verdict']})
    print(json.dumps({'H8':'EXACT_INTEGRAL_BUGAENKO_ISOMETRY','H9':summaries}))


if __name__=='__main__':run()
