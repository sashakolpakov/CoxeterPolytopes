#!/usr/bin/env python3
"""Literal comparison to Chen--Labbe Figure14, bottom-right rank10 graph.

The figure transcription is a human-checked primary-source input. All
matrix/level/polar/cusp-count comparisons below are exact and replayable.
"""
import argparse, hashlib, json
from itertools import combinations
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_all_ideal_seven_simplex_literature_20261017'
BASE='fresh_all_ideal_seven_simplex_realization_20261017'

def inertia(A):
    A=s.Matrix(A);out=[0,0,0]
    while A.rows:
        k=next((i for i in range(A.rows) if A[i,i]),None)
        if k is not None:
            ids=[k]+[i for i in range(A.rows) if i!=k];A=A.extract(ids,ids)
            a=A[0,0];out[0 if a>0 else 1]+=1
            A=A[1:,1:]-A[1:,:1]*A[:1,1:]/a
        else:
            pair=next(((i,j) for i,j in combinations(range(A.rows),2) if A[i,j]),None)
            if pair is None:out[2]+=A.rows;break
            ids=list(pair)+[i for i in range(A.rows) if i not in pair];A=A.extract(ids,ids)
            out[0]+=1;out[1]+=1
            A=A[2:,2:]-A[2:,:2]*A[:2,:2].inv()*A[:2,2:]
    return out

def run():
    d=json.loads((ROOT/f'results/{BASE}.json').read_text())
    # In the printed figure, number the octagon clockwise from its top
    # vertex. Node8 is the interior leaf at the left cycle vertex6;
    # node9 is the exterior leaf at the right cycle vertex2.
    edges=[(i,(i+1)%8) for i in range(8)]+[(6,8),(2,9)]
    G=2*s.eye(10)
    for i,j in edges:G[i,j]=G[j,i]=-1
    figure_to_repository=[1,0,7,6,5,4,3,2,9,8]
    target=s.Matrix(d['original_gram'])
    assert G==target.extract(figure_to_repository,figure_to_repository)
    deletions=[]
    for p in combinations(range(10),2):
        ids=[i for i in range(10) if i not in p];sig=inertia(G.extract(ids,ids))
        assert sig[1]==0
        deletions.append(dict(deleted=list(p),inertia=sig))
    assert inertia(G.extract([i for i in range(10) if i!=8],[i for i in range(10) if i!=8]))[1]==1
    J=target.inv();N=s.Matrix(d['root_coordinates'])
    assert N[10,:].T==-2*J[:,8] and N[11,:].T==-2*J[:,9]
    assert J[8,8]==J[9,9]==s.Rational(1,2) and J[8,9]==-s.Rational(1,2)
    t,u=s.symbols('t u')
    assert s.expand((s.Matrix([t,u]).T*J[8:,8:]*s.Matrix([t,u]))[0]-(t-u)**2/2)==0
    nonsimple=[i for i,r in enumerate(d['vertices']) if len(r['active'])>9]
    assert len(nonsimple)==3 and all(d['vertices'][i]['type']=='ideal' for i in nonsimple)
    assert max(len(f['vertices']) for f in d['facets'])==21<len(d['vertices'])-1
    paths=[f'scripts/{TAG}.py',f'results/{BASE}.json',f'results/{BASE}.md',
           f'data/references/{TAG}_chen_labbe.pdf']
    return dict(status='PASS_LITERAL_PUBLISHED_LEVEL2_DIAGRAM_AND_POLAR_TRUNCATION_IDENTIFICATION',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        primary_url='https://arxiv.org/pdf/1310.8608v6',
        journal_doi='https://doi.org/10.1007/s10711-014-0004-1',
        primary_location='Figure14, bottom-right graph, printed page24 of arXiv v6',
        figure_transcription=dict(cycle=list(range(8)),leaf_edges=[[6,8],[2,9]],all_labels=3,
            vertex_numbering='Octagon clockwise from top;8 interior-left leaf,9 exterior-right leaf'),
        figure_to_repository=figure_to_repository,two_node_deletions=deletions,level=2,
        canonical_polar_root_identity=True,nonsimple_vertex_indices=nonsimple,
        is_pyramid=False,hyperideal_edge_square='(t-u)^2/2; only its midpoint is ideal',
        qualification='Published original ten-root diagram identified literally. Twelve-wall polar completion follows canonically; no located explicit source naming its all-ideal7-face, and no novelty or exhaustive-literature claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(out.read_text())==d
    else:out.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],'saved_replay=',a.verify)
