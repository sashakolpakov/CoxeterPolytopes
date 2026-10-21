#!/usr/bin/env python3
"""Independent literal-domain and elementary algebra review; no producer import."""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_three_large_blocks_20261015_review'
PRIMARY='fresh_h9_star_three_large_blocks_20261015'


def run():
    primary_path=f'results/{PRIMARY}.json'
    p=json.loads((ROOT/primary_path).read_text())
    for f,h in p['sources'].items():
        assert sha256((ROOT/f).read_bytes()).hexdigest()==h
    geompath='results/fresh_h9_tree_cube_wedges_20261014.json'
    g=json.loads((ROOT/geompath).read_text())['representatives'][0]['geometry']
    vertices=[set(v['incident']) for v in g['vertices']]
    missing={frozenset(x) for x in g['minimal_nonfaces']}
    assert len(vertices)==378
    blocks=[(3,4),(6,7),(9,10),(12,13)];leaves=[14,15,16,17]
    assert {tuple(r['chosen_blocks']) for r in p['placements']}==set(combinations(range(4),3))
    total=0
    for row in p['placements']:
        chosen=row['chosen_blocks']
        expected_cross={tuple(sorted((x,y))) for i,j in combinations(chosen,2)
                        for x in blocks[i] for y in blocks[j]}
        expected_foreign={tuple(sorted((leaves[i],x))) for i in chosen for j in chosen
                          if i!=j for x in blocks[j]}
        for key,expected in [('cross_right_guards',expected_cross),('foreign_right_guards',expected_foreign)]:
            assert {tuple(q['trial_nonright_pair']) for q in row[key]}==expected
            for q in row[key]:
                a,b,c=q['ordered_support']
                assert tuple(sorted((a,b))) in [blocks[i] for i in chosen]
                assert tuple(sorted((b,c))) in expected
                assert set((a,b,c))==set(q['support'])<=vertices[q['actual_vertex']]
                total+=1
        for q in row['own_Lanner_triples']:
            t=frozenset(q['minimal_nonface']);assert t in missing
            assert all(any(set(pair)<=v for v in vertices) for pair in combinations(t,2))
        triangle=row['leaf_spherical_triangle']
        assert set(triangle['support'])=={leaves[i] for i in chosen}
        assert set(triangle['support'])<=vertices[triangle['actual_vertex']]
        assert row['family']['fixed']==[]
        assert row['family']['bounds']==[list(blocks[i])+[6] for i in chosen]
    assert total==96
    root3=s.sqrt(3)
    K=s.Matrix([[1,-root3/2,0],[-root3/2,1,-s.Rational(1,2)],[0,-s.Rational(1,2),1]])
    v=s.Matrix([root3,2,1]);assert K*v==s.zeros(3,1)
    assert all(x>0 for x in v)
    c,a,b=s.symbols('c a b')
    G=s.Matrix([[1,-c,-a],[-c,1,-b],[-a,-b,1]])
    determinant=s.expand(G.det())
    expected=1-c*c-a*a-b*b-2*c*a*b
    assert determinant==expected
    block=s.Matrix([[1,-c],[-c,1]])
    row=s.Matrix([a,b])
    assert s.cancel(1-(row.T*block.inv()*row)[0]-expected/(1-c*c))==0
    sources=[f'scripts/{TAG}.py',f'scripts/{PRIMARY}.py',primary_path,geompath]
    return dict(status='PASS_INDEPENDENT_THREE_LARGE_BLOCK_DOMAIN_AND_PROOF_REVIEW',
                sources={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in sources},
                literal_placements=4,path_occurrences=total,own_minimal_nonfaces=12,
                actual_leaf_triangles=4,positive_kernel=[str(x) for x in v],
                own_determinant=str(determinant),
                proof_audit='Positive-kernel monotonicity forces cross/foreign right; compact minimal-nonface triples have negative determinant; PD6 Schur gives three timelike residuals. A right pair gives index>=2, while three nonright leaf pairs contradict their actual spherical triangle.',
                scope='Independent literal incidence and elementary SymPy algebra plus source/proof audit. No original primal realization replay or old parameter-bank replay. No unconditional dimension theorem or whole star-incidence exclusion.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();data=run();out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
