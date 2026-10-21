#!/usr/bin/env python3
"""Dimension-independent three-timelike-leaf obstruction on literal star0."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_three_large_blocks_20261015'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
BATCH='results/fresh_h9_star_pi4_resume_20261015_search_batch_015.json'
BLOCKS=[(3,4),(6,7),(9,10),(12,13)]
LEAVES=[14,15,16,17]


def qm(x,y):return (x[0]*y[0]+3*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
Z=(F(0),F(0));O=(F(1),F(0));H=(F(-1,2),F(0));R=(F(0),F(-1,2))


def add(p,q):
    r=p.copy()
    for k,v in q.items():r[k]=r.get(k,F(0))+v
    return {k:v for k,v in r.items() if v}


def mul(p,q):
    r={}
    for a,x in p.items():
        for b,y in q.items():
            k=tuple(i+j for i,j in zip(a,b));r[k]=r.get(k,F(0))+x*y
    return {k:v for k,v in r.items() if v}


def polynomial_det(M):
    out={};n=len(M)
    for p in permutations(range(n)):
        s=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        term={(0,0,0):F(s)}
        for i,j in enumerate(p):term=mul(term,M[i][j])
        out=add(out,term)
    return out


def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    data=json.loads((ROOT/GEOM).read_text())['representatives'][0]['geometry']
    V=[set(v['incident']) for v in data['vertices']]
    MN={tuple(m) for m in data['minimal_nonfaces']}
    assert len(V)==378
    G=[[O,R,Z],[R,O,H],[Z,H,O]];u=[(F(0),F(1)),(F(2),F(0)),O]
    assert all(sum((qm(a,b)[j] for a,b in zip(row,u)),F(0))==0 for row in G for j in (0,1))
    one={(0,0,0):F(1)};c={(1,0,0):F(-1)};a={(0,1,0):F(-1)};b={(0,0,1):F(-1)}
    polynomial=polynomial_det([[one,c,a],[c,one,b],[a,b,one]])
    expected={(0,0,0):F(1),(2,0,0):F(-1),(0,2,0):F(-1),(0,0,2):F(-1),(1,1,1):F(-2)}
    assert polynomial==expected
    def face(nodes):
        nodes=sorted(nodes);i=next(i for i,v in enumerate(V) if set(nodes)<=v)
        return dict(support=nodes,actual_vertex=i)
    rows=[];families=[]
    for chosen in combinations(range(4),3):
        cross=[];foreign=[]
        for i,j in combinations(chosen,2):
            for x in BLOCKS[i]:
                for y in BLOCKS[j]:
                    p=next(p for p in BLOCKS[i] if p!=x)
                    cross.append(dict(ordered_support=[p,x,y],large_pair=sorted([p,x]),
                        trial_nonright_pair=sorted([x,y]),**face([p,x,y])))
        for i in chosen:
            for j in chosen:
                if i==j:continue
                for x in BLOCKS[j]:
                    p=next(p for p in BLOCKS[j] if p!=x)
                    foreign.append(dict(ordered_support=[p,x,LEAVES[i]],large_pair=sorted([p,x]),
                        trial_nonright_pair=sorted([x,LEAVES[i]]),**face([p,x,LEAVES[i]])))
        own=[]
        for i in chosen:
            triple=tuple(sorted((*BLOCKS[i],LEAVES[i])))
            assert triple in MN
            assert all(any(set(pair)<=v for v in V) for pair in combinations(triple,2))
            own.append(dict(block=list(BLOCKS[i]),leaf=LEAVES[i],minimal_nonface=list(triple)))
        assert len(cross)==len(foreign)==12
        source_roots={x for i in chosen for x in BLOCKS[i]};assert len(source_roots)==6
        f=dict(id='three_large_blocks_'+'_'.join(map(str,chosen)),origin='three_timelike_cut_leaves',
               fixed=[],bounds=[list(BLOCKS[i])+[6] for i in chosen])
        families.append(f)
        rows.append(dict(chosen_blocks=list(chosen),pd_six_root_support=sorted(source_roots),
            cross_right_guards=cross,foreign_right_guards=foreign,own_Lanner_triples=own,
            leaf_spherical_triangle=face([LEAVES[i] for i in chosen]),family=f))
    assert len(rows)==4
    batch=json.loads((ROOT/BATCH).read_text())
    assert batch['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR' and batch['cumulative_additional_nodes']==119877
    labels={tuple(r[:2]):r[2] for r in batch['solution']['pair_labels']}
    hits=[f['id'] for f in families if all(labels[a,b]>=m for a,b,m in f['bounds'])]
    assert len(hits)==4
    sources=[f'scripts/{TAG}.py',GEOM,BATCH]
    return dict(status='PASS_FOUR_THREE_LARGE_BLOCK_LORENTZ_OBSTRUCTIONS',
        sources={p:sha(p) for p in sources},case=0,placements=rows,families=families,
        relaxed_path_Gram=[[[str(t) for t in x] for x in row] for row in G],
        positive_kernel=[[str(t) for t in x] for x in u],kernel_basis='1,sqrt3',
        own_triple_determinant=[dict(powers=list(k),coefficient=str(v)) for k,v in sorted(polynomial.items())],
        exact_analytic_proof={
            'cross_right':'A nonright extra contact gives an actual6--3 path, whose positive null vector has nonpositive quadratic form under every strengthening. Hence all indicated cross/foreign pairs are right.',
            'own_residual':'For ownI2(c) with0<=c<1 and leaf cosinesa,b>=0, rho=(1-c*c-a*a-b*b-2*c*a*b)/(1-c*c)<0 because the own minimal nonface is Lanner.',
            'pair_right_impossible':'The PD6 Schur complement of two cut leaves would be diagonal with two strictly negative entries if their ambient pair were right; that gives two negative directions.',
            'triangle_failure':'Consequently all three actual leaf-triangle pairs have integer labels>=3. Its all-ones quadratic form is3-2*(cosA+cosB+cosC)<=0, contradicting positive definiteness.'},
        path_guard_occurrences=96,own_Lanner_occurrences=12,actual_leaf_triangles=4,
        current_cursor=119877,current_survivor_hits=hits,
        theorem_inputs=['Finite proper faces of a compact Coxeter chamber have positive-definite normal Gram matrices.',
                        'A size-three minimal nonface in a compact Coxeter nerve is a Lanner triple: all proper pairs are spherical and its determinant is strictly negative.',
                        'A Lorentzian normal space has negative index one; Schur complement by a positive-definite support preserves the negative index.'],
        scope='All finite internal labels>=6 on any three indicated blocks; all own labels, other labels and dotted distances unrestricted. Dimension-independent under the literal spherical-face/Lanner premises. Four exact predicates for compact star0, not whole-incidence exclusion or an unconditional dimension bound.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'],out['path_guard_occurrences'],out['current_cursor'])
