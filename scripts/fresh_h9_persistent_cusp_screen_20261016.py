#!/usr/bin/env python3
"""Exact old49 pair inventory for a permanent-cusp obstruction screen.

This initial inventory does not claim complete affine-subdiagram coverage.
"""
import argparse
import json
from collections import Counter
from itertools import combinations
import fresh_nonarithmetic_common_escape_20261016_cap_review as C

B=C.B
TAG='fresh_h9_persistent_cusp_screen_20261016'


def run():
    data=C.geometry()
    N,cov,NI,covI=[data[k] for k in ('N','cov','NI','covI')]
    norms=[B.dot(n,c) for n,c in zip(N,cov)]
    normI=[B.iv(q) for q in norms]
    pair_matrix=[[2 if i==j else None for j in range(49)] for i in range(49)]
    counts=Counter()
    pairs=[]
    for i in range(49):
        for j in range(i):
            gI=B.idot(NI[i],covI[j])
            if gI.hi<0 and (gI*gI-normI[i]*normI[j]).lo>0:
                kind='dotted';label=0
            else:
                g=B.dot(N[i],cov[j])
                if g==0:
                    kind='right';label=2
                elif B.sg(g*g-norms[i]*norms[j])==0:
                    assert B.sg(g)<0
                    kind='parallel';label=-1
                else:
                    relation=B.relation(g,norms[i]*norms[j])
                    assert relation['kind']=='finite',relation
                    kind='finite';label=relation['label']
            pair_matrix[i][j]=pair_matrix[j][i]=label
            counts[kind]+=1
            pairs.append(dict(walls=[j,i],kind=kind,label=label))
        print('exact pair row',i,flush=True)
    forbidden=[]
    for i,j in combinations(range(49),2):
        m=pair_matrix[i][j]
        if m>2 and m not in (3,4,6):
            orth=[k for k in range(49) if k not in (i,j) and pair_matrix[i][k]==pair_matrix[j][k]==2]
            edges=[[a,b,pair_matrix[a][b]] for a,b in combinations(orth,2) if pair_matrix[a][b] not in (0,2)]
            forbidden.append(dict(walls=[i,j],label=m,common_orthogonal=orth,nonright_nondotted_edges=edges))
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'scripts/{B.TAG}.py',C.PRIMARY,
           'results/fresh_nonarithmetic_common_escape_20261016_cap_review.json']
    return dict(status='PASS_EXACT_49_WALL_PAIR_INVENTORY_FOR_CUSP_SCREEN',
        sources={p:B.digest(p) for p in paths},pair_counts=dict(counts),
        pair_matrix=pair_matrix,forbidden_angle_edges=forbidden,
        scope='Exact complete pair inventory and common orthogonality domains. Does not by itself exclude affine circuits, completion or finite volume.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=B.ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
