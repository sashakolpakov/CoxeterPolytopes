#!/usr/bin/env python3
"""Conditional face-coverage reduction and complete local triangle-growth alphabet.

No feasible weight rule or complete white-extension enumeration is asserted.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_marked_discount_census_20261016 as D
import fresh_sqrt5_quad_filters_20261015 as G

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_discounted_face_domain_reduction_20261016'


def run():
    models=[(n,g) for n,g in D.C.models() if len(g)<=21]
    assert len(models)==68
    checks=[];eligible=0
    for name,g in models:
        assert max(map(len,g.values()),default=0)<=3
        branch_nodes=[u for u in g if len(g[u])==3]
        assert all(m==3 for u in branch_nodes for m in g[u].values())
        pairs=[]
        for a,b in combinations(g,2):
            if D.eligible('long_adjacent3',name,g,a,b):
                assert len(g)>=6 and g[a].get(b)==3 and D.C.base(g,a,b)==1
                pairs.append([a,b]);eligible+=1
        checks.append(dict(component=name,rank=len(g),branch_nodes=branch_nodes,eligible_pairs=pairs))
    # A positive null combination creates an ideal point in any coherent
    # acute chamber: the common strict future interior fixes its orientation.
    affine_a2=[[2 if i==j else -1 for j in range(3)] for i in range(3)]
    affine_d4=[[2 if i==j else -int(i==0 or j==0) for j in range(5)] for i in range(5)]
    nulls=[]
    for name,matrix,vector in [('affine_A2',affine_a2,[1,1,1]),('affine_D4',affine_d4,[2,1,1,1,1])]:
        assert all(x>0 for x in vector)
        assert [sum(a*b for a,b in zip(row,vector)) for row in matrix]==[0]*len(vector)
        nulls.append(dict(diagram=name,normalized_gram=matrix,positive_null_coefficients=vector))
    alphabet=[3,4,5,6,10]
    cat=json.loads((ROOT/'results/fresh_sqrt5_lanner_catalogue_20261015.json').read_text())
    assert cat['finite_alphabet']==[2]+alphabet
    triangles=json.loads((ROOT/'results/fresh_sqrt5_triangle_weights_20261015.json').read_text())
    assert len(triangles['markings'])==76
    assert all(sum(c['h4_only_weight'] for c in row['corners'])>=2 for row in triangles['markings'])
    degree2=[list(labels) for labels in product(alphabet,repeat=2)
        if (G.K(4)-sum((G.SQUARES[m] for m in labels),G.K())).sign()>0]
    assert degree2==[[3,3],[3,4],[3,5],[4,3],[5,3]]
    degree3_black=[list(labels) for labels in product(alphabet,repeat=3)
        if all((G.K(4)-G.SQUARES[a]-G.SQUARES[b]).sign()>0 for a,b in combinations(labels,2))]
    assert degree3_black==[[3,3,3],[3,3,4],[3,3,5],[3,4,3],[3,5,3],[4,3,3],[5,3,3]]
    # For a k-sided face, the integral baseline is B and each of its d
    # discounted corners contributes exactly1 to B. Other additions are
    # nonnegative. At delta<=2/3, the guaranteed payment is B-2d/3.
    tables={};unsafe={}
    for sides,target in [(3,2),(4,1)]:
        table=[];bad=[]
        for baseline in range(0,9):
            if sides==3 and baseline<2:continue  # frozen complete triangle theorem
            for discounted in range(min(baseline,sides)+1):
                if sides==3 and discounted==3:continue  # compact affine-A2 exclusion
                lower=F(baseline)-F(2*discounted,3)
                record=dict(baseline=baseline,discounted_corners=discounted,
                    guaranteed_payment=str(lower),safe=lower>=target)
                table.append(record)
                if lower<target:bad.append([baseline,discounted])
        if sides==3:assert bad==[[2,1],[2,2],[3,2]]
        else:assert bad==[[0,0],[1,1],[2,2]]
        tables[str(sides)]=table;unsafe[str(sides)]=bad
    # B>=9 is automatically safe for either face, since d<=4.
    assert F(9)-F(8,3)>2
    paths=[f'scripts/{TAG}.py',f'scripts/{D.TAG}.py',f'scripts/{D.C.TAG}.py',
        f'scripts/{G.TAG}.py','results/fresh_sqrt5_quad_filters_20261015_review.json',
        'results/fresh_sqrt5_lanner_catalogue_20261015.json',
        'results/fresh_sqrt5_triangle_weights_20261015.json',
        'results/fresh_sqrt5_triangle_weights_20261015_review.json']
    return dict(status='PASS_CONDITIONAL_DISCOUNTED_FACE_DOMAIN_REDUCTION',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        assumptions=['distance-two discount exactly0','0<=adjacent discount<=2/3',
            'adjacent discount only on a3-edge in an entire spherical component of rank>=6',
            'all added compensation and cross payments nonnegative',
            'compact arithmetic ground field Q(sqrt5), rank21 finite vertex diagrams'],
        finite_component_checks=checks,eligible_pair_count=eligible,positive_affine_null_witnesses=nulls,
        baseline='Graph-distance1/2 indicator plus the two entire-H4 corrections; all76 minimal triangle markings have total>=2 and the frozen extension-monotonicity proof applies.',
        face_payment_tables=tables,remaining_cases=unsafe,
        triangle_extension_domains=dict(degree1_labels=alphabet,degree2_ordered_labels=degree2,
            degree3_not_all_black=[[3,3,3]],degree3_all_black_ordered_labels=degree3_black,
            degree4_or_more='Impossible: finite corners force anydegree4white root to be thecenter of an induced affineD4; degree>=5 violates a finite corner degreebound.'),
        hereditary_pruning=['baseline>=4','baseline3 with fewer than2 black3edges','baseline2 with no black3edge'],
        growth_scope='Local contact domains for connected white-root growth only. All retained graphs still need allthree spherical corner tests, physical/conjugate Gram tests, exact isomorphism control, rankbound and eventual common-spectator/facet-mark payment evaluation.',
        scope='Conditional structural reduction, not a feasible weight rule, complete triangle/quad enumeration or dimension exclusion. It does not assert any active numerical candidate meets the hypotheses.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['remaining_cases'],flush=True)
