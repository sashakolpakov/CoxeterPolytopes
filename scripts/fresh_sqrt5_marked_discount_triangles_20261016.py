#!/usr/bin/env python3
"""Actual V-arithmetic long-white triangles obstruct two uniform discounts.

These are infinite-volume local chambers, never compact/finite-volume
examples or obstructions to compensated weighting rules.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_sqrt5_d4_two_e8_self_edge_20261016 as E
import fresh_sqrt5_marked_discount_census_20261016 as D

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_marked_discount_triangles_20261016'
CENSUS='results/fresh_sqrt5_marked_discount_census_20261016.json'
TRI='results/fresh_sqrt5_triangle_weights_20261015.json'
Q,phi,C=E.Q,E.phi,D.C


def saved_vector(v):return [x.saved() for x in v]
def saved_matrix(g):return [saved_vector(r) for r in g]
def signature(g):
    p=E.pivots(g);signs=[x.sign() for x in p]
    return dict(positive=signs.count(1),negative=signs.count(-1),zero=signs.count(0),pivots=saved_vector(p))


def run():
    census=json.loads((ROOT/CENSUS).read_text());tri=json.loads((ROOT/TRI).read_text())
    assert census['status']=='PASS_EXHAUSTIVE_374_HIGH_MARKED_TYPES_AND_DISCOUNT_BUDGETS'
    known={C.code(g):name for name,g in C.models() if len(g)<=21}
    records=[]
    for policy in ('long_local33_distance2','long_adjacent3'):
        assert next(r for r in census['policies'] if r['policy']==policy)['exact_marked_maximum']=='35'
        if policy=='long_local33_distance2':
            core=5;black=[0,1,2];white=[3,4]
            labels=dict(zip(combinations(range(core),2),tri['markings'][46]['labels']))
            g=E.graph_gram(22,[(a,b,-phi if m==5 else -1) for (a,b),m in labels.items() if m!=2])
        else:
            core=3;black=[0,1,2];white=[]
            labels={(0,1):2,(0,2):3,(1,2):10}
            g=E.graph_gram(22,[(0,2,-1),(1,2,-2-phi)])
            g[1][1]=2*(2+phi)
        tail=list(range(core,22));white+=tail
        labels[tuple(sorted((0,tail[0])))]=3;g[0][tail[0]]=g[tail[0]][0]=Q(-1)
        for a,b in zip(tail,tail[1:]):
            labels[a,b]=3;g[a][b]=g[b][a]=Q(-1)
        assert len(white)==19 and len(tail)==(17 if core==5 else 19)
        assert all(g[i][i].sign()>0 for i in range(22))
        assert all(g[i][j].sign()<=0 for i in range(22) for j in range(i))
        physical=signature(g);conjugate=signature([[v.conj() for v in r] for r in g])
        assert (physical['positive'],physical['negative'])==(21,1)
        assert (conjugate['positive'],conjugate['negative'])==(22,0)
        cartan=[[2*g[i][j]/g[i][i] for j in range(22)] for i in range(22)]
        assert all(v.a.denominator==v.b.denominator==1 for row in cartan for v in row)
        cycle=next((i,j,cartan[i][j]*cartan[j][i]) for i in range(22) for j in range(i)
                   if (cartan[i][j]*cartan[j][i]).b)
        # A tree has no independent long cycles: all cyclic products are
        # products of edge two-cycles, already integral inQ(phi).
        nonzero_edges=[(i,j) for i in range(22) for j in range(i) if g[i][j]]
        assert len(nonzero_edges)==21
        inverse=E.inverse(g)
        dual=[[-inverse[j][i] for j in range(22)] for i in range(22)]
        vertices=[];points=[];payments=[]
        for omitted in black:
            v=dual[omitted];gv=E.mv(g,v)
            assert gv==[-Q(i==omitted) for i in range(22)]
            norm=E.dot(v,gv);assert norm.sign()<0
            support=[i for i in range(22) if i!=omitted]
            p=E.pivots(E.sub(g,support));assert all(x.sign()>0 for x in p)
            forest=C.pieces(C.graph(support,labels))
            types=sorted(known[C.code(h)] for h in forest)
            pair=[i for i in black if i!=omitted]
            h=next(h for h in forest if pair[0] in h)
            weight=C.base(h,*pair) if pair[1] in h else F(0)
            name=known[C.code(h)]
            changed=pair[1] in h and D.eligible(policy,name,h,*pair)
            baseline=F(0);discounted=F(0);eligible_count=0
            for component in forest:
                typ=known[C.code(component)]
                for a,b in combinations(component,2):
                    if 21 in (a,b):continue
                    w=C.base(component,a,b);flag=D.eligible(policy,typ,component,a,b)
                    baseline+=w;discounted+=w-D.DELTA*flag;eligible_count+=flag
            assert discounted<=35
            payments.append(weight-D.DELTA*changed)
            vertices.append(dict(omitted=omitted,point=saved_vector(v),norm=norm.saved(),
                 active=support,finite_types=types,positive_corner_pivots=saved_vector(p),
                 triangle_pair=pair,triangle_pair_component=name,baseline_corner_payment=str(weight),
                 discounted_corner_payment=str(payments[-1]),changed=changed,
                 facet_baseline=str(baseline),eligible_facet_corners=eligible_count,
                 discounted_facet_total=str(discounted),facet_root=21))
            points.append(v)
        pairings=[[E.dot(v,E.mv(g,w)) for w in points] for v in points]
        assert all(z.sign()<0 for row in pairings for z in row)
        c0=[sum((v[j] for v in points),Q()) for j in range(22)]
        assert E.dot(c0,E.mv(g,c0)).sign()<0
        all_dual=[sum((v[j] for v in dual),Q()) for j in range(22)]
        for power in range(1,121):
            epsilon=F(1,1<<power)
            center=[a+epsilon*b for a,b in zip(c0,all_dual)]
            cp=E.mv(g,center);cn=E.dot(center,cp)
            if cn.sign()<0 and all(E.dot(v,cp).sign()<0 for v in points):break
        else:raise AssertionError('No strict center')
        assert all(x.sign()<0 for x in cp)
        facets=[]
        for i in range(22):
            norm=g[i][i]*g[i][i]*cn-g[i][i]*cp[i]*cp[i]
            future=g[i][i]*cn-cp[i]*cp[i]
            assert norm.sign()<0 and future.sign()<0
            assert all((g[i][i]*cp[j]-cp[i]*g[j][i]).sign()==(0 if i==j else -1) for j in range(22))
            facets.append(dict(root=i,norm=norm.saved(),center_pairing=future.saved()))
        # The last tail dual is weakly feasible and spacelike. Perturbing
        # toward the strict center gives an explicit strict escape.
        weak=dual[21];weaknorm=E.dot(weak,E.mv(g,weak));assert weaknorm.sign()>0
        for epower in range(1,121):
            eta=F(1,1<<epower);escape=[a+eta*b for a,b in zip(weak,center)]
            escape_pairs=E.mv(g,escape);escape_norm=E.dot(escape,escape_pairs)
            if escape_norm.sign()>0 and E.dot(escape,cp).sign()<0:break
        else:raise AssertionError('No strict escape')
        assert all(x.sign()<0 for x in escape_pairs)
        white_pivots=E.pivots(E.sub(g,white));assert all(x.sign()>0 for x in white_pivots)
        projected=E.schur(g,white,black)
        assert signature(projected)['negative']==1
        assert all(all(v.sign()>0 for v in E.pivots(E.sub(projected,list(pair))))
                   for pair in combinations(range(3),2))
        assert sum(payments)==F(35,19)<2
        records.append(dict(policy=policy,dimension=21,facet_dimension=20,facet_root=21,
            core_size=core,tail_length=len(tail),tail_attachment=0,white=white,black=black,
            labels=[dict(pair=list(p),label=m) for p,m in sorted(labels.items()) if m!=2],
            Gram=saved_matrix(g),physical=physical,conjugate=conjugate,
            integral_Cartan=saved_matrix(cartan),groundfield_cycle=dict(pair=list(cycle[:2]),product=cycle[2].saved()),
            finite_triangle_vertices=vertices,vertex_mutual_pairings=saved_matrix(pairings),
            projected_triangle_Gram=saved_matrix(projected),white_positive_pivots=saved_vector(white_pivots),
            strict_center=dict(point=saved_vector(center),norm=cn.saved(),pairings=saved_vector(cp),epsilon=str(epsilon)),
            actual_facets=facets,strict_escape=dict(point=saved_vector(escape),norm=escape_norm.saved(),
                pairings=saved_vector(escape_pairs),center_pairing=E.dot(escape,cp).saved(),eta=str(eta)),
            required_triangle_payment='2',actual_discounted_payment=str(sum(payments)),deficit='3/19'))
    paths=[f'scripts/{TAG}.py',f'scripts/{D.TAG}.py',CENSUS,TRI,
           'scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016.py']
    return dict(status='PASS_TWO_ACTUAL_V_ARITHMETIC_TRIANGLES_OBSTRUCT_UNCOMPENSATED_DISCOUNTS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},discount='3/19',witnesses=records,
        conclusion='Both uniform long-component discounts attain the desired vertex bound35 but fail on an actual compact triangle in a dimension21 V-arithmetic local chamber: their sum is35/19 instead of2.',
        next_compensation='The distance-two witness couples the changed long A21 corner to an H4 corner with an external marked A17 tail. The adjacent witness couples the changed long A21 corner to an I2(10) corner with an external marked A19 tail. A compensated rule must repay at least3/19 somewhere on each triangle; these specific contexts have vertex-budget slack.',
        scope='Explicit local V-arithmetic Gram data with full Lorentz signature,positive-definite conjugate,integral Cartan entries,actual facets,three future finite triangle corners,and strict infinite-volume escapes. This does not realize a compact arithmetic polytope,exclude all compensated discounts,or give a dimension bound. No old paid-corner pruning is used.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(p.read_text())
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
