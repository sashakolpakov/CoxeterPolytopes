#!/usr/bin/env python3
"""Exact six-simplex minimal-diagram partition and direct-white frontier.

Uses the frozen Qsqrt2 arithmetic engine and the frozen complete26 diagrams.
No new quasi-Lanner census or outside-white rank cutoff is performed.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path
import networkx as nx
import fresh_finite_edge_cusp_tradeoff12_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_six_simplex_ideal_frontier_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
Q=A.Q

def graph(L,ids):
    G=nx.Graph();G.add_nodes_from(ids)
    for a,b,m in L['edges']:
        if a in ids and b in ids:G.add_edge(a,b,label=m)
    return G

def affine_profiles(L,h,diagrams):
    n=L['order'];ids=[i for i in range(n) if i!=h];target=graph(L,ids)
    profiles={}
    for M in diagrams:
        if M['order']!=n:continue
        for z in M['ideal_deletions']:
            source=graph(M,[i for i in range(n) if i!=z])
            matcher=nx.algorithms.isomorphism.GraphMatcher(source,target,
                edge_match=lambda x,y:x['label']==y['label'])
            for iso in matcher.isomorphisms_iter():
                labels=[2]*len(ids)
                for x,y,m in M['edges']:
                    if x==z:labels[ids.index(iso[y])]=m
                    if y==z:labels[ids.index(iso[x])]=m
                key=tuple(labels)
                witness=dict(diagram=M['index'],affine_deletion=z,
                    source_to_target=[[i,iso[i]] for i in sorted(iso)])
                if key not in profiles or str(witness)<str(profiles[key]):profiles[key]=witness
    return ids,[dict(labels=list(p),witness=profiles[p]) for p in sorted(profiles)]

def run():
    census=json.loads((ROOT/f'results/{CENSUS}.json').read_text())
    diagrams=census['diagrams'];rows=[];profile_cache={};hard=[]
    for L in diagrams:
        n=L['order']
        if n<7:continue
        G=A.gram(n,L['edges']);assert A.inertia(G)==[n-1,1,0]
        invs={}
        for j in range(n):
            ids=[i for i in range(n) if i!=j];H=A.sub(G,ids)
            if j in L['ideal_deletions']:
                assert A.inertia(H)==[n-2,0,1];A.kernel1(H)
            else:invs[j]=(ids,A.inverse(H))
        for W in combinations(range(n),n-7):
            B=[i for i in range(n) if i not in W]
            I=[i for i in L['ideal_deletions'] if i in B]
            row=dict(diagram=L['index'],order=n,internal_whites=list(W),black=B,
                black_affine_deletions=I)
            if len(I)>=2:
                row.update(method='external_whites_orthogonal_to_L',possible_ideal_counts=[len(I)])
            elif I:
                a=I[0];thresholds=[]
                for j in B:
                    if j==a:continue
                    ids,J=invs[j];t=(1/J[ids.index(a)][ids.index(a)]).rational()
                    assert t>0;thresholds.append([j,str(t)])
                low=min(F(t) for j,t in thresholds);ties=[j for j,t in thresholds if F(t)==low]
                for j in B:
                    if j==a:continue
                    ids,J=invs[j];H=A.sub(G,ids);H[ids.index(a)][ids.index(a)]-=Q(low)
                    assert A.inertia(H)==([n-2,0,1] if j in ties else [n-1,0,0])
                row.update(method='rank_one_diagonal_decrease',decreased_black=a,
                    thresholds=thresholds,minimum=str(low),ties=ties,
                    possible_ideal_counts=sorted(set([1,1+len(ties)])))
            else:
                h=L['ideal_deletions'][0];key=(L['index'],h)
                if key not in profile_cache:profile_cache[key]=affine_profiles(L,h,diagrams)
                ids,profiles=profile_cache[key]
                vectors=[[Q() if m==2 else Q(1) if m==3 else Q(0,1) for m in p['labels']] for p in profiles]
                mandatory=[2 if G[h][i]==Q() else 3 if G[h][i]==Q(-1) else 4 for i in ids]
                mi=next(i for i,p in enumerate(profiles) if p['labels']==mandatory)
                costs={}
                for b in B:
                    remaining=[i for i in ids if i!=b];J=A.inverse(A.sub(G,remaining))
                    loc=[ids.index(i) for i in remaining];r=len(loc)
                    costs[b]=[[sum((vectors[u][loc[i]]*J[i][j]*vectors[v][loc[j]]
                        for i in range(r) for j in range(r)),Q())
                        for v in range(len(vectors))] for u in range(len(vectors))]
                pairs=[];compatible=[]
                for p in range(len(profiles)):
                    ds=[]
                    for b,C in costs.items():
                        det=(Q(2)-C[mi][mi])*(Q(2)-C[p][p])-C[mi][p]*C[mi][p]
                        ds.append([b,det.data()])
                    bad=[b for b,C in costs.items() if
                        ((Q(2)-C[mi][mi])*(Q(2)-C[p][p])-C[mi][p]*C[mi][p]).sign()<0]
                    if not bad:compatible.append(p)
                    pairs.append(dict(other_profile=p,right_native_determinants=ds,excluded_at=bad))
                row.update(method='no_black_affine_deletion',black_inertia=A.inertia(A.sub(G,B)),
                    selected_affine_deletion=h,affine_indices=ids,profiles=profiles,
                    mandatory_profile=mi,profile_cost_matrices=[dict(black=b,matrix=A.matrix_data(C)) for b,C in costs.items()],
                    mandatory_pair_checks=pairs,possible_additional_profile_indices_at_native_right=compatible,
                    scope='Complete direct-white alphabet under minimum-order choice ofL; pair survivors are necessary conditions only. Arbitrary indirect whites remain an analytic obligation.')
                hard.append(row)
            if I:
                allowed=[i for i in row['possible_ideal_counts'] if 7-i>=2]
                assert all(i<=3 for i in allowed)
                row['ideal_counts_with_at_least_two_finite_vertices']=allowed
            rows.append(row)
    counts=Counter(r['method'] for r in rows)
    assert len(rows)==539 and counts==Counter(external_whites_orthogonal_to_L=175,
        rank_one_diagonal_decrease=290,no_black_affine_deletion=74)
    assert sum(r['black_inertia']==[6,0,1] for r in hard)==4
    assert sum(r['black_inertia']==[7,0,0] for r in hard)==70
    paths=[f'scripts/{TAG}.py',f'scripts/{CENSUS}.py',f'results/{CENSUS}.json',
        'scripts/fresh_finite_edge_cusp_tradeoff12_20261017.py',
        'results/fresh_finite_edge_cusp_tradeoff10_20261017.json',
        'results/fresh_finite_edge_cusp_tradeoff10_20261017.md',
        'results/fresh_multi_ideal_triangle_scope_20261017.md']
    return dict(status='PASS_EXACT_SIX_SIMPLEX_539_CASE_PARTITION_AND_DIRECT_WHITE_FRONTIER',
        scope='465 cases proved at most3 ideals for at least2finite vertices;74hardcases remain. Complete direct-white alphabets require minimum-orderL. No bound on arbitrary supporting-white rank, fullsix-simplextheorem or dimensionclaim.',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        cases=rows,counts=dict(counts),proved_easy_cases=465,
        hard_cases=74,hard_black_inertia_counts={'affine_rank6':4,'positive_definite_rank7':70},
        easy_ideal_count_distributions=dict(Counter(','.join(map(str,r['possible_ideal_counts'])) for r in rows if 'possible_ideal_counts'in r)),
        minimum_order_lemma='Choose L with least order among all indefinite principal subdiagrams of the actual star. For affineF=L-h and any outsidewhiteu touchingF, F+u is indefinite and everyproper subset has smallerorder, hencePSD. Therefore its direct contact profile occurs in the frozen same-order quasiLanner catalogue. Whites rightF may be arbitrarily numerous and connected; eliminating them contributes entrywise nonnegative PSD correction on the directwhite block.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(out.read_text())==d
    else:out.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],d['counts'],'saved_replay=',a.verify,flush=True)
