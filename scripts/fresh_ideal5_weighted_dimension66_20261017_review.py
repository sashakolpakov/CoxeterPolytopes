#!/usr/bin/env python3
"""Independent dimension66 scalar reconstruction and local-rule interfaces.

Standard library only: no project arithmetic imports. Exact forest upper
bounds are inherited from the separately checked finite-budget certificate.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = 'fresh_ideal5_weighted_dimension66_20261017'
TAG = PRIMARY+'_review'
BASE = 'fresh_ideal5_weighted_dimension68_20261017'
LOCAL = 'fresh_finite_d17_long_pair_support_20261017'
OVER = 'fresh_distance8_high_overlap_20261017'
FOREST = 'fresh_rank17_reduced_forest_frontier_20261017'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def cost(r):
    special = {11:F(65,12),12:F(75,13),13:F(43,7),14:F(32,5),15:F(107,16),16:F(120,17)}
    return F(r,2) if r <= 10 else special[r] if r <= 16 else F(163,23)


def anchors(total, ranks=(10,11,13,14,15,16,17)):
    if not ranks:
        if total == 0:
            yield ()
        return
    r, *rest = ranks
    for multiplicity in range(total//r+1):
        for tail in anchors(total-r*multiplicity,tuple(rest)):
            yield (r,)*multiplicity+tail


def shapes(d):
    out = set()
    for low in range(10):
        for twelve in (0,1):
            for high in (0,*range(18,d+1)):
                remainder = d-low-12*twelve-high
                if remainder < 0:
                    continue
                optional = ((low,) if low else ())+((12,) if twelve else ())+((high,) if high else ())
                for repeated in anchors(remainder):
                    p = tuple(sorted(optional+repeated))
                    assert p not in out
                    out.add(p)
    return sorted(out)


def scalar(ranks):
    h = [1]
    for r in ranks:
        nxt = []; window = 0
        for j in range(len(h)+r):
            if j < len(h):
                window += h[j]
            if 0 <= j-r-1 < len(h):
                window -= h[j-r-1]
            nxt.append(window)
        h = nxt
    d = sum(ranks); v = prod(r+1 for r in ranks); half = d//2
    assert len(h) == d+1 and h == h[::-1] and sum(h) == v
    g = sum(h[j]*min(j,d-j) for j in range(d+1))
    shifts = {0:1}
    for r in ranks:
        nxt = dict(shifts)
        for j,a in shifts.items():
            nxt[j+r+1] = nxt.get(j+r+1,0)-a
        shifts = {j:a for j,a in nxt.items() if a}
    s = len(ranks)
    cdf_g = half*v-2*sum(a*(comb(half-j+s,s+1) if half-j+s >= s+1 else 0) for j,a in shifts.items())
    assert g == cdf_g
    curvature = F(5*g-d*v,2*v)
    absolute = F(sum(h[j]*abs(2*j-d) for j in range(d+1)),2*v)
    assert curvature == F(3*d,4)-F(5,2)*absolute
    payment = sum(map(cost,ranks),F())
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=v,h_polynomial=h,G=g,
        normalized_curvature=str(curvature),normalized_cost=str(payment),margin=str(curvature-payment))


def arithmetic():
    rows = []; minima = []
    for n,expected_count,expected_minimum,expected_rank,finite in (
        (67,1761,'851/165888',[11]*6,504),
        (68,1906,'411499/3234816',[11]*5+[12],512),
    ):
        group = [scalar(p) for p in shapes(n-1)]
        assert len(group) == expected_count and all(F(r['margin']) > 0 for r in group)
        minimum = min(group,key=lambda r:F(r['margin']))
        assert (minimum['margin'],minimum['ranks']) == (expected_minimum,expected_rank)
        beta = F(n*(n-6),8) if n%2==0 else F((n-1)*(n-5),8)
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=minimum['margin'],
            minimizing_ranks=minimum['ranks'],finite_maximum=str(finite),finite_curvature_coefficient=str(beta),
            finite_slack=str(beta-finite),finite_maximum_witness=[f'D{n}'],strict_global_slack=5*(n-1)))
        rows.extend(group)
    gain = cost(11)+cost(13)-2*cost(12)
    assert gain == F(23,1092) > 0
    assert all(cost(r)==F(r,2) for r in range(1,11))
    assert all(cost(r)==F(163,23) for r in range(17,69))
    boundary = scalar((13,)*5)
    assert boundary['margin'] == '-1565/10976'
    return dict(reduced_profiles=rows,dimension_minima=minima,boundary_at_n66=boundary,
        duplicate_rank12_cost_increase=str(gain),
        formulas=['sliding-window product of uniform polynomials',
            'G=floor(d/2)*V-2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1)',
            'kappa/V=3d/4-(5/2)*E(abs(sum U_i-d/2))'])


def graph(kind,n):
    g = [[0]*n for _ in range(n)]
    edges = [(i,i+1,3) for i in range(n-1)]
    if kind == 'B':
        edges[-1] = (n-2,n-1,4)
    if kind == 'D':
        edges = [(i,i+1,3) for i in range(n-3)]+[(n-3,n-2,3),(n-3,n-1,3)]
    for i,j,label in edges:
        g[i][j] = g[j][i] = label
    return g


def distances(g):
    out = []
    for start in range(len(g)):
        d = {start:0}; todo = [start]
        for i in todo:
            for j,bond in enumerate(g[i]):
                if bond and j not in d:
                    d[j] = d[i]+1; todo.append(j)
        assert len(d) == len(g)
        out.append(d)
    return out


def basic(g):
    dist = distances(g)
    leaves = {i for i,row in enumerate(g) if sum(bool(x) for x in row)==1}
    return {p for p in combinations(range(len(g)),2) if dist[p[0]][p[1]] <= 8 or
            (dist[p[0]][p[1]]==9 and set(p)&leaves)}


def native_pattern(gram,ids):
    g = [[F(x) for x in row] for row in gram]
    p = [[0]*len(ids) for _ in ids]
    for i,j in combinations(range(len(ids)),2):
        a,b=ids[i],ids[j]
        if g[a][b]:
            assert g[a][b]<0
            ratio = 4*g[a][b]**2/(g[a][a]*g[b][b])
            if ratio not in (1,2):
                return None
            p[i][j]=p[j][i]=3 if ratio==1 else 4
    return p


def induced_maps(pattern,target):
    if pattern is None:
        return set()
    n,m=len(pattern),len(target)
    degrees=[sum(bool(x) for x in row) for row in pattern]
    start=max(range(n),key=lambda i:(degrees[i],i));order=[start];seen={start}
    for i in order:
        for j in sorted(range(n),key=lambda j:(-degrees[j],j)):
            if pattern[i][j] and j not in seen:
                seen.add(j);order.append(j)
    assert len(order)==n
    host_degrees=[sum(bool(x) for x in row) for row in target]
    maps=set();current=[-1]*n
    def visit(k,available):
        if k==n:
            maps.add(tuple(current));return
        i=order[k]
        for v in sorted(available):
            if host_degrees[v]<degrees[i]:
                continue
            if all(pattern[i][j]==target[v][current[j]] for j in order[:k]):
                current[i]=v;visit(k+1,available-{v});current[i]=-1
    visit(0,set(range(m)))
    return maps


def local_interfaces(local,over):
    assert local['new_payment_theorem'] and local['cross_component_coupons_unchanged'] and local['ideal_weights_unchanged']
    eligible=[c for c in over['isolation']['cases'] if c.get('component_order_bound')==17]
    assert len(eligible)==140 and len(over['isolation']['cases'])==2614
    out=[];removed_by_type={}
    for kind,expected_cost,expected_maps in (('A',114,48),('B',105,18),('D',127,168)):
        saved=next(r for r in local['records'] if r['kind']==kind)
        target=graph(kind,17);assert target==saved['graph']
        base=basic(target);assert base==set(map(tuple,saved['basic_pairs']))
        kept=set(map(tuple,saved['retained_pairs']));removed=set(map(tuple,saved['removed_pairs']))
        assert len(kept)==saved['new_cost']==expected_cost and base<=kept
        assert kept.isdisjoint(removed) and kept|removed==set(combinations(range(17),2))
        autos=induced_maps(target,target)
        assert autos==set(map(tuple,saved['automorphisms']))
        assert all(tuple(sorted(m[i] for i in p)) in removed for m in autos for p in removed)
        cases={r['index']:r for r in saved['selected_high_contexts']}
        assert set(cases)=={r['index'] for r in eligible}
        pair_union=set();map_count=0;nonempty=0
        for c in eligible:
            ids=c['marked_component'];marks=c['preferred_corner']
            pattern=native_pattern(over['isolation']['cores'][c['core_key']]['gram'],ids)
            maps=induced_maps(pattern,target)
            assert maps=={tuple(r['images']) for r in cases[c['index']]['maps']}
            images={tuple(sorted(m[ids.index(i)] for i in marks)) for m in maps}
            assert images<=kept
            pair_union|=images;map_count+=len(maps);nonempty+=bool(maps)
        assert map_count==expected_maps and pair_union==set(map(tuple,saved['selected_high_pairs']))
        if kind=='A':
            prescribed=base|pair_union|{(0,13),(3,16)}
        elif kind=='B':
            prescribed=base|pair_union
        else:
            deleted={(0,10),(0,11)}|{(1,j) for j in range(10,17)}
            prescribed=set(combinations(range(17),2))-deleted
        assert kept==prescribed
        removed_by_type[kind+'17']=removed
        out.append(dict(kind=kind,basic_pairs=len(base),retained_cost=len(kept),removed_pairs=sorted(removed),
                        selected_contexts=nonempty,induced_maps=map_count,marked_images=sorted(pair_union)))
    # Recompute every saved simultaneous total. The geometric sealing and
    # original core inertia are inherited from the frozen local theorem.
    totals=[]
    for row in local['simultaneous_literal_totals']:
        total=F()
        for c in row['corners']:
            expected_debit=0
            for match in c['whole_rank17_matches']:
                images=set(map(tuple,match['marked_images']))
                removed=removed_by_type[match['type']]
                assert images<=removed or images.isdisjoint(removed)
                expected_debit+=images<=removed
            assert expected_debit==c['debit']
            retained=F(c['old_weight'])-c['debit']
            assert retained==F(c['retained_weight']) and retained>=0
            total+=retained
        assert total==F(row['retained_total']) and total>=1
        totals.append(dict(source=row['source'],case=row['case'],retained_total=str(total)))
    assert len(totals)==174
    return dict(records=out,literal_totals=totals,induced_maps=sum(r['induced_maps'] for r in out),
                original_core_census_and_sealing='Inherited frozen geometric inputs; not rerun.')


def main(prepare=False,compare_only=False,verify=False):
    result=arithmetic()
    if prepare:
        p=Path('/private/tmp/dimension66_independent_scalar_20261017.json')
        p.write_text(json.dumps(result,separators=(',',':'))+'\n')
        print('PASS_INDEPENDENT_DIMENSION66_PREPARATION',result['dimension_minima'],flush=True)
        return
    # Completed after the primary and finite-budget reviewer are frozen.
    stems=(PRIMARY,BASE,BASE+'_review',LOCAL,OVER,FOREST,FOREST+'_review')
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    primary=docs[PRIMARY]
    assert primary['dimension_bound']==66 and primary['inherited_dimension_bound']==68
    assert primary['excluded_dimensions']==[67,68]
    assert primary['unrestricted_finite_volume'] and not primary['arithmetic_hypothesis'] and primary['cusp_restriction'] is None
    by={tuple(r['ranks']):r for r in primary['reduced_profiles']}
    assert len(by)==len(primary['reduced_profiles'])==len(result['reduced_profiles'])==3667
    assert by=={tuple(r['ranks']):r for r in result['reduced_profiles']}
    assert primary['dimension_minima']==result['dimension_minima']
    assert primary['boundary_at_n66']['profile']==result['boundary_at_n66']
    result['local_interfaces']=local_interfaces(docs[LOCAL],docs[OVER])
    finite=docs[FOREST]['variants'][1]
    assert finite['assumed_A17_base']==114 and finite['assumed_D17_base']==127
    assert finite['assumed_B17_base_upper_bound']==127>=105
    assert finite['dimension_maxima']==docs[FOREST+'_review']['variants'][1]['dimension_maxima']
    witness_checks=[]
    for n,wanted in ((67,504),(68,512)):
        row=finite['dimension_maxima'][n]
        assert row['doubled_cost']==2*wanted and row['witness_components']==[f'D{n}']
        count=len(basic(graph('D',n)))
        assert count==wanted==8*n-32
        witness_checks.append(dict(dimension=n,type=f'D{n}',direct_pair_cost=count))
    result['finite_witnesses']=witness_checks
    assert len(basic(graph('D',66)))==496
    assert finite['dimension_maxima'][66]['doubled_cost']==992
    assert primary['boundary_at_n66']==dict(profile=result['boundary_at_n66'],finite_maximum='496',
        finite_curvature_coefficient='495',finite_slack='-1',finite_maximum_witness=['D66'],new_bound65_asserted=False)
    assert primary['ideal_rank_envelope']==docs[BASE]['ideal_rank_envelope']
    assert primary['small_order_costs']==docs[BASE]['small_order_costs']
    assert primary['concentration']==docs[BASE]['concentration']
    assert primary['ideal_weights_unchanged'] and primary['cross_component_coupons_unchanged']
    assert docs[BASE]['dimension_bound']==docs[BASE+'_review']['dimension_bound']==68
    for r in result['local_interfaces']['records']:
        saved=next(s for s in primary['finite_rank17_interfaces'] if s['kind']==r['kind'])
        assert (saved['rank'],saved['new_cost'],saved['selected_contexts'],saved['selected_maps']) == (
            17,r['retained_cost'],r['selected_contexts'],r['induced_maps'])
        assert list(map(tuple,saved['removed_pairs']))==r['removed_pairs']
    paths=[Path(__file__)]
    for stem,doc in docs.items():
        assert doc['status'].startswith('PASS'),stem
        for rel,h in doc['sources'].items():
            assert digest(ROOT/rel)==h,(stem,rel)
        paths.extend(ROOT/f'{folder}/{stem}.{ext}' for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    if compare_only:
        print('PASS_INDEPENDENT_DIMENSION66_PRIMARY_COMPARISON',len(result['reduced_profiles']),
              'profiles;234 independently enumerated finite maps;174 saved totals',flush=True)
        return
    paths.extend(ROOT/f'results/{s}.md' for s in (
        'fresh_cusp_link_shape_curvature_20261017','fresh_cusp_factor_payment_framework_20261017',
        'fresh_distance14_unrestricted_dimension156_20261017'))
    result.update(status='PASS_INDEPENDENT_DIMENSION66_SCALARS_AND_INTERFACES',dimension_bound=66,
        sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
        scope='No project arithmetic imports: complete3667 scalar profiles in two dimensions, three curvature formulas, all234 finite-rank17 selected maps independently enumerated,174 saved literal totals, canonical retained pair rules and direct D67/D68 attaining costs. Exact forest upper certificates, original all-label geometry and sealing remain inherited from separately reviewed packages.')
    result=json.loads(json.dumps(result));p=ROOT/f'results/{TAG}.json'
    if verify:
        assert result==json.loads(p.read_text())
    else:
        p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['dimension_minima'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--compare-only',action='store_true');p.add_argument('--verify',action='store_true')
    args=p.parse_args();assert sum((args.prepare,args.compare_only,args.verify))<=1
    main(args.prepare,args.compare_only,args.verify)
