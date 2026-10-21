#!/usr/bin/env python3
"""Independent Fraction scalar and literal-interface review of dimension69.

No producer arithmetic imports. The scalar reconstruction independently uses
optional-first partitions, sliding convolution, cumulative bounded compositions
and the absolute-deviation identity. The unchanged finite budget is inherited.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
import json
from math import comb,prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_ideal5_weighted_dimension69_20261017'
TAG=PRIMARY+'_review'
BASE='fresh_ideal5_weighted_dimension70_20261017'
LOCAL='fresh_b14_b15_singleton_completion_transfer_20261017'
OVER='fresh_distance8_high_overlap_20261017'
PRIOR_B15='fresh_order15_affine_b_pair_reduction_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
FINITE='fresh_finite8_d12_b10_coupon_budget_20261017'
def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def choose(n,k):return comb(n,k)if n>=k>=0 else 0

def scalar(ranks):
    # Sliding-window multiplication, rather than the producer's nested sums.
    h=[1]
    for r in ranks:
        out=[];window=0
        for j in range(len(h)+r):
            if j<len(h):window+=h[j]
            if 0<=j-r-1<len(h):window-=h[j-r-1]
            out.append(window)
        h=out
    d=sum(ranks);V=prod(r+1 for r in ranks);half=d//2
    assert len(h)==d+1 and h==h[::-1]and sum(h)==V
    g=sum(a*min(j,d-j)for j,a in enumerate(h))
    # Numerator of product(1-z^(r+1)); repeated shifts are combined exactly.
    shifts={0:1}
    for r in ranks:
        nxt=dict(shifts)
        for t,a in shifts.items():nxt[t+r+1]=nxt.get(t+r+1,0)-a
        shifts={t:a for t,a in nxt.items()if a}
    s=len(ranks)
    # G=sum_{t=1}^{floor(d/2)}(V-2*CDF(t-1)); hockey-stick summation.
    gcdf=half*V-2*sum(a*choose(half-t+s,s+1)for t,a in shifts.items())
    assert g==gcdf
    curvature=F(5*g-d*V,2*V)
    absolute=F(sum(a*abs(2*j-d)for j,a in enumerate(h)),2*V)
    assert curvature==F(3*d,4)-F(5,2)*absolute
    cost=sum((budget(r)for r in ranks),F())
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=V,h_polynomial=h,G=g,
        independent_cumulative_G=gcdf,normalized_curvature=str(curvature),
        normalized_cost=str(cost),margin=str(curvature-cost))

def budget(r):
    return F(r,2)if r<=11 else F(75,13)if r==12 else F(43,7)if r==13 else F(32,5)if r==14 else F(107,16)if r==15 else F(120,17)if r==16 else F(163,23)

def affine_feature_maxima():
    # Independent spine counts at order m>=18. Between consecutive breakpoints
    # cost/m is a+b/m, so endpoints and the limiting value5 exhaust the maximum.
    specs={'B':(-6,(8,9,10,11),13,F(160,23)),
           'C':(-15,(6,7,8,9),14,F(75,11)),
           'D':(3,(10,11,12,13),12,F(163,23))}
    rows=[]
    for kind,(constant,shifts,cap,expected)in specs.items():
        breakpoints={s+cap for s in shifts}
        points=sorted({18}|{x+e for x in breakpoints for e in(-1,0,1)if x+e>=18})
        checks=[[m,5*m+constant+sum(min(m-s,cap)for s in shifts)]for m in points]
        maximum=max([F(5)]+[F(c,m)for m,c in checks]);assert maximum==expected
        cuts=sorted({18}|breakpoints);branches=[]
        for i,lo in enumerate(cuts):
            hi=cuts[i+1]if i+1<len(cuts)else None
            sample=F(lo)+F(1,2)
            active=[sample-s<cap for s in shifts]
            branches.append(dict(kind=kind,lower_order=lo,upper_order=hi,
                slope=5+sum(active),intercept=constant+sum(-s if a else cap for s,a in zip(shifts,active))))
        rows.append(dict(kind=kind,constant=constant,shifts=list(shifts),window_cap=cap,
            exact_breakpoint_checks=checks,linear_branches=branches,
            limiting_ratio='5',maximum_ratio=str(maximum)))
    return rows

def anchors(total, allowed=(11,13,14,15,16,17)):
    if len(allowed)==1:
        if total%allowed[0]==0:
            yield [allowed[0]]*(total//allowed[0])
        return
    r=allowed[0]
    for count in range(total//r+1):
        for tail in anchors(total-count*r,allowed[1:]):
            yield [r]*count+tail


def arithmetic(compare_primary=True):
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())if compare_primary else None
    profiles=[];minima=[]
    for d in (69,):
        shapes=set()
        for high in(0,*range(18,d+1)):
            for low in range(11):
                for twelve in(0,1):
                    left=d-high-low-12*twelve
                    if left<0:continue
                    exceptional=([high]if high else[])+([low]if low else[])+[12]*twelve
                    for rest in anchors(left):shapes.add(tuple(sorted(exceptional+rest)))
        group=[scalar(rs)for rs in sorted(shapes)]
        assert all(F(r['margin'])>0 for r in group)
        profiles+=group;n=d+1
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta==8*n
        m=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=m['margin'],
            minimizing_ranks=m['ranks'],finite_slack=str(beta-8*n),strict_global_slack=5*(n-1)))
    assert len(profiles)==1339 and minima==[dict(ambient_dimension=70,profiles=1339,minimum_margin='1819/51840',minimizing_ranks=[11,11,11,11,11,14],finite_slack='0',strict_global_slack=345)]
    assert budget(11)+budget(13)-2*budget(12)==F(19,182)>0
    boundary=scalar((11,11,11,11,11,13))
    assert boundary['margin']=='-2029/13608'
    if primary is not None:
        by={tuple(r['ranks']):r for r in primary['reduced_profiles']}
        assert set(by)=={tuple(r['ranks'])for r in profiles}and len(by)==len(primary['reduced_profiles'])==len(profiles)
        for r in profiles:
            assert {k:v for k,v in r.items()if k!='independent_cumulative_G'}==by[tuple(r['ranks'])]
        assert minima==primary['dimension_minima']
        assert {k:v for k,v in boundary.items()if k!='independent_cumulative_G'}==primary['boundary_at_n69']
    return dict(profiles=profiles,dimension_minima=minima,boundary=boundary,
        minimum_margin=str(min(F(r['margin'])for r in profiles)),
        independent_affine_feature_maxima=affine_feature_maxima(),
        concentration=dict(duplicate_rank12_cost_increase=str(budget(11)+budget(13)-2*budget(12)),
            terminal_repeated_ranks=[11,13,14,15,16,17],optional_rank12=True,
            at_most_one_low_rank=list(range(1,11)),at_most_one_high_rank_minimum=18,
            termination='Each nontrivial concentration increases the integer sum of squared parts, bounded by d^2.'),
        formulas=['G=sum_j h_j min(j,d-j)',
            'G=floor(d/2)*V-2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1)',
            'kappa/V=3d/4-(5/2)*E(abs(sum U_i-d/2))'])


def freeze(g):return tuple(tuple(F(x)for x in row)for row in g)

@lru_cache(None)
def inverse(g):
    n=len(g);a=[list(row)+[F(i==j)for j in range(n)]for i,row in enumerate(g)]
    for j in range(n):
        k=next(k for k in range(j,n)if a[k][j]);a[j],a[k]=a[k],a[j]
        z=a[j][j];a[j]=[x/z for x in a[j]]
        for k in range(n):
            if k!=j and a[k][j]:
                z=a[k][j];a[k]=[x-z*y for x,y in zip(a[k],a[j])]
    return tuple(tuple(row[n:])for row in a)

@lru_cache(None)
def inertia(g):
    # Exact congruence elimination, including a hyperbolic 2x2 pivot when
    # every remaining diagonal is zero. No producer matrix code is imported.
    a=[list(row)for row in g];positive=negative=zero=0
    while a:
        n=len(a);pivot=next((i for i in range(n)if a[i][i]),None)
        if pivot is not None:
            q=a[pivot][pivot];positive+=q>0;negative+=q<0
            keep=[i for i in range(n)if i!=pivot]
            a=[[a[i][j]-a[i][pivot]*a[pivot][j]/q for j in keep]for i in keep]
            continue
        pair=next(((i,j)for i in range(n)for j in range(i+1,n)if a[i][j]),None)
        if pair is None:zero+=n;break
        p,q=pair;t=a[p][q];positive+=1;negative+=1
        keep=[i for i in range(n)if i not in pair]
        a=[[a[i][j]-(a[i][p]*a[q][j]+a[i][q]*a[p][j])/t for j in keep]for i in keep]
    return positive,negative,zero

def sub(g,ids):return tuple(tuple(g[i][j]for j in ids)for i in ids)

def components(g,ids):
    left=set(ids);out=[]
    while left:
        todo=[min(left)];seen=set(todo)
        for i in todo:
            for j in sorted(left-seen):
                if i!=j and g[i][j]:seen.add(j);todo.append(j)
        left-=seen;out.append(sorted(seen))
    return out


def graph_distances(g,starts):
    d={i:0 for i in starts}; todo=list(starts)
    for i in todo:
        for j,x in enumerate(g[i]):
            if i!=j and x and j not in d:
                d[j]=d[i]+1; todo.append(j)
    assert len(d)==len(g)
    return d


def b_target(m):
    g=[[0]*m for _ in range(m)]
    for i in range(m-3):
        g[i][i+1]=g[i+1][i]=4 if i==m-4 else 3
    for i in (m-2,m-1): g[0][i]=g[i][0]=3
    return g


def local_interfaces(local,small,over,prior_b15):
    records=[]; checked_cores=set(); total_maps=0; target_contexts=0; unique_cases=set()
    assert local['added_finite_coupons']==[]
    assert local['triangle_minimum']==2 and local['quadrilateral_minimum']==1
    for rec in local['records']:
        m=rec['order']; assert m in (14,15)
        drop=(0,10) if m==14 else (1,11)
        target=b_target(m); assert target==rec['target']
        ds=[graph_distances(target,[i]) for i in range(m)]
        identity=list(range(m)); flip=identity.copy(); flip[-2],flip[-1]=flip[-1],flip[-2]
        autos=[identity,flip]
        assert set(map(tuple,autos))==set(map(tuple,rec['automorphisms']))
        for a in autos:
            assert all(target[i][j]==target[a[i]][a[j]] for i in range(m) for j in range(m))
            assert tuple(sorted(a[i] for i in drop))==drop
        old=next(r for r in small['records'] if r['order']==m and r['kind']=='B')
        assert old['target']==target
        paid=set(map(tuple,old['retained_paid_pairs'] if m==14 else prior_b15['new_paid_pairs']))
        assert len(paid)==rec['inherited_cost']==(87 if m==14 else 97)
        reduced=paid-{drop}; assert len(reduced)==rec['new_cost']==(86 if m==14 else 96)
        assert set(map(tuple,rec['new_paid_pairs']))==reduced and rec['removed_pair']==list(drop)
        features={0,m-4,m-3}
        short={p for p in combinations(range(m),2) if ds[p[0]][p[1]]<=8}
        triangles={p for p in combinations(range(m),2) if ds[p[0]][p[1]]<=5 or
                   (ds[p[0]][p[1]]<=8 and any(ds[i][f]<=5 for i in p for f in features))}
        basic={p for p in combinations(range(m),2) if ds[p[0]][p[1]]<=5 or
               (ds[p[0]][p[1]]<=9 and any(ds[i][f]<=5 for i in p for f in features))}
        assert triangles==set(map(tuple,rec['triangle_pairs_retained'])) and triangles<=reduced
        assert short==set(map(tuple,rec['short8_pairs_retained'])) and short<=reduced
        assert basic<=reduced and ds[drop[0]][drop[1]]==10
        affected={}; all_saved_maps=0
        for r in old['case_embeddings']:
            old_case=over['isolation']['cases'][r['case']]
            positions=[old_case['marked_component'].index(i) for i in old_case['preferred_corner']]
            all_images={tuple(sorted(mp[i] for i in positions)) for mp in r['embeddings']}
            assert all_images==set(map(tuple,r['marked_images']))
            all_saved_maps+=len(r['embeddings'])
            if drop in all_images: affected[r['case']]=r
        assert sorted(affected)==rec['affected_cases']==[c['index'] for c in rec['cases']]
        counts=Counter(); maps=0; literals=[]
        for case in rec['cases']:
            index=case['index']; unique_cases.add(index); target_contexts+=1
            saved=over['isolation']['cases'][index]
            for key in ('core_key','method','left','right','match','preferred_corner'):
                assert case[key]==saved[key]
            gram=freeze(over['isolation']['cores'][case['core_key']]['gram'])
            assert gram==freeze(case['core_gram']); n=len(gram); checked_cores.add(gram)
            assert all(gram[i][i]==2 for i in range(n))
            assert all(gram[i][j] in (-1,0) for i in range(n) for j in range(n) if i!=j)
            J=inverse(gram)
            assert inertia(gram)==(n-1,1,0) and J==freeze(case['core_inverse'])
            u=case['left'][1]; v=[saved['right_copy_map'][i] for i in case['right'][1]]
            black=set(u+v); white=sorted(set(range(n))-black); marks=case['preferred_corner']
            assert len(black)==4 and inertia(sub(gram,white))==(n-4,0,0)
            for b in black: assert inertia(sub(gram,sorted(set(white)|{b})))==(n-3,0,0)
            corner=sorted(set(white)|set(marks)); omitted=sorted(black-set(marks))
            assert inertia(sub(gram,corner))==(n-2,0,0)
            blocks=components(gram,corner)
            base=next(b for b in blocks if marks[0] in b)
            assert marks[1] in base and base==case['base_nodes']==saved['marked_component']
            allowed=set(omitted)
            for block in blocks:
                ib=inverse(sub(gram,block))
                assert all(x>=0 for row in ib for x in row)
                allowed.update(i for t,i in enumerate(block) if ib[t][t]<=2)
            assert sorted(allowed)==case['contact_domain']==saved['possible_contacts']
            ends=case['end_sets']; proof=case['end_projection_proof']
            assert ends==proof['end_sets'] and sorted(allowed)==proof['allowed_contacts']
            assert set(ends[0]).isdisjoint(ends[1]) and set(ends[0]+ends[1])==allowed
            du=graph_distances(gram,u); dv=graph_distances(gram,v)
            shortcut=max((du[i] for i in ends[0]),default=0)+max((dv[i] for i in ends[1]),default=0)+2
            assert shortcut==proof['both_ends_shortcut_bound']<=8
            assert all(J[i][j]>0 for end in ends for i in end for j in end)
            assert [sub(J,end) for end in ends]==[freeze(g) for g in proof['inverse_end_blocks']]
            others=sorted(set(corner)-set(base))
            assert others==proof['other_old_corner_nodes'] and omitted==proof['old_nodes_absent_from_corner']
            assert all(gram[i][j]==0 for i in base for j in others)
            expected={tuple(mp) for mp in affected[index]['embeddings']
                      if tuple(sorted(mp[base.index(i)] for i in marks))==drop}
            seen=set()
            for row in case['maps']:
                mp=tuple(row['map']); assert mp in expected and mp not in seen; seen.add(mp)
                assert len(mp)==len(set(mp))==len(base)
                assert all((int(2-gram[a][b]) if a!=b and gram[a][b] else 0)==target[mp[i]][mp[j]]
                           for i,a in enumerate(base) for j,b in enumerate(base))
                missing=sorted(set(range(m))-set(mp)); assert missing==row['missing_target_nodes']
                reverse={mp[t]:i for t,i in enumerate(base)}
                terminal=m-3; neighbor=m-4
                assert terminal in missing and neighbor in reverse and target[terminal][neighbor]==4
                i=reverse[neighbor]; assert i in marks and i in allowed
                lower=2*J[i][i]; assert lower>2
                witness=dict(new_target_node=terminal,old_target_node=neighbor,core_node=i,
                             coxeter_label=4,norm_two_contact_square='2',
                             full_core_inverse_diagonal=str(J[i][i]),forced_projection_lower_bound=str(lower))
                assert witness==row['witness'] and witness in row['forced_entries']
                assert row['outcome']=='projection_exclusion'
                counts[(str(J[i][i]),str(lower))]+=1; maps+=1
                literals.append(dict(case=index,map=list(mp),terminal_target=terminal,core_neighbor=i,
                                     inverse_diagonal=str(J[i][i]),norm_two_contact_square='2',projection_lower_bound=str(lower)))
            assert seen==expected
        assert maps==rec['affected_map_count']==(46 if m==14 else 18)
        assert rec['map_outcomes']=={'projection_exclusion':maps}
        total_maps+=maps
        records.append(dict(order=m,affine_type=rec['affine_type'],removed_pair=list(drop),
                            old_cost=len(paid),new_cost=len(reduced),affected_cases=sorted(affected),affected_maps=maps,
                            minimum_forced_projection=str(min(F(b) for a,b in counts)),
                            full_target_saved_map_count=all_saved_maps,
                            projection_histogram=[dict(inverse_diagonal=a,forced_projection=b,maps=k) for (a,b),k in sorted(counts.items())],
                            triangle_units=len(triangles),short8_units=len(short),basic_units=len(basic),
                            literal_witnesses=literals,new_finite_cost='0'))
    assert total_maps==64 and target_contexts==34 and len(unique_cases)==29 and len(checked_cores)==13
    return dict(records=records,total_maps=total_maps,target_contexts=target_contexts,
                distinct_contexts=len(unique_cases),distinct_lorentz_core_grams=len(checked_cores),
                scope='Independent norm-two Gram inverse/inertia, actual preferred-corner contact domain, end block/distance proof, complete saved affected-map comparison and all64 forced doubled-terminal bounds. Original exhaustive core/map censuses remain frozen inputs.')


def interfaces(docs,data):
    primary=docs[PRIMARY]; base=docs[BASE]; review=docs[BASE+'_review']; finite=docs[FINITE]
    assert primary['dimension_bound']==69 and primary['inherited_dimension_bound']==70
    assert base['dimension_bound']==review['dimension_bound']==70 and primary['excluded_dimensions']==[70]
    assert primary['unrestricted_finite_volume'] and not primary['arithmetic_hypothesis'] and primary['cusp_restriction'] is None
    assert primary['finite_vertex_budget']==finite['finite_vertex_budget']=='8n'
    assert primary['finite_budget_dimension_maximum']==finite['dimension_maximum']==75
    previous_budget=review['composition_interfaces']['inherited_finite_budget']
    assert previous_budget['previous_independent_endpoints']==13034 and len(finite['dimension_endpoints'])==13034
    local=local_interfaces(docs[LOCAL],docs[SMALL],docs[OVER],docs[PRIOR_B15])
    costs={k:dict(v) for k,v in base['small_order_costs'].items()}
    assert costs['14']=={'A':70,'B':87,'C':84,'D':82}; costs['14']['B']=86
    assert costs['15']=={'A':75,'B':97,'C':95,'D':96}; costs['15']['B']=96
    assert primary['small_order_costs']==costs
    assert all(F(max(row.values()),int(m))==budget(int(m)-1) for m,row in costs.items())
    env=dict(base['ideal_rank_envelope']); env.update({'13':'43/7','14':'32/5'})
    assert primary['ideal_rank_envelope']==env
    for saved,independent in zip(primary['affine_B_interfaces'],local['records']):
        for key in ('order','old_cost','new_cost','removed_pair','affected_cases','affected_maps','minimum_forced_projection','new_finite_cost'):
            assert saved[key]==independent[key]
    assert data['concentration']['terminal_repeated_ranks']==primary['concentration']['repeated_cap_ranks']
    return dict(affine_B=local,small_order_costs=costs,ideal_rank_envelope=env,
                inherited_finite_budget=dict(budget='8n',maximum_dimension=75,previous_independent_endpoints=13034),
                inherited_ceiling=70,newly_excluded_dimension=70,finite_slack='0',strict_global_slack=345,
                scope='Nonnegative finite coefficient and strictly positive global additive constant suffice, including zero cusps. All ideal coefficients are positive. Original complete core/map censuses and unchanged finite-budget enumeration are inherited.')


def main(prepare=False,compare_only=False,verify=False):
    data=arithmetic(not prepare)
    if prepare:
        Path('/private/tmp/dimension69_independent_scalar.json').write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')
        print('PASS_INDEPENDENT_SCALAR_PREPARATION',len(data['profiles']),data['dimension_minima'],'boundary',data['boundary']['margin'],flush=True)
        return
    stems=[PRIMARY,BASE,BASE+'_review',LOCAL,SMALL,OVER,PRIOR_B15,FINITE]
    docs={t:json.loads((ROOT/f'results/{t}.json').read_text()) for t in stems}
    data['composition_interfaces']=interfaces(docs,data)
    if compare_only:
        print('PASS_INDEPENDENT_PRIMARY_COMPARISON',len(data['profiles']),'profiles;64 exact local maps',flush=True)
        return
    paths=[Path(__file__)]
    for t in stems:
        assert docs[t]['status'].startswith('PASS'),t
        for p,h in docs[t]['sources'].items(): assert digest(ROOT/p)==h,(t,p)
        paths.extend(ROOT/f'{folder}/{t}.{ext}' for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    paths.extend(ROOT/f'results/{t}.md' for t in ('fresh_cusp_link_shape_curvature_20261017','fresh_cusp_factor_payment_framework_20261017'))
    data.update(status='PASS_INDEPENDENT_DIMENSION69_SCALARS_AND_COMPOSITION',dimension_bound=69,
                sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
                scope='No project arithmetic imports:1339 scalar profiles and15 all-rank branches;34 target/context occurrences,64 exact saved maps,13 independently inverted Lorentzian cores and all actual preferred-corner contact domains. Original complete censuses and unchanged13034 finite-budget endpoints remain inherited; full mathematical composition reviewed.')
    data=json.loads(json.dumps(data)); path=ROOT/f'results/{TAG}.json'
    if verify: assert data==json.loads(path.read_text())
    else: path.write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')
    print(data['status'],len(data['profiles']),'profiles;64 exact local maps;zero finite slack and global345',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--prepare',action='store_true'); p.add_argument('--compare-only',action='store_true'); p.add_argument('--verify',action='store_true')
    a=p.parse_args(); assert sum((a.prepare,a.compare_only,a.verify))<=1; main(a.prepare,a.compare_only,a.verify)

