#!/usr/bin/env python3
"""Independent Fraction scalar and literal-interface review of dimension70.

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
PRIMARY='fresh_ideal5_weighted_dimension70_20261017'
TAG=PRIMARY+'_review'
BASE='fresh_ideal5_weighted_dimension71_20261017'
LOCAL='fresh_d16_partial_weight_transfer_20261017'
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
    return F(r,2)if r<=11 else F(75,13)if r==12 else F(87,14)if r==13 else F(97,15)if r==14 else F(107,16)if r==15 else F(120,17)if r==16 else F(163,23)

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
    for d in (70,):
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
        assert beta>8*n
        m=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=m['margin'],
            minimizing_ranks=m['ranks'],finite_slack=str(beta-8*n),strict_global_slack=5*(n-1)))
    assert budget(11)+budget(13)-2*budget(12)==F(16,91)>0
    boundary=scalar((11,11,11,11,11,14))
    assert F(boundary['margin'])<0
    if primary is not None:
        by={tuple(r['ranks']):r for r in primary['reduced_profiles']}
        assert set(by)=={tuple(r['ranks'])for r in profiles}and len(by)==len(primary['reduced_profiles'])==len(profiles)
        for r in profiles:
            assert {k:v for k,v in r.items()if k!='independent_cumulative_G'}==by[tuple(r['ranks'])]
        assert minima==primary['dimension_minima']
        assert {k:v for k,v in boundary.items()if k!='independent_cumulative_G'}==primary['boundary_at_n70']
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

def finite_component(g,nodes):
    adj={i:[j for j in nodes if i!=j and g[i][j]]for i in nodes}
    assert all(g[i][i]==2 for i in nodes)
    assert all(g[i][j]==-1 for i in nodes for j in adj[i])
    assert sum(map(len,adj.values()))==2*(len(nodes)-1)
    forks=[i for i in nodes if len(adj[i])==3]
    if not forks:
        assert max(map(len,adj.values()))<=2
        return dict(nodes=nodes,type='A'+str(len(nodes)),roles={str(i):'endpoint'for i in nodes if len(adj[i])<=1})
    assert len(forks)==1 and max(map(len,adj.values()))==3
    center=forks[0];arms=[]
    for j in adj[center]:
        previous=center;node=j;length=1
        while len(adj[node])==2:
            nxt=next(k for k in adj[node]if k!=previous);previous,node=node,nxt;length+=1
        assert len(adj[node])==1;arms.append(length)
    arms.sort()
    if arms[:2]==[1,1]:
        leaves=[i for i in adj[center]if len(adj[i])==1]
        return dict(nodes=nodes,type='D'+str(len(nodes)),roles={str(center):'fork_center',**{str(i):'short_fork_leaf'for i in leaves}})
    assert arms in([1,2,2],[1,2,3],[1,2,4])
    return dict(nodes=nodes,type='E'+str(len(nodes)),roles={})

def target_graph():
    g=[[0]*16 for _ in range(16)]
    for i,j in [(i,i+1)for i in range(11)]+[(0,12),(0,13),(11,14),(11,15)]:g[i][j]=g[j][i]=3
    return g

def distances(g):
    result=[]
    for root in range(len(g)):
        d={root:0};todo=[root]
        for i in todo:
            for j,label in enumerate(g[i]):
                if i!=j and label and j not in d:d[j]=d[i]+1;todo.append(j)
        assert len(d)==len(g);result.append(d)
    return result

def local_interfaces(local,small):
    target=target_graph();ds=distances(target)
    assert target==local['target']
    old=next(r for r in small['records']if r['order']==16 and r['kind']=='D')
    assert old['target']==target
    paid=set(map(tuple,old['retained_paid_pairs']));assert len(paid)==111
    half={p for p in paid if ds[p[0]][p[1]]>9}
    assert len(half)==10 and half==set(map(tuple,local['half_weight_pairs']))
    weights={p:F(1,2)if p in half else F(p in paid)for p in combinations(range(16),2)}
    assert weights=={tuple(r['pair']):F(r['weight'])for r in local['weighted_pairs']}
    assert sum(weights.values())==106 and len(weights)==120
    assert set(weights)-paid==set(map(tuple,local['inherited_zero_pairs']))
    assert not local['added_finite_coupons']
    features={0,11}
    triangles={p for p in weights if ds[p[0]][p[1]]<=5 or(ds[p[0]][p[1]]<=8 and any(ds[i][f]<=5 for i in p for f in features))}
    short={p for p in weights if ds[p[0]][p[1]]<=8}
    assert triangles==set(map(tuple,local['triangle_pairs_retained_at_weight_one']))
    assert short==set(map(tuple,local['short8_pairs_retained_at_weight_one']))
    assert all(weights[p]==1 for p in triangles|short)
    autos=[]
    for reverse,left,right in product((False,True),repeat=3):
        a=list(range(16))
        if left:a[12],a[13]=a[13],a[12]
        if right:a[14],a[15]=a[15],a[14]
        if reverse:a=[11-x if x<12 else {12:14,13:15,14:12,15:13}[x]for x in a]
        assert all(target[i][j]==target[a[i]][a[j]]for i in range(16)for j in range(16))
        assert {tuple(sorted(a[i]for i in p))for p in half}==half
        autos.append(tuple(a))
    assert set(autos)==set(map(tuple,local['automorphisms']))
    affected={r['case']:r for r in old['case_embeddings']if set(map(tuple,r['marked_images']))&half}
    assert sorted(affected)==local['affected_cases']==[r['index']for r in local['cases']]
    assert len(affected)==54
    totals=Counter();map_counts=Counter();checked_grams=set();corners=0;half_records=[]
    for case in local['cases']:
        core=freeze(case['core_gram']);J=inverse(core)
        assert J==freeze(case['core_inverse'])and inertia(core)==(len(core)-1,1,0)
        con=case['constraints']
        assert all(J[i][i]==F(con['projection_diagonals'][str(i)])for i in con['allowed_contacts'])
        expected={tuple(mp)for mp in affected[case['index']]['embeddings']if tuple(sorted(mp[case['base_nodes'].index(i)]for i in case['preferred_corner']))in half}
        seen=set()
        for row in case['completions']:
            kind=row['outcome'];totals[kind]+=1;map_counts[kind]+=len(row['embeddings'])
            for e in row['embeddings']:
                mp=tuple(e['map']);assert mp in expected and mp not in seen;seen.add(mp)
                assert e['marked_image']==sorted(mp[case['base_nodes'].index(i)]for i in case['preferred_corner'])
                assert all((int(2-core[a][b])if a!=b and core[a][b]else 0)==target[mp[i]][mp[j]]
                    for i,a in enumerate(case['base_nodes'])for j,b in enumerate(case['base_nodes']))
            required=row['required_core_contacts']
            if kind=='excluded_multiple_core_contact_roots':
                assert len(required)>1 and case['method']=='lorentz_residual';continue
            if kind=='excluded_projection_above_norm':
                assert len(required)==1 and len(required[0]['core_nodes'])==1
                node=required[0]['core_nodes'][0]
                assert J[node][node]==F(row['projection_square'])>2;continue
            g=freeze(row['gram']);checked_grams.add(g)
            assert sub(g,range(len(core)))==core
            assert inertia(g)==tuple(row['signature'])
            image=row['target_to_literal'];affine=sub(g,image)
            assert affine==tuple(tuple(F(2 if i==j else -bool(target[i][j]))for j in range(16))for i in range(16))
            assert inertia(affine)==(15,0,1)
            if kind=='excluded_lorentz_index':assert inertia(g)[1]>1;continue
            assert kind=='sealed_completion'and not row['new_roots']and not required
            assert con['candidate_nodes']==[]
            white=row['white'];black=set(range(len(g)))-set(white)
            assert len(black)==4 and inertia(sub(g,white))==(len(white),0,0)
            for b in black:assert inertia(sub(g,sorted(set(white)|{b})))==(len(white)+1,0,0)
            z=list(map(F,row['positive_affine_kernel']));assert len(z)==16 and all(x>0 for x in z)
            assert all(sum(affine[i][j]*z[j]for j in range(16))==0 for i in range(16))
            finite=F();ideal=F()
            for c in row['corners']:
                corners+=1;marks=c['marks'];ids=sorted(set(white)|set(marks));sig=inertia(sub(g,ids))
                assert sig==tuple(c['signature'])
                if c['kind']=='affineD15_order16':
                    assert marks==case['preferred_corner']and sig==(len(ids)-1,0,1)
                    assert sorted(image)in components(g,ids)
                    assert F(c['ideal_weight'])==F(1,2);ideal+=F(1,2);continue
                assert c['kind']=='finite'and sig==(len(ids),0,0)
                blocks=[finite_component(g,ns)for ns in components(g,ids)]
                assert blocks==c['components']
                chosen=[next(i for i,b in enumerate(blocks)if m in b['nodes'])for m in marks]
                assert chosen==c['marked_component_indices']
                unit=int(chosen[0]==chosen[1]and len(blocks[chosen[0]]['nodes'])<=17)
                halfroles=[]
                for i,j in((0,1),(1,0)):
                    a,b=blocks[chosen[i]],blocks[chosen[j]];role=b['roles'].get(str(marks[j]))
                    if len(a['nodes'])<=5 and b['type']in{'D'+str(n)for n in range(11,17)}and role in('fork_center','short_fork_leaf'):
                        halfroles.append(dict(small_type=a['type'],small_mark=marks[i],large_type=b['type'],large_mark=marks[j],large_role=role))
                coupon=F(1,2)if halfroles else F()
                assert unit==c['finite_unit']and coupon==F(c['finite_half'])and halfroles==c['coupon_roles']
                finite+=unit+coupon
                if halfroles:half_records.append(dict(case=case['index'],marks=marks,roles=halfroles))
            assert ideal==F(row['ideal_total'])==F(1,2)
            assert finite==F(row['finite_payment_total'])>=F(1,2)
            assert ideal+finite==F(row['total_payment'])>=1
        assert seen==expected
    assert dict(totals)==local['completion_outcomes']
    assert sum(totals.values())==66 and sum(map_counts.values())==local['affected_map_count']==328
    assert corners==40 and len(half_records)==10
    return dict(target_order=16,affine_type='D15',old_cost=111,new_cost='106',half_weight_pairs=sorted(half),
        invariant_automorphisms=len(autos),affected_cases=sorted(affected),affected_maps=sum(map_counts.values()),
        completion_outcomes=dict(totals),completion_map_counts=dict(map_counts),
        distinct_full_grams_independently_checked=len(checked_grams),literal_sealed_corners=corners,
        finite_half_records=half_records,triangle_units=len(triangles),short8_units=len(short),
        scope='Independent Fraction inverse/inertia and saved induced-map/completion/corner/payment checks. Completeness of the original overlap and induced-map censuses, plus geometric contact and no-bridge lemmas, remains inherited and is mathematically reviewed.')

def interfaces(docs,data):
    primary=docs[PRIMARY];base=docs[BASE];review=docs[BASE+'_review'];finite=docs[FINITE]
    assert primary['dimension_bound']==70 and primary['inherited_dimension_bound']==71
    assert base['dimension_bound']==review['dimension_bound']==71
    assert primary['excluded_dimensions']==[71]
    assert primary['unrestricted_finite_volume']and not primary['arithmetic_hypothesis']and primary['cusp_restriction']is None
    assert primary['finite_vertex_budget']==finite['finite_vertex_budget']=='8n'
    assert primary['finite_budget_dimension_maximum']==finite['dimension_maximum']==75
    assert review['composition_interfaces']['finite_budget']['independently_enumerated_endpoints']==13034
    assert len(finite['dimension_endpoints'])==13034
    local=local_interfaces(docs[LOCAL],docs[SMALL])
    costs={k:dict(v)for k,v in base['small_order_costs'].items()}
    assert costs['16']=={'A':80,'B':107,'C':103,'D':111};costs['16']['D']=106
    assert primary['small_order_costs']==costs
    assert all(F(max(row.values()),int(m))==budget(int(m)-1)for m,row in costs.items())
    env=dict(base['ideal_rank_envelope']);env['15']='107/16';assert primary['ideal_rank_envelope']==env
    assert data['concentration']['terminal_repeated_ranks']==primary['concentration']['repeated_cap_ranks']
    for key in('affected_cases','affected_maps','completion_outcomes'):
        assert primary['weighted_order16'][key]==local[key]
    return dict(weighted_order16=local,small_order_costs=costs,ideal_rank_envelope=env,
        inherited_finite_budget=dict(budget='8n',maximum_dimension=75,previous_independent_endpoints=13034),
        retained_tail=71,finite_slack='19/2',strict_global_slack=350,
        scope='The whole-factor local proof and modified face-payment composition receive a complete mathematical/source read. The unchanged finite-budget enumeration and complete previous theorem are frozen dependencies, not repeated calculations.')

def main(prepare=False,compare_only=False,verify=False):
    data=arithmetic(not prepare)
    if prepare:
        Path('/private/tmp/dimension70_independent_scalar.json').write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')
        print('PASS_INDEPENDENT_SCALAR_PREPARATION',len(data['profiles']),data['dimension_minima'],'boundary',data['boundary']['margin'],flush=True);return
    stems=[PRIMARY,BASE,BASE+'_review',LOCAL,SMALL,FINITE]
    docs={t:json.loads((ROOT/f'results/{t}.json').read_text())for t in stems}
    data['composition_interfaces']=interfaces(docs,data)
    if compare_only:
        print('PASS_INDEPENDENT_PRIMARY_COMPARISON',len(data['profiles']),'profiles;',data['composition_interfaces']['weighted_order16'],flush=True);return
    paths=[Path(__file__)]
    for t in stems:
        assert docs[t]['status'].startswith('PASS'),t
        for p,h in docs[t]['sources'].items():assert digest(ROOT/p)==h,(t,p)
        paths.extend(ROOT/f'{folder}/{t}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths.extend(ROOT/f'results/{t}.md'for t in('fresh_cusp_link_shape_curvature_20261017','fresh_cusp_factor_payment_framework_20261017'))
    data.update(status='PASS_INDEPENDENT_DIMENSION70_SCALARS_AND_COMPOSITION',dimension_bound=70,
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        scope='No producer arithmetic imports:1432 scalar profiles and15 all-rank cost branches,54 saved support cases/328maps,66 completion outcomes and40 sealed finite/ideal corners. Original exhaustive core/map censuses and unchanged13034 finite-budget endpoints remain frozen inputs; full mathematical/source composition reviewed.')
    data=json.loads(json.dumps(data));path=ROOT/f'results/{TAG}.json'
    if verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')
    print(data['status'],len(data['profiles']),'profiles; full saved local interfaces',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--compare-only',action='store_true');p.add_argument('--verify',action='store_true')
    a=p.parse_args();assert sum((a.prepare,a.compare_only,a.verify))<=1;main(a.prepare,a.compare_only,a.verify)
