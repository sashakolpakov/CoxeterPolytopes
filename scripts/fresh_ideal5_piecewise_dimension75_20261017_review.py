#!/usr/bin/env python3
"""Independent stdlib Fraction scalar and composition review of dimension75.

No producer/project arithmetic imports. Sliding-window convolution and a
bounded-composition cumulative-tail identity independently compute curvature.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb,prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_ideal5_piecewise_dimension75_20261017'
TAG=PRIMARY+'_review'
LOCAL='fresh_finite8_ideal5_combined_payments_20261017'
CENSUS='fresh_ideal5_small_factor_payments_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
PAIR='fresh_d17_singleton_pair_coupon_20261017'
FINITE='fresh_finite8_a4_d13_neighbor_coupon_budget_20261017'
BASE='fresh_ideal5_piecewise_dimension77_20261017'

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
    return F(r,2)if r<=12 else F(87,14)if r==13 else F(20,3)if r==14 else F(111,16)if r==15 else F(126,17)if r==16 else F(163,23)

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

def anchors(total, allowed=(12,14,16,17)):
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
    for d in range(75,77):
        shapes=set()
        for high in(0,*range(18,d+1)):
            for low in range(12):
                for thirteen in(0,1):
                    for fifteen in(0,1):
                        left=d-high-low-13*thirteen-15*fifteen
                        if left<0:continue
                        exceptional=([high]if high else[])+([low]if low else[])+[13]*thirteen+[15]*fifteen
                        for rest in anchors(left):shapes.add(tuple(sorted(exceptional+rest)))
        group=[scalar(rs)for rs in sorted(shapes)]
        assert all(F(r['margin'])>0 for r in group)
        profiles+=group;n=d+1
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>8*n
        m=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=m['margin'],
            minimizing_ranks=m['ranks'],finite_slack=str(beta-8*n),strict_global_slack=5*(n-1)))
    assert len(profiles)==2164
    assert [m['profiles']for m in minima]==[1053,1111]
    assert budget(12)+budget(14)-2*budget(13)==F(5,21)>0
    assert budget(14)+budget(16)-2*budget(15)==F(83,408)>0
    assert min(F(r['margin'])for r in profiles)==F(11615,448188)
    boundary=scalar((12,12,12,12,12,14))
    assert F(boundary['margin'])==F(-16873,171366)
    if primary is not None:
        by={tuple(r['ranks']):r for r in primary['reduced_profiles']}
        assert set(by)=={tuple(r['ranks'])for r in profiles}and len(by)==len(primary['reduced_profiles'])==2164
        for r in profiles:
            assert {k:v for k,v in r.items()if k!='independent_cumulative_G'}==by[tuple(r['ranks'])]
        assert minima==primary['dimension_minima']
        assert {k:v for k,v in boundary.items()if k!='independent_cumulative_G'}==primary['boundary_at_n75']
    return dict(profiles=profiles,dimension_minima=minima,boundary=boundary,
        minimum_margin=str(min(F(r['margin'])for r in profiles)),
        independent_affine_feature_maxima=affine_feature_maxima(),
        formulas=['G=sum_j h_j min(j,d-j)',
            'G=floor(d/2)*V-2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1)',
            'kappa/V=3d/4-(5/2)*E(abs(sum U_i-d/2))'])

def inertia(matrix):
    """Rational congruence elimination, independent of all producer helpers."""
    a=[[F(x)for x in row]for row in matrix];answer=[0,0,0]
    while a:
        n=len(a);pivot=next((i for i in range(n)if a[i][i]),None)
        if pivot is None:
            edge=next(((i,j)for i in range(n)for j in range(i+1,n)if a[i][j]),None)
            if edge is None:answer[2]+=n;break
            # Replace basis vector e_i by e_i+e_j, giving a nonzero diagonal.
            i,j=edge
            a[i]=[a[i][k]+a[j][k]for k in range(n)]
            for k in range(n):a[k][i]+=a[k][j]
            pivot=i
        order=[pivot]+[i for i in range(n)if i!=pivot]
        a=[[a[i][j]for j in order]for i in order];p=a[0][0]
        answer[0 if p>0 else 1]+=1
        a=[[a[i][j]-a[i][0]*a[0][j]/p for j in range(1,n)]for i in range(1,n)]
    return answer


def literal_core_check(saved):
    edges=[(0,1),(0,2),(2,3),(0,4),(4,5),(5,6),(6,7),(7,8),(8,9)]
    edges+= [(i+10,j+10 if j<9 else 9)for i,j in edges.copy()]
    g=[[F(2 if i==j else 0)for j in range(19)]for i in range(19)]
    for i,j in edges:g[i][j]=g[j][i]=F(-1)
    assert g==[[F(x)for x in row]for row in saved['gram']]
    principal=lambda nodes:[[g[i][j]for j in nodes]for i in nodes]
    assert inertia(g)==saved['signature']==[17,1,1]
    radical=saved['radical'];assert len(radical)==19
    assert all(sum(g[i][j]*radical[j]for j in range(19))==0 for i in range(19))
    assert radical[3]==-2 and radical[13]==2
    white=sorted(set(range(19))-{3,4,13,14});assert white==saved['white']
    assert inertia(principal(white))==[15,0,0]
    for b in(3,4,13,14):assert inertia(principal(white+[b]))==[16,0,0]
    affine=sorted(set(range(19))-{3,13});assert affine==saved['affine_nodes']
    kernel=saved['affine_positive_kernel'];a=principal(affine)
    assert len(kernel)==17 and min(kernel)>0
    assert all(sum(a[i][j]*kernel[j]for j in range(17))==0 for i in range(17))
    assert inertia(a)==[16,0,1]
    paid=[]
    for corner in saved['corners']:
        nodes=sorted(white+corner['marks']);assert inertia(principal(nodes))==corner['signature']
        if not F(corner['half_coupon']):continue
        assert F(corner['half_coupon'])==F(1,2)and corner['signature']==[17,0,0]
        assert {c['type']for c in corner['components']}=={'A4','D13'}
        for c in corner['components']:
            ns=c['nodes'];mark=c['mark'];assert mark in corner['marks']and mark in ns
            assert len(ns)==int(c['type'][1:])
            adj={i:{j for j in ns if g[i][j]<0}for i in ns}
            assert all(g[i][j]==0 for i in ns for j in nodes if j not in ns)
            seen={ns[0]};todo=list(seen)
            for i in todo:
                for j in adj[i]-seen:seen.add(j);todo.append(j)
            assert seen==set(ns)
            if c['type']=='A4':
                assert sorted(map(len,adj.values()))==[1,1,2,2]and len(adj[mark])==1
                assert c['role']=='endpoint'
            else:
                center=next(i for i in ns if len(adj[i])==3)
                assert sorted(map(len,adj.values()))==[1]*3+[2]*9+[3]
                assert len([i for i in adj[center]if len(adj[i])==1])==2
                assert mark in adj[center]and len(adj[mark])==2
                assert c['role']=='long_arm_neighbor_of_fork'
        paid.append(corner['marks'])
    assert paid==[[3,14],[4,13]]
    return dict(core_signature=[17,1,1],white_signature=[15,0,0],strict_edge_checks=4,
        affine_signature=[16,0,1],radical=radical,positive_kernel=kernel,
        independently_checked_coupon_corners=paid,total_new_payment='1')


def interfaces(inputs,arithmetic_data):
    primary,local,pair,finite,small=(inputs[s]for s in(PRIMARY,LOCAL,PAIR,FINITE,SMALL))
    assert primary['dimension_bound']==75 and primary['inherited_dimension_bound']==77
    assert primary['excluded_dimensions']==[76,77]
    assert primary['unrestricted_finite_volume']is True and primary['arithmetic_hypothesis']is False
    assert primary['cusp_restriction']is None
    assert inputs[BASE]['dimension_bound']==inputs[BASE+'_review']['dimension_bound']==77
    assert local['finite_vertex_budget']==finite['finite_vertex_budget']==primary['finite_vertex_budget']=='8n'
    assert local['finite_budget_dimension_maximum']==finite['dimension_maximum']==80
    assert local['triangle_minimum']=='2'and local['quadrilateral_minimum']=='1'
    assert local['high_overlap_contexts']==2614 and local['selected_small_patterns']==len(small['patterns'])==1153
    expected={'1_to_12':'r/2','13':'87/14','14':'20/3','15':'111/16','16':'126/17','17_and_larger':'163/23'}
    assert primary['ideal_rank_envelope']==expected
    prior=dict(local['ideal_rank_envelope']);assert prior.pop('16')=='127/17'
    assert prior=={k:v for k,v in expected.items()if k!='16'}
    targetrow=next(r for r in small['records']if r['order']==17 and r['kind']=='D')
    old=set(map(tuple,targetrow['retained_paid_pairs']));new=set(map(tuple,pair['new_paid_pairs']))
    assert len(old)==127 and len(new)==126 and old-new=={(1,11)}and new<old
    assert pair['old_cost']==127 and pair['new_cost']==126
    assert pair['removed_pair']==[1,11]and pair['affected_overlap_cases']==[2518]
    assert pair['affine_type']=='D16'and pair['diagram_order']==17
    affected=sorted(r['case']for r in targetrow['case_embeddings']if [1,11]in r['marked_images'])
    assert affected==[2518]
    g=targetrow['target'];adj=[{j for j,x in enumerate(row)if x and i!=j}for i,row in enumerate(g)]
    distance={1:0};todo=[1]
    for i in todo:
        for j in adj[i]:
            if j not in distance:distance[j]=distance[i]+1;todo.append(j)
    assert len(distance)==17 and distance[11]==10
    autos=pair['automorphisms'];assert len(autos)==len({tuple(a)for a in autos})==8
    for perm in autos:
        assert sorted(perm)==list(range(17))and {perm[1],perm[11]}=={1,11}
        assert all(g[i][j]==g[perm[i]][perm[j]]for i in range(17)for j in range(17))
        assert {tuple(sorted((perm[i],perm[j])))for i,j in new}==new
    sets={(r['order'],r['kind']):set(map(tuple,r['retained_paid_pairs']))for r in small['records']}
    sets[17,'D']=new
    assert len(local['triangle_pair_retention'])==16
    for row in local['triangle_pair_retention']:
        assert set(map(tuple,row['required_pairs']))<=sets[row['order'],row['kind']]
    assert set(map(tuple,pair['triangle_pairs_retained']))<=new
    costs={str(m):{kind:len(sets[m,kind])for kind in'ABCD'}for m in range(14,18)}
    assert costs==primary['small_order_costs']
    for m in range(14,18):assert F(max(costs[str(m)].values()),m)==budget(m-1)
    assert next(r for r in local['short8_reduction']if r['order']==17 and r['kind']=='D')['unpaid_short8_pairs']==[]
    core=literal_core_check(pair['literal_core'])
    component=finite['component_count']
    assert component['A4_endpoints']==2 and component['D13_new_marks']==1
    assert F(component['half_weight'])==F(1,2)
    assert F(component['cost_per_A4_D13_component_pair'])==1
    assert F(component['D13_unit_cost'])==F(13*12,2)==78
    assert F(component['D13_base_slack'])==8*13-78==26
    assert F(component['additional_D13_slack'])==26-8==18
    assert component['additional_D13_rank']==2
    assert 80-11-2==67 and F(67,4)-18==F(-5,4)
    branches=[branch for family in arithmetic_data['independent_affine_feature_maxima']for branch in family['linear_branches']]
    keys=('kind','lower_order','upper_order','slope','intercept')
    saved=[{k:r[k]for k in keys}for r in inputs[CENSUS]['basic_affine_cost']['all_rank_branches']]
    order=lambda r:(r['kind'],r['lower_order'])
    assert sorted(branches,key=order)==sorted(saved,key=order)
    assert primary['effective_compact_facet_ambient_bound']==75
    assert primary['unchanged_compact_arithmetic_bound']==21 and primary['unchanged_arbitrary_compact_bound']==29
    return dict(independently_checked_literal_core=core,removed_pair=[1,11],removed_pair_distance=10,
        old_cost=127,new_cost=126,affected_saved_map_cases=affected,checked_saved_automorphisms=8,
        elementary_triangle_inclusions=16,small_factor_costs=costs,
        independent_all_rank_branches=branches,finite_extra_correction_at_max_a='-5/4 per D13',
        inherited_high_overlap_cases=2614,inherited_preferred_patterns=1153,
        scope='Independent2164 scalar reconstruction and small literal/interface checks plus complete written/source composition review. No duplicate2614 overlap or1153 embedding census; saved automorphisms are validated, not independently enumerated.')


def run():
    stems=(PRIMARY,LOCAL,PAIR,FINITE,SMALL,CENSUS,BASE,BASE+'_review')
    inputs={stem:json.loads((ROOT/f'results/{stem}.json').read_text())for stem in stems}
    for stem,doc in inputs.items():
        assert doc['status'].startswith('PASS')
        for name,h in doc['sources'].items():assert digest(ROOT/name)==h,(stem,name)
    out=arithmetic(True);out['composition_review']=interfaces(inputs,out)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths +=[ROOT/f'results/{stem}.md'for stem in('fresh_cusp_link_shape_curvature_20261017',
        'fresh_cusp_factor_payment_framework_20261017','fresh_distance14_unrestricted_dimension156_20261017')]
    out.update(status='PASS_INDEPENDENT_IDEAL5_DIMENSION75_REVIEW',dimension_bound=75,
        inherited_dimension_bound=77,excluded_dimensions=[76,77],
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        external_refereeing=False,proof_assistant_formalization=False)
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--arithmetic-only',action='store_true');p.add_argument('--compare',action='store_true')
    p.add_argument('--verify',action='store_true');a=p.parse_args()
    assert sum((a.arithmetic_only,a.compare,a.verify))<=1
    if a.arithmetic_only or a.compare:
        out=arithmetic(a.compare)
        print('PASS_INDEPENDENT_CONDITIONAL_75_SCALARS',len(out['profiles']),out['minimum_margin'],out['dimension_minima'],out['boundary']['margin'],flush=True)
    else:
        out=run();path=ROOT/f'results/{TAG}.json'
        if a.verify:assert out==json.loads(path.read_text())
        else:path.write_text(json.dumps(out,indent=2)+'\n')
        print(out['status'],len(out['profiles']),out['minimum_margin'],out['boundary']['margin'],flush=True)
