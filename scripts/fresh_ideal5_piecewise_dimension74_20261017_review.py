#!/usr/bin/env python3
"""Independent stdlib Fraction scalar and composition review of dimension74.

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
PRIMARY='fresh_ideal5_piecewise_dimension74_20261017'
TAG=PRIMARY+'_review'
LOCAL='fresh_finite8_ideal5_combined_payments_20261017'
CENSUS='fresh_ideal5_small_factor_payments_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
PAIR='fresh_order15_affine_b_pair_reduction_20261017'
FINITE='fresh_finite8_a4_d13_neighbor_coupon_budget_20261017'
BASE='fresh_ideal5_piecewise_dimension75_20261017'

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
    return F(r,2)if r<=12 else F(87,14)if r==13 else F(97,15)if r==14 else F(111,16)if r==15 else F(126,17)if r==16 else F(163,23)

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
    for d in range(74,75):
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
    assert len(profiles)==999
    assert [m['profiles']for m in minima]==[999]
    assert budget(12)+budget(14)-2*budget(13)==F(4,105)>0
    assert budget(14)+budget(16)-2*budget(15)==F(7,2040)>0
    assert min(F(r['margin'])for r in profiles)==F(63,74698)
    boundary=scalar((1,12,12,12,12,12,12))
    assert F(boundary['margin'])==F(-168143,742586)
    if primary is not None:
        by={tuple(r['ranks']):r for r in primary['reduced_profiles']}
        assert set(by)=={tuple(r['ranks'])for r in profiles}and len(by)==len(primary['reduced_profiles'])==999
        for r in profiles:
            assert {k:v for k,v in r.items()if k!='independent_cumulative_G'}==by[tuple(r['ranks'])]
        assert minima==primary['dimension_minima']
        assert {k:v for k,v in boundary.items()if k!='independent_cumulative_G'}==primary['boundary_at_n74']
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


def inverse(matrix):
    n=len(matrix);a=[[F(x)for x in row]+[F(i==j)for j in range(n)]for i,row in enumerate(matrix)]
    for j in range(n):
        k=next(i for i in range(j,n)if a[i][j]);a[j],a[k]=a[k],a[j]
        scale=a[j][j];a[j]=[x/scale for x in a[j]]
        for i in range(n):
            if i!=j:
                scale=a[i][j];a[i]=[x-scale*y for x,y in zip(a[i],a[j])]
    return [row[n:]for row in a]


def component_type(g,nodes):
    assert inertia([[g[i][j]for j in nodes]for i in nodes])==[len(nodes),0,0]
    adj={i:{j for j in nodes if i!=j and g[i][j]}for i in nodes}
    assert all(g[i][j]==-1 for i in nodes for j in adj[i])
    seen={nodes[0]};todo=list(seen)
    for i in todo:
        for j in adj[i]-seen:seen.add(j);todo.append(j)
    assert seen==set(nodes)and sum(map(len,adj.values()))==2*(len(nodes)-1)
    degrees=sorted(map(len,adj.values()))
    if degrees==[1,1]+[2]*(len(nodes)-2):return 'A'+str(len(nodes))
    assert degrees==[1]*3+[2]*(len(nodes)-4)+[3]
    fork=next(i for i in nodes if len(adj[i])==3)
    assert sum(len(adj[j])==1 for j in adj[fork])==2
    return 'D'+str(len(nodes))


def seven_local_checks(local):
    edges=[(0,1),(0,2),(2,3),(0,4),(4,5),(5,6),(6,7),(7,8),(8,9)]
    rows=[]
    for row in local['cases']:
        assert row['left'][0]==row['right'][0]==25
        right={b:a for a,b in row['match']};n=10
        for i in range(10):
            if i not in right:right[i]=n;n+=1
        g=[[F(2 if i==j else 0)for j in range(n)]for i in range(n)]
        for i,j in edges+[(right[i],right[j])for i,j in edges]:g[i][j]=g[j][i]=F(-1)
        assert g==[[F(x)for x in r]for r in row['gram']]
        assert inertia(g)==row['signature']==[n-1,1,0]
        j=inverse(g)
        black=set(row['left'][1])|{right[i]for i in row['right'][1]}
        white=sorted(set(range(n))-black);marks=row['marks'];corner=sorted(white+marks)
        assert len(black)==4 and len(set(marks))==2 and set(marks)<=black
        assert inertia([[g[i][k]for k in white]for i in white])==[n-4,0,0]
        assert inertia([[g[i][k]for k in corner]for i in corner])==[n-2,0,0]
        nodes=row['marked_component'];assert set(marks)<=set(nodes)<=set(corner)
        assert all(g[i][k]==0 for i in nodes for k in corner if k not in nodes)
        typ=component_type(g,nodes)
        contacts={r['contact_node']:r for r in row['allowed_single_contacts']}
        for i,r in contacts.items():
            assert j[i][i]==F(r['projection_square'])and 2-j[i][i]==F(r['residual_square'])
            assert 2*j[i][i]>2
        if row['index']!=1641:
            assert n==17 and typ=='D14'and set(contacts)=={3,13}
            assert not set(contacts)&set(nodes)
            assert all(j[i][i]==2 for i in contacts)
            rows.append(dict(index=row['index'],signature=[16,1,0],finite_type='D14',allowed_contacts=[3,13],
                direct_contacts_to_marked_component=[],whole_component_order=14))
        else:
            assert n==16 and typ=='A12'and set(contacts)=={1,3,11,13}
            assert set(nodes)&set(contacts)=={1,11}
            products=[]
            for i in(1,11):
                assert j[i][i]==2 and sum(bool(g[i][k])for k in nodes if k!=i)==1
                for k in sorted(contacts):
                    assert j[i][k]>0 or j[i][k]<=-2
                    products.append([i,k,str(j[i][k])])
                extended=[r+[F(-int(k==i))]for k,r in enumerate(g)]
                extended.append([F(-int(k==i))for k in range(n)]+[F(2)])
                assert component_type(extended,nodes+[n])=='A13'
            rows.append(dict(index=1641,signature=[15,1,0],finite_type='A12',allowed_contacts=[1,3,11,13],
                endpoint_extension_types=['A13','A13'],mixed_products=products,maximum_whole_component_order=13))
        assert row['finite_payment']==1
    assert [r['index']for r in rows]==[1284,1623,1641,1646,1666,1933,2326]
    return rows


def interfaces(docs,arithmetic_data):
    primary,base,local,maps,payments=(docs[s]for s in(PRIMARY,BASE,PAIR,SMALL,LOCAL))
    assert primary['dimension_bound']==74 and primary['inherited_dimension_bound']==75
    assert base['dimension_bound']==docs[BASE+'_review']['dimension_bound']==75
    assert primary['excluded_dimensions']==[75]and primary['unrestricted_finite_volume']is True
    assert primary['arithmetic_hypothesis']is False and primary['cusp_restriction']is None
    assert primary['finite_vertex_budget']==base['finite_vertex_budget']=='8n'
    assert primary['finite_budget_dimension_maximum']==base['finite_budget_dimension_maximum']==80
    assert local['affine_type']=='B14'and local['diagram_order']==15
    assert local['old_cost']==100 and local['new_cost']==97 and local['new_finite_coupons']==[]
    rec=next(r for r in maps['records']if r['order']==15 and r['kind']=='B')
    old=set(map(tuple,rec['retained_paid_pairs']));new=set(map(tuple,local['new_paid_pairs']))
    drop={(0,11),(11,13),(11,14)}
    assert len(old)==100 and len(new)==97 and new<old and old-new==drop
    assert set(map(tuple,local['removed_pairs']))==drop
    affected=sorted(r['case']for r in rec['case_embeddings']if set(map(tuple,r['marked_images']))&drop)
    assert affected==local['affected_cases']==[1284,1623,1641,1646,1666,1933,2326]
    g=rec['target'];distances=[]
    for a,b in sorted(drop):
        distance={a:0};todo=[a]
        for i in todo:
            for k,v in enumerate(g[i]):
                if i!=k and v and k not in distance:distance[k]=distance[i]+1;todo.append(k)
        assert len(distance)==15 and distance[b]>9;distances.append([a,b,distance[b]])
    assert distances==[[0,11,11],[11,13,12],[11,14,12]]
    perm=list(range(15));perm[13],perm[14]=14,13
    assert all(g[i][j]==g[perm[i]][perm[j]]for i in range(15)for j in range(15))
    assert {tuple(sorted((perm[i],perm[j])))for i,j in drop}==drop
    retained={ (r['order'],r['kind']):set(map(tuple,r['retained_paid_pairs']))for r in maps['records']}
    retained[15,'B']=new
    # Keep the already proved D17 singleton deletion from the75 theorem.
    priorcosts=base['small_order_costs'];assert priorcosts['17']['D']==126
    oldD=retained[17,'D'];assert len(oldD)==127 and(1,11)in oldD
    retained[17,'D']=oldD-{(1,11)}
    for row in payments['triangle_pair_retention']:
        assert set(map(tuple,row['required_pairs']))<=retained[row['order'],row['kind']]
    assert len(payments['triangle_pair_retention'])==16
    assert next(r for r in payments['short8_reduction']if r['order']==15 and r['kind']=='B')['unpaid_short8_pairs']==[]
    costs={str(m):{k:len(retained[m,k])for k in'ABCD'}for m in range(14,18)}
    assert costs==primary['small_order_costs']
    for m in range(14,18):assert F(max(costs[str(m)].values()),m)==budget(m-1)
    expected=dict(base['ideal_rank_envelope']);expected['14']='97/15'
    assert primary['ideal_rank_envelope']==expected
    branches=[b for r in arithmetic_data['independent_affine_feature_maxima']for b in r['linear_branches']]
    keys=('kind','lower_order','upper_order','slope','intercept')
    saved=[{k:r[k]for k in keys}for r in docs[CENSUS]['basic_affine_cost']['all_rank_branches']]
    sort=lambda r:(r['kind'],r['lower_order'])
    assert sorted(branches,key=sort)==sorted(saved,key=sort)
    assert primary['unchanged_compact_arithmetic_bound']==21 and primary['unchanged_arbitrary_compact_bound']==29
    assert primary['effective_compact_facet_ambient_bound']==74
    return dict(removed_pair_distances=distances,affected_saved_map_cases=affected,
        local_exact_checks=seven_local_checks(local),triangle_inclusions=16,small_order_costs=costs,
        additional_finite_cost='0',inherited_dimension_bound=75,
        independent_all_rank_branches=branches,
        scope='Independent999-row scalar reconstruction and seven literal rational-core/interface checks, with full source/mathematical composition read. The complete all-label supporting-white domains and overlap/map census remain frozen inputs; no duplicate census.')


def run():
    stems=(PRIMARY,BASE,BASE+'_review',PAIR,SMALL,LOCAL,CENSUS)
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in stems}
    for s,d in docs.items():
        assert d['status'].startswith('PASS')
        for path,h in d['sources'].items():assert digest(ROOT/path)==h,(s,path)
    out=arithmetic(True);out['composition_review']=interfaces(docs,out)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths +=[ROOT/f'results/{s}.md'for s in('fresh_cusp_link_shape_curvature_20261017',
        'fresh_cusp_factor_payment_framework_20261017','fresh_distance14_unrestricted_dimension156_20261017')]
    out.update(status='PASS_INDEPENDENT_IDEAL5_DIMENSION74_REVIEW',dimension_bound=74,
        inherited_dimension_bound=75,excluded_dimensions=[75],
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        external_refereeing=False,proof_assistant_formalization=False)
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--arithmetic-only',action='store_true')
    p.add_argument('--compare',action='store_true');p.add_argument('--verify',action='store_true')
    a=p.parse_args();assert sum((a.arithmetic_only,a.compare,a.verify))<=1
    if a.arithmetic_only or a.compare:
        d=arithmetic(a.compare)
        print('PASS_INDEPENDENT_CONDITIONAL_74_SCALARS',len(d['profiles']),d['minimum_margin'],d['boundary']['margin'],flush=True)
    else:
        d=run();path=ROOT/f'results/{TAG}.json'
        if a.verify:assert d==json.loads(path.read_text())
        else:path.write_text(json.dumps(d,indent=2)+'\n')
        print(d['status'],len(d['profiles']),d['minimum_margin'],d['boundary']['margin'],flush=True)


