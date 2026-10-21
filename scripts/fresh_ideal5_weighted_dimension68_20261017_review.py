#!/usr/bin/env python3
"""Independent dimension68 scalar and local-interface review.

No project arithmetic imports. Primary bindings are added only after the
primary source/data/report and all local interfaces have frozen.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from math import comb,prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_ideal5_weighted_dimension68_20261017'
TAG=PRIMARY+'_review'
BASE='fresh_ideal5_weighted_dimension69_20261017'
SUPPORT='fresh_order12_payment_support_20261017'
A='fresh_cycle12_ideal5_transfer_20261017'
B='fresh_order12_affine_b_half_transfer_20261017'
C='fresh_order12_affine_c_endpoint_reduction_20261017'
D='fresh_affine_d11_leaf_pair_transfer_20261017'
A17='fresh_finite_a17_endpoint_payment_support_20261017'
FOREST='fresh_a17_endpoint_d11_upgrade_forest_budget_20261017'
FINITE='fresh_a11_cycle_coupon_domination_20261017'
TRI='fresh_ideal5_finite8_triangle_payments_20261017'
COMBINED='fresh_finite8_ideal5_combined_payments_20261017'
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
    return F(r,2)if r<=10 else F(65,12)if r==11 else F(75,13)if r==12 else F(43,7)if r==13 else F(32,5)if r==14 else F(107,16)if r==15 else F(120,17)if r==16 else F(163,23)

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

def anchors(total, allowed=(10,11,13,14,15,16,17)):
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
    for d in (68,):
        shapes=set()
        for high in(0,*range(18,d+1)):
            for low in range(10):
                for twelve in(0,1):
                    left=d-high-low-12*twelve
                    if left<0:continue
                    exceptional=([high]if high else[])+([low]if low else[])+[12]*twelve
                    for rest in anchors(left):shapes.add(tuple(sorted(exceptional+rest)))
        group=[scalar(rs)for rs in sorted(shapes)]
        assert all(F(r['margin'])>0 for r in group)
        profiles+=group;n=d+1
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta==544
        m=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=m['margin'],
            minimizing_ranks=m['ranks'],finite_slack=str(beta-544),strict_global_slack=5*(n-1)))
    assert len(profiles)==2060 and minima==[dict(ambient_dimension=69,profiles=2060,minimum_margin='3641/13608',minimizing_ranks=[11,11,11,11,11,13],finite_slack='0',strict_global_slack=340)]
    assert budget(11)+budget(13)-2*budget(12)==F(23,1092)>0
    boundary=scalar((11,11,11,11,11,12))
    assert boundary['margin']=='411499/3234816'
    if primary is not None:
        by={tuple(r['ranks']):r for r in primary['reduced_profiles']}
        assert set(by)=={tuple(r['ranks'])for r in profiles}and len(by)==len(primary['reduced_profiles'])==len(profiles)
        for r in profiles:
            assert {k:v for k,v in r.items()if k!='independent_cumulative_G'}==by[tuple(r['ranks'])]
        for actual,saved in zip(minima,primary['dimension_minima']):
            assert all(saved[k]==v for k,v in actual.items())
            assert saved['finite_maximum']==saved['finite_curvature_coefficient']=='544'
        assert len(primary['dimension_minima'])==len(minima)==1
        assert {k:v for k,v in boundary.items()if k!='independent_cumulative_G'}==primary['boundary_at_n68']['profile']
        assert primary['boundary_at_n68']['finite_slack']=='-17'
        assert primary['boundary_at_n68']['finite_maximum']=='544'
        assert primary['boundary_at_n68']['finite_curvature_coefficient']=='527'
    return dict(profiles=profiles,dimension_minima=minima,boundary=boundary,
        minimum_margin=str(min(F(r['margin'])for r in profiles)),
        independent_affine_feature_maxima=affine_feature_maxima(),
        concentration=dict(duplicate_rank12_cost_increase=str(budget(11)+budget(13)-2*budget(12)),
            terminal_repeated_ranks=[10,11,13,14,15,16,17],optional_rank12=True,
            at_most_one_low_rank=list(range(1,10)),at_most_one_high_rank_minimum=18,
            termination='Each nontrivial concentration increases the integer sum of squared parts, bounded by d^2.'),
        formulas=['G=sum_j h_j min(j,d-j)',
            'G=floor(d/2)*V-2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1)',
            'kappa/V=3d/4-(5/2)*E(abs(sum U_i-d/2))'])



def canonical(kind):
    n=12;g=[[0]*n for _ in range(n)]
    if kind=='A':edges=[(i,(i+1)%n,3)for i in range(n)]
    elif kind=='C':edges=[(i,i+1,4 if i in(0,10)else 3)for i in range(11)]
    elif kind=='B':edges=[(i,i+1,4 if i==8 else 3)for i in range(9)]+[(0,10,3),(0,11,3)]
    else:
        assert kind=='D';edges=[(i,i+1,3)for i in range(7)]+[(0,8,3),(0,9,3),(7,10,3),(7,11,3)]
    for i,j,w in edges:g[i][j]=g[j][i]=w
    return g


def distances(g):
    allrows=[]
    for start in range(len(g)):
        ds={start:0};queue=[start]
        for i in queue:
            for j,w in enumerate(g[i]):
                if w and j not in ds:ds[j]=ds[i]+1;queue.append(j)
        assert len(ds)==len(g)
        allrows.append([ds[i]for i in range(len(g))])
    return allrows


def internal_cost(typ):
    r=int(typ[1:])
    if r<=17:return F(r*(r-1),2)-int(typ=='A17')
    assert typ[0]in'ABD'
    return F(8*r-(32 if typ[0]=='D'else 34))


def cross_cost(a,b,cycle_coupons=True):
    ra,rb=int(a[1:]),int(b[1:])
    if ra>rb:return cross_cost(b,a,cycle_coupons)
    if ra>5 or rb<=5:return F()
    value=F(3*ra,2)if b[0]=='D'and 11<=rb<=16 else F(ra,2)if b in('A13','A15')else F(ra)if b=='A17'else F()
    if a=='A1'and b in('D11','D12'):value+=1
    if a=='A2'and b=='B10':value+=1
    if a=='A4'and b=='D13':value+=1
    if cycle_coupons and b=='A11'and a in('A1','A2'):value+=1
    return value


def forest_cost(components):
    return sum((internal_cost(t)for t in components),F())+sum((cross_cost(a,b)for a,b in combinations(components,2)),F())


def finite_interface(docs):
    finite=docs[FINITE];old=docs[FOREST];review=docs[FOREST+'_review']
    types=['A'+str(r)for r in range(1,76)]+['B'+str(r)for r in range(2,76)]+['D'+str(r)for r in range(4,76)]+['E6','E7','E8','F4','H3','H4','I2']
    assert len(types)==228 and internal_cost('A11')==internal_cost('D11')==55
    rows=[]
    for t in types:
        before=cross_cost('A11',t,False);after=cross_cost('A11',t);dominant=cross_cost('D11',t)
        assert before==0 and dominant>=after
        rows.append(dict(partner=t,old_A11_cross_twice=int(2*before),new_A11_cross_twice=int(2*after),
            D11_cross_twice=int(2*dominant),domination_slack_twice=int(2*(dominant-after))))
    assert rows==finite['partner_domination']
    maxima=finite['dimension_maxima'];assert maxima==old['variants'][1]['dimension_maxima']
    assert old['variants'][1]['D11_upgrade'] and review['producer_arithmetic_imported']is False
    assert review['variants'][1]['rank69_cost']=='544'
    witnesses=[]
    for n,row in enumerate(maxima):
        types=row['witness_components'];assert sum(int(t[1:])for t in types)==n and 'A11'not in types
        cost=forest_cost(types);assert 2*cost==row['doubled_cost']
        witnesses.append(dict(rank=n,components=types,cost=str(cost)))
    assert len(witnesses)==76 and F(witnesses[69]['cost'])==544
    assert forest_cost(['A1']+['D17']*4)==forest_cost(['A1']+['A17']*4)==544
    assert F((69-1)*(69-5),8)==544 and F(68*(68-6),8)==527
    assert finite['added_coupons']==docs[A]['finite_coupons']
    return dict(partner_domination=rows,attaining_forests=witnesses,
        exact_rank69_cost='544',finite_curvature_at69='544',finite_slack_at69='0',
        finite_cost_at68='544',finite_curvature_at68='527',finite_slack_at68='-17',
        inherited_upper_bound='Frozen complete large-type/Bellman certificate and its no-producer review; not rerun.',
        new_upper_bound_argument='A11 to D11 preserves rank and base and dominates each partner interaction. Large-large interactions vanish, so sequential replacement is valid.')


def local_interfaces(docs):
    allpairs=set(combinations(range(12),2));rows=[]
    for kind in'ABCD':
        g=canonical(kind);ds=distances(g);saved=next(r for r in docs[SUPPORT]['records']if r['kind']==kind)
        assert g==saved['target']
        features={i for i in range(12)if sum(bool(x)for x in g[i])>=3 or any(x>=4 for x in g[i])}
        elementary={p for p in allpairs if ds[p[0]][p[1]]<=5 or(ds[p[0]][p[1]]<=8 and any(ds[i][f]<=5 for i in p for f in features))}
        short8={p for p in allpairs if ds[p[0]][p[1]]<=8}
        if kind=='A':weights={p:F(ds[p[0]][p[1]]<=5)for p in allpairs}
        elif kind=='B':weights={p:F(1,2)if p in((9,10),(9,11))else F(1)for p in allpairs}
        elif kind=='C':weights={p:F(p!=(0,11))for p in allpairs}
        else:weights={p:F(p not in((8,10),(8,11),(9,10),(9,11)))for p in allpairs}
        assert all(weights[p]==1 for p in elementary)
        changed=sorted(p for p,w in weights.items()if w!=1)
        if kind=='A':assert set(changed)==set(map(tuple,docs[A]['pair_rule']['transferred']))
        elif kind=='B':assert weights=={tuple(r['pair']):F(r['weight'])for r in docs[B]['weighted_pairs']}
        elif kind=='C':assert {p for p,w in weights.items()if w}==set(map(tuple,docs[C]['new_paid_pairs']))
        else:assert {p for p,w in weights.items()if w}==set(map(tuple,docs[D]['retained_pairs']))
        if kind!='A':assert all(weights[p]==1 for p in short8)
        cost=sum(weights.values(),F());assert cost=={'A':60,'B':65,'C':65,'D':62}[kind]
        out=dict(kind=kind,order=12,cost=int(cost),changed_pairs=[list(p)for p in changed],
            changed_distances=[ds[i][j]for i,j in changed],elementary_triangle_pairs=[list(p)for p in sorted(elementary)],
            weighted_pairs=[dict(pair=list(p),weight=str(weights[p]))for p in sorted(allpairs)])
        assert out==next(r for r in docs[PRIMARY]['order12_interfaces']if r['kind']==kind)
        rows.append(out)
    b=docs[B];counts=Counter(r['outcome']for r in b['cases'])
    assert counts==Counter(sealed_radical_core=8,whole_finite_B11_before_any_affine_completion=2,zero_or_one_sealed_unit_contact_root=4)
    assert sum(len(r['embeddings'])for r in b['cases'])==b['affected_maps']==28
    sealed=[]
    for r in b['cases']:
        if r['outcome']=='whole_finite_B11_before_any_affine_completion':
            assert r['finite_unit']==1 and r['whole_finite_component']['type']=='B11'
            continue
        groups=[r['sealed_payment']]if 'sealed_payment'in r else[r['no_extra_root_payment'],r['one_extra_root_payment']]
        for z in groups:
            finite=ideal=F()
            for c in z['corners']:
                f=F(c['finite_unit'])+F(c['finite_half']);w=F(c['ideal_weight'])
                if f:assert c['kind']=='finite'and c['signature'][1:]==[0,0]
                if w:assert w==F(1,2)and c['kind']=='affineB11_order12'
                finite+=f;ideal+=w
            assert finite>=F(1,2)and finite+ideal>=1
            assert str(finite)==z['finite_total']and str(ideal)==z['ideal_total']
            sealed.append(dict(case=r['index'],finite=str(finite),ideal=str(ideal)))
    gaps=b['immediate_D13_strict_gap']['rows'];assert len(gaps)==3
    for r in gaps:
        R=r['R'];far=r['far_maximum_cost']
        assert R in(9,11)and F(1)-F(2+far,R)==F(r['strict_gap'])>0
    assert docs[C]['complete_high_patterns']==758 and len(docs[C]['high_case_embeddings'])==15
    assert sum(len(r['embeddings'])for r in docs[C]['high_case_embeddings'])==30
    assert docs[C]['low_interface']['affine_order12_families']=={'B':2,'D':10}
    assert docs[COMBINED]['high_overlap_contexts']==2614 and docs[COMBINED]['low_overlap_contexts']==66
    assert docs[TRI]['rule']['ideal_whole_factor_order_at_most']==9
    assert not b['added_finite_coupons']and not docs[C]['added_finite_coupons']
    assert docs[A17]['candidate_removed_pair']==[0,16]
    assert F(docs[A17]['legacy_replacements']['remaining_total_at_affected_case'])==2
    return dict(order12_interfaces=rows,B_case_partition=dict(counts),B_sealed_payment_totals=sealed,
        immediate_D13_gaps=gaps,C_high_pattern_count=758,C_successful_contexts=15,C_maps=30,
        complete_high_contexts=2614,complete_low_closures=66,
        scope='Independent canonical graphs/BFS, all264 pair weights and elementary triangle sets; saved B terminal-payment arithmetic and explicit immediate gaps. Complete physical local classifications and supporting-root arguments are inherited and mathematically read, not recounted.')


def interfaces(docs,data):
    p=docs[PRIMARY];base=docs[BASE]
    assert p['dimension_bound']==68 and p['inherited_dimension_bound']==base['dimension_bound']==docs[BASE+'_review']['dimension_bound']==69
    assert p['excluded_dimensions']==[69]and p['unrestricted_finite_volume']and not p['arithmetic_hypothesis']and p['cusp_restriction']is None
    local=local_interfaces(docs);finite=finite_interface(docs)
    costs={k:dict(v)for k,v in base['small_order_costs'].items()};costs['12']={'A':60,'B':65,'C':65,'D':62}
    assert p['small_order_costs']==costs
    assert all(F(max(v.values()),int(k))==budget(int(k)-1)for k,v in costs.items())
    assert p['ideal_rank_envelope']=={str(r):str(budget(r))for r in range(11,17)}
    c=p['concentration'];assert c['repeated_cap_ranks']==[10,11,13,14,15,16,17]
    assert c['at_most_one_low_part']==list(range(1,10))and c['at_most_one_each_middle_part']==[12]and c['at_most_one_high_part_minimum']==18
    assert F(c['double12_budget_increase'])==F(23,1092)
    assert p['finite_vertex_budget_at_excluded_dimension']=='544'
    assert forest_cost(p['finite_maximum_witness'])==544
    return dict(local=local,finite=finite,small_order_costs=costs,
        inherited_ceiling=69,newly_excluded_dimension=69,finite_slack='0',strict_global_slack=340,
        conclusion='The nonnegative finite term, positive ideal terms and strict global constant340 exclude69, including a polytope with no ideal vertices. The completed69 theorem excludes all higher dimensions.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--compare-only',action='store_true');p.add_argument('--verify',action='store_true');args=p.parse_args()
    assert sum((args.prepare,args.compare_only,args.verify))<=1
    data=arithmetic(not args.prepare)
    if args.prepare:
        path=Path('/private/tmp/dimension68_independent_scalar.json')
        path.write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')
        print('PASS_INDEPENDENT_DIMENSION68_SCALAR_PREPARATION',len(data['profiles']),data['dimension_minima'],flush=True)
        return
    stems=[PRIMARY,BASE,BASE+'_review',SUPPORT,A,B,C,D,A17,FOREST,FOREST+'_review',FINITE,TRI,COMBINED]
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in stems}
    data['composition_interfaces']=interfaces(docs,data)
    if args.compare_only:
        print('PASS_INDEPENDENT_DIMENSION68_COMPARISON',len(data['profiles']),'profiles;264 pair weights;228 partner comparisons;76 finite witnesses',flush=True)
        return
    paths=[Path(__file__)]
    for s in stems:
        doc=docs[s];assert doc['status'].startswith('PASS'),s
        for name,h in doc['sources'].items():assert digest(ROOT/name)==h,(s,name)
        paths.extend(ROOT/f'{folder}/{s}.{ext}'for folder,ext in[('scripts','py'),('results','json'),('results','md')])
    paths.extend(ROOT/f'results/{s}.md'for s in('fresh_cusp_link_shape_curvature_20261017','fresh_cusp_factor_payment_framework_20261017'))
    data.update(status='PASS_INDEPENDENT_DIMENSION68_SCALARS_AND_COMPOSITION',dimension_bound=68,
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        scope='No project arithmetic imports.2060 independent scalar profiles with three curvature identities;264 canonical order12 weights and triangle interfaces;228 A11 domination comparisons and76 finite attaining witnesses. All-label classifications, complete local geometry/corner maps and exhaustive Bellman upper bounds are inherited, not rerun. Full mathematical composition reviewed.')
    data=json.loads(json.dumps(data));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')
    print(data['status'],len(data['profiles']),'profiles; finite slack0/global340',flush=True)


if __name__=='__main__':main()
