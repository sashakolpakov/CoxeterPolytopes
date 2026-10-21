#!/usr/bin/env python3
"""Independent rational scalar reconstruction for the proposed ceiling65.

No project arithmetic imports. Preparation alone makes no global claim;
final bindings require the completed geometric and finite-budget packages.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb, prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_ideal5_weighted_dimension65_20261017'
TAG=PRIMARY+'_review'
BASE='fresh_ideal5_weighted_dimension66_20261017'
BC='fresh_order12_bc_unit_coupon_transfer_20261017'
B14='fresh_order14_affine_b_fork_projection_transfer_20261017'
LARGE_D='fresh_large_finite_d_leaf_half_transfer_20261017'
FINITE='fresh_a2_d11_full_coupon_forest_budget_20261017'
SUPPORT='fresh_order12_payment_support_20261017'
COMBINED='fresh_finite8_ideal5_combined_payments_20261017'
REPEATED=(10,11,12,14,15,16,17)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def cost(r):
    assert isinstance(r,int) and r>=1
    special={11:F(64,12),12:F(75,13),13:F(6),14:F(32,5),15:F(107,16),16:F(120,17)}
    return F(r,2) if r<=10 else special[r] if r<=16 else F(163,23)


def repeated(total,position=0):
    if position==len(REPEATED):
        if total==0:
            yield ()
        return
    rank=REPEATED[position]
    for count in range(total//rank+1):
        for tail in repeated(total-count*rank,position+1):
            yield (rank,)*count+tail


def profiles(total):
    """Optional parts first, independently of the primary's rank recursion."""
    seen=set()
    for low in range(10):
        for thirteen in (0,1):
            for high in (0,*range(18,total+1)):
                remain=total-low-13*thirteen-high
                if remain<0:
                    continue
                optional=((low,)if low else())+((13,)if thirteen else())+((high,)if high else())
                for tail in repeated(remain):
                    ranks=tuple(sorted(optional+tail))
                    assert ranks not in seen and sum(ranks)==total
                    assert sum(r<10 for r in ranks)<=1 and ranks.count(13)<=1 and sum(r>=18 for r in ranks)<=1
                    seen.add(ranks)
    return sorted(seen)


def scalar(ranks):
    h=[1]
    for rank in ranks:
        window=0;new=[]
        for j in range(len(h)+rank):
            if j<len(h):window+=h[j]
            if 0<=j-rank-1<len(h):window-=h[j-rank-1]
            new.append(window)
        h=new
    d=sum(ranks);v=prod(r+1 for r in ranks)
    assert len(h)==d+1 and h==h[::-1] and sum(h)==v
    g=sum(a*min(j,d-j)for j,a in enumerate(h))
    # Independent inclusion-exclusion cumulative-composition identity.
    shifts={0:1}
    for rank in ranks:
        out=dict(shifts)
        for shift,a in shifts.items():out[shift+rank+1]=out.get(shift+rank+1,0)-a
        shifts={k:a for k,a in out.items()if a}
    half=d//2;s=len(ranks)
    cdf=half*v-2*sum(a*(comb(half-j+s,s+1)if half-j+s>=s+1 else 0)for j,a in shifts.items())
    assert cdf==g
    curvature=F(5*g-d*v,2*v)
    absolute=F(sum(a*abs(2*j-d)for j,a in enumerate(h)),2*v)
    assert curvature==F(3*d,4)-F(5,2)*absolute
    charge=sum(map(cost,ranks),F())
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=v,h_polynomial=h,G=g,
                normalized_curvature=str(curvature),normalized_cost=str(charge),margin=str(curvature-charge))


def arithmetic():
    rows=[scalar(p)for p in profiles(65)]
    assert len(rows)==1659 and all(F(r['margin'])>0 for r in rows)
    best=min(rows,key=lambda r:F(r['margin']))
    assert(best['margin'],best['ranks'])==('13367/495924',[10,11,12,16,16])
    boundary_rows=[scalar(p)for p in profiles(64)]
    boundary=min(boundary_rows,key=lambda r:F(r['margin']))
    gain=cost(12)+cost(14)-2*cost(13)
    assert gain==F(11,65)>0
    assert all(cost(r)==F(r,2)for r in range(1,11))
    assert all(cost(r)==F(163,23)for r in range(17,76))
    return dict(reduced_profiles=rows,dimension_minima=[dict(ambient_dimension=66,profiles=len(rows),
        minimum_margin=best['margin'],minimizing_ranks=best['ranks'],
        finite_maximum='495',finite_curvature_coefficient='495',finite_slack='0',
        finite_maximum_witness=['D66'],strict_global_slack=325)],
        boundary_at_n65=boundary,boundary_terminal_profiles=len(boundary_rows),
        boundary_negative_profiles=sum(F(r['margin'])<0 for r in boundary_rows),
        duplicate_rank13_cost_increase=str(gain),
        concentration=dict(repeated=list(REPEATED),at_most_one_low=list(range(1,10)),
                           at_most_one_rank13=True,at_most_one_high_lower_bound=18),
        formulas=['sliding-window uniform-polynomial product',
            'G=floor(d/2)*V-2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1)',
            'kappa/V=3*d/4-(5/2)*E(abs(sum U_i-d/2))'])


def affine_graph(kind,n):
    g=[[0]*n for _ in range(n)]
    if kind=='B':
        es=[(i,i+1,4 if i==n-4 else 3)for i in range(n-3)]+[(0,n-2,3),(0,n-1,3)]
    else:
        assert kind=='C'
        es=[(i,i+1,4 if i in(0,n-2)else 3)for i in range(n-1)]
    for i,j,m in es:g[i][j]=g[j][i]=m
    return g


def distances(g):
    rows=[]
    for root in range(len(g)):
        dd={root:0};queue=[root]
        for i in queue:
            for j,m in enumerate(g[i]):
                if m and j not in dd:dd[j]=dd[i]+1;queue.append(j)
        assert len(dd)==len(g);rows.append(dd)
    return rows


def basic(g):
    ds=distances(g)
    features=[i for i,row in enumerate(g)if sum(bool(v)for v in row)>=3 or any(v>=4 for v in row)]
    return {p for p in combinations(range(len(g)),2)if ds[p[0]][p[1]]<=5 or
            ds[p[0]][p[1]]<=9 and min(ds[i][f]for i in p for f in features)<=5}


def local_interfaces(bc,b14,support,combined):
    pair_records=[]
    for kind,wanted in(('B',64),('C',63)):
        saved=bc[kind];g=affine_graph(kind,12);rec=next(z for z in support['records']if z['kind']==kind)
        assert g==rec['target']
        paid=set(map(tuple,saved['retained_unit_pairs']));ds=distances(g)
        assert paid==basic(g)==set(map(tuple,rec['basic_pairs']))and len(paid)==wanted
        assert {p for p in combinations(range(12),2)if ds[p[0]][p[1]]<=8}<=paid
        removed=set(map(tuple,saved['removed_half_pairs']if kind=='B'else saved['removed_unit_pairs']))
        assert all(ds[i][j]==10 for i,j in removed)
        maps=0
        patterns={z['index']:z for z in support['high_patterns']}
        for z in bc['reconstructed_high_support'][kind]:
            pattern=patterns[z['index']];p=pattern['labels'];marks=pattern['marks'];images=set()
            for mp in z['embeddings']:
                assert len(mp)==len(set(mp))==len(p)
                assert all(p[i][j]==g[mp[i]][mp[j]]for i in range(len(p))for j in range(len(p)))
                images.add(tuple(sorted(mp[i]for i in marks)));maps+=1
            assert images==set(map(tuple,z['marked_images']))
        pair_records.append(dict(kind=kind,order=12,basic_unit_pairs=len(paid),
            removed_pairs=sorted(removed),saved_induced_maps_checked=maps,
            elementary_triangle_scope='Retained entire basic ideal5/feature5 rule; inherited elementary triangle theorem.'))
    variants=[];upgrades=Counter()
    for case in bc['B']['cases']:
        for row in case['payment_variants']:
            total=F()
            for c in row['corners']:
                if c['kind']!='finite':
                    assert c['kind']=='affineB11_order12'and F(c['new_weight'])==0
                    continue
                left,right=[next(z for z in c['components']if m in z['nodes'])for m in c['marks']]
                unit=int(left is right)
                assert unit==c['finite_unit']and(not unit or len(left['nodes'])<=13)
                coupon=F(c['old_half']);reason=None
                for small,large,mark in((left,right,c['marks'][1]),(right,left,c['marks'][0])):
                    if small['type']in('A1','A2')and large['type']=='D11'and large['roles'].get(str(mark))=='short_fork_leaf':
                        assert not unit and coupon==F(1,2)
                        coupon=F(1);reason='existing_A1_upgrade'if small['type']=='A1'else'new_A2_upgrade'
                assert coupon==F(c['new_coupon'])and F(c['new_weight'])==unit+coupon
                assert(c['upgrade']is None)==(reason is None)
                if reason:assert c['upgrade']['source']==reason;upgrades[reason]+=1
                total+=unit+coupon
            assert total==F(row['new_finite_total'])and total>=1 and F(row['new_ideal_total'])==0
            variants.append(dict(case=case['index'],variant=row['variant'],finite_total=str(total)))
    assert len(variants)==16 and Counter(z['finite_total']for z in variants)=={'1':8,'2':8}
    assert dict(upgrades)==bc['used_upgrades']==dict(new_A2_upgrade=12,existing_A1_upgrade=4)
    cprojections=[]
    assert [r['index']for r in bc['C']['cases']]==[164,1452]
    for case in bc['C']['cases']:
        assert case['retained_finite_unit']==1
        for row in case['maps']:
            assert len(row['forced_entries'])==1
            q=row['forced_entries'][0]
            value=F(q['norm_two_contact_square'])*F(q['old_gauge_square'])*F(q['inverse_diagonal'])
            assert q['native_label']==4 and value==F(q['forced_projection_lower_bound'])==16>2
            cprojections.append(dict(case=case['index'],projection=str(value)))
    assert len(cprojections)==4
    g=affine_graph('B',14);paid=set(map(tuple,b14['retained_pairs']))
    assert g==b14['target']and len(paid)==84 and basic(g)<=paid
    assert set(map(tuple,b14['basic_pairs_retained']))==basic(g)
    triangle=next(z for z in combined['triangle_pair_retention']if z['kind']=='B'and z['order']==14)
    assert set(map(tuple,triangle['required_pairs']))<=paid
    assert set(map(tuple,b14['triangle_pairs_retained']))==set(map(tuple,triangle['required_pairs']))
    assert set(map(tuple,b14['removed_pairs']))=={(10,12),(10,13)}
    assert all(distances(g)[i][j]==11 for i,j in b14['removed_pairs'])
    assert bc['low_actual_cases']==66 and bc['low_actual_variants']==69
    assert not bc['large_D_D10_coupon_used']and bc['existing_A1_D11_full_upgrade_required']
    assert bc['new_coupon_increment_per_A2_D11_component_pair']=='2'
    assert bc['doubled_finite_budget_increment']=='4*a2*d11'
    assert b14['finite_weights_unchanged']and not b14['added_finite_coupons']
    return dict(order12_pair_rules=pair_records,B_finite_variants=variants,
        C_forced_terminal_projections=cprojections,upgrades=dict(upgrades),
        B14_triangle_pairs=len(triangle['required_pairs']),B14_retained_unit_pairs=84,
        saved_map_scope='Only supplied maps are checked; upstream complete search and arbitrary-support sealing remain bound geometric inputs.')


def finite_certificate(data):
    """Independent count recursion and original-component coupon recount."""
    limit=data['maximum_rank'];assert limit==75
    ranks=tuple(range(1,10))+tuple(range(17,limit+1))
    def component(rank):return('A'if rank<=9 else'D',rank)
    def own(t):
        k,r=t
        if r<=16:return r*(r-1)
        if r==17:return 228 if k=='A'else 254
        assert k=='D'
        return 16*r-(66 if r>=20 else 64)
    def pair(a,b):
        out=0
        for small,large in((a,b),(b,a)):
            if small[1]>5:continue
            k,r=large
            if(k,r)==('D',10):out+=2*small[1]
            elif k=='D'and 11<=r<=16:
                out+=3*small[1]
                if small==('A',1)and r in(11,12):out+=2
                if small==('A',2)and r==11:out+=4
                if small==('A',4)and r==13:out+=2
            elif(k,r)==('A',17):out+=2*small[1]
        return out
    expected={}
    def visit(rank,left,counts):
        if rank==18:
            ts=[('A'if r==17 else'D',r)for r,n in zip(range(10,18),counts)for _ in range(n)]
            expected[tuple(counts)]=ts;return
        for n in range(left//rank+1):visit(rank+1,left-rank*n,counts+[n])
    visit(10,limit,[]);assert len(expected)==1608
    cert=data['certificate'];domain=cert['large_count_vectors']
    assert len(domain)==1608 and {tuple(z['component_counts'])for z in domain}==set(expected)
    tables={tuple(z['interaction_key']):z for z in cert['residual_tables']}
    assert len(tables)==len(cert['residual_tables'])==357
    costs={};domains=[]
    for row in domain:
        cc=row['component_counts'];ts=expected[tuple(cc)];used=sum(r for k,r in ts);base=sum(own(t)for t in ts)
        key=(2*cc[0]+3*sum(cc[1:7])+2*cc[7],cc[1]+cc[2],cc[3],cc[1])
        assert row['used_rank']==used and row['doubled_base_cost']==base and tuple(row['interaction_key'])==key
        actual={r:own(component(r))+sum(pair(component(r),t)for t in ts)for r in ranks}
        if key in costs:assert costs[key]==actual
        costs[key]=actual;domains.append((used,base,key))
    assert set(costs)==set(tables)
    alternatives=0;attainers=0
    for key,table in tables.items():
        values=table['doubled_optima'];last=table['last_rank']
        assert len(values)==len(last)==limit+1 and values[0]==last[0]==0
        for n in range(1,limit+1):
            for r in ranks:
                if r<=n:assert values[n]>=values[n-r]+costs[key][r];alternatives+=1
            r=last[n];assert r in ranks and r<=n and values[n]==values[n-r]+costs[key][r];attainers+=1
    global_checks=0
    for n,row in enumerate(cert['dimension_maxima']):
        assert row['dimension']==n
        for used,base,key in domains:
            if used<=n:assert row['doubled_cost']>=base+tables[key]['doubled_optima'][n-used];global_checks+=1
        ts=[(name[0],int(name[1:]))for name in row['witness_components']]
        assert all(t==component(t[1])or t in[('D',r)for r in range(10,17)]+[('A',17)]for t in ts)
        assert sum(r for k,r in ts)==n
        value=sum(map(own,ts))+sum(pair(a,b)for a,b in combinations(ts,2))
        assert value==row['doubled_cost']
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert F(row['finite_curvature_coefficient'])==beta and F(row['finite_slack'])==beta-F(value,2)
    assert len(cert['dimension_maxima'])==76
    assert(alternatives,attainers,global_checks)==(860013,26775,22504)
    f=cert['dimension_maxima'][66];assert f['doubled_cost']==990 and f['witness_components']==['D66']
    # Direct finite D66 graph, with the two short-leaf-nine units halved.
    n=66;g=[[0]*n for _ in range(n)]
    for a,b in[(i,i+1)for i in range(n-3)]+[(n-3,n-2),(n-3,n-1)]:g[a][b]=g[b][a]=3
    ds=distances(g);leaves={i for i,row in enumerate(g)if sum(bool(v)for v in row)==1}
    total=F();short=leaf9=halves=0
    for i,j in combinations(range(n),2):
        if ds[i][j]<=8:total+=1;short+=1
        elif ds[i][j]==9 and {i,j}&leaves:
            leaf9+=1;half=bool({i,j}&{n-2,n-1});halves+=half;total+=F(1,2)if half else 1
    assert(short,leaf9,halves,total)==(493,3,2,495)
    return dict(large_count_vectors=1608,interaction_tables=357,
        bellman_upper_alternatives=alternatives,attaining_predecessors=attainers,
        global_upper_alternatives=global_checks,typed_attainers=76,
        direct_D66=dict(short_pairs=493,leaf9_pairs=3,halved_pairs=2,cost='495'),
        type_domination_scope='Full source/math read; finite classification and same-rank domination remain inherited inputs.')


def main(prepare=False,compare_only=False,verify=False):
    result=arithmetic()
    if prepare:
        path=Path('/private/tmp/dimension65_independent_scalar_20261017.json')
        path.write_text(json.dumps(result,separators=(',',':'))+'\n')
        print('PASS_CONDITIONAL_DIMENSION65_SCALARS',result['dimension_minima'],
              'boundary',result['boundary_at_n65']['ranks'],result['boundary_at_n65']['margin'],
              'boundary_profiles',result['boundary_terminal_profiles'],flush=True)
        return
    stems=(PRIMARY,BASE,BASE+'_review',BC,B14,LARGE_D,FINITE,SUPPORT,COMBINED)
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in stems}
    p=docs[PRIMARY]
    assert p['dimension_bound']==65 and p['inherited_dimension_bound']==66 and p['excluded_dimensions']==[66]
    assert p['unrestricted_finite_volume']and not p['arithmetic_hypothesis']and p['cusp_restriction']is None
    assert docs[BASE]['dimension_bound']==docs[BASE+'_review']['dimension_bound']==66
    assert {tuple(r['ranks']):r for r in p['reduced_profiles']}=={tuple(r['ranks']):r for r in result['reduced_profiles']}
    assert len(p['reduced_profiles'])==1659 and p['dimension_minima']==result['dimension_minima']
    assert p['boundary_at_n65']==dict(profile=result['boundary_at_n65'],profiles=1528,negative_profiles=16,
        finite_maximum='487',finite_curvature_coefficient='480',finite_slack='-7',
        finite_maximum_witness=['D65'],new_bound64_asserted=False)
    assert result['boundary_terminal_profiles']==1528 and result['boundary_negative_profiles']==16
    assert p['concentration']['repeated_cap_ranks']==list(REPEATED)
    assert p['concentration']['at_most_one_low_part']==list(range(1,10))
    assert p['concentration']['at_most_one_each_middle_part']==[13]
    assert p['concentration']['at_most_one_high_part_minimum']==18
    assert F(p['concentration']['double13_budget_increase'])==F(11,65)
    assert F(p['concentration']['invalid_double12_budget_increase'])==cost(11)+cost(13)-2*cost(12)==-F(8,39)
    assert docs[BC]['new_order12_costs']==dict(A=60,B=64,C=63,D=62)
    assert docs[B14]['new_order14_costs']==dict(A=70,B=84,C=84,D=82)
    small=dict(docs[BASE]['small_order_costs']);small['12']=docs[BC]['new_order12_costs'];small['14']=docs[B14]['new_order14_costs']
    assert p['small_order_costs']==small
    assert p['ideal_rank_envelope']=={str(r):str(cost(r))for r in range(11,17)}
    assert all(F(max(small[str(r+1)].values()),r+1)==cost(r)for r in range(11,17))
    result['local_interfaces']=local_interfaces(docs[BC],docs[B14],docs[SUPPORT],docs[COMBINED])
    result['finite_certificate_check']=finite_certificate(docs[FINITE])
    assert docs[LARGE_D]['routing']==dict(whole_component_at_most17=2531,adaptive_simple_leaf_neighbor_degree2=9,
        protected_simple_leaf_neighbor_degree2=18,two_protected_leaf9_corners=40,new_small_component_D10_half=16)
    assert all(docs[s]['triangle_minimum']==2 and docs[s]['quadrilateral_minimum']==1 for s in(BC,LARGE_D))
    result['composition_review']=dict(
        finite_budget_hypothesis_supplied_by=BC,
        preserved_interfaces=['elementary triangles and short8','all2614 high-overlap cases and separate low/disjoint routes',
            'universal and later cycle/tree replacements','joint large-D/D10 and A2/D11 coupon budget',
            'actual cusp fiber multiplicity','corrected duplicate13 concentration','completed66 higher tail'],
        global_contradiction='Finite slack0, every ideal margin positive, strict additive constant325.',
        inherited_scope='Original geometric censuses, arbitrary-support sealing and all-rank convex-order/fiber/shape theorems are reviewed bound inputs, not independently rerun here.')
    paths=[Path(__file__)]
    for stem,doc in docs.items():
        assert doc['status'].startswith('PASS'),stem
        for rel,h in doc['sources'].items():assert digest(ROOT/rel)==h,(stem,rel)
        paths.extend(ROOT/f'{folder}/{stem}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths.extend(ROOT/f'results/{stem}.md'for stem in('fresh_cusp_factor_payment_framework_20261017',
        'fresh_cusp_link_shape_curvature_20261017','fresh_distance14_unrestricted_dimension156_20261017'))
    if compare_only:
        print('PASS_INDEPENDENT_DIMENSION65_COMPARISON',len(result['reduced_profiles']),result['finite_certificate_check'],flush=True)
        return
    result.update(status='PASS_INDEPENDENT_DIMENSION65_SCALARS_AND_INTERFACES',dimension_bound=65,
        sources={str(q.relative_to(ROOT)):digest(q)for q in paths},
        scope='No producer arithmetic imports: all1659 primary scalar rows and1528 boundary profiles by optional-first enumeration, sliding/CDF/absolute-deviation arithmetic; saved local pair/map/finite-payment interfaces; independent1608-vector/357-table finite upper-and-attainment certificate check. Original geometric census and all-rank structural lemmas remain separately reviewed inputs.')
    result=json.loads(json.dumps(result));path=ROOT/f'results/{TAG}.json'
    if verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(result['status'],result['dimension_minima'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--prepare',action='store_true');ap.add_argument('--compare-only',action='store_true');ap.add_argument('--verify',action='store_true')
    args=ap.parse_args();assert sum((args.prepare,args.compare_only,args.verify))<=1
    main(args.prepare,args.compare_only,args.verify)
