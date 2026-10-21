#!/usr/bin/env python3
"""Residual cost below one improves the reviewed165 charge bound to164.

Uses frozen independently reviewed grounded arithmetic as a read-only helper.
All newly required thresholds, inequalities and small-cost cases are recomputed.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_nikulin_height_charge_20261015_review as H

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_defect_frontier_20261015'
OLD='fresh_nikulin_height_charge_20261015'
B=F(21,20)

def bound(a2,b2):return max(min(a2,b2),max(a2,b2)/B**2)
def literal_models():
    return {'A1':([[2,-2],[-2,2]],[1,1]),'E6':(H.star([2,2,2]),[1]*7),
      'E7':(H.star([1,3,3]),[1]*8),'E8':(H.star([1,2,5]),[1]*9),
      'F4':([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-2,0],[0,0,-2,4,-2],[0,0,0,-2,4]],[1,1,1,2,2]),
      'G2':([[2,-1,0],[-1,2,-3],[0,-3,6]],[1,1,3])}
def single(ds,R,om,side,k):
    i,j=om[side],om[1-side]
    co=[R(i,j),-R(i,j)+R(k,j)-R(i,k),R(i,k)]
    minimum=H.minimum(co)
    height0=bound(ds[k],ds[j])
    height1=bound(ds[i]+ds[k]+2*H.sqrtq(ds[i]*ds[k]),ds[j])
    return height0*R(k,j),height1*minimum

def widened_window(old):
    records=[];singles=0;special=0
    for row in old['exceptional_cases']:
        ds,diags,R,cross=H.model(*literal_models()[row['type']]);T=F(row['target']);om=row['omissions'];out=[]
        for c in row['single_cases']:
            side,k=c['side'],c['node'];cost=diags[om[side]][k]
            assert cost==F(c['cost'])
            lo=single(ds,R,om,side,k);assert all(z+cost>=T for z in lo)
            out.append(dict(side=side,node=k,cost=str(cost),energy_lower=list(map(str,lo)),unit_margins=[str(z+cost-T)for z in lo]))
        for c in row['multiple_cases']:
            if c['status']=='UNIT_CHARGE_VIA_E7_HEIGHT':
                assert row['type']=='E7'and row['forbidden']is None
                assert min(F(5,2),F(9,2)/B**2+1)>T;special+=1
            else:assert F(c['cost'])>=F(3,2)and F(c['defect_upper'])<=F(2,9)
        singles+=len(out);records.append(dict(type=row['type'],forbidden=row['forbidden'],singletons=out))
    assert singles==71 and special==1
    cases=[]
    for row in old['classical_A']['small_cycles']:
        T=F(row['target'])
        for c in row['cases']:
            h=bound(F(len(c['left_support'])**2),F(len(c['right_support'])**2))
            cost=F(c['left_min_cost'])+F(c['right_min_cost'])
            lo=h*F(c['convex_dual']['minimum_energy'])+cost
            assert lo>=T
            cases.append(dict(cycle_size=row['cycle_size'],left=c['left_support'],right=c['right_support'],unit_lower=str(lo),margin=str(lo-T)))
    assert len(cases)==580
    return dict(exceptional_records=records,cycle_cases=cases,unchanged='Complete support sets, KKT certificates, zero-cost bounds, A large-rank tails and BCD estimates inherited unchanged from frozen independently reviewed165 proof.')

def cheap_exceptions(old):
    result=[]
    for row in old['exceptional_cases']:
        name=row['type'];ds,diags,R,cross=H.model(*literal_models()[name]);om=row['omissions'];T=F(row['target'])
        # The entire residual product has fewer than two retained coordinates;
        # this list includes every possible position of its unique coordinate.
        allowed=[(side,k,c)for side,i in enumerate(om)for k,c in diags[i].items()if c<1]
        cases=[]
        for side,k,cost in allowed:
            if T<=F(2,3):reason='target<=2/3';lo=F()
            else:
                bounds=single(ds,R,om,side,k);lo=min(bounds);assert lo>=T
                reason='W>=target by complete singleton energy/height bound'
            assert T-lo<=F(2,3)
            cases.append(dict(side=side,node=k,diagonal=str(cost),reason=reason,energy_lower=str(lo),loss_upper=str(T-lo)))
        result.append(dict(type=name,forbidden=row['forbidden'],target=str(T),all_diagonals_below_one=cases))
    assert len(result)==40
    return result

# Tiny rational polynomial operations, ascending coefficients.
def add(a,b):return[(a[i]if i<len(a)else F())+(b[i]if i<len(b)else F())for i in range(max(len(a),len(b)))]
def mul(a,b):
    c=[F()]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c
def scale(a,b):return[x*b for x in a]
def cheap_A():
    out=[]
    for start,eps in((2,0),(3,1)):
        m=[F(start),F(2)]
        # Point row: R(h-1)-2(m-1)/9+4/9.
        num=add(add(mul(m,m),scale(m,24)),[-36 if not eps else-81]);den=scale(m,36)
        assert all(x>=0 for x in num)and all(x>=0 for x in den)and den[0]>0
        out.append(dict(parity=eps,cycle_size_coefficients=list(map(str,m)),point_loss_at_most='4/9',
            cleared_margin_numerator=list(map(str,num)),positive_denominator=list(map(str,den))))
    # If omitted entry is positive, m2/3 have target<=4/9. For m>=4,
    # E_min >= (m-3)/4, and height squared >=4/B².
    # Exact differences between E_min and (m-3)/4 are positive:
    # even (2m-3)/(4(m-1)); odd (m-3)/(2(m-1)).
    for m,eps in((4,0),(5,1)):
        exact=F(m*(m-2),4*(m-1))if not eps else F((m+1)*(m-3),4*(m-1))
        difference=F(2*m-3,4*(m-1))if not eps else F(m-3,2*(m-1))
        assert exact-F(m-3,4)==difference>0
    margin=[F(106,441),F(302,441)] # m=4+t
    assert (F(4-3)/B**2-F(2*(4-1),9))==margin[0]
    assert 1/B**2-F(2,9)==margin[1]>0
    return dict(all_diagonals='d(m-d)/m<1 permits only d=1; d>=2 has minimum2(m-2)/m>=1 for m>=4.',
        point_rows=out,positive_omitted_energy_min_even='m*(m-2)/(4*(m-1))',
        positive_omitted_energy_min_odd='(m+1)*(m-3)/(4*(m-1))',
        positive_omitted_uniform_energy='(m-3)/4 for m>=4',
        unit_energy_margin_at_m4_plus_t=list(map(str,margin)),
        conclusion='W>=T-4/9; hence residual loss<=2/3, in every ordinary/avoid-node A rank.')

def cheap_BCD():
    # Every low-diagonal possibility from the complete frozen24 prefix classes.
    # Remaining classes have integer diagonal>=1 or rank/2>=1.
    rows=[
      dict(type='B',modified=False,side='r',nodes='0 or1',cost='r/4',rank_range=[2,3],max_target='2/3',reason='small target'),
      dict(type='B',modified=True,side='r-1',nodes='0 or1',cost='(r-1)/4',rank_range=[3,4],max_target='2/3',reason='small target'),
      dict(type='B',modified=True,side='r-1',nodes='r',cost='1/2',rank_range='all r>=3',reason='linear loss0; point energy>=target'),
      dict(type='B',modified=True,side='separateB2',nodes='all',rank_range=[2],max_target='2/9',reason='small target'),
      dict(type='C',modified=True,side='1',nodes='0',cost='1/2',rank_range='all r>=2',reason='linear loss0; point energy>=target'),
      dict(type='C',modified=True,side='1',nodes='r',cost='(r-1)/2',rank_range=[2],max_target='2/9',reason='small target')]
    # Recompute all24 independent affine-polynomial identities. Their coverage
    # was proved/reviewed in the frozen charge package; this class split is new.
    all_classes=H.classical_BCD()
    # Necessary finite rank bounds are exact: r/4<1 givesr<=3;
    # (r-1)/4<1 givesr<=4; (r-1)/2<1 givesr<=2.
    assert max(F(2*r,9)for r in(2,3))==F(2,3)
    assert max(F(2*(r-1),9)for r in(3,4))==F(2,3)
    assert F(2*(2-1),9)==F(2,9)
    # Audit every coordinate class at the cone's coefficient level. A cost
    # polynomial with constant>=1 admits no low-diagonal point. Otherwise
    # either its loss polynomial is zero or its exact parameter inequality
    # forces target<=2/3. Finite leftover integer parameters are enumerated
    # only when a positive coefficient bounds them; zero-coefficient variables
    # cannot change the corresponding target (checked below).
    coverage=[]
    for c in all_classes:
        cost=list(map(F,c['cost_coefficients']));loss=list(map(F,c['loss_coefficients']));rank=list(map(F,c['rank_coefficients']))
        target=scale(rank,F(2,9));target[0]-=F(2*int(c['modified']),9)
        if cost[0]>=1:reason='diagonal>=1'
        elif all(z==0 for z in loss):reason='zero linear loss'
        else:
            assert all(z>=0 for z in cost)
            domains=[]
            for a,b in zip(cost[1:],target[1:]):
                if not a:assert b==0;domains.append([0])
                else:domains.append([k for k in range(10)if cost[0]+a*k<1]);assert cost[0]+a*10>=1
            targets=[]
            for params in product(*domains):
                if cost[0]+sum(a*k for a,k in zip(cost[1:],params))<1:
                    val=target[0]+sum(a*k for a,k in zip(target[1:],params));assert val<=F(2,3);targets.append(val)
            assert targets;reason='low diagonal forces target<=2/3'
        coverage.append(dict(type=c['type'],modified=c['modified'],omissions=c['omissions'],node=c['node'],reason=reason))
    return dict(low_diagonal_classes=rows,all24_class_coverage=coverage,conclusion='Every single-contact factor with cost<1 has T-W<=2/3, including all modified choices.')

def sharp_residual_example():
    # Affine B3 has right-omission inverse diagonal3/4 atnode0.
    # Both rows atnode0 have the same normalized potential, so E=W=0.
    A=[[2,0,-1,0],[0,2,-1,0],[-1,-1,2,-2],[0,0,-2,4]]
    ds,diags,R,cross=H.model(A,[1,1,1,2])
    assert ds[0]==1 and diags[3][0]==F(3,4)and R(0,0)==0
    return dict(type='affineB3',omissions=[0,3],both_rows_node=0,contacts=[1,1],
        cost='3/4',height_ratio=1,energy=0,target='2/3',loss='2/3')

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old=json.loads((ROOT/f'results/{OLD}.json').read_text())
    review=json.loads((ROOT/f'results/{OLD}_review.json').read_text())
    assert review['status']=='PASS_INDEPENDENT_EXACT_HEIGHT_CHARGE_DEFECT_BUDGET_AND_165_REVIEW'
    window=widened_window(old);ex=cheap_exceptions(old);a=cheap_A();bcd=cheap_BCD();sharp=sharp_residual_example()
    cases=[]
    # With<=1 defective factor, old unit charge+one2/9defect suffices.
    cases.append(dict(defective_factors='0 or1',loss_bound=str(4+F(2,9)),proof='total cost4 plus at most one2/9defect'))
    for e8 in range(3):
        minimum=e8*F(14,9)+(2-e8)*F(3,2)
        residual=4-minimum;target=e8*F(16,9)+(2-e8)*F(14,9)
        if e8:
            assert residual<1;bound=target+F(2,3);reason='subunit residual-product lemma'
        else:bound=target+residual;reason='remaining factors satisfy unit charge'
        assert bound<=F(38,9)
        cases.append(dict(defective_factors=2,ordinary_E8_count=e8,minimum_bad_cost=str(minimum),
            residual_cost_budget=str(residual),bad_target_sum=str(target),loss_bound=str(bound),proof=reason))
    near=(F(2*164,9)-F(2,9)-F(38,9))/2-2
    far=B*F(164-43,9)
    assert near==14 and far==14+F(7,60)
    paths=[f'scripts/{TAG}.py']
    for stem in(OLD,OLD+'_review','fresh_nikulin_unit_charge_20261015_dimension169',
                 'fresh_nikulin_coupled_charge_20261015_dimension166'):
        paths += [f'scripts/{stem}.py',f'results/{stem}.json']
    for p in paths:
        if p.endswith('.json'):
            for q,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():assert sha(q)==h,(p,q)
    return json.loads(json.dumps(dict(status='PASS_EXACT_RESIDUAL_COST_ENVELOPE_CANDIDATE164_PENDING_INDEPENDENT_REVIEW',
        sources={p:sha(p)for p in paths},height_ratio_window=str(B),widened_165_checks=window,
        residual_exceptional_cases=ex,residual_A_all_ranks=a,residual_BCD_all_ranks=bcd,
        sharp_single_contact_residual_example=sharp,
        inverse_diagonal_lemma='For a norm-two positive definite principal Gram, J_kk=1/(2-v^T B^-1 v)>=1/2. Entrywise inverse positivity and contact gap imply global cost<1 allows at most one retained coordinate.',
        aggregate_cases=cases,total_loss_at_most='38/9',near_pairing_lower='M*(N-38)/9',far_pairing_lower='M*(N-43)/9',
        threshold_checks=dict(N=164,near_lower=str(near),far_lower=str(far),near_increment='t/9',far_increment='7t/60',
            ideal_branch_bound=164,compact_facet_branch_inherited_bound=143,combined_candidate_bound=164),
        scope='New candidate164 from complete residual-cost envelope and reviewed165 chain. All affine ranks and preset avoid-node choices covered; arbitrary real Coxeter-gap contacts. Independent review pending. No novelty/current-best or formal verification claim.')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    data=run();out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
