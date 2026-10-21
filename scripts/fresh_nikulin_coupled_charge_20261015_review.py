#!/usr/bin/env python3
"""Independent grounded-energy review of the ratio-dependent7/6 estimate."""
import argparse,json
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
import fresh_nikulin_unit_charge_20261015_review as U

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_coupled_charge_20261015_review'
PRIMARY='fresh_nikulin_coupled_charge_20261015'
BASE='fresh_nikulin_unit_charge_20261015'
ALPHA=Q(7,6);B=Q(21,20)

def height_from_squares(a,b):return max(min(a,b),max(a,b)/B**2)

def run():
    saved=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    old=json.loads((ROOT/f'results/{BASE}.json').read_text())
    reviewed=json.loads((ROOT/f'results/{BASE}_review.json').read_text())
    literals={
      'A1':([[2,-2],[-2,2]],[Q(1)]*2),
      'E6':(U.H.star([2,2,2]),[Q(1)]*7),'E7':(U.H.star([1,3,3]),[Q(1)]*8),
      'E8':(U.H.star([1,2,5]),[Q(1)]*9),
      'F4':([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-2,0],[0,0,-2,4,-2],[0,0,0,-2,4]],
            [Q(1),Q(1),Q(1),Q(2),Q(2)]),
      'G2':([[2,-1,0],[-1,2,-3],[0,-3,6]],[Q(1),Q(1),Q(3)])}
    models={name:U.model(*data) for name,data in literals.items()}
    records=[];single_count=two_count=0
    for previous in old['exceptional_cases']:
        name=previous['type'];f=previous['forbidden'];omits=previous['omissions']
        row=next(x for x in saved['exceptional_cases'] if x['type']==name and x['forbidden']==f)
        for key in ('type','rank','forbidden','omissions','target','gram','kernel'):
            assert row[key]==previous[key]
        delta2,diags,R,cross=models[name];target=Q(row['target']);limit=target/ALPHA
        U.physical_gram(*literals[name],row['gram'])
        assert [{str(k):str(v) for k,v in diags[i].items()} for i in omits]==row['inverse_diagonals']
        assert min(delta2[i] for i in omits)*R(*omits)>=target
        cheap=[(side,k,c) for side,i in enumerate(omits) for k,c in diags[i].items() if c<limit]
        assert [(x['side'],x['node'],Q(x['cost'])) for x in row['single_contact_cases']]==cheap
        singles=[]
        for side,k,cost in cheap:
            i,j=omits[side],omits[1-side]
            triple=[R(i,j),R(k,j),R(i,k)]
            coeff=[triple[0],-triple[0]+triple[1]-triple[2],triple[2]]
            minimum=U.minimum_quadratic(coeff)
            hzero=height_from_squares(delta2[k],delta2[j])
            positive_square=s.simplify(delta2[i]+delta2[k]+2*s.sqrt(delta2[i]*delta2[k]))
            hpositive=height_from_squares(positive_square,s.Rational(delta2[j]))
            lower0=hzero*triple[1]+ALPHA*cost
            lower1=s.simplify(hpositive*minimum+ALPHA*cost)
            assert lower0>=target and lower1>=target
            p=next(x for x in row['single_contact_cases'] if (x['side'],x['node'])==(side,k))
            assert list(map(Q,p['resistances']))==triple and Q(p['minimum_energy'])==minimum
            for key,value in [('zero_omitted_height_square',hzero),('positive_omitted_height_square',hpositive),
                              ('zero_omitted_lower',lower0),('positive_omitted_lower',lower1)]:
                assert s.simplify(s.sympify(p[key])-value)==0
            singles.append(dict(side=side,node=k,energy_coefficients=list(map(str,coeff)),
                                height_lower_squares=list(map(str,(hzero,hpositive))),lower_bounds=list(map(str,(lower0,lower1)))))
        surviving=[]
        for size in range(2,len(cheap)+1):
            for support in combinations(cheap,size):
                cost=sum((c for _,_,c in support),Q(0))
                for (a,k,_),(b,l,_) in combinations(support,2):
                    if a==b:cost+=2*cross(omits[a],k,l)
                if cost<limit:surviving.append((support,cost))
        if surviving:
            assert name=='E7' and f is None
            assert surviving==[(((0,4,Q(1,2)),(1,7,Q(1,2))),Q(1))]
            assert Q(1,2)+ALPHA==Q(5,3)>target
            # This identical support and exact polynomial were independently
            # polarized in the frozen4/3 review; no new two-contact identity.
            assert next(x for x in reviewed['exceptional_cases'] if x['type']=='E7' and x['forbidden'] is None)['nontrivial_multiple_supports']==1
            two_count+=1
        assert bool(surviving)==(row['two_contact_case'] is not None)
        single_count+=len(singles)
        records.append(dict(type=name,forbidden=f,omissions=omits,target=str(target),
                            single_contacts=singles,surviving_multiple_supports=len(surviving)))
    assert len(records)==40 and single_count==54 and two_count==1

    # Derive small-cycle resistances from an independent inverse of the path
    # obtained by grounding one node, rather than using the producer formula.
    cycles=[]
    for m in range(2,20):
        J=U.H.inv([[Q(2 if i==j else -1 if abs(i-j)==1 else 0) for j in range(m-1)] for i in range(m-1)])
        R=[Q(0)]+[J[k][k] for k in range(m-1)]
        assert R==[Q(k*(m-k),m) for k in range(m)]
        h=m//2;target=Q(2*(m-1),9);limit=target/ALPHA
        possible=[k for k in range(1,m) if R[k]<=2 and R[k]<limit]
        losses={k:R[h]+R[k]-R[(k-h)%m] for k in possible}
        lam=max([ALPHA]+[losses[k]/R[k] for k in possible])
        margin=R[h]-target-(lam-ALPHA)*limit;assert margin>=0
        oldrow=next(x for x in saved['classical_A']['finite_cases'] if x['cycle_size']==m)
        assert [x['node'] for x in oldrow['allowed_retained_nodes']]==possible
        assert Q(oldrow['uniform_margin'])==margin and Q(oldrow['maximum_loss_per_cost'])==lam
        cycles.append(dict(order=m,allowed_nodes=possible,loss_per_cost=str(lam),point_surplus_margin=str(margin)))
    tail=[]
    for m,eps in ((20+2*U.t,0),(21+2*U.t,1)):
        assert U.nonnegative(m-9) and (m-9).subs(U.t,0)>0
        for d in (1,2):
            margin=s.expand(d*(m-7*d-6*eps));assert U.nonnegative(margin)
            tail.append(dict(order=str(m),distance=d,cleared_margin=str(margin)))

    bcd=[];variables={str(x):x for x in (U.t,U.u,U.v)}
    for row in reviewed['all_rank_coordinate_proof']['prefix_coordinate_classes']:
        cost=s.sympify(row['diagonal'],locals=variables)
        loss=s.sympify(row['linear_loss'],locals=variables)
        rank=s.sympify(row['rank_domain']['r'],locals=variables)
        target=s.Rational(2,9)*(rank-int(row['modified']))
        point=s.sympify(row['point_distance'],locals=variables)
        assert U.nonnegative(point-target)
        direct=s.factor(cost-loss);distant=s.factor(cost-target/s.Rational(ALPHA))
        if U.nonnegative(direct):kind='loss<=cost';margin=direct
        else:
            assert U.nonnegative(distant) and distant.subs({U.t:0,U.u:0,U.v:0})>0
            kind='cost exceeds new threshold';margin=distant
        bcd.append(dict(type=row['type'],modified=row['modified'],omissions=row['omissions'],
                        contact=row['contact'],reason=kind,margin=str(margin)))
    assert len(bcd)==24 and ALPHA*Q(1,2)>Q(2,9)
    assert (Q(2*166,9)-ALPHA*4-Q(2,9))/2-2==14
    assert B*Q(166-43,9)==Q(287,20)>14
    # The same near-balanced scalar argument cannot exclude n166 at M=1.
    assert (Q(2*165,9)-ALPHA*4-Q(2,9))/2-2<14
    files=[f'scripts/{TAG}.py']
    for tag in (PRIMARY,BASE+'_review',BASE+'_dimension169'):
        files.extend([f'scripts/{tag}.py',f'results/{tag}.json',f'results/{tag}.md'])
    digest=lambda p:sha256((ROOT/p).read_bytes()).hexdigest()
    for p in files:
        if p.endswith('.json'):
            for f,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():assert digest(f)==h,(p,f)
    return dict(status='PASS_INDEPENDENT_GROUNDED_RATIO_DEPENDENT_CHARGE_AND_166_REVIEW',
                sources={p:digest(p) for p in files},exceptional_cases=records,small_cycles=cycles,
                unbounded_cycle_classes=tail,BCD_new_thresholds=bcd,
                charge='7/6',height_ratio_limit='21/20',exceptional_single_contacts=54,
                exceptional_multiple_supports=1,first_excluded_dimension=167,
                ideal_branch_bound=166,compact_facet_branch_bound=143,combined_candidate_bound=166,
                scope='Independent grounded Fraction energy, all exceptional support subsets, path-inverse cycle checks and analytic review. Frozen coordinate identities, geometric common-height/ridge proof and compact143 branch explicitly inherited. Internal research proof candidate, not external validation or novelty certification.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
