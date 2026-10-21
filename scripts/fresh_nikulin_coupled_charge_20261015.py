#!/usr/bin/env python3
"""Exact 7/6 charge for nearly balanced affine heights; candidate bound166.

All exceptional supports are covered. Classical A uses finite small-rank
resistance formulas and symbolic tails; B/C/D coordinate-class identities are
explicitly inherited from, and rebound to, the frozen reviewed4/3 proof.
"""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_coupled_charge_20261015'
OLD='fresh_nikulin_unit_charge_20261015'
ALPHA=s.Rational(7,6)
BALANCE=s.Rational(21,20)


def star(arms):
    G=2*s.eye(1+sum(arms));k=1
    for arm in arms:
        previous=0
        for _ in range(arm):
            G[k,previous]=G[previous,k]=-1
            previous,k=k,k+1
    return G


def quadratic_minimum(E,u):
    p=s.Poly(E,u);a,b,c=p.nth(2),p.nth(1),p.nth(0)
    assert a>=0
    trials=[(s.S(0),c),(s.S(1),p.eval(1))]
    if a and 0<=-b/(2*a)<=1:
        trials.append((-b/(2*a),c-b*b/(4*a)))
    return s.simplify(min(v for _,v in trials)),trials


def height_square(h1,h2):
    # Each height>=its h_i. Their ratio is at most BALANCE, hence the smaller
    # height is also at least max(h1,h2)/BALANCE.
    return s.simplify(max(min(h1,h2)**2,max(h1,h2)**2/BALANCE**2))


def exceptional(name,G,omissions,target,forbidden):
    delta=G.nullspace()[0]
    if delta[0]<0:delta=-delta
    delta=s.simplify(delta/min(delta));assert all(x>=1 for x in delta)
    P=(G+delta*delta.T).inv()
    def resistance(i,j):
        z=s.zeros(G.rows,1);z[i]+=1/delta[i];z[j]-=1/delta[j]
        return s.simplify((z.T*P*z)[0])
    inverses=[];diagonals=[]
    for i in omissions:
        kept=[k for k in range(G.rows)if k!=i]
        B=G.extract(kept,kept);J=B.inv()
        assert all(B[:k,:k].det()>0 for k in range(1,B.rows+1))
        assert all(v>=0 for v in J)
        inverses.append((kept,J));diagonals.append({k:J[a,a]for a,k in enumerate(kept)})
    i,j=omissions
    point=s.simplify(min(delta[i],delta[j])**2*resistance(i,j))
    assert point>=target
    limit=target/ALPHA
    cheap=[(side,k,c)for side,dd in enumerate(diagonals)for k,c in dd.items()if c<limit]
    singles=[];u=s.symbols('u',real=True)
    for side,k,cost in cheap:
        i,j=omissions[side],omissions[1-side]
        triple=[resistance(i,j),resistance(k,j),resistance(i,k)]
        E=s.expand((1-u)*triple[0]+u*triple[1]-u*(1-u)*triple[2])
        emin,trials=quadratic_minimum(E,u)
        hzero=height_square(delta[k],delta[j])
        hpositive=height_square(delta[i]+delta[k],delta[j])
        lower0=s.simplify(hzero*triple[1]+ALPHA*cost)
        lower1=s.simplify(hpositive*emin+ALPHA*cost)
        assert lower0>=target and lower1>=target,(name,omissions,k,lower0,lower1,target)
        singles.append(dict(side=side,node=k,cost=str(cost),
                            kernel_coefficients=list(map(str,[delta[i],delta[j],delta[k]])),
                            resistances=list(map(str,triple)),energy_polynomial=str(E),
                            minimum_trials=[[str(x),str(y)]for x,y in trials],minimum_energy=str(emin),
                            zero_omitted_height_square=str(hzero),positive_omitted_height_square=str(hpositive),
                            zero_omitted_lower=str(lower0),positive_omitted_lower=str(lower1),
                            zero_omitted_margin=str(lower0-target),positive_omitted_margin=str(lower1-target)))
    surviving_pairs=[]
    for x,y in combinations(cheap,2):
        side,k,cost=x;side2,k2,cost2=y
        total=cost+cost2
        if side==side2:
            kept,J=inverses[side]
            total+=2*J[kept.index(k),kept.index(k2)]
        if total<limit:surviving_pairs.append((x,y,total))
    if cheap:
        assert 3*min(c for _,_,c in cheap)>=limit
    two=None
    if surviving_pairs:
        assert name=='E7'and forbidden is None and tuple(omissions)==(3,6)
        assert surviving_pairs==[((0,4,s.Rational(1,2)),(1,7,s.Rational(1,2)),s.S(1))]
        v=s.symbols('v',real=True);z=s.zeros(G.rows,1)
        z[3]=(1-u)/delta[3];z[4]=u/delta[4]
        z[6]=-(1-v)/delta[6];z[7]=-v/delta[7]
        E=s.expand((z.T*P*z)[0]);assert E==(1+u*u+v*v)/2
        lower=s.Rational(1,2)+ALPHA
        assert lower==s.Rational(5,3)>target
        two=dict(nodes=[[0,4],[1,7]],energy=str(E),cost_lower='1',
                 combined_lower=str(lower),margin=str(lower-target))
    return dict(type=name,rank=G.rows-1,forbidden=forbidden,omissions=list(omissions),target=str(target),
                gram=[[str(x)for x in row]for row in G.tolist()],kernel=list(map(str,delta)),
                point_energy=str(point),inverse_diagonals=[{str(k):str(c)for k,c in dd.items()}for dd in diagonals],
                cost_threshold=str(limit),single_contact_cases=singles,two_contact_case=two,
                three_contacts_excluded=bool(not cheap or 3*min(c for _,_,c in cheap)>=limit))


def cycle_checks():
    records=[]
    for m in range(2,20):
        r=m-1;T=s.Rational(2*r,9);limit=T/ALPHA;h=m//2
        def R(k):
            k=k%m
            return s.Rational(k*(m-k),m)
        Rij=R(h);allowed=[]
        for k in range(1,m):
            cost=R(k)
            if cost<=2 and cost<limit:
                ell=Rij+cost-R(k-h)
                allowed.append(dict(node=k,cost=str(cost),linear_loss=str(ell),ratio=str(ell/cost)))
        lam=max([ALPHA]+[s.Rational(row['ratio'])for row in allowed])
        margin=s.factor(Rij-T-(lam-ALPHA)*limit)
        assert margin>=0,(m,margin)
        records.append(dict(cycle_size=m,rank=r,target=str(T),point_resistance=str(Rij),
                            cost_threshold=str(limit),allowed_retained_nodes=allowed,
                            maximum_loss_per_cost=str(lam),uniform_margin=str(margin)))
    t=s.symbols('t',nonnegative=True,integer=True)
    tails=[]
    for m,eps in [(20+2*t,0),(21+2*t,1)]:
        assert s.expand(3*(m-3)-2*m)==m-9
        for d in (1,2):
            # 6m*(ALPHA*c-ell_upper)=d*(m-7d-6eps).
            margin=s.expand(d*(m-7*d-6*eps))
            assert all(c>=0 for c in s.Poly(margin,t).all_coeffs())
            tails.append(dict(cycle_size=str(m),distance=d,cleared_loss_margin=str(margin),
                              distance_three_budget_excess_numerator=str(m-9)))
    m=s.symbols('m',positive=True,integer=True)
    even=s.factor(m/4-s.Rational(2,9)*(m-1))
    odd=s.factor((m*m-1)/(4*m)-s.Rational(2,9)*(m-1))
    assert s.cancel(even-(m+8)/36)==0 and s.cancel(odd-(m-1)*(m+9)/(36*m))==0
    return dict(finite_cases=records,all_rank_tails=tails,
                even_point_slack=str(even),odd_point_slack=str(odd),
                mechanism='If C<T/alpha, E+alpha*C >=Rij-(lambda-alpha)C. '
                          'The listed positive margin covers every small rank. '
                          'For the tails each loss is directly<=alpha times its diagonal cost.')


def classical_bcd(old):
    # These complete coordinate index classes and their resistance identities
    # were independently reviewed in the frozen4/3 package. Recompute the
    # stronger threshold comparison, rather than inferring it from a rank scan.
    r,k=s.symbols('r k',positive=True,integer=True)
    t,u,v=s.symbols('t u v',nonnegative=True,integer=True)
    variables={'r':r,'k':k,'t':t,'u':u,'v':v}
    rows=[]
    for row in old['classical_all_rank_checks']['prefix_coordinate_cases']:
        cost=s.sympify(row['inverse_diagonal'],locals=variables)
        ell=s.sympify(row['linear_loss'],locals=variables)
        T=s.Rational(2,9)*(r-1 if row['modified']else r)
        far=row['proof']=='cost exceeds trivial-charge threshold'
        margin=cost-T/ALPHA if far else cost-ell
        substitution={variables[a]:s.sympify(b,locals=variables)for a,b in row['domain_substitution'].items()}
        value=s.factor(margin.subs(substitution,simultaneous=True))
        numerator,denominator=s.fraction(value)
        assert all(c>=0 for c in s.Poly(s.expand(numerator),t,u,v).coeffs())
        assert all(c>=0 for c in s.Poly(s.expand(denominator),t,u,v).coeffs())
        assert denominator.subs({t:0,u:0,v:0})>0
        if far:assert value.subs({t:0,u:0,v:0})>0
        rows.append(dict(type=row['type'],modified=row['modified'],side=row['side'],
                         node_class=row['node_class'],inverse_diagonal=str(cost),linear_loss=str(ell),
                         proof=row['proof'],domain_substitution=row['domain_substitution'],
                         new_cleared_margin=str(value)))
    assert len(rows)==24
    assert s.cancel(r/4-s.Rational(2,9)*r-r/36)==0
    assert s.cancel((r-1)/4-s.Rational(2,9)*(r-1)-(r-1)/36)==0
    assert s.simplify(r/2-s.Rational(2,9)*r)==5*r/18
    assert s.cancel((r-1)/2-s.Rational(2,9)*(r-1)-5*(r-1)/18)==0
    assert ALPHA*s.Rational(1,2)>s.Rational(2,9)
    return dict(inherited_complete_coordinate_classes=rows,modified_B2=dict(point='1',
                cost_minimum='1/2',charged_minimum=str(ALPHA/2),target='2/9'),
                scope='All ranks; coordinate-class coverage inherited from frozen reviewed4/3 source. '
                      'Every far-node exclusion is rechecked at the new T/(7/6) threshold.')


def run():
    old=json.loads((ROOT/f'results/{OLD}.json').read_text())
    assert old['status']=='PASS_EXACT_UNIT_CHARGE_OBSTRUCTION_AND_FOUR_THIRDS_UNIFORM_LEMMA'
    mats={'A1':s.Matrix([[2,-2],[-2,2]]),'E6':star([2,2,2]),'E7':star([1,3,3]),'E8':star([1,2,5]),
          'F4':s.Matrix([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-s.sqrt(2),0],
                         [0,0,-s.sqrt(2),2,-1],[0,0,0,-1,2]]),
          'G2':s.Matrix([[2,-1,0],[-1,2,-s.sqrt(3)],[0,-s.sqrt(3),2]])}
    ordinary={'A1':(0,1),'E6':(2,4),'E7':(3,6),'E8':(1,6),'F4':(0,4),'G2':(0,2)}
    cases=[]
    for name,G in mats.items():
        rank=G.rows-1
        cases.append(exceptional(name,G,ordinary[name],s.Rational(2*rank,9),None))
        for forbidden in range(G.rows):
            if name=='A1':ij=(1-forbidden,1-forbidden)
            elif name=='E6':ij=tuple(x for x in(2,4,6)if x!=forbidden)[:2]
            elif name=='E7':ij=(2,5)if forbidden in(4,7)else(4,7)
            elif name=='E8':ij=(1,6)if forbidden in(3,7)else(3,7)
            elif name=='F4':ij=(0,4)if forbidden in(1,3)else(1,3)
            else:ij=tuple(x for x in range(3)if x!=forbidden)
            assert forbidden not in ij
            cases.append(exceptional(name,G,ij,s.Rational(2*(rank-1),9),forbidden))
    assert len(cases)==40
    for record,previous in zip(cases,old['exceptional_cases']):
        for key in ['type','rank','forbidden','omissions','target','gram','kernel']:
            assert record[key]==previous[key]
    classical_A=cycle_checks()
    classical_other=classical_bcd(old)
    t=s.symbols('t',nonnegative=True,integer=True)
    small=s.expand((s.Rational(2,9)*(166+t)-ALPHA*4-s.Rational(2,9))/2-2)
    large=s.expand(BALANCE*(166+t-43)/9)
    assert small==14+t/9
    assert large==s.Rational(287,20)+s.Rational(7,60)*t
    paths=[f'scripts/{TAG}.py']
    for stem in [OLD,OLD+'_review',OLD+'_dimension169','fresh_nikulin_avoid_node_20261015',
                 'fresh_nikulin_avoid_node_20261015_review','fresh_nikulin_narrow_20261015_products',
                 'fresh_nikulin_narrow_20261015_compact_facet']:
        paths.extend([f'scripts/{stem}.py',f'results/{stem}.json'])
        data=json.loads((ROOT/f'results/{stem}.json').read_text())
        for p,h in data.get('sources',{}).items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    paths.extend([f'results/{OLD}.md',f'results/{OLD}_review.md',f'results/{OLD}_dimension169.md'])
    return dict(status='PASS_EXACT_BALANCED_SEVEN_SIXTHS_CHARGE_CANDIDATE166',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                charge=str(ALPHA),common_height_ratio_limit=str(BALANCE),
                exceptional_cases=cases,classical_A=classical_A,classical_BCD=classical_other,
                counts=dict(exceptional_cases=40,single_contact_occurrences=sum(len(x['single_contact_cases'])for x in cases),
                            two_contact_occurrences=sum(x['two_contact_case']is not None for x in cases),
                            cycle_small_cases=18,cycle_tail_classes=4,BCD_index_classes=24),
                local_statement='If M<=21/20, W+(7/6)C>=2r/9, or2(r-1)/9 for the prescribed avoid-node factor.',
                global_statement=dict(total_endpoint_cost_budget=4,avoid_penalty='2/9',
                                      near_balanced_pairing_lower='M*(N-40)/9',
                                      other_pairing_lower='M*(N-43)/9',
                                      N166_plus_t_small=str(small),N166_plus_t_large=str(large),
                                      ideal_branch_at_most=166,compact_facet_branch_at_most=143,
                                      combined_candidate_at_most=166),
                scope='Research proof candidate pending independent review. No label cutoff. '
                      'Exact exceptional support coverage and all-rank classical proof; '
                      'same omissions in both height branches. Inherits the reviewed169 geometric '
                      'framing/endpoint, common-height and compact-facet arguments. '
                      'No claim of optimality, formal verification, or publication priority.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
