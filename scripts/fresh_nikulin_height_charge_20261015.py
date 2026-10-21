#!/usr/bin/env python3
"""Near-balanced unit charge with at most two exceptional defects: candidate165.

No producer arithmetic imports. Literal exceptional Grams, complete low-cost
supports, exact convex-hull dual witnesses, and analytic classical tails.
"""
import argparse
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_height_charge_20261015'
OLD='fresh_nikulin_unit_charge_20261015'
B=s.Rational(31,30)


def star(arms):
    G=2*s.eye(1+sum(arms));k=1
    for arm in arms:
        prev=0
        for _ in range(arm):G[k,prev]=G[prev,k]=-1;prev,k=k,k+1
    return G


def subsets(nodes,allow_empty=False):
    for k in range(0 if allow_empty else 1,len(nodes)+1):yield from combinations(nodes,k)


def height_bound(h,k):
    return s.simplify(max(min(h,k)**2,max(h,k)**2/B**2))


def quadratic_min(E,u):
    P=s.Poly(s.expand(E),u);a,b,c=[P.nth(k)for k in(2,1,0)]
    assert a>=0;trials=[(s.S(0),c),(s.S(1),P.eval(1))]
    if a and 0<=-b/(2*a)<=1:trials.append((-b/(2*a),c-b*b/(4*a)))
    return min(v for _,v in trials),trials


def exceptional(name,G,omissions,target,forbidden):
    delta=G.nullspace()[0]
    if delta[0]<0:delta=-delta
    delta=s.simplify(delta/min(delta));assert all(x>=1 for x in delta)
    P=(G+delta*delta.T).inv()
    def R(i,j):
        v=s.zeros(G.rows,1);v[i]+=1/delta[i];v[j]-=1/delta[j]
        return s.simplify((v.T*P*v)[0])
    inverses=[];diagonals=[]
    for i in omissions:
        kept=[k for k in range(G.rows)if k!=i];Q=G.extract(kept,kept);J=Q.inv()
        assert all(Q[:k,:k].det()>0 for k in range(1,Q.rows+1))and all(c>=0 for c in J)
        inverses.append((kept,J));diagonals.append({k:J[a,a]for a,k in enumerate(kept)})
    point=s.simplify(min(delta[omissions[0]],delta[omissions[1]])**2*R(*omissions))
    assert point>=target
    # If C>=target no defect is possible. Otherwise every contact occurs
    # in this exact cheap list. Enumerating every subset has no support cutoff.
    cheap=[(side,k,c)for side,d in enumerate(diagonals)for k,c in d.items()if c<target]
    supports=[]
    for support in subsets(cheap):
        cost=sum((c for _,_,c in support),s.S(0))
        for (side,k,_),(side2,l,_)in combinations(support,2):
            if side==side2:
                kept,J=inverses[side];cost+=2*J[kept.index(k),kept.index(l)]
        if cost<target:supports.append((support,s.simplify(cost)))
    singles=[];multiple=[];u,v=s.symbols('u v',real=True)
    for support,cost in supports:
        if len(support)==1:
            side,k,_=support[0];i,j=omissions[side],omissions[1-side]
            E=s.expand((1-u)*R(i,j)+u*R(k,j)-u*(1-u)*R(i,k));emin,trials=quadratic_min(E,u)
            H0=height_bound(delta[k],delta[j]);H1=height_bound(delta[i]+delta[k],delta[j])
            L0=s.simplify(H0*R(k,j)+cost);L1=s.simplify(H1*emin+cost)
            assert L0>=target and L1>=target,(name,forbidden,side,k,L0,L1,target)
            singles.append(dict(side=side,node=k,cost=str(cost),energy=str(E),energy_min=str(emin),
                                minimum_trials=[[str(a),str(b)]for a,b in trials],
                                zero_omitted_lower=str(L0),positive_omitted_lower=str(L1),
                                zero_omitted_margin=str(L0-target),positive_omitted_margin=str(L1-target)))
        else:
            assert len(support)==2
            nodes=[(side,k)for side,k,_ in support]
            special=name=='E7'and forbidden is None and nodes==[(0,4),(1,7)]
            if special:
                z=s.zeros(G.rows,1);z[3]=(1-u)/delta[3];z[4]=u/delta[4]
                z[6]=-(1-v)/delta[6];z[7]=-v/delta[7]
                E=s.expand((z.T*P*z)[0]);assert E==(1+u*u+v*v)/2 and cost==1
                # If neither omitted entry is present, u=v=1. Otherwise at
                # least one actual height is >=3, and both are >=3/B.
                no_omitted=s.Rational(3,2)+cost;some_omitted=s.Rational(9,2)/B**2+cost
                assert min(no_omitted,some_omitted)>target
                multiple.append(dict(nodes=list(map(list,nodes)),cost=str(cost),status='UNIT_CHARGE_VIA_E7_HEIGHT',energy=str(E),
                                     neither_omitted_lower=str(no_omitted),some_omitted_lower=str(some_omitted)))
            else:
                defect=s.simplify(target-cost)
                assert cost>=s.Rational(3,2)and defect<=s.Rational(2,9)
                assert name in('E7','E8')
                multiple.append(dict(nodes=list(map(list,nodes)),cost=str(cost),status='POSSIBLE_BUDGETED_DEFECT',defect_upper=str(defect)))
    return dict(type=name,rank=G.rows-1,forbidden=forbidden,omissions=list(omissions),target=str(target),
                gram=[[str(x)for x in r]for r in G.tolist()],kernel=list(map(str,delta)),point_energy=str(point),
                inverse_diagonals=[{str(k):str(c)for k,c in d.items()}for d in diagonals],
                all_cheap_contact_count=len(cheap),all_cheap_subsets_checked=2**len(cheap)-1,
                low_cost_supports=len(supports),single_cases=singles,multiple_cases=multiple)


def convex_dual(m,U,V):
    """Exact closest-pair witness; the returned full KKT signs prove the bound."""
    R=lambda k:s.Rational((k%m)*(m-k%m),m)
    K=lambda i,j:(R(i)+R(j)-R(i-j))/2
    assert all(K(i,j)==min(i,j)-s.Rational(i*j,m)for i in set(U+V)for j in set(U+V))
    nodes=tuple(U)+tuple(V);n=len(nodes);na=len(U);sign=[1]*na+[-1]*len(V)
    H=s.Matrix(n,n,lambda i,j:sign[i]*sign[j]*K(nodes[i],nodes[j]));A=s.zeros(2,n)
    for i in range(n):A[0 if i<na else 1,i]=1
    candidates=[]
    for left in subsets(tuple(range(na))):
        for right in subsets(tuple(range(na,n))):
            active=left+right;Q=H.extract(active,active);C=A[:,active]
            M=Q.row_join(C.T).col_join(C.row_join(s.zeros(2,2)))
            rhs=s.Matrix([0]*len(active)+[1,1]);sol=s.linsolve((M,rhs))
            if sol is s.EmptySet:continue
            row=next(iter(sol))
            if any(x.free_symbols for x in row):continue
            x=s.zeros(n,1)
            for i,val in zip(active,row[:-2]):x[i]=val
            if any(val<0 for val in x):continue
            lag=s.Matrix(row[-2:]);res=H*x+A.T*lag
            if any(val<0 for val in res):continue
            assert A*x==s.ones(2,1)and all(x[i]*res[i]==0 for i in range(n))
            energy=s.simplify((x.T*H*x)[0]);assert energy>=0
            candidates.append((energy,tuple(x),tuple(lag),tuple(res)))
    assert candidates,(m,U,V)
    # One feasible dual certificate is sufficient; no claim depends on
    # enumerating all degenerate minimizers or solving a floating QP.
    best=min(candidates,key=lambda z:(z[0],z[1]));energy,x,lag,res=best
    return dict(minimum_energy=str(energy),barycentric_weights=list(map(str,x)),
                equality_multipliers=list(map(str,lag)),dual_slacks=list(map(str,res)))


def classical_A():
    rows=[];count=0
    for m in range(2,18):
        h=m//2;T=s.Rational(2*(m-1),9);R=lambda k:s.Rational((k%m)*(m-k%m),m)
        domains=[]
        for omitted in(0,h):
            nodes=[k for k in range(m)if k!=omitted and R(k-omitted)<=2];domain=[]
            # All subsets of the complete per-contact budget domain.
            for support in subsets(nodes,True):
                cost=sum((R(k-omitted)for k in support),s.S(0))
                cost+=sum((R(k-omitted)+R(l-omitted)-R(k-l)for k,l in combinations(support,2)),s.S(0))
                if cost<=2:domain.append((support,cost))
            domains.append(domain)
        cases=[]
        for (left,cp),(right,cq)in product(*domains):
            if cp+cq>=T:continue
            for a,b in product((False,True),repeat=2):
                U=left+((0,)if a else());V=right+((h,)if b else())
                if not U or not V:continue
                witness=convex_dual(m,U,V);H=height_bound(len(U),len(V))
                lower=s.simplify(H*s.Rational(witness['minimum_energy'])+cp+cq)
                assert lower>=T,(m,U,V,lower,T)
                cases.append(dict(left_support=list(U),right_support=list(V),left_min_cost=str(cp),right_min_cost=str(cq),
                                  minimum_height_square=str(H),convex_dual=witness,unit_charge_lower=str(lower),margin=str(lower-T)))
        rows.append(dict(cycle_size=m,rank=m-1,omissions=[0,h],target=str(T),
                         retained_domains=[[[list(z),str(c)]for z,c in d]for d in domains],cases=cases))
        count+=len(cases)
    assert count==580
    t=s.symbols('t',nonnegative=True,integer=True);tails=[]
    for m,eps in[(18+2*t,0),(19+2*t,1)]:
        # Each endpoint cost<=2 excludes distance>=3. For d<=2,
        # ell/c <= (m+eps)/(m-2). The entire two-row cost is <=4.
        R=m/4 if not eps else(m*m-1)/(4*m);T=s.Rational(2,9)*(m-1)
        margin=s.factor(R-T-4*(2+eps)/(m-2));num,den=s.fraction(margin)
        assert all(c>=0 for c in s.Poly(num,t).all_coeffs())and num.subs(t,0)>0
        assert all(c>=0 for c in s.Poly(den,t).all_coeffs())and den.subs(t,0)>0
        assert s.expand(3*(m-3)-2*m)==m-9
        tails.append(dict(cycle_size=str(m),parity=eps,unit_charge_margin=str(margin),
                          cleared_numerator=str(s.expand(num)),positive_denominator=str(s.expand(den))))
    return dict(small_cycles=rows,exact_convex_dual_cases=count,all_rank_tail_classes=tails,
                avoid_node='Rotation of opposite nodes proves the stronger ordinary target whenever possible; affineA1 coincident omissions have target0.')


def classical_BCD(old):
    r,k=s.symbols('r k',positive=True,integer=True);t,u,v=s.symbols('t u v',nonnegative=True,integer=True)
    variables=dict(r=r,k=k,t=t,u=u,v=v);rows=[]
    for row in old['classical_all_rank_checks']['prefix_coordinate_cases']:
        cost=s.sympify(row['inverse_diagonal'],locals=variables);ell=s.sympify(row['linear_loss'],locals=variables)
        T=s.Rational(2,9)*(r-1 if row['modified']else r)
        far=row['proof']=='cost exceeds trivial-charge threshold';margin=cost-T if far else cost-ell
        subs={variables[a]:s.sympify(b,locals=variables)for a,b in row['domain_substitution'].items()}
        val=s.factor(margin.subs(subs,simultaneous=True));num,den=s.fraction(val)
        assert all(c>=0 for c in s.Poly(num,t,u,v).coeffs())and all(c>=0 for c in s.Poly(den,t,u,v).coeffs())
        assert den.subs({t:0,u:0,v:0})>0
        if far:assert val.subs({t:0,u:0,v:0})>0
        rows.append(dict(type=row['type'],modified=row['modified'],side=row['side'],node_class=row['node_class'],
                         inverse_diagonal=str(cost),linear_loss=str(ell),unit_margin=str(val),proof=row['proof']))
    assert len(rows)==24
    assert s.Rational(1,2)>s.Rational(2,9)
    return dict(all_rank_classes=rows,modified_B2='Any nonzero retained cost>=1/2>target2/9; zero-cost point energy1.',
                inherited_input='Complete prefix-coordinate identities and index-class coverage from frozen independently reviewed4/3 proof.')


def run():
    old=json.loads((ROOT/f'results/{OLD}.json').read_text());previous=json.loads((ROOT/'results/fresh_nikulin_coupled_charge_20261015.json').read_text())
    mats={'A1':s.Matrix([[2,-2],[-2,2]]),'E6':star([2,2,2]),'E7':star([1,3,3]),'E8':star([1,2,5]),
          'F4':s.Matrix([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-s.sqrt(2),0],[0,0,-s.sqrt(2),2,-1],[0,0,0,-1,2]]),
          'G2':s.Matrix([[2,-1,0],[-1,2,-s.sqrt(3)],[0,-s.sqrt(3),2]])}
    records=[]
    for frozen in previous['exceptional_cases']:
        name=frozen['type'];rec=exceptional(name,mats[name],tuple(frozen['omissions']),s.Rational(frozen['target']),frozen['forbidden'])
        for key in('type','rank','forbidden','omissions','target','gram','kernel'):assert rec[key]==frozen[key]
        records.append(rec)
    assert len(records)==40
    singles=sum(len(r['single_cases'])for r in records);multiple=sum(len(r['multiple_cases'])for r in records)
    possible=[c for r in records for c in r['multiple_cases']if c['status']=='POSSIBLE_BUDGETED_DEFECT']
    assert singles==71 and multiple==16 and len(possible)==15
    A=classical_A();BCD=classical_BCD(old)
    t=s.symbols('t',nonnegative=True,integer=True)
    # Three bad factors already cost>=9/2, beyond the global budget4.
    assert 3*s.Rational(3,2)>4
    loss=4+2*s.Rational(2,9);assert loss==s.Rational(40,9)
    near=s.expand((s.Rational(2,9)*(165+t)-s.Rational(2,9)-loss)/2-2)
    far=s.expand(B*(165+t-43)/9)
    assert near==14+t/9 and far==14+s.Rational(1,135)+s.Rational(31,270)*t
    paths=[f'scripts/{TAG}.py']
    for stem in(OLD,OLD+'_review',OLD+'_dimension169','fresh_nikulin_coupled_charge_20261015',
                'fresh_nikulin_coupled_charge_20261015_review','fresh_nikulin_coupled_charge_20261015_dimension166'):
        paths += [f'scripts/{stem}.py',f'results/{stem}.json']
        data=json.loads((ROOT/f'results/{stem}.json').read_text())
        for name,h in data.get('sources',{}).items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    return dict(status='PASS_EXACT_HEIGHT_CHARGE_AND_TWO_DEFECT_AGGREGATION_CANDIDATE165_PENDING_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},height_window=str(B),
                exceptional_cases=records,classical_A=A,classical_BCD=BCD,
                counts=dict(exceptional_records=40,single_cases=singles,multiple_cases=multiple,
                            possible_defect_occurrences=15,classical_convex_duals=580,classical_symbolic_classes=26),
                aggregate=dict(unit_charge_exception='Only retained support patterns of cost>=3/2; defect<=2/9.',
                               total_cost_budget=4,maximum_defective_factors=2,total_loss=str(loss),
                               near_pairing_lower='M*(N-39)/9',far_pairing_lower='M*(N-43)/9',
                               near_N165_plus_t=str(near),far_N165_plus_t=str(far),
                               ideal_branch_at_most=165,compact_facet_branch_at_most=143,combined_candidate=165),
                scope='Research proof candidate pending independent review. Same preset omissions; arbitrary real nonzero contacts>=1; all affine ranks and every avoid-node choice. No novelty/current-best claim; inherited geometric169/166 chain remains unchanged.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'])
