#!/usr/bin/env python3
"""Exact obstruction to unit charge and a uniform optimal 4/3 local charge.

Exceptional support budgets are complete; classical ranks use closed formulas.
No finite-label cutoff and no producer arithmetic import.
"""
import argparse
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_unit_charge_20261015'
ALPHA=s.Rational(4,3)


def star(arms):
    G=2*s.eye(1+sum(arms));k=1
    for length in arms:
        last=0
        for _ in range(length):G[k,last]=G[last,k]=-1;last=k;k+=1
    return G


def data(G,omissions):
    null=G.nullspace();assert len(null)==1
    delta=null[0]
    if delta[0]<0:delta=-delta
    delta=s.simplify(delta/min(delta));assert all(x>=1 for x in delta)
    P=(G+delta*delta.T).inv()
    def resistance(i,j):
        v=s.zeros(G.rows,1);v[i]+=1/delta[i];v[j]-=1/delta[j]
        return s.simplify((v.T*P*v)[0])
    inverses=[];diags=[]
    for i in omissions:
        kept=[k for k in range(G.rows)if k!=i];B=G.extract(kept,kept)
        J=B.inv();assert all(x>=0 for x in J)
        assert all(B[:k,:k].det()>0 for k in range(1,B.rows+1))
        inverses.append((kept,J));diags.append({k:J[a,a]for a,k in enumerate(kept)})
    i,j=omissions
    L=s.simplify(min(delta[i]**2,delta[j]**2)*resistance(i,j))
    return delta,P,resistance,inverses,diags,L


def qminimum(q,u):
    q=s.Poly(q,u);c,b,a=(q.nth(k)for k in range(3))
    assert a>=0
    vals=[q.eval(0),q.eval(1)]
    if a and 0<=-b/(2*a)<=1:vals.append(c-b*b/(4*a))
    return s.simplify(min(vals))


def exceptional_case(name,G,omissions,target,forbidden=None):
    delta,P,R,inverses,diags,L=data(G,omissions);assert L>=target
    limit=target/ALPHA
    cheap=[(side,k,c)for side,dd in enumerate(diags)for k,c in dd.items()if c<limit]
    records=[];u=s.symbols('u',real=True)
    for side,k,cost in cheap:
        i,j=omissions[side],omissions[1-side]
        triple=[R(i,j),R(k,j),R(i,k)]
        E=s.expand((1-u)*triple[0]+u*triple[1]-u*(1-u)*triple[2])
        minimum=qminimum(E,u)
        zero=s.simplify(min(delta[k]**2,delta[j]**2)*triple[1]+ALPHA*cost)
        positive=s.simplify(min(delta[i]**2,delta[j]**2)*minimum+ALPHA*cost)
        assert zero>=target and positive>=target,(name,omissions,k,zero,positive,target)
        records.append(dict(side=side,node=k,cost=str(cost),null_coefficients=list(map(str,[delta[i],delta[j],delta[k]])),
                            resistances=list(map(str,triple)),energy_polynomial=str(E),minimum_energy=str(minimum),
                            zero_omitted_total_lower=str(zero),positive_omitted_total_lower=str(positive)))
    # A diagonal lower bound excludes all two-contact supports except the two
    # isolated cheap leaves of the non-narrow E7 pair.
    low_pairs=[]
    for x,y in combinations(cheap,2):
        side,k,cost=x;side2,k2,cost2=y;value=cost+cost2
        if side==side2:
            kept,J=inverses[side];value+=2*J[kept.index(k),kept.index(k2)]
        if value<limit:low_pairs.append((x,y,value))
    two_case=None
    if low_pairs:
        assert name=='E7'and forbidden is None and tuple(omissions)==(3,6)
        assert low_pairs==[((0,4,s.Rational(1,2)),(1,7,s.Rational(1,2)),s.S(1))]
        assert 3*min(c for _,_,c in cheap)>=limit
        v=s.symbols('v',real=True);row=s.zeros(G.rows,1)
        row[3]=(1-u)/delta[3];row[4]=u/delta[4]
        row[6]=-(1-v)/delta[6];row[7]=-v/delta[7]
        E=s.expand((row.T*P*row)[0]);assert E==(1+u*u+v*v)/2
        assert s.Rational(1,2)+ALPHA>=target
        two_case=dict(nodes=[[0,4],[1,7]],energy_polynomial=str(E),cost_lower='1',
                      combined_lower=str(s.Rational(1,2)+ALPHA))
    return dict(type=name,rank=G.rows-1,forbidden=forbidden,omissions=list(omissions),target=str(target),
                gram=[[str(x)for x in row]for row in G.tolist()],kernel=list(map(str,delta)),
                point_energy=str(L),inverse_diagonals=[{str(k):str(v)for k,v in dd.items()}for dd in diags],
                trivial_charge_threshold=str(limit),single_contact_cases=records,two_contact_case=two_case)


def omission_countercensus(name,G):
    r=G.rows-1;target=s.Rational(2*r,9);eligible=[];excluded=[]
    for i,j in combinations(range(G.rows),2):
        delta,P,R,inverses,diags,L=data(G,(i,j))
        if L<target:excluded.append(dict(omissions=[i,j],zero_cost_point_energy=str(L)));continue
        costs=[{i:s.S(0),**diags[0]},{j:s.S(0),**diags[1]}]
        worst=s.S(0);witness=None
        for a,b in product(range(G.rows),repeat=2):
            c,d=costs[0][a],costs[1][b]
            if c>2 or d>2 or not c+d:continue
            W=s.simplify(min(delta[a]**2,delta[b]**2)*R(a,b))
            bound=s.factor((target-W)/(c+d))
            if bound>worst:
                worst=bound;witness=dict(row_nodes=[a,b],heights=[str(delta[a]),str(delta[b])],
                                        endpoint_costs=[str(c),str(d)],weighted_energy=str(W),
                                        coefficient_lower_bound=str(bound))
        eligible.append(dict(omissions=[i,j],zero_cost_point_energy=str(L),required_charge_lower=str(worst),witness=witness))
    assert len(eligible)+len(excluded)==G.rows*(G.rows-1)//2
    return dict(type=name,rank=r,target=str(target),all_unordered_omission_pairs=G.rows*(G.rows-1)//2,
                ineligible_zero_cost_pairs=excluded,eligible_pairs=eligible,
                best_point_row_lower_bound=str(min(s.Rational(x['required_charge_lower'])for x in eligible)))


def classical_checks():
    # A_m cycle: m=r+1, opposite distance floor(m/2), epsilon=m mod2.
    # Ell <= (m+epsilon)d/m, cost=d(m-d)/m. Positivity after multiplying
    # by3m is d*(m-4d-3epsilon).
    t,u,v=s.symbols('t u v',nonnegative=True,integer=True);cycle=[]
    for m in range(2,7):assert s.Rational(m-1,m)>=s.Rational(m-1,6)
    for m in (7,8,9):
        assert s.Rational(2*(m-2),m)>=s.Rational(m-1,6)
        epsilon=m%2;assert m-4-3*epsilon>=0
        cycle.append(dict(m=str(m),allowed_distance=[1],cleared_charge_margin=str(m-4-3*epsilon)))
    for m,epsilon in ((10+2*t,0),(11+2*t,1)):
        # d=3 is already above the endpoint budget2; diagonals increase
        # from d=0 to the antipode, so every allowed retained node has d<=2.
        assert s.simplify(3*(m-3)-2*m)==m-9
        for d in (1,2):
            margin=s.expand(d*(m-4*d-3*epsilon));assert all(c>=0 for c in s.Poly(margin,t).all_coeffs())
            cycle.append(dict(m=str(m),allowed_distance=[d],cleared_charge_margin=str(margin)))
    m=s.symbols('m',positive=True,integer=True)
    assert s.cancel(m/4-s.Rational(2,9)*(m-1)-(m+8)/36)==0
    assert s.cancel((m*m-1)/(4*m)-s.Rational(2,9)*(m-1)-(m-1)*(m+9)/(36*m))==0
    # Coordinate table: R values are squared distances between the prefix
    # potentials described in the report. The table covers every index class.
    r,k=s.symbols('r k',integer=True,positive=True)
    rows=[]
    def add(typ,mod,side,nodes,delta2,Rij,Rik,Rkj,sub,far=False):
        cost=s.factor(delta2*Rik);ell=s.factor(Rij+Rik-Rkj)
        target=s.Rational(2,9)*(r-1 if mod else r)
        margin=cost-target/ALPHA if far else cost-ell
        # Two independent nonnegative parameters cover every allowed k and
        # every excess rank r-k; do not sample the diagonal k=r/2.
        if k in sub:
            kval=sub[k].subs(t,u)
            sub={k:kval,r:sub[r].subs(t,v).subs(k,kval)}
        expr=s.factor(margin.subs(sub,simultaneous=True));num,den=s.fraction(expr)
        assert all(c>=0 for c in s.Poly(s.expand(num),t,u,v).coeffs())
        assert all(c>=0 for c in s.Poly(s.expand(den),t,u,v).coeffs())and den.subs({t:0,u:0,v:0})>0
        if far:assert expr.subs({t:0,u:0,v:0})>0
        rows.append(dict(type=typ,modified=mod,side=side,node_class=nodes,null_coefficient_squared=str(delta2),
                         point_resistances=list(map(str,[Rij,Rik,Rkj])),inverse_diagonal=str(cost),linear_loss=str(ell),
                         proof='cost exceeds trivial-charge threshold'if far else'linear loss at most diagonal',
                         domain_substitution={str(a):str(b)for a,b in sub.items()},cleared_margin=str(expr)))
    # B: v0=0,v1=e1,vk=(e1+...+ek)/2 for2<=k<=r.
    # Generic k-cases use r=k+1+t (or k+2+t); residual margins are then
    # independent of rank and positive for every k>=its indicated minimum.
    add('B',False,0,'1',1,r/4,1,r/4,{r:2+t})
    add('B',False,0,'2<=k<=r-1',4,r/4,k/4,(r-k)/4,{r:k+1+t,k:2+t})
    add('B',False,0,'r',2,r/4,r/4,0,{r:2+t})
    add('B',False,'r','0',1,r/4,r/4,0,{r:2+t},True)
    add('B',False,'r','1',1,r/4,r/4,1,{r:2+t},True)
    add('B',False,'r','2<=k<=r-1',4,r/4,(r-k)/4,k/4,{r:k+1+t,k:2+t})
    # Modified B: omit r-1 and0, r>=3. r2 is checked separately.
    add('B',True,'r-1','0',1,(r-1)/4,(r-1)/4,0,{r:3+t},True)
    add('B',True,'r-1','1',1,(r-1)/4,(r-1)/4,1,{r:3+t},True)
    add('B',True,'r-1','2<=k<=r-2',4,(r-1)/4,(r-1-k)/4,k/4,{r:k+2+t,k:2+t})
    add('B',True,'r-1','r',2,(r-1)/4,s.Rational(1,4),r/4,{r:3+t})
    add('B',True,0,'1',1,(r-1)/4,1,(r-1)/4,{r:3+t})
    add('B',True,0,'2<=k<=r-1',4,(r-1)/4,k/4,(r-1-k)/4,{r:k+1+t,k:2+t})
    add('B',True,0,'r',2,(r-1)/4,r/4,s.Rational(1,4),{r:3+t})
    # C: vk=(e1+...+ek)/sqrt2. The other ordinary side is reversal.
    add('C',False,0,'1<=k<=r-1',2,r/2,k/2,(r-k)/2,{r:k+1+t,k:1+t})
    add('C',False,0,'r',1,r/2,r/2,0,{r:2+t},True)
    add('C',True,1,'0',1,(r-1)/2,s.Rational(1,2),r/2,{r:2+t})
    add('C',True,1,'2<=k<=r-1',2,(r-1)/2,(k-1)/2,(r-k)/2,{r:k+1+t,k:2+t})
    add('C',True,1,'r',1,(r-1)/2,(r-1)/2,0,{r:2+t},True)
    add('C',True,'r','0',1,(r-1)/2,r/2,s.Rational(1,2),{r:2+t},True)
    add('C',True,'r','1<=k<=r-1',2,(r-1)/2,(r-k)/2,(k-1)/2,{r:k+1+t,k:1+t})
    # D: v_{r-1}=(e1+...+e_{r-1}-er)/2; vr=sum(e)/2.
    # Fork reversal gives the other omitted-side formulas.
    add('D',False,0,'1',1,r/4,1,r/4,{r:4+t})
    add('D',False,0,'2<=k<=r-2',4,r/4,k/4,(r-k)/4,{r:k+2+t,k:2+t})
    add('D',False,0,'r-1',1,r/4,r/4,1,{r:4+t},True)
    add('D',False,0,'r',1,r/4,r/4,0,{r:4+t},True)
    assert s.Rational(4,3)*s.Rational(1,2)>=s.Rational(2,9) # modified B2
    return dict(cycle_cases=cycle,prefix_coordinate_cases=rows,modified_B2=dict(point_energy='1',minimum_cost='1/2',target='2/9'))


def run():
    mats={'A1':s.Matrix([[2,-2],[-2,2]]),'E6':star([2,2,2]),'E7':star([1,3,3]),'E8':star([1,2,5]),
          'F4':s.Matrix([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-s.sqrt(2),0],[0,0,-s.sqrt(2),2,-1],[0,0,0,-1,2]]),
          'G2':s.Matrix([[2,-1,0],[-1,2,-s.sqrt(3)],[0,-s.sqrt(3),2]])}
    omissions={'A1':(0,1),'E6':(2,4),'E7':(3,6),'E8':(1,6),'F4':(0,4),'G2':(0,2)}
    records=[]
    for name,G in mats.items():
        r=G.rows-1
        records.append(exceptional_case(name,G,omissions[name],s.Rational(2*r,9)))
        for f in range(G.rows):
            if name=='A1':ij=(1-f,1-f)
            elif name=='E6':ij=tuple(x for x in (2,4,6)if x!=f)[:2]
            elif name=='E7':ij=(2,5)if f in (4,7)else(4,7)
            elif name=='E8':ij=(1,6)if f in (3,7)else(3,7)
            elif name=='F4':ij=(0,4)if f in (1,3)else(1,3)
            else:ij=tuple(x for x in range(3)if x!=f)
            assert f not in ij
            records.append(exceptional_case(name,G,ij,s.Rational(2*(r-1),9),f))
    assert len(records)==40
    census=[omission_countercensus(name,mats[name])for name in ('E6','E7','E8')]
    assert [len(x['eligible_pairs'])for x in census]==[6,1,4]
    assert [x['best_point_row_lower_bound']for x in census]==['1','10/9','4/3']
    classical=classical_checks()
    t=s.symbols('t',nonnegative=True,integer=True)
    pairing=s.expand((s.Rational(2,9)*(169+t)-ALPHA*4-s.Rational(2,9))/2-2)
    assert pairing==14+t/9
    paths=[f'scripts/{TAG}.py']
    for dep in ['fresh_nikulin_avoid_node_20261015','fresh_nikulin_avoid_node_20261015_review',
                'fresh_nikulin_dimension_181_candidate_20261015','fresh_nikulin_narrow_20261015',
                'fresh_nikulin_narrow_20261015_products','fresh_nikulin_energy_charge_20261015_review']:
        paths.extend([f'scripts/{dep}.py',f'results/{dep}.json',f'results/{dep}.md'])
        doc=json.loads((ROOT/'results'/f'{dep}.json').read_text())
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    return dict(status='PASS_EXACT_UNIT_CHARGE_OBSTRUCTION_AND_FOUR_THIRDS_UNIFORM_LEMMA',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                charge_coefficient='4/3',unit_charge_impossible=True,exceptional_cases=records,
                complete_point_row_omission_censuses=census,classical_all_rank_checks=classical,
                local_nonnarrow_statement='energy/M +(4/3)C >=2r/9',
                local_avoid_node_statement='energy/M +(4/3)C >=2(r-1)/9',
                endpoint_cost_budget=4,avoid_penalty='2/9',negative_pairing_lower_bound='M*(N-43)/9',
                threshold_N169_plus_t=str(pairing),ideal_branch_dimension_at_most=169,
                compact_facet_branch_dimension_at_most=143,combined_candidate_dimension_at_most=169,
                scope='Exact exceptional cases and complete point-row omission obstructions, plus symbolic classical rank formulas. Uniform all-real contacts with every nonzero entry>=1 and each endpoint Schur budget<=2; no label cutoff. The 4/3 coefficient is optimal for this fixed-slope preset-omission local lemma. Global consequence inherits frozen framing/endpoint and compact-facet proofs; independent review is still required. No novelty/current-best claim.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(out.read_text())==result
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['combined_candidate_dimension_at_most'])
