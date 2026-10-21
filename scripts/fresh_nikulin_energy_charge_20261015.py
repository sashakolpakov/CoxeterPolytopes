#!/usr/bin/env python3
"""Exact contact-cost improvement of the frozen product-cusp estimate.

No label scan: finite diagonal/energy identities plus symbolic rank inequalities.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_energy_charge_20261015'


def star(arms):
    n=1+sum(arms);G=2*s.eye(n);node=1
    for length in arms:
        previous=0
        for _ in range(length):
            G[previous,node]=G[node,previous]=-1
            previous=node;node+=1
    return G


def matrix_data(G, omissions):
    kernel=G.nullspace();assert len(kernel)==1
    delta=kernel[0]
    if delta[0]<0:delta=-delta
    delta=s.simplify(delta/min(delta))
    assert all(x>0 for x in delta)
    J=(G+delta*delta.T).inv()
    def R(i,j):
        v=s.zeros(G.rows,1);v[i]+=1/delta[i];v[j]-=1/delta[j]
        return s.simplify((v.T*J*v)[0])
    diagonals=[]
    for omitted in omissions:
        kept=[j for j in range(G.rows) if j!=omitted]
        B=G.extract(kept,kept);inverse=B.inv()
        assert all(B[:k,:k].det()>0 for k in range(1,B.rows+1))
        assert all(x>=0 for x in inverse)
        diagonals.append({j:s.simplify(inverse[a,a]) for a,j in enumerate(kept)})
    i,j=omissions
    L=s.simplify(min(delta[i]**2,delta[j]**2)*R(i,j))
    return delta,R,diagonals,L


def run():
    exceptional=[]
    matrices={
        'A1':(s.Matrix([[2,-2],[-2,2]]),(0,1)),
        'E6':(star([2,2,2]),(2,4)),
        'E7':(star([1,3,3]),(3,6)),
        'E8':(star([1,2,5]),(3,7)),
        'F4':(s.Matrix([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-s.sqrt(2),0],
                        [0,0,-s.sqrt(2),2,-1],[0,0,0,-1,2]]),(0,4)),
        'G2':(s.Matrix([[2,-1,0],[-1,2,-s.sqrt(3)],[0,-s.sqrt(3),2]]),(0,2)),
    }
    for name,(G,omissions) in matrices.items():
        delta,R,diagonals,L=matrix_data(G,omissions);rank=G.rows-1
        baseline=s.Rational(2*rank,9)
        assert L>=baseline
        row=dict(type=name,rank=rank,omissions=list(omissions),
                 gram=[[str(x) for x in line] for line in G.tolist()],
                 kernel=[str(x) for x in delta],weighted_point_energy=str(L),
                 retained_inverse_diagonals=[{str(k):str(v) for k,v in d.items()} for d in diagonals])
        if name not in ('E7','E8'):
            minimum=min(v for d in diagonals for v in d.values())
            assert 2*minimum>=baseline
            row['direct_cost_charge_margin']=str(s.simplify(2*minimum-baseline))
        else:
            expected=[{4},{7}] if name=='E7' else [set(),{8}]
            cheap_records=[]
            for side,d in enumerate(diagonals):
                cheap={j for j,value in d.items() if value<1}
                assert cheap==expected[side]
                for leaf in cheap:
                    assert d[leaf]==s.Rational(1,2)
                    i,j=omissions[side],omissions[1-side]
                    assert delta[i]==delta[j]==2 and delta[leaf]==1
                    values=[R(i,j),R(leaf,j),R(i,leaf)]
                    assert values==[s.Rational(1,2),1,s.Rational(1,2)]
                    u=s.symbols('u',real=True)
                    mixture=s.expand((1-u)*values[0]+u*values[1]-u*(1-u)*values[2])
                    assert mixture==(1+u*u)/2
                    cheap_records.append(dict(side=side,omitted=i,opposite_omitted=j,cheap_leaf=leaf,
                                              resistance_values=[str(x) for x in values],mixture=str(mixture)))
            assert baseline<2
            row['cheap_contact_cases']=cheap_records
            row['cheap_case_combined_energy_and_charge_lower_bound']='2'
        exceptional.append(row)

    # The closed inverse formulas are inherited from the frozen classical proof.
    # The following comparisons are exact for every indicated rank, not samples.
    r=s.symbols('r',positive=True)
    small_A=s.factor(2*r/(r+1)-2*r/9)
    large_A=s.factor(2*r/(r+1)-(2-r/36))
    assert s.cancel(small_A-2*r*(8-r)/(9*(r+1)))==0
    assert s.cancel(large_A-(r-8)*(r+9)/(36*(r+1)))==0
    for rank in range(1,9):assert small_A.subs(r,rank)>=0
    t=s.symbols('t',integer=True,nonnegative=True)
    assert s.cancel(large_A.subs(r,9+t)-(t+1)*(t+18)/(36*(t+10)))==0
    # B2/D2 and B3/D3 low-rank presentations have min diagonal r/4;
    # their loss allowance is already <=2*(r/4).
    low_B=[]
    for rank in (2,3):
        assert s.Rational(rank,2)>=s.Rational(2*rank,9)
        vectors=[]
        v=s.zeros(rank,1);v[0]=v[1]=-1;vectors.append(v)
        for a in range(rank-1):
            v=s.zeros(rank,1);v[a]=1;v[a+1]=-1;vectors.append(v)
        v=s.zeros(rank,1);v[rank-1]=s.sqrt(2);vectors.append(v)
        gram=s.Matrix([[(a.T*b)[0] for b in vectors] for a in vectors])
        _,_,diags,L=matrix_data(gram,(0,rank))
        minimum=min(x for dd in diags for x in dd.values())
        assert minimum==s.Rational(rank,4)
        assert L>=s.Rational(2*rank,9)
        low_B.append(dict(rank=rank,cost_lower_bound=str(minimum),weighted_point_energy=str(L)))
    for rank in range(4,10):assert 2>=s.Rational(2*rank,9)
    B_D_margin=s.expand(2-(s.Rational(5,4)-(10+t)/36))
    C_margin=s.expand(2-(s.Rational(5,2)-5*(10+t)/18))
    assert B_D_margin==t/36+s.Rational(37,36)
    assert C_margin==5*t/18+s.Rational(41,18)
    # Cost total <=4, factor containing the narrow wall keeps its old loss<=2.
    energy=s.Rational(2,9)*(189+t)-2*4-2
    pairing=s.expand(energy/2-2)
    assert pairing==14+t/9

    files=[f'scripts/{TAG}.py']
    for dep in ['fresh_nikulin_narrow_20261015','fresh_nikulin_narrow_20261015_products',
                'fresh_nikulin_narrow_20261015_compact_facet','fresh_dimension225_hostile_20261015']:
        files.extend([f'scripts/{dep}.py',f'results/{dep}.json',f'results/{dep}.md'])
        obj=json.loads((ROOT/'results'/f'{dep}.json').read_text())
        for f,h in obj['sources'].items():assert sha256((ROOT/f).read_bytes()).hexdigest()==h
    return {
        'status':'PASS_EXACT_CONTACT_COST_CHARGE_FOR_DIMENSION_189_CANDIDATE',
        'sources':{f:sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        'exceptional_components':exceptional,
        'classical_charge_margins':{'A_r_le8':str(small_A),'A_r_ge9':str(large_A),
                                    'B_D_r10_plus_t':str(B_D_margin),'C_r10_plus_t':str(C_margin)},
        'B_low_rank_costs':low_B,
        'non_narrow_loss_per_unit_schur_cost':2,
        'total_two_endpoint_schur_cost_at_most':4,
        'separate_narrow_component_loss_at_most':2,
        'total_energy_lower_bound_divided_by_M':'2N/9-10',
        'negative_pairing_lower_bound_divided_by_M':'N/9-7',
        'N189_plus_t_pairing':str(pairing),
        'ideal_narrow_facet_dimension_at_most':189,
        'compact_facet_dimension_at_most':143,
        'combined_candidate_dimension_at_most':189,
        'first_excluded_dimension':190,
        'scope':'Exact exceptional matrices, inverse diagonals and resistance identities plus all-rank scalar checks; inherited classical localization and geometric framing arguments explicitly bound. No finite-angle cutoff or label enumeration. Research proof candidate, no novelty/current-best assertion.'
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();data=run();destination=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert data==json.loads(destination.read_text())
    else:destination.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
