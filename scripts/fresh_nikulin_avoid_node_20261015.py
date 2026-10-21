#!/usr/bin/env python3
"""Exact fixed-omission/contact-cost certificate for the narrow affine factor.

No finite-label enumeration. Frozen geometric and all-rank localization inputs
are explicitly bound; new exceptional and low-rank data are rebuilt here.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_avoid_node_20261015'


def star(arms):
    G=2*s.eye(1+sum(arms));node=1
    for length in arms:
        previous=0
        for _ in range(length):
            G[node,previous]=G[previous,node]=-1
            previous=node;node+=1
    return G


def factor(G, omissions):
    null=G.nullspace();assert len(null)==1
    delta=null[0]
    if delta[0]<0:delta=-delta
    delta=s.simplify(delta/min(delta));assert all(v>0 for v in delta)
    regular=(G+delta*delta.T).inv()
    def R(i,j):
        v=s.zeros(G.rows,1);v[i]+=1/delta[i];v[j]-=1/delta[j]
        return s.simplify((v.T*regular*v)[0])
    diagonals=[]
    for omitted in omissions:
        kept=[j for j in range(G.rows) if j!=omitted]
        B=G.extract(kept,kept)
        assert all(B[:k,:k].det()>0 for k in range(1,B.rows+1))
        J=B.inv();assert all(v>=0 for v in J)
        diagonals.append({j:s.simplify(J[a,a]) for a,j in enumerate(kept)})
    i,j=omissions
    L=s.simplify(min(delta[i]**2,delta[j]**2)*R(i,j))
    return delta,R,diagonals,L


def classical(typ,r):
    roots=[];v=s.zeros(r,1)
    if typ=='B':v[0]=v[1]=-1
    else:v[0]=-s.sqrt(2)
    roots.append(v)
    for i in range(r-1):
        v=s.zeros(r,1);v[i]=1;v[i+1]=-1;roots.append(v)
    v=s.zeros(r,1);v[-1]=s.sqrt(2);roots.append(v)
    return s.Matrix([[(a.T*b)[0] for b in roots] for a in roots])


def run():
    matrices={
        'A1':s.Matrix([[2,-2],[-2,2]]),
        'E6':star([2,2,2]),'E7':star([1,3,3]),'E8':star([1,2,5]),
        'F4':s.Matrix([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-s.sqrt(2),0],
                       [0,0,-s.sqrt(2),2,-1],[0,0,0,-1,2]]),
        'G2':s.Matrix([[2,-1,0],[-1,2,-s.sqrt(3)],[0,-s.sqrt(3),2]])}
    records=[];cheap=[]
    for name,G in matrices.items():
        r=G.rows-1;target=s.Rational(2*(r-1),9)
        for forbidden in range(G.rows):
            if name=='A1':omissions=(1-forbidden,1-forbidden)
            elif name=='E6':omissions=tuple(j for j in (2,4,6) if j!=forbidden)[:2]
            elif name=='E7':omissions=(2,5) if forbidden in (4,7) else (4,7)
            elif name=='E8':omissions=(1,6) if forbidden in (3,7) else (3,7)
            elif name=='F4':omissions=(0,4) if forbidden in (1,3) else (1,3)
            else:omissions=tuple(j for j in range(3) if j!=forbidden)
            assert forbidden not in omissions
            delta,R,diags,L=factor(G,omissions);assert L>=target
            row=dict(type=name,rank=r,forbidden=forbidden,omissions=list(omissions),
                     gram=[[str(v) for v in a] for a in G.tolist()],
                     kernel=list(map(str,delta)),point_energy=str(L),target=str(target),
                     inverse_diagonals=[{str(k):str(v) for k,v in a.items()} for a in diags])
            minimum=min(v for a in diags for v in a.values())
            if name!='E8':
                assert 2*minimum>=target
                row['proof_case']='point energy if C=0; inverse-diagonal charge if C>0'
                row['charge_margin']=str(s.simplify(2*minimum-target))
            else:
                # If C>=target/2 the nonnegative energy suffices. Otherwise
                # there is exactly one contact on the following cheap list.
                limit=target/2
                cheap_nodes=[(side,j) for side,dd in enumerate(diags) for j,v in dd.items() if v<limit]
                expected=[(1,7),(1,8)] if omissions==(1,6) else [(1,8)]
                assert cheap_nodes==expected
                assert 2*minimum>=limit  # two retained contacts are impossible
                details=[]
                for side,leaf in cheap_nodes:
                    i,j=omissions[side],omissions[1-side]
                    vals=[R(i,j),R(leaf,j),R(i,leaf)]
                    u=s.symbols('u',real=True)
                    mixture=s.expand((1-u)*vals[0]+u*vals[1]-u*(1-u)*vals[2])
                    if omissions==(3,7):
                        assert list(delta[k] for k in (i,j,leaf))==[2,2,1]
                        assert vals==[s.Rational(1,2),1,s.Rational(1,2)]
                        assert mixture==(1+u*u)/2
                        zero_a_energy=1;positive_a_energy=2
                    else:
                        assert delta[i]==delta[j]==3 and delta[leaf]==9-leaf
                        assert vals==([s.Rational(2,9),s.Rational(7,18),s.Rational(1,6)] if leaf==7
                                      else [s.Rational(2,9),s.Rational(8,9),s.Rational(2,3)])
                        assert mixture==(s.Rational(2,9)+u*u/6 if leaf==7
                                         else s.Rational(2,9)+2*u*u/3)
                        zero_a_energy=s.simplify(min(delta[leaf]**2,delta[j]**2)*vals[1])
                        positive_a_energy=2
                    cost_min=diags[side][leaf]
                    assert zero_a_energy+2*cost_min>=target and positive_a_energy>=target
                    rec=dict(omissions=list(omissions),side=side,leaf=leaf,
                             kernel_values=[str(delta[k]) for k in (i,j,leaf)],
                             cost_coefficient=str(cost_min),resistances=list(map(str,vals)),
                             mixture=str(mixture),zero_omitted_contact_energy_lower=str(zero_a_energy),
                             positive_omitted_contact_energy_lower=str(positive_a_energy))
                    details.append(rec)
                    if rec not in cheap:cheap.append(rec)
                row['proof_case']='C>=7/9 or exactly one listed cheap contact'
                row['cheap_cases']=details
            records.append(row)
    assert len(records)==34 and len(cheap)==3

    # Only endpoint-avoid modifications are new in the classical families.
    # Low ranks are literal matrices. Larger ranks use the bound closed forms.
    low=[]
    for typ in ('B','C'):
        for r in range(2,10):
            G=classical(typ,r)
            omissions=(r-1,0) if typ=='B' else (1,r)
            delta,R,diags,L=factor(G,omissions);target=s.Rational(2*(r-1),9)
            assert L>=target
            minimum=min(v for dd in diags for v in dd.values())
            assert minimum>=s.Rational(1,2)
            rec=dict(type=typ,rank=r,omissions=list(omissions),
                     gram=[[str(v) for v in a] for a in G.tolist()],kernel=list(map(str,delta)),
                     point_energy=str(L),target=str(target),
                     inverse_diagonals=[{str(k):str(v) for k,v in dd.items()} for dd in diags])
            if r<=5:
                assert 2*minimum>=target
                rec['proof_case']='direct minimum-diagonal charge'
            else:
                i,j=omissions;leaf=r if typ=='B' else 0
                assert [(side,k) for side,dd in enumerate(diags) for k,v in dd.items() if v<1]==[(0,leaf)]
                assert diags[0][leaf]==s.Rational(1,2)
                denom=4 if typ=='B' else 2
                vals=[R(i,j),R(leaf,j),R(i,leaf)]
                assert vals==[s.Rational(r-1,denom),s.Rational(r,denom),s.Rational(1,denom)]
                u=s.symbols('u',real=True)
                mixture=s.expand((1-u)*vals[0]+u*vals[1]-u*(1-u)*vals[2])
                assert mixture==(r-1+u*u)/denom
                assert s.Rational(r-1,denom)>=target and min(delta)>=1
                rec['proof_case']='C>=target/2 or sole isolated A1 contact'
                rec['cheap_leaf']=leaf;rec['resistances']=list(map(str,vals));rec['mixture']=str(mixture)
            low.append(rec)
    assert len(low)==16

    # For r>=10, the modified retained diagram has an isolated A1 and
    # the remaining finite B/D inverse diagonals are >=1, hence C>0=>2C>=1.
    t=s.symbols('t',nonnegative=True,integer=True);r=10+t
    B=s.expand((r-5)/4+1-s.Rational(2,9)*(r-1))
    C=s.expand((r-5)/2+1-s.Rational(2,9)*(r-1))
    assert B==t/36+s.Rational(1,4)
    assert C==5*t/18+s.Rational(3,2)
    pointB=s.expand((r-1)/4-s.Rational(2,9)*(r-1))
    pointC=s.expand((r-1)/2-s.Rational(2,9)*(r-1))
    assert pointB==t/36+s.Rational(1,4)
    assert pointC==5*t/18+s.Rational(5,2)
    energy=s.Rational(2,9)*(181+t)-2*4-s.Rational(2,9)
    pairing=s.expand(energy/2-2)
    assert pairing==14+t/9

    # The unavoidable local A1 loss is exhibited with both point contacts
    # on the only node other than the forbidden one: E=C=0 but r=1.
    assert records[0]['point_energy']=='0' and records[1]['point_energy']=='0'
    files=[f'scripts/{TAG}.py']
    for dep in ['fresh_nikulin_dimension_189_candidate_20261015',
                'fresh_nikulin_energy_charge_20261015','fresh_nikulin_energy_charge_20261015_review',
                'fresh_nikulin_narrow_20261015_products','fresh_nikulin_narrow_20261015_compact_facet']:
        files.extend([f'scripts/{dep}.py',f'results/{dep}.json',f'results/{dep}.md'])
        obj=json.loads((ROOT/f'results/{dep}.json').read_text())
        for name,h in obj['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    return dict(status='PASS_EXACT_FIXED_AVOID_NODE_CHARGE_FOR_DIMENSION_181_CANDIDATE',
                sources={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
                exceptional_avoid_node_records=records,E8_distinct_cheap_cases=cheap,
                modified_classical_low_ranks=low,
                all_rank_margins={'B_modified_r10_plus_t':str(B),'C_modified_r10_plus_t':str(C),
                                  'B_modified_point_r10_plus_t':str(pointB),'C_modified_point_r10_plus_t':str(pointC)},
                narrow_component_lemma='energy/M + 2C >= 2(r-1)/9',
                nonnarrow_component_lemma='energy/M + 2C >= 2r/9',
                total_two_endpoint_schur_cost_at_most=4,
                unavoidable_uniform_avoid_node_penalty='2/9',
                total_energy_lower_bound_divided_by_M='2N/9-8-2/9',
                negative_pairing_lower_bound_divided_by_M='(N-55)/9',
                N181_plus_t_pairing=str(pairing),ideal_narrow_facet_dimension_at_most=181,
                compact_facet_dimension_at_most=143,combined_candidate_dimension_at_most=181,
                first_excluded_dimension=182,
                scope='All affine types and every forbidden node. Omissions preset independently of contact rows. No label cutoff; inherited classical localization and reviewed geometric/source inputs are bound. New lemma and resulting research proof candidate await separate independent review; no novelty/current-best claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();destination=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(destination.read_text())
    else:destination.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
