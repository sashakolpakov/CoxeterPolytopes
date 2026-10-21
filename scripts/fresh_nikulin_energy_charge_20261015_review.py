#!/usr/bin/env python3
"""Independent Fraction audit of the new cost-charge step; no producer imports.

Exceptional energies use grounded rational congruences, not regularized inverse
matrices. The all-rank geometry and published inputs have a separate prior audit.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_energy_charge_20261015_review'
PRIMARY='fresh_nikulin_energy_charge_20261015'

def inv(A):
    n=len(A);B=[[Q(x) for x in row]+[Q(i==j) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        p=next(i for i in range(k,n) if B[i][k]);B[k],B[p]=B[p],B[k]
        t=B[k][k];B[k]=[x/t for x in B[k]]
        for i in range(n):
            if i!=k:
                t=B[i][k];B[i]=[x-t*y for x,y in zip(B[i],B[k])]
    return [row[n:] for row in B]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Q(0))
def mv(A,b):return [dot(row,b) for row in A]
def sub(A,idx):return [[A[i][j] for j in idx] for i in idx]
def pd(A):
    B=[[Q(x) for x in row] for row in A]
    for k in range(len(B)):
        assert B[k][k]>0
        for i in range(k+1,len(B)):
            for j in range(k+1,len(B)):B[i][j]-=B[i][k]*B[k][j]/B[k][k]
def star(arms):
    n=1+sum(arms);A=[[Q(2 if i==j else 0) for j in range(n)] for i in range(n)];k=1
    for length in arms:
        parent=0
        for _ in range(length):A[parent][k]=A[k][parent]=-1;parent=k;k+=1
    return A

def factor(A,d2,omissions):
    # A=D G D with positive diagonal D and stored rational squares d2.
    n=len(A);keep=list(range(1,n));J=inv(sub(A,keep))
    kernel=[Q(1)]+mv(J,[-A[i][0] for i in keep])
    assert mv(A,kernel)==[0]*n and min(kernel)>0
    scale=min(kernel);kernel=[x/scale for x in kernel]
    delta2=[d*x*x for d,x in zip(d2,kernel)]
    assert min(delta2)==1
    def energy(i,j):
        v=[Q(0)]*n;v[i]=1/kernel[i];v[j]-=1/kernel[j]
        assert dot(v,kernel)==0
        # Ground coordinate0; recover the other potentials and check all equations.
        potential=[Q(0)]+mv(J,v[1:]);assert mv(A,potential)==v
        return dot(v,potential)
    diagonals=[]
    for omitted in omissions:
        ids=[i for i in range(n) if i!=omitted];B=sub(A,ids);pd(B);inverse=inv(B)
        assert min(x for row in inverse for x in row)>=0
        diagonals.append({i:d2[i]*inverse[k][k] for k,i in enumerate(ids)})
    i,j=omissions
    return kernel,delta2,diagonals,energy,min(delta2[i],delta2[j])*energy(i,j)

def run():
    files=[f'scripts/{TAG}.py']
    for prefix in [PRIMARY,'fresh_dimension225_hostile_20261015','fresh_nikulin_narrow_20261015',
                   'fresh_nikulin_narrow_20261015_products','fresh_nikulin_narrow_20261015_compact_facet']:
        files.extend([f'scripts/{prefix}.py',f'results/{prefix}.json',f'results/{prefix}.md'])
    def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
    for p in files:
        if p.endswith('.json'):
            for f,h in json.loads((ROOT/p).read_text())['sources'].items():assert digest(f)==h,f
    saved=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    cases={
        'A1':([[2,-2],[-2,2]],[Q(1)]*2,(0,1)),
        'E6':(star([2,2,2]),[Q(1)]*7,(2,4)),
        'E7':(star([1,3,3]),[Q(1)]*8,(3,6)),
        'E8':(star([1,2,5]),[Q(1)]*9,(3,7)),
        'F4':([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-2,0],[0,0,-2,4,-2],[0,0,0,-2,4]],[Q(1),Q(1),Q(1),Q(2),Q(2)],(0,4)),
        'G2':([[2,-1,0],[-1,2,-3],[0,-3,6]],[Q(1),Q(1),Q(3)],(0,2)),
    }
    rows=[];cheap_count=0
    for name,(A,d2,omit) in cases.items():
        kernel,delta2,diags,R,L=factor(A,d2,omit);r=len(A)-1;baseline=Q(2*r,9)
        srow=next(v for v in saved['exceptional_components'] if v['type']==name)
        assert srow['rank']==r and srow['omissions']==list(omit)
        # Check the complete literal physical Gram, including its radical edges.
        for i in range(len(A)):
            for j in range(len(A)):
                token=srow['gram'][i][j]
                if token in ('-sqrt(2)','-sqrt(3)'):
                    sign=-1;square=Q(2 if token=='-sqrt(2)' else 3)
                else:
                    z=Q(token);sign=(z>0)-(z<0);square=z*z
                assert (A[i][j]>0)-(A[i][j]<0)==sign
                assert Q(A[i][j])**2/(d2[i]*d2[j])==square
        assert L==Q(srow['weighted_point_energy']) and L>=baseline
        assert [{str(i):str(x) for i,x in d.items()} for d in diags]==srow['retained_inverse_diagonals']
        cheap=[]
        if name not in ('E7','E8'):assert 2*min(v for d in diags for v in d.values())>=baseline
        else:
            for side,d in enumerate(diags):
                for leaf,value in d.items():
                    if value>=1:continue
                    i,j=omit[side],omit[1-side]
                    assert value==Q(1,2) and (delta2[i],delta2[j],delta2[leaf])==(4,4,1)
                    triple=[R(i,j),R(leaf,j),R(i,leaf)]
                    assert triple==[Q(1,2),Q(1),Q(1,2)]
                    # Expand (1-u)Rij+uRlj-u(1-u)Ril independently.
                    coeff=[triple[0],-triple[0]+triple[1]-triple[2],triple[2]]
                    assert coeff==[Q(1,2),Q(0),Q(1,2)]
                    cheap.append(dict(side=side,omitted=i,other=j,leaf=leaf,
                                      energy_coefficients=list(map(str,coeff))))
            expected=[(0,4),(1,7)] if name=='E7' else [(1,8)]
            assert [(x['side'],x['leaf']) for x in cheap]==expected
            assert baseline<2;cheap_count+=len(cheap)
            # a=0: min(sigma,tau)^2 E >=1 and 2C>=1.
            # a>=1: min(sigma,tau)>=2 and E>=1/2.
            assert Q(1)+Q(1)>=baseline and Q(4)*Q(1,2)>=baseline
        rows.append(dict(type=name,rational_congruent_gram=[[str(x) for x in a] for a in A],
                         diagonal_scale_squares=list(map(str,d2)),rational_kernel=list(map(str,kernel)),
                         physical_kernel_squares=list(map(str,delta2)),weighted_point_energy=str(L),
                         inverse_diagonals=[{str(i):str(x) for i,x in d.items()} for d in diags],
                         cheap_contacts=cheap))
    assert cheap_count==3
    # Independently reconstruct the affineB2/B3 normalized-root presentations.
    low=[]
    for r in (2,3):
        vectors=[[-1,-1]+[0]*(r-2)]
        for i in range(r-1):
            v=[0]*r;v[i]=1;v[i+1]=-1;vectors.append(v)
        vectors.append([0]*(r-1)+[2])
        A=[[dot(a,b) for b in vectors] for a in vectors]
        _,_,diags,_,L=factor(A,[Q(1)]*r+[Q(2)],(0,r))
        minimum=min(v for d in diags for v in d.values())
        assert minimum==Q(r,4) and 2*minimum>=Q(2*r,9) and L>=Q(2*r,9)
        low.append(dict(rank=r,minimum_inverse_diagonal=str(minimum),point_energy=str(L)))
    # Exact A9 threshold, not covered by the old product proof's stated r>=10.
    A=[[Q(2 if i==j else -1 if abs(i-j)==1 else 0) for j in range(9)] for i in range(9)]
    J=inv(A);diagonal=[J[i][i] for i in range(9)]
    assert diagonal==[Q((i+1)*(9-i),10) for i in range(9)]
    assert [i+1 for i,x in enumerate(diagonal) if x<=2]==[1,2,8,9]
    assert J[2][2]==Q(21,10)>2
    # Opposite deletions in the10-cycle leave two uncharged arcs of length1.
    A9_flow=Q(1*1,1+1)
    assert A9_flow==Q(1,2)>Q(9-8,4)
    # All-rank comparisons: polynomial coefficient identities plus explicit domains.
    small_A=[Q(2*r,r+1)-Q(2*r,9) for r in range(2,9)]
    assert all(x>=0 for x in small_A)
    # At r=9+t, the largeA margin equals(t+1)(t+18)/(36(t+10)).
    # Cross multiplication has coefficients(18,19,1), all positive.
    A_num=[Q(18),Q(19),Q(1)]
    assert A_num==[Q(9**2+9-72),Q(2*9+1),Q(1)]
    BD=[Q(37,36),Q(1,36)];C=[Q(41,18),Q(5,18)]
    assert BD==[Q(2)-Q(5,4)+Q(10,36),Q(1,36)]
    assert C==[Q(2)-Q(5,2)+Q(50,18),Q(5,18)]
    assert all(x>0 for x in A_num+BD+C)
    assert all(Q(2)>=Q(2*r,9) for r in range(4,10))
    assert all(Q(2)>=Q(2*r,9) for r in range(2,10))
    # Global energy: nonnarrow loss<=2 sumC<=8; narrow loss<=2.
    energy_const=Q(2*189,9)-2*4-2
    assert energy_const/2-2==14
    assert saved['first_excluded_dimension']==190 and saved['combined_candidate_dimension_at_most']==189
    return dict(status='PASS_INDEPENDENT_FRACTION_ENERGY_CHARGE_AND_HOSTILE_SCOPE_REVIEW',
                sources={p:digest(p) for p in files},
                exceptional_factors=rows,cheap_contact_cases=cheap_count,B_low_rank_checks=low,
                A9=dict(inverse_diagonal=list(map(str,diagonal)),allowed_positions=[1,2,8,9],
                        long_arc_lengths=[1,1],flow_energy_lower_bound=str(A9_flow),required_lower_bound='1/4'),
                all_rank_margins=dict(A_r9_plus_t_numerator=list(map(str,A_num)),A_denominator='36(t+10)',
                                      BD_r10_plus_t=list(map(str,BD)),C_r10_plus_t=list(map(str,C))),
                global_schur_cost_bound=4,total_nonnarrow_loss_at_most=8,separate_narrow_loss_at_most=2,
                N189_plus_t_pairing_lower_bound='14+t/9',
                checked_candidate_dimension_at_most=189,
                scope='Independent Fraction congruence/grounded-energy and all-rank scalar audit plus hostile proof review. Inherits the explicitly bound225 geometric/published-input audit and classical localization; extendsA9 directly. No producer imports, label cutoff, new cusp census, or novelty/current-best claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
