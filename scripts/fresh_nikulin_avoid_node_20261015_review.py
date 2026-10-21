#!/usr/bin/env python3
"""Independent rational-congruence/grounded-energy review of avoid-node accounting."""
import argparse,json
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import fresh_nikulin_energy_charge_20261015_review as H
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_avoid_node_20261015_review'
PRIMARY='fresh_nikulin_avoid_node_20261015'

def factor(A,d2,omissions):
    n=len(A);ids=list(range(1,n));J=H.inv(H.sub(A,ids))
    ker=[Q(1)]+H.mv(J,[-A[i][0]for i in ids]);assert H.mv(A,ker)==[0]*n and min(ker)>0
    raw=[d*k*k for d,k in zip(d2,ker)];scale=min(raw);delta2=[x/scale for x in raw]
    def R(i,j):
        v=[Q(0)]*n;v[i]+=1/ker[i];v[j]-=1/ker[j]
        potential=[Q(0)]+H.mv(J,v[1:]);assert H.mv(A,potential)==v
        return scale*H.dot(v,potential)
    diags=[]
    for omit in omissions:
        kept=[i for i in range(n)if i!=omit];B=H.sub(A,kept);H.pd(B);K=H.inv(B)
        assert all(v>=0 for row in K for v in row)
        diags.append({i:d2[i]*K[j][j]for j,i in enumerate(kept)})
    return delta2,diags,R,min(delta2[i]for i in omissions)*R(*omissions)

def check_saved(A,d2,row,diags,L):
    for i in range(len(A)):
        for j in range(len(A)):
            t=row['gram'][i][j]
            if t in ('-sqrt(2)','-sqrt(3)'):sign=-1;sq=Q(2 if t=='-sqrt(2)'else 3)
            else:z=Q(t);sign=(z>0)-(z<0);sq=z*z
            assert sign==(A[i][j]>0)-(A[i][j]<0)
            assert A[i][j]**2/(d2[i]*d2[j])==sq
    assert [{str(k):str(v)for k,v in a.items()}for a in diags]==row['inverse_diagonals']
    assert L==Q(row['point_energy'])

def classical(typ,r):
    vectors=[([-1,-1]+[0]*(r-2))if typ=='B'else([-2]+[0]*(r-1))]
    for i in range(r-1):
        v=[0]*r;v[i]=1;v[i+1]=-1;vectors.append(v)
    vectors.append([0]*(r-1)+[2]);A=[[H.dot(a,b)for b in vectors]for a in vectors]
    return A,([Q(1)]*r+[Q(2)])if typ=='B'else([Q(2)]+[Q(1)]*(r-1)+[Q(2)])

def run():
    saved=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    cases={
      'A1':([[2,-2],[-2,2]],[Q(1)]*2),
      'E6':(H.star([2,2,2]),[Q(1)]*7),'E7':(H.star([1,3,3]),[Q(1)]*8),'E8':(H.star([1,2,5]),[Q(1)]*9),
      'F4':([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-2,0],[0,0,-2,4,-2],[0,0,0,-2,4]],[Q(1),Q(1),Q(1),Q(2),Q(2)]),
      'G2':([[2,-1,0],[-1,2,-3],[0,-3,6]],[Q(1),Q(1),Q(3)])}
    records=[]
    for row in saved['exceptional_avoid_node_records']:
        typ=row['type'];A,d2=cases[typ];omits=row['omissions'];s=row['forbidden'];r=len(A)-1;target=Q(2*(r-1),9)
        assert s not in omits and row['rank']==r and target==Q(row['target'])
        delta2,diags,R,L=factor(A,d2,omits);check_saved(A,d2,row,diags,L);assert L>=target
        minimum=min(v for d in diags for v in d.values());cheap=[]
        if typ!='E8':assert 2*minimum>=target
        else:
            # Above target/2 charge alone pays. Below it there is at most one contact.
            assert 2*minimum>=target/2
            for side,d in enumerate(diags):
                for leaf,value in d.items():
                    if value>=target/2:continue
                    i,j=omits[side],omits[1-side];triple=[R(i,j),R(leaf,j),R(i,leaf)]
                    co=[triple[0],-triple[0]+triple[1]-triple[2],triple[2]]
                    assert co[1]==0 and min(co)>=0
                    # Positive omitted contact => both row heights are at least the omitted-node coefficients.
                    positive=min(delta2[i],delta2[j])*co[0]
                    zero=min(delta2[leaf],delta2[j])*triple[1]
                    assert positive>=target and zero+2*value>=target
                    cheap.append(dict(side=side,leaf=leaf,energy_coefficients=list(map(str,co)),
                                      positive_contact_energy_lower=str(positive),zero_contact_energy_lower=str(zero),
                                      charge_lower=str(2*value)))
            assert {(x['side'],x['leaf'])for x in cheap}==({(1,7),(1,8)}if omits==[1,6]else{(1,8)})
        records.append(dict(type=typ,forbidden=s,omissions=omits,point_energy=str(L),
                            minimum_diagonal=str(minimum),physical_kernel_squares=list(map(str,delta2)),cheap_cases=cheap))
    assert len(records)==34
    low=[]
    for typ in ('B','C'):
        for r in range(2,13):
            A,d2=classical(typ,r);omit=(r-1,0)if typ=='B'else(1,r)
            delta2,diags,R,L=factor(A,d2,omit);target=Q(2*(r-1),9);assert L>=target
            if r<10:
                row=next(x for x in saved['modified_classical_low_ranks']if x['type']==typ and x['rank']==r)
                assert row['omissions']==list(omit);check_saved(A,d2,row,diags,L)
            minimum=min(v for d in diags for v in d.values());assert minimum>=Q(1,2)
            if r<=5:assert 2*minimum>=target
            else:
                leaf=r if typ=='B'else 0
                assert [(side,k)for side,d in enumerate(diags)for k,v in d.items()if v<1]==[(0,leaf)]
                assert diags[0][leaf]==Q(1,2)
                i,j=omit;den=4 if typ=='B'else 2
                vals=[R(i,j),R(leaf,j),R(i,leaf)]
                assert vals==[Q(r-1,den),Q(r,den),Q(1,den)]
                assert -vals[0]+vals[1]-vals[2]==0 and vals[0]>=target
                assert min(delta2)>=1
            low.append(dict(type=typ,rank=r,minimum_diagonal=str(minimum),point_energy=str(L)))
    # Unbounded ranks: inherited modified-support localization is explicit in product proof.
    # These are coefficient identities in t=r-10, not a finite-rank extrapolation.
    B=[Q(5,4)+1-Q(18,9),Q(1,4)-Q(2,9)]
    C=[Q(5,2)+1-Q(18,9),Q(1,2)-Q(2,9)]
    assert B==[Q(1,4),Q(1,36)] and C==[Q(3,2),Q(5,18)]
    assert all(x>0 for x in B+C)
    N=181;energy=Q(2*N,9)-8-Q(2,9);assert energy/2-2==14
    assert saved['combined_candidate_dimension_at_most']==181
    assert all(x['point_energy']=='0'for x in records if x['type']=='A1')
    files=[f'scripts/{TAG}.py','scripts/fresh_nikulin_energy_charge_20261015_review.py']
    for dep in (PRIMARY,'fresh_nikulin_dimension_189_candidate_20261015','fresh_nikulin_narrow_20261015_products','fresh_nikulin_energy_charge_20261015_review'):
        files.extend([f'scripts/{dep}.py',f'results/{dep}.json',f'results/{dep}.md'])
    digest=lambda p:sha256((ROOT/p).read_bytes()).hexdigest()
    for p in files:
        if p.endswith('.json'):
            for f,h in json.loads((ROOT/p).read_text())['sources'].items():assert digest(f)==h,f
    return dict(status='PASS_INDEPENDENT_AVOID_NODE_RATIONAL_ENERGY_AND_PROOF_REVIEW',sources={p:digest(p)for p in files},
       exceptional_cases=records,modified_classical_checks=low,
       all_rank_margins={'B_r10_plus_t':list(map(str,B)),'C_r10_plus_t':list(map(str,C))},
       total_two_endpoint_schur_cost_at_most=4,uniform_avoid_node_loss='2/9',first_excluded_dimension=182,
       combined_candidate_dimension_at_most=181,
       scope='Independent rational-congruence grounded solves, complete inverse-diagonal comparisons and analytic case review. Prior reviewed all-rank localization, actual-ridge framing, common-height identity and compact143 branch are inherited explicitly. Internal research review, not external validation or novelty certification.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
