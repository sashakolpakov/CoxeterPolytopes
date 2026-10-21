#!/usr/bin/env python3
"""Independent Golden-class replay of the complete finite leaf-angle census."""
import copy
import json
from itertools import combinations,combinations_with_replacement
from pathlib import Path
from fractions import Fraction as F
from coxeter_polytopes.exact.golden import Golden as K

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/'data/enumerations/fresh_h9_symmetric_finite_caps_20261013.json'
OUT=ROOT/'results/fresh_h9_symmetric_finite_caps_verify_20261013.json'
Z,O,T=K(0),K(1),K(0,1)
dec=lambda x:K(*map(F,x))
vec=lambda x:tuple(map(dec,x))
dot=lambda x,y:sum((a*b for a,b in zip(x,y,strict=True)),Z)
mv=lambda g,x:tuple(dot(row,x)for row in g)

def radical_sign(a,b,d):
    assert d>=0
    if not d or not b:return a.sign()
    if not a:return b.sign()
    if a.sign()==b.sign():return a.sign()
    difference=a*a-b*b*d
    return a.sign()*difference.sign()

def check(data):
    source=json.loads((ROOT/data['source']).read_text());G=tuple(vec(r)for r in source['gram'])
    b=F(1,2)+T;p=1+T;w=2+b;c=2*p
    for i in range(18):
        for j in range(18):
            expected=K(2)if i==j else -c if max(i,j)<2 else (
                -p if min(i,j)<2 and min(i,j)==(max(i,j)-2)%2 else Z if min(i,j)<2 else
                -w if (i-2)//2==(j-2)//2 else -b if i%2!=j%2 else Z)
            assert G[i][j]==expected
    basis=(0,1,*range(2,18,2));assert list(basis)==data['basis_indices']
    H=tuple(tuple(G[i][j]for j in basis)for i in basis);I=tuple(vec(r)for r in data['inverse_basis_gram'])
    assert H==tuple(vec(r)for r in data['basis_gram'])
    assert all(sum((H[i][k]*I[k][j]for k in range(10)),Z)==(O if i==j else Z)for i in range(10)for j in range(10))
    for i in range(8):
        assert all(G[2+2*i][j]+G[3+2*i][j]==G[2][j]+G[3][j]for j in range(18))
    assert all(G[0][j]+G[1][j]==p/b*(G[2][j]+G[3][j])for j in range(18))
    assert dec(data['center_sum_to_one_leaf_pair_sum_ratio'])==p/b
    expected_orbits={(l,r,k)for l,r in combinations_with_replacement((Z,O,T),2)for k in(range(9)if l!=r else(8,))}
    actual={(dec(r['left_leaf_magnitude']),dec(r['right_leaf_magnitude']),r['left_count'])for r in data['orbits']}
    assert actual==expected_orbits and len(actual)==len(data['orbits'])==data['orbit_count']==30
    assert len({l+r for l,r in combinations_with_replacement((Z,O,T),2)})==6
    E=(O,-O,*(Z for _ in range(8)));A=dot(E,mv(I,E));assert A>0
    counts=dict(positive_center=0,forbidden_small_center=0,accepted=0,nonreal=0);accepted=[]
    for index,row in enumerate(data['orbits']):
        left,right,k=dec(row['left_leaf_magnitude']),dec(row['right_leaf_magnitude']),row['left_count']
        total=left+right;cs=-p*total/b;fixed=(Z,cs,*[-(left if i<k else right)for i in range(8)])
        B=2*dot(E,mv(I,fixed));C=dot(fixed,mv(I,fixed))-2;d=B*B-4*A*C
        assert total==dec(row['pair_sum'])and cs==dec(row['center_sum'])
        assert (A,B,C)==vec(row['quadratic'])and d==dec(row['discriminant'])
        if d<0:
            counts['nonreal']+=1;assert not row['solutions'];continue
        branches=[1]if not d else[-1,1]
        assert [r['branch']for r in row['solutions']]==branches
        for solution,branch in zip(row['solutions'],branches,strict=True):
            a=-B/(2*A);bq=K(branch)/(2*A)
            assert vec(solution['first_center'])==(a,bq)
            assert vec(solution['second_center'])==(cs-a,-bq)
            # Check the quadratic identity coefficient by coefficient, not numerically.
            assert A*(a*a+bq*bq*d)+B*a+C==0 and 2*A*a*bq+B*bq==0
            centers=((a,bq),(cs-a,-bq));signs=[radical_sign(u,v,d)for u,v in centers]
            shifted=[radical_sign(u+1,v,d)for u,v in centers]
            assert signs==solution['center_signs']and shifted==solution['center_plus_one_signs']
            if any(s>0 for s in signs):status='positive_center'
            elif any(s<0 and t>0 for s,t in zip(signs,shifted)):status='forbidden_small_center'
            else:
                status='accepted';accepted.append(index)
                assert (left,right,k)==(Z,T,4)and not d and a==-1 and cs-a==-1
                coeff=mv(I,tuple(x+a*y for x,y in zip(fixed,E)))
                assert dot(coeff,mv(H,coeff))==2
            assert status==solution['status'];counts[status]+=1
    assert counts==data['solution_counts']and accepted==data['accepted_orbit_indices']
    assert counts['accepted']==1 and data['new_caps']==0
    patterns=set()
    for row in source['caps']:
        r=vec(row['coefficients']);cov=mv(G,r)
        assert dot(r,cov)==2 and cov[:2]==(-O,-O)
        assert all(cov[2+2*i:4+2*i]in((Z,-T),(-T,Z))for i in range(8))
        patterns.add(tuple(i for i in range(8)if cov[2+2*i]==Z))
    assert patterns==set(combinations(range(8),4))
    assert data['all_accepted_mirrors_already_in_source']==list(range(70))
    return {'rank':10,'exact_orbits':30,'real_orbit_solutions':sum(counts.values()),
        'solution_counts':counts,'all_accepted_caps_already_in_old_bank':70,'new_caps':0}

def main():
    data=json.loads(PATH.read_text());checks=check(data);bad=copy.deepcopy(data)
    bad['orbits'][0]['quadratic'][2][0]=str(F(bad['orbits'][0]['quadratic'][2][0])+1)
    try:check(bad)
    except AssertionError:forged='REJECTED'
    else:raise AssertionError('Corrupted norm equation accepted')
    out={'status':'PASS','arithmetic':'Independent repository Golden implementation; verified inverse certificate and reconstructed every quadratic and radical sign.',
        'checks':checks,'forged_norm_equation':forged}
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
