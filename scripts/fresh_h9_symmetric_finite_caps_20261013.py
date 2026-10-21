#!/usr/bin/env python3
"""Complete finite golden leaf-pairing census on the new symmetric seed."""
import argparse
from itertools import combinations_with_replacement
import json
from fractions import Fraction as Q
import fresh_h9_symmetric_balanced_caps_verify_20261013 as f

ROOT=f.ROOT
OUT=ROOT/'data/enumerations/fresh_h9_symmetric_finite_caps_20261013.json'
Z,O,T=f.Z,f.O,f.T

def inv(x):
    c=f.conj(x);n=f.mul(x,c);assert n[1]==0 and n[0]
    return f.sc(c,1/n[0])

def inverse(g):
    n=len(g);a=[list(row)+[O if i==j else Z for j in range(n)]for i,row in enumerate(g)]
    for i in range(n):
        p=next(j for j in range(i,n)if a[j][i]!=Z);a[p],a[i]=a[i],a[p]
        s=inv(a[i][i]);a[i]=[f.mul(x,s)for x in a[i]]
        for j in range(n):
            if j!=i:
                s=a[j][i];a[j]=[f.sub(x,f.mul(s,y))for x,y in zip(a[j],a[i])]
    return tuple(tuple(row[n:])for row in a)

def rsign(a,b,d):
    """Exact sign of a+b*sqrt(d); works also for square or zero d."""
    assert f.sign(d)>=0
    sa,sb=f.sign(a),f.sign(b)
    if d==Z or not sb:return sa
    if not sa:return sb
    if sa==sb:return sa
    return sa*f.sign(f.sub(f.mul(a,a),f.mul(f.mul(b,b),d)))

def enc(x):return list(map(str,x))
def matenc(m):return[[enc(x)for x in row]for row in m]

def run(verify=False):
    source=json.loads(f.PATH.read_text());G=tuple(f.vec(row)for row in source['gram'])
    basis=(0,1,*range(2,18,2));H=tuple(tuple(G[i][j]for j in basis)for i in basis);I=inverse(H)
    assert all(f.dot(H[i],tuple(I[k][j]for k in range(10)))==(O if i==j else Z)for i in range(10)for j in range(10))
    b=f.add(f.sc(O,Q(1,2)),T);p=f.add(O,T);ratio=f.mul(p,inv(b))
    # The eight leaf pair sums coincide; the center sum equals p/(8b)*S.
    for i in range(8):
        assert tuple(f.add(G[2+2*i][j],G[3+2*i][j])for j in range(18))==tuple(f.add(G[2][j],G[3][j])for j in range(18))
    assert all(f.add(G[0][j],G[1][j])==f.mul(ratio,f.add(G[2][j],G[3][j]))for j in range(18))
    E=(O,f.sc(O,-1),*(Z for _ in range(8)));A=f.dot(E,f.mv(I,E));assert f.sign(A)>0
    rows=[];accepted=[];counts={'positive_center':0,'forbidden_small_center':0,'accepted':0,'nonreal':0}
    sums=[]
    for left,right in combinations_with_replacement((Z,O,T),2):
        total=f.add(left,right);assert total not in sums;sums.append(total)
        center_sum=f.sc(f.mul(ratio,total),-1)
        for k in(range(9)if left!=right else(8,)):
            fixed=(Z,center_sum,*(f.sc(left if i<k else right,-1)for i in range(8)))
            B=f.sc(f.dot(E,f.mv(I,fixed)),2);C=f.sub(f.dot(fixed,f.mv(I,fixed)),f.sc(O,2))
            disc=f.sub(f.mul(B,B),f.sc(f.mul(A,C),4));solutions=[]
            if f.sign(disc)<0:counts['nonreal']+=1
            else:
                a=f.mul(f.sc(B,-1),inv(f.sc(A,2)))
                for branch in((1,)if disc==Z else(-1,1)):
                    c=f.sc(inv(f.sc(A,2)),branch)
                    centers=((a,c),(f.sub(center_sum,a),f.sc(c,-1)))
                    signs=[rsign(x,y,disc)for x,y in centers]
                    if any(v>0 for v in signs):status='positive_center'
                    elif any(v<0 and rsign(f.add(x,O),y,disc)>0 for v,(x,y)in zip(signs,centers)):
                        status='forbidden_small_center'
                    else:
                        # The census leaves only the old balanced orbit.
                        assert left==Z and right==T and k==4 and disc==Z
                        assert centers[0][0]==centers[1][0]==f.sc(O,-1)
                        status='accepted';accepted.append(len(rows))
                    counts[status]+=1
                    solutions.append({'branch':branch,'first_center':[enc(a),enc(c)],
                        'second_center':[enc(centers[1][0]),enc(centers[1][1])],
                        'center_signs':signs,'center_plus_one_signs':[rsign(f.add(x,O),y,disc)for x,y in centers],
                        'status':status})
            rows.append({'left_leaf_magnitude':enc(left),'right_leaf_magnitude':enc(right),
                'left_count':k,'pair_sum':enc(total),'center_sum':enc(center_sum),
                'quadratic':[enc(A),enc(B),enc(C)],'discriminant':enc(disc),'solutions':solutions})
    assert len(rows)==30 and len(accepted)==1 and counts['accepted']==1
    # Every assignment in the accepted orbit is already in the old bank.
    old_matches=[];old_patterns=set()
    for index,cap in enumerate(source['caps']):
        c=f.vec(cap['coefficients']);cov=f.mv(G,c)
        assert cov[:2]==(f.sc(O,-1),f.sc(O,-1)) and f.dot(c,cov)==f.sc(O,2)
        assert sum(cov[2+2*i]==Z for i in range(8))==4
        assert all(cov[2+2*i:4+2*i]in((Z,f.sc(T,-1)),(f.sc(T,-1),Z))for i in range(8))
        old_patterns.add(tuple(cov[2+2*i]==Z for i in range(8)))
        old_matches.append(index)
    assert len(old_matches)==len(old_patterns)==70
    out={'status':'PASS_COMPLETE_DECLARED_LEAF_ALPHABET','source':str(f.PATH.relative_to(ROOT)),
        'leaf_pairing_alphabet':[enc(Z),enc(f.sc(O,-1)),enc(f.sc(T,-1))],
        'basis_indices':list(basis),'basis_gram':matenc(H),'inverse_basis_gram':matenc(I),
        'center_sum_to_one_leaf_pair_sum_ratio':enc(ratio),'orbit_count':len(rows),
        'orbits':rows,'solution_counts':counts,'accepted_orbit_indices':accepted,
        'all_accepted_mirrors_already_in_source':old_matches,'new_caps':0,
        'scope':'Every norm-two real normal with all sixteen leaf pairings in {0,-1,-phi} is covered, allowing arbitrary real center pairings. A necessary Coxeter gate excludes positive center pairings and pairings strictly between -1 and 0. The only remaining orbit is exactly the old 70 balanced caps. Other leaf labels and dotted leaf pairs are not excluded.'}
    if verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':out['status'],'orbits':len(rows),'solution_counts':counts,'new_caps':0}));return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
