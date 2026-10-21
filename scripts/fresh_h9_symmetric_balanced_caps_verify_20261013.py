#!/usr/bin/env python3
"""Independent exact arithmetic and geometry replay of the golden cap seed."""
import copy
import json
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/'data/enumerations/fresh_h9_symmetric_balanced_caps_20261013.json'
OUT=ROOT/'results/fresh_h9_symmetric_balanced_caps_verify_20261013.json'
Z=(Q(0),Q(0));O=(Q(1),Q(0));T=(Q(0),Q(1))
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):return(x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def sc(x,c):return(x[0]*c,x[1]*c)
def sub(x,y):return add(x,sc(y,-1))
def conj(x):return(x[0]+x[1],-x[1])
def sign(x):
    # Twice x = A+B sqrt5; all comparisons are rational.
    a,b=2*x[0]+x[1],x[1]
    if not b:return(a>0)-(a<0)
    if not a:return(b>0)-(b<0)
    if a*b>0:return(a>0)-(a<0)
    d=a*a-5*b*b
    assert d!=0
    return ((a>0)-(a<0))*((d>0)-(d<0))
def dec(x):return tuple(Q(t)for t in x)
def vec(x):return tuple(map(dec,x))
def dot(x,y):
    a=Z
    for b,c in zip(x,y):a=add(a,mul(b,c))
    return a
def mv(g,x):return tuple(dot(row,x)for row in g)


def check(d):
    G=tuple(vec(row)for row in d['gram']);assert len(G)==18
    b=sc(add(O,sc(T,2)),Q(1,2));w=add(sc(O,2),b);c=sc(add(O,T),2);p=add(O,T)
    for i in range(18):
        for j in range(18):
            if i==j:value=sc(O,2)
            elif max(i,j)<2:value=sc(c,-1)
            elif min(i,j)<2:
                center=min(i,j);leaf=max(i,j)-2
                value=sc(p,-1)if center==leaf%2 else Z
            else:
                ii,jj=i-2,j-2
                value=sc(w,-1)if ii//2==jj//2 else(sc(b,-1)if ii%2!=jj%2 else Z)
            assert G[i][j]==value
    # Seven positive standard eigenvalues4. Remaining Gram blocks are P,Q.
    lp=sc(b,-8);lm=add(sc(O,4),sc(b,8))
    cp=sub(sc(O,2),c);cm=add(sc(O,2),c)
    assert sub(mul(cp,lp),sc(mul(p,p),8))==Z
    assert sign(add(cp,lp))<0
    qdet=sub(mul(cm,lm),sc(mul(p,p),8))
    assert qdet==add(sc(O,32),sc(T,40)) and sign(qdet)>0 and sign(cm)>0
    # The conjugate P rank1 and Q rank2 are both positive definite on their images.
    assert sign(conj(add(cp,lp)))>0 and sign(conj(qdet))>0 and sign(conj(cm))>0
    S=vec(d['strict_center_coefficients']);assert S==(Z,Z)+tuple(O for _ in range(16))
    GS=mv(G,S);assert sign(dot(S,GS))<0 and all(sign(x)<0 for x in GS)
    caps=d['caps'];assert len(caps)==70
    assert {tuple(row['plus_indices'])for row in caps}==set(combinations(range(8),4))
    normals=[];covectors=[]
    for row in caps:
        subset=set(row['plus_indices']);r=vec(row['coefficients'])
        expected=[Z,Z]
        for i in range(8):
            e=1 if i in subset else -1
            expected.extend([add(sc(sub(sc(O,2),T),Q(1,8)),sc(T,Q(e,8))),
                             sub(sc(sub(sc(O,2),T),Q(1,8)),sc(T,Q(e,8)))])
        assert r==tuple(expected)
        gr=mv(G,r);assert dot(r,gr)==sc(O,2) and sign(dot(S,gr))<0
        assert gr[:2]==(sc(O,-1),sc(O,-1))
        for i in range(8):
            assert gr[2+2*i:4+2*i]==((Z,sc(T,-1))if i in subset else(sc(T,-1),Z))
        normals.append(r);covectors.append(gr)
    compatible=[]
    for i,j in combinations(range(70),2):
        h=len(set(caps[i]['plus_indices'])^set(caps[j]['plus_indices']))
        pairing=dot(normals[i],covectors[j])
        assert pairing==sub(sc(O,2),sc(mul(T,T),Q(h,4)))
        # Equal norm2 and pairing in Q(phi): finite Coxeter possibilities
        # are0,-1,-phi; the cyclotomic degree formula rules out all others.
        okay=pairing in(Z,sc(O,-1),sc(T,-1))or sign(add(pairing,sc(O,2)))<=0
        assert okay==(h==8)
        if okay:compatible.append((i,j))
    assert len(compatible)==35 and tuple(d['selected_caps'])in compatible
    Y=vec(d['escape_coefficients']);GY=mv(G,Y)
    assert sign(dot(Y,GY))>0 and sign(dot(S,GY))<0
    assert all(sign(x)<=0 for x in GY)
    selected=d['selected_caps']
    assert all(sign(dot(Y,covectors[i]))<0 for i in selected)
    # A small center perturbation turns the feasible spacelike escape strict.
    eps=Q(1,1024);strict=tuple(add(x,sc(y,eps))for x,y in zip(Y,S))
    assert sign(dot(strict,mv(G,strict)))>0
    assert all(sign(x)<0 for x in mv(G,strict))
    assert all(sign(dot(strict,covectors[i]))<0 for i in selected)
    # Replay its S8 transport separately for every complementary pair.
    for i,j in compatible:
        order=sorted(caps[i]['plus_indices'])+sorted(set(range(8))-set(caps[i]['plus_indices']))
        transported=list(strict)
        for source,target in enumerate(order):
            transported[2+2*target]=strict[2+2*source]
            transported[3+2*target]=strict[3+2*source]
        gt=mv(G,transported)
        assert sign(dot(transported,gt))>0 and sign(dot(S,gt))<0
        assert all(sign(x)<0 for x in gt)
        assert sign(dot(transported,covectors[i]))<0 and sign(dot(transported,covectors[j]))<0
    return {'seed_signature':[9,1,8],'conjugate_signature':[10,0,8],
            'individually_compatible_caps':70,'compatible_cap_pairs':35,
            'maximum_compatible_cap_count':2,'selected_cone_walls':20,
            'strict_infinite_volume_escape':True,
            'strict_escapes_for_all_complementary_pairs':35,
            'all_subsets_of_this_70_cap_dictionary_excluded':True}


def main():
    d=json.loads(PATH.read_text());checks=check(d);bad=copy.deepcopy(d)
    bad['caps'][0]['coefficients'][2][0]=str(Q(bad['caps'][0]['coefficients'][2][0])+1)
    try:check(bad)
    except AssertionError:mutation='REJECTED'
    else:raise AssertionError('Forged cap accepted')
    result={'status':'PASS','arithmetic':'stdlib Fraction over Q(phi)',
            'checks':checks,'forged_cap':mutation}
    OUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
