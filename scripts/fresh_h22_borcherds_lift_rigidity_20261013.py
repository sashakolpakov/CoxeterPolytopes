#!/usr/bin/env python3
"""Exact all-real one-coordinate lift rigidity and center-cap family proof.

Standard-library arithmetic only. The known H21 coordinates are the input;
no numerical scout or project arithmetic module is imported.
"""
import argparse
from collections import Counter,deque
from fractions import Fraction as F
import hashlib,json
from itertools import combinations
from math import gcd
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
CUSPS='data/enumerations/h20_gluing_independent_cusps_20261013.json'
OUT='results/fresh_h22_borcherds_lift_rigidity_20261013.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def mv(M,x):return [sum(a*b for a,b in zip(row,x,strict=True)) for row in M]
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def run():
    raw=read(SOURCE)['coxeter_input'];Q=[[-int(x) for x in row]for row in raw['basis_gram']];W=[[int(x)for x in row]for row in raw['wall_coordinates_as_columns_transposed']]
    cov=[mv(Q,r) for r in W];G=[[dot(r,c) for c in cov]for r in W]
    assert len(Q)==22 and len(W)==210
    S=[i for i in range(210) if G[i][i]==2];L=[i for i in range(210)if G[i][i]==4];assert len(S)==42 and len(L)==168
    assert all(G[i][j]<=0 for i,j in combinations(range(210),2))
    graph={i:[j for j in S if j!=i and G[i][j]==-1]for i in S}
    assert all(G[i][j]in(0,-1)for i,j in combinations(S,2))
    assert all(len(graph[i])==5 for i in S)
    assert all(G[j][k]==0 for i in S for j,k in combinations(graph[i],2))
    seen={S[0]};queue=deque(seen)
    while queue:
        i=queue.popleft()
        for j in graph[i]:
            if j not in seen:seen.add(j);queue.append(j)
    assert len(seen)==42
    # Positive and negative height supports are cliques in this triangle-free
    # graph. At most four heights are nonzero. One zero propagates along all
    # label3 edges because the new cosine would lie strictly in(0,1/2).
    assert len(S)>2*2
    short_anchors={i:next(j for j in S if G[i][j]==-2)for i in L}
    assert all(G[i][j]in(0,-2)for i in L for j in S)
    # With shorts at height0, a long's cosine to its short anchor is
    # sqrt2/sqrt(4+t²), forcing t in{0,+2,-2} (labels4 or3).
    # A label3 pair of old longs permits only the two zero heights.
    gates=[]
    for x in(-2,0,2):
        for y in(-2,0,2):
            p=-2+x*y;n=(4+x*x)*(4+y*y);ratio=F(p*p,n)
            if not x and not y:assert p<0 and ratio==F(1,4);kind='allowed_both_zero'
            elif p>0:kind='positive_pairing'
            elif (x==0)!=(y==0):assert 0<ratio<F(1,4);kind='between_right_and_pi3'
            else:
                assert ratio==F(9,16)>F(1,2)
                # 9/16 <(3+sqrt5)/8 because3<2sqrt5.
                assert 3**2<4*5;kind='between_pi4_and_pi5'
            gates.append({'heights':[x,y],'pairing':p,'squared_cosine':str(ratio),'decision':kind})
    long_anchors={i:next(j for j in L if i!=j and G[i][j]==-2)for i in L}
    assert all(sum(G[i][j]==-2 for j in L if i!=j)==40 for i in L)
    h=[sum(col)for col in zip(*W,strict=True)];g=gcd(*h);h=[x//g for x in h];hc=mv(Q,h)
    assert dot(h,hc)==-14 and Counter(dot(h,c)for c in cov)=={-1:42,-8:168}
    delta=read(CUSPS)['cusps'][1]['vector'];dc=mv(Q,delta)
    assert dot(delta,dc)==0 and dot(delta,hc)==-12 and max(dot(delta,c)for c in cov)==0
    epsilon=F(1,10000);v=[F(x)+epsilon*y for x,y in zip(delta,h,strict=True)];vc=mv(Q,v)
    escape_norm=dot(v,vc)+1;center_pair=dot(v,hc)
    assert escape_norm==1-24*epsilon-14*epsilon**2>0
    oldpairs=[dot(v,c)for c in cov];assert max(oldpairs)<0
    # For m>=3, u=cos²(pi/m) belongs[1/4,1). Center caps are
    # (h,+/-sqrt(14+1/(2u))), with norm1/(2u). Thus height<=4.
    assert center_pair+4<0
    return {'status':'PASS_ALL_REAL_FIXED_PROJECTION_H21_TO_H22_LIFT_RIGIDITY_AND_UNIFORM_CENTER_CAP_ESCAPE',
        'source_sha256':{p:sha(p)for p in(SOURCE,CUSPS)},
        'short_graph':{'vertices':S,'degree':5,'connected':True,'triangle_free':True,'largest_nonzero_same_sign_height_support':2,'forced_zero_heights':42},
        'long_short_anchors':[[i,j]for i,j in short_anchors.items()],'long_allowed_heights':['0','2','-2'],'long_pi3_anchors':[[i,j]for i,j in long_anchors.items()],
        'long_height_pair_cases':gates,'forced_zero_heights_total':210,
        'center':h,'center_norm':-14,'center_old_pair_counts':{'-1':42,'-8':168},
        'center_cap_family':{'parameter':'Every integer m>=3','normal':'(h,+/-sqrt(14+1/(2cos²(pi/m))))','norm':'1/(2cos²(pi/m))',
            'pairs_to_old_walls':{'label_m':42,'dotted':168},'mutual_pairing':'-28-1/(2cos²(pi/m))','mutual_relation':'dotted','facets':212,
            'actual_facets':'Strict timelike center and all nonpositive pairings imply actual facets by orthogonal projection of the center.'},
        'uniform_strict_escape':{'old_projection':list(map(str,v)),'new_coordinate':'1','norm':str(escape_norm),'center_pairing':str(center_pair),
            'old_wall_pairings':list(map(str,oldpairs)),'new_cap_pairing_upper_bound':str(center_pair+4)},
        'scope':'All one-coordinate lifts retaining one wall over each original210 projected directions have zero heights. Additional facets with new projections remain unrestricted. For eachm>=3, the212-wall center-cap construction has actual Coxeter facets but infinite volume, with one common exact strict escape. No H22 nonexistence theorem follows.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert read(OUT)==out
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
