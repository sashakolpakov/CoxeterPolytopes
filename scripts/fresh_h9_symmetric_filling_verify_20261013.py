#!/usr/bin/env python3
"""Independent rational-polynomial replay of the S8 filling obstruction.

No producer, symbolic library, numerical package, or previous helper imports.
The geometric rank/sign argument is recorded in the companion report.
"""
import copy
import json
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/'data/enumerations/fresh_h9_symmetric_filling_proof_20261013.json'
OUT=ROOT/'results/fresh_h9_symmetric_filling_verify_20261013.json'
ZERO={}; E=(0,)*5


def const(n):return {E:F(n)}if n else {}
def variable(i):
    e=[0]*5;e[i]=1;return {tuple(e):F(1)}
def add(*items):
    d={}
    for item in items:
        for m,c in item.items():d[m]=d.get(m,F(0))+c
    return {m:c for m,c in d.items()if c}
def mul(a,b):
    d={}
    for m,c in a.items():
        for n,e in b.items():
            k=tuple(x+y for x,y in zip(m,n));d[k]=d.get(k,F(0))+c*e
    return {m:c for m,c in d.items()if c}
def scale(a,c):return mul(a,const(c))
def decode(rows):
    d={tuple(m):F(c)for m,c in rows}
    assert len(d)==len(rows) and all(len(m)==5 for m in d)
    return {m:c for m,c in d.items()if c}
def dot(a,b):return add(*(mul(x,y)for x,y in zip(a,b)))
def mv(g,v):return [dot(row,v)for row in g]


def check(data):
    assert data['variables']==['a','b','w','p','q']
    a,b,w,p,q=[variable(i)for i in range(5)]
    one,two=const(1),const(2);delta=add(two,a);total=add(a,b)
    def gram(matched):
        g=[[two if i==j else {} for j in range(17)]for i in range(17)]
        for i in range(8):
            u,v=2*i,2*i+1
            g[u][v]=g[v][u]=scale(matched,-1)
            g[16][u]=g[u][16]=scale(p,-1)
            g[16][v]=g[v][16]=scale(q,-1)
            for j in range(i):
                h,k=2*j,2*j+1
                g[u][h]=g[h][u]=g[v][k]=g[k][v]=scale(a,-1)
                g[u][k]=g[k][u]=g[v][h]=g[h][v]=scale(b,-1)
        return g
    G=gram(w)
    C=[]
    for parity in (0,1):
        for i in range(7):
            v=[{}for _ in range(17)];v[2*i+parity]=one;v[14+parity]=const(-1)
            C.append(v)
    # Seven independent contrasts of each species: the Kronecker factor
    # I7+J7 is positive definite, eigenvalues1 (multiplicity6),8 (multiplicity1).
    for i in range(14):
        for j in range(14):
            factor=2 if i%7==j%7 else 1
            weight=delta if i//7==j//7 else add(b,scale(w,-1))
            assert dot(C[i],mv(G,C[j]))==scale(weight,factor)
    determinant=add(mul(delta,delta),scale(mul(add(b,scale(w,-1)),add(b,scale(w,-1))),-1))
    assert determinant==decode(data['standard_block_determinant'])
    assert data['contrast_rank']==14>data['ambient_rank_bound']==10
    assert {r['branch']for r in data['branches']}=={'plus','minus'}
    S=[one]*16+[{}]
    for row in data['branches']:
        plus=row['branch']=='plus'
        matched=add(two,a,b)if plus else add(b,scale(two,-1),scale(a,-1))
        assert matched==decode(row['matched_weight'])
        # This verifies each root of the standard-block determinant separately.
        assert mul(delta,delta)==mul(add(b,scale(matched,-1)),add(b,scale(matched,-1)))
        H=gram(matched)
        lp=scale(total,-8)if plus else add(const(4),scale(a,-6),scale(b,-8))
        if plus:
            signs=[s for i in range(8)for s in ((1,-1)if i<4 else(-1,1))]+[0]
            Y=[add(scale(delta,int(i<16)),scale(total,4*signs[i]))for i in range(17)]
            norm=scale(mul(mul(total,delta),add(scale(a,3),scale(b,4),const(-2))),128)
            sy=scale(mul(total,delta),-128)
            expected=[]
            for i in range(8):
                sigma=1 if i<4 else -1
                expected.extend([scale(mul(total,delta),8*(-1+sigma)),
                                 scale(mul(total,delta),8*(-1-sigma))])
            expected.append(scale(mul(delta,add(p,q)),-8))
        else:
            Y=[add(total,lp if i in (14,15)else{})if i<16 else{}for i in range(17)]
            norm=scale(mul(mul(lp,delta),add(lp,total)),4)
            sy=scale(mul(lp,delta),4)
            expected=[{}for _ in range(14)]+[scale(mul(lp,delta),2)]*2
            expected.append(scale(mul(delta,add(p,q)),-2))
        assert Y==[decode(x)for x in row['witness_coefficients']]
        pairs=mv(H,Y)
        assert pairs==expected==[decode(x)for x in row['wall_pairings']]
        assert dot(Y,pairs)==norm==decode(row['witness_norm'])
        assert dot(S,pairs)==sy==decode(row['sum_witness_pairing'])
        assert dot(S,mv(H,S))==scale(lp,16)==decode(row['sum_norm'])
    # Exact inequalities used in the report, reduced to nonnegative variables.
    # Plus: a>0 implies a>=1, giving3a+4b-2>=1;
    # b>0 implies b>=1, giving3a+4b-2>=2. No bounded label alphabet.
    # Minus: b>=2+a implies lp<=-12-14a and lp+s<=-10-12a.
    lp=add(const(4),scale(a,-6),scale(b,-8))
    gap=add(b,const(-2),scale(a,-1))
    assert lp==add(const(-12),scale(a,-14),scale(gap,-8))
    assert add(lp,total)==add(const(-10),scale(a,-12),scale(gap,-7))
    return {'symbolic_gram_branches':2,'generic_extra_wall_pairings':True,
            'rank_factorization':'(I7+J7) tensor [[2+a,b-w],[b-w,2+a]]',
            'explicit_escape_norm_identities':2,'bounded_angle_alphabet':False}


def main():
    data=json.loads(ARCHIVE.read_text());checks=check(data);bad=copy.deepcopy(data)
    bad['branches'][0]['witness_norm'][0][1]=str(F(bad['branches'][0]['witness_norm'][0][1])+1)
    try:check(bad)
    except AssertionError:mutation='REJECTED'
    else:raise AssertionError('Corrupted norm was accepted')
    result={'status':'PASS','arithmetic':'stdlib Fraction polynomial ring',
            'checks':checks,'wrong_norm_mutation':mutation,
            'scope':'All-parameter S8 symmetric filling obstruction under the '
                    'rank, Coxeter-weight, acuteness and strict-interior hypotheses.'}
    OUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
