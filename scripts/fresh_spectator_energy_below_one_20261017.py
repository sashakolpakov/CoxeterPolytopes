#!/usr/bin/env python3
"""Exact identities for the all-rank strict spectator spectrum D<1."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_spectator_energy_below_one_20261017'
BASE='fresh_small_spectator_difference_gap_20261017'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def matrix(n,edges,last_scale=False):
    g=[[F(2 if i==j else 0)for j in range(n)]for i in range(n)]
    if last_scale:g[-1][-1]=4
    for i,j,w in edges:g[i][j]=g[j][i]=-F(w)
    return g


def tree_pivots(g):
    diagonal=[g[i][i]for i in range(len(g))]
    adj=[{j for j,x in enumerate(row)if x and i!=j}for i,row in enumerate(g)]
    active=set(range(len(g)));out=[]
    while active:
        leaf=min(i for i in active if len(adj[i]&active)<=1)
        p=diagonal[leaf];assert p>0
        out.append(str(p))
        for j in adj[leaf]&active:diagonal[j]-=g[leaf][j]**2/p
        active.remove(leaf)
    return out


def check_column(g,x,k):
    assert all(sum((a*b for a,b in zip(row,x)),F())==int(i==k)
        for i,row in enumerate(g))


def run():
    old=json.loads((ROOT/f'results/{BASE}.json').read_text())
    assert old['status'].startswith('PASS')
    for name,h in old['sources'].items():assert digest(ROOT/name)==h,name
    ar=[];templates=[];adjacent=0
    for n in range(1,65):
        g=matrix(n,[(i,i+1,1)for i in range(n-1)])
        inv=[[F((min(i,j)+1)*(n-max(i,j)),n+1)for j in range(n)]for i in range(n)]
        # Every column is verified against the sparse native Gram, not only
        # its proposed diagonal. The proof gives the formula for arbitrary n.
        for j in range(n):
            assert all(2*inv[i][j]-(inv[i-1][j]if i else 0)
                -(inv[i+1][j]if i+1<n else 0)==int(i==j)for i in range(n))
        value=F(n,n+1);assert inv[0][0]==value
        assert all(inv[i][i]>=1 for i in range(1,n-1))
        pairs=[]
        for i in range(n-1):
            energy=inv[i][i]+inv[i+1][i+1]-2*inv[i][i+1]
            assert energy==value
            v=[inv[j][i]-inv[j][i+1]for j in range(n)]
            expected=[F(j+1,n+1)if j<=i else-F(n-j,n+1)for j in range(n)]
            assert v==expected and all(v)
            pairs.append(dict(marks=[i,i+1],separate_edge_extensions_PD=inv[i][i]<2 and inv[i+1][i+1]<2))
            adjacent+=1
        ar.append(dict(rank=n,energy=str(value),inverse_leaf_column=[str(x[0])for x in inv],
            adjacent_contact_pairs=pairs))
    for kind,start in(('B',2),('D',4)):
        for n in range(start,65):
            if kind=='B':
                edges=[(i,i+1,2 if i==n-2 else 1)for i in range(n-1)]
                g=matrix(n,edges,True)
                first=[F(1)]*(n-1)+[F(1,2)]
                second=([F(1)]+[F(2)]*(n-2)+[F(1)])if n>=3 else None
            else:
                edges=[(i,i+1,1)for i in range(n-3)]+[(n-3,n-2,1),(n-3,n-1,1)]
                g=matrix(n,edges)
                first=[F(1)]*(n-2)+[F(1,2)]*2
                second=[F(1)]+[F(2)]*(n-3)+[F(1)]*2
            pivots=tree_pivots(g);check_column(g,first,0)
            if second is not None:check_column(g,second,1)
            assert first[0]==1 and(second is None or second[1]==2)
            templates.append(dict(kind=kind,rank=n,positive_pivots=pivots,
                first_inverse_column=list(map(str,first)),second_inverse_column=None if second is None else list(map(str,second)),
                first_normalized_inverse_diagonal='1',second_normalized_inverse_diagonal='2'if second is not None else None))
    assert adjacent==2016 and len(templates)==124
    # At p=sqrt2,q=t=1 the two-node energy is2-(2/3)sqrt2.
    # Its excess over1 is positive because1^2>2*(2/3)^2.
    assert F(1)>2*F(2,3)**2
    examples=[]
    for value in(F(4,7),F(11,14),F(2,3),F(3,4),F(4,5),F(9,16),F(8,17)):
        rank=value/(1-value)
        allowed=rank.denominator==1 and rank>=1
        examples.append(dict(rational_energy=str(value),forced_rank=str(rank),allowed=allowed))
    assert [r['allowed']for r in examples]==[False,False,True,True,True,False,False]
    paths=[Path(__file__)]+[ROOT/f'{folder}/{BASE}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ALL_RANK_SPECTATOR_ENERGY_BELOW_ONE',
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        hypotheses=dict(G='positive-definite norm-two acute Coxeter Gram',contacts='native, nonnegative, disjoint supports',
            separate_augmented_Grams='strictly positive definite',energy='D<1'),
        spectrum=['0','r/(r+1) for integer r>=1','2-sqrt(2)','1-sqrt(5)/5'],
        rational_equality='One unit contact at a leaf of a whole A_r, or two unit contacts at adjacent nodes of a whole A_r; both separate extensions must still be PD. All other contacts vanish.',
        A_identity_checks=ar,BD_obstruction_checks=templates,
        nonunit_A2_energy=['2','-2/3'],nonunit_A2_excess_positive_square_difference=str(1-2*F(2,3)**2),
        rational_energy_applications=examples,counts=dict(A_ranks=64,A_adjacent_pairs=adjacent,BD_templates=len(templates)),
        arbitrary_rank_proof='Written forest, support, native-gap and B/D-principal-subtree argument. A inverse columns and B/D elimination formulas hold for all ranks; finite samples only corroborate them.',
        scope='Universal strict-edge spectator lemma only. No complete new corner-transfer theorem, ambient dimension improvement or finite-volume construction.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],flush=True)
