#!/usr/bin/env python3
"""All-rank short-path alternative at spectator energy five quarters."""
import argparse
import json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
import fresh_spectator_energy_below_one_20261017 as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_spectator_five_quarters_transfer_20261017'

def gram(n,edges):
    g=2*s.eye(n)
    for i,j,t in edges:g[i,j]=g[j,i]=-t
    return g

def distances(g,i):
    out={i:0};todo=[i]
    for j in todo:
        for k in range(g.rows):
            if j!=k and g[j,k] and k not in out:out[k]=out[j]+1;todo.append(k)
    return out

def run():
    old=json.loads((ROOT/f'results/{B.TAG}.json').read_text())
    assert old['status']=='PASS_ALL_RANK_SPECTATOR_ENERGY_BELOW_ONE'
    for name,h in old['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    d,N=s.symbols('d N',positive=True)
    E=d*(N-d)/N
    assert s.factor(E-2*d/(d+2)-d**2*(N-d-2)/(N*(d+2)))==0
    assert s.factor(2*d/(d+2)-s.Rational(4,3)-2*(d-4)/(3*(d+2)))==0
    ar=[]
    for n in range(5,33):
        J=s.Matrix(n,n,lambda i,j:s.Rational((min(i,j)+1)*(n-max(i,j)),n+1))
        g=gram(n,[(i,i+1,1)for i in range(n-1)])
        assert g*J==s.eye(n)
        count=0;least=None
        for i,j in combinations(range(n),2):
            if j-i<4:continue
            a,b,c=J[i,i],J[j,j],J[i,j]
            assert a-c>0 and b-c>0
            e=a+b-2*c
            assert e==s.Rational((j-i)*(n+1-j+i),n+1) and e>=s.Rational(4,3)
            least=e if least is None else min(e,least);count+=1
        ar.append(dict(rank=n,far_pairs=count,minimum=str(least)))
    classical=[]
    for kind,start in(('B',2),('D',4)):
        for n in range(start,25):
            if kind=='B':
                edges=[(i,i+1,s.sqrt(2)if i==n-2 else 1)for i in range(n-1)]
                J=s.Matrix(n,n,lambda i,j:s.Rational(n,2)if i==j==n-1 else
                    (min(i,j)+1)/s.sqrt(2)if max(i,j)==n-1 else min(i,j)+1)
            else:
                edges=[(i,i+1,1)for i in range(n-3)]+[(n-3,n-2,1),(n-3,n-1,1)]
                J=s.Matrix(n,n,lambda i,j:s.Rational(n if i==j else n-2,4)if min(i,j)>=n-2 else
                    s.Rational(min(i,j)+1,2)if max(i,j)>=n-2 else min(i,j)+1)
            g=gram(n,edges);assert (g*J).applyfunc(s.simplify)==s.eye(n)
            allowed=[i for i in range(n)if J[i,i]<2]
            assert allowed==([0]+([n-1]if n<=3 else[])if kind=='B'else[0]+([n-2,n-1]if n<=7 else[]))
            far=[]
            for i,j in combinations(allowed,2):
                dist=distances(g,i)[j]
                if dist<4:continue
                a,b,c=J[i,i],J[j,j],J[i,j]
                assert a-c>=0 and b-c>=0
                e=s.simplify(a+b-2*c)
                assert kind=='D' and n in(6,7) and e==s.Rational(n,4) and e>=s.Rational(3,2)
                far.append(dict(marks=[i,j],distance=dist,unit_energy=str(e),
                    row_sums=list(map(str,(a-c,b-c)))))
            classical.append(dict(kind=kind,rank=n,diagonals=list(map(str,J.diagonal())),allowed=allowed,far_pairs=far))
    exceptional=[]
    for n in(6,7,8):
        edges=[];node=1
        for length in(n-4,2,1):
            previous=0
            for _ in range(length):edges.append((previous,node,1));previous=node;node+=1
        assert node==n
        g=gram(n,edges);J=g.inv();allowed=[i for i in range(n)if J[i,i]<2]
        assert len(allowed)=={6:2,7:1,8:0}[n]
        pairs=[]
        for i,j in combinations(allowed,2):
            a,b,c=J[i,i],J[j,j],J[i,j]
            assert n==6 and distances(g,i)[j]==4
            assert a==b==s.Rational(4,3) and c==s.Rational(2,3)
            pairs.append(dict(marks=[i,j],distance=4,unit_energy='4/3',row_sums=['2/3','2/3']))
        exceptional.append(dict(kind=f'E{n}',gram=[list(map(str,row))for row in g.tolist()],
            inverse=[list(map(str,row))for row in J.tolist()],allowed=allowed,pairs=pairs))
    # Once no opposite contacts lie in one component, all component energies
    # are one-sided. At total5/4 each of two components is in[1/2,3/4].
    candidates=[s.Rational(r,r+1)for r in range(1,4)]
    rank_pairs=[(r,t)for r in range(1,4)for t in range(r,4)
                if s.Rational(r,r+1)+s.Rational(t,t+1)==s.Rational(5,4)]
    assert rank_pairs==[(1,3)]
    paths=[Path(__file__)]+[ROOT/f'{folder}/{B.TAG}.{ext}'for folder,ext in
        (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ALL_RANK_SPECTATOR_FIVE_QUARTERS_ALTERNATIVE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        hypotheses=old['hypotheses']|dict(energy='D=5/4',both_supports='nonempty'),
        A_far_pair_identities=ar,BD_inverse_identities=classical,E_exact_cases=exceptional,
        same_component_far_lower_bound='4/3',
        no_mixed_component_energies=list(map(str,candidates)),separated_rank_pairs=rank_pairs,
        theorem='Either some x-contact and y-contact have distance<=3 in a whole spectator component, giving a new-black path of length<=5; or the only contacted components are whole A1 and A3, with one unit leaf contact on each and opposite contact sides.',
        arbitrary_rank_proof='Strict augmented forests permit at most one contact per component on each side. The all-rank A/B/D inverse formulas and finite E table exclude distance>=4 at energy<=5/4. F4,H3,H4,I2 have diameter<=3. Without mixed components, the one-sided spectrum below1 forces A1+A3.',
        scope='All-rank local spectator alternative. No finite coupon, complete face transfer, ambient dimension bound or construction is asserted.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
