#!/usr/bin/env python3
"""Exact finite diagram checks for the persistent odd-neighbor doubling gate."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_dimension_alternative_routes_20261014_doubling_gate'
BASE='results/fresh_dimension_alternative_routes_20261014.json'
def norm(x):return s.expand(x)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    base=json.loads((ROOT/BASE).read_text())
    for p,h in base['sources'].items():assert sha(p)==h
    cos={2:s.Integer(0),3:s.Rational(1,2),4:s.sqrt(2)/2,5:(1+s.sqrt(5))/4}
    G=s.eye(10)
    for a,b,m in base['source_nonright_edges']:G[a-1,b-1]=G[b-1,a-1]=-cos[m]
    z=2*s.sqrt(2)+s.sqrt(10);w=17+8*s.sqrt(5)
    for (a,b),q in { (6,7):z,(7,8):w,(8,9):z }.items():G[a-1,b-1]=G[b-1,a-1]=-q
    def certify(M):
        n=M.rows;relations={};counts=Counter()
        assert all(norm(M[i,i]-1)==0 for i in range(n))
        for i,j in combinations(range(n),2):
            x=norm(M[i,j]);m=next((k for k,c in cos.items() if norm(x+c)==0),None)
            if m is None:
                assert (-x-1).is_positive,(i,j,x);m='dotted'
            relations[i,j]=m;counts[str(m)]+=1
        odd={i for (a,b),m in relations.items() if isinstance(m,int) and m>=3 and m%2 for i in (a,b)}
        witnesses=[]
        for i in range(n):
            if i in odd:
                j=next(j for j in range(n) if j!=i and isinstance(relations[tuple(sorted([i,j]))],int) and relations[tuple(sorted([i,j]))]>=3 and relations[tuple(sorted([i,j]))]%2)
                witnesses.append(dict(facet=i,kind='odd',neighbor=j,label=relations[tuple(sorted([i,j]))]))
            else:
                j=next(j for j in odd if isinstance(relations[tuple(sorted([i,j]))],int) and relations[tuple(sorted([i,j]))]>=3)
                witnesses.append(dict(facet=i,kind='nonright_to_odd',neighbor=j,label=relations[tuple(sorted([i,j]))]))
        eligible=[i for i in range(n) if not any(isinstance(m,int) and m%2 for p,m in relations.items() if i in p)]
        return dict(facets=n,pair_counts=dict(sorted(counts.items())),odd_facets=sorted(odd),doubling_facets=eligible,witnesses=witnesses)
    original=certify(G);assert original['doubling_facets']==[6,7]
    doubles=[]
    for k in [6,7]:
        roots=[(i,0) for i in range(10) if i!=k]+[(i,1) for i in range(10) if i!=k and G[i,k]!=0]
        M=s.zeros(len(roots))
        for a,(i,si) in enumerate(roots):
            for b,(j,sj) in enumerate(roots):
                M[a,b]=norm(G[i,j] if si==sj else G[i,j]-2*G[i,k]*G[j,k])
        result=certify(M);assert result['facets']==12
        doubles.append(dict(doubling_original_facet=k+1,facets=[dict(source=i+1,reflected=bool(v)) for i,v in roots],
                            exact_gram=[[str(x) for x in row] for row in M.tolist()],checks=result))
    paths=[BASE,f'scripts/{TAG}.py']
    return dict(status='PASS_NO_ALL_RIGHT_FACET_AFTER_ANY_FINITE_LEGITIMATE_DOUBLINGS',sources={p:sha(p) for p in paths},
        initial_diagram=original,first_doubles=doubles,
        invariant='Every facet has an odd finite-angle neighbor, or a nonright finite-angle neighbor that itself has an odd finite-angle neighbor.',
        theorem='This invariant is preserved by doubling across any wall all of whose meeting labels are even. No all-right facet can therefore occur after a finite sequence.',
        scope='Known compact H6 ten-facet source and any finite sequence of legitimate single-wall Coxeter doublings. No claim about general gluing, non-doubling reflection subgroups, rolling, source angle surgery, or all compact H9 constructions.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
