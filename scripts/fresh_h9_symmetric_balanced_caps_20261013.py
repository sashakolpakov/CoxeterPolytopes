#!/usr/bin/env python3
"""Produce a new golden 18-wall seed and 70 nonsymmetric balanced caps."""
import json
from itertools import combinations
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/enumerations/fresh_h9_symmetric_balanced_caps_20261013.json'
t=s.symbols('t');P=s.Poly(t*t-t-1,t)
def red(x):return s.rem(s.Poly(s.expand(x),t),P).as_expr()
def enc(x):
    q=s.Poly(red(x),t);return [str(q.nth(0)),str(q.nth(1))]


def main():
    b=(1+2*t)/2;w=2+b;c=2+2*t
    G=2*s.eye(18);G[0,1]=G[1,0]=-c
    for i in range(8):
        u,v=2+2*i,3+2*i
        G[0,u]=G[u,0]=G[1,v]=G[v,1]=-1-t
        G[u,v]=G[v,u]=-w
        for j in range(i):
            h,k=2+2*j,3+2*j
            G[u,k]=G[k,u]=G[v,h]=G[h,v]=-b
    S=s.Matrix([0,0]+[1]*16);caps=[]
    for subset in combinations(range(8),4):
        eps=[1 if i in subset else -1 for i in range(8)]
        Z=s.Matrix([0,0]+[x for e in eps for x in (e,-e)])
        r=(2-t)*S/8+t*Z/8
        caps.append({'plus_indices':list(subset),'coefficients':[enc(x)for x in r]})
    selected=[0,69];eps=[1 if i in (0,1,4,5)else-1 for i in range(8)]
    Z=s.Matrix([0,0]+[x for e in eps for x in (e,-e)])
    witness=2*S+4*b*Z
    data={'schema':'golden-symmetric-seed-balanced-caps-v1',
          'field':'phi^2=phi+1, phi>0',
          'gram':[[enc(G[i,j])for j in range(18)]for i in range(18)],
          'seed_parameters':{'a':enc(0),'b':enc(b),'w':enc(w),
                             'p':enc(1+t),'q':enc(0),'c':enc(c)},
          'strict_center_coefficients':[enc(x)for x in S],
          'caps':caps,'selected_caps':selected,
          'escape_coefficients':[enc(x)for x in witness],
          'escape_balanced_plus_indices':[0,1,4,5],
          'scope':'Exact 20-wall Coxeter cone, infinite volume. All70 caps '
                  'individually meet the seed angle gate; only complementary '
                  'balanced subsets are mutually Coxeter compatible. '
                  'Arbitrary further mirrors remain open.'}
    OUT.write_text(json.dumps(data,indent=2)+'\n');print({'caps':len(caps),'selected':selected})


if __name__=='__main__':main()
