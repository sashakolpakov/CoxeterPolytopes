#!/usr/bin/env python3
"""Produce symbolic identities for an all-parameter S8 filling obstruction."""
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/enumerations/fresh_h9_symmetric_filling_proof_20261013.json'
a,b,w,p,q = s.symbols('a b w p q')
VARIABLES = (a,b,w,p,q)


def encode(f):
    return [[list(m),str(c)] for m,c in s.Poly(s.expand(f),*VARIABLES).terms() if c]


def gram():
    g = 2*s.eye(17)
    for i in range(8):
        u,v = 2*i,2*i+1
        g[u,v] = g[v,u] = -w
        g[16,u] = g[u,16] = -p
        g[16,v] = g[v,16] = -q
        for j in range(i):
            h,k = 2*j,2*j+1
            g[u,h] = g[h,u] = g[v,k] = g[k,v] = -a
            g[u,k] = g[k,u] = g[v,h] = g[h,v] = -b
    return g


def main():
    G = gram(); S=s.Matrix([1]*16+[0]); delta=2+a; total=a+b
    contrast=s.Matrix([t for i in range(8) for t in ((1,-1) if i<4 else (-1,1))]+[0])
    C=s.zeros(17,14)
    for i in range(7):
        C[2*i,i]=1; C[14,i]=-1
        C[2*i+1,7+i]=1; C[15,7+i]=-1
    K=s.ones(7)+s.eye(7)
    B=s.Matrix([[2+a,b-w],[b-w,2+a]])
    assert C.T*G*C == s.kronecker_product(B,K)
    rows=[]
    for branch in ('plus','minus'):
        matched=2+a+b if branch=='plus' else b-2-a
        H=G.subs(w,matched); lp=-8*total if branch=='plus' else 4-6*a-8*b
        if branch=='plus':
            Y=delta*S+4*total*contrast
        else:
            last=s.zeros(17,1);last[14]=last[15]=1
            Y=total*S+lp*last
        pairs=H*Y
        rows.append({'branch':branch,'matched_weight':encode(matched),
                     'witness_coefficients':[encode(x)for x in Y],
                     'wall_pairings':[encode(x)for x in pairs],
                     'witness_norm':encode((Y.T*H*Y)[0]),
                     'sum_norm':encode((S.T*H*S)[0]),
                     'sum_witness_pairing':encode((S.T*H*Y)[0])})
    data={'schema':'S8-H9-continuous-filling-obstruction-v1',
          'variables':[str(x)for x in VARIABLES],
          'standard_block_determinant':encode(B.det()),
          'contrast_rank':14,'ambient_rank_bound':10,
          'coxeter_orbit_weight_condition':'a,b each either 0 or at least 1',
          'branches':rows,
          'scope':'Sixteen S8 leaf normals of square2, with weights a within '
                  'each species, b between different pairs, w at matched pairs. '
                  'Any extra walls have constant pairings -p,-q with the two '
                  'species, separately for each extra wall. All weights '
                  'nonnegative; strict timelike interior; Lorentz rank10. '
                  'No compact Coxeter completion. Every changed branch except '
                  'a=b=0,w=2 has an explicit infinite-volume escape.'}
    OUT.write_text(json.dumps(data,indent=2)+'\n')
    print({'branches':2,'standard_block':s.factor(B.det()),'output':str(OUT.relative_to(ROOT))})


if __name__=='__main__':main()
