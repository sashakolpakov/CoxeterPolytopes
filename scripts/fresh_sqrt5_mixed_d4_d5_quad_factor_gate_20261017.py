#!/usr/bin/env python3
"""Exact factor gate for free opposite products in mixed D4/D5 quadrilaterals."""
import argparse, hashlib, json
from collections import Counter
from itertools import combinations
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_mixed_d4_d5_quad_factor_gate_20261017'
INPUT='results/fresh_sqrt5_mixed_d4_d5_quad_patterns_20261017.json'
u,v=s.symbols('u v')

def decode(data):
    return s.Poly(sum(s.Rational(c)*u**a*v**b for (a,b),c in data),u,v)

def encode(poly):
    return [[list(m),str(c)] for m,c in sorted(s.Poly(poly,u,v).terms())]

def run():
    data=json.loads((ROOT/INPUT).read_text())
    for p,h in data['sources'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    targets=[]; roots=set(); factor_count=0
    for target in data['targets']:
        records=[]; retained=set(); discarded=set()
        for rec in target['records']:
            row=rec['row']; n=row['nodes']; gram=2*s.eye(n)
            for pair,label in zip(combinations(range(n),2),row['labels']):
                i,j=pair
                if pair==(0,2): q=u; assert label==0
                elif pair==(1,3): q=v; assert label==0
                else: assert label in(2,3); q=2-label
                gram[i,j]=gram[j,i]=q
            # Reconstruct the Schur matrix directly, without producer arithmetic.
            white=list(range(4,n)); black=list(range(4))
            W=gram.extract(white,white); D=gram.extract(black,white)
            schur=gram.extract(black,black)-D*W.inv()*D.T
            for i in black:
                for j in black:
                    assert s.expand(schur[i,j]-decode(rec['symbolic']['residual'][i][j]).as_expr())==0
            determinant=s.Poly(schur.det(method='domain-ge'),u,v)
            assert determinant==decode(rec['symbolic']['determinant']) and not determinant.is_zero
            coefficient,factors=s.factor_list(determinant.as_expr(),u,v)
            rebuilt=coefficient; branches=[]; survives=False
            for factor,multiplicity in factors:
                p=s.Poly(factor,u,v); rebuilt*=factor**multiplicity; factor_count+=1
                du,dv=p.degree(u),p.degree(v)
                assert max(du,dv)>0
                branch=dict(polynomial=encode(p),multiplicity=multiplicity,degrees=[du,dv])
                if not du or not dv:
                    variable=u if du else v
                    q=s.Poly(factor,variable); assert q.degree()==1
                    root=-q.nth(0)/q.nth(1); assert root.is_Rational
                    roots.add(root)
                    branch.update(status='EXCLUDED_RATIONAL_EDGE_GRAM',variable=str(variable),root=str(root))
                else:
                    survives=True
                    branch.update(status='COUPLED_PARAMETER_BRANCH_UNRESOLVED')
                branches.append(branch)
            assert s.Poly(s.expand(rebuilt),u,v)==determinant
            (retained if survives else discarded).add(rec['index'])
            records.append(dict(index=rec['index'],coefficient=str(coefficient),branches=branches,
                status='HAS_UNRESOLVED_COUPLED_BRANCH' if survives else 'EXCLUDED_ALL_BRANCHES_RATIONAL_EDGE'))
        assert retained.isdisjoint(discarded) and len(retained|discarded)==len(target['records'])
        maps=[x for x in target['maps'] if x['target'] in retained]
        corners=sorted({tuple(x['ports']) for x in maps})
        targets.append(dict(name=target['name'],forest=target['forest'],dimension=target['dimension'],
            original_colored_patterns=target['colored_patterns'],retained_colored_patterns=len(retained),
            retained_indices=sorted(retained),excluded_indices=sorted(discarded),retained_literal_patterns=len(maps),
            retained_possible_corners=[list(p) for p in corners],retained_corner_count=len(corners),records=records))
    paths=[f'scripts/{TAG}.py',INPUT,'scripts/fresh_sqrt5_mixed_d4_d5_quad_patterns_20261017.py']
    return dict(status='PASS_EXACT_MIXED_D_QUADRILATERAL_RATIONAL_FACTOR_GATE',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        rational_roots=sorted(map(str,roots)),factor_count=factor_count,targets=targets,
        rational_edge_lemma='At an actual compact Coxeter edge, the d finite vertex normals plus the new endpoint normal have full Lorentz rank d+1. Here all products except that edge framing product are rational (0,-1,2). If the framing product is rational, this entire normalized Gram is rational. Its connected blocks have a rational norm-two gauge; under the nonidentity embedding of Q(sqrt5), any sign gauge is diagonally congruent to that same rational Gram. Its negative direction contradicts arithmetic conjugate positivity.',
        scope='All exact rational univariate determinant factors are excluded for actual arithmetic compact quadrilaterals with the stated four equal forests. Coupled factors remain only necessary rank equations, not certified solutions. No remaining branch has been excluded by an arbitrary parameter bound. No mixed incidence inequality or dimension bound is asserted.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if a.verify: assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],[(x['name'],x['retained_colored_patterns'],x['retained_corner_count']) for x in data['targets']],flush=True)
