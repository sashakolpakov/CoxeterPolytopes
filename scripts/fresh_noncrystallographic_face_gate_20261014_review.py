#!/usr/bin/env python3
"""Separate SymPy check of the weighted I2(5) projection identities."""
import argparse
from hashlib import sha256
from pathlib import Path
import json
import sympy as S

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_noncrystallographic_face_gate_20261014_review'
def run():
    ph=(1+S.sqrt(5))/2;k=S.simplify(4-ph**2)
    G=S.Matrix([[4,-ph],[-ph,1]]);I=S.simplify(G.inv())
    def projection_square(v):return S.simplify(1-(S.Matrix([v])*I*S.Matrix(v))[0])
    squares=[projection_square([-1,0]),projection_square([-S.sqrt(2),0]),projection_square([-ph,0]),projection_square([-ph,-1])]
    assert squares[0]>0 and all(x<0 for x in squares[1:])
    first=S.simplify(I[0,0]/squares[0]);opp=S.simplify((1+I[0,0])/squares[0])
    assert S.simplify(first-ph**2)==0 and S.simplify(opp-(3+2*ph))==0
    assert first>1 and opp>1
    assert 0<squares[0]<1
    # Each unit wall really gives a facet of the I2(5) face: the
    # three-root Gram is PD and its inverse is entrywise nonnegative.
    H=S.Matrix([[4,-ph,-1],[-ph,1,0],[-1,0,1]])
    assert all(S.simplify(H[:n,:n].det())>0 for n in (1,2,3))
    HI=S.simplify(H.inv());assert all(v>=0 for v in HI)
    assert S.simplify(ph**2-4)<0
    assert S.simplify(20+8*ph-4)>0
    srcs=[f'scripts/{TAG}.py','scripts/fresh_noncrystallographic_face_gate_20261014.py','results/fresh_noncrystallographic_face_gate_20261014.json']
    return dict(status='PASS_INDEPENDENT_PROJECTION_IDENTITIES_AND_THEOREM_SOURCE_REVIEW',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in srcs},
        projected_squares=list(map(str,squares)),projected_unit_cosines=[str(first),str(opp)],
        actual_unit_facet_gram_pivotal_determinants=[str(S.simplify(H[:n,:n].det()))for n in (1,2,3)],
        mathematical_dependencies=[
            dict(source='Allcock, Infinitely many hyperbolic Coxeter groups through dimension19',url='https://arxiv.org/pdf/0903.0138',location='Theorems2.1/2.2, printedpp743–745',checked='I2(m), m outside2,3,4,6, has no An/D5 component; its face is Coxeter.'),
            dict(source='Potyagailo–Vinberg, On right-angled reflection groups in hyperbolic spaces',url='https://ems.press/content/serial-article-files/43032',location='§2, printedp65',checked='Compact right-angled hyperbolic polytopes have dimension at most4.')],
        scope='Separate SymPy arithmetic for all seed projection identities and actual unit-facet positivity. Source/proof review of the compact-face cusp argument and the restricted right-angled completion obstruction. The published general theorems are dependencies, not computationally re-proved; no unrestricted H20 nonexistence claim.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run();dest=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(dest.read_text())==out
    else:dest.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
