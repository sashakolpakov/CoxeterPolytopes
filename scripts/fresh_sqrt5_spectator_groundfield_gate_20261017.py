#!/usr/bin/env python3
"""Independent exact membership/sign/integrality gate; no producer imports."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sympy as s
from sympy.polys.polyerrors import IsomorphismFailed
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_spectator_groundfield_gate_20261017'
OLD='fresh_sqrt5_multispectator_cap_capacity_20261017'
TRIPLE='fresh_sqrt5_triple_a1_cap_capacity_20261017'
SQ=s.sqrt(5)
X=s.Symbol('X')

def inK(x):
    try:
        a=s.to_number_field(x,SQ).as_expr()
        assert s.simplify(x-a)==0
        return True,str(a)
    except IsomorphismFailed:
        return False,None

def run():
    d=json.loads((ROOT/'results'/f'{OLD}.json').read_text())
    t=json.loads((ROOT/'results'/f'{TRIPLE}.json').read_text())
    assert t['arithmetic_ground_field']=='Q(sqrt5)'
    for obj in (d,t):
        for p,h in obj['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    out=[];changes=[];summary=[]
    cases=[('two',d,s.Rational(3,2),[(0,0),(1,0),(0,1),(0,2)],['disjoint','common_head','chain','swapped']),
           ('three',t,s.Integer(1),[(0,0),(1,0),(2,0),(0,1),(1,1),(0,2),(1,2)],
            ['j0_c0','j1_c0','j2_c0','j0_c1','j1_c1','j0_c2','j1_c2'])]
    for kind,obj,constant,pairs,names in cases:
        assert len(obj['tau_pair_tests'])==(11 if kind=='two' else 5)
        for z in obj['tau_pair_tests']:
            a,b=map(s.Rational,z['tau']);tau=a+b*SQ
            rootK,rootexpr=inK(s.sqrt(tau))
            h=constant-tau/2
            mask=[];oldmask=[]
            for (j,c),name in zip(pairs,names):
                H=h+s.Rational(j,2);C=s.Rational(c,2)
                p=H+C*s.sqrt(tau);w=s.expand(p*p)
                sign=int(s.sign(p));assert sign in (-1,0,1)
                poly=s.Poly(s.minimal_polynomial(w,X),X).monic()
                integral=all(a.q==1 for a in poly.all_coeffs())
                belongs,expr=inK(w)
                # Independent field embedding agrees with the elementary radical criterion.
                assert belongs==(rootK or H*C==0)
                weak=sign<=0 and integral;strong=weak and belongs
                saved=z['pairs'][name]
                assert saved['product_sign']==sign
                assert saved['two_cycle_integral']==integral
                if kind=='two':assert saved['allowed_by_necessary_tests']==weak
                else:
                    assert saved['two_cycle_in_Qsqrt5']==belongs
                    assert saved['passes_acuteness_and_integrality']==weak
                    assert saved['allowed_by_necessary_tests']==strong
                mask.append(strong);oldmask.append(weak)
                row=dict(family=kind,tau=z['tau'],pair=name,j=j,c=c,
                    square_root_in_K=rootK,square_root_expression=rootexpr,
                    product_sign=sign,monic_minimal_polynomial_descending=list(map(str,poly.all_coeffs())),
                    two_cycle_integral=integral,two_cycle_in_K=belongs,two_cycle_expression=expr,
                    allowed_weak=weak,allowed_fixed_groundfield=strong)
                out.append(row)
                if weak and not strong:changes.append(row)
            summary.append(dict(family=kind,tau=z['tau'],weak_pattern=oldmask,fixed_groundfield_pattern=mask))
    assert len(out)==79
    assert sum(r['family']=='three' for r in changes)==2
    paths=[Path(__file__),ROOT/'scripts'/f'{OLD}.py',ROOT/'results'/f'{OLD}.json',
           ROOT/'scripts'/f'{TRIPLE}.py',ROOT/'results'/f'{TRIPLE}.json']
    return dict(status='PASS_INDEPENDENT_FIXED_GROUND_FIELD_PAIR_GATE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        ground_field='Q(sqrt5)',normalization='All compared cap roots have square2, so the Vinberg two-cycle is the squared native product.',
        criterion='For H,C in K and positive tau in K, (H+C*sqrt(tau))^2 belongs to K iff sqrt(tau) belongs to K or H*C=0.',
        counts=dict(total_pair_tests=len(out),two_spectator_tests=44,three_spectator_tests=35,
                    newly_forbidden=len(changes),two_spectator_new=sum(z['family']=='two'for z in changes),
                    three_spectator_new=sum(z['family']=='three'for z in changes)),
        pair_checks=out,pattern_summary=summary,newly_forbidden=changes,
        scope='Independent SymPy field embeddings, minimal polynomials and exact signs; no producer arithmetic imports, no duplicate physical cap or colored-orbit reconstruction, no new optimization or capacity-model assertion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
