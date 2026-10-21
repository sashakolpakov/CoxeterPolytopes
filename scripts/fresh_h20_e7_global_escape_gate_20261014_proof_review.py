#!/usr/bin/env python3
"""Narrow independent symbolic proof checks; no ambient reconstruction.

The separate arithmetic review reconstructs the target coefficients. Here we
check the rational corner, saved strict signs with another Fraction engine,
and the exact symbolic substitutions/minimizer identities.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as sp
from fresh_h20_e7_face_20261014_review import sign

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_global_escape_gate_20261014_proof_review'
SOURCE='results/fresh_h20_e7_global_escape_gate_20261014.json'
FACE='results/fresh_h20_e7_face_20261014.json'
FACE_REVIEW='results/fresh_h20_e7_face_20261014_review.json'


def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def load(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d.get('sources',{}).items():assert sha(f)==h,f
    return d


def sgn(x):return sign(tuple(map(F,x)))


def run():
    d,f,review=[load(p) for p in [SOURCE,FACE,FACE_REVIEW]]
    B=sp.Matrix([[F(x) for x in row] for row in f['gram']])
    G=sp.Matrix([[F(x) for x in row] for row in f['face_gram']])
    corner=sp.diag(B,G);_,D0=corner.LDLdecomposition(hermitian=False)
    assert all(D0[i,i]>0 for i in range(20))
    inverse=corner.inv();assert all(x>=0 for x in inverse)
    assert inverse.tolist()==[[F(x) for x in row] for row in d['finite_corner']['inverse']]
    assert d['finite_corner']['walls']==f['E7_support']+f['face_walls']
    J=G.inv();assert J[0,0]==sp.Rational(11,12) and J[11,11]==sp.Rational(2,3)
    assert all(J[i,11]==0 for i in range(11))
    assert review['intrinsic_form_determinant']==-72
    assert review['unique_positive_residual_nonright_profile']==dict(wall=33,label=3,residual_fraction='1/4')
    tc=d['target_coefficients'];assert all(sgn(x)>0 for x in tc['a']+[tc['B'],tc['w']])
    right=d['right_lift'];assert all(F(v)>0 for v in right['multipliers'].values())
    assert len(right['residual_coefficients'])==14 and all(sgn(x)<0 for x in right['residual_coefficients'])
    n=d['nonright_lift'];assert len(n['mixed_coefficients'])==13
    assert all(sgn(x)>0 for x in n['mixed_coefficients']) and sgn(n['mixed_rhs'])>0
    def sub(x,y):return [F(a)-F(b) for a,b in zip(x,y)]
    assert sgn(sub(n['mixed_coefficients'][12],n['mixed_rhs']))>0
    row=n['cut_boundary_rows']['16'];rhs=n['cut_boundary_rhs']['16']
    assert sgn(row[0])<0 and sgn(sub(row[0],rhs))>0
    assert all(sgn(row[i])>0 for i in range(1,12))
    assert all(sgn(n[k])>0 for k in ['rank_one_positive_gap','denominator','C','minimum_norm_minus_half_times_D'])
    assert len(n['positive_gradient_numerators'])==11 and all(sgn(x)>0 for x in n['positive_gradient_numerators'])
    # B*(r.oldwall) = (linear cut-boundary row minus rhs) + 2*z*(r.E).
    aa,cc,ww,BB,zz,oo,tt=sp.symbols('aa cc ww BB zz oo tt')
    escape=aa+ww-BB*tt;pair=-cc-2*zz*tt-oo
    boundary=-BB*cc-2*zz*aa-BB*oo-2*zz*ww
    assert sp.expand(BB*pair-boundary-2*zz*escape)==0
    # One free coordinate remains after d0=sqrt2,d12=0 and the other
    # A11 coordinates vanish. M denotes the positive scalar B^2.
    M,a,C,z=sp.symbols('M a C z');D=M-3*a*a
    optimum=3*a*C/D
    objective=sp.Rational(11,6)+sp.Rational(2,3)*z*z-2*(C+a*z)**2/M
    assert sp.factor(sp.diff(objective,z).subs(z,optimum))==0
    assert sp.factor(objective.subs(z,optimum)-(sp.Rational(11,6)-2*C*C/D))==0
    assert sp.factor(((sp.Rational(11,6)-2*C*C/D)-sp.Rational(1,2))*D-
                     (sp.Rational(4,3)*D-2*C*C))==0
    ai,Ji=sp.symbols('ai Ji')
    half_gradient=Ji*sp.sqrt(2)-2*ai*(C+a*optimum)/M
    assert sp.factor(D*half_gradient-(D*Ji*sp.sqrt(2)-2*ai*C))==0
    paths=[SOURCE,FACE,FACE_REVIEW,f'scripts/{TAG}.py',
           'scripts/fresh_h20_e7_global_escape_gate_20261014.py',
           'scripts/fresh_h20_e7_face_20261014_review.py']
    return dict(status='PASS_E7_GLOBAL_GATE_SOURCE_AND_SYMBOLIC_PROOF_REVIEW',
                sources={p:sha(p) for p in paths},
                rational_corner_dimension=20,nonnegative_corner_inverse_entries=400,
                right_residual_strict_signs=14,mixed_positive_coefficients=13,
                constrained_positive_gradient_signs=11,symbolic_identities_checked=4,
                time_orientation='Any strict future center and nonnegative full-corner inverse force t>0 for every valid acute-compatible addition; the old finite apex is automatically retained.',
                endpoint_strengthening='The same necessary inequalities hold under r.E>=0, so the strict minimum gap excludes equality as well: all walls within the stated two-lift scope have r.E<0.',
                scope='Source/proof review with separate rational and symbolic algebra plus exact signs of saved coefficient rows. Ambient normals and target coefficients are not independently reconstructed here; that arithmetic is covered by the separate coefficient replay. No arbitrary-completion obstruction is claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
