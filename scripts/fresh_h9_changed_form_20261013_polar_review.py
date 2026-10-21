#!/usr/bin/env python3
"""Small independent derivation of the two concrete final cap formulas."""
import argparse
import json
from pathlib import Path
from coxeter_polytopes.exact.golden import Golden as K
from scripts.fresh_h9_quartic_completion_20261013_verify import load_old,matrix

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_h9_changed_form_20261013_polar_review.json'

def run(verify=False):
    base=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_3_6_20261013.json').read_text())
    I=tuple(tuple(K.from_json(x)for x in row)for row in base['seed_gram_inverse'])
    phi=K(0,1);delta2=2+phi
    coefficients=(I[4][4]+2*I[4][8]+I[8][8],2*(I[4][8]+I[8][8]),I[8][8]-2/delta2)
    assert coefficients==(K(-20,'25/2'),K(-25,15),K('-42/5','24/5'))
    A,B,C=coefficients;D=B*B-4*A*C
    assert D==K(-62,39)and D>0
    assert -B/(2*A)==phi and 1/(2*A)==K('3/5',1)
    old=load_old(matrix(base['seed_gram']));polar=old[348]
    a=tuple(K(*x)for x in polar['a']);b=tuple(K(*x)for x in polar['b']);q=K(*polar['q'])
    assert all(not x for x in b)and a[4]+a[8]==0
    factor=-a[8]
    assert factor<0 and delta2*factor**2==q/2
    assert I[6][9]==0
    d9=(2-2*I[6][6])/I[9][9]
    assert d9==12+20*phi and d9>4 and 2*I[6][6]+d9*I[9][9]==2
    out={'status':'PASS_SOURCE_REVIEW_AND_INDEPENDENT_SCALAR_IDENTITIES',
        'polar_norm_quadratic_after_division_by_2_plus_phi':[x.to_json()for x in coefficients],
        'polar_discriminant':D.to_json(),'positive_parameter_form':'phi + (3/5+phi)*sqrt(-62+39phi)',
        'polar_pair_factor_times_sqrt_2_plus_phi':factor.to_json(),
        'polar_pi3_normalized_square':True,'opposite_radicand':d9.to_json(),'opposite_norm':2,
        'tower_embedding_review':'Old quadratic coefficients and radicands are both transported into the common positive-radical tower before pairing. Dependent radical extensions are handled by exact recursive signs without division in the larger compositum.',
        'scope':'Source review of both concrete cap producers and their explicit compositum embeddings, plus a separate Golden-arithmetic derivation of the two norm identities and polar pi/3 angle. This is not a second replay of all 44 or 45 pairs or the volume certificates.'}
    if verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
