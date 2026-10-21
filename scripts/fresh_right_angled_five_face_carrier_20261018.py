#!/usr/bin/env python3
"""Exact identity coefficients for a universal RA five-face carrier lemma."""
from hashlib import sha256
from math import comb
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_right_angled_five_face_carrier_20261018.json'
PINS={
 'results/ideal_ra4_second_volume_20261018.md':'9dff11c375bddcf62f16f3f01ddabd7f5b5c137756cf02ebff275076e6869c25',
 'results/ideal_ra4_second_volume_equality_20261018.md':'4d227e9f7f4577f0a61b30a3ebf2790c763feeb4f1acddf7799a83f96fc1d7e2'}

def build():
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    # Variable order v,e,r,g,f,D,A,K,1. Every displayed row equals zero.
    premises=[[-8,1,0,0,0,0,0,0,0],[-24,0,3,-2,0,0,0,0,0],
      [1,-1,1,-1,1,0,0,0,-2],[-8,0,0,2,0,1,0,0,0],
      [0,0,0,-2,24,0,1,0,0],[0,0,0,0,0,-3,1,1,0]]
    targets={
      'g=3v+3f-6':([-3,1,-3,0,0,0],[-3,0,0,1,-3,0,0,0,6]),
      'D=2v-6f+12':([6,-2,6,1,0,0],[-2,0,0,0,6,1,0,0,-12]),
      'A=3D-48':([-24,8,-24,-3,1,0],[0,0,0,0,0,-3,1,0,48]),
      'K=48':([24,-8,24,3,-1,1],[0,0,0,0,0,0,0,1,-48])}
    checks=[]
    for name,(multipliers,row) in targets.items():
        assert [sum(c*p[j] for c,p in zip(multipliers,premises)) for j in range(9)]==row
        checks.append({'identity':name,'multipliers':multipliers,'row':row})
    models=[{'type':'ideal24cell','vertices':24,'facets':24,'delta':0,'kappa':0},
      {'type':'double24cell','vertices':42,'facets':38,'delta':4,'kappa':-2}]
    for m in models:
        assert m['delta']==m['vertices']-m['facets']
        assert m['kappa']==3*m['delta']-(m['facets']-24)
    dimensions=[{'dimension':n,'finite_vertex_five_face_incidence':comb(n,5),
      'four_face_five_cofaces':n-4,'weighted_carrier_coefficient':48} for n in range(5,13)]
    assert [r['finite_vertex_five_face_incidence'] for r in dimensions[-4:]]==[126,252,462,792]
    return {'status':'EXACT_IDENTITIES_FOR_UNIVERSAL_FIVE_FACE_CARRIER_LEMMA',
      'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'sources':PINS,
      'variables':['v','e','r','g','f','D','A','K','1'],'premise_rows':premises,
      'identity_checks':checks,'known_four_facet_types':models,'dimensions':dimensions,
      'conclusions':{'ideal5_sum_kappa':48,'positive_kappa_facet_required':True,
        'required_facet_volume_pi_squared_at_least':'3',
        'weighted_mixed_bound':'48*(f5-binomial(n,5)*v_finite) <= (n-4)*sum_ideal4 max(kappa,0)',
        'new_dimension_exclusion':False,'new_example':False,'arithmetic_field_restriction':False},
      'scope':'identity arithmetic only; universal face and cusp incidences proved in accompanying report'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=build()
    if a.verify:assert json.loads(OUT.read_text())==d
    else:
        assert not OUT.exists();OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({'status':'PASS','verify':a.verify,'identity_checks':4,'dimensions':8,
      'ideal5_kappa_sum':48,'new_dimension_exclusion':False}))
