#!/usr/bin/env python3
"""Exact derived density bounds; geometric/external premises are in the report."""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_right_angled_mixed_five_face_density_20261018.json'
PINS={
 'results/right_angled_five_face_carrier_20261018.md':'ce90d9a3fa3c217e6811163a48f84d8ccadaafbe9c99eb4d62ef86de36dc1e0d',
 'results/ideal_ra4_second_volume_equality_20261018.md':'4d227e9f7f4577f0a61b30a3ebf2790c763feeb4f1acddf7799a83f96fc1d7e2'}
def nikulin(n,k,l):
    a,b=n//2,(n+1)//2
    assert 1<=l<k<=b
    return Q(comb(n-l,n-k)*(comb(a,l)+comb(b,l)),comb(a,k)+comb(b,k))
def build():
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    rows=[]
    for n in range(9,13):
        u=nikulin(n,5,4);edge_u=nikulin(n,5,1)
        fraction=(u-16)/40;c=1-fraction;C=comb(n,5)
        # Global edge identity 2e=n*v_fin+2^(n-1)*v_ideal.
        multiplier=Q(comb(n-1,4),1)/edge_u
        aa=multiplier*Q(n,2);bb=multiplier*2**(n-2)
        ratio=c*bb/(C-c*aa)
        assert 0<fraction<1 and C-c*aa>0
        rows.append({'n':n,'average_four_facets_in_five_faces_strict_upper':str(u),
          'average_edges_in_five_faces_strict_upper':str(edge_u),
          'ideal_five_face_fraction_strict_upper':str(fraction),
          'nonideal_five_face_fraction_strict_lower':str(c),
          'finite_vertex_five_face_incidence':C,
          'five_face_count_strict_lower_coefficients':[str(aa),str(bb)],
          'finite_to_ideal_vertex_ratio_strict_lower':str(ratio)})
    assert [r['ideal_five_face_fraction_strict_upper'] for r in rows]==['7/20','7/20','1/10','1/10']
    assert rows[-1]['finite_to_ideal_vertex_ratio_strict_lower']=='512/437'
    return {'status':'NECESSARY_MIXED_FIVE_FACE_DENSITIES_NOT_DIMENSION_EXCLUSION',
      'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'sources':PINS,
      'external_premises':{
        'ideal_RA4_min_facets':24,'equality_only_regular_24cell':True,
        'general_finite_volume_RA5_min_facets':16,
        'nikulin_khovanskii_scope':'simple at edges; used only l=1,4 k=5 n=9..12',
        'primary_formula_reference':'https://arxiv.org/html/2303.09533v1#S6',
        'reported_Dufour_constants':'https://arxiv.org/html/2603.29925v1#S2'},
      'internal_positive_carrier_minima':{'facets':25,'delta':5,'vertices':30,'facets_plus_vertices':55},
      'ideal_RA5_min_facets':56,'density_rows':rows,
      'new_dimension_bound':False,'optimality_asserted':False,
      'Dufour_original_full_proof_reviewed':False,
      'scope':'exact arithmetic consequences of explicitly stated universal and published premises'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=build()
    if a.verify:assert json.loads(OUT.read_text())==d
    else:
        assert not OUT.exists();OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({'status':'PASS','verify':a.verify,'ideal5_min_facets':56,
      'fractions':[x['ideal_five_face_fraction_strict_upper'] for x in d['density_rows']],
      'dimension12_vertex_ratio':d['density_rows'][-1]['finite_to_ideal_vertex_ratio_strict_lower']}))
